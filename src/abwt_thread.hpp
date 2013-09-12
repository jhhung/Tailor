#ifndef ABWT_THREAD_HPP_
#define ABWT_THREAD_HPP_

#include <vector>
#include <sstream>
#include "abwt_format.hpp"
#include "abwt_search.hpp"
#include "boost/thread.hpp"
#include "boost/thread/mutex.hpp"

template <typename T>
class ABWT_threads : private ABWT_search<T>
{
private:
	static boost::mutex _io_mutex;
private:
	static const int _poolSize = 5000 ;

	std::vector <Fastq> _queryBuffer {};
	std::stringstream _resultBuffer {};

	std::istream* _in {nullptr};
	std::ostream* _out {nullptr};

	int _minLen {18};

public:
	ABWT_threads () {}
	ABWT_threads (T& table, std::istream* in, std::ostream* out, int minLen) :
		ABWT_search<T> {table},
		_in {in},
		_out {out},
		_minLen {minLen}
	{}

	ABWT_threads (ABWT_threads&& other):
		ABWT_search<T> {other},
		_queryBuffer {std::move (other._queryBuffer)},
		_in {other._in},
		_out {other._out}
	{

	}

	ABWT_threads& operator=(const ABWT_threads&) = delete;

	void operator () () {
		while (_in->good ()) {

			{ 	/// reading from input
				boost::mutex::scoped_lock lock(_io_mutex);
				for (int i = 0 ; i < _poolSize && _in->good (); ++i)
					_queryBuffer.emplace_back (*_in);
			}

			{	/// do searching
				for (const auto & _query : _queryBuffer)
					this->start_tailing_match_Dual(_query, &_resultBuffer, _minLen);
				this->start_end_pos_ = {0,0};
			}

			{	/// writting
				boost::mutex::scoped_lock lock(_io_mutex);
				*_out << _resultBuffer.rdbuf ();
				_resultBuffer.str(std::string());
			}
			_queryBuffer.clear ();
		}
	}
};

template <typename T>
boost::mutex ABWT_threads<T>::_io_mutex;


#endif /* ABWT_THREAD_HPP_ */
