/*
# Tailor, a BWT-based aligner for non-templated RNA tailing
# Copyright (C) 2014 Min-Te Chou, Bo W Han, Jui-Hung Hung
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef ABWT_THREAD_HPP_
#define ABWT_THREAD_HPP_

#include <vector>
#include <sstream>
#include "abwt_format.hpp"
#include "abwt_search.hpp"
#include "boost/thread.hpp"
#include "boost/thread/mutex.hpp"
#include <type_traits>

template <typename T>
class ABWT_threads : private ABWT_search<T>
{
private:
	static boost::mutex _io_mutex;
private:
	static const int _poolSize = 5000 ;

	std::vector <Fastq> _queryBuffer {};
	std::stringstream _resultBuffer {};
	std::vector<INTTYPE> _resultBuffer2 {};

	std::istream* _in {nullptr};
	std::ostream* _out {nullptr};

	int _minLen;
	bool _allowMismatch;

public:
	ABWT_threads () {}
	ABWT_threads (T& table, std::istream* in, std::ostream* out, int minLen, bool allow_mismatch) :
		ABWT_search<T> {table},
		_in {in},
		_out {out},
		_minLen {minLen},
		_allowMismatch {allow_mismatch}

	{ }

	ABWT_threads (ABWT_threads&& other):
		ABWT_search<T> {other},
		_queryBuffer {std::move (other._queryBuffer)},
		_in {other._in},
		_out {other._out},
		_minLen {other._minLen},
		_allowMismatch {other._allowMismatch}
	{ }

	ABWT_threads& operator=(const ABWT_threads&) = delete;

	void operator () () {
		while (1) {
			{ 	/// reading from input
				boost::mutex::scoped_lock lock(this->_io_mutex);
				if (!_in->good ())
					break;
				for (int i = 0 ; i < _poolSize && _in->good (); ++i)
					_queryBuffer.emplace_back (*_in);
			}

			{	/// do searching
				for (const auto & _query : _queryBuffer)
				{
					if(_query.seq_size() < 12)
					{
						continue;
					}
					typedef std::integral_constant<bool, true> true_type;
					typedef std::integral_constant<bool, false> false_type;
                    if (_allowMismatch)
                        this->template start_tailing_match_Dual <true_type>(_query, &_resultBuffer, _minLen);
                    else
                        this->template start_tailing_match_Dual <false_type>(_query, &_resultBuffer, _minLen);
					this->start_end_pos_ = {0,0};
				}
			}

			{	/// writting
				boost::mutex::scoped_lock lock(this->_io_mutex);
				if(_resultBuffer.tellp() != 0)
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
