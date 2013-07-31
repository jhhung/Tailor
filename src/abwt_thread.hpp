#ifndef ABWT_THREAD_HPP_
#define ABWT_THREAD_HPP_

#include <vector>
#include "abwt_format.hpp"
#include "abwt_search.hpp"
#include "boost/thread.hpp"

template <typename T>
class ABWT_thread : private ABWT_search<T>
{
private:
	Fastq _query;
	std::ostream* _out {nullptr};
	int _minLen {18};
public:
	ABWT_thread (T& table, Fastq&& q, std::ostream* out, int minLen) :
		ABWT_search<T> {table},
		_query {std::move (q)},
		_out {out},
		_minLen {minLen}
		{ }

	ABWT_thread (const ABWT_thread& other):
		ABWT_search<T> {other},
		_query {other._query},
		_out {other._out},
		_minLen {other._minLen}
	{}

	ABWT_thread (ABWT_thread&& other):
		ABWT_search<T> {other},
		_query {std::move (other._query)},
		_out {other._out},
		_minLen {other._minLen}
	{}

	void operator () () {
		this->start_tailing_match_Dual(_query, _out);
	}
};

template <typename T>
class ABWT_threads : private ABWT_search<T>
{
private:
	std::vector <Fastq> _queryPool {};
	std::ostream* _out {nullptr};
	int _minLen {18};
public:
	ABWT_threads (T& table, std::vector<Fastq>&& pool, std::ostream* out, int minLen) :
		ABWT_search<T> {table},
		_queryPool {std::move (pool)},
		_out {out},
		_minLen {minLen}
	{}
	ABWT_threads (const ABWT_threads& other):
			ABWT_search<T> {other},
			_queryPool {other._queryPool},
			_out {other._out},
			_minLen {other._minLen}
	{}

	ABWT_threads (ABWT_threads&& other):
			ABWT_search<T> {other},
			_queryPool {std::move (other._queryPool)},
			_out {other._out},
			_minLen {other._minLen}
	{}

	void operator () () {
		for (const auto & _query : _queryPool)
		this->start_tailing_match_Dual(_query, _out, _minLen);
		this->start_end_pos_ = {0,0};
	}
};


#endif /* ABWT_THREAD_HPP_ */
