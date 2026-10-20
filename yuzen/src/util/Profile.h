/**************************************************************************
*
*  Profile.h : wall-clock time spent in named parts of a run
*
*  With the environment variable YUZEN_PROFILE set (to anything but 0), each
*  YUZEN_PROF("name") scope adds the wall-clock time it took, and one call,
*  to the total of that name; at exit the totals go to standard error, one
*  line each:
*
*      [profile] threads 8
*      [profile] h.linsolve        1060      12.345678
*
*  (name, calls, seconds). tools/profile_threads.py reads them. Without the
*  variable a scope costs one test of a static flag, and nothing is timed.
*  A scope inside a parallel region times nothing: the parts measured are
*  the serial stretches between parallel regions and the regions as a whole,
*  as the thread that runs the time loop sees them. Times of nested scopes
*  are each inclusive. Nothing here changes what the solver computes.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#ifndef YUZEN_PROFILE_H
#define YUZEN_PROFILE_H

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace yuzen_prof {

struct Entry
{
	const char *name;
	double seconds;
	long long calls;
};

inline std::vector<Entry>& Table()
{
	static std::vector<Entry> t;
	return t;
}

//  The totals, merged by name (a template instantiated for 1-D and 2-D
//  registers its scopes once for each), in the order first registered.
inline void Report()
{
	std::vector<Entry> merged;
	for(const Entry &e : Table())
	{
		bool found = false;
		for(Entry &m : merged)
			if(std::strcmp(m.name, e.name) == 0) { m.seconds += e.seconds; m.calls += e.calls; found = true; break; }
		if(!found) merged.push_back(e);
	}
#ifdef _OPENMP
	std::fprintf(stderr, "[profile] threads %d\n", omp_get_max_threads());
#else
	std::fprintf(stderr, "[profile] threads 1\n");
#endif
	for(const Entry &m : merged)
		std::fprintf(stderr, "[profile] %-20s %12lld %14.6f\n", m.name, m.calls, m.seconds);
}

inline bool Enabled()
{
	static const bool on = []() {
		const char *e = std::getenv("YUZEN_PROFILE");
		const bool b = (e != 0 && *e != 0 && std::strcmp(e, "0") != 0);
		if(b) std::atexit(Report);
		return b;
	}();
	return on;
}

inline int Register(const char *name)
{
	Table().push_back(Entry{name, 0.0, 0});
	return (int)Table().size() - 1;
}

inline double Now()
{
	return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

class Scope
{
public:
	explicit Scope(int id) : _id(id), _on(false), _t0(0.0)
	{
		if(!Enabled()) return;
#ifdef _OPENMP
		if(omp_in_parallel()) return;
#endif
		_on = true;
		_t0 = Now();
	}
	~Scope()
	{
		if(!_on) return;
		Entry &e = Table()[(size_t)_id];
		e.seconds += Now() - _t0;
		++e.calls;
	}
private:
	int _id;
	bool _on;
	double _t0;
	Scope(const Scope&);
	void operator=(const Scope&);
};

} // namespace yuzen_prof

#define YUZEN_PROF_CAT2(a, b) a##b
#define YUZEN_PROF_CAT(a, b) YUZEN_PROF_CAT2(a, b)
//  Times the rest of the enclosing block under the given name.
#define YUZEN_PROF(name) \
	static const int YUZEN_PROF_CAT(yz_prof_id_, __LINE__) = ::yuzen_prof::Register(name); \
	::yuzen_prof::Scope YUZEN_PROF_CAT(yz_prof_scope_, __LINE__)(YUZEN_PROF_CAT(yz_prof_id_, __LINE__))

#endif
