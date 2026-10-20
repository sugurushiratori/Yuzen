/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Grain.h : the thread count of a linear solve, by its size
*
*  A parallel region costs a fork and a join whatever the work inside it, and
*  a linear solve of a few thousand unknowns makes a thousand regions a time
*  step (an AMG cycle, a BiCGStab iteration). On a run of T threads each of
*  them uses T, and below some hundreds of unknowns a thread that cost is
*  more than the work: the time of a small run then grows with T. So a
*  linear solve of n unknowns runs on
*
*      threads = min(T, max(1, ceil(n / grain)))
*
*  and the rest of the step (the assembly, the explicit terms) on T.
*
*  That pays only where changing the size of a team is cheap. The LLVM and
*  Intel OpenMP runtimes (libomp, libiomp5) keep the threads a smaller team
*  leaves idle at hand when KMP_HOT_TEAMS_MODE=1, which PrepareRuntime() sets
*  (unless it is set already) before the runtime starts. GCC's libgomp parks
*  them, and waking them again for the next large region costs more than the
*  smaller team saves; there the default is no grain. Measured on a Xeon Max
*  9480 socket (56 cores), 300 to 10^6 cells: with libomp and grain 500 the
*  time no longer grows with T anywhere (by at most 7 %), and the fastest
*  time is 1.15 times that of libgomp without it (tools/thread_variants.py).
*
*  YUZEN_GRAIN overrides the default: "lin=G" sets the grain, "off" (or
*  "lin=0") turns it off. With Eigen the result does not depend on the thread
*  count; with AMGCL it can, as with any change of thread count.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#ifndef YUZEN_GRAIN_H
#define YUZEN_GRAIN_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(_OPENMP) && (defined(__GNUC__) || defined(__clang__))
//  Exported by the LLVM and Intel runtimes only: a null address means another
//  runtime (libgomp)
extern "C" void kmp_set_defaults(char const *) __attribute__((weak));
#endif

namespace yuzen { namespace grain {

const long long kGrainKmp = 500;	// the default grain with libomp / libiomp5

/** @brief Whether the OpenMP runtime is LLVM's or Intel's (hot teams). */
inline bool KmpRuntime()
{
#if defined(_OPENMP) && (defined(__GNUC__) || defined(__clang__))
	return &kmp_set_defaults != 0;
#else
	return false;
#endif
}

/** @brief The OpenMP runtime, for the log and JobInfo. */
inline const char* RuntimeName()
{
#if defined(_OPENMP)
	return KmpRuntime() ? "LLVM/Intel (libomp, libiomp5)" : "other (libgomp)";
#else
	return "none";
#endif
}

/**
 * @brief Before the first OpenMP call (the start of main): hot teams with
 *        the LLVM and Intel runtimes, unless KMP_HOT_TEAMS_MODE is set.
 */
inline void PrepareRuntime()
{
	if(!KmpRuntime() || std::getenv("KMP_HOT_TEAMS_MODE")) return;
#if defined(_WIN32)
	_putenv_s("KMP_HOT_TEAMS_MODE", "1");
#else
	setenv("KMP_HOT_TEAMS_MODE", "1", 0);
#endif
}

/**
 * @brief The grain a YUZEN_GRAIN value sets: `deflt` for a null pointer or
 *        an empty value, 0 for "off". `bad` collects what is not recognized.
 */
inline long long Parse(const char *s, long long deflt, std::string *bad = 0)
{
	if(s == 0 || *s == '\0') return deflt;
	long long g = deflt;
	std::string cur;
	for(const char *p = s; ; ++p)
	{
		const char ch = *p;
		if(ch == ',' || ch == ' ' || ch == ';' || ch == '\0')
		{
			if(!cur.empty())
			{
				bool ok = false;
				if(cur == "off") { g = 0; ok = true; }
				else if(cur.compare(0, 4, "lin=") == 0)
				{
					const std::string val = cur.substr(4);
					char *end = 0;
					const long long v = std::strtoll(val.c_str(), &end, 10);
					if(!val.empty() && *end == '\0' && v >= 0) { g = v; ok = true; }
				}
				if(!ok && bad) { if(!bad->empty()) *bad += " "; *bad += cur; }
				cur.clear();
			}
			if(ch == '\0') break;
		}
		else cur += ch;
	}
	return g;
}

/** @brief The grain of this run (YUZEN_GRAIN and the runtime, read once). */
inline long long Get()
{
	static const long long g = Parse(std::getenv("YUZEN_GRAIN"), KmpRuntime() ? kGrainKmp : 0);
	return g;
}

/** @brief "lin=G" or "off". */
inline std::string Names(long long g)
{
	if(g <= 0) return "off";
	std::ostringstream o;
	o<<"lin="<<g;
	return o.str();
}

/** @brief min(T, max(1, ceil(n / g))); T when g is 0. */
inline int Threads(long long n, long long g, int T)
{
	if(g <= 0 || T <= 1) return T;
	long long k = (n + g - 1) / g;
	if(k < 1) k = 1;
	return k < T ? (int)k : T;
}

/** @brief The thread count of the run, as set before the model is built. */
inline int& RunThreads() { static int t = 0; return t; }

/** @brief Called once the thread count of the run is final. */
inline void Start()
{
#if defined(_OPENMP)
	RunThreads() = omp_get_max_threads();
#else
	RunThreads() = 1;
#endif
}

/**
 * @brief A linear solve of n unknowns: for its lifetime, the thread count of
 *        its size.
 */
class SolveScope
{
public:
	explicit SolveScope(long long n) : _saved(0), _set(false)
	{
#if defined(_OPENMP)
		const long long g = Get();
		if(g <= 0) return;
		const int T  = RunThreads() > 0 ? RunThreads() : omp_get_max_threads();
		const int nt = Threads(n, g, T);
		_saved = omp_get_max_threads();
		if(nt == _saved) return;
		omp_set_num_threads(nt);
		_set = true;
#else
		(void)n;
#endif
	}
	~SolveScope()
	{
#if defined(_OPENMP)
		if(_set) omp_set_num_threads(_saved);
#endif
	}
private:
	int _saved;
	bool _set;
	SolveScope(const SolveScope&);
	SolveScope& operator=(const SolveScope&);
};

} } // namespace yuzen::grain

#endif // YUZEN_GRAIN_H
