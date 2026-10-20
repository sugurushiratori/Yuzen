/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Opt.h : switches for the algorithms of the parallel-efficiency work
*
*  Every one of these is on by default. The environment variable YUZEN_OPT
*  turns them off, one by one or all at once, to compare a run against the
*  algorithm it replaced (tools/opt_effects.py measures each on its own):
*
*      YUZEN_OPT=-ew,-extrap   these two off, the others on
*      YUZEN_OPT=none          all off: the algorithms they replaced
*      YUZEN_OPT=none,bilu     only bilu on
*
*  Tokens are comma or space separated and applied in order, from the
*  default set; "default" (or "all") restores it.
*
*    regions   light loops run serially below 64 iterations per thread;
*              fewer parallel regions in the assembly, and its serial
*              loops made parallel                          (bit-identical)
*    assembly  the surface tension and viscosity of each face computed
*              once per Newton pass and shared by A and its Jacobian; a
*              product A h that was thrown away is not formed (bit-identical)
*    numa      the large arrays placed by first touch from all threads, as
*              the parallel loops later read them           (bit-identical)
*    ew        Eisenstat-Walker: the linear tolerance of each Newton pass
*              follows the progress of the Newton iteration
*    extrap    the Newton iteration starts from the film extrapolated from
*              the last two steps instead of the last one
*    concpc    a light preconditioner for the concentration solves (SPAI0
*              with AMGCL, Jacobi with Eigen) instead of AMG / ILUT
*    bilu      Eigen: block-Jacobi ILUT over a number of blocks fixed by
*              the size, factorized and applied in parallel (the result
*              does not depend on the thread count)
*    amgcoarse AMGCL: coarsen down to 500 unknowns (instead of 3000)
*              before the sequential direct solve of the coarsest level,
*              on a system of more than 3000 unknowns (one of 3000 or
*              fewer is solved directly, as without the switch)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#ifndef YUZEN_OPT_H
#define YUZEN_OPT_H

#include <cstdlib>
#include <cstring>
#include <string>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace yuzen { namespace opt {

enum Flag
{
	Regions   = 1u << 0,
	Assembly  = 1u << 1,
	EW        = 1u << 2,
	ConcPC    = 1u << 3,
	BlockILU  = 1u << 4,
	AmgCoarse = 1u << 6,
	Numa      = 1u << 7,
	Extrap    = 1u << 8
};

const unsigned kDefault = Regions | Assembly | Numa | EW | Extrap | ConcPC | BlockILU | AmgCoarse;

struct Token { const char *name; unsigned flag; };

inline const Token* Tokens(int &n)
{
	static const Token t[] = {
		{"regions", Regions}, {"assembly", Assembly}, {"numa", Numa}, {"ew", EW},
		{"extrap", Extrap}, {"concpc", ConcPC}, {"bilu", BlockILU}, {"amgcoarse", AmgCoarse}
	};
	n = (int)(sizeof(t) / sizeof(t[0]));
	return t;
}

/**
 * @brief The flags a YUZEN_OPT value leaves on, starting from kDefault
 *        (a null pointer: kDefault). `bad` collects the tokens that are
 *        not recognized (space separated), for the caller to report.
 */
inline unsigned Parse(const char *s, std::string *bad = 0)
{
	unsigned f = kDefault;
	if(s == 0) return f;
	int n; const Token *t = Tokens(n);
	std::string cur;
	for(const char *p = s; ; ++p)
	{
		const char ch = *p;
		if(ch == ',' || ch == ' ' || ch == '\t' || ch == '\0')
		{
			if(!cur.empty())
			{
				bool off = false;
				std::string name = cur;
				if(name[0] == '-') { off = true; name = name.substr(1); }
				bool found = true;
				if(!off && name == "none") f = 0;
				else if(!off && (name == "all" || name == "default")) f = kDefault;
				else
				{
					found = false;
					for(int i = 0; i < n && !found; ++i)
						if(name == t[i].name) { f = off ? (f & ~t[i].flag) : (f | t[i].flag); found = true; }
				}
				if(!found && bad) { if(!bad->empty()) *bad += " "; *bad += cur; }
				cur.clear();
			}
			if(ch == '\0') break;
		}
		else cur += ch;
	}
	return f;
}

/** @brief The flags of this run (YUZEN_OPT, read once). */
inline unsigned Flags()
{
	static const unsigned f = Parse(std::getenv("YUZEN_OPT"));
	return f;
}

inline bool On(unsigned flag) { return (Flags() & flag) != 0; }

/**
 * @brief Whether a light loop of n iterations (a copy, a clear, a boundary
 *        row) is worth a parallel region: with regions on (the default),
 *        only from 64 iterations per thread; with it off, always.
 *        Its fork and join cost about 0.2 us per thread (1.7 us at 8
 *        threads, 24 us at 112 on a Xeon Max 9480), a light iteration a few
 *        ns. Used as the if() clause of such loops; a loop run serially
 *        computes every element as the parallel one does.
 */
inline bool Par(long long n)
{
	if(!On(Regions)) return true;
#if defined(_OPENMP)
	return n >= 64LL * omp_get_max_threads();
#else
	return false;
#endif
}

/**
 * @brief For a loop that ran serially before regions: parallel only with
 *        regions on, and then as Par().
 */
inline bool ParNew(long long n) { return On(Regions) && Par(n); }

/** @brief The tokens of `f`, comma separated ("none" for 0). */
inline std::string Names(unsigned f)
{
	int n; const Token *t = Tokens(n);
	std::string r;
	for(int i = 0; i < n; ++i)
		if(f & t[i].flag) { if(!r.empty()) r += ","; r += t[i].name; }
	return r.empty() ? "none" : r;
}

} } // namespace yuzen::opt

#endif // YUZEN_OPT_H
