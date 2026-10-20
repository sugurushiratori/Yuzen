/**************************************************************************
*
*  ThreadCap.h : how many threads a run of a given size may use
*
*  Every parallel loop costs a fork and a join, whatever the work inside it;
*  a time step makes thousands of loops, and on a small grid each one has
*  too little work to pay for many threads. And past one socket the
*  reductions and the memory traffic cost more than the extra cores return.
*  On a two-socket Xeon Max 9480 node (56 cores a socket), leveling_1d_reg
*  (300 cells) ran 3.8 times slower on 56 threads than on 2 with AMGCL, and
*  every size measured, up to 1024^2, was slower on 112 threads than on 56
*  (by up to 45 %). So the thread count of a run is at most
*
*      T = min(cores of one socket, max(kMinThreads, ceil(cells / cellsPerThread)))
*
*  With the LLVM or Intel OpenMP runtime (libomp, libiomp5) the size term is
*  left out, unless YUZEN_THREAD_CAP sets it: linear solves then take threads
*  by their size (util/Grain.h), and on that node no case was more than 7 %
*  slower on 56 threads than at its best count (SizeTermApplies).
*
*  with cellsPerThread = kCellsPerThread by default. "cells" is the number
*  of grid cells the run solves on: Nx for a 1-D film, Nx x Ny for a 2-D
*  one, times Nz for ConcZ. On that node T was within 10 % of the fastest
*  thread count for every size measured (manual, sec:parallel:cap). The
*  socket term needs the processor topology (Linux: /sys/devices/system/
*  cpu); where it cannot be read, or where a socket has fewer than
*  kMinThreads cores, only the size term applies. A run asking for
*  kMinThreads threads or fewer is never affected.
*
*  What happens to a request above T is the caller's (Yuzen.cpp): asked for
*  explicitly (OMP_NUM_THREADS set), the run stops and says what to ask
*  for; left to the OpenMP default, the run uses T.
*
*  YUZEN_THREAD_CAP in the environment overrides it: 0 (or "off") turns the
*  cap off, a positive integer replaces cellsPerThread. Nothing else is read.
*
*  With the Eigen backend the answer does not depend on the thread count, so
*  the cap changes only the run time; with AMGCL it can move the last digits,
*  as any change of thread count can.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#ifndef YUZEN_THREADCAP_H
#define YUZEN_THREADCAP_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <map>
#include <set>
#include <utility>

namespace yuzen { namespace util {

const int kMinThreads     = 8;		// never capped below this
const int kCellsPerThread = 1000;	// cells each thread above kMinThreads must have

/**
 * @brief The thread limit of a run of `cells` grid cells.
 * @param cellsPerThread cells per thread; 0 or less means no cap
 * @param socketCores    physical cores of one socket; 0 or less: unknown
 * @return the limit T of the formula at the top of this file, or 0 when
 *         the cap is off
 */
inline int ThreadLimit(long long cells, int cellsPerThread = kCellsPerThread, int socketCores = 0,
                       bool sizeTerm = true)
{
	if(cellsPerThread <= 0) return 0;
	if(!sizeTerm) return socketCores >= kMinThreads ? socketCores : 0;
	const long long want = (cells + cellsPerThread - 1) / cellsPerThread;
	long long cap = want < kMinThreads ? kMinThreads : want;
	if(socketCores >= kMinThreads && cap > socketCores) cap = socketCores;
	return cap > 1000000 ? 1000000 : (int)cap;
}

/**
 * @brief The thread count a run of `cells` grid cells uses when `requested`
 *        threads are available: `requested` when it is kMinThreads or less
 *        or the cap is off, else the smaller of it and ThreadLimit().
 */
inline int CappedThreads(long long cells, int requested, int cellsPerThread = kCellsPerThread,
                         int socketCores = 0, bool sizeTerm = true)
{
	if(requested <= kMinThreads || cellsPerThread <= 0) return requested;
	const int cap = ThreadLimit(cells, cellsPerThread, socketCores, sizeTerm);
	return (cap > 0 && cap < requested) ? cap : requested;
}

/**
 * @brief Whether the size term applies: YUZEN_THREAD_CAP set (to a number of
 *        cells per thread), or an OpenMP runtime other than LLVM's and
 *        Intel's. With those two, linear solves take threads by their size
 *        (util/Grain.h) and the time of a small run no longer grows with the
 *        thread count, so only the socket term is left.
 * @param env       the value of YUZEN_THREAD_CAP, or a null pointer
 * @param kmpRuntime whether the runtime is LLVM's or Intel's
 */
inline bool SizeTermApplies(const char *env, bool kmpRuntime)
{
	return (env != 0 && *env != '\0') || !kmpRuntime;
}

/**
 * @brief cellsPerThread as YUZEN_THREAD_CAP sets it: kCellsPerThread when the
 *        variable is unset or empty, 0 (no cap) for "0", "off" or anything
 *        that is not a positive integer, else its value.
 * @param env the value of the variable, or a null pointer when it is unset
 */
inline int CellsPerThreadFromEnv(const char *env)
{
	if(env == 0 || *env == '\0') return kCellsPerThread;
	char *end = 0;
	const long v = std::strtol(env, &end, 10);
	if(end == env || *end != '\0' || v <= 0 || v > 1000000000L) return 0;
	return (int)v;
}

/**
 * @brief The physical cores of the smallest socket, from the processor
 *        topology Linux gives in /sys/devices/system/cpu (core_id and
 *        physical_package_id of each logical CPU; hyper-threads of one
 *        core counted once). 0 where it cannot be read. Read once.
 */
inline int CoresPerSocket()
{
	static int cached = -1;
	if(cached >= 0) return cached;
	std::map<int, std::set<int> > cores;		// package -> core ids
	for(int cpu = 0; cpu < 65536; ++cpu)
	{
		char dir[96], path[160];
		std::snprintf(dir, sizeof(dir), "/sys/devices/system/cpu/cpu%d", cpu);
		std::error_code ec;
		if(!std::filesystem::exists(dir, ec)) break;		// the CPUs are numbered from 0 without gaps
		int pkg = -1, core = -1;
		std::snprintf(path, sizeof(path), "%s/topology/physical_package_id", dir);
		if(FILE *f = std::fopen(path, "r")) { if(std::fscanf(f, "%d", &pkg) != 1) pkg = -1; std::fclose(f); }
		std::snprintf(path, sizeof(path), "%s/topology/core_id", dir);
		if(FILE *f = std::fopen(path, "r")) { if(std::fscanf(f, "%d", &core) != 1) core = -1; std::fclose(f); }
		if(pkg >= 0 && core >= 0) cores[pkg].insert(core);		// an offline CPU has no topology
	}
	int n = 0;
	for(const auto &p : cores)
		if(n == 0 || (int)p.second.size() < n) n = (int)p.second.size();
	cached = n;
	return n;
}

} } // namespace yuzen::util

#endif // YUZEN_THREADCAP_H
