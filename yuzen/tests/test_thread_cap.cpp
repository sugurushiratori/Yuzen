/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_thread_cap.cpp : the thread count a run of a given size uses
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "doctest/doctest.h"
#include "ThreadCap.h"

using yuzen::util::CappedThreads;
using yuzen::util::CellsPerThreadFromEnv;
using yuzen::util::kCellsPerThread;
using yuzen::util::kMinThreads;
using yuzen::util::ThreadLimit;
using yuzen::util::CoresPerSocket;

TEST_CASE("ThreadCap: a request of kMinThreads or fewer is never changed")
{
	for(int t = 1; t <= kMinThreads; ++t)
	{
		CHECK(CappedThreads(1, t) == t);
		CHECK(CappedThreads(300, t) == t);
		CHECK(CappedThreads(100000000LL, t) == t);
	}
}

TEST_CASE("ThreadCap: min(socket, max(kMinThreads, ceil(cells/cellsPerThread))), never above the request")
{
	CHECK(kMinThreads == 8);
	CHECK(kCellsPerThread == 1000);
	CHECK(CappedThreads(300, 112) == 8);			// leveling_1d_reg
	CHECK(CappedThreads(3600, 112) == 8);			// leveling_2d_reg
	CHECK(CappedThreads(8000, 112) == 8);
	CHECK(CappedThreads(8001, 112) == 9);
	CHECK(CappedThreads(16384, 112) == 17);			// 128^2
	CHECK(CappedThreads(65536, 112) == 66);			// 256^2, the socket unknown
	CHECK(CappedThreads(65536, 112, kCellsPerThread, 56) == 56);
	CHECK(CappedThreads(108000, 112, kCellsPerThread, 56) == 56);	// concz_2d
	CHECK(CappedThreads(262144, 112, kCellsPerThread, 56) == 56);	// 512^2
	CHECK(CappedThreads(262144, 112) == 112);		// the request
	CHECK(CappedThreads(262144, 48, kCellsPerThread, 56) == 48);
	CHECK(CappedThreads(3000000000LL, 112) == 112);	// no overflow
	CHECK(CappedThreads(3600, 112, kCellsPerThread, 4) == 8);		// a socket below kMinThreads is not used
	CHECK(CappedThreads(65536, 112, 1500) == 44);
	CHECK(CappedThreads(300, 112, 0) == 112);		// cap off
	CHECK(CappedThreads(300, 112, -1) == 112);
	CHECK(ThreadLimit(300, 0) == 0);
	CHECK(ThreadLimit(300) == 8);
	CHECK(ThreadLimit(16384, kCellsPerThread, 56) == 17);
}

TEST_CASE("ThreadCap: the cores of a socket are read once, never negative")
{
	const int a = CoresPerSocket();
	CHECK(a >= 0);
	CHECK(CoresPerSocket() == a);
}

TEST_CASE("ThreadCap: YUZEN_THREAD_CAP")
{
	CHECK(CellsPerThreadFromEnv(0) == kCellsPerThread);
	CHECK(CellsPerThreadFromEnv("") == kCellsPerThread);
	CHECK(CellsPerThreadFromEnv("0") == 0);
	CHECK(CellsPerThreadFromEnv("off") == 0);
	CHECK(CellsPerThreadFromEnv("-5") == 0);
	CHECK(CellsPerThreadFromEnv("1500x") == 0);
	CHECK(CellsPerThreadFromEnv("1000") == 1000);
	CHECK(CellsPerThreadFromEnv("1") == 1);
}

TEST_CASE("ThreadCap: socket only with the LLVM and Intel runtimes")
{
	using yuzen::util::SizeTermApplies;
	CHECK(SizeTermApplies(0, false));			// libgomp: the size term
	CHECK(!SizeTermApplies(0, true));			// libomp: the socket only
	CHECK(!SizeTermApplies("", true));
	CHECK(SizeTermApplies("1000", true));		// asked for
	CHECK(CappedThreads(300, 112, kCellsPerThread, 56, false) == 56);
	CHECK(CappedThreads(300, 56, kCellsPerThread, 56, false) == 56);
	CHECK(CappedThreads(300, 112, kCellsPerThread, 0, false) == 112);	// the socket unknown: no cap
	CHECK(CappedThreads(300, 112, 0, 56, false) == 112);			// cap off
	CHECK(ThreadLimit(300, kCellsPerThread, 56, false) == 56);
	CHECK(ThreadLimit(300, kCellsPerThread, 4, false) == 0);
}
