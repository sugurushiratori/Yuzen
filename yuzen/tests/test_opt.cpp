/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_opt.cpp : the parsing of YUZEN_OPT (util/Opt.h) and YUZEN_GRAIN
*                 (util/Grain.h)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "doctest/doctest.h"
#include "Opt.h"
#include "Grain.h"

using namespace yuzen::opt;

TEST_CASE("Opt: the default set, tokens on and off, none, all and unknown ones")
{
	CHECK(Parse(0) == kDefault);
	CHECK(Parse("") == kDefault);
	CHECK(kDefault == (unsigned)(Regions | Assembly | Numa | EW | Extrap | ConcPC | BlockILU | AmgCoarse));
	CHECK(Parse("none") == 0u);
	CHECK(Parse("none,regions") == (unsigned)Regions);
	CHECK(Parse("-ew,-extrap") == (kDefault & ~(unsigned)(EW | Extrap)));
	CHECK(Parse(" none  bilu numa ") == (unsigned)(Numa | BlockILU));
	CHECK(Parse("none,default") == kDefault);
	CHECK(Parse("none,all") == kDefault);
	std::string bad;
	CHECK(Parse("none,amgcoarse,foo,-bar", &bad) == (unsigned)AmgCoarse);
	CHECK(bad == "foo -bar");
	CHECK(Names(0) == "none");
	CHECK(Names(Regions | Extrap) == "regions,extrap");
	CHECK(Parse(("none," + Names(kDefault)).c_str()) == kDefault);
}

TEST_CASE("Grain: parsing YUZEN_GRAIN and the thread count of a linear solve")
{
	namespace G = yuzen::grain;
	CHECK(G::Parse(0, 500) == 500);			// unset: the default of the runtime
	CHECK(G::Parse("", 0) == 0);
	CHECK(G::Parse("lin=1000", 500) == 1000);
	CHECK(G::Parse("off", 500) == 0);
	CHECK(G::Parse("lin=0", 500) == 0);
	std::string bad;
	CHECK(G::Parse("foo=3,lin=x,lin=-2,amg=5,lin=250", 500, &bad) == 250);
	CHECK(bad == "foo=3 lin=x lin=-2 amg=5");
	CHECK(G::Names(0) == "off");
	CHECK(G::Names(500) == "lin=500");
	CHECK(G::Threads(3600, 0, 56) == 56);		// no grain: the run's count
	CHECK(G::Threads(3600, 1000, 56) == 4);
	CHECK(G::Threads(3601, 1200, 56) == 4);	// ceil
	CHECK(G::Threads(10, 500, 56) == 1);
	CHECK(G::Threads(1000000, 500, 56) == 56);
	CHECK(G::Threads(1000000, 500, 1) == 1);
}
