/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_config.cpp : test_config.cpp  --  regression tests for IO/ConfigBlock  *
*  Pins down how the .inp format is interpreted.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "doctest/doctest.h"
#include "ConfigBlock.h"

#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace yuzen::config;

namespace {

// Helper that creates a temporary file, has it parsed, and cleans up
// afterwards
struct TempFile
{
	std::string path;
	TempFile(const std::string &name, const std::string &body) : path(name)
	{
		std::ofstream f(path.c_str());
		f << body;
	}
	~TempFile() { std::remove(path.c_str()); }
};

} // namespace

TEST_CASE("ConfigBlock: scalar variables in .inp")
{
	TempFile f("test_cfg_scalar.inp",
		"// a comment line\n"
		"SolverName = \"Yuzen-1.0\"\n"
		"Count  = 42\n"
		"Ratio  = 0.25\n"
		"Big    = 1.5e-3\n"
		"Flag   = true\n"
		"Off    = false\n"
		"# the other comment marker\n");

	Block *b = Block::newread(f.path);
	REQUIRE(b != 0);

	std::string s; b->getin("SolverName", s); CHECK(s == "Yuzen-1.0");
	int    i = 0;   b->getin("Count",  i);    CHECK(i == 42);
	double d = 0;   b->getin("Ratio",  d);    CHECK(d == doctest::Approx(0.25));
	double e = 0;   b->getin("Big",    e);    CHECK(e == doctest::Approx(1.5e-3));
	bool   t = false; b->getin("Flag", t);    CHECK(t == true);
	bool   u = true;  b->getin("Off",  u);    CHECK(u == false);

	// With a default: a variable that is absent keeps the default
	int missing = 7; b->getin("NoSuchKey", missing, 7); CHECK(missing == 7);
	// A variable that is absent and has no default throws
	CHECK_THROWS_AS(b->getin("NoSuchKey", i), VariableNotFound);

	delete b;
}

//  An integer key given a value outside int, or no number at all, must come
//  out the same on every platform: (int) of such a double is undefined, and
//  x86 and ARM give different answers for NaN. Out of range saturates; NaN
//  and inf are not whole numbers and are refused.
TEST_CASE("ConfigBlock: integers out of range saturate, NaN and inf are refused")
{
	TempFile f("test_cfg_intsat.toml",
		"Big   = 4000000000\n"
		"Small = -4.0e9\n"
		"NaN   = nan\n"
		"Inf   = inf\n");
	Block *b = Block::newread(f.path);
	REQUIRE(b != 0);
	int i = 0;
	b->getin("Big",   i); CHECK(i == std::numeric_limits<int>::max());
	b->getin("Small", i); CHECK(i == std::numeric_limits<int>::min());
	CHECK_THROWS(b->getin("NaN", i));
	CHECK_THROWS(b->getin("Inf", i));
	delete b;
}

TEST_CASE("ConfigBlock: vector variables in .inp")
{
	TempFile f("test_cfg_vec.inp",
		"Center = (0.5, 0)\n"
		"Width  = (1.0, 0.5)\n"
		"Names  = (\"a\", \"b\", \"c\")\n");

	Block *b = Block::newread(f.path);
	std::vector<double> c; b->getin("Center", c);
	REQUIRE(c.size() == 2);
	CHECK(c[0] == doctest::Approx(0.5));
	CHECK(c[1] == doctest::Approx(0.0));

	std::vector<std::string> n; b->getin("Names", n);
	REQUIRE(n.size() == 3);
	CHECK(n[1] == "b");

	std::vector<double> def; def.push_back(9.0);
	std::vector<double> miss; b->getin("NoSuchVec", miss, def);
	REQUIRE(miss.size() == 1);
	CHECK(miss[0] == doctest::Approx(9.0));

	delete b;
}

TEST_CASE("ConfigBlock: nested and repeated blocks in .inp")
{
	TempFile f("test_cfg_block.inp",
		"block Grid\n"
		"  Nx = 300\n"
		"  block Sub\n"
		"    Depth = 2\n"
		"  end\n"
		"end\n"
		"block Substrate\n"
		"  block RectBump\n    Order = 1\n  end\n"
		"  block RectBump\n    Order = 2\n  end\n"
		"  block RectBump\n    Order = 3\n  end\n"
		"end\n");

	Block *b = Block::newread(f.path);

	CHECK(b->getBlockSize("Grid") == 1);
	int nx = 0; b->getBlock("Grid")->getin("Nx", nx); CHECK(nx == 300);

	int depth = 0; b->getBlock("Grid")->getBlock("Sub")->getin("Depth", depth);
	CHECK(depth == 2);

	CHECK(b->getBlock("Substrate")->getBlockSize("RectBump") == 3);
	for(int k = 0; k < 3; ++k)
	{
		int order = 0;
		b->getBlock("Substrate")->getBlock("RectBump", k)->getin("Order", order);
		CAPTURE(k);
		CHECK(order == k + 1);
	}

	// A block that is absent throws (callers use this to test for its
	// presence)
	CHECK_THROWS_AS(b->getBlockSize("NoSuchBlock"), BlockNotFound);
	CHECK_THROWS_AS(b->getBlock("NoSuchBlock"), BlockNotFound);

	delete b;
}

TEST_CASE("ConfigBlock: a .toml input gives the same tree as the block-structured .inp")
{
	TempFile a("test_cfg_same.inp",
		"Name = \"x\"\n"
		"block Grid\n  Nx = 300\n  Vec = (1.0, 2.0)\n end\n"
		"block Rep\n  V = 1\n end\n"
		"block Rep\n  V = 2\n end\n");
	TempFile t("test_cfg_same.toml",
		"Name = \"x\"\n"
		"[Grid]\nNx = 300\nVec = [1.0, 2.0]\n"
		"[[Rep]]\nV = 1\n"
		"[[Rep]]\nV = 2\n");

	Block *bi = Block::newread(a.path);
	Block *bt = Block::newread(t.path);

	std::string s1, s2;
	bi->getin("Name", s1); bt->getin("Name", s2);
	CHECK(s1 == s2);

	int n1 = 0, n2 = 0;
	bi->getBlock("Grid")->getin("Nx", n1);
	bt->getBlock("Grid")->getin("Nx", n2);
	CHECK(n1 == n2);

	std::vector<double> v1, v2;
	bi->getBlock("Grid")->getin("Vec", v1);
	bt->getBlock("Grid")->getin("Vec", v2);
	REQUIRE(v1.size() == v2.size());
	for(size_t k = 0; k < v1.size(); ++k) CHECK(v1[k] == doctest::Approx(v2[k]));

	CHECK(bi->getBlockSize("Rep") == bt->getBlockSize("Rep"));
	for(int k = 0; k < 2; ++k)
	{
		int x1 = 0, x2 = 0;
		bi->getBlock("Rep", k)->getin("V", x1);
		bt->getBlock("Rep", k)->getin("V", x2);
		CAPTURE(k);
		CHECK(x1 == x2);
	}

	delete bi; delete bt;
}

TEST_CASE("ConfigBlock: a file that does not exist throws BlockNotFound")
{
	CHECK_THROWS_AS(Block::newread("no_such_file_12345.inp"), BlockNotFound);
}

TEST_CASE("ConfigBlock: entries nobody asked for are reported, with the likely intent")
{
	//  collectUnknown() is what Input::Read() stops on. Before it existed an
	//  entry no reader asked for was simply never read: a misspelled key kept
	//  its default, and a key in the wrong table did too, without a word.
	TempFile t("test_cfg_unknown.toml",
		"Name = \"x\"\n"
		"Stray = 1\n"
		"[Initial]\n"
		"PrecursorRatio = 0.02\n"
		"[Initial.SurfaceShape]\n"
		"TotalHeight = 1e-5\n"
		"[Initial.SurfaceShape.EdgeBump]\n"
		"alpha = 2.0\n"
		"InnerLength = 1e-4\n"
		"[Physics]\n"
		"[Physics.Evaporaton]\n"
		"Coef = 1e-7\n");
	Block *b = Block::newread(t.path);
	REQUIRE(b != 0);

	//  What a reader would ask for.
	std::string s; b->getin("Name", s);
	Block *ini = b->getBlock("Initial");
	Block *shp = ini->getBlock("SurfaceShape");
	double d = 0; shp->getin("TotalHeight", d); shp->getin("PrecursorRatio", d, 0.0);
	Block *eb = shp->getBlock("EdgeBump");
	eb->getin("Alpha", d, 1.0); eb->getin("InnerLength", d, 1.0);
	Block *ph = b->getBlock("Physics");
	CHECK_THROWS_AS(ph->getBlockSize("Evaporation"), BlockNotFound);

	std::vector<std::string> u;
	b->collectUnknown("", u);
	std::string all;
	for(size_t i = 0; i < u.size(); ++i) all += u[i] + "\n";
	INFO(all);
	REQUIRE(u.size() == 4);
	CHECK(all.find("key   Stray") != std::string::npos);
	CHECK(all.find("key   Initial.PrecursorRatio   (this key is read in [Initial.SurfaceShape])") != std::string::npos);
	CHECK(all.find("key   Initial.SurfaceShape.EdgeBump.alpha   (did you mean Initial.SurfaceShape.EdgeBump.Alpha?)") != std::string::npos);
	CHECK(all.find("block [Physics.Evaporaton]   (did you mean [Physics.Evaporation]?)") != std::string::npos);

	//  An input whose every entry is asked for reports nothing.
	TempFile t2("test_cfg_known.toml", "Name = \"x\"\n[Grid]\nNx = 3\n");
	Block *c = Block::newread(t2.path);
	c->getin("Name", s); int nx = 0; c->getBlock("Grid")->getin("Nx", nx);
	std::vector<std::string> v; c->collectUnknown("", v);
	CHECK(v.empty());

	delete b; delete c;
}

//------------------------------------------------------------------
//  A value is converted whole or refused, with the key named (the fuzzer
//  of the input reader, yuzen/tests/fuzz/fuzz_input.cpp, feeds such texts).
//  What strtod makes of a text is not enough: "fast" would be 0, "1e-3x"
//  1e-3, a NaN would pass the range checks and be cast to int, "2.7" for an
//  int would be 3, and any word that is not a spelling of true false.
//------------------------------------------------------------------
TEST_CASE("ConfigBlock: a value that is not of the type asked for is refused")
{
	TempFile f("test_cfg_types.toml",
		"Word = \"fast\"\n"
		"Trail = \"1e-3x\"\n"
		"NotANumber = nan\n"
		"Infinite = -inf\n"
		"Fraction = 2.7\n"
		"Whole = 3.0\n"
		"Exp = 1e3\n"
		"Big = 4000000000\n"
		"Maybe = \"maybe\"\n"
		"No = \"no\"\n"
		"Zero = 0\n"
		"Empty = []\n"
		"Two = [1, 2]\n"
		"One = [5]\n"
		"Pair = [0.5, \"x\"]\n"
		"[Grid]\n"
		"Nx = \"many\"\n");
	Block *b = Block::newread(f.path);
	REQUIRE(b != 0);
	double d = 0; int i = 0; bool t = false; std::string s;
	std::vector<double> v;

	CHECK_THROWS_AS(b->getin("Word", d), VariableNotFound);
	CHECK_THROWS_AS(b->getin("Word", d, 1.0), VariableNotFound);	// present, so not the default
	CHECK_THROWS_AS(b->getin("Trail", d), VariableNotFound);
	CHECK_THROWS_AS(b->getin("NotANumber", d), VariableNotFound);
	CHECK_THROWS_AS(b->getin("NotANumber", i), VariableNotFound);
	CHECK_THROWS_AS(b->getin("Infinite", d), VariableNotFound);
	CHECK_THROWS_AS(b->getin("Fraction", i), VariableNotFound);
	b->getin("Whole", i);  CHECK(i == 3);
	b->getin("Exp", i);    CHECK(i == 1000);
	b->getin("Big", i);    CHECK(i == 2147483647);		// saturates
	CHECK_THROWS_AS(b->getin("Maybe", t), VariableNotFound);
	t = true;  b->getin("No", t);   CHECK(t == false);
	t = true;  b->getin("Zero", t); CHECK(t == false);
	CHECK_THROWS_AS(b->getin("Empty", d, 1.0), VariableNotFound);	// not the default
	CHECK_THROWS_AS(b->getin("Empty", s), VariableNotFound);
	CHECK_THROWS_AS(b->getin("Two", d), VariableNotFound);			// not the first
	b->getin("One", d);    CHECK(d == 5.0);
	CHECK_THROWS_AS(b->getin("Pair", v), VariableNotFound);
	b->getin("Two", v);    CHECK(v.size() == 2);

	//  The message names the key, with its block
	try { b->getBlock("Grid")->getin("Nx", i); FAIL("not refused"); }
	catch(const VariableNotFound &e) { CHECK(std::string(e.what()).find("Grid.Nx") != std::string::npos); }
	delete b;
}

TEST_CASE("ConfigBlock: TOML values the input has no use for are refused, not dropped")
{
	//  A date and an array inside an array were dropped without a word: the
	//  key was neither read nor reported as unknown.
	TempFile f1("test_cfg_date.toml", "[Numerics]\nEndTime = 1979-05-27\n");
	CHECK_THROWS_AS(Block::newread(f1.path), BlockNotFound);
	TempFile f2("test_cfg_nested.toml", "Center = [0.5, [1.0], 0.0]\n");
	CHECK_THROWS_AS(Block::newread(f2.path), BlockNotFound);
	TempFile f3("test_cfg_tablearr.toml", "Center = [0.5, {a = 1}]\n");
	CHECK_THROWS_AS(Block::newread(f3.path), BlockNotFound);
}

TEST_CASE("ConfigBlock: blocks nested without end are refused, not recursed into")
{
	//  The parsers recurse once per level: without a limit, 200000 nested
	//  blocks would run the .inp parser out of stack.
	std::string inp, toml = "[a";
	for(int k = 0; k < 200000; ++k) inp += "block a\n";
	for(int k = 0; k < 100; ++k) toml += ".a";
	TempFile f1("test_cfg_deep.inp", inp);
	CHECK_THROWS_AS(Block::newread(f1.path), BlockNotFound);
	TempFile f2("test_cfg_deep.toml", toml + "]\nx = 1\n");
	CHECK_THROWS_AS(Block::newread(f2.path), BlockNotFound);
}

TEST_CASE("ConfigBlock: a non-ASCII character in a TOML key is a parse error")
{
	//  toml++ 3.4.0 reached TOML_UNREACHABLE (undefined behavior) for the
	//  characters U+00A1..U+0499, U+2C5E..U+3057 (but U+3000) and
	//  U+FB26..U+FEFE in a bare key or a table header (patched in
	//  third_party/toml++/toml.hpp). U+3042 is one.
	TempFile f1("test_cfg_unicode1.toml", "[Physics.Soluto\xc4\x9e" "pillarity]\nEnable = true\n");
	CHECK_THROWS_AS(Block::newread(f1.path), BlockNotFound);
	TempFile f2("test_cfg_unicode2.toml", "\xe3\x81\x82 = 1\n");
	CHECK_THROWS_AS(Block::newread(f2.path), BlockNotFound);
	//  ... and a comment or a string may hold any of them
	TempFile f3("test_cfg_unicode3.toml", "# \xe3\x81\x82\xc3\xa9\nName = \"\xe3\x81\x82\"\n");
	Block *b = Block::newread(f3.path);
	std::string s; b->getin("Name", s); CHECK(s == "\xe3\x81\x82");
	delete b;
}
