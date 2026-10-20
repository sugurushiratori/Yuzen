/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_restart.cpp : pins down the indexing of the ASCII restart
*
*  The Neutral format (field*****.dat) is both an output and the input of
*  an ASCII restart. The "line number that was written" and the "grid point
*  that is read" must therefore agree. Two ways for them not to:
*
*    - 1D: reading with Of<Xedge>(ix+1). Of() takes a **linear index**, and
*          in 1-D linear index = x index + 1, so the value would go in one
*          cell too early, and the shape move by one cell every time the
*          output is read back.
*    - 2D: dividing the line number by Ny instead of Nx when splitting it
*          into ix / jy. On a grid with Nx != Ny the read would be
*          completely broken.
*
*  Both are caught by checking that "writing and reading back the same
*  shape returns it to the same places".
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "doctest/doctest.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "LubFields.h"
#include "Restart.h"

using namespace yuzen;
using namespace yuzen::grid;

namespace {

Input::I_Grid MakeGrid(const char *dim, int nx, int ny)
{
	Input::I_Grid g;
	g.Dimension = dim;
	g.Xmin = -1.0; g.Xmax = 1.0;
	g.Ymin = -1.0; g.Ymax = 1.0;
	g.Nx = nx; g.Ny = ny; g.Nz = 1;
	return g;
}

Input::I_Substrate FlatSubstrate()
{
	Input::I_Substrate s;
	s.TotalHeight = 0.0;
	return s;
}

//  A temporary file. It is used from doctest, so a fixed name is fine.
struct TempFile
{
	std::string path;
	explicit TempFile(const char *name) : path(std::string("./") + name) {}
	~TempFile() { std::remove(path.c_str()); }
};

} // namespace


//==================================================================
//  1D
//==================================================================
TEST_CASE("ASCII restart (1D): line i lands at x index i+1")
{
	const int NX = 16;

	Field_WM1D field;
	field.Build(MakeGrid("xh", NX, 1), FlatSubstrate(), 1.0);
	C_Grid1D *G = field.GetGridPtr();

	//  The same number of lines as the output (= nx()-2). The two ends are
	//  boundaries, and in between there are NX cells.
	const int nrow = G->nx() - 2;
	REQUIRE(nrow == NX + 2);

	TempFile tf("test_restart_1d.dat");
	{
		std::ofstream fo(tf.path.c_str());
		fo << "xx \t HS \t Ux \t CC \t HH \t SS\n";
		for(int r = 0; r < nrow; r++)
			//  Put the line number itself into HH. Fill SS and CC with
			//  different values as well.
			fo << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t"
			   << (0.001 * r) << "\t" << (100.0 + r) << "\t" << (1000.0 + r) << "\n";
	}

	for(Xedge xe = G->First<Xedge>(); xe != G->Last<Xedge>(); ++xe)
	{
		field.HH(xe) = -1; field.SS(xe) = -1; field.CC(xe) = -1;
	}

	Input::I_Restart inp;
	inp.Enable = true; inp.ASCII = true; inp.FileName = tf.path;
	inp.RestartTime = 0; inp.FileIndex = 0;

	REQUIRE(C_Restart::ReadAscii_WM1D(inp, field, 1.0) == true);

	for(int r = 0; r < nrow; r++)
	{
		Xedge xe = G->At<Xedge>(r + 1);
		REQUIRE(field.HH(xe) == doctest::Approx(100.0 + r));
		REQUIRE(field.SS(xe) == doctest::Approx(1000.0 + r));
		REQUIRE(field.CC(xe) == doctest::Approx(0.001 * r));
	}
}


TEST_CASE("ASCII restart (1D): writing out and reading back returns the same places")
{
	const int NX = 16;

	Field_WM1D field;
	field.Build(MakeGrid("xh", NX, 1), FlatSubstrate(), 1.0);
	C_Grid1D *G = field.GetGridPtr();

	//  Put in a clearly asymmetric shape
	for(Xedge xe = G->First<Xedge>(); xe != G->Last<Xedge>(); ++xe)
	{
		const int ix = G->ix(xe);
		field.HH(xe) = 1.0 + 0.125 * ix;
		field.SS(xe) = 0.0;
		field.CC(xe) = 0.5;
	}

	//  Build the ordering used on the output side (Field_WM1D_Out packs it
	//  into lines)
	Field_WM1D_Out out;
	out.Build(field, 1.0, 1.0, 1.0, false);
	out.ConvertFromStg(field);

	const int nrow = out.nx();
	TempFile tf("test_restart_1d_rt.dat");
	{
		std::ofstream fo(tf.path.c_str());
		fo.precision(17);
		fo << "xx \t HS \t Ux \t CC \t HH \t SS\n";
		for(int r = 0; r < nrow; r++)
			fo << out.xx(r) << "\t" << out.HS(r) << "\t" << out.Ux(r) << "\t"
			   << out.CC(r) << "\t" << out.HH(r) << "\t" << out.SS(r) << "\n";
	}

	Field_WM1D back;
	back.Build(MakeGrid("xh", NX, 1), FlatSubstrate(), 1.0);

	Input::I_Restart inp;
	inp.Enable = true; inp.ASCII = true; inp.FileName = tf.path;
	inp.RestartTime = 0; inp.FileIndex = 0;
	REQUIRE(C_Restart::ReadAscii_WM1D(inp, back, 1.0) == true);

	//  The outermost ghost layer is not written out, so we compare the
	//  interior only
	for(int ix = 2; ix <= G->nx() - 3; ix++)
	{
		Xedge xe = G->At<Xedge>(ix);
		REQUIRE(back.HH(xe) == doctest::Approx(field.HH(xe)).epsilon(1e-13));
	}
}


//==================================================================
//  2D
//==================================================================
TEST_CASE("ASCII restart (2D): the line decomposition is right when Nx != Ny")
{
	//  This is the crucial point. Unless Nx and Ny are made different, a
	//  division by Ny instead of Nx stays hidden.
	const int NX = 12, NY = 5;

	Field_WM2D field;
	field.Build(MakeGrid("xyh", NX, NY), FlatSubstrate(), 1.0);
	C_Grid2D *G = field.GetGridPtr();

	for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf)
	{
		field.HH(zf) = 1.0 + G->ix(zf) + 0.01 * G->iy(zf);
		field.SS(zf) = 0.0;
		field.CC(zf) = 0.5;
	}

	Field_WM2D_Out out;
	out.Build(field, 1.0, 1.0, 1.0, false);
	out.ConvertFromStg(field);

	TempFile tf("test_restart_2d_rt.dat");
	{
		std::ofstream fo(tf.path.c_str());
		fo.precision(17);
		fo << "xx \t yy \t HS \t Ux \t Uy \t CC \t HH \t SS\n";
		for(int jy = 0; jy < out.ny(); jy++)
		for(int ix = 0; ix < out.nx(); ix++)
			fo << out.xx(ix) << "\t" << out.yy(jy) << "\t"
			   << out.HS(ix,jy) << "\t" << out.Ux(ix,jy) << "\t" << out.Uy(ix,jy) << "\t"
			   << out.CC(ix,jy) << "\t" << out.HH(ix,jy) << "\t" << out.SS(ix,jy) << "\n";
	}

	Field_WM2D back;
	back.Build(MakeGrid("xyh", NX, NY), FlatSubstrate(), 1.0);

	Input::I_Restart inp;
	inp.Enable = true; inp.ASCII = true; inp.FileName = tf.path;
	inp.RestartTime = 0; inp.FileIndex = 0;
	REQUIRE(C_Restart::ReadAscii_WM2D(inp, back, 1.0) == true);

	for(int ix = 2; ix <= G->nx() - 3; ix++)
	for(int jy = 2; jy <= G->ny() - 3; jy++)
	{
		Zface zf = G->At<Zface>(ix, jy);
		REQUIRE(back.HH(zf) == doctest::Approx(field.HH(zf)).epsilon(1e-13));
	}
}


//==================================================================
//  A malformed restart file is refused, not half read
//==================================================================
//  A row with fewer numbers than a field*.dat has, or a .vts block whose
//  byte count is not the size of the grid, is refused. Read anyway, the
//  numbers not in the file would be taken from uninitialized memory (so the
//  restart would depend on the compiler and the optimization level), and a
//  block longer than the grid would be written past the end of its array.
TEST_CASE("Restart: a short row or a wrong-sized block is refused")
{
	const int NX = 16;
	Input::I_Restart inp;
	inp.Enable = true; inp.RestartTime = 0; inp.FileIndex = 0;

	SUBCASE("ASCII 1D, a row of five numbers")
	{
		Field_WM1D field;
		field.Build(MakeGrid("xh", NX, 1), FlatSubstrate(), 1.0);
		const int nrow = field.GetGrid().nx() - 2;
		TempFile tf("test_restart_1d_short.dat");
		{
			std::ofstream fo(tf.path.c_str());
			fo << "xx \t HS \t Ux \t CC \t HH \t SS\n";
			for(int r = 0; r < nrow; r++)
			{
				fo << 0.0 << "\t" << 1.0 << "\t" << 0.0 << "\t" << 0.5 << "\t" << 1.0;
				if(r != nrow / 2) fo << "\t" << 0.0;
				fo << "\n";
			}
		}
		inp.ASCII = true; inp.FileName = tf.path;
		CHECK(C_Restart::ReadAscii_WM1D(inp, field, 1.0) == false);
	}

	SUBCASE("ASCII 2D, a row of seven numbers")
	{
		const int NY = 5;
		Field_WM2D field;
		field.Build(MakeGrid("xyh", NX, NY), FlatSubstrate(), 1.0);
		const int nrow = (field.GetGrid().nx() - 2) * (field.GetGrid().ny() - 2);
		TempFile tf("test_restart_2d_short.dat");
		{
			std::ofstream fo(tf.path.c_str());
			fo << "xx \t yy \t HS \t Ux \t Uy \t CC \t HH \t SS\n";
			for(int r = 0; r < nrow; r++)
			{
				fo << "0\t0\t1\t0\t0\t0.5\t1";
				if(r != nrow - 1) fo << "\t0";
				fo << "\n";
			}
		}
		inp.ASCII = true; inp.FileName = tf.path;
		CHECK(C_Restart::ReadAscii_WM2D(inp, field, 1.0) == false);
	}

	SUBCASE("binary 1D, blocks of the wrong size or cut short")
	{
		Field_WM1D field;
		field.Build(MakeGrid("xh", NX, 1), FlatSubstrate(), 1.0);
		const int nx = field.GetGrid().nx() - 2;
		const size_t np = (size_t)nx * 2;
		const std::vector<double> v(2 * np, 0.5);
		TempFile tf("test_restart_1d_block.vts");
		inp.ASCII = false; inp.FileName = tf.path;
		const auto write = [&](unsigned int bytes, size_t nblocks, size_t ndata, bool whole = true)
		{
			std::ofstream fo(tf.path.c_str(), std::ios::binary);
			fo << "<YuzenData Nx=\"" << nx << "\" Nz=\"2\" />\n"
			   << "<AppendedData encoding=\"raw\">\n   _";
			for(size_t b = 0; b < nblocks; b++)
			{
				fo.write((const char*)&bytes, sizeof(bytes));
				fo.write((const char*)v.data(), (std::streamsize)(ndata * sizeof(double)));
			}
			if(whole) fo << "\n</AppendedData>\n</VTKFile>\n";
		};
		//  The right file is read ...
		write((unsigned int)(np * sizeof(double)), 3, np);
		CHECK(C_Restart::Read(inp, field, 1.0) == true);
		//  ... a block that claims more than the grid holds is not ...
		write((unsigned int)(2 * np * sizeof(double)), 3, 2 * np);
		CHECK(C_Restart::Read(inp, field, 1.0) == false);
		//  ... and neither is a file that ends inside its last block.
		write((unsigned int)(np * sizeof(double)), 3, np - 1, false);
		CHECK(C_Restart::Read(inp, field, 1.0) == false);
		//  ... nor one cut short after the blocks the reader takes
		write((unsigned int)(np * sizeof(double)), 3, np, false);
		CHECK(C_Restart::Read(inp, field, 1.0) == false);
	}
}


//==================================================================
//  The marker search of the binary reader
//==================================================================
//  Extruct(fin, st, ed) is what the binary readers find the <YuzenData .../>
//  header and the start of the appended data with. No restart can tell a
//  wrong result from a right one -- the header's attributes are found with
//  find() whatever trails them, and the text before "   _" is discarded --
//  so it is pinned down here: the text between the markers, the end marker
//  itself not included, the stream left just past it, and a missing start
//  marker reported as a failed stream.
std::string Extruct(std::ifstream &fin, const std::string &st, const std::string &ed,
                    std::string *skipped, size_t limit);

TEST_CASE("Binary restart: the marker search returns what lies between the markers")
{
	TempFile tf("test_restart_extruct.txt");
	{
		std::ofstream fo(tf.path.c_str(), std::ios::binary);
		fo << "junk<YuzenData Nx=\"7\" Nz=\"3\" />rest<AppendedData encoding=\"raw\">\n   _BIN";
	}

	std::ifstream fin(tf.path.c_str(), std::ios::binary);
	std::string skipped;
	CHECK(Extruct(fin, "<YuzenData", "/>", &skipped, 1u << 20) == " Nx=\"7\" Nz=\"3\" ");
	CHECK(skipped == "junk<YuzenData");
	REQUIRE(fin.good());

	//  The stream is just past "/>": the next search starts at "rest".
	CHECK(Extruct(fin, "<AppendedData encoding=\"raw\">", "   _", &skipped, 1u << 20) == "\n");
	CHECK(skipped == "rest<AppendedData encoding=\"raw\">");
	char next[4] = {0, 0, 0, 0};
	fin.read(next, 3);
	CHECK(std::string(next) == "BIN");

	//  A start marker that is not there fails the stream at once.
	std::ifstream fin2(tf.path.c_str(), std::ios::binary);
	CHECK(Extruct(fin2, "<Missing", "/>", &skipped, 1u << 20) == "");
	CHECK(fin2.fail());

	//  ... and one beyond the search limit is not found either.
	std::ifstream fin3(tf.path.c_str(), std::ios::binary);
	CHECK(Extruct(fin3, "<YuzenData", "/>", &skipped, 8) == "");
	CHECK(fin3.fail());
}


//------------------------------------------------------------------
//  field*.state: what a restart needs beyond the field file. Every double
//  must come back to the bit, scalars and arrays alike, and a state file
//  is only taken for the model and the time it was written for.
//------------------------------------------------------------------
#include "RestartState.h"
#include "TimeSteps.h"
#include <cmath>
#include <vector>

TEST_CASE("C_RestartState: written and read back to the bit")
{
	TempFile tf("test_restart_state.state");
	C_RestartState st;
	st.SetText("model", "Well-Mixed-1D");
	st.Set("evaporated", 1.0 / 3.0);
	std::vector<double> two(2);
	two[0] = std::nextafter(1.0, 2.0); two[1] = -1e-300;
	st.Set("outflowed", two);
	std::vector<double> v(1000);
	for(size_t i = 0; i < v.size(); i++) v[i] = std::sin(0.1 * (double)i) * std::exp(-(double)i);
	st.SetVector("height.old", v);
	st.SetVector("empty", std::vector<double>());
	REQUIRE(st.Write(tf.path));

	C_RestartState r;
	REQUIRE(r.Read(tf.path));
	CHECK(r.Text("model") == "Well-Mixed-1D");
	CHECK(r.Get("evaporated") == 1.0 / 3.0);
	CHECK(r.Get("outflowed", 0) == two[0]);
	CHECK(r.Get("outflowed", 1) == two[1]);
	CHECK(r.Get("missing", 0, 7.0) == 7.0);
	REQUIRE(r.Vector("height.old").size() == v.size());
	CHECK(r.Vector("height.old") == v);
	CHECK(r.HasVector("empty"));

	//  A field file is not a state file.
	TempFile bad("test_restart_state_bad.state");
	{ std::ofstream f(bad.path.c_str()); f << "# x y\n1 2 3\n"; }
	CHECK_FALSE(r.Read(bad.path));
}

TEST_CASE("C_RestartState: the state file beside a field file")
{
	CHECK(C_RestartState::PathFor("run/field00012.vts") == "run/field00012.state");
	CHECK(C_RestartState::PathFor("run/field00012.dat") == "run/field00012.state");
	CHECK(C_RestartState::PathFor("run.d/field") == "run.d/field.state");
}

TEST_CASE("C_RestartState: taken only for its model and its time")
{
	TempFile tf("test_restart_load.state");
	C_TimeSteps ts;
	ts.Build(1e-3, 0.0, 1.0, true);
	ts.Add_OutputItem("Field", 0.1);
	ts.Set(0.25);
	C_RestartState st;
	st.SetText("model", "ConcZ-2D");
	st.PutClock(ts, std::vector<std::string>(1, "Field"));
	REQUIRE(st.Write(tf.path));

	//  Load() looks for the state file beside the field file it is given.
	const std::string field = tf.path.substr(0, tf.path.size() - 6) + ".vts";
	std::string msg;
	C_RestartState r;
	CHECK(r.Load(field, "ConcZ-2D", 0.25, msg));
	CHECK(r.Load(field, "ConcZ-2D", 0.250001, msg));		// six digits of a .pvd
	CHECK_FALSE(r.Load(field, "ConcZ-2D", 0.3, msg));		// a fresh start from that state
	CHECK_FALSE(r.Load(field, "Well-Mixed-1D", 0.25, msg));
	CHECK_FALSE(r.Load("no_such_file.vts", "ConcZ-2D", 0.25, msg));
}


//------------------------------------------------------------------
//  Damaged restart files are refused, not trusted (the fuzzer of the
//  readers, yuzen/tests/fuzz/fuzz_restart.cpp, feeds such files). Trusted,
//  each of these would read out of bounds, leave values uninitialized, put
//  a NaN into the field, or throw out of the run.
//------------------------------------------------------------------
namespace {

void WriteBytes(const std::string &path, const std::string &bytes)
{
	std::ofstream f(path.c_str(), std::ios::binary);
	f.write(bytes.data(), (std::streamsize)bytes.size());
}

//  How a whole .vts ends
const std::string kVtsEnd = "\n</AppendedData>\n</VTKFile>\n";

std::string U32(unsigned int n) { return std::string((const char*)&n, 4); }
std::string Doubles(size_t n, double v)
{
	std::vector<double> d(n, v);
	return std::string((const char*)d.data(), n * sizeof(double));
}

//  A well-mixed .vts of nx x ny x 2 points: lens[b] is the byte count
//  written before block b, each block holding ndata values val.
std::string WmVts(int nx, int ny, const unsigned int lens[3], double val, size_t ndata)
{
	std::string s = "<VTKFile>\n<YuzenData Time=\"0\" Model=\"WM2D\" Nx=\"" + std::to_string(nx)
	              + "\" Ny=\"" + std::to_string(ny) + "\" Nz=\"2\" />\n"
	              + "<AppendedData encoding=\"raw\">\n   _";
	for(int b = 0; b < 3; b++) s += U32(lens[b]) + Doubles(ndata, val);
	return s + kVtsEnd;
}

Input::I_Restart RestartOf(const std::string &path, bool ascii)
{
	Input::I_Restart r;
	r.Enable = true; r.ASCII = ascii; r.FileName = path; r.RestartTime = 0; r.FileIndex = 0;
	return r;
}

} // namespace

TEST_CASE("Binary restart (WM): a damaged .vts is refused")
{
	Field_WM1D field;
	field.Build(MakeGrid("xh", 8, 1), FlatSubstrate(), 1.0);
	const int nx = field.GetGrid().nx() - 2;
	const unsigned int bytes = (unsigned int)(nx * 2 * sizeof(double));
	TempFile tf("test_restart_bad.vts");
	Input::I_Restart r = RestartOf(tf.path, false);

	const unsigned int good[3] = {bytes, bytes, bytes};
	WriteBytes(tf.path, WmVts(nx, 1, good, 0.5, nx * 2));
	CHECK(C_Restart::Read(r, field, 1.0));

	//  A block that claims more than the grid holds: read into the array
	//  sized from the grid, it would go past its end.
	const unsigned int big[3] = {bytes * 64, bytes, bytes};
	WriteBytes(tf.path, WmVts(nx, 1, big, 0.5, nx * 2 * 64));
	CHECK_FALSE(C_Restart::Read(r, field, 1.0));

	//  Cut short inside its data
	const std::string cut = WmVts(nx, 1, good, 0.5, nx * 2);
	WriteBytes(tf.path, cut.substr(0, cut.size() - kVtsEnd.size() - 12));
	CHECK_FALSE(C_Restart::Read(r, field, 1.0));

	//  Not a number
	WriteBytes(tf.path, WmVts(nx, 1, good, std::nan(""), nx * 2));
	CHECK_FALSE(C_Restart::Read(r, field, 1.0));

	//  No header, and a grid size that is not a number
	WriteBytes(tf.path, "<VTKFile>\n");
	CHECK_FALSE(C_Restart::Read(r, field, 1.0));
	WriteBytes(tf.path, "<YuzenData Nx=\"8x\" Ny=\"1\" Nz=\"2\" />");
	CHECK_FALSE(C_Restart::Read(r, field, 1.0));
}

TEST_CASE("Binary restart (CZ): a block longer than the grid is refused")
{
	Input::I_Grid g = MakeGrid("xh", 6, 1);
	g.Nz = 4;
	Field_CZ field;
	field.Build(g, FlatSubstrate(), 1.0);
	const int nx = field.GetGrid3().nx(), ny = field.GetGrid3().ny(), nz = field.GetGrid3().nz();
	const size_t np = (size_t)nx * ny * nz;
	TempFile tf("test_restart_bad_cz.vts");
	Input::I_Restart r = RestartOf(tf.path, false);

	const auto vts = [&](unsigned int concBytes) {
		std::string s = "<VTKFile>\n<YuzenData Time=\"0\" Model=\"ConcZ\" Nx=\"" + std::to_string(nx)
		              + "\" Ny=\"" + std::to_string(ny) + "\" Nz=\"" + std::to_string(nz) + "\" />\n"
		              + "<AppendedData encoding=\"raw\">\n   _";
		s += U32(concBytes) + Doubles(concBytes / 8, 0.5);
		s += U32((unsigned int)(3 * np * 8)) + Doubles(3 * np, 0.0);
		s += U32((unsigned int)(3 * np * 8)) + Doubles(3 * np, 1.0);
		return s + kVtsEnd;
	};
	WriteBytes(tf.path, vts((unsigned int)(np * 8)));
	CHECK(C_Restart::Read(r, field, 1.0));
	WriteBytes(tf.path, vts((unsigned int)(np * 8 * 16)));
	CHECK_FALSE(C_Restart::Read(r, field, 1.0));
}

TEST_CASE("ASCII restart (WM): a row that is not six finite numbers is refused")
{
	Field_WM1D field;
	field.Build(MakeGrid("xh", 4, 1), FlatSubstrate(), 1.0);
	const int nrow = field.GetGrid().nx() - 2;
	TempFile tf("test_restart_bad.dat");
	Input::I_Restart r = RestartOf(tf.path, true);

	const auto file = [&](const std::string &row) {
		std::ofstream f(tf.path.c_str());
		f << "xx\tHS\tUx\tCC\tHH\tSS\n";
		for(int i = 0; i < nrow; i++) f << (i == 1 ? row : std::string("0\t1\t0\t0.5\t1\t0")) << "\n";
	};
	file("0\t1\t0\t0.5\t1\t0");
	CHECK(C_Restart::ReadAscii_WM1D(r, field, 1.0));
	file("0\t1\t0\tnan\t1\t0");				// must not go into the field
	CHECK_FALSE(C_Restart::ReadAscii_WM1D(r, field, 1.0));
	file("0\t1\t0\t0.5");					// would leave h and s uninitialized
	CHECK_FALSE(C_Restart::ReadAscii_WM1D(r, field, 1.0));
	file(std::string(3000, '1'));			// a very long line
	CHECK_FALSE(C_Restart::ReadAscii_WM1D(r, field, 1.0));
}

TEST_CASE("C_RestartState: a damaged state file is not continued from")
{
	TempFile tf("test_restart_damaged.state");
	const std::string field = tf.path.substr(0, tf.path.size() - 6) + ".vts";
	const std::string head = "# Yuzen restart state 1\ntext model M\ntime 0\n";
	std::string msg;
	C_RestartState r;

	//  An array longer than the file: allocated as it stands, a length near
	//  2^64 would throw std::length_error out of the run.
	WriteBytes(tf.path, head + "vector field.h 18446744073709551615\nbinary\n");
	CHECK_FALSE(r.Read(tf.path));
	WriteBytes(tf.path, head + "vector field.h 3\nbinary\n" + Doubles(2, 1.0));
	CHECK_FALSE(r.Read(tf.path));

	//  The counters are cast to integers by the models: a NaN, a negative
	//  or a fractional value, or one out of range, was undefined behavior
	//  there; a word was read as 0.
	const char *bad[] = {"step nan\n", "step -1\n", "steadyCount 1e300\n", "adjustCount.height 2.5\n",
	                     "out.Field 0 0 -3 0.1 0\n", "dt 0 0\n", "evaporated inf\n", "volume0 x\n"};
	for(const char *b : bad)
	{
		WriteBytes(tf.path, head + b + "binary\n");
		CHECK_FALSE(r.Load(field, "M", 0.0, msg));
	}
	WriteBytes(tf.path, head + "step 12\ndt 0.1 0.1\nbinary\n");
	CHECK(r.Load(field, "M", 0.0, msg));
}
