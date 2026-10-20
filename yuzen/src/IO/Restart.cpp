/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Restart.cpp : when raw data is written in the AppendedData section of a VTK
*  file, it is written as blocks of [#Bytes][DATA].
*
* #Bytes states the byte count of DATA as an integer. The byte count of that
* integer type itself is given by the "header-type" attribute of the <VTKFile>
* element. "header-type" is either UInt32 or UInt64 (VTK-4.2); the default is
* UInt32.
* As it stands, this code is written on the assumption that sizeof(int) == 4
* and sizeof(double) == 8.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
using namespace yuzen::util;
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
using namespace std;

#include "Restart.h"
#include "CZAsciiLayout.h"

#include "LubFields.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif


//------------------------------------------------------------------
//  Reads fin up to and including the marker st, then returns what follows
//  up to (not including) ed, leaving fin just past ed.
//
//  Each character is compared against the end of what has been read, so the
//  cost is linear in the distance read. Searching the whole growing string
//  after every character would be quadratic: harmless on a real .vts, whose
//  markers sit in the first kilobyte, and hours on a file that does not
//  contain them (a 6 MB field*.dat handed to the binary reader). The search
//  for st also gives up after `limit` bytes, so a file of the wrong kind
//  fails at once.
//------------------------------------------------------------------
static bool EndsWith(const std::string &s, const std::string &tail)
{
	return s.size() >= tail.size()
	    && s.compare(s.size() - tail.size(), tail.size(), tail) == 0;
}

std::string Extruct(std::ifstream &fin, const std::string &st, const std::string &ed,
                    std::string *skipped = 0, size_t limit = 1u << 20)
{
	std::string str = "";
	char ch;
	bool found = false;
	while(str.size() < limit && fin.get(ch))
	{
		str += ch;
		if(EndsWith(str, st)) { found = true; break; }
	}
	if(skipped) *skipped = str;
	if(!found) { fin.setstate(std::ios::failbit); return ""; }

	str = "";
	while(fin.get(ch))
	{
		str += ch;
		if(EndsWith(str, ed)) return str.substr(0, str.size() - ed.size());
	}
	return str;
}


std::string Extruct(std::string &str_in, const std::string &st, const std::string &ed)
{
	size_t p1, p2;

	p1 = str_in.find(st);
	if(p1 == string::npos) return "";
	p2 = str_in.find(ed, p1 + st.length());
	if(p2 == string::npos) return "";

	return str_in.substr(p1 + st.length(), p2-(p1 + st.length()));
}

//------------------------------------------------------------------
//  Extracts the header element <YuzenData .../> of a binary restart.
//------------------------------------------------------------------
static std::string ExtructHeader(std::ifstream &fin)
{
	return Extruct(fin, "<YuzenData", "/>");
}

//------------------------------------------------------------------
//  The guards of the readers. A restart file is an input like any other:
//  it may be cut short, be of another grid or model, or have been edited.
//  No reader trusts it: an appended block whose byte count is longer than
//  the array the grid sizes (which would overflow the heap), a file cut
//  short, a row of a field*.dat that does not parse (which would leave
//  values uninitialized) and a NaN each refuse the file with the reason, on
//  the terminal and in FatalErrors.txt, and the run stops (yuzen/tests/fuzz
//  feeds the readers such files).
//------------------------------------------------------------------

//  A .vts ends with </VTKFile>. One that does not was cut short -- written
//  by a run killed or out of disk space -- even when the blocks the reader
//  takes are whole, and it is not used.
static bool EndsWhole(const FilePath &rfile, const char *tag)
{
	std::ifstream f(rfile.c_str(), std::ios::binary);
	f.seekg(0, std::ios::end);
	const std::streamoff size = f.tellg();
	const std::streamoff n = size < 64 ? size : 64;
	std::string tail((size_t)(n > 0 ? n : 0), '\0');
	if(n > 0)
	{
		f.seekg(size - n);
		f.read(&tail[0], n);
	}
	if(f && tail.find("</VTKFile>") != std::string::npos) return true;
	Global::FatalError<<"Restart("<<tag<<"): "<<rfile<<" does not end with </VTKFile>: it was cut"
	                  <<" short (by a run killed or out of disk space while writing it) and is not used."<<nl;
	return false;
}

//  A non-negative integer attribute of the <YuzenData> header; -1 when it
//  is missing or not a number.
static long HeaderInt(std::string &header, const std::string &name)
{
	const std::string v = Extruct(header, name + "=\"", "\"");
	if(v.empty()) return -1;
	char *end = 0;
	const long n = std::strtol(v.c_str(), &end, 10);
	if(end == v.c_str() || *end != '\0' || n < 0 || n > 100000000L) return -1;
	return n;
}

//  One block of the appended data: its length (4 bytes), which must be
//  `bytes`, then the data into dst (skipped when dst is null).
static bool ReadBlock(std::ifstream &fin, void *dst, size_t bytes, const char *what,
                      const FilePath &file)
{
	unsigned int n = 0;
	if(!fin.read((char*)&n, sizeof(n)))
	{
		Global::FatalError<<"Restart: "<<file<<" ends before its "<<what<<" block."<<nl;
		return false;
	}
	if((size_t)n != bytes)
	{
		Global::FatalError<<"Restart: the "<<what<<" block of "<<file<<" holds "<<n
		                  <<" bytes; the grid of the input needs "<<bytes<<"."<<nl;
		return false;
	}
	if(dst) fin.read((char*)dst, static_cast<std::streamsize>(bytes));
	else    fin.ignore(static_cast<std::streamsize>(bytes));
	if(!fin)
	{
		Global::FatalError<<"Restart: "<<file<<" ends inside its "<<what<<" block."<<nl;
		return false;
	}
	return true;
}

//  Every value finite. `what` names the array in the message.
static bool AllFinite(const std::vector<double> &v, const char *what, const FilePath &file)
{
	for(size_t k = 0; k < v.size(); k++)
		if(!std::isfinite(v[k]))
		{
			Global::FatalError<<"Restart: "<<file<<" holds a "<<what<<" that is not a finite"
			                  <<" number ("<<v[k]<<", entry "<<k<<")."<<nl;
			return false;
		}
	return true;
}

//  The numbers of one row of a field*.dat: exactly ncol of them, all finite.
static bool ParseRow(const std::string &line, double *v, int ncol)
{
	const char *p = line.c_str();
	char *e;
	for(int k = 0; k < ncol; k++)
	{
		v[k] = std::strtod(p, &e);
		if(e == p || !std::isfinite(v[k])) return false;
		p = e;
	}
	while(*p == ' ' || *p == '\t' || *p == '\r') p++;
	return *p == '\0';
}

static bool RowError(const FilePath &file, long row, int ncol)
{
	Global::FatalError<<"Restart: row "<<row<<" of "<<file<<" is not "<<ncol
	                  <<" finite numbers."<<nl;
	return false;
}


//------------------------------------------------------------------
//  The rows of a well-mixed field*.dat after its header line, each one
//  checked to be ncol finite numbers.
//------------------------------------------------------------------
static bool ReadAsciiRows(const FilePath &rfile, int ncol, std::vector<std::vector<double> > &rows)
{
	std::ifstream fin(rfile.c_str());
	if(!fin.good())
	{
		Global::FatalError<<"Restart: cannot open "<<rfile<<"."<<nl;
		return false;
	}
	std::string line;
	std::getline(fin, line);					// the header
	//  A file written with Output.FieldData.WritePotential has the
	//  potential as one more column, last: read past
	const int nread = (line.find("PHI[") != std::string::npos) ? ncol + 1 : ncol;
	long row = 0;
	while(std::getline(fin, line))
	{
		row++;
		std::vector<double> v((size_t)nread);
		if(!ParseRow(line, v.data(), nread)) return RowError(rfile, row, nread);
		v.resize((size_t)ncol);
		rows.push_back(v);
	}
	return true;
}

bool C_Restart::ReadAscii_WM1D(Input::I_Restart &inp, Field_WM1D &field, double H_coef)
{
	FilePath rfile = inp.FileName;

	//  x  h+s  u  c  h  s
	std::vector<std::vector<double> > rows;
	if(!ReadAsciiRows(rfile, 6, rows)) return false;
	const int Nx = (int)rows.size();

	Global::Info << "Grid size: Nx = " << Nx << ", Grid.nx -2 = " << field.GetGrid().nx()-2 << nl;

	if(field.GetGrid().nx() - 2 != Nx)
	{
		Global::FatalError<<"Restart(WM1D): grid size mismatch: "<<rfile<<" has "<<Nx
		                  <<" rows, the grid of the input "<<field.GetGrid().nx() - 2<<"."<<nl;
		return false;
	}

	//------------------------------------------------------------------
	//  Row ix corresponds to the grid x index ix+1.
	//
	//  OutputNeutral::Write / Field_WM1D_Out::ConvertScalar "write the grid x
	//  index g to row g-1". There are two ghost layers on each side, so
	//  row 0 = the Xmin boundary, rows 1..Nx-2 = the cell centers, and
	//  row Nx-1 = the Xmax boundary.
	//
	//  ** Do not use Of<Xedge>() here. **
	//  Of() takes the **linear index**, and in 1D the linear index is the x
	//  index + 1, so Of<Xedge>(ix+1) lands on x index ix and reading and
	//  writing end up one cell apart. Use At<>(), which is given grid
	//  coordinates.
	//------------------------------------------------------------------
	for(int ix = 0; ix < Nx; ix ++)
	{
		Xedge xe = field.GetGrid().At<Xedge>(ix + 1);
		const std::vector<double> &v = rows[(size_t)ix];
		field.HH(xe) = v[4]*H_coef;
		field.SS(xe) = v[5]*H_coef;
		field.CC(xe) = v[3];
	}
	return true;
}


bool C_Restart::ReadAscii_WM2D(Input::I_Restart &inp, Field_WM2D &field, double H_coef)
{
	FilePath rfile = inp.FileName;

	//  x  y  h+s  u  v  c  h  s
	std::vector<std::vector<double> > rows;
	if(!ReadAsciiRows(rfile, 8, rows)) return false;
	const int Nsum = (int)rows.size();

	int Nx = field.GetGrid().nx() - 2;
	int Ny = field.GetGrid().ny() - 2;

	Global::Info << "Grid size: Nx*Ny = " << Nsum << ", (Grid.nx-2)*(*)Grid.ny-2) = " << Nx*Ny << nl;

	if(Nsum != Nx*Ny)
	{
		Global::FatalError<<"Restart(WM2D): grid size mismatch: "<<rfile<<" has "<<Nsum
		                  <<" rows, the grid of the input "<<Nx*Ny<<"."<<nl;
		return false;
	}

	//------------------------------------------------------------------
	//  OutputNeutral::Write(Field_WM2D_Out&) writes with jy as the outer loop
	//  and ix as the inner one (IL = nx*jy + ix). So from the row number kk,
	//      ix = kk % Nx,  jy = kk / Nx
	//  and the grid indices are each of those plus 1.
	//  ** Divide by Nx, not by Ny. ** Dividing by Ny breaks on any grid with
	//  Nx != Ny.
	//------------------------------------------------------------------
	for(int kk = 0; kk < Nsum; kk++)
	{
		int ix_g2 = kk % Nx + 1;
		int jy_g2 = kk / Nx + 1;
		Zface zf = field.GetGrid().At<Zface>(ix_g2, jy_g2);
		const std::vector<double> &v = rows[(size_t)kk];
		field.HH(zf) = v[6]*H_coef;
		field.SS(zf) = v[7]*H_coef;
		field.CC(zf) = v[5];
	}
	return true;
}

//------------------------------------------------------------------
//  The well-mixed .vts: the <YuzenData> header, then the appended blocks
//  Conc, Height (h+s) and Substrate of nx*ny*nz values each.
//------------------------------------------------------------------
static bool ReadVtsWM(const FilePath &rfile, int nx, int ny, const char *tag,
                      std::vector<double> &cc, std::vector<double> &hs, std::vector<double> &ss)
{
	std::ifstream fin;
	fin.open(rfile.c_str(), ios::binary);
	if(!fin.good())
	{
		Global::FatalError<<"Restart("<<tag<<"): cannot open "<<rfile<<"."<<nl;
		return false;
	}

	std::string str = ExtructHeader(fin);
	if(!fin.good())
	{
		Global::FatalError<<"Restart("<<tag<<"): "<<rfile<<" is not a .vts written by Yuzen"
		                  <<" (no <YuzenData> header). To restart from a field*.dat, set"
		                  <<" Restart.ASCII = true."<<nl;
		return false;
	}
	const long Nx = HeaderInt(str, "Nx");
	const long Ny = (ny == 1 && Extruct(str, "Ny=\"", "\"").empty()) ? 1 : HeaderInt(str, "Ny");
	const long Nz = HeaderInt(str, "Nz");

	if(Nx != nx || Ny != ny || Nz != 2)
	{
		Global::FatalError<<"Restart("<<tag<<"): grid size mismatch: "<<rfile<<" has "
		                  <<Nx<<" x "<<Ny<<" x "<<Nz<<" points, the input needs "
		                  <<nx<<" x "<<ny<<" x 2."<<nl;
		return false;
	}

	str = Extruct(fin, 	"<AppendedData encoding=\"raw\">", "   _");
	if(!fin.good())
	{
		Global::FatalError<<"Restart("<<tag<<"): "<<rfile<<" has no appended data block."<<nl;
		return false;
	}
	if(!EndsWhole(rfile, tag)) return false;

	const size_t np = (size_t)nx*ny*2;
	cc.assign(np, 0); hs.assign(np, 0); ss.assign(np, 0);
	return ReadBlock(fin, cc.data(), np*sizeof(double), "Conc", rfile)
	    && ReadBlock(fin, hs.data(), np*sizeof(double), "Height", rfile)
	    && ReadBlock(fin, ss.data(), np*sizeof(double), "Substrate", rfile)
	    && AllFinite(cc, "concentration", rfile)
	    && AllFinite(hs, "height", rfile)
	    && AllFinite(ss, "substrate height", rfile);
}

bool C_Restart::Read(Input::I_Restart &inp, Field_WM1D &field, double H_coef)
{
	FilePath rfile = inp.FileName;
	const int Nx = field.GetGrid().nx()-2;
	std::vector<double> cc, hs, ss;
	if(!ReadVtsWM(rfile, Nx, 1, "WM1D", cc, hs, ss)) return false;

	Xedge xe;
	int ix, ind;
	for(int ii = field.GetGrid().imin(); ii <= field.GetGrid().imax(); ++ii)
	{
		xe = field.GetGrid().Of<Xedge>(ii);
		ix = field.GetGrid().ix(xe);

		if(field.GetGrid().CheckGhost(xe) == C_Grid1D::Outside) continue;

		ind = (Nx) + (ix-1);
		field.CC(xe) = cc[ind];
		field.SS(xe) = ss[ind]*H_coef;
		field.HH(xe) = (hs[ind] - ss[ind])*H_coef;
	}
	return true;
}

bool C_Restart::Read(Input::I_Restart &inp, Field_WM2D &field, double H_coef)
{
	FilePath rfile = inp.FileName;
	const int Nx = field.GetGrid().nx()-2, Ny = field.GetGrid().ny()-2;
	std::vector<double> cc, hs, ss;
	if(!ReadVtsWM(rfile, Nx, Ny, "WM2D", cc, hs, ss)) return false;

	Zface zf;
	int ix, jy, ind;
	for(int ii = field.GetGrid().imin(); ii <= field.GetGrid().imax(); ++ii)
	{
		zf = field.GetGrid().Of<Zface>(ii);
		ix = field.GetGrid().ix(zf); jy = field.GetGrid().iy(zf);

		if(field.GetGrid().CheckGhost(zf) == C_Grid2D::Outside) continue;

		ind = (Nx*Ny) + (Nx * (jy-1) ) + (ix-1);
		field.CC(zf) = cc[ind];
		field.SS(zf) = ss[ind]*H_coef;
		field.HH(zf) = (hs[ind] - ss[ind])*H_coef;
	}
	return true;
}

bool C_Restart::Read(Input::I_Restart &inp, Field_CZ &field, double H_coef)
{
	FilePath rfile = inp.FileName;

	std::ifstream fin;
	fin.open(rfile.c_str(), ios::binary);
	if(!fin.good())
	{
		Global::FatalError<<"Restart(CZ): cannot open " << rfile << "." << nl;
		return false;
	}

	std::string str = ExtructHeader(fin);
	if(!fin.good() || str.find("Model=\"ConcZ\"") == std::string::npos)
	{
		Global::FatalError<<"Restart(CZ): " << rfile << " is not a ConcZ .vts written by"
		                  <<" Yuzen (no <YuzenData Model=\"ConcZ\"> header). To restart"
		                  <<" from a field*.dat, set Restart.ASCII = true." << nl;
		return false;
	}

	const long Nx = HeaderInt(str, "Nx");
	const long Ny = HeaderInt(str, "Ny");
	const long Nz = HeaderInt(str, "Nz");

	if(field.GetGrid3().nx() != Nx || field.GetGrid3().ny() != Ny || field.GetGrid3().nz() != Nz)
	{
		Global::FatalError<<"Restart(CZ): grid size mismatch: the file has "
		                  <<Nx<<" x "<<Ny<<" x "<<Nz<<" points including ghost layers, the input "
		                  <<field.GetGrid3().nx()<<" x "<<field.GetGrid3().ny()<<" x "<<field.GetGrid3().nz()<<"."<<nl;
		return false;
	}

	//  Files written before the thickness arrays were added carry the film
	//  only in their point coordinates, laid out as z = s + h*kz/nz; they are
	//  still read, by inverting exactly that.
	std::string xml;
	str = Extruct(fin, "<AppendedData encoding=\"raw\">", "   _", &xml);
	if(!fin.good())
	{
		Global::FatalError<<"Restart(CZ): " << rfile << " has no appended data block." << nl;
		return false;
	}
	if(!EndsWhole(rfile, "CZ")) return false;
	const bool hasThickness = xml.find("Name=\"FilmThickness\"") != std::string::npos;

	const size_t np = (size_t)Nx*Ny*Nz;
	std::vector<double> cc(np), coord(3*np), thick, subst;
	const size_t bytes = np*sizeof(double);
	if(!ReadBlock(fin, cc.data(), bytes, "Conc", rfile)
	   || !ReadBlock(fin, 0, 3*bytes, "Velocity", rfile)
	   || !ReadBlock(fin, coord.data(), 3*bytes, "Coordinates", rfile))
		return false;
	if(hasThickness)
	{
		thick.resize(np); subst.resize(np);
		if(!ReadBlock(fin, thick.data(), bytes, "FilmThickness", rfile)
		   || !ReadBlock(fin, subst.data(), bytes, "Substrate", rfile))
			return false;
	}
	if(!AllFinite(cc, "concentration", rfile) || !AllFinite(coord, "coordinate", rfile)
	   || !AllFinite(thick, "film thickness", rfile) || !AllFinite(subst, "substrate height", rfile))
		return false;
	fin.close();

	int ix, jy, kz;
	for(Zface zf = field.GetGrid2().First<Zface>(); zf != field.GetGrid2().Last<Zface>(); ++zf)
	{
		ix = field.GetGrid2().ix(zf);
		jy = field.GetGrid2().iy(zf);
		const size_t col = (size_t)Nx*jy + ix;
		double hh, ss;
		if(hasThickness)
		{
			hh = thick[col];
			ss = subst[col];
		}
		else
		{
			const double z0   = coord[3*col + 2];
			const double ztop = coord[3*((size_t)Nx*Ny*(Nz-1) + col) + 2];
			ss = z0;
			hh = (ztop - z0) * (double)Nz / (double)(Nz - 1);
		}
		field.SS(zf) = ss * H_coef;
		field.HH(zf) = hh * H_coef;
	}
	for(Cell cl = field.GetGrid3().First<Cell>(); cl != field.GetGrid3().Last<Cell>(); ++cl)
	{
		ix = field.GetGrid3().ix(cl);
		jy = field.GetGrid3().iy(cl);
		kz = field.GetGrid3().iz(cl);
		field.CC(cl) = cc[((size_t)Ny*Nx)*kz + (size_t)Nx*jy + ix];
	}
	return true;
}

//------------------------------------------------------------------
//  ASCII restart of a CZ field, from a field*.dat.
//
//  OutputNeutral::Write(Field_CZ&) writes the points of cz_ascii::Points --
//  wall face, cell centers, wall face along each axis -- and this reader walks
//  the same points, so each row comes back to the cell it came from. The file
//  is written with 17 significant digits, so the state comes back as written:
//  the real cells exactly, the ghost layer next to each wall to rounding
//  (it is recovered from the mean written on the face). The outer ghost
//  layer is not in the file; it is left as the initial field had it and set
//  by the boundary conditions at the first solve.
//------------------------------------------------------------------
bool C_Restart::ReadAscii_CZ(Input::I_Restart &inp, Field_CZ &field, double H_coef)
{
	FilePath rfile = inp.FileName;
	std::ifstream fin(rfile.c_str());
	if(!fin.good())
	{
		Global::FatalError<<"Restart(CZ): cannot open " << rfile << "." << nl;
		return false;
	}

	const C_Grid2D &G2 = field.GetGrid2();
	const C_Grid3D &G3 = field.GetGrid3();
	const bool twoD = (G3.ny() != 1);
	const int  ncol = twoD ? 10 : 8;

	std::string line;
	std::getline(fin, line);
	if(line.empty() || line[0] != '#')
	{
		Global::FatalError<<"Restart(CZ): " << rfile << " does not start with the header of a"
		                  <<" ConcZ field*.dat." << nl;
		return false;
	}

	//  The rows are the points of cz_ascii::Points along x, y and z, x the
	//  fastest, exactly as OutputNeutral::Write(Field_CZ&) walks them.
	const std::vector<cz_ascii::Pt> px = cz_ascii::Points(G3.IminReal(), G3.ImaxReal(),
		[&](int i){ return G3.XOf<Cell>(i); }, G3.XminReal(), G3.XmaxReal());
	const std::vector<cz_ascii::Pt> py = !twoD
		? std::vector<cz_ascii::Pt>(1, cz_ascii::Pt{0, 0, 0.0})
		: cz_ascii::Points(G3.JminReal(), G3.JmaxReal(),
			[&](int j){ return G3.YOf<Cell>(0, j); }, G3.YminReal(), G3.YmaxReal());
	const std::vector<cz_ascii::Pt> pz = cz_ascii::Points(G3.KminReal(), G3.KmaxReal(),
		[&](int k){ return G3.ZOf<Cell>(0, 0, k); }, G3.ZminReal(), G3.ZmaxReal());
	const size_t nx = px.size(), ny = py.size(), nz = pz.size();

	std::vector<double> rc, rh, rs;		// c, h and s of every row, in row order
	rc.reserve(nx*ny*nz);
	long nrow = 0;
	const long expected = (long)(nx*ny*nz);
	while(nrow < expected && std::getline(fin, line))
	{
		double v[10];
		int k = 0;
		const char *p = line.c_str();
		char *e;
		for(; k < ncol; k++)
		{
			v[k] = strtod(p, &e);
			if(e == p || !std::isfinite(v[k])) break;
			p = e;
		}
		if(k != ncol) break;
		//  1-D: x z c u w h+s h s     2-D: x y z c u v w h+s h s
		rc.push_back(v[twoD ? 3 : 2]);
		rh.push_back(v[twoD ? 8 : 6]);
		rs.push_back(v[twoD ? 9 : 7]);
		nrow++;
	}
	if(nrow != expected || (std::getline(fin, line) && !line.empty()))
	{
		Global::FatalError<<"Restart(CZ): " << rfile << " does not match the grid: expected "
		                  <<expected<<" rows of "<<ncol<<" finite numbers, row "<<nrow + 1
		                  <<" is not one of them." << nl;
		return false;
	}

	//  A row at a cell center is that cell. A row on a wall face, the
	//  substrate or the free surface is the mean of the cells either side,
	//  one of them a ghost; the rows are taken in order of how many such faces
	//  they sit on, so that when a row is reached every cell of its mean but
	//  one is already known, and that one -- the ghost next to the wall -- is
	//  recovered from it. The outer ghost layer is not in the file; it is left
	//  as the initial field had it and set at the first solve.
	const auto row = [&](size_t i, size_t j, size_t k){ return (k*ny + j)*nx + i; };
	for(int nface = 0; nface <= 3; nface++)
	for(size_t k = 0; k < nz; k++)
	for(size_t j = 0; j < ny; j++)
	for(size_t i = 0; i < nx; i++)
	{
		const cz_ascii::Pt &X = px[i], &Y = py[j], &Z = pz[k];
		if(!cz_ascii::IsCenter(X) + !cz_ascii::IsCenter(Y) + !cz_ascii::IsCenter(Z) != nface) continue;
		//  The unknown cell is the ghost on each face axis.
		const int gx = (X.a < G3.IminReal() ? X.a : X.b);
		const int gy = (Y.a < (twoD ? G3.JminReal() : 0) ? Y.a : Y.b);
		const int gz = (Z.a < G3.KminReal() ? Z.a : Z.b);
		const int xs[2] = {X.a, X.b}, ys[2] = {Y.a, Y.b}, zs[2] = {Z.a, Z.b};

		//  c: the mean over up to 8 cells (each counted with its weight).
		double known = 0;
		for(int a = 0; a < 2; a++) for(int b = 0; b < 2; b++) for(int c = 0; c < 2; c++)
		{
			const bool ghost = (xs[a] == gx || cz_ascii::IsCenter(X))
			                && (ys[b] == gy || cz_ascii::IsCenter(Y))
			                && (zs[c] == gz || cz_ascii::IsCenter(Z));
			if(!ghost) known += field.CC(G3.At<Cell>(xs[a], ys[b], zs[c]));
		}
		const int nsame = (1 << (3 - nface));	// copies of the unknown in the 8
		field.CC(G3.At<Cell>(gx, gy, gz)) = (8.0*rc[row(i,j,k)] - known) / nsame;

		//  h and s: the same, over the in-plane cells, once per column.
		if(k != 1 || !cz_ascii::IsCenter(Z)) continue;
		const int nf2 = !cz_ascii::IsCenter(X) + !cz_ascii::IsCenter(Y);
		double kh = 0, ks = 0;
		for(int a = 0; a < 2; a++) for(int b = 0; b < 2; b++)
		{
			const bool ghost = (xs[a] == gx || cz_ascii::IsCenter(X))
			                && (ys[b] == gy || cz_ascii::IsCenter(Y));
			if(ghost) continue;
			const Zface z2 = G2.At<Zface>(xs[a], ys[b]);
			kh += field.HH(z2);
			ks += field.SS(z2);
		}
		const int n2 = (1 << (2 - nf2));
		const Zface zg = G2.At<Zface>(gx, gy);
		field.HH(zg) = (4.0*rh[row(i,j,k)]*H_coef - kh) / n2;
		field.SS(zg) = (4.0*rs[row(i,j,k)]*H_coef - ks) / n2;
	}
	return true;
}
