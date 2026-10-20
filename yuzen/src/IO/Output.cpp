/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Output.cpp : when raw data is written in the AppendedData section of a VTK
*  file, it is written as blocks of [#Bytes][DATA].
*
* #Bytes states the byte count of DATA as an integer. The byte count of that
* integer type itself is given by the "header-type" attribute of the <VTKFile>
* element. "header-type" is either UInt32 or UInt64 (VTK-4.2); the default is
* UInt32.
* As it stands, this code is written on the assumption that sizeof(int) == 4
* and sizeof(double) == 8
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#include "YuzenUtil.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <system_error>
using namespace std;

#include "MyMath.h"
#include "Output.h"
#include "CZAsciiLayout.h"
#include "NonDimensionalNumbers.h"

using namespace yuzen::util;
using namespace Global;

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif


void OutVelocity(ofstream& fout, double *vx, int nx)
{
	unsigned int size_v;
	unsigned int d_size = 8;

	int nz = 2;
	size_v = (nx * nz) * d_size * 3;

	double ux, uy, uz;
	int ix, jy, kz;
	fout.write((char*)&size_v, sizeof(int));
	for(kz=0; kz < nz ; kz++){
		for(jy=0; jy < 1 ; jy++){
			for(ix = 0; ix < nx; ix++){
				ux = vx[ix];
				uy = 0;
				uz = 0;
				fout.write((char*)&ux, static_cast<std::streamsize>(d_size));
				fout.write((char*)&uy, static_cast<std::streamsize>(d_size));
				fout.write((char*)&uz, static_cast<std::streamsize>(d_size));
			}
		}
	}
}

void OutVelocity(ofstream& fout, double *vx, double *vy, int nx, int ny)
{
	unsigned int size_v;
	unsigned int d_size = 8;

	int nz = 2;
	size_v = (nx * ny * nz) * d_size * 3;

	double ux, uy, uz;
	int ix, jy, kz,s1;
	fout.write((char*)&size_v, sizeof(int));
	for(kz=0; kz < nz ; kz++){
		for(jy=0; jy < ny ; jy++){
			for(ix = 0; ix < nx; ix++){
				s1 = nx*jy +ix;
				ux = vx[s1];
				uy = vy[s1];
				uz = 0;
				fout.write((char*)&ux, static_cast<std::streamsize>(d_size));
				fout.write((char*)&uy, static_cast<std::streamsize>(d_size));
				fout.write((char*)&uz, static_cast<std::streamsize>(d_size));
			}
		}
	}
}

void OutVelocity(ofstream& fout, const VectorField2D<> &vv)
{
	unsigned int  size_v;
	unsigned int  d_size = 8;
	int nx = vv.Grid().nx();
	int ny = vv.Grid().ny();
	int nz = 2;
	size_v = (nx * ny * nz) * d_size * 3;

	std::vector<double> ux(nx*ny*nz, 0.0);
	std::vector<double> uy(nx*ny*nz, 0.0);
	std::vector<double> uz(nx*ny*nz, 0.0);

	int ix, jy, kz, s1, s2;
	Xedge xe; Yedge ye;
	for(Zface zf = vv.Grid().First<Zface>(); zf != vv.Grid().Last<Zface>(); ++zf)
	{
		ix = vv.Grid().ix(zf);
		jy = vv.Grid().iy(zf);
		xe = vv.Grid().At<Xedge>(ix,jy);
		ye = vv.Grid().At<Yedge>(ix,jy);

		kz = 0;
		ux[(ny*nx)*kz + nx*jy +ix] = vv.x(ye);		uy[(ny*nx)*kz + nx*jy +ix] = vv.y(xe);
		uz[(ny*nx)*kz + nx*jy +ix] = 0;

		kz = 1;
		s1 = (ny*nx)*(kz  ) + nx*jy +ix;
		s2 = (ny*nx)*(kz-1) + nx*jy +ix;
		ux[s1] = ux[s2];
		uy[s1] = uy[s2];
		uz[s1] = 0;
	}

	fout.write((char*)&size_v, sizeof(int));
	for(kz=0; kz < nz ; kz++){
		for(jy=0; jy < ny ; jy++){
			for(ix = 0; ix < nx; ix++){
				s1 = (ny*nx)*kz + nx*jy +ix;
				fout.write((char*)&ux[s1], static_cast<std::streamsize>(d_size));
				fout.write((char*)&uy[s1], static_cast<std::streamsize>(d_size));
				fout.write((char*)&uz[s1], static_cast<std::streamsize>(d_size));
			}
		}
	}
}


//----------------------------
//   Writing a file safely
//----------------------------

//  Reports a file that could not be written, on standard error and in
//  FatalErrors.txt. The caller stops the run (StopOnOutputError).
static bool WriteFailed(const std::string &path, const std::string &why)
{
	Global::FatalError<<"Could not write "<<path<<": "<<why<<"."<<nl;
	return false;
}

//  Closes a file written under tmp and, if every write reached the disk,
//  renames it to path. A field file, a .pvd or a state file is therefore
//  either the whole file or not there: a run killed while writing, or out of
//  disk space, leaves at worst the tmp file beside it, which no reader opens.
static bool Commit(std::ofstream &fout, const std::string &tmp, const std::string &path)
{
	fout.close();		// the last of the buffer is written here, and can fail here
	if(fout.fail())
	{
		const std::string why = System::LastIOError();
		std::remove(tmp.c_str());
		return WriteFailed(path, why);
	}
	std::string why;
	if(!System::ReplaceFile(tmp, path, why))
	{
		std::remove(tmp.c_str());
		return WriteFailed(path, why);
	}
	return true;
}

#ifdef _DEF_YUZEN_
namespace yuzen {
#endif
void StopOnOutputError()
{
	Global::FatalError<<"The run is stopped: its output could not be written into "<<JobInfo.ResultDir()
	                  <<" (full disk, quota, file size limit or permissions). The field files already"
	                  <<" there are complete, and the last one with its .state can be restarted from."<<nl;
	Global::FatalError.Abort();
}
#ifdef _DEF_YUZEN_
}
#endif


//----------------------------
//   C_PVD_File
//----------------------------


//  The lines of a file, empty if there is none. A last line without its
//  newline is the row a run was writing when it was killed or ran out of
//  space, and is left out. Only a regular file is read: a device in its
//  place (/dev/full, /dev/zero) never ends and has no newline, and reading
//  it took every byte of memory there was before the run could report that
//  it cannot write there (found by error_paths under AddressSanitizer).
static std::vector<std::string> ReadLines(const char *path)
{
	std::vector<std::string> lines;
	std::error_code ec;
	if(!std::filesystem::is_regular_file(std::filesystem::path(path), ec)) return lines;
	std::ifstream fin(path, std::ios::binary);
	std::string l;
	while(std::getline(fin, l))
	{
		if(fin.eof()) break;		// no newline after it
		if(!l.empty() && l[l.size() - 1] == '\r') l.erase(l.size() - 1);
		lines.push_back(l);
	}
	return lines;
}

//  The time a row of a timeline file starts with; false for a header line.
static bool RowTime(const std::string &l, double &t)
{
	const char *p = l.c_str();
	char *e;
	t = strtod(p, &e);
	return e != p;
}

void C_PVD_File::Build(C_TimeSteps &tt, NonDimensionalNumbers &nd, bool dimensioned, const FilePath &resultDir, const std::string &fname, const std::string &head, double st_time, double dt, int fst)
{
	C_Output_Interface::Build(tt, nd, dimensioned, resultDir, fname, ".pvd");

	_header = head;
	_delta_t = dt;
	_fnum_st = fst;
	_start_time = st_time;

	_fnum = _fnum_st;
	_time = _start_time;

	//  The file is named after the last part of fname: a JobName that is a
	//  path ("runs/case1") names a result directory runs/case1 and a file
	//  case1.pvd inside it, not runs/case1/runs/case1.pvd
	std::string base = fname;
	while(base.size() > 1 && (base[base.size() - 1] == '/' || base[base.size() - 1] == '\\'))
		base.erase(base.size() - 1);
	const size_t slash = base.find_last_of("/\\");
	if(slash != std::string::npos && slash + 1 < base.size()) base = base.substr(slash + 1);
	_fpath =  resultDir / base;
	_fpath += _ext;

	_previous.clear();
	{
		const std::vector<std::string> old = ReadLines(_fpath.c_str());
		for(size_t k = 0; k < old.size(); k++)
			if(old[k].find("<DataSet") != std::string::npos) _previous.push_back(old[k]);
	}
	_entries.clear();

	_constructed = true;
	if(!Rewrite()) StopOnOutputError();
}

//  The whole file, the closing tags included, so that it is valid XML
//  after every field file -- ParaView opens the .pvd of a run that is still
//  going or that was killed. It is written under a temporary name and
//  renamed, as the field files are.
bool C_PVD_File::Rewrite()
{
	if(!_constructed) return false;
	const std::string fpath = _fpath, ftmp = fpath + ".tmp";
	errno = 0;
	ofstream fout(ftmp.c_str());
	if(!fout.good()) return WriteFailed(fpath, System::LastIOError());
	//  One write of the whole text: the file is rewritten after every field
	//  output, so a flush per entry (std::endl) made a run with N outputs
	//  cost N^2/2 system calls -- 20 ms per output after 30,000 of them.
	std::string text = "<VTKFile type=\"Collection\">\n  <Collection>\n";
	for(size_t k = 0; k < _entries.size(); k++) { text += _entries[k]; text += '\n'; }
	text += "  </Collection>\n</VTKFile>\n";
	fout.write(text.data(), static_cast<std::streamsize>(text.size()));
	return Commit(fout, ftmp, fpath);
}

bool C_PVD_File::operator++()
{
	return Add(_time);
}

bool C_PVD_File::Add(double timeNow)
{
	if(!_constructed) return false;

	char fname[5000];
	sprintf_s(fname, "%s%05d%s", _header.c_str(), _fnum, ".vts");
	//  12 significant digits: at the default 6, outputs 1e-4 s apart after
	//  t = 100 s would carry the same time stamp.
	std::ostringstream line;
	line<<std::setprecision(12)<<"    <DataSet timestep=\""<<timeNow<<"\" file=\""<<fname<<"\"/>";
	_entries.push_back(line.str());
	if(!Rewrite()) return false;

	_fnum++;
	_time += _delta_t;
	return true;
}

bool C_PVD_File::Resume(double tmax)
{
	if(!_constructed) return false;
	//  An entry is kept when it points at a field file numbered below the
	//  one this run writes next (SetFNum() has been given Restart.FileIndex).
	//  Its time stamp alone cannot tell: it is written to 12 digits, which can
	//  round the time restarted from up past tmax,
	//  and the entry of the very file restarted from would go missing
	//  (tools/proptest.py seed 2283). The time decides only for an entry
	//  whose file name is not one of ours.
	std::vector<std::string> kept;
	const std::string head = "file=\"" + _header;
	for(size_t k = 0; k < _previous.size(); k++)
	{
		const std::string &line = _previous[k];
		const size_t f = line.find(head);
		if(f != std::string::npos)
		{
			const char *s = line.c_str() + f + head.size();
			char *e;
			const long n = strtol(s, &e, 10);
			if(e != s && strncmp(e, ".vts\"", 5) == 0)
			{
				if(n < _fnum) kept.push_back(line);
				continue;
			}
		}
		const size_t p = line.find("timestep=\"");
		if(p == std::string::npos) continue;
		if(strtod(line.c_str() + p + 10, 0) <= tmax) kept.push_back(line);
	}
	_previous.clear();
	_entries.insert(_entries.begin(), kept.begin(), kept.end());
	return Rewrite();
}

bool C_PVD_File::Close()
{
	if(!_constructed) return false;
	//  Valid after every entry already; written once more for a run that
	//  wrote no field file
	return Rewrite();
}


//----------------------------
//   OutputVTK
//----------------------------
void OutputVTK::Build(C_TimeSteps &tt, NonDimensionalNumbers &nd, bool dimensioned, const FilePath &resultDir, const std::string &fhead, int fst)
{
	C_Output_Interface::Build(tt, nd, dimensioned, resultDir, fhead, ".vts", fst);

	_constructed = true;
}

std::string OutputVTK::FileName() const
{
	char wfile_head[2000];

	sprintf_s(wfile_head, "%s%05d", _head.c_str(), _fnum);
	FilePath fpath = _dir / wfile_head;
	fpath += _ext;
	return fpath.c_str();
}

bool OutputVTK::Write(Field_WM1D_Out &ff)
{
	if(!_constructed) return false;

	const std::string fpath = FileName(), ftmp = fpath + ".tmp";
	errno = 0;
	ofstream fout(ftmp.c_str(), ios::binary);
	if(!fout.good()) return WriteFailed(fpath, System::LastIOError());


	unsigned int size_s,size_v;
	char extent[3000];
	unsigned int offset = 0;
	unsigned int i_size = 4;
	unsigned int d_size = 8;

	int nx = ff.nx();
	int ny = 1;
	int nz = ff.nz();

	size_s = (nx * ny * nz) * d_size;
	size_v = (nx * ny * nz) * d_size * 3;

	sprintf_s(extent,"%d %d %d %d %d %d",1,nx,1,ny,1,nz);

//	fout<<"<VTKFile type=\"StructuredGrid\" GhostLevel=\"0\" byte_order=\"LittleEndian\" version=\"1.0\" header_type=\"UInt32\">"<<endl;	//[2013.3.25] not supported by ParaView yet
	fout<<"<VTKFile type=\"StructuredGrid\" GhostLevel=\"0\" byte_order=\"LittleEndian\">"<<endl;

	//Adds a tag of our own to the VTK file
	double time = _timeStep->NowTime();
	if(this->_dimensioned) time *= _NonDim->TimeScale();
	fout<<"<YuzenData Time=\""<<std::scientific<<time<<"\" Model=\"WM2D\" Nx=\""<<nx<<"\" Ny=\""<<ny<<"\" Nz=\""<<nz<<"\" ";
	fout<<"/>"<<endl;

	fout<<"<StructuredGrid WholeExtent=\""<<extent<<"\">"<<endl;
	fout<<"<Piece Extent=\""<<extent<<"\">"<<endl;
	fout<<"<PointData>"<<endl;
	fout<<"<DataArray Name=\"Conc\" type=\"Float64\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_s + i_size;
	fout<<"<DataArray Name=\"Height\" type=\"Float64\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_s + i_size;
	fout<<"<DataArray Name=\"Substrate\" type=\"Float64\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_s + i_size;

	fout<<"<DataArray Name=\"Velocity\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;

	fout<<"<DataArray Name=\"V_Ca\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	fout<<"<DataArray Name=\"V_Ma\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	fout<<"<DataArray Name=\"V_Bo\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	fout<<"<DataArray Name=\"V_Pg\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	fout<<"<DataArray Name=\"V_Sp\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	//  The driving potential (WritePotential): its data after the
	//  coordinates, so that the blocks before them stay where a reader
	//  taking them in order (the restart) finds them
	if(ff.HasPotential())
		fout<<"<DataArray Name=\"Potential\" type=\"Float64\" format=\"appended\" offset=\""
		    <<offset + size_v + i_size<<"\"/>"<<endl;


	fout<<"</PointData>"<<endl;
	fout<<"<Points>"<<endl;
	fout<<"<DataArray Name=\"Coordinate\" type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	fout<<"</Points>"<<endl;
	fout<<"</Piece>"<<endl;
	fout<<"</StructuredGrid>"<<endl;

	fout<<"<AppendedData encoding=\"raw\">"<<endl;

	fout<<"   _";

	double xx, yy, zz;

	int ix, jy, kz;

	fout.write((char*)&size_s, sizeof(int));	//First writes the length of the data that follows, in 4 bytes (sizeof(int))
	fout.write((char*)ff.C_ptr(), static_cast<std::streamsize>(size_s));	//Writes the data

	fout.write((char*)&size_s, sizeof(int));
	fout.write((char*)ff.HS_ptr(), static_cast<std::streamsize>(size_s));

	fout.write((char*)&size_s, sizeof(int));
	fout.write((char*)ff.S_ptr(), static_cast<std::streamsize>(size_s));

	OutVelocity(fout, ff.Ux_ptr(), nx);

	OutVelocity(fout, ff.Ux_Ca_ptr(), nx);
	OutVelocity(fout, ff.Ux_Ma_ptr(), nx);
	OutVelocity(fout, ff.Ux_Bo_ptr(), nx);
	OutVelocity(fout, ff.Ux_Pg_ptr(), nx);
	OutVelocity(fout, ff.Ux_Sp_ptr(), nx);

	// Coordinates
	double zeta;
	double hs, ss;
	fout.write((char*)&size_v, sizeof(int));
	for(kz=0; kz < nz ; kz++){
		zeta = (double)kz;
		for(jy=0; jy < ny ; jy++){
			for(ix = 0; ix < nx; ix++){
				xx = ff.xx(ix);
				yy = 0;
				hs = ff.HS(ix);
				ss = ff.SS(ix);
				fout.write((char*)&xx, static_cast<std::streamsize>(d_size));
				fout.write((char*)&yy, static_cast<std::streamsize>(d_size));
				zz = ss + (hs-ss)*zeta;
				fout.write((char*)&zz, static_cast<std::streamsize>(d_size));
			}
		}
	}

	if(ff.HasPotential())
	{
		fout.write((char*)&size_s, sizeof(int));
		fout.write((char*)ff.PHI_ptr(), static_cast<std::streamsize>(size_s));
	}

	fout<<endl;
	fout<<"</AppendedData>"<<endl;
	fout<<"</VTKFile>"<<endl;

	if(!Commit(fout, ftmp, fpath)) return false;
	_fnum++;
	return true;
}

bool OutputVTK::Write(Field_WM2D_Out &ff)
{
	if(!_constructed) return false;

	const std::string fpath = FileName(), ftmp = fpath + ".tmp";
	errno = 0;
	ofstream fout(ftmp.c_str(), ios::binary);
	if(!fout.good()) return WriteFailed(fpath, System::LastIOError());

	unsigned int size_s,size_v;
	char extent[3000];
	unsigned int offset = 0;
	unsigned int i_size = 4;
	unsigned int d_size = 8;

	int nx = ff.nx();
	int ny = ff.ny();
	int nz = ff.nz();

	size_s = (nx * ny * nz) * d_size;
	size_v = (nx * ny * nz) * d_size * 3;

	sprintf_s(extent,"%d %d %d %d %d %d",1,nx,1,ny,1,nz);

//	fout<<"<VTKFile type=\"StructuredGrid\" GhostLevel=\"0\" byte_order=\"LittleEndian\" version=\"1.0\" header_type=\"UInt32\">"<<endl;	//[2013.3.25] not supported by ParaView yet
	fout<<"<VTKFile type=\"StructuredGrid\" GhostLevel=\"0\" byte_order=\"LittleEndian\">"<<endl;

	//Adds a tag of our own to the VTK file
	double time = _timeStep->NowTime();
	if(this->_dimensioned) time *= _NonDim->TimeScale();
	fout<<"<YuzenData Time=\""<<std::scientific<<time<<"\" Model=\"WM2D\" Nx=\""<<nx<<"\" Ny=\""<<ny<<"\" Nz=\""<<nz<<"\" ";
	fout<<"/>"<<endl;

	fout<<"<StructuredGrid WholeExtent=\""<<extent<<"\">"<<endl;
	fout<<"<Piece Extent=\""<<extent<<"\">"<<endl;
	fout<<"<PointData>"<<endl;
	fout<<"<DataArray Name=\"Conc\" type=\"Float64\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_s + i_size;
	fout<<"<DataArray Name=\"Height\" type=\"Float64\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_s + i_size;
	fout<<"<DataArray Name=\"Substrate\" type=\"Float64\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_s + i_size;

	fout<<"<DataArray Name=\"Velocity\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;

	fout<<"<DataArray Name=\"V_Ca\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	fout<<"<DataArray Name=\"V_Ma\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	fout<<"<DataArray Name=\"V_Bo\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	fout<<"<DataArray Name=\"V_Pg\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	fout<<"<DataArray Name=\"V_Sp\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	//  The driving potential (WritePotential): its data after the
	//  coordinates, so that the blocks before them stay where a reader
	//  taking them in order (the restart) finds them
	if(ff.HasPotential())
		fout<<"<DataArray Name=\"Potential\" type=\"Float64\" format=\"appended\" offset=\""
		    <<offset + size_v + i_size<<"\"/>"<<endl;


	fout<<"</PointData>"<<endl;
	fout<<"<Points>"<<endl;
	fout<<"<DataArray Name=\"Coordinate\" type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	fout<<"</Points>"<<endl;
	fout<<"</Piece>"<<endl;
	fout<<"</StructuredGrid>"<<endl;

	fout<<"<AppendedData encoding=\"raw\">"<<endl;

	fout<<"   _";

	double xx, yy, zz;

	int ix, jy, kz;

	fout.write((char*)&size_s, sizeof(int));	//First writes the length of the data that follows, in 4 bytes (sizeof(int))
	fout.write((char*)ff.C_ptr(), static_cast<std::streamsize>(size_s));	//Writes the data

	fout.write((char*)&size_s, sizeof(int));
	fout.write((char*)ff.HS_ptr(), static_cast<std::streamsize>(size_s));

	fout.write((char*)&size_s, sizeof(int));
	fout.write((char*)ff.S_ptr(), static_cast<std::streamsize>(size_s));

	OutVelocity(fout, ff.Ux_ptr(), ff.Uy_ptr(), nx, ny);

	OutVelocity(fout, ff.Ux_Ca_ptr(), ff.Uy_Ca_ptr(), nx, ny);
	OutVelocity(fout, ff.Ux_Ma_ptr(), ff.Uy_Ma_ptr(), nx, ny);
	OutVelocity(fout, ff.Ux_Bo_ptr(), ff.Uy_Bo_ptr(), nx, ny);
	OutVelocity(fout, ff.Ux_Pg_ptr(), ff.Uy_Pg_ptr(), nx, ny);
	OutVelocity(fout, ff.Ux_Sp_ptr(), ff.Uy_Sp_ptr(), nx, ny);

	// Coordinates
	double zeta;
	double hs, ss;
	fout.write((char*)&size_v, sizeof(int));
	for(kz=0; kz < nz ; kz++){
		zeta = (double)kz;
		for(jy=0; jy < ny ; jy++){
			for(ix = 0; ix < nx; ix++){
				xx = ff.xx(ix);
				yy = ff.yy(jy);
				hs = ff.HS(ix,jy);
				ss = ff.SS(ix,jy);
				fout.write((char*)&xx, static_cast<std::streamsize>(d_size));
				fout.write((char*)&yy, static_cast<std::streamsize>(d_size));
				zz = ss + (hs-ss)*zeta;
				fout.write((char*)&zz, static_cast<std::streamsize>(d_size));
			}
		}
	}

	if(ff.HasPotential())
	{
		fout.write((char*)&size_s, sizeof(int));
		fout.write((char*)ff.PHI_ptr(), static_cast<std::streamsize>(size_s));
	}

	fout<<endl;
	fout<<"</AppendedData>"<<endl;
	fout<<"</VTKFile>"<<endl;

	if(!Commit(fout, ftmp, fpath)) return false;
	_fnum++;
	return true;
}

bool OutputVTK::Write(Field_CZ &ff)
{
	if(!_constructed) return false;

	const std::string fpath = FileName(), ftmp = fpath + ".tmp";
	errno = 0;
	ofstream fout(ftmp.c_str(), ios::binary);
	if(!fout.good()) return WriteFailed(fpath, System::LastIOError());

	unsigned int  size_s,size_v;
	char extent[3000];
	unsigned int offset = 0;
	unsigned int i_size = 4;
	unsigned int d_size = 8;

	int nx = ff.GetGrid3().nx();
	int ny = ff.GetGrid3().ny();
	int nz = ff.GetGrid3().nz();

	size_s = (nx * ny * nz) * d_size;
	size_v = (nx * ny * nz) * d_size * 3;

	sprintf_s(extent,"%d %d %d %d %d %d",1,nx,1,ny,1,nz);

//	fout<<"<VTKFile type=\"StructuredGrid\" GhostLevel=\"0\" byte_order=\"LittleEndian\" version=\"1.0\" header_type=\"UInt32\">"<<endl;	//[2013.3.25] not supported by ParaView yet
	fout<<"<VTKFile type=\"StructuredGrid\" GhostLevel=\"0\" byte_order=\"LittleEndian\">"<<endl;

	//Adds a tag of our own to the VTK file
	double time = _timeStep->NowTime();
	if(this->_dimensioned) time *= _NonDim->TimeScale();
	fout<<"<YuzenData Time=\""<<std::scientific<<time<<"\" Model=\"ConcZ\" Nx=\""<<nx<<"\" Ny=\""<<ny<<"\" Nz=\""<<nz<<"\" ";
	fout<<"/>"<<endl;

	fout<<"<StructuredGrid WholeExtent=\""<<extent<<"\">"<<endl;
	fout<<"<Piece Extent=\""<<extent<<"\">"<<endl;
	fout<<"<PointData>"<<endl;
	fout<<"<DataArray Name=\"Conc\" type=\"Float64\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_s + i_size;
	fout<<"<DataArray Name=\"Velocity\" type=\"Float64\"  NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset<<"\"/>"<<endl;
	offset += size_v + i_size;
	//  The film thickness and the substrate height, repeated down each column.
	//  They come after the coordinates in the appended block, so their
	//  offsets are past it. The restart reads them back directly: taking them
	//  from the point coordinates instead means undoing the zeta mapping,
	//  which the ghost layers make easy to get wrong.
	const unsigned int offCoord = offset;
	const unsigned int offThick = offCoord + size_v + i_size;
	const unsigned int offSubst = offThick + size_s + i_size;
	fout<<"<DataArray Name=\"FilmThickness\" type=\"Float64\" format=\"appended\" offset=\""<<offThick<<"\"/>"<<endl;
	fout<<"<DataArray Name=\"Substrate\" type=\"Float64\" format=\"appended\" offset=\""<<offSubst<<"\"/>"<<endl;

	fout<<"</PointData>"<<endl;
	fout<<"<Points>"<<endl;
	fout<<"<DataArray Name=\"Coordinate\" type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offCoord<<"\"/>"<<endl;
	fout<<"</Points>"<<endl;
	fout<<"</Piece>"<<endl;
	fout<<"</StructuredGrid>"<<endl;

	fout<<"<AppendedData encoding=\"raw\">"<<endl;

	fout<<"   _";

	double xx, yy, zz;
	std::vector<double> hs(nx*ny, 0.0);
	std::vector<double> ss(nx*ny, 0.0);
	std::vector<double> hh(nx*ny, 0.0);
	std::vector<double> cc(nx*ny*nz, 0.0);
	std::vector<double> ux(nx*ny*nz, 0.0);
	std::vector<double> uy(nx*ny*nz, 0.0);
	std::vector<double> uz(nx*ny*nz, 0.0);

	double L0 = 1, H0 = 1, U0 = 1;
	if(_dimensioned)
	{
		L0 = _NonDim->ExtentScale() * 1e3;
		H0 = _NonDim->HeightScale() * 1e6;
		U0 = _NonDim->VelocityScale();
	}
	//  The through-thickness velocity is made dimensionless with H0/T0, not
	//  L0/T0: in m/s it is U0 * H0/L0.
	const double W0 = _dimensioned ? U0 * _NonDim->HeightScale() / _NonDim->ExtentScale() : 1.0;

	int ix, jy, kz;
	for(Zface zf = ff.GetGrid2().First<Zface>(); zf != ff.GetGrid2().Last<Zface>(); ++zf)
	{
		ix = ff.GetGrid2().ix(zf);
		jy = ff.GetGrid2().iy(zf);
		hs[nx*jy+ix] = (ff.HH(zf) + ff.SS(zf)) * H0;
		ss[nx*jy+ix] = (ff.SS(zf)) * H0;
		hh[nx*jy+ix] = (ff.HH(zf)) * H0;
	}
	for(Cell cl = ff.GetGrid3().First<Cell>(); cl != ff.GetGrid3().Last<Cell>(); ++cl)
	{
		ix = ff.GetGrid3().ix(cl);
		jy = ff.GetGrid3().iy(cl);
		kz = ff.GetGrid3().iz(cl);
		const C_Grid3D &G3 = ff.GetGrid3();
		cc[(ny*nx)*kz + nx*jy +ix] = ff.CC(cl);
		//  The velocity components live on the faces; the point is the cell
		//  center, so each is the mean of the two faces of the cell, or the
		//  one face there is in the outermost ghost layer.
		ux[(ny*nx)*kz + nx*jy +ix] = (ix+1 < nx ? 0.5*(ff.Ux(G3.HalfXm(cl)) + ff.Ux(G3.HalfXp(cl)))
		                                         : ff.Ux(G3.HalfXm(cl))) * U0;
		uy[(ny*nx)*kz + nx*jy +ix] = (ny == 1 ? ff.Uy(G3.HalfYm(cl))
		                              : jy+1 < ny ? 0.5*(ff.Uy(G3.HalfYm(cl)) + ff.Uy(G3.HalfYp(cl)))
		                                          : ff.Uy(G3.HalfYm(cl))) * U0;
		uz[(ny*nx)*kz + nx*jy +ix] = (kz+1 < nz ? 0.5*(ff.Uz(G3.HalfZm(cl)) + ff.Uz(G3.HalfZp(cl)))
		                                         : ff.Uz(G3.HalfZm(cl))) * W0;
	}

	fout.write((char*)&size_s, sizeof(int));	//First writes the length of the data that follows, in 4 bytes (sizeof(int))
	fout.write((char*)cc.data(), static_cast<std::streamsize>(size_s));	//Writes the data

	// Velocity vector
	fout.write((char*)&size_v, sizeof(int));
	for(kz=0; kz < nz ; kz++){
		for(jy=0; jy < ny ; jy++){
			for(ix = 0; ix < nx; ix++){
				fout.write((char*)&ux[(ny*nx)*kz + nx*jy +ix], static_cast<std::streamsize>(d_size));
				fout.write((char*)&uy[(ny*nx)*kz + nx*jy +ix], static_cast<std::streamsize>(d_size));
				fout.write((char*)&uz[(ny*nx)*kz + nx*jy +ix], static_cast<std::streamsize>(d_size));
			}
		}
	}

	// Coordinates
	//
	//  The file carries every cell, ghost layers included, so each layer is
	//  placed at the grid's own zeta for that cell index -- below 0 and above 1
	//  for the ghosts, which lie outside the film. kz/nz over all nz layers
	//  would put the real cells at the wrong heights and the top of the drawn
	//  film at (nz-1)/nz of the true free surface: 93 % of it for
	//  concz_1d_reg in ParaView, and the same 93 % of the thickness lost on
	//  every restart that reads the coordinates back.
	double zeta;
	fout.write((char*)&size_v, sizeof(int));
	for(kz=0; kz < nz ; kz++){
		zeta = ff.GetGrid3().ZOf<Cell>(0, 0, kz);
		for(jy=0; jy < ny ; jy++){
			for(ix = 0; ix < nx; ix++){
				xx = ff.GetGrid2().XOf<Zface>(ix, jy) * L0;
				yy = ff.GetGrid2().YOf<Zface>(ix, jy) * L0;
				fout.write((char*)&xx, static_cast<std::streamsize>(d_size));
				fout.write((char*)&yy, static_cast<std::streamsize>(d_size));
				zz = ss[nx*jy+ix] + (hs[nx*jy+ix] -ss[nx*jy+ix]) * zeta;
				fout.write((char*)&zz, static_cast<std::streamsize>(d_size));
			}
		}
	}

	// Film thickness, then substrate, one value per point
	fout.write((char*)&size_s, sizeof(int));
	for(kz=0; kz < nz ; kz++)
		for(jy=0; jy < ny ; jy++)
			for(ix = 0; ix < nx; ix++)
				fout.write((char*)&hh[nx*jy+ix], static_cast<std::streamsize>(d_size));
	fout.write((char*)&size_s, sizeof(int));
	for(kz=0; kz < nz ; kz++)
		for(jy=0; jy < ny ; jy++)
			for(ix = 0; ix < nx; ix++)
				fout.write((char*)&ss[nx*jy+ix], static_cast<std::streamsize>(d_size));

	fout<<endl;
	fout<<"</AppendedData>"<<endl;
	fout<<"</VTKFile>"<<endl;

	if(!Commit(fout, ftmp, fpath)) return false;
	_fnum++;
	return true;
}


//----------------------------
//   OutputNeutral
//----------------------------
void OutputNeutral::Build(C_TimeSteps &tt, NonDimensionalNumbers &nd, bool dimensioned, const FilePath &resultDir, const std::string &fhead, int fst)
{
	C_Output_Interface::Build(tt, nd, dimensioned, resultDir, fhead, ".dat", fst);

	_constructed = true;
}

std::string OutputNeutral::FileName() const
{
	char wfile_head[2000];

	sprintf_s(wfile_head, "%s%05d", _head.c_str(), _fnum);
	FilePath fpath = _dir / wfile_head;
	fpath += _ext;
	return fpath.c_str();
}

bool OutputNeutral::Write(Field_WM1D_Out &ff)
{
	if(!_constructed) return false;

	const std::string fpath = FileName(), ftmp = fpath + ".tmp";
	errno = 0;
	ofstream fout(ftmp.c_str());
	if(!fout.good()) return WriteFailed(fpath, System::LastIOError());

	//----------------------------------------------------------------------
	//  Write every digit.
	//
	//  This format is also the input of the ASCII restart. At the default 6
	//  significant digits the pressure is destroyed once the film thickness is
	//  put through the second difference (dividing by dx^2). On a grid fine
	//  enough to resolve the contact line (dx = 0.3 um), the output precision
	//  alone threw the pressure off by several percent, and the velocity that
	//  comes out of it by three orders of magnitude.
	//----------------------------------------------------------------------
	fout<<std::setprecision(17);

	int ix;
	//  Units in the header, so the file says what it holds. The scales are
	//  applied by Field_WM1D_Out::ConvertFromStg (model/LubFields.cpp).
	//  The driving potential (WritePotential) last, so that the columns
	//  before it keep their places
	const bool phi = ff.HasPotential();
	if(_dimensioned)
		fout<<"xx[mm] \t HS[um] \t Ux[m/s] \t CC[-] \t HH[um] \t SS[um]"<<(phi ? " \t PHI[Pa]" : "")<<endl;
	else
		fout<<"xx[-] \t HS[-] \t Ux[-] \t CC[-] \t HH[-] \t SS[-]"<<(phi ? " \t PHI[-]" : "")<<endl;
	for(ix = 0; ix < ff.nx(); ix++)
	{
		fout<<scientific
			<<ff.xx(ix)
			<<"\t" << ff.HH(ix) + ff.SS(ix)
			<<"\t" << ff.Ux(ix)
			<<"\t" << ff.CC(ix)
		    	<<"\t" << ff.HH(ix)
		   	<<"\t" << ff.SS(ix);
		if(phi) fout<<"\t" << ff.PHI(ix);
		fout<<endl;
	}
	if(!Commit(fout, ftmp, fpath)) return false;
	_fnum++;
	return true;
}

bool OutputNeutral::Write(Field_WM2D_Out &ff)
{
	if(!_constructed) return false;

	const std::string fpath = FileName(), ftmp = fpath + ".tmp";
	errno = 0;
	ofstream fout(ftmp.c_str());
	if(!fout.good()) return WriteFailed(fpath, System::LastIOError());

	int ix, jy;
	//----------------------------------------------------------------------
	//  Write every digit.
	//
	//  This format is also the input of the ASCII restart. At the default 6
	//  significant digits the pressure is destroyed once the film thickness is
	//  put through the second difference (dividing by dx^2). On a grid fine
	//  enough to resolve the contact line (dx = 0.3 um), the output precision
	//  alone threw the pressure off by several percent, and the velocity that
	//  comes out of it by three orders of magnitude.
	//----------------------------------------------------------------------
	fout<<std::setprecision(17);

	const bool phi = ff.HasPotential();
	if(_dimensioned)
		fout<<"xx[mm] \t yy[mm] \t HS[um] \t Ux[m/s] \t Uy[m/s] \t CC[-] \t HH[um] \t SS[um]"<<(phi ? " \t PHI[Pa]" : "")<<endl;
	else
		fout<<"xx[-] \t yy[-] \t HS[-] \t Ux[-] \t Uy[-] \t CC[-] \t HH[-] \t SS[-]"<<(phi ? " \t PHI[-]" : "")<<endl;
	for(jy = 0; jy < ff.ny(); jy++)
	{
		for(ix = 0; ix < ff.nx(); ix++)
		{
			fout<<scientific
				<<ff.xx(ix) <<"\t"<<ff.yy(jy)
				<<"\t"<< ff.HH(ix,jy) + ff.SS(ix,jy)
				<<"\t"<< ff.Ux(ix,jy) << "\t" << ff.Uy(ix,jy)
				<<"\t" << ff.CC(ix,jy) << "\t" << ff.HH(ix,jy)  <<"\t"<< ff.SS(ix,jy);
			if(phi) fout<<"\t" << ff.PHI(ix,jy);
			fout<<endl;
		}
	}
	if(!Commit(fout, ftmp, fpath)) return false;
	_fnum++;
	return true;
}
bool OutputNeutral::Write(Field_CZ &ff)
{
	if(!_constructed) return false;

	const std::string fpath = FileName(), ftmp = fpath + ".tmp";
	errno = 0;
	ofstream fout(ftmp.c_str());
	if(!fout.good()) return WriteFailed(fpath, System::LastIOError());
	C_Grid2D *G2 = ff.GetGridPtr2();
	C_Grid3D *G3 = ff.GetGridPtr3();

	//----------------------------------------------------------------------
	//  Write every digit.
	//
	//  This format is also the input of the ASCII restart. At the default 6
	//  significant digits the pressure is destroyed once the film thickness is
	//  put through the second difference (dividing by dx^2). On a grid fine
	//  enough to resolve the contact line (dx = 0.3 um), the output precision
	//  alone threw the pressure off by several percent, and the velocity that
	//  comes out of it by three orders of magnitude.
	//----------------------------------------------------------------------
	fout<<std::setprecision(17);

	//  Same convention as every other field writer: in-plane lengths in mm,
	//  film thickness in um, velocities in m/s, following Output.Dimensioned
	//  like every other file. It is also the ASCII restart input for ConcZ
	//  (C_Restart::ReadAscii_CZ, which walks the cells in the order of the
	//  loop below), so the columns and their order are part of that contract.
	double L0 = 1, H0 = 1, U0 = 1;
	if(_dimensioned)
	{
		L0 = _NonDim->ExtentScale() * 1e3;
		H0 = _NonDim->HeightScale() * 1e6;
		U0 = _NonDim->VelocityScale();
	}
	//  Uz is made dimensionless with H0/T0, not L0/T0 (see the .vts writer).
	const double W0 = _dimensioned ? U0 * _NonDim->HeightScale() / _NonDim->ExtentScale() : 1.0;

	//  A header line, so the columns say what they are and in which unit.
	//  The well-mixed files have always had one; this file had none.
	fout<<"# "<<(_dimensioned ? "xx[mm]" : "xx[-]");
	if(G3->ny() != 1) fout<<"\t"<<(_dimensioned ? "yy[mm]" : "yy[-]");
	fout<<"\t"<<(_dimensioned ? "zz[um]" : "zz[-]")
	    <<"\tCC[-]"
	    <<"\t"<<(_dimensioned ? "Ux[m/s]" : "Ux[-]");
	if(G3->ny() != 1) fout<<"\t"<<(_dimensioned ? "Uy[m/s]" : "Uy[-]");
	fout<<"\t"<<(_dimensioned ? "Uz[m/s]" : "Uz[-]")
	    <<"\t"<<(_dimensioned ? "HS[um]" : "HS[-]")
	    <<"\t"<<(_dimensioned ? "HH[um]" : "HH[-]")
	    <<"\t"<<(_dimensioned ? "SS[um]" : "SS[-]")
	    <<endl;

	//  The points written, axis by axis, in the layout of the well-mixed files
	//  (Field_WM2D_Out): the wall face, the cell centers, the wall face -- and,
	//  through the thickness, the substrate, the cell centers, the free surface.
	//  A value on a face is the mean of the two cells either side of it (a
	//  ghost and the first real cell), so the top row is the concentration at
	//  the free surface. Along y of a 1-D run there is the one cell, jy = 0.
	const std::vector<cz_ascii::Pt> px = cz_ascii::Points(G3->IminReal(), G3->ImaxReal(),
		[&](int i){ return G3->XOf<Cell>(i); }, G3->XminReal(), G3->XmaxReal());
	const std::vector<cz_ascii::Pt> py = (G3->ny() == 1)
		? std::vector<cz_ascii::Pt>(1, cz_ascii::Pt{0, 0, 0.0})
		: cz_ascii::Points(G3->JminReal(), G3->JmaxReal(),
			[&](int j){ return G3->YOf<Cell>(0, j); }, G3->YminReal(), G3->YmaxReal());
	const std::vector<cz_ascii::Pt> pz = cz_ascii::Points(G3->KminReal(), G3->KmaxReal(),
		[&](int k){ return G3->ZOf<Cell>(0, 0, k); }, G3->ZminReal(), G3->ZmaxReal());

	for(const cz_ascii::Pt &Z : pz)
	for(const cz_ascii::Pt &Y : py)
	for(const cz_ascii::Pt &X : px)
	{
		const int xs[2] = {X.a, X.b}, ys[2] = {Y.a, Y.b}, zs[2] = {Z.a, Z.b};
		double hh = 0, ss = 0, cc = 0, uu = 0, vv = 0, ww = 0;
		for(int a = 0; a < 2; a++) for(int b = 0; b < 2; b++)
		{
			const Zface c2 = G2->At<Zface>(xs[a], ys[b]);
			hh += 0.25 * ff.HH(c2);
			ss += 0.25 * ff.SS(c2);
			for(int c = 0; c < 2; c++)
				cc += 0.125 * ff.CC(G3->At<Cell>(xs[a], ys[b], zs[c]));
		}
		//  A velocity component on a face is the one stored there; at a cell
		//  center, the mean of the two faces of that cell. Across the other
		//  two axes it is averaged like a scalar.
		for(int b = 0; b < 2; b++) for(int c = 0; c < 2; c++)
		{
			const Cell cl = G3->At<Cell>(X.a, ys[b], zs[c]);
			uu += 0.25 * (X.a == X.b ? 0.5*(ff.Ux(G3->HalfXm(cl)) + ff.Ux(G3->HalfXp(cl)))
			                          : ff.Ux(G3->HalfXp(cl)));
		}
		if(G3->ny() != 1)
			for(int a = 0; a < 2; a++) for(int c = 0; c < 2; c++)
			{
				const Cell cl = G3->At<Cell>(xs[a], Y.a, zs[c]);
				vv += 0.25 * (Y.a == Y.b ? 0.5*(ff.Uy(G3->HalfYm(cl)) + ff.Uy(G3->HalfYp(cl)))
				                          : ff.Uy(G3->HalfYp(cl)));
			}
		for(int a = 0; a < 2; a++) for(int b = 0; b < 2; b++)
		{
			const Cell cl = G3->At<Cell>(xs[a], ys[b], Z.a);
			ww += 0.25 * (Z.a == Z.b ? 0.5*(ff.Uz(G3->HalfZm(cl)) + ff.Uz(G3->HalfZp(cl)))
			                          : ff.Uz(G3->HalfZp(cl)));
		}

		//  The height of the point, which is what the matching .vts carries in
		//  its point coordinates. Writing zeta here instead would give a
		//  column that runs to the height scale everywhere, including where
		//  the film is thinner than that, under a header saying [um].
		const double zz = ss + hh * Z.x;

		//  1-D: 8 columns, x z c u w h+s h s.
		//  2-D: 10 columns, x y z c u v w h+s h s.
		fout<<scientific<<X.x * L0 <<"\t";
		if(G3->ny() != 1) fout<<Y.x * L0 <<"\t";
		fout<< zz * H0 <<"\t"<< cc <<"\t"<< uu * U0 <<"\t";
		if(G3->ny() != 1) fout<< vv * U0 <<"\t";
		fout<< ww * W0
			<<"\t" << (hh+ss) * H0 <<"\t"<< hh * H0 << "\t" << ss * H0
			<<endl;
	}
	if(!Commit(fout, ftmp, fpath)) return false;
	_fnum++;
	return true;
}



//----------------------------
//   OutputTimeline
//----------------------------

void OutputTimeline::Build(C_TimeSteps &tt, NonDimensionalNumbers &nd, bool dimensioned,
	const FilePath &resultDir, const std::string &fhead, int num, std::vector<std::string> &labels, std::string pre_str)
{
	C_Output_Interface::Build(tt, nd, dimensioned, resultDir, fhead, ".dat");

	_num_vals = num;
	_labels.resize(_num_vals);
	_labels = labels;

	_constructed = true;

	_fpath =  resultDir / _head;
	_fpath += _ext;

	_previous = ReadLines(_fpath.c_str());

	errno = 0;
	ofstream fout(_fpath.c_str());
	if(!fout.good())
	{
		WriteFailed(_fpath, System::LastIOError());
		StopOnOutputError();
		return;
	}

	if(!pre_str.empty())
		fout<<pre_str<<nl;

	fout<<"Time ";
	if(_dimensioned) fout<<"[s]";
	for(int ii = 0; ii < _num_vals; ii++)
	{
		fout<<"\t"<<_labels[ii];
	}
	fout<<nl;
	fout.close();
	if(fout.fail())
	{
		WriteFailed(_fpath, System::LastIOError());
		StopOnOutputError();
	}
}
bool OutputTimeline::Resume(double tmax)
{
	if(!_constructed) return false;
	errno = 0;
	ofstream fout(_fpath.c_str(), ios::app);
	double t;
	for(size_t k = 0; k < _previous.size(); k++)
		if(RowTime(_previous[k], t) && t <= tmax) fout<<_previous[k]<<nl;
	_previous.clear();
	fout.close();
	if(fout.fail()) return WriteFailed(_fpath, System::LastIOError());
	return true;
}

//  One row, appended. The file is opened and closed for every row, so that
//  what is on disk is up to date while the run goes on; a run killed in the
//  middle of a row leaves that row without its newline, and a restart into
//  the same directory drops it (ReadLines).
bool OutputTimeline::Write(const ostringstream& oss)
{
	errno = 0;
	ofstream fout(_fpath.c_str(), ios::app);
	fout<<oss.str() << nl;
	fout.close();
	if(fout.fail()) return WriteFailed(_fpath, System::LastIOError());
	return true;
}
