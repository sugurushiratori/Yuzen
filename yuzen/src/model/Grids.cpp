/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Grids.cpp : construction of the computational grids
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <fstream>

#include "Grids.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif


const int C_Grid1D::_n_stencil = 5;
const int C_Grid2D::_n_stencil = 13;

////////////////////////////////////////////
C_Grid1D::C_Grid1D(const RectGrid<1> &g) : _num_ghostX(2)
{
	double xmin, xmax, dx;
	dx = g.dx();

	_Xmin_real = g.xmin();
	_Xmax_real = g.xmax();

	xmin = g.xmin() - _num_ghostX*dx;
	xmax = g.xmax() + _num_ghostX*dx;
	int nx = g.nx() + _num_ghostX*2;
	RectGrid<1>::Build(xmin, xmax, nx);
}
void C_Grid1D::Build(const Input::I_Grid &g)
{
	double xmin, xmax, dx;
	dx = (g.Xmax - g.Xmin) / (double)g.Nx;

	_Xmin_real = g.Xmin;
	_Xmax_real = g.Xmax;

	xmin = g.Xmin - _num_ghostX*dx;
	xmax = g.Xmax + _num_ghostX*dx;
	int nx = g.Nx + _num_ghostX*2;
	RectGrid<1>::Build(xmin, xmax, nx);
}

void C_Grid1D::operator=(const C_Grid1D &g)
{
	_Xmin_real = g._Xmin_real;
	_Xmax_real = g._Xmax_real;
	_num_ghostX = g._num_ghostX;
	RectGrid<1>::operator=(g);
}
////////////////////////////////////////////
C_Grid2D::C_Grid2D(const RectGrid<2> &g) : _num_ghostX(2), _num_ghostY(2)
{
	double xmin, xmax, dx;
	double ymin, ymax, dy;
	dx = g.dx();
	dy = g.dy();

	_Xmin_real = g.xmin();
	_Xmax_real = g.xmax();
	_Ymin_real = g.ymin();
	_Ymax_real = g.ymax();

	xmin = _Xmin_real - _num_ghostX*dx;
	xmax = _Xmax_real + _num_ghostX*dx;
	ymin = _Ymin_real - _num_ghostY*dy;
	ymax = _Ymax_real + _num_ghostY*dy;
	int nx = g.nx() + _num_ghostX*2;
	int ny = g.ny() + _num_ghostY*2;
	RectGrid<2>::Build(xmin, xmax, nx, ymin, ymax, ny);
}
void C_Grid2D::operator=(const C_Grid2D &g)
{
	_Xmin_real = g._Xmin_real;
	_Xmax_real = g._Xmax_real;
	_Ymin_real = g._Ymin_real;
	_Ymax_real = g._Ymax_real;
	_num_ghostX = g._num_ghostX;
	_num_ghostY = g._num_ghostY;
	RectGrid<2>::operator=(g);
}
void C_Grid2D::Build(const Input::I_Grid &g)
{
	if(g.Ny == 1) _num_ghostY = 0;

	// Build the grid with the ghost cells added on.
	_Xmin_real = g.Xmin;
	_Xmax_real = g.Xmax;
	_Ymin_real = g.Ymin;
	_Ymax_real = g.Ymax;

	double xmin, xmax, dx;
	double ymin, ymax, dy;
	dx = (g.Xmax - g.Xmin) / (double)g.Nx;
	dy = (g.Ymax - g.Ymin) / (double)g.Ny;

	xmin = g.Xmin - _num_ghostX*dx;
	xmax = g.Xmax + _num_ghostX*dx;
	ymin = g.Ymin - _num_ghostY*dy;
	ymax = g.Ymax + _num_ghostY*dy;
	int nx = g.Nx + _num_ghostX*2;
	int ny = g.Ny + _num_ghostY*2;

	RectGrid<2>::Build(xmin, xmax, nx, ymin, ymax, ny);
}


////////////////////////////////////////////
C_Grid3D::C_Grid3D(const RectGrid<3> &g)   : _num_ghostX(2), _num_ghostY(2), _num_ghostZ(2)
{
	double xmin, xmax, dx;
	double ymin, ymax, dy;
	double zmin, zmax, dz;
	dx = g.dx();
	dy = g.dy();
	dz = g.dz();

	_Xmin_real = g.xmin();
	_Xmax_real = g.xmax();
	_Ymin_real = g.ymin();
	_Ymax_real = g.ymax();
	_Zmin_real = g.zmin();
	_Zmax_real = g.zmax();

	xmin = _Xmin_real - _num_ghostX*dx;
	xmax = _Xmax_real + _num_ghostX*dx;
	ymin = _Ymin_real - _num_ghostY*dy;
	ymax = _Ymax_real + _num_ghostY*dy;
	zmin = _Zmin_real - _num_ghostZ*dz;
	zmax = _Zmax_real + _num_ghostZ*dz;
	int nx = g.nx() + _num_ghostX*2;
	int ny = g.ny() + _num_ghostY*2;
	int nz = g.nz() + _num_ghostZ*2;
	RectGrid<3>::Build(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz);
}
void C_Grid3D::operator=(const C_Grid3D &g)
{
	_Xmin_real = g._Xmin_real;
	_Xmax_real = g._Xmax_real;
	_Ymin_real = g._Ymin_real;
	_Ymax_real = g._Ymax_real;
	_Zmin_real = g._Zmin_real;
	_Zmax_real = g._Zmax_real;
	_num_ghostX = g._num_ghostX;
	_num_ghostY = g._num_ghostY;
	_num_ghostZ = g._num_ghostZ;
	RectGrid<3>::operator=(g);
}
void C_Grid3D::Build(const Input::I_Grid &g)
{
	if(g.Ny == 1) _num_ghostY = 0;

	// Build the grid with the ghost cells added on.
	_Xmin_real = g.Xmin;
	_Xmax_real = g.Xmax;
	_Ymin_real = g.Ymin;
	_Ymax_real = g.Ymax;
	_Zmin_real =  0.0;
	_Zmax_real =  1.0;

	double xmin, xmax, dx;
	double ymin, ymax, dy;
	double zmin, zmax, dz;
	dx = (g.Xmax - g.Xmin) / (double)g.Nx;
	dy = (g.Ymax - g.Ymin) / (double)g.Ny;
	dz = (  1.0  -  0.0  ) / (double)g.Nz;
	xmin = g.Xmin - _num_ghostX*dx;
	xmax = g.Xmax + _num_ghostX*dx;
	ymin = g.Ymin - _num_ghostY*dy;
	ymax = g.Ymax + _num_ghostY*dy;
	zmin = 0.0    - _num_ghostZ*dz;
	zmax = 1.0    + _num_ghostZ*dz;
	int nx = g.Nx + _num_ghostX*2;
	int ny = g.Ny + _num_ghostY*2;
	int nz = g.Nz + _num_ghostZ*2;

	RectGrid<3>::Build(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz);
}
