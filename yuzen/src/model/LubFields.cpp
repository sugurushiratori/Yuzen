/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LubFields.cpp : implementation of the lubrication-approximation fields
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "Opt.h"
#include <iostream>
#include <fstream>

#include "LubFields.h"
#include "NonDimensionalNumbers.h"
#include <cmath>
#include "FieldVisitor.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif


void Velocity_WM1D::Build(C_Grid1D &g)
{
	_G = &g;
	_mixed.Build(g);
	_LaplacePressure.Build(g);
	_SolutoCapillarity.Build(g);
	_DisjoiningPressure.Build(g);
	_VerticalGravity.Build(g);
	_HorizontalGravity.Build(g);
	_Centrifugal.Build(g);

	_constructed = true;
}
void Velocity_WM1D::Build(const Velocity_WM1D& vv)
{
	_G = vv._G;
	_mixed.Build(*_G);
	_LaplacePressure.Build(*_G);
	_SolutoCapillarity.Build(*_G);
	_DisjoiningPressure.Build(*_G);
	_VerticalGravity.Build(*_G);
	_HorizontalGravity.Build(*_G);
	_Centrifugal.Build(*_G);

	_constructed = true;
	operator=(vv);
}
void Velocity_WM1D::operator=(const Velocity_WM1D& vv)
{
	if(!_constructed) Build(vv);

	_mixed = vv._mixed;
	_LaplacePressure = vv._LaplacePressure;
	_SolutoCapillarity = vv._SolutoCapillarity;
	_DisjoiningPressure = vv._DisjoiningPressure;
	_VerticalGravity = vv._VerticalGravity;
	_HorizontalGravity = vv._HorizontalGravity;
	_Centrifugal = vv._Centrifugal;
}
void Velocity_WM1D::Cal_MixedVelocity()
{
	Node nd;
	double ux;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(nd, ux) if(yuzen::opt::Par((_G->imax() - _G->imin() + 1)))
#endif
	for(int ii = _G->imin(); ii <= _G->imax(); ++ii)
	{
		nd = _G->Of<Node>(ii);

		//--------------------------------------------------------------
		//  Plain sum of the per-term velocities.
		//
		//  The governing equation is h_t = -div(Q); when Q is expanded inside
		//  that -div, only the gravity term comes out with a + sign and the
		//  rest with a -. The _Magnitude of each physical term is given
		//  following that convention, so the velocity field each term writes
		//  already has its physical direction. A plain sum is therefore the
		//  right way to combine them.
		//
		//  *** Do not write this as a chain of differences (Math::Diff and
		//  the like). *** Every further term added would flip the sign of all
		//  the terms after it, and no new physical term could be added. The
		//  mixed velocity is used to advect the concentration equation, so
		//  getting the sign wrong reverses the direction in which the solute
		//  is transported.
		//--------------------------------------------------------------
		ux = _LaplacePressure.x(nd)
		   + _SolutoCapillarity.x(nd)
		   + _VerticalGravity.x(nd)
		   + _HorizontalGravity.x(nd)
		   + _Centrifugal.x(nd)
		   + _DisjoiningPressure.x(nd);

		_mixed.x(nd) = ux;
	}
}

////////////////////////////////////
Field_WM1D::Field_WM1D()
{
}

void Field_WM1D::Build(const Input::I_Grid &g, const Input::I_Substrate &i_subs, double scale)
{
	_G.Build(g);
	_height.Build(_G);
	_substrate.Build(_G);
	_conc.Build(_G);
	_velocity.Build(_G);

	C_Substrate::MakeSubstrateProfile(_G, _substrate, i_subs, scale);

	_constructed = true;
}

void Field_WM1D::Build(const Field_WM1D& ff)
{
	_G = ff.GetGrid();
	_height.Build(_G);
	_substrate.Build(_G);
	_conc.Build(_G);
	_velocity.Build(_G);

	_constructed = true;
	operator=(ff);
}
void Field_WM1D::operator=(const Field_WM1D& ff)
{
	if(!_constructed) Build(ff);
	else
	{
		_height = ff.HH();
		_substrate = ff.SS();
		_conc = ff.CC();
		_velocity = ff.UU();
	}
}

void Field_WM1D::Update_H(Vector<double> &hh)
{
	Xedge xe;
	int irow;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(xe,irow) if(yuzen::opt::Par((_G.imax() - _G.imin() + 1)))
#endif
	for(int ii = _G.imin(); ii <= _G.imax(); ii++)
	{
		xe = _G.Of<Xedge>(ii);
		irow = _G.ix(xe);
		_height(xe) = hh[irow];
	}
}

void Field_WM1D::Update_C(Vector<double> &cc)
{
	Xedge xe;
	int irow;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(xe,irow) if(yuzen::opt::Par((_G.imax() - _G.imin() + 1)))
#endif
	for(int ii = _G.imin(); ii <= _G.imax(); ii++)
	{
		xe = _G.Of<Xedge>(ii);
		irow = _G.ix(xe);
		_conc(xe) = cc[irow];
	}
}
void Field_WM1D::Update_C_fromHr(Vector<double> &Hr)
{
	Xedge xe;
	double hh;
	int irow;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(xe,irow,hh) if(yuzen::opt::Par((_G.imax() - _G.imin() + 1)))
#endif
	for(int ii = _G.imin(); ii <= _G.imax(); ii++)
	{
		xe = _G.Of<Xedge>(ii);
		irow = _G.ix(xe);
		hh = _height(xe);
		if(Math::IsZero(hh)) hh = 1;
		_conc(xe) = Hr[irow] / hh;
	}
}
void Field_WM1D::Accept(FieldVisitor &v)
{
	v.Visit(*this);
}


//------------------------
//   Output-side field
//------------------------

Field_WM1D_Out::Field_WM1D_Out()
{
	_dimensioned = false;
}

void Field_WM1D_Out::Build(const Field_WM1D& ff, double L0, double H0, double U0, bool dimensioned)
{
	_G = ff.GetGrid();
	_dimensioned = dimensioned;
	_L0 = L0;
	_H0 = H0;
	_U0 = U0;

	_nx = _G.nx();
	if(_nx != 1) _nx -=2;

	_nz = 2;

	_hh.assign(_nx, 0.0);
	_ss.assign(_nx*_nz, 0.0);
	_cc.assign(_nx*_nz, 0.0);
	_hs.assign(_nx*_nz, 0.0);
	_xx.assign(_nx, 0.0);

	_ux.assign(_nx, 0.0);
	_ux_Ca.assign(_nx, 0.0);
	_ux_Ma.assign(_nx, 0.0);
	_ux_Bo.assign(_nx, 0.0);
	_ux_Pg.assign(_nx, 0.0);
	_ux_Sp.assign(_nx, 0.0);

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_nx))
#endif
	for(int ii = 0; ii< _nx; ii++)
	{
		_hh[ii] = 0;
		_ux[ii] = 0;
		_ux_Ca[ii] = 0;
		_ux_Ma[ii] = 0;
		_ux_Bo[ii] = 0;
		_ux_Pg[ii] = 0;
		_ux_Sp[ii] = 0;
	}
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_nx*_nz))
#endif
	for(int ii = 0; ii< _nx*_nz; ii++)
	{
		_ss[ii] = 0;
		_cc[ii] = 0;
		_hs[ii] = 0;
	}

	double xmin, xmax;
	double dx;
	dx = _G.dx();
	xmin = _G.XminReal();
	xmax = _G.XmaxReal();

	_xx[0]		= xmin;
	if(_nx != 1)
	{
		_xx[1]		= xmin + 0.5*dx;
		_xx[_nx-1]	= xmax;
		_xx[_nx-2]	= xmax - 0.5*dx;
		for(int ii = 2; ii<_nx-2; ii++) _xx[ii] = _xx[ii-1] + dx;
	}


	if(_dimensioned)
	{
		for(int ii = 0; ii<_nx; ii++) _xx[ii] *= _L0 *1.0e3;	//the in-plane scale is mm
	}


}
//  The driving potential at a site from its film, substrate, concentration
//  and the second difference of h + s (PotentialParams)
double yuzen::PotentialAt(double h, double hs, double c, double lap, const PotentialParams &pp)
{
	double phi = 0;
	if(pp.capillary) phi += (pp.nd ? pp.nd->SurfaceTension(c) : 1.0) * lap;
	if(pp.dp != 0 && h > 0)
	{
		const double r = pp.hstar / h;
		phi += pp.dp * (std::pow(r, pp.n) - std::pow(r, pp.m));
	}
	if(pp.bo != 0) phi -= pp.bo * hs;
	return phi;
}

void Field_WM1D_Out::ConvertPotential(Field_WM1D &ff, const PotentialParams &pp)
{
	//  The potential on the sites of the field, ghosts as their inner
	//  neighbor (a boundary row is then that value), then laid out as the
	//  other scalars
	ScalarField1D<> phi(ff.HH());
	const double dx = _G.dx(), dx2 = dx * dx;
	for(int ii = _G.imin(); ii <= _G.imax(); ++ii)
	{
		Xedge xe = _G.Of<Xedge>(ii);
		if(_G.CheckGhost(xe) != C_Grid1D::Inside) continue;
		const Xedge xp = _G.ShiftXp(xe), xm = _G.ShiftXm(xe);
		const double vo = ff.HH(xe) + ff.SS(xe);
		const double lap = ((ff.HH(xp) + ff.SS(xp)) - 2.0 * vo + (ff.HH(xm) + ff.SS(xm))) / dx2;
		phi(xe) = PotentialAt(ff.HH(xe), vo, ff.CC(xe), lap, pp);
	}
	for(int ii = _G.imin(); ii <= _G.imax(); ++ii)
	{
		Xedge xe = _G.Of<Xedge>(ii);
		const int g = _G.CheckGhost(xe);
		if(g == C_Grid1D::bXmin) phi(xe) = phi(_G.ShiftXp(xe));
		else if(g == C_Grid1D::bXmax) phi(xe) = phi(_G.ShiftXm(xe));
	}
	_phi.assign((size_t)_nx * _nz, 0.0);
	ConvertScalar(_phi, phi);
	if(_dimensioned) for(int ii = 0; ii < _nx; ii++) _phi[ii] *= pp.scale;
	for(int kz = 1; kz < _nz; kz++)			// the same at each point of the column
		for(int ii = 0; ii < _nx; ii++) _phi[(size_t)_nx * kz + ii] = _phi[ii];
	_hasPhi = true;
}

void Field_WM1D_Out::ConvertScalar(std::vector<double> &ff, ScalarField1D<> &ss)
{
	int ix;
	Xedge xe;
	for(int ii = _G.imin(); ii <= _G.imax(); ++ii)
	{
		xe = _G.Of<Xedge>(ii);
		ix = _G.ix(xe);
		switch(_G.CheckGhost(xe))
		{
		case C_Grid1D::Inside:
			ff[ix-1]	= ss(xe);
			break;
		case C_Grid1D::bXmin:
			ff[ix-1]	= 0.5*(ss(xe) +ss(_G.ShiftXp(xe)));
			break;
		case C_Grid1D::bXmax:
			ff[ix-1]	= 0.5*(ss(xe) +ss(_G.ShiftXm(xe)));
			break;
		default: break;
		}
	}

}
void Field_WM1D_Out::ConvertVelocity(std::vector<double> &ux, VectorField1D<> &vv)
{
	int ix;
	Xedge xe;
	for(int ii = _G.imin(); ii <= _G.imax(); ++ii)
	{
		xe = _G.Of<Xedge>(ii);
		ix = _G.ix(xe);
		switch(_G.CheckGhost(xe))
		{
		case C_Grid1D::Inside:
			ux[ix-1]	= 0.5*(vv.x(_G.HalfXp(xe)) +vv.x(_G.HalfXm(xe)));
			break;
		case C_Grid1D::bXmin:
			ux[ix-1]	= vv.x(_G.HalfXp(xe));
			break;
		case C_Grid1D::bXmax:
			ux[ix-1]	= vv.x(_G.HalfXm(xe));
			break;
		default: break;
		}
	}
}

void Field_WM1D_Out::ConvertFromStg(Field_WM1D &ff)
{
	assert(_G.nx() == _nx + _G.NumGhostX());

	ConvertScalar(_hh, ff.HH());
	ConvertScalar(_ss, ff.SS());
	ConvertScalar(_cc, ff.CC());

	ConvertVelocity(_ux, ff.UU().UU());
	ConvertVelocity(_ux_Ca, ff.UU().U_LaplacePressure());
	ConvertVelocity(_ux_Ma, ff.UU().U_SolutoCapillarity());
	ConvertVelocity(_ux_Bo, ff.UU().U_VerticalGravity());
	ConvertVelocity(_ux_Pg, ff.UU().U_HorizontalGravity());
	ConvertVelocity(_ux_Sp, ff.UU().U_Centrifugal());

	int ix;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_nx))
#endif
	for(ix = 0; ix < _nx; ix++)
	{
		_hs[ix] = _hh[ix] + _ss[ix];
		_hs[_nx +ix] = _hh[ix] + _ss[ix];
		_cc[_nx +ix] = _cc[ix];
		_ss[_nx +ix] = _ss[ix];
	}

	if(_dimensioned == false) return;

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_nx))
#endif
	for(int ii = 0; ii< _nx; ii++)
	{
		_hh[ii] *= _H0  * 1e6;
		_ux[ii] *= _U0;
		_ux_Ca[ii] *= _U0;
		_ux_Ma[ii] *= _U0;
		_ux_Bo[ii] *= _U0;
		_ux_Pg[ii] *= _U0;
		_ux_Sp[ii] *= _U0;
	}
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_nx*_nz))
#endif
	for(int ii = 0; ii< _nx*_nz; ii++)
	{
		_ss[ii] *= _H0 * 1e6;
		_hs[ii] *= _H0 * 1e6;
	}

}



/////////////////////////////////////////
void Velocity_WM2D::Build(C_Grid2D &g)
{
	_G = &g;
	_mixed.Build(g);
	_LaplacePressure.Build(g);
	_SolutoCapillarity.Build(g);
	_DisjoiningPressure.Build(g);
	_VerticalGravity.Build(g);
	_HorizontalGravity.Build(g);
	_Centrifugal.Build(g);

	_constructed = true;
}
void Velocity_WM2D::Build(const Velocity_WM2D& vv)
{
	_G = vv._G;
	_mixed.Build(*_G);
	_LaplacePressure.Build(*_G);
	_SolutoCapillarity.Build(*_G);
	_DisjoiningPressure.Build(*_G);
	_VerticalGravity.Build(*_G);
	_HorizontalGravity.Build(*_G);
	_Centrifugal.Build(*_G);

	_constructed = true;
	operator=(vv);
}
void Velocity_WM2D::operator=(const Velocity_WM2D& vv)
{
	if(!_constructed) Build(vv);
	else
	{
		_mixed = vv._mixed;
		_LaplacePressure = vv._LaplacePressure;
		_SolutoCapillarity = vv._SolutoCapillarity;
		_DisjoiningPressure = vv._DisjoiningPressure;
		_VerticalGravity = vv._VerticalGravity;
		_HorizontalGravity = vv._HorizontalGravity;
		_Centrifugal = vv._Centrifugal;
	}
}
void Velocity_WM2D::Cal_MixedVelocity()
{
	Xedge xe;
	Yedge ye;
	double ux, uy;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(xe,ye,ux,uy) if(yuzen::opt::Par((_G->imax() - _G->imin() + 1)))
#endif
	for(int ii = _G->imin(); ii <= _G->imax(); ++ii)
	{
		xe = _G->Of<Xedge>(ii);		ye = _G->Of<Yedge>(ii);

		//  Plain sum of the per-term velocities. See the note in the 1-D version for why.
		ux = _LaplacePressure.x(ye)
		   + _SolutoCapillarity.x(ye)
		   + _VerticalGravity.x(ye)
		   + _HorizontalGravity.x(ye)
		   + _Centrifugal.x(ye)
		   + _DisjoiningPressure.x(ye);

		uy = _LaplacePressure.y(xe)
		   + _SolutoCapillarity.y(xe)
		   + _VerticalGravity.y(xe)
		   + _HorizontalGravity.y(xe)
		   + _Centrifugal.y(xe)
		   + _DisjoiningPressure.y(xe);

		_mixed.x(ye) = ux;		_mixed.y(xe) = uy;
	}
}

////////////////////////////////////
Field_WM2D::Field_WM2D()
{
}

void Field_WM2D::Build(const Input::I_Grid &g, const Input::I_Substrate &i_subs, double scale)
{
	_G.Build(g);
	_height.Build(_G);
	_substrate.Build(_G);
	_conc.Build(_G);
	_velocity.Build(_G);

	C_Substrate::MakeSubstrateProfile(_G, _substrate, i_subs, scale);

	_constructed = true;
}
void Field_WM2D::Build(const Field_WM2D& ff)
{
	_G = ff.GetGrid();
	_height.Build(_G);
	_substrate.Build(_G);
	_conc.Build(_G);
	_velocity.Build(_G);

	_constructed = true;
	operator=(ff);
}
void Field_WM2D::operator=(const Field_WM2D& ff)
{
	if(!_constructed) Build(ff);
	else
	{
		_height = ff.HH();
		_substrate = ff.SS();
		_conc = ff.CC();
		_velocity = ff.UU();
	}
}

void Field_WM2D::Update_H(Vector<double> &hh)
{
	Zface zf;
	int ix, jy;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(zf,ix,jy) if(yuzen::opt::Par((_G.imax() - _G.imin() + 1)))
#endif
	for(int ii = _G.imin(); ii <= _G.imax(); ii++)
	{
		zf = _G.Of<Zface>(ii);
		ix = _G.ix(zf);
		jy = _G.iy(zf);
		_height(zf) = hh[jy*_G.nx() + ix];
	}
}

void Field_WM2D::Update_C(Vector<double> &cc)
{
	Zface zf;
	int ix, jy;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(zf,ix,jy) if(yuzen::opt::Par((_G.imax() - _G.imin() + 1)))
#endif
	for(int ii = _G.imin(); ii <= _G.imax(); ii++)
	{
		zf = _G.Of<Zface>(ii);
		ix = _G.ix(zf);
		jy = _G.iy(zf);
		_conc(zf) = cc[jy*_G.nx() + ix];
	}
}
void Field_WM2D::Update_C_fromHr(Vector<double> &Hr)
{
	Zface zf;
	int ix, jy;
	double hh;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(zf,ix,jy,hh) if(yuzen::opt::Par((_G.imax() - _G.imin() + 1)))
#endif
	for(int ii = _G.imin(); ii <= _G.imax(); ii++)
	{
		zf = _G.Of<Zface>(ii);
		ix = _G.ix(zf);
		jy = _G.iy(zf);
		hh = _height(zf);
		if(Math::IsZero(hh)) hh = 1;
		_conc(zf) = Hr[jy*_G.nx() + ix] / hh;
	}
}
void Field_WM2D::Accept(FieldVisitor &v)
{
	v.Visit(*this);
}


//------------------------
//   Output-side field
//------------------------

Field_WM2D_Out::Field_WM2D_Out()
{
	_dimensioned = false;
}

void Field_WM2D_Out::Build(const Field_WM2D& ff, double L0, double H0, double U0, bool dimensioned)
{
	_G = ff.GetGrid();
	_dimensioned = dimensioned;
	_L0 = L0;
	_H0 = H0;
	_U0 = U0;

	_nx = _G.nx();
	_ny = _G.ny();
	if(_nx != 1) _nx -=2;
	if(_ny != 1) _ny -=2;

	_nz = 2;

	_hh.assign(_nx*_ny, 0.0);
	_ss.assign(_nx*_ny*_nz, 0.0);
	_cc.assign(_nx*_ny*_nz, 0.0);
	_hs.assign(_nx*_ny*_nz, 0.0);
	_xx.assign(_nx, 0.0);
	_yy.assign(_ny, 0.0);

	_ux.assign(_nx*_ny, 0.0);
	_uy.assign(_nx*_ny, 0.0);
	_ux_Ca.assign(_nx*_ny, 0.0);
	_uy_Ca.assign(_nx*_ny, 0.0);
	_ux_Ma.assign(_nx*_ny, 0.0);
	_uy_Ma.assign(_nx*_ny, 0.0);
	_ux_Bo.assign(_nx*_ny, 0.0);
	_uy_Bo.assign(_nx*_ny, 0.0);
	_ux_Pg.assign(_nx*_ny, 0.0);
	_uy_Pg.assign(_nx*_ny, 0.0);
	_ux_Sp.assign(_nx*_ny, 0.0);
	_uy_Sp.assign(_nx*_ny, 0.0);

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_nx*_ny))
#endif
	for(int ii = 0; ii< _nx*_ny; ii++)
	{
		_hh[ii] = 0;
		_ux[ii] = 0;		_uy[ii] = 0;
		_ux_Ca[ii] = 0;		_uy_Ca[ii] = 0;
		_ux_Ma[ii] = 0;		_uy_Ma[ii] = 0;
		_ux_Bo[ii] = 0;		_uy_Bo[ii] = 0;
		_ux_Pg[ii] = 0;		_uy_Pg[ii] = 0;
		_ux_Sp[ii] = 0;		_uy_Sp[ii] = 0;
	}
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_nx*_ny*_nz))
#endif
	for(int ii = 0; ii< _nx*_ny*_nz; ii++)
	{
		_ss[ii] = 0;
		_cc[ii] = 0;
		_hs[ii] = 0;
	}

	double xmin, xmax, ymin, ymax;
	double dx, dy;
	dx = _G.dx();
	dy = _G.dy();
	xmin = _G.XminReal();
	xmax = _G.XmaxReal();
	ymin = _G.YminReal();
	ymax = _G.YmaxReal();

	_xx[0]		= xmin;
	if(_nx != 1)
	{
		_xx[1]		= xmin + 0.5*dx;
		_xx[_nx-1]	= xmax;
		_xx[_nx-2]	= xmax - 0.5*dx;
		for(int ii = 2; ii<_nx-2; ii++) _xx[ii] = _xx[ii-1] + dx;
	}

	_yy[0]		= ymin;
	if(_ny != 1)
	{
		_yy[1]		= ymin + 0.5*dy;
		_yy[_ny-1]	= ymax;
		_yy[_ny-2]	= ymax - 0.5*dy;
		for(int jj = 2; jj<_ny-2; jj++)	_yy[jj] = _yy[jj-1] + dy;
	}

	if(_dimensioned)
	{
		for(int ii = 0; ii<_nx; ii++) _xx[ii] *= _L0 *1.0e3;	//the in-plane scale is mm
		for(int jj = 0; jj<_ny; jj++) _yy[jj] *= _L0 *1.0e3;
	}


}
void Field_WM2D_Out::ConvertPotential(Field_WM2D &ff, const PotentialParams &pp)
{
	ScalarField2D<> phi(ff.HH());
	const double dx = _G.dx(), dy = _G.dy(), dx2 = dx * dx, dy2 = dy * dy;
	for(int ii = _G.imin(); ii <= _G.imax(); ++ii)
	{
		Zface zf = _G.Of<Zface>(ii);
		if(_G.CheckGhost(zf) != C_Grid2D::Inside) continue;
		auto v = [&](const Zface &z) { return ff.HH(z) + ff.SS(z); };
		const double vo = v(zf);
		const double lap = (v(_G.ShiftXp(zf)) - 2.0 * vo + v(_G.ShiftXm(zf))) / dx2
		                 + (v(_G.ShiftYp(zf)) - 2.0 * vo + v(_G.ShiftYm(zf))) / dy2;
		phi(zf) = PotentialAt(ff.HH(zf), vo, ff.CC(zf), lap, pp);
	}
	//  Ghosts as their inner neighbor: the edges first, then the corners
	//  from the edges
	for(int pass = 0; pass < 2; ++pass)
		for(int ii = _G.imin(); ii <= _G.imax(); ++ii)
		{
			Zface zf = _G.Of<Zface>(ii);
			switch(_G.CheckGhost(zf))
			{
			case C_Grid2D::bXmin: phi(zf) = phi(_G.ShiftXp(zf)); break;
			case C_Grid2D::bXmax: phi(zf) = phi(_G.ShiftXm(zf)); break;
			case C_Grid2D::bYmin: phi(zf) = phi(_G.ShiftYp(zf)); break;
			case C_Grid2D::bYmax: phi(zf) = phi(_G.ShiftYm(zf)); break;
			case C_Grid2D::CornerMM: if(pass) phi(zf) = phi(_G.ShiftXpYp(zf)); break;
			case C_Grid2D::CornerPM: if(pass) phi(zf) = phi(_G.ShiftXmYp(zf)); break;
			case C_Grid2D::CornerMP: if(pass) phi(zf) = phi(_G.ShiftXpYm(zf)); break;
			case C_Grid2D::CornerPP: if(pass) phi(zf) = phi(_G.ShiftXmYm(zf)); break;
			default: break;
			}
		}
	const size_t layer = (size_t)_nx * _ny;
	_phi.assign(layer * _nz, 0.0);
	ConvertScalar(_phi, phi);
	if(_dimensioned) for(size_t k = 0; k < layer; k++) _phi[k] *= pp.scale;
	for(int kz = 1; kz < _nz; kz++)
		for(size_t k = 0; k < layer; k++) _phi[layer * kz + k] = _phi[k];
	_hasPhi = true;
}

void Field_WM2D_Out::ConvertScalar(std::vector<double> &ff, ScalarField2D<> &ss)
{
	int ix, jy;
	Zface zf;
	for(int ii = _G.imin(); ii <= _G.imax(); ++ii)
	{
		zf = _G.Of<Zface>(ii);
		ix = _G.ix(zf); jy = _G.iy(zf);
		if(jy == 0) jy = 1;	////////
		switch(_G.CheckGhost(zf))
		{
		case C_Grid2D::Inside:
			ff[IL(ix-1,jy-1)]	= ss(zf);
			break;
		case C_Grid2D::bXmin:
			ff[IL(ix-1,jy-1)]	= 0.5*(ss(zf) +ss(_G.ShiftXp(zf)));
			break;
		case C_Grid2D::bXmax:
			ff[IL(ix-1,jy-1)]	= 0.5*(ss(zf) +ss(_G.ShiftXm(zf)));
			break;
		case C_Grid2D::bYmin:
			ff[IL(ix-1,jy-1)]	= 0.5*(ss(zf) +ss(_G.ShiftYp(zf)));
			break;
		case C_Grid2D::bYmax:
			ff[IL(ix-1,jy-1)]	= 0.5*(ss(zf) +ss(_G.ShiftYm(zf)));
			break;
		case C_Grid2D::CornerMM:
			ff[IL(ix-1,jy-1)]	= 0.25*(ss(zf) +ss(_G.ShiftXp(zf)) +ss(_G.ShiftYp(zf)) +ss(_G.ShiftXpYp(zf)));
			break;
		case C_Grid2D::CornerPM:
			ff[IL(ix-1,jy-1)]	= 0.25*(ss(zf) +ss(_G.ShiftXm(zf)) +ss(_G.ShiftYp(zf)) +ss(_G.ShiftXmYp(zf)));
			break;
		case C_Grid2D::CornerMP:
			ff[IL(ix-1,jy-1)]	= 0.25*(ss(zf) +ss(_G.ShiftXp(zf)) +ss(_G.ShiftYm(zf)) +ss(_G.ShiftXpYm(zf)));
			break;
		case C_Grid2D::CornerPP:
			ff[IL(ix-1,jy-1)]	= 0.25*(ss(zf) +ss(_G.ShiftXm(zf)) +ss(_G.ShiftYm(zf)) +ss(_G.ShiftXmYm(zf)));
			break;
		default: break;
		}
	}

}
void Field_WM2D_Out::ConvertVelocity(std::vector<double> &ux, std::vector<double> &uy, VectorField2D<> &vv)
{
	int ix, jy;
	Zface zf;
	for(int ii = _G.imin(); ii <= _G.imax(); ++ii)
	{
		zf = _G.Of<Zface>(ii);
		ix = _G.ix(zf); jy = _G.iy(zf);
		if(jy == 0) jy = 1;	////////
		switch(_G.CheckGhost(zf))
		{
		case C_Grid2D::Inside:
			ux[IL(ix-1,jy-1)]	= 0.5*(vv.x(_G.HalfXp(zf)) +vv.x(_G.HalfXm(zf)));
			uy[IL(ix-1,jy-1)]	= 0.5*(vv.y(_G.HalfYp(zf)) +vv.y(_G.HalfYm(zf)));
			break;
		case C_Grid2D::bXmin:
			ux[IL(ix-1,jy-1)]	= vv.x(_G.HalfXp(zf));
			uy[IL(ix-1,jy-1)]	= 0.25*(vv.y(_G.HalfYp(zf)) +vv.y(_G.HalfYm(zf)) +vv.y(_G.HalfYp(_G.ShiftXp(zf))) +vv.y(_G.HalfYm(_G.ShiftXp(zf))));
			break;
		case C_Grid2D::bXmax:
			ux[IL(ix-1,jy-1)]	= vv.x(_G.HalfXm(zf));
			uy[IL(ix-1,jy-1)]	= 0.25*(vv.y(_G.HalfYp(zf)) +vv.y(_G.HalfYm(zf)) +vv.y(_G.HalfYp(_G.ShiftXm(zf))) +vv.y(_G.HalfYm(_G.ShiftXm(zf))));
			break;
		case C_Grid2D::bYmin:
			ux[IL(ix-1,jy-1)]	= 0.25*(vv.x(_G.HalfXp(zf)) +vv.x(_G.HalfXm(zf)) +vv.x(_G.HalfXp(_G.ShiftYp(zf))) +vv.x(_G.HalfXm(_G.ShiftYp(zf))));
			uy[IL(ix-1,jy-1)]	= vv.y(_G.HalfYp(zf));
			break;
		case C_Grid2D::bYmax:
			ux[IL(ix-1,jy-1)]	= 0.25*(vv.x(_G.HalfXp(zf)) +vv.x(_G.HalfXm(zf)) +vv.x(_G.HalfXp(_G.ShiftYm(zf))) +vv.x(_G.HalfXm(_G.ShiftYm(zf))));
			uy[IL(ix-1,jy-1)]	= vv.y(_G.HalfYm(zf));
			break;
		case C_Grid2D::CornerMM:
			ux[IL(ix-1,jy-1)]	= 0.5*(vv.x(_G.HalfXp(zf)) +vv.x(_G.HalfXp(_G.ShiftYp(zf))));
			uy[IL(ix-1,jy-1)]	= 0.5*(vv.y(_G.HalfYp(zf)) +vv.y(_G.HalfYp(_G.ShiftXp(zf))));
			break;
		case C_Grid2D::CornerPM:
			ux[IL(ix-1,jy-1)]	= 0.5*(vv.x(_G.HalfXm(zf)) +vv.x(_G.HalfXm(_G.ShiftYp(zf))));
			uy[IL(ix-1,jy-1)]	= 0.5*(vv.y(_G.HalfYp(zf)) +vv.y(_G.HalfYp(_G.ShiftXm(zf))));
			break;
		case C_Grid2D::CornerMP:
			ux[IL(ix-1,jy-1)]	= 0.5*(vv.x(_G.HalfXp(zf)) +vv.x(_G.HalfXp(_G.ShiftYm(zf))));
			uy[IL(ix-1,jy-1)]	= 0.5*(vv.y(_G.HalfYm(zf)) +vv.y(_G.HalfYm(_G.ShiftXp(zf))));
			break;
		case C_Grid2D::CornerPP:
			ux[IL(ix-1,jy-1)]	= 0.5*(vv.x(_G.HalfXm(zf)) +vv.x(_G.HalfXm(_G.ShiftYm(zf))));
			uy[IL(ix-1,jy-1)]	= 0.5*(vv.y(_G.HalfYm(zf)) +vv.y(_G.HalfYm(_G.ShiftXm(zf))));
			break;
		default: break;
		}
	}
}

void Field_WM2D_Out::ConvertFromStg(Field_WM2D &ff)
{
	assert(_G.nx() == _nx + _G.NumGhostX() && _G.ny() == _ny + _G.NumGhostY());

	ConvertScalar(_hh, ff.HH());
	ConvertScalar(_ss, ff.SS());
	ConvertScalar(_cc, ff.CC());

	ConvertVelocity(_ux,_uy, ff.UU().UU());
	ConvertVelocity(_ux_Ca,_uy_Ca, ff.UU().U_LaplacePressure());
	ConvertVelocity(_ux_Ma,_uy_Ma, ff.UU().U_SolutoCapillarity());
	ConvertVelocity(_ux_Bo,_uy_Bo, ff.UU().U_VerticalGravity());
	ConvertVelocity(_ux_Pg,_uy_Pg, ff.UU().U_HorizontalGravity());
	ConvertVelocity(_ux_Sp,_uy_Sp, ff.UU().U_Centrifugal());

	int ix, jy;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(jy) if(yuzen::opt::Par(_nx*_ny))
#endif
	for(ix = 0; ix < _nx; ix++)
	{
		for(jy = 0; jy < _ny; jy++)
		{
			_hs[IL(ix, jy, 0)] = _hh[IL(ix,jy)] + _ss[IL(ix,jy)];
			_hs[IL(ix, jy, 1)] = _hh[IL(ix,jy)] + _ss[IL(ix,jy)];
			_cc[IL(ix, jy, 1)] = _cc[IL(ix,jy,0)];
			_ss[IL(ix, jy, 1)] = _ss[IL(ix,jy,0)];
		}
	}

	if(_dimensioned == false) return;

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_nx*_ny))
#endif
	for(int ii = 0; ii< _nx*_ny; ii++)
	{
		_hh[ii] *= _H0  * 1e6;
		_ux[ii] *= _U0;			_uy[ii] *= _U0;
		_ux_Ca[ii] *= _U0;		_uy_Ca[ii] *= _U0;
		_ux_Ma[ii] *= _U0;		_uy_Ma[ii] *= _U0;
		_ux_Bo[ii] *= _U0;		_uy_Bo[ii] *= _U0;
		_ux_Pg[ii] *= _U0;		_uy_Pg[ii] *= _U0;
		_ux_Sp[ii] *= _U0;		_uy_Sp[ii] *= _U0;
	}
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_nx*_ny*_nz))
#endif
	for(int ii = 0; ii< _nx*_ny*_nz; ii++)
	{
		_ss[ii] *= _H0 * 1e6;
		_hs[ii] *= _H0 * 1e6;
	}

}


//------------------------
//   Field for the ConcZ model
//------------------------

Field_CZ::Field_CZ()
{
}
void Field_CZ::Build(const Input::I_Grid &g, const Input::I_Substrate &i_subs, double scale)
{
	// Build a grid with Nz = 1 for H and S.
	_G2.Build(g);
	_G3.Build(g);
	_height.Build(_G2);
	_substrate.Build(_G2);
	_conc.Build(_G3);
	_velocity.Build(_G3);

	C_Substrate::MakeSubstrateProfile(_G2,_substrate, i_subs, scale);

	_constructed = true;
}

void Field_CZ::Build(const Field_CZ& ff)
{
	_G2 = ff.GetGrid2();
	_G3 = ff.GetGrid3();

	_height.Build(ff.HH().Grid());
	_substrate.Build(ff.SS().Grid());
	_conc.Build(ff.CC().Grid());
	_velocity.Build(ff.UU().Grid());

	_constructed = true;
	operator=(ff);
}

void Field_CZ::operator=(const Field_CZ& ff)
{
	if(!_constructed) Build(ff);
	else
	{
		_height = ff.HH();
		_substrate = ff.SS();
		_conc = ff.CC();
		_velocity = ff.UU();
	}
}

void Field_CZ::Update_H(Vector<double> &hh)
{
	int ix, jy;

	for(Zface zf = _G2.First<Zface>(); zf != _G2.Last<Zface>(); ++zf)
	{
		ix = _G2.ix(zf);
		jy = _G2.iy(zf);
		_height(zf) = hh[_G2.nx()*jy +ix];
	}
}

void Field_CZ::Update_C(Vector<double> &cc)
{
	int ix, jy, kz;

	for(Cell cl = _G3.First<Cell>(); cl != _G3.Last<Cell>(); ++cl)
	{
		ix = _G3.ix(cl);
		jy = _G3.iy(cl);
		kz = _G3.iz(cl);
		_conc(cl) = cc[_G3.nx()*_G3.ny()*kz + _G3.nx()*jy +ix];
	}
}

 void Field_CZ::Accept(FieldVisitor &v)
{
	v.Visit(*this);
}
