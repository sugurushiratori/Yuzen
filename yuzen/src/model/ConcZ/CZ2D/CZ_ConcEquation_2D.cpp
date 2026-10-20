/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CZ_ConcEquation_2D.cpp : implementation of the concentration equation of the ConcZ model
*
*  ConcEquation.cpp
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "Opt.h"
#include "TimeSteps.h"
#include "Profile.h"
#include "CZ_ConcEquation_2D.h"
#include <vector>
#include <cmath>
#include <algorithm>

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

using namespace yuzen::util;
using namespace Global;


//##################################
//           Conc Equation
//##################################

void CZ_ConcEquation_2D::SetFlagBC()
{
	for(Node nd = _G3->First<Node>(); nd != _G3->Last<Node>(); ++nd)
	{
		if(_G3->CheckGhost(nd) != C_Grid3D::Inside)
			this->SetAsBoundary(_G3->IL(nd));
	}
}
void CZ_ConcEquation_2D::Build(Input &inp, C_Fluid *prop, C_TimeSteps *ts, Field_CZ *ff, NonDimensionalNumbers *nonDim, BoundaryCondition2D *bc)
{
	//Set up the pointers.
	_timeStep = ts;
	_field = ff;
	_prop = prop;
	_boundary = bc;
	_NonDim = nonDim;

	//Set up everything to do with the grid.
	_G2 = ff->GetGridPtr2();
	_G3 = ff->GetGridPtr3();
	_Nx = _G3->nx();
	_Ny = _G3->ny();
	_Nz = _G3->nz();

	_dtdh.resize(_Nx*_Ny);
	_mms.Build(inp, *nonDim, _Ny != 1, true);

	_fdm.Build(C_Interpolation_base4::Central_2nd);

	_VecB = 0;
	//  The face interpolation of the advection (Numerics.Convection.Concentration),
	//  in x, y and zeta alike
	_scheme.Build(C_Interpolation_base4::FromKeyOr(inp.Numerics.Convection.Concentration,
	                                               C_Interpolation_base4::Upwind_3rd));

	//This function calls SetMatrix, so it must not run until every model has been built.
	C_ImplicitSolver::Build(ff->CC(), *ts, C_ImplicitSolver::CrankNicholson, false);

	//------------------------
	//  Everything to do with the linear solver
	//------------------------
	//  Set the preconditioner, the scaling and the convergence criterion
	//  together. What decides the shape of the linear system is the in-plane
	//  dimension of the film. ConcZ also has a grid in the film-thickness
	//  direction, but that has no bearing on the "auto" choice.
	this->_bcg.Configure(inp.Numerics.ConcEquation.LinearSolver, (_Ny != 1) ? 2 : 1);
	this->_bcg.SetConcRole(true);
}

//-----------------------------------------------------------------------
//  The concentration equation in conservative form on the moving zeta grid.
//
//  With z = s + zeta h, a cell (ix, jy, kz) holds the solute h c dx dy dzeta,
//  and the equation is
//
//      (h c)_t + d_x(Qx c) + d_y(Qy c) + d_zeta(W c)
//          = d_x(h Fx) + d_y(h Fy) + d_zeta(G),
//
//  Qx = h u and Qy = h v the volume fluxes through the x and y faces (per
//  unit zeta), W = h Omega the one through the zeta faces, Fx = Df D (c_x -
//  A c_zeta) and its y counterpart the diffusive flux at fixed z, A = (zeta
//  h_x + s_x)/h, and
//
//      G = Df D [ As c_zeta / h - h A (c_x - A c_zeta) - h Ay (c_y - Ay c_zeta) ]
//
//  minus the diffusive flux through a zeta face (the conormal one: what
//  crosses the surface zeta = const, not only its vertical part). The volume
//  fluxes satisfy the discrete geometric conservation law cell by cell,
//
//      dzeta h_t + dzeta (Qx_R - Qx_L)/dx + dzeta (Qy_N - Qy_S)/dy + W_T - W_B = 0,
//
//  because W is built from it, from W = 0 on the substrate up
//  (BuildVolumeFluxes): what goes through the top face of the top cell is
//  then what leaves the column through its free surface, the evaporation.
//  Subtracting c times that from the equation, a row reads
//
//      h c_t = - sum_f F_f (c_f - c_P) / delta_f + (diffusion),
//
//  the sum over the faces f of the cell with F_f its outward volume flux and
//  c_f the face value of the advection scheme: a uniform c is a solution
//  exactly, whatever the film does. At the free surface no solute crosses:
//  the total flux W c - G through the top face is zero, so that face
//  contributes + W_T c_P / dzeta (the solvent leaves, the solute stays) and
//  no diffusive flux; at the substrate the total flux is zero as well. With
//  h c_t = h_mid (c_n+1 - c_n)/dt and c_P at the midpoint, the change of h c
//  over a step is (h c)_n+1 - (h c)_n = h_mid dc + c_mid dh, and the sum of
//  the rows is the sum of the face fluxes: the solute of a closed film is
//  conserved to rounding (with a limited scheme, to the tolerance its face
//  values are settled to, SettleLimiter).
//
//  (The non-conservative form c_t + u c_x + Omega c_zeta = L c, with the
//  surface condition D c_zeta = Pe E h c, does not do: the conormal parts of
//  the diffusive flux, h A D c_x at a sloping surface and s_x D c_x on a
//  sloping substrate, go through the boundary, and a drying film that is not
//  flat loses 3 to 25 % of the solute its evaporated fraction should have
//  concentrated, whatever dt, Nx and Nz.)
//
//  The advection is the scheme of Numerics.Convection.Concentration on every
//  face. A linear scheme is implicit with its own weights, one linear
//  function of c for the two cells of a face, so that the flux they exchange
//  is one number. A limited one is frozen at the field and settled to the
//  midpoint of the step (AdvWeights, SettleLimiter).
//-----------------------------------------------------------------------
void CZ_ConcEquation_2D::BuildVolumeFluxes()
{
	const bool has_y = (_Ny != 1);
	const double dx = _G3->dx();
	const double dy = has_y ? _G3->dy() : 1.0;
	const double dzeta = _G3->dz();
	const int jlo = has_y ? 2      : 0;
	const int jhi = has_y ? _Ny-2  : 1;
	const size_t n3 = (size_t)(_Nx*_Ny*_Nz);
	if(_Qx.size() != n3) { _Qx.assign(n3, 0.0); _Qy.assign(n3, 0.0); _W.assign(n3, 0.0); }

	//  x faces ix = 2 .. _Nx-2 (between ix-1 and ix): the film at the face
	//  times the mean of the velocity profile over each cell (layerMean),
	//  and the difference between the flux of the film-height equation and
	//  the sum of the layers spread evenly through the depth (the two are
	//  integrals of the same profile, the one with the positivity preserving
	//  h factors and the depth means of 1/mu, the other from the profile at
	//  the cell centers; they differ by the truncation error in x).
	//
	//  The mean, not the value at the center: a cell exchanges what crosses
	//  its whole face. The flux the profile at the centers adds up to falls
	//  short of the film's by u'' dzeta^2/24 a cell, which, spread back over
	//  the layers, put an error of that order into the velocity of every
	//  layer, whatever the grid in x (with 4 cells through the film, the
	//  solute of tools/mms.py cz_conc_1d converged at 0.3 in dx).
	auto layerMean = [&](const double *u, double *um)
	{
		//  u[2 .. _Nz-2]: the profile at the cell centers, and above the
		//  surface its continuation. In the cell next to the substrate the
		//  mean of the parabola through u = 0 there and the two lowest
		//  centers, elsewhere u + u'' dzeta^2/24 with the centered second
		//  difference (exact for a cubic); the profile of a uniform
		//  viscosity is a parabola.
		const int k0 = 2, k1 = _Nz-3;
		for(int kz = k0; kz <= k1; kz++)
		{
			if(kz == k0) um[kz] = u[kz] + (u[kz+1] - 3.0*u[kz])/18.0;
			else         um[kz] = u[kz] + (u[kz+1] - 2.0*u[kz] + u[kz-1])/24.0;
		}
	};
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for
#endif
	for(int ix = 2; ix <= _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		std::vector<double> u(_Nz, 0.0), um(_Nz, 0.0);
		for(int kz = 2; kz <= _Nz-2; kz++) u[kz] = _field->Ux(_G3->At<Xface>(ix, jy, kz));
		layerMean(u.data(), um.data());
		const double hf = 0.5*(_field->HH(_G2->At<Zface>(ix-1, jy)) + _field->HH(_G2->At<Zface>(ix, jy)));
		double sum = 0.0;
		for(int kz = 2; kz < _Nz-2; kz++)
		{
			const double q = hf * um[kz];
			_Qx[IL(ix, jy, kz)] = q;
			sum += q * dzeta;
		}
		if(_qFilmX)
		{
			const double d = (*_qFilmX)[IL(ix, jy)] - sum;
			for(int kz = 2; kz < _Nz-2; kz++) _Qx[IL(ix, jy, kz)] += d;
		}
	}
	if(has_y)
	{
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for
#endif
		for(int ix = 2; ix < _Nx-2; ix++)
		for(int jy = 2; jy <= _Ny-2; jy++)
		{
			std::vector<double> u(_Nz, 0.0), um(_Nz, 0.0);
			for(int kz = 2; kz <= _Nz-2; kz++) u[kz] = _field->Uy(_G3->At<Yface>(ix, jy, kz));
			layerMean(u.data(), um.data());
			const double hf = 0.5*(_field->HH(_G2->At<Zface>(ix, jy-1)) + _field->HH(_G2->At<Zface>(ix, jy)));
			double sum = 0.0;
			for(int kz = 2; kz < _Nz-2; kz++)
			{
				const double q = hf * um[kz];
				_Qy[IL(ix, jy, kz)] = q;
				sum += q * dzeta;
			}
			if(_qFilmY)
			{
				const double d = (*_qFilmY)[IL(ix, jy)] - sum;
				for(int kz = 2; kz < _Nz-2; kz++) _Qy[IL(ix, jy, kz)] += d;
			}
		}
	}

	//  The zeta faces of every column, from the substrate up. _W at kz is the
	//  bottom face of cell kz; kz = _Nz-2 is the free surface.
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for
#endif
	for(int ix = 2; ix < _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		const double ht = _dtdh[IL(ix, jy)];
		//  The sideways outflow of 1-D (the centrifugal term), spread through
		//  the column as the velocity profile goes, adding up to what the
		//  film-height equation takes out
		const double side = _sideSink ? (*_sideSink)[IL(ix, jy)] : 0.0;
		double wsum = 0.0;
		if(side != 0.0 && _ppSide && !_ppSide->empty())
			for(int kz = 2; kz < _Nz-2; kz++)
				wsum += 0.5*((*_ppSide)[_G3->At<Xface>(ix, jy, kz).Index()]
				           + (*_ppSide)[_G3->At<Xface>(ix+1, jy, kz).Index()]) * dzeta;
		double w = 0.0;
		_W[IL(ix, jy, 2)] = 0.0;
		for(int kz = 2; kz < _Nz-2; kz++)
		{
			double div = ht + (_Qx[IL(ix+1, jy, kz)] - _Qx[IL(ix, jy, kz)])/dx;
			if(has_y) div += (_Qy[IL(ix, jy+1, kz)] - _Qy[IL(ix, jy, kz)])/dy;
			if(side != 0.0)
			{
				double o = 1.0;
				if(wsum > 0.0)
					o = 0.5*((*_ppSide)[_G3->At<Xface>(ix, jy, kz).Index()]
					       + (*_ppSide)[_G3->At<Xface>(ix+1, jy, kz).Index()]) / wsum;
				div += side * o;
			}
			w -= dzeta * div;
			_W[IL(ix, jy, kz+1)] = w;
		}
		//  The free surface: the recursion arrives at what leaves the column
		//  there (evaporation, less the manufactured source) only to the
		//  rounding of the sum of the layers' divergences, which can be large
		//  terms (a strong shear over a sloping substrate); taken from the
		//  film-height equation instead, it is that exactly -- zero without
		//  evaporation, so that a uniform concentration is not stirred by
		//  that rounding, which a stiff step amplifies (to 1e-12 in
		//  tools/equivalence.py cz_uniform). The top cell's volume balance
		//  then closes to that rounding, unamplified.
		if(_colSrc)
			_W[IL(ix, jy, _Nz-2)] = -((*_colSrc)[IL(ix, jy)] + side);
	}
}

double CZ_ConcEquation_2D::CourantRate() const
{
	const size_t n3 = (size_t)(_Nx*_Ny*_Nz);
	if(_Qx.size() != n3) return 0.0;
	const bool has_y = (_Ny != 1);
	const double dx = _G3->dx();
	const double dy = has_y ? _G3->dy() : 1.0;
	const double dzeta = _G3->dz();
	const int jlo = has_y ? 2      : 0;
	const int jhi = has_y ? _Ny-2  : 1;
	double rate = 0.0;
	for(int ix = 2; ix < _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		const double h = _field->HH(_G2->At<Zface>(ix, jy));
		if(!(h > 0.0)) continue;
		for(int kz = 2; kz < _Nz-2; kz++)
		{
			double r = std::max(std::fabs(_Qx[IL(ix, jy, kz)]), std::fabs(_Qx[IL(ix+1, jy, kz)])) / dx;
			if(has_y) r += std::max(std::fabs(_Qy[IL(ix, jy, kz)]), std::fabs(_Qy[IL(ix, jy+1, kz)])) / dy;
			r += std::max(std::fabs(_W[IL(ix, jy, kz)]), std::fabs(_W[IL(ix, jy, kz+1)])) / dzeta;
			rate = std::max(rate, r / h);
		}
	}
	return rate;
}

void CZ_ConcEquation_2D::SetMatrixA()
{
	_crsA.ClearData();
	_VecB = 0;
	BuildVolumeFluxes();
	if(_react.size() != (size_t)_crsA.RowSize()) _react.assign(_crsA.RowSize(), 0.0);

	const double dx = _G3->dx();
	const double dzeta = _G3->dz();
	const double As = _NonDim->As();
	const double Df = _NonDim->Df();

	const bool has_y = (_Ny != 1);
	const double dy   = has_y ? _G3->dy() : 1.0;
	const int jlo = has_y ? 2      : 0;
	const int jhi = has_y ? _Ny-2  : 1;
	const double tmid = _timeStep->NowTime() + 0.5*_timeStep->DT().now;

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for
#endif
	for(int ix = 2; ix < _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		auto C  = [&](int i, int j, int k) { return _field->CC(_G3->At<Cell>(i, j, k)); };
		auto Hh = [&](int i, int j) { return _field->HH(_G2->At<Zface>(i, j)); };
		auto Ss = [&](int i, int j) { return _field->SS(_G2->At<Zface>(i, j)); };

		//  The manufactured source of the column, and the conormal flux of
		//  the manufactured solution through its substrate and its surface
		std::vector<double> mmsSrc;
		double mmsGb = 0.0, mmsGt = 0.0, mmsCs = 0.0;
		const double xc = _G2->X(_G2->At<Zface>(ix, jy));
		const double yc = has_y ? _G2->Y(_G2->At<Zface>(ix, jy)) : 0.0;
		if(_mms.Enabled())
		{
			//  The mean of the source over each cell (three Gauss points), as
			//  the rows are cell balances
			const int nc = _Nz-4;
			const double g = 0.5*dzeta*std::sqrt(0.6);
			std::vector<double> zc(3*nc), sc(3*nc);
			for(int k = 0; k < nc; k++)
			{
				const double z0 = dzeta*(0.5 + (double)k);
				zc[3*k] = z0 - g;	zc[3*k+1] = z0;	zc[3*k+2] = z0 + g;
			}
			_mms.SCZColumnNd(xc, yc, tmid, zc.data(), 3*nc, sc.data());
			mmsSrc.resize(nc);
			for(int k = 0; k < nc; k++)
				mmsSrc[k] = (5.0*sc[3*k] + 8.0*sc[3*k+1] + 5.0*sc[3*k+2])/18.0;
			mmsGb = _mms.CZFluxGNd(xc, yc, 0.0, tmid);
			mmsGt = _mms.CZFluxGNd(xc, yc, 1.0, tmid);
			mmsCs = _mms.CZnd(xc, yc, 1.0, tmid);
		}

		//  Whether a face of this column is that of a wall (Symmetry)
		const bool wallXmin = (ix == 2) && BoundaryCondition::IsWall(_boundary->TypeXmin());
		const bool wallXmax = (ix == _Nx-3) && BoundaryCondition::IsWall(_boundary->TypeXmax());
		const bool wallYmin = has_y && (jy == 2) && BoundaryCondition::IsWall(_boundary->TypeYmin());
		const bool wallYmax = has_y && (jy == _Ny-3) && BoundaryCondition::IsWall(_boundary->TypeYmax());

		const double hoo = Hh(ix, jy);
		const double hpo = Hh(ix+1, jy), hPo = Hh(ix+2, jy), hmo = Hh(ix-1, jy), hMo = Hh(ix-2, jy);
		const double spo = Ss(ix+1, jy), sPo = Ss(ix+2, jy), smo = Ss(ix-1, jy), sMo = Ss(ix-2, jy);
		const double soo = Ss(ix, jy);
		const double hL = 0.5*(hoo+hmo), hR = 0.5*(hoo+hpo);
		//  The slopes of the mapping, at the column center and on the x faces
		const double hx0 = (hMo -8.0*hmo +8.0*hpo-hPo)/(12.0*dx);
		const double sx0 = (sMo -8.0*smo +8.0*spo-sPo)/(12.0*dx);
		const double hxL = (hMo-27.0*hmo+27.0*hoo-hpo)/(24.0*dx);
		const double sxL = (sMo-27.0*smo+27.0*soo-spo)/(24.0*dx);
		const double hxR = (hmo-27.0*hoo+27.0*hpo-hPo)/(24.0*dx);
		const double sxR = (smo-27.0*soo+27.0*spo-sPo)/(24.0*dx);
		double hy0 = 0, sy0 = 0, hyS = 0, syS = 0, hyN = 0, syN = 0, hS = 0, hN = 0;
		if(has_y)
		{
			const double hoP = Hh(ix, jy+2), hop = Hh(ix, jy+1), hom = Hh(ix, jy-1), hoM = Hh(ix, jy-2);
			const double soP = Ss(ix, jy+2), sop = Ss(ix, jy+1), som = Ss(ix, jy-1), soM = Ss(ix, jy-2);
			hS = 0.5*(hoo+hom);  hN = 0.5*(hoo+hop);
			hy0 = (hoM -8.0*hom +8.0*hop-hoP)/(12.0*dy);
			sy0 = (soM -8.0*som +8.0*sop-soP)/(12.0*dy);
			hyS = (hoM-27.0*hom+27.0*hoo-hop)/(24.0*dy);
			syS = (soM-27.0*som+27.0*soo-sop)/(24.0*dy);
			hyN = (hom-27.0*hoo+27.0*hop-hoP)/(24.0*dy);
			syN = (som-27.0*soo+27.0*sop-soP)/(24.0*dy);
		}

		for(int kz = 2; kz < _Nz-2; kz++)
		{
			const int irow = IL(ix, jy, kz);
			const double zeta = dzeta*(-1.5 +(double)kz);
			const bool bottom = (kz == 2), top = (kz == _Nz-3);

			//  The coefficients of h c_t = R(c): xz[di+2][dk+2] of the cell
			//  (ix+di, jy, kz+dk), yz[dj+2][dk+1] of (ix, jy+dj, kz+dk), dj != 0;
			//  and the part of R that does not depend on c.
			double xz[5][5] = {{0}}, yz[5][3] = {{0}};
			double src = 0.0, react = 0.0;
			auto X = [&](int di, int dk, double v) { xz[di+2][dk+2] += v; };
			auto Y = [&](int dj, int dk, double v) { if(dj == 0) xz[2][dk+2] += v; else yz[dj+2][dk+1] += v; };

			//---------------- advection ----------------
			//  Each face adds -F (c_f - c_P) to R, F its outward volume
			//  flux; the c_P parts are in the diagonal (below).
			//  x faces
			{
				const double QL = _Qx[IL(ix, jy, kz)], QR = _Qx[IL(ix+1, jy, kz)];
				const double fL[4] = {C(ix-2,jy,kz), C(ix-1,jy,kz), C(ix,jy,kz), C(ix+1,jy,kz)};
				const double fR[4] = {C(ix-1,jy,kz), C(ix,jy,kz), C(ix+1,jy,kz), C(ix+2,jy,kz)};
				double wL[4], wR[4], kL = 0, kR = 0;
				AdvWeights(fL, QL, wL, kL, C_Interpolation_base4::R);
				AdvWeights(fR, QR, wR, kR, C_Interpolation_base4::L);
				//  face L: inflow +QL (c_L - c_P), face R: outflow -QR (c_R - c_P)
				for(int t = 0; t < 4; t++)
				{
					X(t-2, 0, +QL*wL[t]/dx);
					X(t-1, 0, -QR*wR[t]/dx);
				}
				src += (QL*kL - QR*kR)/dx;
			}
			if(has_y)
			{
				const double QS = _Qy[IL(ix, jy, kz)], QN = _Qy[IL(ix, jy+1, kz)];
				const double fS[4] = {C(ix,jy-2,kz), C(ix,jy-1,kz), C(ix,jy,kz), C(ix,jy+1,kz)};
				const double fN[4] = {C(ix,jy-1,kz), C(ix,jy,kz), C(ix,jy+1,kz), C(ix,jy+2,kz)};
				double wS[4], wN[4], kS = 0, kN = 0;
				AdvWeights(fS, QS, wS, kS, C_Interpolation_base4::R);
				AdvWeights(fN, QN, wN, kN, C_Interpolation_base4::L);
				for(int t = 0; t < 4; t++)
				{
					Y(t-2, 0, +QS*wS[t]/dy);
					Y(t-1, 0, -QN*wN[t]/dy);
				}
				src += (QS*kS - QN*kN)/dy;
			}
			//  zeta faces
			{
				const double WB = _W[IL(ix, jy, kz)], WT = _W[IL(ix, jy, kz+1)];
				if(!bottom)
				{
					const double fB[4] = {C(ix,jy,kz-2), C(ix,jy,kz-1), C(ix,jy,kz), C(ix,jy,kz+1)};
					double wB[4], kB = 0;
					AdvWeights(fB, WB, wB, kB, C_Interpolation_base4::R);
					for(int t = 0; t < 4; t++) X(0, t-2, +WB*wB[t]/dzeta);
					src += WB*kB/dzeta;
				}
				if(!top)
				{
					const double fT[4] = {C(ix,jy,kz-1), C(ix,jy,kz), C(ix,jy,kz+1), C(ix,jy,kz+2)};
					double wT[4], kT = 0;
					AdvWeights(fT, WT, wT, kT, C_Interpolation_base4::L);
					for(int t = 0; t < 4; t++) X(0, t-1, -WT*wT[t]/dzeta);
					src -= WT*kT/dzeta;
				}
				else
				{
					//  The free surface: no solute crosses it, while W_T of
					//  solvent does. (With the manufactured solution, the
					//  solute flux W_T c_s - G of that solution crosses it.)
					react += WT/dzeta;
					if(_mms.Enabled()) src -= (WT*mmsCs - mmsGt)/dzeta;
				}
				//  The substrate: nothing crosses it (the manufactured
				//  solution's -G, which it carries diffusively).
				if(bottom && _mms.Enabled()) src -= mmsGb/dzeta;
			}

			//---------------- diffusion ----------------
			const double Coo = C(ix, jy, kz);
			//  x faces: + (h Fx)_R / dx - (h Fx)_L / dx
			{
				const double DL = Df * _NonDim->Diffusivity(0.5*(C(ix-1,jy,kz) + Coo));
				const double DR = Df * _NonDim->Diffusivity(0.5*(Coo + C(ix+1,jy,kz)));
				const double AL = (hxL*zeta + sxL) / hL;
				const double AR = (hxR*zeta + sxR) / hR;
				//  No solute through a wall face. The mirrored ghosts make the
				//  normal part of the diffusive flux zero there and the slope
				//  of the mapping (AL, AR) too, but the latter only to the
				//  rounding of its four-point difference, which put that much
				//  solute through the wall every step.
				const double aL = wallXmin ? 0.0 : hL*DL/dx, aR = wallXmax ? 0.0 : hR*DR/dx;
				//  face L: (c_P - c_W)/dx - AL (c_W,k+1 + c_P,k+1 - c_W,k-1 - c_P,k-1)/(4 dzeta)
				X(-1, 0, +aL/dx);
				X( 0, 0, -aL/dx);
				X(-1, 1, +aL*AL/(4.0*dzeta));	X(0, 1, +aL*AL/(4.0*dzeta));
				X(-1,-1, -aL*AL/(4.0*dzeta));	X(0,-1, -aL*AL/(4.0*dzeta));
				//  face R
				X( 1, 0, +aR/dx);
				X( 0, 0, -aR/dx);
				X( 0, 1, -aR*AR/(4.0*dzeta));	X(1, 1, -aR*AR/(4.0*dzeta));
				X( 0,-1, +aR*AR/(4.0*dzeta));	X(1,-1, +aR*AR/(4.0*dzeta));
			}
			if(has_y)
			{
				const double DS = Df * _NonDim->Diffusivity(0.5*(C(ix,jy-1,kz) + Coo));
				const double DN = Df * _NonDim->Diffusivity(0.5*(Coo + C(ix,jy+1,kz)));
				const double AS = (hyS*zeta + syS) / hS;
				const double AN = (hyN*zeta + syN) / hN;
				const double aS = wallYmin ? 0.0 : hS*DS/dy, aN = wallYmax ? 0.0 : hN*DN/dy;
				Y(-1, 0, +aS/dy);
				Y( 0, 0, -aS/dy);
				Y(-1, 1, +aS*AS/(4.0*dzeta));	Y(0, 1, +aS*AS/(4.0*dzeta));
				Y(-1,-1, -aS*AS/(4.0*dzeta));	Y(0,-1, -aS*AS/(4.0*dzeta));
				Y( 1, 0, +aN/dy);
				Y( 0, 0, -aN/dy);
				Y( 0, 1, -aN*AN/(4.0*dzeta));	Y(1, 1, -aN*AN/(4.0*dzeta));
				Y( 0,-1, +aN*AN/(4.0*dzeta));	Y(1,-1, +aN*AN/(4.0*dzeta));
			}
			//  zeta faces: + G_T / dzeta - G_B / dzeta, the interior ones only
			//  (through the substrate and the surface the total flux is set)
			for(int side = 0; side < 2; side++)
			{
				if(side == 0 && bottom) continue;
				if(side == 1 && top) continue;
				const int k0 = (side == 0) ? kz-1 : kz;		// the face lies between k0 and k0+1
				const double zf = zeta + (side == 0 ? -0.5 : 0.5)*dzeta;
				const double Dz = Df * _NonDim->Diffusivity(0.5*(C(ix,jy,k0) + C(ix,jy,k0+1)));
				const double hA  = hx0*zf + sx0;		// h A on the face
				const double hAy = hy0*zf + sy0;
				const double g = (side == 1) ? +1.0 : -1.0;	// + G_T, - G_B
				const double kzz = g * Dz * (As/hoo + (hA*hA + hAy*hAy)/hoo) / (dzeta*dzeta);
				X(0, k0+1-kz, +kzz);
				X(0, k0  -kz, -kzz);
				//  - h A c_x on the face, c_x from the columns either side at k0, k0+1
				const double kx = -g * Dz * hA / (4.0*dx) / dzeta;
				X( 1, k0-kz, +kx);	X( 1, k0+1-kz, +kx);
				X(-1, k0-kz, -kx);	X(-1, k0+1-kz, -kx);
				if(has_y)
				{
					const double ky = -g * Dz * hAy / (4.0*dy) / dzeta;
					Y( 1, k0-kz, +ky);	Y( 1, k0+1-kz, +ky);
					Y(-1, k0-kz, -ky);	Y(-1, k0+1-kz, -ky);
				}
			}

			//---------------- the row ----------------
			//  c_t = -A c - B with A = -R/h, B = -src/h - (manufactured source).
			//  The diagonal is minus the sum of the rest, plus the free-surface
			//  term: the row of a uniform c is that term alone, exactly.
			const double rh = 1.0/hoo;
			double offsum = 0.0;
			auto put = [&](int col, double v) { _crsA(irow, col) = -v*rh; offsum += -v*rh; };
			put(IL(ix  , jy, kz-2), xz[2][0]);
			put(IL(ix-1, jy, kz-1), xz[1][1]);
			put(IL(ix  , jy, kz-1), xz[2][1]);
			put(IL(ix+1, jy, kz-1), xz[3][1]);
			put(IL(ix-2, jy, kz  ), xz[0][2]);
			put(IL(ix-1, jy, kz  ), xz[1][2]);
			put(IL(ix+1, jy, kz  ), xz[3][2]);
			put(IL(ix+2, jy, kz  ), xz[4][2]);
			put(IL(ix-1, jy, kz+1), xz[1][3]);
			put(IL(ix  , jy, kz+1), xz[2][3]);
			put(IL(ix+1, jy, kz+1), xz[3][3]);
			put(IL(ix  , jy, kz+2), xz[2][4]);
			if(has_y)
			{
				put(IL(ix, jy-2, kz  ), yz[0][1]);
				put(IL(ix, jy-1, kz  ), yz[1][1]);
				put(IL(ix, jy+1, kz  ), yz[3][1]);
				put(IL(ix, jy+2, kz  ), yz[4][1]);
				put(IL(ix, jy-1, kz-1), yz[1][0]);
				put(IL(ix, jy+1, kz-1), yz[3][0]);
				put(IL(ix, jy-1, kz+1), yz[1][2]);
				put(IL(ix, jy+1, kz+1), yz[3][2]);
			}
			_react[irow] = -react*rh;
			_crsA.Diag(irow) = -offsum + _react[irow];
			_VecB(irow) = -src*rh - (_mms.Enabled() ? mmsSrc[kz-2] : 0.0);
		}
	}

	SetBoundaryCondition();
}

//  The weights w[0..3] of the four stencil points in the face value of the
//  advection, for a flow F through the face (toward stencil point 2 when
//  F >= 0), seen from the row of stencil point rowPos, and the part k of the
//  face value that does not depend on c. A linear scheme: its own weights,
//  the same for both rows. A limited one: the limiter frozen at the field,
//  in the upwind form for the upwind cell and the downwind form for the
//  other (FaceWeights), which makes the frozen matrix an M-matrix and keeps
//  c bounded; the two forms give the face the same value only at the field
//  they were frozen at, so the step is solved again until its midpoint is
//  that field (SettleLimiter). (A single form for both rows -- the upwind
//  value implicitly and the limited correction from the field, or the face
//  value linearized about the field -- is conservative at any c but not
//  bounded: a front carried at a Courant number of 0.6 overshoots by 2e-3,
//  tools/advection.py front_cz.)
void CZ_ConcEquation_2D::AdvWeights(const double f[4], double F, double w[4], double &k, int rowPos) const
{
	k = 0.0;
	_scheme.FaceWeights(f, F, _scheme.IsLimited() ? rowPos : -1, w);
}

//  The ghost rows. The two ghost layers at each end of a column are the
//  images of the two cells inside, c(1) = c(2) - dzeta c_zeta(0), c(0) =
//  c(3) - 3 dzeta c_zeta(0), and at the top c(_Nz-2) = c(_Nz-3) + dzeta
//  c_zeta(1), c(_Nz-1) = c(_Nz-4) + 3 dzeta c_zeta(1), with the zeta gradient
//  that makes the conormal flux what it is there:
//
//      G = Df D [ (As + (h A)^2 + (h Ay)^2) c_zeta / h - h A c_x - h Ay c_y ],
//
//  G = 0 on the substrate (h A = s_x) and G = W_T c_s at the free surface
//  (h A = h_x + s_x), W_T the solvent leaving through it (with the
//  manufactured solution, G = its own G there, and W_T c_s - (W_T c_s - G)
//  at the top). c_x, c_y, c_s and D are those of the field (frozen, as the
//  other coefficients). Without a slope this is D c_zeta = Pe E h c at the
//  surface and c_zeta = 0 at the substrate. The ghosts are read by the
//  advection stencils next to the ends and by the diffusive flux of the x
//  faces of the end cells; the flux through the ends themselves is set in
//  the rows. At both ends in x and y they are copied from the opposite side
//  if periodic, and mirrored otherwise.
void CZ_ConcEquation_2D::SetBoundaryCondition(bool atEnd)
{
	const double dzeta = _G3->dz();
	const double dx = _G3->dx();
	const double As = _NonDim->As();
	const double Df = _NonDim->Df();

	const bool has_y = (_Ny != 1);
	const double dy = has_y ? _G3->dy() : 1.0;
	const int jlo = has_y ? 2      : 0;
	const int jhi = has_y ? _Ny-2  : 1;
	//  The conditions of the step are taken at its midpoint, with the
	//  coefficients of the stage; atEnd: of the state the field holds at
	//  the end of the step (RefreshGhosts).
	const double tmid = _timeStep->NowTime() + (atEnd ? 1.0 : 0.5)*_timeStep->DT().now;

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for
#endif
	for(int ix = 2; ix < _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		auto C  = [&](int i, int j, int k) { return _field->CC(_G3->At<Cell>(i, j, k)); };
		auto Hh = [&](int i, int j) { return _field->HH(_G2->At<Zface>(i, j)); };
		auto Ss = [&](int i, int j) { return _field->SS(_G2->At<Zface>(i, j)); };
		const double hh = Hh(ix, jy);
		const double xc = _G2->X(_G2->At<Zface>(ix, jy));
		const double yc = has_y ? _G2->Y(_G2->At<Zface>(ix, jy)) : 0.0;
		const double hx0 = (Hh(ix-2,jy) -8.0*Hh(ix-1,jy) +8.0*Hh(ix+1,jy) -Hh(ix+2,jy))/(12.0*dx);
		const double sx0 = (Ss(ix-2,jy) -8.0*Ss(ix-1,jy) +8.0*Ss(ix+1,jy) -Ss(ix+2,jy))/(12.0*dx);
		double hy0 = 0, sy0 = 0;
		if(has_y)
		{
			hy0 = (Hh(ix,jy-2) -8.0*Hh(ix,jy-1) +8.0*Hh(ix,jy+1) -Hh(ix,jy+2))/(12.0*dy);
			sy0 = (Ss(ix,jy-2) -8.0*Ss(ix,jy-1) +8.0*Ss(ix,jy+1) -Ss(ix,jy+2))/(12.0*dy);
		}
		double gB = 0.0, gT = 0.0, csT = 0.0;
		if(_mms.Enabled())
		{
			gB = _mms.CZFluxGNd(xc, yc, 0.0, tmid);
			gT = _mms.CZFluxGNd(xc, yc, 1.0, tmid);
			csT = _mms.CZnd(xc, yc, 1.0, tmid);
		}

		//  The substrate
		{
			const int k = 2;
			const double cb = 0.5*(C(ix,jy,k-1) + C(ix,jy,k));
			const double D  = Df * _NonDim->Diffusivity(cb);
			const double cx = (C(ix+1,jy,k) - C(ix-1,jy,k))/(2.0*dx);
			const double cy = has_y ? (C(ix,jy+1,k) - C(ix,jy-1,k))/(2.0*dy) : 0.0;
			const double hA = sx0, hAy = sy0;
			const double cz = hh*(gB/D + hA*cx + hAy*cy)/(As + hA*hA + hAy*hAy);
			_crsA.Diag(IL(ix, jy, 0)) = 1.0;		_crsA(IL(ix, jy, 0), IL(ix, jy, 3)) = -1.0;
			_Rhs(IL(ix, jy, 0)) = -3.0*dzeta*cz;
			_crsA.Diag(IL(ix, jy, 1)) = 1.0;		_crsA(IL(ix, jy, 1), IL(ix, jy, 2)) = -1.0;
			_Rhs(IL(ix, jy, 1)) = -dzeta*cz;
		}
		//  The free surface
		{
			const int k = _Nz-3;
			const double cs = 0.5*(C(ix,jy,k) + C(ix,jy,k+1));
			const double D  = Df * _NonDim->Diffusivity(cs);
			const double cx = (C(ix+1,jy,k) - C(ix-1,jy,k))/(2.0*dx);
			const double cy = has_y ? (C(ix,jy+1,k) - C(ix,jy-1,k))/(2.0*dy) : 0.0;
			const double hA = hx0 + sx0, hAy = hy0 + sy0;
			//  The solvent leaving through the surface: that of the step
			//  (the top face of the column's volume balance), or at the end
			//  of the step the evaporation rate of the end state (less the
			//  manufactured source of the film, which the step's flux has
			//  in it too)
			double WT = _W.empty() ? 0.0 : _W[IL(ix, jy, _Nz-2)];
			if(atEnd)
			{
				WT = _NonDim->Ev() * _NonDim->EvapRate(*_field, _G2->At<Zface>(ix, jy));
				//  limited to the solvent the column holds, as the
				//  film-height equation limits it
				double cm = 0.0;
				for(int kk = 2; kk < _Nz-2; kk++) cm += C(ix, jy, kk);
				cm /= (double)(_Nz-4);
				WT = C_Evaporation::SolventLimit(WT, hh, cm, _timeStep->DT().now);
				if(_mms.Enabled()) WT -= _mms.Snd(xc, yc, tmid);
			}
			const double G  = _mms.Enabled() ? WT*(cs - csT) + gT : WT*cs;
			const double cz = hh*(G/D + hA*cx + hAy*cy)/(As + hA*hA + hAy*hAy);
			const double grad = dzeta*cz;
			//  The two ghosts are the images of the two top cells across the
			//  surface, as at the substrate (0 <- 3, 1 <- 2): c(_Nz-2) - c(_Nz-3)
			//  and c(_Nz-1) - c(_Nz-4) are the gradient over one and over three
			//  cells (the outer one taken as c(_Nz-2) + grad was right only for
			//  a linear profile, an O(dzeta) error in the advection of the top
			//  cell).
			_crsA.Diag(IL(ix, jy, _Nz-2)) = 1.0;	_crsA(IL(ix, jy, _Nz-2), IL(ix, jy, _Nz-3)) = -1.0;
			_Rhs(IL(ix, jy, _Nz-2)) = grad;
			_crsA.Diag(IL(ix, jy, _Nz-1)) = 1.0;	_crsA(IL(ix, jy, _Nz-1), IL(ix, jy, _Nz-4)) = -1.0;
			_Rhs(IL(ix, jy, _Nz-1)) = 3.0*grad;
		}
	}
	int ix, jy, kz;
	//  The x ghost layers. In two dimensions they are placed on every y row.
	const bool perX = (_boundary->TypeXmin()==BoundaryCondition::Periodic
	                && _boundary->TypeXmax()==BoundaryCondition::Periodic);
	const int jAll_lo = has_y ? 0     : 0;
	const int jAll_hi = has_y ? _Ny   : 1;

	for(kz = 0; kz < _Nz; kz++)
	for(jy = jAll_lo; jy < jAll_hi; jy++)
	{
		if(perX)
		{
			_crsA.Diag(IL(0    , jy, kz)) = 1.0;		_crsA(IL(0    , jy, kz), IL(_Nx-4, jy, kz)) = -1.0;
			_crsA.Diag(IL(1    , jy, kz)) = 1.0;		_crsA(IL(1    , jy, kz), IL(_Nx-3, jy, kz)) = -1.0;
			_crsA.Diag(IL(_Nx-2, jy, kz)) = 1.0;		_crsA(IL(_Nx-2, jy, kz), IL(2    , jy, kz)) = -1.0;
			_crsA.Diag(IL(_Nx-1, jy, kz)) = 1.0;		_crsA(IL(_Nx-1, jy, kz), IL(3    , jy, kz)) = -1.0;
		}
		else
		{
			_crsA.Diag(IL(0    , jy, kz)) = 1.0;		_crsA(IL(0    , jy, kz), IL(3    , jy, kz)) = -1.0;
			_crsA.Diag(IL(1    , jy, kz)) = 1.0;		_crsA(IL(1    , jy, kz), IL(2    , jy, kz)) = -1.0;
			_crsA.Diag(IL(_Nx-2, jy, kz)) = 1.0;		_crsA(IL(_Nx-2, jy, kz), IL(_Nx-3, jy, kz)) = -1.0;
			_crsA.Diag(IL(_Nx-1, jy, kz)) = 1.0;		_crsA(IL(_Nx-1, jy, kz), IL(_Nx-4, jy, kz)) = -1.0;
		}
	}

	if(!has_y) return;

	//  The y ghost layers. The corners have already been written on the x side,
	//  so only the interior x values are looped over.
	const bool perY = (_boundary->TypeYmin()==BoundaryCondition::Periodic
	                && _boundary->TypeYmax()==BoundaryCondition::Periodic);

	for(kz = 0; kz < _Nz; kz++)
	for(ix = 2; ix < _Nx-2; ix++)
	{
		if(perY)
		{
			_crsA.Diag(IL(ix, 0    , kz)) = 1.0;		_crsA(IL(ix, 0    , kz), IL(ix, _Ny-4, kz)) = -1.0;
			_crsA.Diag(IL(ix, 1    , kz)) = 1.0;		_crsA(IL(ix, 1    , kz), IL(ix, _Ny-3, kz)) = -1.0;
			_crsA.Diag(IL(ix, _Ny-2, kz)) = 1.0;		_crsA(IL(ix, _Ny-2, kz), IL(ix, 2    , kz)) = -1.0;
			_crsA.Diag(IL(ix, _Ny-1, kz)) = 1.0;		_crsA(IL(ix, _Ny-1, kz), IL(ix, 3    , kz)) = -1.0;
		}
		else
		{
			_crsA.Diag(IL(ix, 0    , kz)) = 1.0;		_crsA(IL(ix, 0    , kz), IL(ix, 3    , kz)) = -1.0;
			_crsA.Diag(IL(ix, 1    , kz)) = 1.0;		_crsA(IL(ix, 1    , kz), IL(ix, 2    , kz)) = -1.0;
			_crsA.Diag(IL(ix, _Ny-2, kz)) = 1.0;		_crsA(IL(ix, _Ny-2, kz), IL(ix, _Ny-3, kz)) = -1.0;
			_crsA.Diag(IL(ix, _Ny-1, kz)) = 1.0;		_crsA(IL(ix, _Ny-1, kz), IL(ix, _Ny-4, kz)) = -1.0;
		}
	}
}




//  One step of the concentration, second order in dt.
//
//  The operator is formed with the film halfway through the step (the height
//  equation sets it before calling this) and its velocities at the midpoint.
//  Its concentration-dependent parts -- the diffusivity, the limiter of a
//  limited scheme and the ghosts of the conormal conditions -- are formed from
//  the field's concentration; solved once, that is the concentration at the
//  start of the step and the step is first order. So a predictor is solved
//  with it, and then the same system again from the same start with the
//  concentration of the field set halfway to the prediction.
//
//  Two things the base class does would each keep it first order: it forms
//  the explicit part A c_old with the matrix of the step before, and the
//  ghost cells of c_old carry the surface condition of the step before (and
//  at the first step none at all, which loses half a step of solute). Here
//  A c_old is formed with this stage's matrix, after the ghosts of c_old
//  have been set from this stage's boundary rows.
int CZ_ConcEquation_2D::SolveStage(int &iter, double &resid)
{
	SetMatrixA();
	ApplyBoundaryRows(_phi_old);
	ExplicitPart(_phi_old, _Ax_old);
	const double dt = _timeStep->DT().now;
	const int n = _crsA.RowSize();
	if(Math::IsZero(_theta))
	{
		_phi_new = _phi_old;
		for(int r = 0; r < n; r++)
			if(!_crsA.IsBC(r)) _phi_new[r] = _phi_old[r] - dt*(_Ax_old[r] + _VecB[r]);
		ApplyBoundaryRows(_phi_new);
		iter = 0; resid = 0;
		return 0;
	}
	MakeMatrix_IA();
	//  Solved for the change over the step, (I + theta dt A) d = -dt (A c_old
	//  + B), d = 0 in the ghost rows (c_old already satisfies them): a right
	//  side that is zero -- a uniform concentration that nothing changes --
	//  gives no change at all, and the solver's tolerance is relative to the
	//  change rather than to c.
	Vector<double> rhs = _phi_old, d = _phi_old;
	rhs = 0; d = 0;
	bool any = false;
	for(int r = 0; r < n; r++)
	{
		if(_crsA.IsBC(r)) continue;
		rhs[r] = -dt*(_Ax_old[r] + _VecB[r]);
		if(rhs[r] != 0.0) any = true;
	}
	int ch = 0;
	iter = 0; resid = 0;
	if(any) { YUZEN_PROF("c.linsolve"); ch = _bcg.Solve(_crsIA, rhs, d, iter, resid); }
	if(ch != 0)
	{
		Global::Warning<<"Semi-Implicit: BiCGStab did not converge: ret = "<< ch << ", iter = "<<iter<< ", resid = "<< std::scientific <<resid<< nl;
		bool finite = true;
		for(int i = 0; i < d.size() && finite; i++) finite = std::isfinite(d[i]);
		if(!finite) d = 0;
	}
	//  The solute of the step is sum h d over the film, which the rows make
	//  the sum of face fluxes that cancel -- for the exact solution. The
	//  solver stops at a residual r, and sum h r, which the stiffness of the
	//  vertical diffusion makes 1e-14 of the solute a step, is lost or
	//  gained. That sum is put back evenly over the film: d + alpha, with r
	//  = -dt (A c_theta + B) - d, c_theta = c_old + theta d (the flux
	//  balance of the rows, in the difference form). The shift moves
	//  c_theta too, and A of a uniform shift is the surface term alone (the
	//  solvent leaving the top cells, _react), so r goes to r - alpha (1 +
	//  theta dt react), and alpha = sum h r / sum h (1 + theta dt react)
	//  takes sum h r to zero. (sum h r / sum h would leave the part
	//  theta dt react of the shift in r: in a step that evaporates a sizeable
	//  part of the film (0.3 of it) a quarter of sum h r stays -- 1e-10 of
	//  the solute in one step, tools/proptest.py seeds 1054 and 2105.)
	//  (Taking that flux balance itself as d made the solute exact, but it
	//  is the solution plus r, and r is not small pointwise: the mirror
	//  image of a run differed by 1e-10 in c and, through the Marangoni
	//  stress, by 1e-8 in h. The error of d itself is r smoothed by the
	//  inverse of the matrix.) A uniform concentration that nothing changes
	//  has d = 0 and r = 0.
	if(any && ch == 0)
	{
		//  The ghosts of d as the rows give them from its interior (the
		//  solver holds them only to its residual): a flux across a periodic
		//  seam is formed on each side from a ghost there, and the two agree
		//  only if the ghosts are the copies.
		{
			Vector<double> cn = _phi_old;
			cn += d;
			ApplyBoundaryRows(cn);
			d = cn;
			d -= _phi_old;
		}
		Vector<double> ct = d;
		ct *= _theta;
		ct += _phi_old;
		Vector<double> Ac = ct;
		ExplicitPart(ct, Ac);
		const bool has_y = (_Ny != 1);
		const int jlo = has_y ? 2      : 0;
		const int jhi = has_y ? _Ny-2  : 1;
		double shr = 0.0, sh = 0.0;
		for(int jy = jlo; jy < jhi; jy++)
		for(int ix = 2; ix < _Nx-2; ix++)
		{
			const double h = _field->HH(_G2->At<Zface>(ix, jy));
			for(int kz = 2; kz < _Nz-2; kz++)
			{
				const int r = IL(ix, jy, kz);
				shr += h * (-dt*(Ac[r] + _VecB[r]) - d[r]);
				//  A of a uniform shift is the surface term alone (_react)
				sh  += h * (1.0 + _theta*dt*_react[r]);
			}
		}
		//  Every row, the ghosts too: a uniform shift satisfies the ghost
		//  rows, so they keep their copies
		const double alpha = (sh > 0.0) ? shr / sh : 0.0;
		for(int r = 0; r < n; r++) d[r] += alpha;
	}
	_phi_new = _phi_old;
	_phi_new += d;
	ApplyBoundaryRows(_phi_new);
	return ch;
}

//  A x for the interior rows as sum_j a_ij (x_j - x_i) + r_i x_i (the
//  diagonal is -sum_j a_ij + r_i, see SetMatrixA), so that a uniform x gives
//  r_i x_i exactly and not the rounding of a sum that cancels; zero in the
//  ghost rows.
void CZ_ConcEquation_2D::ExplicitPart(const Vector<double> &x, Vector<double> &Ax) const
{
	const int n = _crsA.RowSize();
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(n))
#endif
	for(int r = 0; r < n; r++)
	{
		if(_crsA.IsBC(r)) { Ax[r] = 0.0; continue; }
		const double xr = x[r];
		double s = _react[r]*xr;
		for(int j = 1; j < _crsA.NonZero(r); j++)
			s += _crsA.Element(r, j) * (x[_crsA.Jcol(r, j)] - xr);
		Ax[r] = s;
	}
}

int CZ_ConcEquation_2D::Solve(int &iter, double &resid)
{
	assert(_initialized);
	int ch = SolveStage(iter, resid);
	if(ch != 0) return ch;

	Vector<double> mid = _phi_old;
	mid += _phi_new;
	mid *= 0.5;
	_field->Update_C(mid);

	int iter2 = 0;
	double resid2 = 0;
	ch = SolveStage(iter2, resid2);
	iter += iter2;
	resid = resid2;
	if(ch == 0) ch = SettleLimiter(iter, resid);
	_field->Update_C(_phi_old);
	return ch;
}

int CZ_ConcEquation_2D::SolveAtMidpoint(const Vector<double> &cmid, int &iter, double &resid)
{
	assert(_initialized);
	Vector<double> c = cmid;
	_field->Update_C(c);
	int ch = SolveStage(iter, resid);
	if(ch == 0) ch = SettleLimiter(iter, resid);
	_field->Update_C(_phi_old);
	return ch;
}

//  A limited scheme: the step solved again with the limiter frozen at the
//  midpoint of what came out, until that midpoint settles (to 1e-12 of the
//  change of the step, kLimiterPasses at most). The face values of the
//  upwind and the downwind form then agree (AdvWeights), and with them the
//  flux the two cells of a face exchange: the solute is conserved to that
//  tolerance, and the step is the midpoint rule of the limited scheme.
int CZ_ConcEquation_2D::SettleLimiter(int &iter, double &resid)
{
	if(!_scheme.IsLimited()) return 0;
	for(int pass = 0; pass < kLimiterPasses; pass++)
	{
		Vector<double> mid = _phi_old;
		mid += _phi_new;
		mid *= 0.5;
		const Vector<double> prev = _phi_new;
		_field->Update_C(mid);
		int it = 0; double rs = 0;
		const int ch = SolveStage(it, rs);
		iter += it; resid = rs;
		if(ch != 0) return ch;
		double dmax = 0, cmax = 0;
		for(int i = 0; i < _phi_new.size(); i++)
		{
			dmax = std::max(dmax, std::fabs(_phi_new[i] - prev[i]));
			cmax = std::max(cmax, std::fabs(_phi_new[i] - _phi_old[i]));
		}
		if(!(dmax > 1e-12*cmax)) break;
	}
	return 0;
}

//  The ghosts of the concentration at the end of the step, from the
//  conditions of the end state (the field must hold it: h and c of the end
//  of the step). Those the step was solved with are the conditions of its
//  midpoint, which the ghosts of the end state differ from by O(dt); the
//  surface concentration the next step starts from (surface tension,
//  Marangoni stress, evaporation) is read from them. The rows are left as
//  they are; nothing but the ghosts changes.
void CZ_ConcEquation_2D::RefreshGhosts()
{
	SetBoundaryCondition(true);
	ApplyBoundaryRows(_phi_new);
}

//  Make x satisfy the boundary rows of the matrix: each ghost row reads
//  a_rr x_r + sum a_rj x_j = Rhs_r. A ghost can depend on another ghost (the
//  outer layer on the inner one, a corner on an edge), so the rows are swept
//  until nothing changes, which takes at most as many sweeps as the chain is
//  long.
void CZ_ConcEquation_2D::ApplyBoundaryRows(Vector<double> &x) const
{
	for(int sweep = 0; sweep < 4; sweep++)
	{
		for(int r = 0; r < _crsA.RowSize(); r++)
		{
			if(!_crsA.IsBC(r)) continue;
			double sum = _Rhs[r];
			for(int j = 1; j < _crsA.NonZero(r); j++)
				sum -= _crsA.Element(r, j) * x[_crsA.Jcol(r, j)];
			x[r] = sum / _crsA.Diag(r);
		}
	}
}
