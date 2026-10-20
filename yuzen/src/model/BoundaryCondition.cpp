/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  BoundaryCondition.cpp : implementation of the application of the boundary conditions
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "Opt.h"
#include "YuzenUtil.h"
using namespace yuzen::util;
#include <iostream>
#include <fstream>

#include "BoundaryCondition.h"
#include "LinearSolver.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

void BoundaryCondition::Dirichlet(CRS_Matrix &crs, Vector<double>& rhs, int i_diag, int col_st, int step, double b)
{
	double cM, cm, cp, cP;
	cM = -1.0/16.0;
	cm = +9.0/16.0;
	cp = +9.0/16.0;
	cP = -1.0/16.0;

	int jcol = col_st;
	if(jcol != i_diag) crs(i_diag, jcol) = cM; else crs.Diag(i_diag) = cM; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cm; else crs.Diag(i_diag) = cm; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cp; else crs.Diag(i_diag) = cp; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cP; else crs.Diag(i_diag) = cP; jcol+=step;
	rhs[i_diag] = b;
}
void BoundaryCondition::Neumann(CRS_Matrix &crs, Vector<double>& rhs,   int i_diag, int col_st, int step, double b, double dx)
{
	double cM, cm, cp, cP;
	cM = +1.0/24.0;
	cm = -27.0/24.0;
	cp = +27.0/24.0;
	cP = -1.0/24.0;

	int jcol = col_st;
	if(jcol != i_diag) crs(i_diag, jcol) = cM; else crs.Diag(i_diag) = cM; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cm; else crs.Diag(i_diag) = cm; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cp; else crs.Diag(i_diag) = cp; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cP; else crs.Diag(i_diag) = cP; jcol+=step;
	rhs[i_diag] = dx * b;
}
void BoundaryCondition::Derive2(CRS_Matrix &crs, Vector<double>& rhs,    int i_diag, int col_st, int step, double b, double dx)
{
	double cM, cm, cp, cP;
	cM = +1.0/2.0;
	cm = -1.0/2.0;
	cp = -1.0/2.0;
	cP = +1.0/2.0;
	int jcol = col_st;
	if(jcol != i_diag) crs(i_diag, jcol) = cM; else crs.Diag(i_diag) = cM; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cm; else crs.Diag(i_diag) = cm; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cp; else crs.Diag(i_diag) = cp; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cP; else crs.Diag(i_diag) = cP; jcol+=step;
	rhs[i_diag] = dx*dx * b;
}
void BoundaryCondition::Derive3(CRS_Matrix &crs, Vector<double>& rhs,    int i_diag, int col_st, int step, double b, double dx)
{
	double cM, cm, cp, cP;
	cM = -1.0;
	cm = +3.0;
	cp = -3.0;
	cP = +1.0;
	int jcol = col_st;
	if(jcol != i_diag) crs(i_diag, jcol) = cM; else crs.Diag(i_diag) = cM; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cm; else crs.Diag(i_diag) = cm; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cp; else crs.Diag(i_diag) = cp; jcol+=step;
	if(jcol != i_diag) crs(i_diag, jcol) = cP; else crs.Diag(i_diag) = cP; jcol+=step;
	rhs[i_diag] = dx*dx*dx * b;
}
void BoundaryCondition::EqualVal(CRS_Matrix &crs, Vector<double>& rhs,    int i_diag, int jcol)
{
	crs.Diag(i_diag)	=  1.0;
	crs(i_diag, jcol)	= -1.0;
	rhs[i_diag]			=  0.0;
}

//==================================================================
//  Set up the matrix for one boundary.
//
//  Xmin / Xmax / Ymin / Ymax (only two directions in 1D) differ in three
//  respects only:
//    - whether the row numbers of the ghost layers are IL(0,jy), IL(_Nx-1,jy),
//      IL(ix,0), and so on
//    - whether the 4-point stencil starts from the inside or from the outside
//    - whether the column stride is 1 or _Nx
//  Those three are collected into EdgeSlot, so the body is written only once.
//
//  *** Do not change the order of insertion into the matrix. ***
//  If the storage order of CRS_Matrix changes, the rounding of the sums in
//  ScalingDiag / Residual changes with it, and the results no longer agree
//  bit for bit.
//==================================================================
void BoundaryCondition::ApplyEdge(CRS_Matrix &crs, Vector<double>& rhs, tField tf,
                                  tBoundaryCondition tbc, const EdgeSlot &e, tEdgeForm form)
{
	const int *L  = e.L;
	const int  cs = e.colStart;
	const int  st = e.step;

	switch(tbc)
	{
	default:
	case Symmetry:
		//  The ghosts are the mirror image of the first two layers, in 1-D
		//  as in 2-D, so that dh/dn = 0 and d3h/dn3 = 0 at the face. Written
		//  as those two derivatives instead (as 1-D was), the two ghost rows
		//  are coupled and a solve, or CRS_Matrix::SatisfyBCRows, leaves the
		//  ghosts off the mirror by rounding of h: the capillary flux of the
		//  edge face, a third difference over them, is then a leak of volume
		//  of eps h D / dx^4 (tools/analytic.py nozzle_balance).
		EqualVal(crs, rhs, L[0], L[3]);
		EqualVal(crs, rhs, L[1], L[2]);
		break;

	case Periodic:
		// The outermost ghost takes the innermost point on the opposite side, and
		// the inner ghost the one just outside that.
		EqualVal(crs, rhs, L[0], e.opp[3]);
		EqualVal(crs, rhs, L[1], e.opp[2]);
		break;

	case Pinning:
		switch(tf)
		{
		case Height:
		case HeightJacobi:
			// 0.5*(h[layer 1] + h[layer 2]) = 0
			crs.Diag(L[1]) = 0.5;	crs(L[1], L[2]) = 0.5;

			// Ux = 0
			crs.Diag(L[0]) = 1.0;
			crs(L[0], L[1]) = -3.0;
			crs(L[0], L[2]) = +3.0;
			crs(L[0], L[3]) = -1.0;
			break;
		case Conc:
			if(form == Form2D)
			{
				EqualVal(crs, rhs, L[1], L[2]);
				EqualVal(crs, rhs, L[0], L[1]);
			}
			else
			{
				Neumann(crs, rhs, L[0], cs, st, 0, 0);
				Derive2(crs, rhs, L[1], cs, st, 0, 0);
			}
			break;
		}
		break;

	case FreshInlet:
		switch(tf)
		{
		case Height:
			Dirichlet(crs, rhs, L[1], cs, st, e.hb);	// the height itself is held fixed
			Derive3(  crs, rhs, L[0], cs, st, 0, 0);
			break;
		case HeightJacobi:
			Dirichlet(crs, rhs, L[1], cs, st, 0.0);	// the correction is zero
			Derive3(  crs, rhs, L[0], cs, st, 0, 0);
			break;
		case Conc:
			Dirichlet(crs, rhs, L[1], cs, st, e.cb);
			crs.Diag(L[0])  = +1.0;
			crs(L[0], L[1]) = -2.0;
			crs(L[0], L[2]) = +1.0;
			break;
		}
		break;

	case Outlet:
		Derive2(crs, rhs, L[0], cs, st, 0, 0);
		Derive3(crs, rhs, L[1], cs, st, 0, 0);
		break;
	}
}

void BoundaryCondition1D::ApplyToField(Field_WM1D& field)
{
	C_Grid1D *G = field.GetGridPtr();
	int nx = G->nx();

	// Temporary arrays
	Vector<double> tmpH(nx), tmpS(nx), tmpC(nx), SS(nx), HH(nx), CC(nx);
	tmpH = 0; tmpS = 0; tmpC = 0;
	//Set up the matrices for the boundary conditions.
	CRS_Matrix crsBC_H, crsBC_S, crsBC_C;
	crsBC_H.Build(nx);crsBC_C.Build(nx);
	for(Xedge xe = G->First<Xedge>(); xe != G->Last<Xedge>(); ++xe)
	{
		if(G->CheckGhost(xe) != G->Inside)
		{	//Flag the boundary-condition part.
			crsBC_H.SetAsBC(G->ix(xe));
			crsBC_C.SetAsBC(G->ix(xe));
		}
		else
		{
			//Everything other than the boundary becomes the identity, and the right-hand side is assigned as it stands.
			crsBC_H.Diag(G->ix(xe)) = 1.0;
			crsBC_C.Diag(G->ix(xe)) = 1.0;
			tmpH[G->ix(xe)] = field.HH(xe);
			tmpC[G->ix(xe)] = field.CC(xe);
			tmpS[G->ix(xe)] = field.SS(xe);
		}
		HH[G->ix(xe)] = field.HH(xe);
		CC[G->ix(xe)] = field.CC(xe);
		SS[G->ix(xe)] = field.SS(xe);
	}
	crsBC_S = crsBC_H;	//Copy only the structure over to the substrate-shape matrix, before the height boundary conditions are set.

	//Set up the matrix for the boundary part.
	Apply(crsBC_H, tmpH, BoundaryCondition::Height, field);
	Apply(crsBC_C, tmpC, BoundaryCondition::Conc, field);
	//Solve the matrix so that everything outside the boundary-condition part is left as it is, while the given boundary conditions are forced onto the boundary part.
	LinearSolver bcg;
	int nz = nx + 4*nx, iter = 0;
	double resid = 0;
	bcg.Build(nx, nz, nz, nz);
	bcg.MaxIteration() = nx;
	bcg.Solve(crsBC_H, tmpH, HH, iter, resid);
	crsBC_H.SatisfyBCRows(tmpH, HH);
	bcg.Solve(crsBC_C, tmpC, CC, iter, resid);
	crsBC_C.SatisfyBCRows(tmpC, CC);

	//  The substrate beyond an edge is its mirror image, except across a
	//  periodic pair, where it is the substrate at the other end. Mirrored
	//  there too, it would put a kink in the free surface at the seam, and the
	//  capillary pressure of that kink would drive a flow out of the edge that
	//  a film moved round the domain does not have.
	BoundaryCondition1D bcS; bcS.Build(nx, SubstrateBC(_tXmin), SubstrateBC(_tXmax));
	bcS.Apply(crsBC_S, tmpS, BoundaryCondition::Height, field);
	bcg.Solve(crsBC_S, tmpS, SS, iter, resid);
	crsBC_S.SatisfyBCRows(tmpS, SS);
	for(Xedge xe = G->First<Xedge>(); xe != G->Last<Xedge>(); ++xe)
	{
		//  The whole initial field is re-solved as a linear system in which the
		//  interior is the identity and only the boundary carries the boundary
		//  conditions, and the result is written back. The interior values do
		//  not change, so in effect this just fills the ghost cells.
		//
		//  Unresolved: forcing the boundary conditions on through a linear
		//  solver like this costs more than filling the ghost cells directly,
		//  and the one benefit - that the definition of the boundary conditions
		//  lives in one place, as the shape of the matrix in Apply() - does not
		//  count for much, since the time integration calls Apply() separately
		//  anyway. It has been left alone because it is only reached once, at
		//  initialization.
		field.HH(xe) = HH[G->ix(xe)];
		field.CC(xe) = CC[G->ix(xe)];
		field.SS(xe) = SS[G->ix(xe)];
	}
}

void BoundaryCondition1D::Apply(CRS_Matrix &crs, Vector<double>& rhs, tField tf, Field_WM1D& ini_Field)
{
	crs.SetAsBC(0);
	crs.SetAsBC(1);
	crs.SetAsBC(_Nx-2);
	crs.SetAsBC(_Nx-1);

	C_Grid1D *G = ini_Field.GetGridPtr();
	Xedge xe0, xe1;

	//---- Xmin ----
	EdgeSlot xm;
	xm.L[0] = 0;  xm.L[1] = 1;  xm.L[2] = 2;  xm.L[3] = 3;
	xm.opp[0] = _Nx-1; xm.opp[1] = _Nx-2; xm.opp[2] = _Nx-3; xm.opp[3] = _Nx-4;
	xm.colStart = xm.L[0];	xm.step = 1;

	xe0 = G->At<Xedge>(1);
	xe1 = G->At<Xedge>(2);
	xm.hb = 0.5 * (ini_Field.HH(xe0) + ini_Field.HH(xe1));
	xm.cb = 0.5 * (ini_Field.CC(xe0) + ini_Field.CC(xe1));

	//---- Xmax ----
	EdgeSlot xp;
	xp.L[0] = _Nx-1; xp.L[1] = _Nx-2; xp.L[2] = _Nx-3; xp.L[3] = _Nx-4;
	xp.opp[0] = 0;  xp.opp[1] = 1;  xp.opp[2] = 2;  xp.opp[3] = 3;
	//  On all four sides the 4-point stencil starts its columns from the
	//  outermost ghost (starting from _Nx-3 would put column _Nx outside the
	//  range of the matrix).
	xp.colStart = xp.L[3];	xp.step = 1;

	//  The inlet value is the mean of the two points either side of the face,
	//  layers 1 and 2, as at Xmin, not layers 0 and 1 -- the two ghosts --
	//  which would hold an inlet at Xmax at a value half a cell outside the
	//  domain, and one whose initial profile is not flat at another value
	//  than the same inlet at Xmin.
	xe0 = G->At<Xedge>(_Nx-2);
	xe1 = G->At<Xedge>(_Nx-3);
	xp.hb = 0.5 * (ini_Field.HH(xe0) + ini_Field.HH(xe1));
	xp.cb = 0.5 * (ini_Field.CC(xe0) + ini_Field.CC(xe1));

	ApplyEdge(crs, rhs, tf, _tXmin, xm, Form1D);
	ApplyEdge(crs, rhs, tf, _tXmax, xp, Form1D);
}

void BoundaryCondition2D::ApplyToField(Field_WM2D& field)
{
	C_Grid2D *G = field.GetGridPtr();
	int nx = G->nx();
	int ny = G->ny();

	// Temporary arrays
	Vector<double> tmpH(nx*ny), tmpS(nx*ny), tmpC(nx*ny), SS(nx*ny), HH(nx*ny), CC(nx*ny);
	tmpH = 0; tmpS = 0; tmpC = 0;
	//Set up the matrices for the boundary conditions.
	CRS_Matrix crsBC_H, crsBC_S, crsBC_C;
	crsBC_H.Build(nx*ny);crsBC_C.Build(nx*ny);
	for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf)
	{
		if(G->CheckGhost(zf) != G->Inside)
		{	//Flag the boundary-condition part.
			crsBC_H.SetAsBC(G->IL(zf));
			crsBC_C.SetAsBC(G->IL(zf));
		}
		else
		{
			//Everything other than the boundary becomes the identity, and the right-hand side is assigned as it stands.
			crsBC_H.Diag(G->IL(zf)) = 1.0;
			crsBC_C.Diag(G->IL(zf)) = 1.0;
			tmpH[G->IL(zf)] = field.HH(zf);
			tmpC[G->IL(zf)] = field.CC(zf);
			tmpS[G->IL(zf)] = field.SS(zf);
		}
		HH[G->IL(zf)] = field.HH(zf);
		CC[G->IL(zf)] = field.CC(zf);
		SS[G->IL(zf)] = field.SS(zf);
	}
	crsBC_S = crsBC_H;	//Copy only the structure over to the substrate-shape matrix, before the height boundary conditions are set.

	//Set up the matrix for the boundary part.
	Apply(crsBC_H, tmpH, BoundaryCondition::Height, field);
	Apply(crsBC_C, tmpC, BoundaryCondition::Conc, field);

	//Solve the matrix so that everything outside the boundary-condition part is left as it is, while the given boundary conditions are forced onto the boundary part.
	LinearSolver bcg;
	int nz = nx*ny + 4*nx +4*ny, iter = 0;
	double resid = 0;
	bcg.Build(nx*ny, nz, nz, nz);
	bcg.MaxIteration() = nx*ny;
	bcg.Solve(crsBC_H, tmpH, HH, iter, resid);
	crsBC_H.SatisfyBCRows(tmpH, HH);
	bcg.Solve(crsBC_C, tmpC, CC, iter, resid);
	crsBC_C.SatisfyBCRows(tmpC, CC);

	//  The substrate beyond an edge, as in 1-D: mirrored, or taken from the
	//  other end across a periodic pair.
	BoundaryCondition2D bcS; bcS.Build(nx, ny, SubstrateBC(_tXmin), SubstrateBC(_tXmax),
	                                   SubstrateBC(_tYmin), SubstrateBC(_tYmax));
	bcS.Apply(crsBC_S, tmpS, BoundaryCondition::Height, field);
	bcg.Solve(crsBC_S, tmpS, SS, iter, resid);
	crsBC_S.SatisfyBCRows(tmpS, SS);
	for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf)
	{
		//  The same as in the 1-D case. The ghost cells of the initial field are
		//  filled by solving a linear system. It costs more, but it is only
		//  reached once, at initialization. See the note at the same place above.
		field.HH(zf) = HH[G->IL(zf)];
		field.CC(zf) = CC[G->IL(zf)];
		field.SS(zf) = SS[G->IL(zf)];
	}
}

void BoundaryCondition2D::Apply(CRS_Matrix &crs, Vector<double>& rhs, tField tf, Field_WM2D& ini_Field)
{
	C_Grid2D *G = ini_Field.GetGridPtr();

	//------------------------------------------------------------------
	//  Build the EdgeSlot of all four edges first.
	//  Only three things differ from edge to edge: the row numbers of the ghost
	//  layers, the first column of the stencil, and the column stride. The
	//  boundary values for FreshInlet are filled in here too.
	//------------------------------------------------------------------
	std::vector<EdgeSlot> Xm(_Ny), Xp(_Ny), Ym(_Nx), Yp(_Nx);

	Zface zf0, zf1;
	for(int jy = 2; jy < _Ny-2; jy++)
	{
		EdgeSlot &m = Xm[jy];
		EdgeSlot &p = Xp[jy];
		for(int k = 0; k < 4; k++)
		{
			m.L[k] = IL(k, jy);			m.opp[k] = IL(_Nx-1-k, jy);
			p.L[k] = IL(_Nx-1-k, jy);	p.opp[k] = IL(k, jy);
		}
		m.colStart = m.L[0];	m.step = 1;		// on the min side, from the outside inwards
		p.colStart = p.L[3];	p.step = 1;		// on the max side, from the inside outwards

		zf0 = G->At<Zface>(1,jy);
		zf1 = G->At<Zface>(2,jy);
		m.hb = 0.5 * (ini_Field.HH(zf0) + ini_Field.HH(zf1));
		m.cb = 0.5 * (ini_Field.CC(zf0) + ini_Field.CC(zf1));

		//  Layers 1 and 2 on the max side too (see the 1-D Apply above), as
		//  the corners already had them.
		zf0 = G->At<Zface>(_Nx-2,jy);
		zf1 = G->At<Zface>(_Nx-3,jy);
		p.hb = 0.5 * (ini_Field.HH(zf0) + ini_Field.HH(zf1));
		p.cb = 0.5 * (ini_Field.CC(zf0) + ini_Field.CC(zf1));
	}
	for(int ix = 2; ix < _Nx-2; ix++)
	{
		EdgeSlot &m = Ym[ix];
		EdgeSlot &p = Yp[ix];
		for(int k = 0; k < 4; k++)
		{
			m.L[k] = IL(ix, k);			m.opp[k] = IL(ix, _Ny-1-k);
			p.L[k] = IL(ix, _Ny-1-k);	p.opp[k] = IL(ix, k);
		}
		m.colStart = m.L[0];	m.step = _Nx;
		p.colStart = p.L[3];	p.step = _Nx;

		zf0 = G->At<Zface>(ix,1);
		zf1 = G->At<Zface>(ix,2);
		m.hb = 0.5 * (ini_Field.HH(zf0) + ini_Field.HH(zf1));
		m.cb = 0.5 * (ini_Field.CC(zf0) + ini_Field.CC(zf1));

		zf0 = G->At<Zface>(ix, _Ny-2);
		zf1 = G->At<Zface>(ix, _Ny-3);
		p.hb = 0.5 * (ini_Field.HH(zf0) + ini_Field.HH(zf1));
		p.cb = 0.5 * (ini_Field.CC(zf0) + ini_Field.CC(zf1));
	}

	//------------------------------------------------------------------
	//  Apply the same body to all four edges. The rows of different edges are
	//  disjoint, so parallelizing this is safe.
	//------------------------------------------------------------------
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_Ny))
#endif
	for(int jy = 2; jy < _Ny-2; jy++)
	{
		ApplyEdge(crs, rhs, tf, _tXmin, Xm[jy], Form2D);
		ApplyEdge(crs, rhs, tf, _tXmax, Xp[jy], Form2D);
	}
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_Nx))
#endif
	for(int ix = 2; ix < _Nx-2; ix++)
	{
		ApplyEdge(crs, rhs, tf, _tYmin, Ym[ix], Form2D);
		ApplyEdge(crs, rhs, tf, _tYmax, Yp[ix], Form2D);
	}

	ApplyCornersAll(crs, rhs, tf, ini_Field);
}

//-------------------------------------------
//  Set up the index tables of the corners.
//-------------------------------------------
void BoundaryCondition2D::SetCornerIndex()
{
	ci_mm.OIoi = IL(0,3);	ci_mm.OIii = IL(1,3);	ci_mm.IIoi = IL(2,3);	ci_mm.IIii = IL(3,3);
	ci_mm.OIoo = IL(0,2);	ci_mm.OIio = IL(1,2);	ci_mm.IIoo = IL(2,2);	ci_mm.IIio = IL(3,2);
	ci_mm.OOoi = IL(0,1);	ci_mm.OOii = IL(1,1);	ci_mm.IOoi = IL(2,1);	ci_mm.IOii = IL(3,1);
	ci_mm.OOoo = IL(0,0);	ci_mm.OOio = IL(1,0);	ci_mm.IOoo = IL(2,0);	ci_mm.IOio = IL(3,0);

	ci_mp.OOoo = IL(0,_Ny-1);	ci_mp.OOio = IL(1,_Ny-1);	ci_mp.IOoo = IL(2,_Ny-1);	ci_mp.IOio = IL(3,_Ny-1);
	ci_mp.OOoi = IL(0,_Ny-2);	ci_mp.OOii = IL(1,_Ny-2);	ci_mp.IOoi = IL(2,_Ny-2);	ci_mp.IOii = IL(3,_Ny-2);
	ci_mp.OIoo = IL(0,_Ny-3);	ci_mp.OIio = IL(1,_Ny-3);	ci_mp.IIoo = IL(2,_Ny-3);	ci_mp.IIio = IL(3,_Ny-3);
	ci_mp.OIoi = IL(0,_Ny-4);	ci_mp.OIii = IL(1,_Ny-4);	ci_mp.IIoi = IL(2,_Ny-4);	ci_mp.IIii = IL(3,_Ny-4);

	ci_pm.IIii = IL(_Nx-4,3);	ci_pm.IIoi = IL(_Nx-3,3);	ci_pm.OIii = IL(_Nx-2,3);	ci_pm.OIoi = IL(_Nx-1,3);
	ci_pm.IIio = IL(_Nx-4,2);	ci_pm.IIoo = IL(_Nx-3,2);	ci_pm.OIio = IL(_Nx-2,2);	ci_pm.OIoo = IL(_Nx-1,2);
	ci_pm.IOii = IL(_Nx-4,1);	ci_pm.IOoi = IL(_Nx-3,1);	ci_pm.OOii = IL(_Nx-2,1);	ci_pm.OOoi = IL(_Nx-1,1);
	ci_pm.IOio = IL(_Nx-4,0);	ci_pm.IOoo = IL(_Nx-3,0);	ci_pm.OOio = IL(_Nx-2,0);	ci_pm.OOoo = IL(_Nx-1,0);

	ci_pp.IOio = IL(_Nx-4,_Ny-1);	ci_pp.IOoo = IL(_Nx-3,_Ny-1);	ci_pp.OOio = IL(_Nx-2,_Ny-1);	ci_pp.OOoo = IL(_Nx-1,_Ny-1);
	ci_pp.IOii = IL(_Nx-4,_Ny-2);	ci_pp.IOoi = IL(_Nx-3,_Ny-2);	ci_pp.OOii = IL(_Nx-2,_Ny-2);	ci_pp.OOoi = IL(_Nx-1,_Ny-2);
	ci_pp.IIio = IL(_Nx-4,_Ny-3);	ci_pp.IIoo = IL(_Nx-3,_Ny-3);	ci_pp.OIio = IL(_Nx-2,_Ny-3);	ci_pp.OIoo = IL(_Nx-1,_Ny-3);
	ci_pp.IIii = IL(_Nx-4,_Ny-4);	ci_pp.IIoi = IL(_Nx-3,_Ny-4);	ci_pp.OIii = IL(_Nx-2,_Ny-4);	ci_pp.OIoi = IL(_Nx-1,_Ny-4);
}

void BoundaryCondition2D::ApplyCornersAll(CRS_Matrix &crs, Vector<double>& rhs, tField tf, Field_WM2D& ini_Field)
{
	C_Grid2D *G = ini_Field.GetGridPtr();

	double ValMM = 0, ValMP = 0, ValPM = 0, ValPP = 0;
	switch(tf)
	{
	case Height:
		ValMM = 0.25*(ini_Field.HH(G->At<Zface>(1,1)) + ini_Field.HH(G->At<Zface>(2,1)) + ini_Field.HH(G->At<Zface>(1,2)) + ini_Field.HH(G->At<Zface>(2,2)));
		ValPM = 0.25*(ini_Field.HH(G->At<Zface>(_Nx-2,1)) + ini_Field.HH(G->At<Zface>(_Nx-3,1)) + ini_Field.HH(G->At<Zface>(_Nx-2,2)) + ini_Field.HH(G->At<Zface>(_Nx-3,2)));
		ValMP = 0.25*(ini_Field.HH(G->At<Zface>(1,_Ny-2)) + ini_Field.HH(G->At<Zface>(2,_Ny-2)) + ini_Field.HH(G->At<Zface>(1,_Ny-3)) + ini_Field.HH(G->At<Zface>(2,_Ny-3)));
		ValPP = 0.25*(ini_Field.HH(G->At<Zface>(_Nx-2,_Ny-2)) + ini_Field.HH(G->At<Zface>(_Nx-3,_Ny-2)) + ini_Field.HH(G->At<Zface>(_Nx-2,_Ny-3)) + ini_Field.HH(G->At<Zface>(_Nx-3,_Ny-3)));
		break;
	case HeightJacobi:
		ValMM = ValMP = ValPM = ValPP = 0;
		break;
	case Conc:
		ValMM = 0.25*(ini_Field.CC(G->At<Zface>(1,1)) + ini_Field.CC(G->At<Zface>(2,1)) + ini_Field.CC(G->At<Zface>(1,2)) + ini_Field.CC(G->At<Zface>(2,2)));
		ValPM = 0.25*(ini_Field.CC(G->At<Zface>(_Nx-2,1)) + ini_Field.CC(G->At<Zface>(_Nx-3,1)) + ini_Field.CC(G->At<Zface>(_Nx-2,2)) + ini_Field.CC(G->At<Zface>(_Nx-3,2)));
		ValMP = 0.25*(ini_Field.CC(G->At<Zface>(1,_Ny-2)) + ini_Field.CC(G->At<Zface>(2,_Ny-2)) + ini_Field.CC(G->At<Zface>(1,_Ny-3)) + ini_Field.CC(G->At<Zface>(2,_Ny-3)));
		ValPP = 0.25*(ini_Field.CC(G->At<Zface>(_Nx-2,_Ny-2)) + ini_Field.CC(G->At<Zface>(_Nx-3,_Ny-2)) + ini_Field.CC(G->At<Zface>(_Nx-2,_Ny-3)) + ini_Field.CC(G->At<Zface>(_Nx-3,_Ny-3)));
		break;
	}
	// Xmin-Ymin
	ApplyCorners(crs,rhs, tf, ValMM, _tXmin, _tYmin, ci_mm, ci_pm, ci_mp, ci_pp);
	// Xmin-Ymax
	ApplyCorners(crs,rhs, tf, ValMP, _tXmin, _tYmax, ci_mp, ci_pp, ci_mm, ci_pm);
	// Xmax-Ymin
	ApplyCorners(crs,rhs, tf, ValPM, _tXmax, _tYmin, ci_pm, ci_mm, ci_pp, ci_mp);
	// Xmax-Ymax
	ApplyCorners(crs,rhs, tf, ValPP, _tXmax, _tYmax, ci_pp, ci_mp, ci_pm, ci_mm);
}

//==================================================================
//  Shared parts of the corner handling
//==================================================================

//  Make each of the four ghost rows a copy of the point given for it.
void BoundaryCondition2D::CopyCorner(CRS_Matrix &crs, Vector<double>& rhs, const stCornerIndex &ci,
                                     int s_ii, int s_io, int s_oi, int s_oo)
{
	EqualVal(crs, rhs, ci.OOii, s_ii);
	EqualVal(crs, rhs, ci.OOio, s_io);
	EqualVal(crs, rhs, ci.OOoi, s_oi);
	EqualVal(crs, rhs, ci.OOoo, s_oo);
}

//  Place a 5-point stencil on the four ghost rows.
//  The order of the rows, and the order in which the columns within a row are
//  inserted, are fixed (diagonal -> IIii -> IIio -> IIoi -> IIoo).
//  Reordering them makes the results no longer agree bit for bit.
void BoundaryCondition2D::CornerStencil(CRS_Matrix &crs, Vector<double>& rhs, const stCornerIndex &ci,
                                        const stCornerRow row[4], double cornerVal, bool setRhs)
{
	const int r[4] = { ci.OOoo, ci.OOio, ci.OOoi, ci.OOii };
	for(int k = 0; k < 4; k++)
	{
		crs.Diag(r[k])     = row[k].d;
		crs(r[k], ci.IIii) = row[k].cII;
		crs(r[k], ci.IIio) = row[k].cIo;
		crs(r[k], ci.IIoi) = row[k].coI;
		crs(r[k], ci.IIoo) = row[k].coo;
		if(setRhs) rhs[r[k]] = row[k].r * cornerVal;
	}
}

//------------------------------------------------------------------
//  Coefficient tables of the stencils. The rows come in the order
//  OOoo, OOio, OOoi, OOii.
//  The coefficients on the right-hand side are written out unfolded, so that
//  the shape of the formula stays visible.
//------------------------------------------------------------------
namespace {

typedef stCornerRow Row;

//  Inlet-Inlet : height (HeightJacobi uses the same coefficients, with a right-hand side of zero)
const Row kInletInletH[4] = {
	{ -1.0, +4.0, -18.0, -18.0, +81.0, -(16.0-72.0+ 8.0) },
	{ -3.0, +2.0,  -9.0, -12.0, +54.0, -( 8.0-48.0+ 8.0) },
	{ -3.0, +2.0, -12.0,  -9.0, +54.0, -(16.0-72.0+24.0) },
	{ -9.0, +1.0,  -6.0,  -6.0, +36.0, -( 8.0-48.0+24.0) }
};
//  Inlet-Inlet : concentration
const Row kInletInletC[4] = {
	{ -49.0, +4.0, -54.0, -54.0, +729.0, -( 64.0-864.0+224.0) },
	{ -49.0, +2.0, -27.0, -20.0, +270.0, -( 32.0-320.0+112.0) },
	{ -49.0, +2.0, -20.0, -27.0, +270.0, -( 32.0-432.0+224.0) },
	{ -49.0, +1.0, -10.0, -10.0, +100.0, -( 16.0-160.0+112.0) }
};
//  Inlet-Outlet
const Row kInletOutlet[4] = {
	{ -1.0, -4.0,  +6.0, +18.0, -27.0, -8.0 },
	{ -3.0, -2.0,  +3.0, +12.0, -18.0, -8.0 },
	{ -1.0, -2.0,  +4.0,  +9.0, -18.0, -8.0 },
	{ -3.0, -1.0,  +2.0,  +6.0, -12.0, -8.0 }
};
//  Outlet-Inlet
const Row kOutletInlet[4] = {
	{ -1.0, -4.0, +18.0,  +6.0, -27.0, -(-16.0+24.0) },
	{ -1.0, -2.0,  +9.0,  +4.0, -18.0, -( -8.0+16.0) },
	{ -3.0, -2.0, +12.0,  +3.0, -18.0, -(-16.0+24.0) },
	{ -3.0, -1.0,  +6.0,  +2.0, -12.0, -( -8.0+16.0) }
};
//  Outlet-Outlet (the right-hand side is not set)
const Row kOutletOutlet[4] = {
	{ -1.0, +4.0, -6.0, -6.0, +9.0, 0.0 },
	{ -1.0, +2.0, -3.0, -4.0, +6.0, 0.0 },
	{ -1.0, +2.0, -4.0, -3.0, +6.0, 0.0 },
	{ -1.0, +1.0, -2.0, -2.0, +4.0, 0.0 }
};

} // namespace


//==================================================================
//  The boundary condition for one corner.
//
//  A 5x5 table of the combinations of the X-direction and Y-direction
//  conditions. There are only two things to be done:
//    - make the rows copies of four points somewhere (Mirror*)
//    - place a 5-point stencil                       (CornerStencil)
//
//  How the copy is taken is decided by whether the X direction is folded back
//  and whether the Y direction is. Symmetry folds back on its own side,
//  Periodic on the opposite side. Pinning / FreshInlet / Outlet do not fold
//  back in that direction at all.
//==================================================================
void BoundaryCondition2D::ApplyCorners(CRS_Matrix &crs, Vector<double>& rhs, tField t_field, double cornerVal,
		tBoundaryCondition tbcX, tBoundaryCondition tbcY,
		const stCornerIndex &ci, const stCornerIndex &ciRO, const stCornerIndex &ciOR, const stCornerIndex &ciRR)
{
	/********************************
	 tbcX: boundary condition in the X direction
	 tbcY: boundary condition in the Y direction
	 ci:   index table of the corner in question
	 ciRO: opposite side in the X direction
	 ciOR: opposite side in the Y direction
	 ciRR: opposite side in both directions
	********************************/
	switch(tbcX)
	{
	default:
	case Symmetry:
		switch(tbcY)
		{
		default:
		case Symmetry:		MirrorXY(crs,rhs, ci, ci);		break;	//Symmetry-Symmetry
		case Periodic:		MirrorXY(crs,rhs, ci, ciOR);	break;	//Symmetry-Periodic
		case Pinning:												//Symmetry-Pinning
		case FreshInlet:											//Symmetry-Inlet
		case Outlet:		MirrorX (crs,rhs, ci, ci);		break;	//Symmetry-Outlet
		}
		break;

	case Periodic:
		switch(tbcY)
		{
		default:
		case Symmetry:		MirrorXY(crs,rhs, ci, ciRO);	break;	//Periodic-Symmetry
		case Periodic:		MirrorXY(crs,rhs, ci, ciRR);	break;	//Periodic-Periodic
		case Pinning:												//Periodic-Pinning
		case FreshInlet:											//Periodic-Inlet
		case Outlet:		MirrorX (crs,rhs, ci, ciRO);	break;	//Periodic-Outlet
		}
		break;

	case Pinning:
		switch(tbcY)
		{
		default:
		case Symmetry:		//Pinning-Symmetry
			switch(t_field)
			{
			case Height:
			case HeightJacobi:
				crs.Diag(ci.OOii) = 0.5; crs(ci.OOii, ci.IIoo) =  0.5;
				crs.Diag(ci.OOio) = 0.5; crs(ci.OOio, ci.IIoi) =  0.5;
				crs.Diag(ci.OOoi) = 1.0; crs(ci.OOoi, ci.OIoo) = -1.0;
				crs.Diag(ci.OOoo) = 1.0; crs(ci.OOoo, ci.OIoi) = -1.0;
				break;
			case Conc:
				//  On an edge, Pinning-Conc copies inwards, but at a corner the mirror image is taken in both directions.
				MirrorXY(crs,rhs, ci, ci);
				break;
			}
			break;

		case Periodic:		//Pinning-Periodic
			MirrorY(crs,rhs, ci, ciOR);
			break;

		case Pinning:		//Pinning-Pinning
			switch(t_field)
			{
			case Height:
			case HeightJacobi:
				crs.Diag(ci.OOii) = 0.25; crs(ci.OOii, ci.IIoo) = 0.25; crs(ci.OOii,ci.OIio) = 0.25; crs(ci.OOii,ci.IOoi) = 0.25;
				break;
			case Conc:
				CopyCorner(crs,rhs, ci, ci.IIoo, ci.OOii, ci.OOii, ci.OOii);
				break;
			}
			break;

		//  Combinations in which Pinning meets FreshInlet or Outlet at a 2-D
		//  corner. Which of the two conditions to impose on the single corner
		//  point is not determined (Pinning holds the height fixed, FreshInlet
		//  prescribes the height and the concentration, and Outlet sets the
		//  third derivative to zero), so these are left unimplemented and abort
		//  the run. Quietly picking one of the two would leave no way of
		//  noticing that the boundary is not doing its job.
		//  If they are ever needed, a condition specific to the corner has to be
		//  decided on and written in the same form as CopyCorner.
		case FreshInlet:	//Pinning-Inlet
			FatalErrorIn("BC: Pinning-Inlet: not yet implemented.");
			Global::FatalError.Abort();
			break;
		case Outlet:		//Pinning-Outlet
			FatalErrorIn("BC: Pinning-Outlet: not yet implemented.");
			Global::FatalError.Abort();
			break;
		}
		break;

	case FreshInlet:
		switch(tbcY)
		{
		default:
		case Symmetry:		MirrorY(crs,rhs, ci, ci);		break;	//Inlet-Symmetry
		case Periodic:		MirrorY(crs,rhs, ci, ciOR);		break;	//Inlet-Periodic

		case Pinning:		//Inlet-Pinning
			FatalErrorIn("BC: Inlet-Pinning: not yet implemented.");
			Global::FatalError.Abort();
			break;

		case FreshInlet:	//Inlet-Inlet
			switch(t_field)
			{
			//  HeightJacobi uses the same coefficients as Height, with the
			//  right-hand side of a correction that keeps the corner as it
			//  is: zero (cornerVal is 0 for HeightJacobi), as the Inlet-Outlet
			//  corners and the FreshInlet edges have it. Left untouched, it
			//  would keep the Height value r*h_corner written there first, and
			//  every Newton iteration would correct the corner by that much
			//  again: the iteration stagnates whatever the step, the step is
			//  cut down to MinDt and the run aborts (a film with inlets on two
			//  adjacent edges, tools/proptest.py seed 1171).
			case Height:		CornerStencil(crs,rhs, ci, kInletInletH, cornerVal, true);	break;
			case HeightJacobi:	CornerStencil(crs,rhs, ci, kInletInletH, cornerVal, true);	break;
			case Conc:			CornerStencil(crs,rhs, ci, kInletInletC, cornerVal, true);	break;
			}
			break;

		case Outlet:		//Inlet-Outlet
			CornerStencil(crs,rhs, ci, kInletOutlet, cornerVal, true);
			break;
		}
		break;

	case Outlet:
		switch(tbcY)
		{
		default:
		case Symmetry:		MirrorY(crs,rhs, ci, ci);		break;	//Outlet-Symmetry
		case Periodic:		MirrorY(crs,rhs, ci, ciOR);		break;	//Outlet-Periodic

		case Pinning:		//Outlet-Pinning
			FatalErrorIn("BC: Outlet-Pinning: not yet implemented.");
			Global::FatalError.Abort();
			break;

		case FreshInlet:	//Outlet-Inlet
			CornerStencil(crs,rhs, ci, kOutletInlet, cornerVal, true);
			break;

		case Outlet:		//Outlet-Outlet (the right-hand side is not set)
			CornerStencil(crs,rhs, ci, kOutletOutlet, cornerVal, false);
			break;
		}
		break;
	}
}
