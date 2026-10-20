/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_jacobian.cpp : each term's Jacobian against a finite difference of
*                      the term itself
*
*  The Newton iteration of the film-height equation solves with
*
*      J(h) = d/dh [ A(h) h ],
*
*  assembled as A(h) (SetMatrixA) plus what each term adds for the
*  dependence of its coefficients on h (SetJacobiMatrix). A mistake in that
*  second part does not change what the iteration converges to, only how
*  fast, so none of the solver tests can see it: they converge anyway, a
*  little slower. It shows up instead as a time step that collapses on a
*  hard case. Here the product J(h) v is compared with the central
*  difference [R(h + e v) - R(h - e v)] / 2e of R(h) = A(h) h, term by
*  term, in 1-D and 2-D, on a film that is not flat, over a substrate that is
*  not flat, with a concentration that varies -- so that every coefficient
*  of every term depends on position.
*
*  The kind of error it catches: the derivative of the face mobility taken
*  at the neighbor for the center column too (vertical gravity, disjoining
*  pressure); all of Pi' added to a Jacobian that already holds Pi/h from A
*  (disjoining pressure); a factor h missing from the sideways outflow of
*  the 1-D centrifugal term.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "doctest/doctest.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "NonDimensionalNumbers.h"
#include "LubFields.h"
#include "BoundaryCondition.h"
#include "CRS_Matrix.h"
#include "WellMixedTraits.h"
#include "FiniteDifferences.h"

using namespace yuzen;
using namespace yuzen::grid;

namespace {

const double PI = 3.14159265358979323846;

enum Term { tLaplace, tVertical, tParallel, tCentrifugal, tSoluto, tDisjoining };

const char *Name(Term t)
{
	switch(t)
	{
	case tLaplace:     return "capillary pressure";
	case tVertical:    return "vertical gravity";
	case tParallel:    return "in-plane gravity";
	case tCentrifugal: return "centrifugal";
	case tSoluto:      return "Marangoni";
	default:           return "disjoining pressure";
	}
}

//  The Navier slip length [m] MakeInput gives (Physics.Slip); 0: none. Set
//  by the slip test cases through SlipScope, for the duration of one case.
double gSlipLength = 0.0;
struct SlipScope
{
	explicit SlipScope(double len) { gSlipLength = len; }
	~SlipScope() { gSlipLength = 0.0; }
};

//  One term on; the fluid properties depend on the concentration.
Input MakeInput(Term t)
{
	Input inp;
	//  Input() leaves the flags unset (the parser sets every one), and a
	//  copy of an indeterminate bool is undefined (UBSan stops on it):
	//  clear those this test does not set itself.
	inp.Physics.ConcField.EnableConvection = inp.Physics.ConcField.EnableDiffusion = false;
	inp.Physics.Source.Enable = false;
	inp.Output.Overwrite = inp.Output.Dimensioned = false;
	inp.Output.FieldData.NeutralFile = inp.Output.FieldData.VTKFile = false;
	inp.Numerics.TimeStepControl.Enable = false;
	inp.Numerics.HeightEquation.LinearSolver.RowScaling = false;
	inp.Numerics.ConcEquation.LinearSolver.RowScaling = false;
	inp.Numerics.ManufacturedSolution.Enable = false;
	inp.Restart.Enable = inp.Restart.ASCII = false;
	inp.LengthScaleType = "capillary";
	inp.LengthScale     = 0;

	inp.Fluid.Density.Value            = 1000.0;
	inp.Fluid.SurfaceTension.BaseValue = 30.0e-3;
	inp.Fluid.SurfaceTension.BaseConc  = 0.5;
	inp.Fluid.SurfaceTension.Coef      = 10.0e-3;
	inp.Fluid.Viscosity.BaseValue      = 100.0e-3;
	inp.Fluid.Viscosity.BaseConc       = 0.5;
	inp.Fluid.Viscosity.Exponent       = 3.0;
	inp.Fluid.Diffusivity.BaseValue    = 1.0e-9;
	inp.Fluid.Diffusivity.BaseConc     = 0.5;
	inp.Fluid.Diffusivity.Exponent     = 0.0;
	inp.Fluid.MolarMass.Solvent        = 100.0;
	inp.Fluid.MolarMass.Resin          = 2000.0;

	inp.Initial.SurfaceShape.TotalHeight    = 10e-6;
	inp.Initial.SurfaceShape.PrecursorRatio = 0.1;

	inp.Physics.ConcField.Model            = "WM";
	inp.Physics.LaplacePressure.Enable     = (t == tLaplace);
	inp.Physics.SolutoCapillarity.Enable   = (t == tSoluto);
	inp.Physics.VerticalGravity.Enable     = (t == tVertical);
	inp.Physics.VerticalGravity.Acceleration     = 9.80665;
	inp.Physics.HorizontalGravity.Enable   = (t == tParallel);
	inp.Physics.HorizontalGravity.Acceleration.assign(2, 0.0);
	inp.Physics.HorizontalGravity.Acceleration[0] = 4.903;
	inp.Physics.HorizontalGravity.Acceleration[1] = -2.0;
	inp.Physics.Centrifugal.Enable         = (t == tCentrifugal);
	inp.Physics.Centrifugal.RotationRate   = 40.0;
	inp.Physics.Centrifugal.RotationAxis.assign(2, 0.0);
	inp.Physics.Centrifugal.RotationAxis[0] = 0.13e-3;
	inp.Physics.Centrifugal.RotationAxis[1] = -0.07e-3;
	inp.Physics.Evaporation.Enable         = false;
	inp.Physics.Evaporation.Coef           = 0.0;
	inp.Physics.DisjoiningPressure.Enable       = (t == tDisjoining);
	inp.Physics.DisjoiningPressure.ContactAngle = 10.0;
	inp.Physics.DisjoiningPressure.ExponentN    = 3.0;
	inp.Physics.DisjoiningPressure.ExponentM    = 2.0;
	inp.Physics.Slip.Enable = gSlipLength > 0;
	inp.Physics.Slip.Length = gSlipLength;
	return inp;
}

Input::I_Grid MakeGrid(bool twoD, int nx, int ny, double scale)
{
	Input::I_Grid g;
	g.Dimension = twoD ? "xyh" : "xh";
	g.Xmin = -0.5e-3 / scale; g.Xmax = 0.5e-3 / scale;
	g.Ymin = -0.4e-3 / scale; g.Ymax = 0.4e-3 / scale;
	if(!twoD) { g.Ymin = -1.0; g.Ymax = 1.0; }
	g.Nx = nx; g.Ny = twoD ? ny : 1; g.Nz = 1;
	return g;
}

Input::I_Substrate FlatSubstrate()
{
	Input::I_Substrate s;
	s.TotalHeight = 0.0;
	return s;
}

//  Smooth fields over the whole grid, ghosts included. x and y are
//  dimensionless; X and Y are the half-widths of the domain.
struct Shapes
{
	double X, Y, hmin;
	//  front: a film rising through a steep front in x and in y instead of
	//  the waves, for the limited schemes. It has no extremum, where the
	//  limiter has a kink at a zero difference (the ramp keeps every
	//  difference above 0.01), and its fronts are steep enough that the
	//  ratio r of successive differences runs from below 0.4 to above 4
	//  through every branch of Koren's limiter, on the grids used here
	//  (28 cells in 1-D, 18 x 14 in 2-D), while no r comes within 0.05 of
	//  a kink (0.4 or 4).
	bool front;
	double H(double x, double y) const
	{
		if(front)
			return hmin + 0.3 * (1.0 + 0.4 * x / X + 0.9 * std::tanh(16.5 * x / X + 0.3))
			                  * (1.0 + 0.4 * y / Y + 0.9 * std::tanh(9.0 * y / Y + 0.3));
		return hmin + 0.6 * (1.0 + std::sin(PI * x / X + 0.3)) * (1.1 + 0.4 * std::cos(PI * y / Y - 0.2));
	}
	double S(double x, double y) const { return 0.3 * std::cos(0.7 * PI * x / X) * std::sin(0.5 * PI * y / Y + 0.4); }
	double C(double x, double y) const { return 0.5 + 0.15 * std::sin(0.8 * PI * x / X - 0.5) + 0.1 * std::cos(0.6 * PI * y / Y); }
	//  The direction of the perturbation: anything that is not smooth
	double V(int k) const { return std::sin(12.9898 * k + 78.233) * 0.5 + std::cos(3.1 * k * k); }
};

//------------------------------------------------------------------
//  The comparison, for one term. Grid, field and term are those of the
//  1-D or the 2-D model (the traits Tr); returns the largest error of
//  J v against the central difference, relative to the largest entry of
//  the central difference, over the rows of the interior.
//------------------------------------------------------------------
template<class Tr, class FieldT, class TermT>
double JacobianError(TermT &term, FieldT &field, const Shapes &sh, double eps)
{
	typedef typename Tr::Grid Grid;
	typedef typename Tr::Site Site;
	Grid *G = field.GetGridPtr();

	int n = 0;
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
		n = std::max(n, Tr::Index(G, s) + 1);

	std::vector<double> h0(n, 0.0), v(n, 0.0);
	std::vector<char> inside(n, 0);
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
	{
		double x, y;
		Tr::SiteXY(G, s, x, y);
		const int i = Tr::Index(G, s);
		h0[i] = sh.H(x, y);
		field.HH(s) = h0[i];
		field.SS(s) = sh.S(x, y);
		field.CC(s) = sh.C(x, y);
		inside[i] = (G->CheckGhost(s) == Grid::Inside);
		if(inside[i]) v[i] = sh.V(i);
	}

	auto setH = [&](double e)
	{
		for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
		{
			const int i = Tr::Index(G, s);
			field.HH(s) = h0[i] + e * v[i];
		}
	};
	auto assembleA = [&](CRS_Matrix &A)
	{
		A.Build(n);
		for(int i = 0; i < n; ++i) { A(i, i) = 0.0; if(!inside[i]) A.SetAsBC(i); }
		term.SetMatrixA(A, field);
	};
	auto residual = [&](double e, Vector<double> &R)
	{
		setH(e);
		CRS_Matrix A;
		assembleA(A);
		Vector<double> h(n);
		for(int i = 0; i < n; ++i) h[i] = h0[i] + e * v[i];
		R.resize(n);
		A.MatVec(R, h);
	};

	setH(0.0);
	CRS_Matrix J;
	assembleA(J);
	term.SetJacobiMatrix(J, field);
	Vector<double> vv(n), Jv(n);
	for(int i = 0; i < n; ++i) vv[i] = v[i];
	J.MatVec(Jv, vv);

	Vector<double> Rp(n), Rm(n);
	residual(+eps, Rp);
	residual(-eps, Rm);
	setH(0.0);

	double worst = 0, scale = 0;
	for(int i = 0; i < n; ++i)
	{
		if(!inside[i]) continue;
		const double fd = (Rp[i] - Rm[i]) / (2.0 * eps);
		scale = std::max(scale, std::fabs(fd));
		worst = std::max(worst, std::fabs(Jv[i] - fd));
	}
	return scale > 0 ? worst / scale : 0.0;
}

//------------------------------------------------------------------
//  What the term puts through the edges: |sum of A(h) h over the interior
//  rows| relative to sum |A(h) h| (the grids are uniform, so the sum is
//  the volume balance up to the cell measure). The ghosts hold the smooth
//  shapes, not what the boundary condition would put there: on a wall the
//  flux of an upwinded term must be zero whatever they hold.
//------------------------------------------------------------------
template<class Tr, class FieldT, class TermT>
double WallLeak(TermT &term, FieldT &field, const Shapes &sh)
{
	typedef typename Tr::Grid Grid;
	typedef typename Tr::Site Site;
	Grid *G = field.GetGridPtr();
	int n = 0;
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
		n = std::max(n, Tr::Index(G, s) + 1);
	std::vector<char> inside(n, 0);
	Vector<double> h(n);
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
	{
		double x, y;
		Tr::SiteXY(G, s, x, y);
		const int i = Tr::Index(G, s);
		h[i] = field.HH(s) = sh.H(x, y);
		field.SS(s) = sh.S(x, y);
		field.CC(s) = sh.C(x, y);
		inside[i] = (G->CheckGhost(s) == Grid::Inside);
	}
	CRS_Matrix A;
	A.Build(n);
	for(int i = 0; i < n; ++i) { A(i, i) = 0.0; if(!inside[i]) A.SetAsBC(i); }
	term.SetMatrixA(A, field);
	Vector<double> R(n);
	A.MatVec(R, h);
	double sum = 0, mag = 0;
	for(int i = 0; i < n; ++i)
		if(inside[i]) { sum += R[i]; mag += std::fabs(R[i]); }
	return mag > 0 ? std::fabs(sum) / mag : 0.0;
}

//------------------------------------------------------------------
//  The conservative form of A h (PhysicsTermBase::AddConservative) of the
//  capillary pressure and the vertical gravity, against the product of
//  their entries of A with h. Returns, over the interior rows,
//    agree: max |A h - AddConservative| / max_i sum_j |A_ij h_j|
//    sum:   |sum of AddConservative| / sum of |AddConservative|
//  With mirror = true the ghosts of h, s and c are the mirror image of the
//  interior (what Symmetry edges give), so the box is closed: every edge
//  face then carries exactly no flux and each interior face's flux is taken
//  out of one cell and put into the other, so sum is at rounding of the
//  fluxes, not of the entries of A times h.
//------------------------------------------------------------------
template<class Tr, class FieldT, class TermT>
void ConservativeForm(TermT &term, FieldT &field, const Shapes &sh, bool mirror,
                      double &agree, double &sum)
{
	typedef typename Tr::Grid Grid;
	typedef typename Tr::Site Site;
	Grid *G = field.GetGridPtr();
	int n = 0;
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
		n = std::max(n, Tr::Index(G, s) + 1);

	//  The rows are ix + nx iy (1-D: ix). The mirror of an index along one
	//  axis: 0 <-> 3, 1 <-> 2 at either end.
	const int nx = G->nx(), ny = n / nx;
	auto mir = [](int i, int nn) { return i < 2 ? 3 - i : (i >= nn - 2 ? 2 * nn - 5 - i : i); };
	std::vector<char> inside(n, 0);
	std::vector<double> hv(n), sv(n), cv(n);
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
	{
		const int i = Tr::Index(G, s);
		inside[i] = (G->CheckGhost(s) == Grid::Inside);
		double x, y;
		Tr::SiteXY(G, s, x, y);
		hv[i] = sh.H(x, y); sv[i] = sh.S(x, y); cv[i] = sh.C(x, y);
	}
	Vector<double> h(n);
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
	{
		const int i = Tr::Index(G, s);
		const int j = mirror ? mir(i % nx, nx) + nx * (ny > 1 ? mir(i / nx, ny) : 0) : i;
		h[i] = field.HH(s) = hv[j];
		field.SS(s) = sv[j];
		field.CC(s) = cv[j];
	}
	CRS_Matrix A;
	A.Build(n);
	for(int i = 0; i < n; ++i) { A(i, i) = 0.0; if(!inside[i]) A.SetAsBC(i); }
	term.SetMatrixA(A, field);
	REQUIRE(term.ConservativeForm());

	Vector<double> Ah(n), y(n);
	A.MatVec(Ah, h);
	y = 0;
	term.AddConservative(y, h, A, field);

	double worst = 0, scale = 0, s = 0, m = 0;
	bool ghostsUntouched = true;
	for(int i = 0; i < n; ++i)
	{
		if(!inside[i]) { if(y[i] != 0.0) ghostsUntouched = false; continue; }
		double row = 0;
		for(int k = 0; k < A.NonZero(i); ++k) row += std::fabs(A.Element(i, k) * h[A.Jcol(i, k)]);
		scale = std::max(scale, row);
		worst = std::max(worst, std::fabs(Ah[i] - y[i]));
		s += y[i]; m += std::fabs(y[i]);
	}
	CHECK(ghostsUntouched);	// the boundary rows are left to the conditions
	agree = scale > 0 ? worst / scale : 0.0;
	sum   = m > 0 ? std::fabs(s) / m : 0.0;
}

//  ConservativeForm of the capillary pressure (vertical = false) or the
//  vertical gravity, in 1-D or 2-D, with edges of kind e
void Conservative(bool twoD, bool vertical, BoundaryCondition::tBoundaryCondition e, bool mirror,
                  double &agree, double &sum)
{
	const Term t = vertical ? tVertical : tLaplace;
	Input inp = MakeInput(t);
	NonDimensionalNumbers nd;
	nd.Build(inp, inp.Initial.SurfaceShape.TotalHeight);
	Input::I_Grid g = MakeGrid(twoD, twoD ? 18 : 28, 14, nd.ExtentScale());
	Shapes sh = {g.Xmax, twoD ? g.Ymax : 1.0, 0.3, false};
	if(twoD)
	{
		Field_WM2D field;
		field.Build(g, FlatSubstrate(), nd.ExtentScale());
		BoundaryCondition2D bc(field.GetGridPtr()->nx(), field.GetGridPtr()->ny(), e, e, e, e);
		if(vertical) { WM2D_Traits::VerticalGravity k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); ConservativeForm<WM2D_Grid>(k, field, sh, mirror, agree, sum); }
		else         { WM2D_Traits::LaplacePressure k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); ConservativeForm<WM2D_Grid>(k, field, sh, mirror, agree, sum); }
		return;
	}
	Field_WM1D field;
	field.Build(g, FlatSubstrate(), nd.ExtentScale());
	BoundaryCondition1D bc(field.GetGridPtr()->nx(), e, e);
	if(vertical) { WM1D_Traits::VerticalGravity k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); ConservativeForm<WM1D_Grid>(k, field, sh, mirror, agree, sum); }
	else         { WM1D_Traits::LaplacePressure k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); ConservativeForm<WM1D_Grid>(k, field, sh, mirror, agree, sum); }
}

//  The 1-D and the 2-D set-up around JacobianError
//  scheme: the face interpolation of the three upwinded terms (in-plane
//  gravity, centrifugal, Marangoni; Numerics.Convection.Height); the other
//  terms have none.
typedef C_Interpolation_base4 IB;
bool LimitedScheme(IB::t_Scheme t) { IB b; b.Build(t); return b.IsLimited(); }

double Error1D(Term t, bool pinned = false, double eps = 1e-6, IB::t_Scheme scheme = IB::KK)
{
	Input inp = MakeInput(t);
	NonDimensionalNumbers nd;
	nd.Build(inp, inp.Initial.SurfaceShape.TotalHeight);
	Input::I_Grid g = MakeGrid(false, 28, 1, nd.ExtentScale());

	Field_WM1D field;
	field.Build(g, FlatSubstrate(), nd.ExtentScale());
	const BoundaryCondition::tBoundaryCondition e = pinned ? BoundaryCondition::Pinning : BoundaryCondition::Symmetry;
	BoundaryCondition1D bc(field.GetGridPtr()->nx(), e, e);
	Shapes sh = {g.Xmax, 1.0, t == tDisjoining ? 0.12 : 0.3, LimitedScheme(scheme)};

	switch(t)
	{
	case tLaplace:     { WM1D_Traits::LaplacePressure   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); return JacobianError<WM1D_Grid>(k, field, sh, eps); }
	case tVertical:    { WM1D_Traits::VerticalGravity   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); return JacobianError<WM1D_Grid>(k, field, sh, eps); }
	case tParallel:    { WM1D_Traits::ParallelGravity   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return JacobianError<WM1D_Grid>(k, field, sh, eps); }
	case tCentrifugal: { WM1D_Traits::Centrifugal       k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return JacobianError<WM1D_Grid>(k, field, sh, eps); }
	case tSoluto:      { WM1D_Traits::SolutoCapillarity k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return JacobianError<WM1D_Grid>(k, field, sh, eps); }
	default:           { WM1D_Traits::DisjoiningPressure k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); return JacobianError<WM1D_Grid>(k, field, sh, eps); }
	}
}

double Error2D(Term t, bool pinned = false, double eps = 1e-6, IB::t_Scheme scheme = IB::KK)
{
	Input inp = MakeInput(t);
	NonDimensionalNumbers nd;
	nd.Build(inp, inp.Initial.SurfaceShape.TotalHeight);
	Input::I_Grid g = MakeGrid(true, 18, 14, nd.ExtentScale());

	Field_WM2D field;
	field.Build(g, FlatSubstrate(), nd.ExtentScale());
	const BoundaryCondition::tBoundaryCondition e = pinned ? BoundaryCondition::Pinning : BoundaryCondition::Symmetry;
	BoundaryCondition2D bc(field.GetGridPtr()->nx(), field.GetGridPtr()->ny(), e, e, e, e);
	Shapes sh = {g.Xmax, g.Ymax, t == tDisjoining ? 0.12 : 0.3, LimitedScheme(scheme)};

	switch(t)
	{
	case tLaplace:     { WM2D_Traits::LaplacePressure   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); return JacobianError<WM2D_Grid>(k, field, sh, eps); }
	case tVertical:    { WM2D_Traits::VerticalGravity   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); return JacobianError<WM2D_Grid>(k, field, sh, eps); }
	case tParallel:    { WM2D_Traits::ParallelGravity   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return JacobianError<WM2D_Grid>(k, field, sh, eps); }
	case tCentrifugal: { WM2D_Traits::Centrifugal       k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return JacobianError<WM2D_Grid>(k, field, sh, eps); }
	case tSoluto:      { WM2D_Traits::SolutoCapillarity k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return JacobianError<WM2D_Grid>(k, field, sh, eps); }
	default:           { WM2D_Traits::DisjoiningPressure k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); return JacobianError<WM2D_Grid>(k, field, sh, eps); }
	}
}

//  WallLeak of an upwinded term (in-plane gravity, centrifugal, Marangoni)
//  on a 1-D or a 2-D grid whose four edges are of kind e
double Leak(bool twoD, Term t, BoundaryCondition::tBoundaryCondition e, IB::t_Scheme scheme)
{
	Input inp = MakeInput(t);
	NonDimensionalNumbers nd;
	nd.Build(inp, inp.Initial.SurfaceShape.TotalHeight);
	Input::I_Grid g = MakeGrid(twoD, twoD ? 18 : 28, 14, nd.ExtentScale());
	Shapes sh = {g.Xmax, twoD ? g.Ymax : 1.0, 0.3, LimitedScheme(scheme)};
	if(twoD)
	{
		Field_WM2D field;
		field.Build(g, FlatSubstrate(), nd.ExtentScale());
		BoundaryCondition2D bc(field.GetGridPtr()->nx(), field.GetGridPtr()->ny(), e, e, e, e);
		switch(t)
		{
		case tParallel:    { WM2D_Traits::ParallelGravity   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return WallLeak<WM2D_Grid>(k, field, sh); }
		case tCentrifugal: { WM2D_Traits::Centrifugal       k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return WallLeak<WM2D_Grid>(k, field, sh); }
		default:           { WM2D_Traits::SolutoCapillarity k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return WallLeak<WM2D_Grid>(k, field, sh); }
		}
	}
	Field_WM1D field;
	field.Build(g, FlatSubstrate(), nd.ExtentScale());
	BoundaryCondition1D bc(field.GetGridPtr()->nx(), e, e);
	switch(t)
	{
	case tParallel:    { WM1D_Traits::ParallelGravity   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return WallLeak<WM1D_Grid>(k, field, sh); }
	//  The 1-D centrifugal term also flings the film out in y (SideOutflow),
	//  by design: not a closed problem
	case tCentrifugal: return 0.0;
	default:           { WM1D_Traits::SolutoCapillarity k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); return WallLeak<WM1D_Grid>(k, field, sh); }
	}
}


//------------------------------------------------------------------
//  The face fluxes of a term (PhysicsTermBase::FaceFlux, SideRate)
//  against its A h: the well-mixed concentration equation carries the
//  solute with them, and the film-height step is A h, so they must be the
//  same operator. Returns, over the interior rows,
//    agree: max |A h - (sum_d +-FaceFlux_d + SideRate)| / max_i sum_j |A_ij h_j|
//    share: max |FaceFlux(s, high face) - FaceFlux(neighbor, low face)|
//           over the interior faces, relative to the largest face flux:
//           the two cells that share a face see the same flux
//------------------------------------------------------------------
template<class Tr, class FieldT, class TermT>
void FaceFluxError(TermT &term, FieldT &field, const Shapes &sh, double &agree, double &share)
{
	typedef typename Tr::Grid Grid;
	typedef typename Tr::Site Site;
	Grid *G = field.GetGridPtr();
	int n = 0;
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
		n = std::max(n, Tr::Index(G, s) + 1);
	std::vector<char> inside(n, 0);
	Vector<double> h(n);
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
	{
		double x, y;
		Tr::SiteXY(G, s, x, y);
		const int i = Tr::Index(G, s);
		h[i] = field.HH(s) = sh.H(x, y);
		field.SS(s) = sh.S(x, y);
		field.CC(s) = sh.C(x, y);
		inside[i] = (G->CheckGhost(s) == Grid::Inside);
	}
	CRS_Matrix A;
	A.Build(n);
	for(int i = 0; i < n; ++i) { A(i, i) = 0.0; if(!inside[i]) A.SetAsBC(i); }
	term.SetMatrixA(A, field);
	Vector<double> Ah(n);
	A.MatVec(Ah, h);
	double worst = 0, scale = 0, fmax = 0, dmax = 0;
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
	{
		const int i = Tr::Index(G, s);
		if(!inside[i]) continue;
		double r = term.SideRate(s, field);
		std::vector<Site> sv(Tr::nStencil);
		G->SetSite(sv, s);
		for(int d = 0; d < Tr::nDir; ++d)
		{
			const double f = term.FaceFlux(s, d, field, sv);
			r += (d % 2 == 0) ? -f : f;
			fmax = std::max(fmax, std::fabs(f));
			if(d % 2 == 1)
			{
				Site a = Tr::Adj(G, s, Tr::Dir(d));
				if(inside[Tr::Index(G, a)])
				{
					std::vector<Site> sa(Tr::nStencil);
					G->SetSite(sa, a);
					dmax = std::max(dmax, std::fabs(f - term.FaceFlux(a, d - 1, field, sa)));
				}
			}
		}
		double row = 0;
		for(int k = 0; k < A.NonZero(i); ++k) row += std::fabs(A.Element(i, k) * h[A.Jcol(i, k)]);
		scale = std::max(scale, row);
		worst = std::max(worst, std::fabs(Ah[i] - r));
	}
	agree = scale > 0 ? worst / scale : 0.0;
	share = fmax > 0 ? dmax / fmax : 0.0;
}

void FaceFlux(bool twoD, Term t, IB::t_Scheme scheme, double &agree, double &share)
{
	Input inp = MakeInput(t);
	NonDimensionalNumbers nd;
	nd.Build(inp, inp.Initial.SurfaceShape.TotalHeight);
	Input::I_Grid g = MakeGrid(twoD, twoD ? 18 : 28, 14, nd.ExtentScale());
	Shapes sh = {g.Xmax, twoD ? g.Ymax : 1.0, t == tDisjoining ? 0.12 : 0.3, LimitedScheme(scheme)};
	const BoundaryCondition::tBoundaryCondition e = BoundaryCondition::Symmetry;
	if(twoD)
	{
		Field_WM2D field;
		field.Build(g, FlatSubstrate(), nd.ExtentScale());
		BoundaryCondition2D bc(field.GetGridPtr()->nx(), field.GetGridPtr()->ny(), e, e, e, e);
		switch(t)
		{
		case tLaplace:     { WM2D_Traits::LaplacePressure   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); FaceFluxError<WM2D_Grid>(k, field, sh, agree, share); return; }
		case tVertical:    { WM2D_Traits::VerticalGravity   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); FaceFluxError<WM2D_Grid>(k, field, sh, agree, share); return; }
		case tParallel:    { WM2D_Traits::ParallelGravity   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); FaceFluxError<WM2D_Grid>(k, field, sh, agree, share); return; }
		case tCentrifugal: { WM2D_Traits::Centrifugal       k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); FaceFluxError<WM2D_Grid>(k, field, sh, agree, share); return; }
		case tSoluto:      { WM2D_Traits::SolutoCapillarity k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); FaceFluxError<WM2D_Grid>(k, field, sh, agree, share); return; }
		default:           { WM2D_Traits::DisjoiningPressure k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); FaceFluxError<WM2D_Grid>(k, field, sh, agree, share); return; }
		}
	}
	Field_WM1D field;
	field.Build(g, FlatSubstrate(), nd.ExtentScale());
	BoundaryCondition1D bc(field.GetGridPtr()->nx(), e, e);
	switch(t)
	{
	case tLaplace:     { WM1D_Traits::LaplacePressure   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); FaceFluxError<WM1D_Grid>(k, field, sh, agree, share); return; }
	case tVertical:    { WM1D_Traits::VerticalGravity   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); FaceFluxError<WM1D_Grid>(k, field, sh, agree, share); return; }
	case tParallel:    { WM1D_Traits::ParallelGravity   k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); FaceFluxError<WM1D_Grid>(k, field, sh, agree, share); return; }
	case tCentrifugal: { WM1D_Traits::Centrifugal       k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); FaceFluxError<WM1D_Grid>(k, field, sh, agree, share); return; }
	case tSoluto:      { WM1D_Traits::SolutoCapillarity k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); k.SetConvectionScheme(scheme); FaceFluxError<WM1D_Grid>(k, field, sh, agree, share); return; }
	default:           { WM1D_Traits::DisjoiningPressure k; k.Build(inp.Physics, field.GetGridPtr(), &nd, &bc); FaceFluxError<WM1D_Grid>(k, field, sh, agree, share); return; }
	}
}

}	// namespace

//  The central difference is exact to O(eps^2) ~ 1e-12 here, and rounding in
//  R(h +- eps v) / 2 eps costs about 1e-10. A term whose Jacobian leaves
//  anything out is off by the size of what it leaves out: the three faults
//  this test found were 2e-2, 3e-2 and 0.5.
const double kTol = 1e-7;

TEST_CASE("Jacobian: every term, 1-D, against the finite difference")
{
	const Term all[] = {tLaplace, tVertical, tParallel, tCentrifugal, tSoluto, tDisjoining};
	for(Term t : all)
	{
		CAPTURE(Name(t));
		CHECK(Error1D(t) < kTol);
	}
}

TEST_CASE("Jacobian: every term, 2-D, against the finite difference")
{
	const Term all[] = {tLaplace, tVertical, tParallel, tCentrifugal, tSoluto, tDisjoining};
	for(Term t : all)
	{
		CAPTURE(Name(t));
		CHECK(Error2D(t) < kTol);
	}
}

TEST_CASE("Jacobian: the capillary term next to a pinned edge")
{
	//  The faces on a pinned edge carry no flux, in A and in J alike
	CHECK(Error1D(tLaplace, true) < kTol);
	CHECK(Error2D(tLaplace, true) < kTol);
}

TEST_CASE("Jacobian: the upwinded terms with every convection scheme")
{
	//  Numerics.Convection.Height. The limited schemes carry the exact
	//  derivative of the limiter in the Jacobian. They are tested on the
	//  film with fronts (Shapes::front), where no face sits within the
	//  finite difference's step of one of the limiter's kinks; the step is
	//  1e-7, since van Leer's curvature, 1/(a + b) at a face with small
	//  differences a and b, makes the central difference itself 4e-7 off
	//  at 1e-6 in 2-D (its error falls as the step squared: 4e-5, 4e-7,
	//  4e-9 at 1e-5, 1e-6, 1e-7).
	const Term up[] = {tParallel, tCentrifugal, tSoluto};
	const IB::t_Scheme all[] = {IB::KK, IB::Upwind_3rd, IB::Central_4th, IB::QUICK,
	                            IB::Upwind_1st, IB::Koren, IB::VanLeer};
	for(IB::t_Scheme sc : all)
		for(Term t : up)
		{
			CAPTURE(IB::Key(sc));
			CAPTURE(Name(t));
			const double eps = LimitedScheme(sc) ? 1e-7 : 1e-6;
			CHECK(Error1D(t, false, eps, sc) < kTol);
			CHECK(Error2D(t, false, eps, sc) < kTol);
			//  and with the wall faces, which carry no flux, on pinned edges
			CHECK(Error1D(t, true, eps, sc) < kTol);
			CHECK(Error2D(t, true, eps, sc) < kTol);
		}
}

TEST_CASE("Upwinded terms: no flux through a Symmetry or a Pinning wall")
{
	//  Each upwinded term gives a wall face zero mobility. With the face
	//  value of the film taken from the ghosts at a wall -- a Symmetry ghost
	//  mirrors the film, so the face value is the film at the wall -- the
	//  in-plane body forces would put liquid through it (VolErr -5.9e-5 over
	//  0.2 s, tools/proptest.py). On a grid open at an edge
	//  (Periodic stands in for it here: the shapes are not periodic) the
	//  sum is of the order of the terms themselves, which shows the test
	//  can see a leak.
	const Term up[] = {tParallel, tCentrifugal, tSoluto};
	const IB::t_Scheme all[] = {IB::KK, IB::Upwind_3rd, IB::Central_4th, IB::QUICK, IB::Upwind_2nd,
	                            IB::Central_2nd, IB::Upwind_1st, IB::Koren, IB::VanLeer};
	const BoundaryCondition::tBoundaryCondition walls[] = {BoundaryCondition::Symmetry, BoundaryCondition::Pinning};
	for(IB::t_Scheme sc : all)
		for(Term t : up)
			for(BoundaryCondition::tBoundaryCondition e : walls)
			{
				CAPTURE(IB::Key(sc));
				CAPTURE(Name(t));
				CAPTURE((int)e);
				CHECK(Leak(false, t, e, sc) < 1e-13);
				CHECK(Leak(true,  t, e, sc) < 1e-13);
			}
	CHECK(Leak(false, tParallel, BoundaryCondition::Periodic, IB::KK) > 1e-3);
	CHECK(Leak(true,  tParallel, BoundaryCondition::Periodic, IB::KK) > 1e-3);
}

TEST_CASE("Capillary and vertical gravity: A h face by face is A h, and closed boxes keep their volume")
{
	//  The Newton iteration takes A h of these two terms from
	//  AddConservative (HeightEquationBase::ApplyA), the Jacobian from the
	//  entries of A. The two must be the same operator: agreement to
	//  rounding of the entries times h, on a film and a substrate that are
	//  not flat with a concentration that varies, in 1-D and 2-D (with the
	//  cross term of the capillary flux), next to Symmetry, Periodic and
	//  Pinning edges.
	const BoundaryCondition::tBoundaryCondition edges[] = {
		BoundaryCondition::Symmetry, BoundaryCondition::Periodic, BoundaryCondition::Pinning};
	for(int twoD = 0; twoD < 2; ++twoD)
	for(int vert = 0; vert < 2; ++vert)
	for(BoundaryCondition::tBoundaryCondition e : edges)
	{
		CAPTURE(twoD); CAPTURE(vert); CAPTURE((int)e);
		double agree = 1, sum = 1;
		Conservative(twoD == 1, vert == 1, e, false, agree, sum);
		CHECK(agree < 1e-13);
		//  With mirror ghosts (what Symmetry gives) the box is closed: the
		//  sum over the cells is the rounding of the face fluxes. A pinned
		//  edge closes the capillary flux on its own.
		Conservative(twoD == 1, vert == 1, e, true, agree, sum);
		CHECK(agree < 1e-13);
		CHECK(sum < 1e-14);
	}
	//  The capillary term with Pinning edges and ghosts that are not a
	//  mirror: the pinned faces carry nothing
	for(int twoD = 0; twoD < 2; ++twoD)
	{
		double agree = 1, sum = 1;
		Conservative(twoD == 1, false, BoundaryCondition::Pinning, false, agree, sum);
		CHECK(sum < 1e-14);
	}
	//  ... while an open edge, with the same ghosts, shows a flux
	{
		double agree = 1, sum = 0;
		Conservative(true, false, BoundaryCondition::Symmetry, false, agree, sum);
		CHECK(sum > 1e-6);
	}
}

TEST_CASE("Face fluxes: every term's FaceFlux is its A h, one flux per face")
{
	//  The well-mixed concentration equation carries the solute with the
	//  volume flux of each face (HeightEquationBase::FaceVolumeFluxes); a
	//  uniform concentration stays uniform only if those are the fluxes of
	//  the film-height equation. Every term, every face interpolation of the
	//  upwinded ones, 1-D (with the sideways outflow of the centrifugal
	//  term) and 2-D, on a film, a substrate and a concentration that vary.
	const Term all[] = {tLaplace, tVertical, tParallel, tCentrifugal, tSoluto, tDisjoining};
	const IB::t_Scheme schemes[] = {IB::KK, IB::Upwind_3rd, IB::Central_4th, IB::QUICK,
	                                IB::Upwind_1st, IB::Koren, IB::VanLeer};
	for(Term t : all)
		for(IB::t_Scheme sc : schemes)
		{
			if(!(t == tParallel || t == tCentrifugal || t == tSoluto) && sc != IB::KK) continue;
			for(int twoD = 0; twoD < 2; ++twoD)
			{
				CAPTURE(Name(t)); CAPTURE(IB::Key(sc)); CAPTURE(twoD);
				double agree = 1, share = 1;
				FaceFlux(twoD == 1, t, sc, agree, share);
				CHECK(agree < 1e-13);
				CHECK(share < 1e-13);
			}
		}
}


//--------------------------------------------------------------------------
//  The same with Navier slip (Mobility.h): the mobility h^3 + 3 Lambda h^2,
//  the shear mobility h^2 + 2 Lambda h, the face mobility of the capillary
//  term from g(h) = int dh / D(h). Two slip lengths against the film scale
//  h0 = 10 um: 3 um (Lambda = 0.3, a / h about 1 to 10: g as written) and
//  10 nm (Lambda = 1e-3, a / h below 0.1: g from its series). The series
//  of D_f for two nearly equal heights, which these films do not reach, is
//  tested on its own in test_positivity_mobility.cpp.
//--------------------------------------------------------------------------
namespace {
const double kSlips[] = {3e-6, 1e-8};
}

TEST_CASE("Slip: the Jacobian of every term, 1-D and 2-D, and next to a pinned edge")
{
	const Term all[] = {tLaplace, tVertical, tParallel, tCentrifugal, tSoluto, tDisjoining};
	for(double len : kSlips)
	{
		SlipScope slip(len);
		for(Term t : all)
		{
			CAPTURE(len); CAPTURE(Name(t));
			CHECK(Error1D(t) < kTol);
			CHECK(Error2D(t) < kTol);
		}
		CHECK(Error1D(tLaplace, true) < kTol);
		CHECK(Error2D(tLaplace, true) < kTol);
	}
}

TEST_CASE("Slip: the upwinded terms with every convection scheme")
{
	const Term up[] = {tParallel, tCentrifugal, tSoluto};
	const IB::t_Scheme all[] = {IB::KK, IB::Upwind_3rd, IB::Central_4th, IB::QUICK,
	                            IB::Upwind_1st, IB::Koren, IB::VanLeer};
	SlipScope slip(kSlips[0]);
	for(IB::t_Scheme sc : all)
		for(Term t : up)
		{
			CAPTURE(IB::Key(sc)); CAPTURE(Name(t));
			const double eps = LimitedScheme(sc) ? 1e-7 : 1e-6;
			CHECK(Error1D(t, false, eps, sc) < kTol);
			CHECK(Error2D(t, false, eps, sc) < kTol);
			CHECK(Error1D(t, true, eps, sc) < kTol);
			CHECK(Error2D(t, true, eps, sc) < kTol);
		}
}

TEST_CASE("Slip: no flux through a wall, A h face by face, one flux per face")
{
	const Term up[] = {tParallel, tCentrifugal, tSoluto};
	const IB::t_Scheme upw[] = {IB::KK, IB::Koren, IB::VanLeer, IB::Upwind_1st};
	const BoundaryCondition::tBoundaryCondition walls[] = {BoundaryCondition::Symmetry, BoundaryCondition::Pinning};
	const Term all[] = {tLaplace, tVertical, tParallel, tCentrifugal, tSoluto, tDisjoining};
	for(double len : kSlips)
	{
		SlipScope slip(len);
		CAPTURE(len);
		for(IB::t_Scheme sc : upw)
			for(Term t : up)
				for(BoundaryCondition::tBoundaryCondition e : walls)
				{
					CAPTURE(IB::Key(sc)); CAPTURE(Name(t)); CAPTURE((int)e);
					CHECK(Leak(false, t, e, sc) < 1e-13);
					CHECK(Leak(true,  t, e, sc) < 1e-13);
				}
		const BoundaryCondition::tBoundaryCondition edges[] = {
			BoundaryCondition::Symmetry, BoundaryCondition::Periodic, BoundaryCondition::Pinning};
		for(int twoD = 0; twoD < 2; ++twoD)
		for(int vert = 0; vert < 2; ++vert)
		for(BoundaryCondition::tBoundaryCondition e : edges)
		{
			CAPTURE(twoD); CAPTURE(vert); CAPTURE((int)e);
			double agree = 1, sum = 1;
			Conservative(twoD == 1, vert == 1, e, true, agree, sum);
			CHECK(agree < 1e-13);
			CHECK(sum < 1e-14);
		}
		for(Term t : all)
			for(IB::t_Scheme sc : upw)
			{
				if(!(t == tParallel || t == tCentrifugal || t == tSoluto) && sc != IB::KK) continue;
				for(int twoD = 0; twoD < 2; ++twoD)
				{
					CAPTURE(Name(t)); CAPTURE(IB::Key(sc)); CAPTURE(twoD);
					double agree = 1, share = 1;
					FaceFlux(twoD == 1, t, sc, agree, share);
					CHECK(agree < 1e-13);
					CHECK(share < 1e-13);
				}
			}
	}
}
