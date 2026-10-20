/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_disjoining.cpp : verification of the disjoining pressure
*                        (wettability / contact line)
*
*  There are two things we want to pin down here.
*
*  (1) The coefficient Dp is consistent with the contact angle that was
*      asked for
*
*      The equilibrium shape is the one with "uniform pressure", that is
*
*          H'' + Dp * Pihat(H) = K
*
*      When the precursor film extends to infinity (K = 0), this can be
*      integrated once after multiplying by H',
*
*          (1/2) H'^2 = -Dp * int_{H*}^{H} Pihat ds
*
*      For H >> H* the right-hand side converges to
*      Dp*H*(n-m)/((n-1)(m-1)), so
*
*          (H')^2 -> As * theta_e^2
*
*      The physical slope is the dimensionless slope H' multiplied by
*      eps = h0/L0, so the slope on the outside of the film is exactly
*      theta_e. In other words, whether "the formula for Dp is consistent
*      with the contact angle that was asked for" can be checked in closed
*      form.
*
*  (2) The discretized term drives no flow on that equilibrium shape
*
*      We build the equilibrium shape (the one satisfying
*      H'' + Dp*Pihat = K) with Runge-Kutta, compute the capillary
*      pressure + disjoining pressure velocity, and confirm that the
*      remaining velocity approaches 0 at second order as the mesh is
*      refined. If the coefficient, the sign or the discretization is
*      wrong in any way, this residual does not vanish.
*
*  If (1) holds but (2) does not, the discretization is wrong; if (2)
*  holds but (1) does not, the problem is "solved cleanly, but at the
*  wrong angle".
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "doctest/doctest.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "NonDimensionalNumbers.h"
#include "LubFields.h"
#include "BoundaryCondition.h"
#include "LaplacePressure_1D.h"
#include "DisjoiningPressureKernel.h"

using namespace yuzen;
using namespace yuzen::grid;

namespace {

const double PI = 3.14159265358979323846;

//------------------------------------------------------------------
//  Assemble the minimal Input that uses the disjoining pressure
//------------------------------------------------------------------
Input MakeInput(double theta_deg, double hstar, double n, double m,
                double h0 = 10e-6)
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
	inp.Fluid.SurfaceTension.Coef      = 0.0;
	inp.Fluid.Viscosity.BaseValue      = 100.0e-3;
	inp.Fluid.Viscosity.BaseConc       = 0.5;
	inp.Fluid.Viscosity.Exponent       = 0.0;	//  keep the viscosity constant
	inp.Fluid.Diffusivity.BaseValue    = 1.0e-9;
	inp.Fluid.Diffusivity.BaseConc     = 0.5;
	inp.Fluid.Diffusivity.Exponent     = 0.0;
	inp.Fluid.MolarMass.Solvent        = 100.0;
	inp.Fluid.MolarMass.Resin          = 2000.0;

	inp.Initial.SurfaceShape.TotalHeight    = h0;
	inp.Initial.SurfaceShape.PrecursorRatio = hstar;

	inp.Physics.ConcField.Model            = "WM";
	inp.Physics.LaplacePressure.Enable     = true;
	inp.Physics.SolutoCapillarity.Enable   = false;
	inp.Physics.VerticalGravity.Enable     = false;
	inp.Physics.VerticalGravity.Acceleration     = 9.80665;
	inp.Physics.HorizontalGravity.Enable   = false;
	inp.Physics.HorizontalGravity.Acceleration.assign(2, 0.0);
	inp.Physics.Centrifugal.Enable         = false;
	inp.Physics.Centrifugal.RotationRate   = 0.0;
	inp.Physics.Evaporation.Enable         = false;
	inp.Physics.Evaporation.Coef           = 0.0;

	inp.Physics.DisjoiningPressure.Enable       = true;
	inp.Physics.DisjoiningPressure.ContactAngle = theta_deg;
	inp.Physics.DisjoiningPressure.ExponentN    = n;
	inp.Physics.DisjoiningPressure.ExponentM    = m;

	return inp;
}

Input::I_Grid MakeGrid1D(int nx, double xmin, double xmax)
{
	Input::I_Grid g;
	g.Dimension = "xh";
	g.Xmin = xmin; g.Xmax = xmax;
	g.Ymin = -1.0; g.Ymax = 1.0;
	g.Nx = nx; g.Ny = 1; g.Nz = 1;
	return g;
}

Input::I_Substrate FlatSubstrate()
{
	Input::I_Substrate s;
	s.TotalHeight = 0.0;
	return s;
}

//------------------------------------------------------------------
//  Pihat and its antiderivative
//------------------------------------------------------------------
struct Disjoin
{
	double hstar, n, m;
	double Pi(double H) const
	{
		const double r = hstar / H;
		return std::pow(r, n) - std::pow(r, m);
	}
	//  int_{H*}^{H} Pihat ds
	double W(double H) const
	{
		const double a = std::pow(hstar, n) / (1.0 - n)
		               * (std::pow(H, 1.0 - n) - std::pow(hstar, 1.0 - n));
		const double b = std::pow(hstar, m) / (1.0 - m)
		               * (std::pow(H, 1.0 - m) - std::pow(hstar, 1.0 - m));
		return a - b;
	}
};

//------------------------------------------------------------------
//  An equilibrium droplet. The surrounding flat film Hf is determined from
//  the apex height Hmax by bisection, then the shape is integrated with
//  Runge-Kutta from the foot up to the apex.
//
//      H'' = K - Dp*Pihat(H),   K = Dp*Pihat(Hf)
//------------------------------------------------------------------
class EquilibriumDroplet
{
public:
	EquilibriumDroplet(double Hmax, double hstar, double n, double m, double Dp)
		: _Dp(Dp)
	{
		_d.hstar = hstar; _d.n = n; _d.m = m;

		//  Hf: G(Hf) = K(Hmax-Hf) - Dp[W(Hmax)-W(Hf)] = 0
		double lo = hstar * (1.0 + 1e-12);
		double hi = hstar * std::pow(n / m, 1.0 / (n - m));
		for(int i = 0; i < 200; i++)
		{
			const double mid = 0.5 * (lo + hi);
			if(G(lo, Hmax) * G(mid, Hmax) <= 0) hi = mid;
			else                                lo = mid;
		}
		_Hf = 0.5 * (lo + hi);
		_K  = _Dp * _d.Pi(_Hf);

		//  Exponential growth rate at the foot
		const double lam = std::sqrt(_Dp * (n - m) / hstar
		                             * std::pow(hstar / _Hf, m) * (1.0));
		//  (Strictly, lam^2 = -Dp*Pihat'(Hf). Since Hf ~ H* the
		//   approximation above is good enough, but we replace it with the
		//   exact expression to be safe.)
		const double r  = hstar / _Hf;
		const double dPi = -(n * std::pow(r, n) - m * std::pow(r, m)) / _Hf;
		_lam = std::sqrt(-_Dp * dPi);
		(void)lam;

		Integrate();
	}

	double Hf() const {return _Hf;}
	double K()  const {return _K;}
	double HalfWidth() const {return _Xapex;}
	double ApexSlope() const {return _apexSlope;}

	//  The symmetric shape, with the apex at X = 0.
	//
	//  Linear interpolation is not good enough here, so this uses cubic
	//  Hermite interpolation. Taking a second difference on this shape
	//  divides by dX^2 (the square of the mesh spacing), so the
	//  interpolation error turns directly into an error in the pressure.
	//  The error of linear interpolation is O(_dX^2) and does not depend
	//  on the mesh spacing, so refining the mesh would actually increase
	//  the "pressure error". With Hermite interpolation it is O(_dX^4),
	//  which does not interfere with the mesh-refinement test.
	double H(double X) const
	{
		const double r = std::fabs(X);
		if(r >= _Xapex) return _Hf;
		//  _prof[i] is the position i*_dX counted from the foot.
		//  X = 0 (the apex) is the last entry.
		const double s = (_Xapex - r) / _dX;
		int i = (int)s;
		if(i < 0) i = 0;
		if(i >= (int)_prof.size() - 1) return _prof.back();
		const double t  = s - i;
		const double t2 = t * t, t3 = t2 * t;
		const double h00 =  2*t3 - 3*t2 + 1.0;
		const double h10 =      t3 - 2*t2 + t;
		const double h01 = -2*t3 + 3*t2;
		const double h11 =      t3 -   t2;
		return h00 * _prof[i]   + h10 * _dX * _dprof[i]
		     + h01 * _prof[i+1] + h11 * _dX * _dprof[i+1];
	}

private:
	Disjoin _d;
	double _Dp, _Hf, _K, _lam, _Xapex, _dX, _apexSlope;
	std::vector<double> _prof, _dprof;

	double G(double Hf, double Hmax) const
	{
		const double K = _Dp * _d.Pi(Hf);
		return K * (Hmax - Hf) - _Dp * (_d.W(Hmax) - _d.W(Hf));
	}

	//  Changing the initial value eta0 at the foot into eta0*c shifts the
	//  shape by -ln(c)/lambda along X. We use this to tune eta0 so that
	//  "the last grid point lands exactly on the apex (H' = 0)".
	//
	//  Without this, the last point is off the apex by up to _dX, and
	//  mirroring the shape left-right creates a kink at the apex. The
	//  resulting height error is constant and does not depend on _dX, so
	//  once it passes through a second difference (dividing by dx^2) the
	//  pressure error grows as the mesh is refined, and the order of
	//  convergence can no longer be measured.
	void Integrate()
	{
		_dX = 0.005 / _lam;
		double eta0 = 1e-12;

		//  First determine the number of steps N (until we pass the apex)
		int N = 0;
		{
			double y0 = _Hf + eta0, y1 = _lam * eta0;
			for(int k = 0; k < 4000000; k++)
			{
				double n0, n1;
				Step(y0, y1, n0, n1);
				if(n1 <= 0) break;
				y0 = n0; y1 = n1; N++;
			}
			//  Shift the shape nearer by the remaining distance from the
			//  last point to the apex
			eta0 *= std::exp(_lam * Gap(y0, y1));
		}

		//  Tune eta0 so that N steps land exactly on the apex (H' = 0)
		for(int pass = 0; pass < 8; pass++)
		{
			double y0 = _Hf + eta0, y1 = _lam * eta0;
			for(int k = 0; k < N; k++)
			{
				double n0, n1;
				Step(y0, y1, n0, n1);
				y0 = n0; y1 = n1;
			}
			const double g = Gap(y0, y1);
			eta0 *= std::exp(_lam * g);
			if(std::fabs(g) < 1e-14 * _dX) break;
		}

		//  Store the shape using the eta0 we settled on
		{
			double y0 = _Hf + eta0, y1 = _lam * eta0;
			_prof.clear();  _dprof.clear();
			_prof.push_back(y0);  _dprof.push_back(y1);
			for(int k = 0; k < N; k++)
			{
				double n0, n1;
				Step(y0, y1, n0, n1);
				y0 = n0; y1 = n1;
				_prof.push_back(y0);  _dprof.push_back(y1);
			}
			_apexSlope = y1;
		}
		_Xapex = _dX * N;
	}

	//  One RK4 step
	void Step(double p0, double p1, double &n0, double &n1) const
	{
		double k0[4], k1[4];
		k0[0] = p1;                       k1[0] = _K - _Dp * _d.Pi(p0);
		k0[1] = p1 + 0.5 * _dX * k1[0];   k1[1] = _K - _Dp * _d.Pi(p0 + 0.5 * _dX * k0[0]);
		k0[2] = p1 + 0.5 * _dX * k1[1];   k1[2] = _K - _Dp * _d.Pi(p0 + 0.5 * _dX * k0[1]);
		k0[3] = p1 + _dX * k1[2];         k1[3] = _K - _Dp * _d.Pi(p0 + _dX * k0[2]);
		n0 = p0 + _dX / 6.0 * (k0[0] + 2 * k0[1] + 2 * k0[2] + k0[3]);
		n1 = p1 + _dX / 6.0 * (k1[0] + 2 * k1[1] + 2 * k1[2] + k1[3]);
	}

	//  Distance from the current point to the apex (H' = 0), one Newton step
	double Gap(double y0, double y1) const
	{
		const double acc = _K - _Dp * _d.Pi(y0);	// = H''
		return (acc < 0) ? (y1 / (-acc)) : 0.0;
	}
};


//------------------------------------------------------------------
//  Put the equilibrium shape on the grid and return the maximum of the
//  capillary pressure + disjoining pressure velocity
//------------------------------------------------------------------
double ResidualVelocity(int nx, double theta_deg, double hstar,
                        double n, double m, double Hmax,
                        double xmin, double xmax,
                        double *halfWidth = 0)
{
	Input inp = MakeInput(theta_deg, hstar, n, m);

	NonDimensionalNumbers nd;
	nd.Build(inp, inp.Initial.SurfaceShape.TotalHeight);

	Input::I_Grid g = MakeGrid1D(nx, xmin, xmax);
	g.Xmin /= nd.ExtentScale();
	g.Xmax /= nd.ExtentScale();

	Field_WM1D field;
	field.Build(g, FlatSubstrate(), nd.ExtentScale());
	C_Grid1D *G = field.GetGridPtr();

	EquilibriumDroplet drop(Hmax, hstar, n, m, nd.Dp());
	if(halfWidth) *halfWidth = drop.HalfWidth();

	for(Xedge xe = G->First<Xedge>(); xe != G->Last<Xedge>(); ++xe)
	{
		field.HH(xe) = drop.H(G->X(xe));
		field.SS(xe) = 0.0;
		field.CC(xe) = inp.Fluid.Viscosity.BaseConc;
		field.UU().U_LaplacePressure().x(G->HalfXm(xe))    = 0.0;
		field.UU().U_DisjoiningPressure().x(G->HalfXm(xe)) = 0.0;
	}

	BoundaryCondition1D bc;

	LaplacePressure_WM1D    lap;
	DisjoiningPressureKernel<WM1D_Grid> dis;
	lap.Build(inp.Physics, G, &nd, &bc);
	dis.Build(inp.Physics, G, &nd, &bc);

	lap.CalVelocity(field);
	dis.CalVelocity(field);

	//  The capillary and disjoining pressures should cancel each other.
	//  What remains is the truncation error.
	double umax = 0; int imax = 0;
	for(int ix = 3; ix <= G->nx() - 4; ix++)
	{
		Xedge xe = G->Of<Xedge>(ix);
		Node nd_ = G->HalfXm(xe);
		const double u = field.UU().U_LaplacePressure().x(nd_)
		               + field.UU().U_DisjoiningPressure().x(nd_);
		if(std::fabs(u) > umax) {umax = std::fabs(u); imax = ix;}
	}
	if(std::getenv("YUZEN_DBG_EQ"))
	{
		Xedge xe = G->Of<Xedge>(imax);
		std::fprintf(stderr, "[eq] nx=%d umax=%.4e at ix=%d X=%.6e H=%.6e  (lap=%.4e dis=%.4e)\n",
		  nx, umax, imax, G->X(xe), field.HH(xe),
		  field.UU().U_LaplacePressure().x(G->HalfXm(xe)),
		  field.UU().U_DisjoiningPressure().x(G->HalfXm(xe)));
	}
	return umax;
}

} // namespace


//==================================================================
//  (1) The calibration of Dp = the contact angle
//==================================================================
TEST_CASE("Disjoining pressure: Dp agrees with the equilibrium contact angle asked for")
{
	const double thetas[] = {2.0, 5.0, 10.0, 20.0};
	const double hstars[] = {0.01, 0.05, 0.2};
	struct { double n, m; } exps[] = {{3.0, 2.0}, {4.0, 3.0}, {9.0, 3.0}};

	for(int it = 0; it < 4; it++)
	for(int ih = 0; ih < 3; ih++)
	for(int ie = 0; ie < 3; ie++)
	{
		const double th = thetas[it], hs = hstars[ih];
		const double n = exps[ie].n, m = exps[ie].m;

		Input inp = MakeInput(th, hs, n, m);
		NonDimensionalNumbers nd;
		nd.Build(inp, inp.Initial.SurfaceShape.TotalHeight);

		const double As  = nd.As();
		const double eps = 1.0 / std::sqrt(As);
		const double rad = th * PI / 180.0;

		//  The formula exactly as implemented
		const double Dp_expect = As * (n - 1.0) * (m - 1.0)
		                       / (2.0 * hs * (n - m)) * rad * rad;
		REQUIRE(nd.Dp() == doctest::Approx(Dp_expect).epsilon(1e-14));

		//  The integral of Pihat from the precursor film to infinity
		//  (analytic)
		//    int_{H*}^{inf} Pihat = H*(m-n) / ((n-1)(m-1))   (< 0)
		const double integ = hs * (m - n) / ((n - 1.0) * (m - 1.0));

		//  From the first integral of the equilibrium,
		//  (H')^2 = -2 Dp * integ
		const double slope2 = -2.0 * nd.Dp() * integ;
		REQUIRE(slope2 > 0);

		//  The physical slope eps*H' must agree with theta_e
		const double slope_phys = eps * std::sqrt(slope2);
		REQUIRE(slope_phys == doctest::Approx(rad).epsilon(1e-12));
	}
}


TEST_CASE("Disjoining pressure: Dp is exactly zero when it is disabled")
{
	Input inp = MakeInput(10.0, 0.2, 3.0, 2.0);
	inp.Physics.DisjoiningPressure.Enable = false;

	NonDimensionalNumbers nd;
	nd.Build(inp, inp.Initial.SurfaceShape.TotalHeight);
	REQUIRE(nd.Dp() == 0.0);
}


//==================================================================
//  (2) No flow arises on the equilibrium shape (second-order convergence)
//==================================================================
TEST_CASE("Disjoining pressure: the residual velocity on the equilibrium shape falls at second order")
{
	const double theta = 10.0;
	const double hstar = 0.02;
	const double n = 3.0, m = 2.0;
	const double Hmax = 1.0;
	const double xmin = -0.3e-3, xmax = 0.3e-3;

	double hw = 0;
	const double u1 = ResidualVelocity(1000, theta, hstar, n, m, Hmax, xmin, xmax, &hw);
	const double u2 = ResidualVelocity(2000, theta, hstar, n, m, Hmax, xmin, xmax);
	const double u4 = ResidualVelocity(4000, theta, hstar, n, m, Hmax, xmin, xmax);

	//  The droplet must fit inside the domain (otherwise the checks below
	//  are meaningless)
	REQUIRE(hw > 0);
	REQUIRE(hw < 0.3e-3 / 1.749042e-03);

	INFO("u(1000) = " << u1 << ", u(2000) = " << u2 << ", u(4000) = " << u4);

	//  Refining the mesh must always make it smaller
	REQUIRE(u2 < u1);
	REQUIRE(u4 < u2);

	//  The ratio must be close to 4 (second order). The interpolation error
	//  from building the shape rides on top of it, so we allow a band of
	//  3.0 - 5.0. If a sign or a coefficient is wrong, the ratio comes out
	//  near 1.
	REQUIRE(u1 / u2 > 3.0);
	REQUIRE(u1 / u2 < 5.0);
	REQUIRE(u2 / u4 > 3.0);
	REQUIRE(u2 / u4 < 5.0);
}


//==================================================================
//  A flat precursor film is itself an equilibrium (Pihat(H*) = 0)
//==================================================================
TEST_CASE("Disjoining pressure: a uniform precursor film produces no velocity")
{
	const double hstar = 0.02;
	Input inp = MakeInput(10.0, hstar, 3.0, 2.0);

	NonDimensionalNumbers nd;
	nd.Build(inp, inp.Initial.SurfaceShape.TotalHeight);

	Input::I_Grid g = MakeGrid1D(64, -0.3e-3, 0.3e-3);
	g.Xmin /= nd.ExtentScale();
	g.Xmax /= nd.ExtentScale();

	Field_WM1D field;
	field.Build(g, FlatSubstrate(), nd.ExtentScale());
	C_Grid1D *G = field.GetGridPtr();

	for(Xedge xe = G->First<Xedge>(); xe != G->Last<Xedge>(); ++xe)
	{
		field.HH(xe) = hstar;
		field.SS(xe) = 0.0;
		field.CC(xe) = 0.5;
		field.UU().U_DisjoiningPressure().x(G->HalfXm(xe)) = 1.0;	// dirty it first
	}

	BoundaryCondition1D bc;
	DisjoiningPressureKernel<WM1D_Grid> dis;
	dis.Build(inp.Physics, G, &nd, &bc);
	dis.CalVelocity(field);

	for(int ix = 3; ix <= G->nx() - 4; ix++)
	{
		Xedge xe = G->Of<Xedge>(ix);
		REQUIRE(std::fabs(field.UU().U_DisjoiningPressure().x(G->HalfXm(xe))) < 1e-12);
	}
}
