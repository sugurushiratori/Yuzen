/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_finite_differences.cpp : regression tests for numerics/FiniteDifferences.h
*
*  These check the consistency of the interpolation coefficients and that
*  the gradient agrees with the analytic solution. The solver does not call
*  GradX/GradY/GradZ, so their sign (left and right stencils the right way
*  round in every overload) is pinned down here.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "doctest/doctest.h"
#include "grid/YuzenGrid.h"
using namespace yuzen::grid;
#include "FiniteDifferences.h"

#include <cmath>

using namespace yuzen;

namespace {

typedef C_Interpolation_base4 IB;

const IB::t_Scheme kSchemes[9] = {
	IB::Upwind_1st, IB::Upwind_2nd, IB::Central_2nd,
	IB::Upwind_3rd, IB::QUICK, IB::KK, IB::Central_4th,
	IB::Koren, IB::VanLeer };

} // namespace

TEST_CASE("C_Interpolation_base4: the interpolation coefficients sum to one (a constant is preserved)")
{
	for(int k = 0; k < 9; ++k)
	{
		C_FiniteDifferences fd;
		fd.Build(kSchemes[k]);
		CAPTURE(fd.SchemeName());

		double s1 = 0.0;
		for(int i = 0; i < fd.NumStencils(); ++i) s1 += fd.c1(i);
		CHECK(s1 == doctest::Approx(1.0));		// the sign-independent part is 1

		double s2 = 0.0;
		for(int i = 0; i < fd.NumStencils(); ++i) s2 += fd.c2(i);
		CHECK(s2 == doctest::Approx(0.0));		// the sign-dependent part cancels out
	}
}

TEST_CASE("C_FiniteDifferences: the gradient of a linear function is reproduced exactly")
{
	// f(x) = a*x + b  ->  df/dx = a  (exact for every scheme)
	const double a = 3.0, b = 1.0;
	RectGrid<1> G;
	G.Build(0.0, 1.0, 100);
	ScalarField1D<> s(G);
	for(int i = G.imin(); i <= G.imax(); ++i)
	{
		Xedge xe = G.Of<Xedge>(i);
		s(xe) = a * G.X(xe) + b;
	}

	for(int k = 0; k < 7; ++k)
	{
		C_FiniteDifferences fd;
		fd.Build(kSchemes[k]);
		CAPTURE(fd.SchemeName());

		// The 3 points at each end are affected by the periodic wrap, so
		// we look at the interior only
		for(int i = G.imin() + 3; i <= G.imax() - 3; ++i)
		{
			Xedge xe = G.Of<Xedge>(i);
			CAPTURE(i);
			CHECK(fd.GradX(s, xe) == doctest::Approx(a).epsilon(1e-9));
		}
	}
}

TEST_CASE("C_FiniteDifferences: Central_2nd gives the gradient of a quadratic exactly")
{
	// f(x) = x^2  ->  df/dx = 2x
	RectGrid<1> G;
	G.Build(0.0, 1.0, 400);
	ScalarField1D<> s(G);
	for(int i = G.imin(); i <= G.imax(); ++i)
	{
		Xedge xe = G.Of<Xedge>(i);
		const double x = G.X(xe);
		s(xe) = x * x;
	}

	C_FiniteDifferences fd;
	fd.Build(IB::Central_2nd);

	double err = 0.0;
	for(int i = G.imin() + 3; i <= G.imax() - 3; ++i)
	{
		Xedge xe = G.Of<Xedge>(i);
		err = std::max(err, std::fabs(fd.GradX(s, xe) - 2.0 * G.X(xe)));
	}
	CAPTURE(err);
	CHECK(err < 1e-10);		// the 2nd-order central difference is exact for a quadratic
}

TEST_CASE("C_FiniteDifferences: the Node variant reproduces a linear gradient exactly")
{
	const double a = -2.5;
	RectGrid<1> G;
	G.Build(0.0, 1.0, 100);
	ScalarField1D<> s(G);
	for(int i = G.imin(); i <= G.imax(); ++i)
	{
		Xedge xe = G.Of<Xedge>(i);
		s(xe) = a * G.X(xe);
	}

	C_FiniteDifferences fd;
	fd.Build(IB::Central_2nd);

	for(int i = G.imin() + 4; i <= G.imax() - 4; ++i)
	{
		Node nd = G.Of<Node>(i);
		CAPTURE(i);
		CHECK(fd.GradX(s, nd) == doctest::Approx(a).epsilon(1e-9));
	}
}


//==================================================================
//  The face values of the convection schemes (Numerics.Convection)
//==================================================================
namespace {

//  The advection derivative a row forms from its two faces, as the
//  concentration equations do: (f_face(i+1/2) - f_face(i-1/2)) / dx, on
//  the points x_i + k dx, k = -2..2, of f.
template<class F>
double RowDerivative(const IB &sc, F f, double x0, double dx, double sign)
{
	double fl[4], fr[4];
	for(int j = 0; j < 4; ++j)
	{
		fl[j] = f(x0 + (j - 2) * dx);
		fr[j] = f(x0 + (j - 1) * dx);
	}
	return (sc.FaceValue(fr, sign) - sc.FaceValue(fl, sign)) / dx;
}

//  Polynomials the derivative must be exact for: degree 3 for the
//  third-order schemes (KK, UTOPIA, Koren where it is not limited), 2 for
//  the second-order ones, 1 for van Leer (nonlinear even where smooth), 0
//  for first-order upwind.
int ExactDegree(IB::t_Scheme t)
{
	switch(t)
	{
	case IB::Upwind_1st: return 0;
	case IB::VanLeer:    return 1;
	case IB::KK: case IB::Upwind_3rd: case IB::Koren: return 3;
	default:             return 2;
	}
}

//  Pseudo-random numbers, reproducible
double Rand(unsigned &st) { st = st * 1103515245u + 12345u; return ((st >> 8) & 0xffff) / 65535.0; }

}	// namespace

TEST_CASE("Convection schemes: the input keys map one to one onto the schemes")
{
	for(int k = 0; k < 9; ++k)
	{
		IB::t_Scheme t;
		const std::string key = IB::Key(kSchemes[k]);
		CAPTURE(key);
		REQUIRE(!key.empty());
		REQUIRE(IB::FromKey(key, t));
		CHECK(t == kSchemes[k]);
		//  and each has a name for the log
		IB sc; sc.Build(kSchemes[k]);
		CHECK(!sc.SchemeName().empty());
	}
	IB::t_Scheme t;
	CHECK(!IB::FromKey("kk", t));
	CHECK(!IB::FromKey("WENO5", t));
}

TEST_CASE("Convection schemes: the advection derivative is exact on polynomials of the scheme's degree")
{
	//  Monotone data whose ratio of successive differences stays inside
	//  0.4 < r < 4, where Koren is not limited
	const double c[4] = {0.3, 2.0, 0.4, 0.15};
	for(int k = 0; k < 9; ++k)
	{
		IB sc; sc.Build(kSchemes[k]);
		const int deg = ExactDegree(kSchemes[k]);
		CAPTURE(sc.SchemeName());
		for(int p = 0; p <= deg; ++p)
		{
			auto f = [&](double x) { double v = 0, xp = 1; for(int q = 0; q <= p; ++q) { v += c[q] * xp; xp *= x; } return v; };
			auto df = [&](double x) { double v = 0, xp = 1; for(int q = 1; q <= p; ++q) { v += q * c[q] * xp; xp *= x; } return v; };
			for(double sg : {+1.0, -1.0})
			{
				CAPTURE(p); CAPTURE(sg);
				const double x0 = 0.37, dx = 0.05;
				CHECK(RowDerivative(sc, f, x0, dx, sg) == doctest::Approx(df(x0)).epsilon(1e-10));
			}
		}
	}
}

TEST_CASE("Convection schemes: order of the advection derivative on a smooth monotone profile")
{
	//  f = tanh(x): the error falls as dx^p, p the degree the derivative is
	//  exact for (1 for first-order upwind), and 2 for van Leer, whose
	//  limiter is nonlinear even where the profile is smooth. (Central4's
	//  face value is a fourth-order interpolant, but the difference of two
	//  such face values is a second-order derivative.)
	auto f  = [](double x) { return std::tanh(x); };
	auto df = [](double x) { const double t = std::tanh(x); return 1.0 - t * t; };
	for(int k = 0; k < 9; ++k)
	{
		IB sc; sc.Build(kSchemes[k]);
		CAPTURE(sc.SchemeName());
		const int deg = ExactDegree(kSchemes[k]);
		const int expect = (kSchemes[k] == IB::VanLeer) ? 2 : std::max(deg, 1);
		double e[2];
		for(int m = 0; m < 2; ++m)
		{
			const double dx = 0.02 / (1 << m);
			e[m] = 0;
			for(double x0 = -0.8; x0 <= 0.8; x0 += 0.1)
				e[m] = std::max(e[m], std::fabs(RowDerivative(sc, f, x0, dx, 1.0) - df(x0)));
		}
		const double order = std::log(e[0] / e[1]) / std::log(2.0);
		CAPTURE(order);
		CHECK(order > expect - 0.15);
		CHECK(order < expect + 0.3);
	}
}

TEST_CASE("Convection schemes: the limited face value is bounded by its two neighbors (TVD)")
{
	unsigned st = 7;
	for(IB::t_Scheme t : {IB::Koren, IB::VanLeer})
	{
		IB sc; sc.Build(t);
		CAPTURE(sc.SchemeName());
		for(int n = 0; n < 2000; ++n)
		{
			double f[4];
			for(int j = 0; j < 4; ++j) f[j] = Rand(st) < 0.3 ? std::floor(Rand(st) * 2.0) : Rand(st);
			for(double sg : {+1.0, -1.0})
			{
				const int u = sg > 0 ? IB::L : IB::R, d = sg > 0 ? IB::R : IB::L, uu = sg > 0 ? IB::LL : IB::RR;
				const double v = sc.FaceValue(f, sg);
				//  between f_U and f_D, and within 2 (f_U - f_UU) of f_U
				CHECK(v >= std::min(f[u], f[d]) - 1e-15);
				CHECK(v <= std::max(f[u], f[d]) + 1e-15);
				const double a = f[u] - f[uu];
				CHECK(std::fabs(v - f[u]) <= std::fabs(a) + 1e-15);
				//  first order at an extremum
				if((f[u] - f[uu]) * (f[d] - f[u]) <= 0) CHECK(v == f[u]);
			}
		}
		//  A step: the face takes the upwind value, whichever side the step is
		const double step[3][4] = {{0, 0, 1, 1}, {0, 1, 1, 1}, {0, 0, 0, 1}};
		for(int i = 0; i < 3; ++i)
		{
			CHECK(sc.FaceValue(step[i], +1.0) >= 0.0);
			CHECK(sc.FaceValue(step[i], +1.0) <= 1.0);
			CHECK(sc.FaceValue(step[i], -1.0) >= 0.0);
			CHECK(sc.FaceValue(step[i], -1.0) <= 1.0);
		}
		CHECK(sc.FaceValue(step[0], +1.0) == 0.0);
		CHECK(sc.FaceValue(step[0], -1.0) == 1.0);
	}
	//  ... which the linear schemes are not: the face of the step overshoots
	for(IB::t_Scheme t : {IB::KK, IB::Upwind_3rd, IB::Central_4th, IB::QUICK})
	{
		IB sc; sc.Build(t);
		CAPTURE(sc.SchemeName());
		const double s1[4] = {0, 0, 0, 1}, s2[4] = {0, 1, 1, 1};
		CHECK((sc.FaceValue(s1, +1.0) < 0.0 || sc.FaceValue(s2, +1.0) > 1.0 || sc.FaceValue(s1, -1.0) < 0.0));
	}
}

TEST_CASE("Convection schemes: mirror and sign symmetry of the face value")
{
	unsigned st = 11;
	for(int k = 0; k < 9; ++k)
	{
		IB sc; sc.Build(kSchemes[k]);
		CAPTURE(sc.SchemeName());
		for(int n = 0; n < 500; ++n)
		{
			double f[4], fm[4], fn[4];
			for(int j = 0; j < 4; ++j) f[j] = Rand(st) - 0.3;
			for(int j = 0; j < 4; ++j) { fm[j] = f[3 - j]; fn[j] = -f[j]; }
			for(double sg : {+1.0, -1.0})
			{
				if(sc.IsLimited())
				{	//  to the last bit
					CHECK(sc.FaceValue(fm, -sg) == sc.FaceValue(f, sg));
					CHECK(sc.FaceValue(fn, sg) == -sc.FaceValue(f, sg));
				}
				else
				{
					CHECK(sc.FaceValue(fm, -sg) == doctest::Approx(sc.FaceValue(f, sg)).epsilon(1e-14));
					CHECK(sc.FaceValue(fn, sg) == doctest::Approx(-sc.FaceValue(f, sg)).epsilon(1e-14));
				}
			}
		}
	}
}

TEST_CASE("Convection schemes: frozen weights reproduce the face value; the derivative matches a finite difference")
{
	unsigned st = 5;
	for(int k = 0; k < 9; ++k)
	{
		IB sc; sc.Build(kSchemes[k]);
		CAPTURE(sc.SchemeName());
		for(int n = 0; n < 500; ++n)
		{
			double f[4];
			for(int j = 0; j < 4; ++j) f[j] = Rand(st);
			for(double sg : {+1.0, -1.0})
			{
				const double v = sc.FaceValue(f, sg);
				for(int rowPos : {-1, (int)IB::L, (int)IB::R})
				{
					double w[4], sum = 0, wsum = 0;
					sc.FaceWeights(f, sg, rowPos, w);
					for(int j = 0; j < 4; ++j) { sum += w[j] * f[j]; wsum += w[j]; }
					CHECK(sum == doctest::Approx(v).epsilon(1e-13));
					CHECK(wsum == doctest::Approx(1.0).epsilon(1e-13));
					if(sc.IsLimited())
					{	//  the signs that make the frozen advection an M-matrix
						const int u = sg > 0 ? IB::L : IB::R, d = sg > 0 ? IB::R : IB::L, uu = sg > 0 ? IB::LL : IB::RR;
						if(rowPos == u) { CHECK(w[uu] <= 0.0); CHECK(w[u] >= 1.0); CHECK(w[u] <= 2.0); CHECK(w[d] == 0.0); }
						else            { CHECK(w[d] >= 0.0); CHECK(w[d] <= 1.0); CHECK(w[u] >= 0.0); CHECK(w[uu] == 0.0); }
					}
				}
				double g[4];
				CHECK(sc.FaceDerivative(f, sg, g) == v);
				//  Away from the kinks (a = 0, b = 0, and Koren's r = 0.4, 4)
				const int u = sg > 0 ? IB::L : IB::R, d = sg > 0 ? IB::R : IB::L, uu = sg > 0 ? IB::LL : IB::RR;
				const double a = f[u] - f[uu], b = f[d] - f[u];
				const double eps = 1e-7;
				if(std::fabs(a) < 1e-3 || std::fabs(b) < 1e-3) continue;
				if(std::fabs(a / b - 0.4) < 1e-3 || std::fabs(a / b - 4.0) < 1e-2) continue;
				for(int j = 0; j < 4; ++j)
				{
					double fp[4], fm[4];
					for(int q = 0; q < 4; ++q) { fp[q] = f[q]; fm[q] = f[q]; }
					fp[j] += eps; fm[j] -= eps;
					const double fd = (sc.FaceValue(fp, sg) - sc.FaceValue(fm, sg)) / (2 * eps);
					CAPTURE(j);
					CHECK(g[j] == doctest::Approx(fd).epsilon(1e-6).scale(1.0));
				}
			}
		}
	}
}

TEST_CASE("Convection schemes: the limiters stay finite and right at extreme differences")
{
	//  van Leer's 2ab/(a+b) and its derivatives b^2/(a+b)^2 were formed from
	//  the products: below |a|, |b| ~ 1e-162 they underflowed (Phi = 0, and
	//  the derivatives 0/0 = NaN in the Jacobian), above ~1e154 they
	//  overflowed. (Relative comparisons: scale 0.)
	for(IB::t_Scheme t : {IB::Koren, IB::VanLeer})
	{
		IB sc; sc.Build(t);
		CAPTURE(sc.SchemeName());
		for(double e : {1e-300, 1e-200, 1e-170, 1e-10, 1.0, 1e100, 1e200})
		{
			CAPTURE(e);
			for(double s : {+1.0, -1.0})
			{
				double pa = -1, pb = -1;
				//  a = b (r = 1): Phi = b for both limiters
				CHECK(sc.Limiter(s*e, s*e, &pa, &pb) == doctest::Approx(s*e).epsilon(1e-15).scale(0.0));
				CHECK(std::isfinite(pa));
				CHECK(std::isfinite(pb));
				if(t == IB::VanLeer) { CHECK(pa == doctest::Approx(0.5)); CHECK(pb == doctest::Approx(0.5)); }
				//  a = b/2: Koren min(2a, (2b + a)/3, 2b) = 5/6 b, van Leer 2ab/(a+b) = 2/3 b
				const double v = sc.Limiter(0.5*s*e, s*e, &pa, &pb);
				CHECK(v == doctest::Approx((t == IB::VanLeer ? 2.0/3.0 : 5.0/6.0)*s*e).epsilon(1e-15).scale(0.0));
				CHECK(std::isfinite(pa));
				CHECK(std::isfinite(pb));
			}
			//  The face value and its derivative on a film that varies by e
			const double f[4] = {0, e, 2*e, 3*e};
			double g[4];
			CHECK(sc.FaceDerivative(f, 1.0, g) == doctest::Approx(1.5*e).epsilon(1e-15).scale(0.0));
			for(int j = 0; j < 4; ++j) CHECK(std::isfinite(g[j]));
		}
	}
}
