/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_positivity_mobility.cpp : pins down that the face mobility is
*                                 positivity preserving
*
*  The well-mixed thickness equation uses
*
*      D = 2 h1^2 h2^2 / (h1 + h2)
*
*  as the mobility at a face (model/well_mixed/LaplacePressureKernel.h).
*  This is Eq. (23) of Diez & Kondic, J. Comput. Phys. 183 (2002) 274-306,
*  obtained by substituting D(h) = h^s with s = 3 into Eqs. (21)(22) of the
*  same paper,
*
*      D_i = (h_{i+1} - h_i) / (g_{i+1} - g_i),   g(h) = ∫ dh/D(h)
*
*  The reason it must take this form is positivity. As the thickness on one
*  side approaches 0, D goes to 0 like O(h^2), so the flux vanishes at a
*  face on the dry side and the numerical solution never steps into h < 0.
*  If this is "tidied up" into an arithmetic mean, the mobility stays finite
*  even when one side is 0, liquid gets pulled out of the dry region, and
*  negative film thicknesses appear.
*
*  The four things pinned down here are:
*    1. Eq. (23) agrees exactly with the general form of Eqs. (21)(22)
*    2. D -> 0 as one side -> 0 (the key to positivity)
*    3. that decay is O(h^2)
*    4. the arithmetic mean does not have this property (so that replacing
*       the formula with it would be noticed)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "doctest/doctest.h"

#include <cmath>

namespace {

//  Exactly the formula used by NonlinearTerm() in LaplacePressureKernel.h.
//  That one needs a Field and a Site, so here we extract just the formula
//  and pin it down. If one of the two is changed, the other must be
//  changed as well.
double MobilityPPS(double h1, double h2)
{
	return 2.0*(h1*h1)*(h2*h2) / (h1 + h2);
}

//  The general form of Eqs. (21)(22) of Diez & Kondic, with D(h) = h^s.
double G(double h, double s)
{
	//  We only use the s != 1 branch (this equation has s = 3)
	return std::pow(h, 1.0 - s) / (1.0 - s);
}

double MobilityGeneral(double h1, double h2, double s)
{
	if(h1 == h2) return std::pow(h1, s);
	return (h2 - h1) / (G(h2, s) - G(h1, s));
}

//  An example of the replacement that must not be made. Kept here for
//  comparison.
double MobilityArithmetic(double h1, double h2)
{
	return 0.5*(h1*h1*h1 + h2*h2*h2);
}

} // namespace


TEST_CASE("PPS mobility matches the general form of Diez & Kondic Eq.(21)(22)")
{
	//  Eq. (23) is the general form with s = 3 substituted in, so the two
	//  agree exactly.
	const double s = 3.0;
	const double pairs[][2] = {
		{1.0, 1.0}, {1.0, 2.0}, {0.7, 1.3}, {2.5, 0.4},
		{1.0, 1.0e-3}, {1.0, 1.0e-6}, {5.0, 5.0},
	};
	for(unsigned i = 0; i < sizeof(pairs)/sizeof(pairs[0]); ++i)
	{
		const double h1 = pairs[i][0], h2 = pairs[i][1];
		const double a = MobilityPPS(h1, h2);
		const double b = MobilityGeneral(h1, h2, s);
		CHECK(a == doctest::Approx(b).epsilon(1.0e-12));
	}
}

TEST_CASE("PPS mobility equals h^3 when the two thicknesses agree")
{
	//  For a uniform h it falls back to D(h) = h^3.
	const double hs[] = {0.25, 1.0, 3.0, 10.0};
	for(unsigned i = 0; i < sizeof(hs)/sizeof(hs[0]); ++i)
	{
		const double h = hs[i];
		CHECK(MobilityPPS(h, h) == doctest::Approx(h*h*h).epsilon(1.0e-14));
	}
}

TEST_CASE("PPS mobility vanishes as one side dries out")
{
	//  This is the key to positivity. The flux vanishes at a face on the
	//  dry side.
	const double h1 = 1.0;
	double prev = 1.0e30;
	const double h2s[] = {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-6, 1.0e-8};
	for(unsigned i = 0; i < sizeof(h2s)/sizeof(h2s[0]); ++i)
	{
		const double d = MobilityPPS(h1, h2s[i]);
		CHECK(d > 0.0);			// it stays positive
		CHECK(d < prev);		// it decreases monotonically
		prev = d;
	}
	//  On a side that is exactly 0, the mobility is 0 too.
	CHECK(MobilityPPS(1.0, 0.0) == doctest::Approx(0.0));
}

TEST_CASE("PPS mobility vanishes at second order in the thin side")
{
	//  As h2 -> 0 we have D ~ 2 h1 h2^2, so dividing h2 by 10 divides D
	//  by 100.
	const double h1 = 1.0;
	const double h2 = 1.0e-3;
	const double d1 = MobilityPPS(h1, h2);
	const double d2 = MobilityPPS(h1, 0.1*h2);
	CHECK(d1/d2 == doctest::Approx(100.0).epsilon(1.0e-3));

	//  We also check the agreement with the asymptotic form itself.
	CHECK(d1 == doctest::Approx(2.0*h1*h2*h2).epsilon(1.0e-5));
}

TEST_CASE("the arithmetic mean does not preserve positivity")
{
	//  This is here so that "tidying up" the mobility into an arithmetic
	//  mean would be noticed. Even when one side is completely dry, the
	//  arithmetic mean leaves a finite flux behind.
	const double h1 = 1.0;
	CHECK(MobilityArithmetic(h1, 0.0) > 0.4);		// h1^3/2 = 0.5
	CHECK(MobilityPPS(h1, 0.0) == doctest::Approx(0.0));

	//  Even on the thin side the two are orders of magnitude apart.
	const double h2 = 1.0e-4;
	CHECK(MobilityArithmetic(h1, h2) > 1.0e-1);
	CHECK(MobilityPPS(h1, h2) < 1.0e-7);
}


//==========================================================================
//  With Navier slip (Physics.Slip; model/Mobility.h): D(h) = h^3 + a h^2,
//  a = 3 Lambda, and the face mobility of the capillary term is the general
//  form above with g(h) = -(x - ln(1 + x)) / a^2, x = a / h. It is computed
//  as written except where it cancels: g from its series for x < 0.1, D_f
//  from its series about the mean height for |h2 - h1| < 1e-3 of it. What
//  is pinned down:
//    5. D_f against a reference without cancellation, to 1e-12, across both
//       switches, for slip lengths from 1e-6 to 3 of the film scale
//    6. it is continuous where it switches
//    7. D_f -> 0 as one side -> 0, and D_f = 0 on a dry side; with slip
//       the decay is O(h) (D ~ a h^2 near the wall), not O(h^2)
//    8. Lambda -> 0 gives back 2 h1^2 h2^2 / (h1 + h2)
//    9. its derivative (the Jacobian) against the finite difference
//==========================================================================
#include "Mobility.h"
#include <algorithm>
#include <vector>

namespace {

typedef long double LD;

//  The 20-point Gauss-Legendre rule in long double, by Newton on P_20
struct GL20
{
	LD x[20], w[20];
	GL20()
	{
		const int n = 20;
		for(int i = 0; i < n; ++i)
		{
			LD z = std::cos(3.14159265358979323846264338327950288L * (i + 0.75L) / (n + 0.5L)), dp = 0;
			for(int it = 0; it < 100; ++it)
			{
				LD p0 = 1, p1 = z;
				for(int k = 2; k <= n; ++k) { LD p2 = ((2*k - 1)*z*p1 - (k - 1)*p0) / k; p0 = p1; p1 = p2; }
				dp = n*(z*p1 - p0) / (z*z - 1);
				const LD dz = p1 / dp;
				z -= dz;
				if(std::fabs(dz) < 1e-30L) break;
			}
			x[i] = z; w[i] = 2 / ((1 - z*z)*dp*dp);
		}
	}
};
const GL20 &Rule() { static const GL20 r; return r; }

//  The reference D_f: delta over the integral of 1/D, the integral by
//  Gauss-Legendre over [h1, h2] when the interval is short against the
//  heights (no cancellation; the integrand is smooth there), from g in
//  long double otherwise (the series of x - ln(1+x) for small x)
LD RefFace(LD h1, LD h2, LD a)
{
	if(h1 == h2) return h1*h1*(h1 + a);
	const LD lo = std::min(h1, h2), hi = std::max(h1, h2);
	LD I;
	if((hi - lo) < 0.2L*lo)
	{
		const LD m = (hi + lo)/2, r = (hi - lo)/2;
		I = 0;
		for(int i = 0; i < 20; ++i) { const LD h = m + r*Rule().x[i]; I += Rule().w[i] / (h*h*(h + a)); }
		I *= r;
	}
	else
	{
		auto g = [a](LD h) {
			const LD x = a / h;
			LD f;
			if(x < 0.1L) { LD t = x*x, s = 0; for(int n = 2; n <= 40; ++n) { s += ((n % 2 == 0) ? t : -t) / n; t *= x; } f = s; }
			else f = x - std::log1p(x);
			return -f / (a*a);
		};
		I = g(hi) - g(lo);
	}
	return (hi - lo) / I;
}

yuzen::Mobility Mob(double lambda) { yuzen::Mobility m; m.Set(lambda); return m; }

const double kLambdas[] = {1e-6, 1e-3, 3e-2, 0.3, 3.0};

} // namespace


TEST_CASE("slip: the face mobility against a reference without cancellation")
{
	const double hs[] = {1e-3, 0.1, 1.0, 7.0};
	//  relative separations (h2 - h1) / h1: none, inside the series of D_f,
	//  at its switch on both sides, and well apart
	const double rel[] = {0.0, 1e-9, 1e-6, 1e-4, 9.99e-4, 1.001e-3, 1e-2, 0.15, 0.5, 3.0, 100.0};
	double worst = 0;
	for(double lam : kLambdas)
	{
		const yuzen::Mobility m = Mob(lam);
		for(double h1 : hs)
			for(double r : rel)
			{
				const double h2 = h1*(1.0 + r);
				const LD ref = RefFace(h1, h2, 3.0L*lam);
				const double got = m.FaceDK(h1, h2), sym = m.FaceDK(h2, h1);
				const double e = (double)(std::fabs(got - ref) / ref);
				worst = std::max(worst, e);
				CAPTURE(lam); CAPTURE(h1); CAPTURE(r);
				CHECK(e < 1e-12);
				CHECK(sym == doctest::Approx(got).epsilon(1e-15));
			}
	}
	MESSAGE("slip face mobility: largest relative error against the reference ", worst);
}

TEST_CASE("slip: the face mobility is continuous where it switches to its series")
{
	const double k = yuzen::Mobility::kSeriesRel;
	for(double lam : kLambdas)
	{
		const yuzen::Mobility m = Mob(lam);
		for(double hm : {1e-3, 0.5, 4.0})
		{
			//  |dl| just below and just above k hm about the same mean
			const double d0 = k*hm*(1.0 - 1e-9), d1 = k*hm*(1.0 + 1e-9);
			const double a = m.FaceDK(hm - d0/2, hm + d0/2), b = m.FaceDK(hm - d1/2, hm + d1/2);
			CAPTURE(lam); CAPTURE(hm);
			CHECK(std::fabs(a - b) / a < 1e-12);
			const double da = m.FaceDK_d1(hm - d0/2, hm + d0/2), db = m.FaceDK_d1(hm - d1/2, hm + d1/2);
			CHECK(std::fabs(da - db) / std::fabs(da) < 1e-9);
		}
	}
}

TEST_CASE("slip: the face mobility vanishes on a dry side, at first order")
{
	for(double lam : kLambdas)
	{
		const yuzen::Mobility m = Mob(lam);
		CAPTURE(lam);
		CHECK(m.FaceDK(1.0, 0.0) == 0.0);
		CHECK(m.FaceDK(0.0, 1.0) == 0.0);
		CHECK(m.FaceDK(1.0, -1e-3) == 0.0);
		CHECK(m.FaceDK_d1(0.0, 1.0) == 0.0);
		double prev = 1e30;
		for(double h2 : {1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8})
		{
			const double d = m.FaceDK(1.0, h2);
			CHECK(d > 0.0);
			CHECK(d < prev);
			prev = d;
		}
		//  Near the dry side g ~ -1/(a h): D_f ~ a h2 (h1 - h2) / (1 + ...),
		//  so dividing h2 by 10 divides D_f by 10 once a h2 dominates h2^2,
		//  i.e. for h2 well below a
		const double a = 3.0*lam, h2 = 1e-4*std::min(a, 1.0);
		CHECK(m.FaceDK(1.0, h2) / m.FaceDK(1.0, 0.1*h2) == doctest::Approx(10.0).epsilon(2e-3));
	}
}

TEST_CASE("slip: a vanishing slip length gives back the no-slip face mobility")
{
	const yuzen::Mobility m = Mob(1e-10);
	const double pairs[][2] = {{1.0, 2.0}, {0.7, 1.3}, {2.5, 0.4}, {1.0, 1.0e-3}, {1.0, 1.0 + 1e-5}, {5.0, 5.0}};
	for(const auto &p : pairs)
		CHECK(m.FaceDK(p[0], p[1]) == doctest::Approx(MobilityPPS(p[0], p[1])).epsilon(1e-7));
}

TEST_CASE("slip: the derivative of the face mobility against the finite difference")
{
	const double rel[] = {-0.6, -1e-2, -5e-4, -1e-6, 0.0, 1e-6, 5e-4, 1e-2, 0.6, 4.0};
	double worst = 0;
	for(double lam : kLambdas)
	{
		const yuzen::Mobility m = Mob(lam);
		for(double h2 : {1e-3, 0.2, 3.0})
			for(double r : rel)
			{
				const double h1 = h2*(1.0 + r);
				//  a step that keeps both points on the same side of the switch
				const double eps = 1e-4*h1*std::max(std::min(std::fabs(r), 1.0), 1e-3);
				const double fd = (m.FaceDK(h1 + eps, h2) - m.FaceDK(h1 - eps, h2)) / (2*eps);
				const double an = m.FaceDK_d1(h1, h2);
				const double e = std::fabs(an - fd) / std::max(std::fabs(fd), m.D(h2)/h2);
				worst = std::max(worst, e);
				CAPTURE(lam); CAPTURE(h2); CAPTURE(r);
				CHECK(e < 1e-6);
			}
	}
	MESSAGE("slip face mobility: largest error of the derivative against the finite difference ", worst);
}

TEST_CASE("the face mobility of the potential-driven terms vanishes against a Pinning ghost")
{
	//  A Pinning edge holds the film at zero with the ghost -h. The
	//  disjoining pressure and the vertical gravity take Mobility::Face,
	//  which must give that face nothing, with and without slip (a mean of
	//  h^3 did it by symmetry; one of D with slip did not, VolErr 2.5e-6 in
	//  a closed box, tools/analytic.py wall_balance_slip); D itself is
	//  extended to h < 0 as an odd function, as h^3 is.
	yuzen::Mobility none;
	none.Set(0.0);
	for(double h : {1e-3, 0.3, 2.0})
	{
		CHECK(none.Face(h, -h) == 0.0);
		CHECK(none.Face(-h, h) == 0.0);
		CHECK(none.Face_d1(h, -h) == 0.0);
		//  the same bits as the capillary kernel's 2 h1^2 h2^2 / (h1 + h2)
		CHECK(none.Face(h, 1.7*h) == MobilityPPS(h, 1.7*h));
	}
	for(double lam : kLambdas)
	{
		const yuzen::Mobility m = Mob(lam);
		for(double h : {1e-3, 0.3, 2.0})
		{
			CAPTURE(lam); CAPTURE(h);
			CHECK(m.Face(h, -h) == 0.0);
			CHECK(m.DperH(-h) == m.DperH(h));
			CHECK(m.dD(-h) == m.dD(h));
			CHECK(m.D(-h) == -m.D(h));
			CHECK(m.dDperH(-h) == -m.dDperH(h));
			CHECK(m.D(h) == h*h*(h + 3.0*lam));
			CHECK(m.Face(h, 1.3*h) == m.FaceDK(h, 1.3*h));
		}
	}
}
