/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Mobility.h : the mobilities of the film, with and without Navier slip
*
*  With a Navier slip length lambda at the substrate (Physics.Slip; Lambda =
*  lambda / h0 dimensionless) the flux of a film of thickness h is
*
*      Q = -(h^3/3 + Lambda h^2) / mu  grad p  +  (h^2/2 + Lambda h) / mu  grad sigma,
*
*  so every term driven by a pressure or a body force (capillary,
*  disjoining, vertical and in-plane gravity, centrifugal) has the mobility
*
*      D(h) = h^3 + a h^2,   a = 3 Lambda       (the 1/3 is in _Magnitude),
*
*  and the shear-driven Marangoni term h^2 + 2 Lambda h (the 1/2 in
*  _Magnitude). Without slip (Lambda = 0) these are h^3 and h^2. The
*  capillary, the in-plane body forces and the Marangoni term keep their
*  own no-slip expressions and ask this class only when Slip() is true; the
*  disjoining pressure and the vertical gravity take Face() always.
*
*  The face mobility of the terms driven by the gradient of the potential
*  -- the capillary pressure, the disjoining pressure and the vertical
*  gravity, all three the same (Face) -- is the one that keeps the film
*  positive (Diez and Kondic; manual, sec:disc:positivity):
*
*      D_f = (h_2 - h_1) / (g(h_2) - g(h_1)),   g(h) = int dh / D(h),
*
*  which vanishes when either side is dry (g -> -infinity). Without slip it
*  is 2 h_1^2 h_2^2 / (h_1 + h_2). With slip
*
*      g(h) = -1/(a h) + ln((h + a)/h) / a^2 = -(x - ln(1 + x)) / a^2,  x = a/h,
*
*  evaluated as written for x >= 0.1 and by its series below (where the
*  two terms cancel); and D_f is taken from a series about the mean height
*  where the two heights differ by less than kSeriesRel of it (where the
*  difference of g cancels). See FaceDK for the errors of each.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#ifndef YUZEN_MOBILITY_H
#define YUZEN_MOBILITY_H

#include <cmath>

namespace yuzen {

class Mobility
{
public:
	//  The relative difference of the two heights of a face below which the
	//  face mobility is taken from its series (FaceDK)
	static constexpr double kSeriesRel = 1e-3;

	/** @brief Set the slip length (dimensionless, >= 0). */
	void Set(double lambda)
	{
		_lambda = lambda > 0 ? lambda : 0.0;
		_a = 3.0 * _lambda;
		_slip = _lambda > 0;
	}
	/** @brief Whether there is slip; the kernels ask this class only then. */
	bool Slip() const {return _slip;}
	double Lambda() const {return _lambda;}

	//--- the pressure mobility D(h) = h^3 + a h^2 and its derivatives ---
	//  Extended to h < 0 as the odd function h |h| (|h| + a), as h^3 is
	//  odd, so that a value below zero (the ghost -h of a Pinning edge, a
	//  Newton guess) is taken as h^3 would take it. For h >= 0 these are the
	//  expressions as written.
	double D(double h) const   {const double u = std::fabs(h); return h*u*(u + _a);}
	double dD(double h) const  {const double u = std::fabs(h); return u*(3.0*u + 2.0*_a);}
	/** @brief D(h)/h = h^2 + a h (|h| (|h| + a), even): the share a cell's
	 *         own height carries when the flux is written as (D/h) h. */
	double DperH(double h) const {const double u = std::fabs(h); return u*(u + _a);}
	double dDperH(double h) const {return h < 0 ? -(2.0*(-h) + _a) : 2.0*h + _a;}

	/**
	 * @brief The face mobility of every term driven by the gradient of the
	 *        potential -- the capillary pressure, the disjoining pressure and
	 *        the vertical gravity -- with or without slip: the
	 *        positive form of Diez and Kondic, 2 h1^2 h2^2 / (h1 + h2)
	 *        without slip (written as LaplacePressureKernel writes it, so
	 *        that the three terms take the same bits), FaceDK with it; 0
	 *        when either side is dry (0 or below, as the ghost -h of a
	 *        Pinning edge is). One mobility for the three makes the flux of
	 *        a face the mobility times the difference of the potential, so
	 *        that the dissipation is a sum of squares and the steady state of
	 *        the scheme is the minimum of its free energy (EnergyBudget.h).
	 */
	double Face(double h1, double h2) const
	{
		if(_slip) return FaceDK(h1, h2);
		if(!(h1 > 0) || !(h2 > 0)) return 0;
		return 2.0*(h1*h1)*(h2*h2) / (h1 + h2);
	}
	/** @brief dFace / dh1. */
	double Face_d1(double h1, double h2) const
	{
		if(_slip) return FaceDK_d1(h1, h2);
		if(!(h1 > 0) || !(h2 > 0)) return 0;
		return 2.0*(h1*h2*h2)*(h1+h2+h2) / ((h1 + h2)*(h1 + h2));
	}

	/** @brief The upwind face height h_f of a body force, as h_f^2 + a h_f
	 *         (times the face height again in A: D(h_f)), and its
	 *         derivative with respect to h_f; extended as DperH. */
	double Upwind(double hf) const {return DperH(hf);}
	double dUpwind(double hf) const {return dDperH(hf);}

	/** @brief The shear mobility over h of the Marangoni term: h + 2 Lambda
	 *         (times the face height in A: h^2 + 2 Lambda h). */
	double Shear(double hf) const {return hf + 2.0*_lambda;}

	/**
	 * @brief The face mobility of the capillary term (header comment),
	 *        positive-preserving, with slip.
	 *
	 * Where |h2 - h1| < kSeriesRel hm (hm the mean) it is the series
	 *
	 *     D_f = D(hm) / (1 + delta^2 K(hm) / 24),  K = D (1/D)'',
	 *
	 * whose truncation error is O((delta/hm)^4), below 1e-12 relative
	 * there; above it, delta / (g(h2) - g(h1)) loses about
	 * eps hm / delta < 1e-12 relative to the cancellation. Either side
	 * dry (0 or below): 0.
	 */
	double FaceDK(double h1, double h2) const
	{
		if(!(h1 > 0) || !(h2 > 0)) return 0;
		const double hm = 0.5*(h1 + h2), dl = h2 - h1;
		if(std::fabs(dl) < kSeriesRel*hm)
			return D(hm) / (1.0 + dl*dl*K(hm)/24.0);
		return dl / (G(h2) - G(h1));
	}
	/** @brief dFaceDK / dh1 (the derivative at the first argument). */
	double FaceDK_d1(double h1, double h2) const
	{
		if(!(h1 > 0) || !(h2 > 0)) return 0;
		const double hm = 0.5*(h1 + h2), dl = h2 - h1;
		if(std::fabs(dl) < kSeriesRel*hm)
		{
			//  F(hm, dl) = D / (1 + dl^2 K / 24); d/dh1 = F_hm / 2 - F_dl
			const double k = K(hm), q = 1.0 + dl*dl*k/24.0;
			const double Fhm = (dD(hm)*q - D(hm)*dl*dl*dK(hm)/24.0) / (q*q);
			const double Fdl = -D(hm)*(dl*k/12.0) / (q*q);
			return 0.5*Fhm - Fdl;
		}
		const double dg = G(h2) - G(h1);
		//  d/dh1 [dl / dg] = (-dg + dl / D(h1)) / dg^2
		return (-dg + dl / D(h1)) / (dg*dg);
	}

private:
	double _lambda = 0, _a = 0;
	bool _slip = false;

	//  g(h) = -(x - ln(1+x)) / a^2, x = a/h
	double G(double h) const
	{
		const double x = _a / h;
		double f;
		if(x < 0.1)
		{
			//  x - ln(1+x) = x^2/2 - x^3/3 + x^4/4 - ...; to x^18 the
			//  remainder is below 0.1^17 / 19 relative to x^2/2
			double t = x*x, s = 0;
			for(int n = 2; n <= 18; ++n)
			{
				s += ((n % 2 == 0) ? t : -t) / n;
				t *= x;
			}
			f = s;
		}
		else f = x - std::log1p(x);
		return -f / (_a*_a);
	}
	//  K = D (1/D)'' = 2 (D'/D)^2 - D''/D, and its derivative
	double K(double h) const
	{
		const double u = dD(h) / D(h);
		return 2.0*u*u - (6.0*h + 2.0*_a) / D(h);
	}
	double dK(double h) const
	{
		const double d = D(h), d1 = dD(h), d2 = 6.0*h + 2.0*_a, d3 = 6.0;
		const double u = d1 / d, du = d2 / d - u*u;
		return 4.0*u*du - (d3*d - d2*d1) / (d*d);
	}
};

}	// namespace yuzen

#endif // YUZEN_MOBILITY_H
