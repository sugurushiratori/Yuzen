/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ManufacturedSource.cpp : the source terms of the method of manufactured
*                           solutions (see the header)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "ManufacturedSource.h"
#include "MyMath.h"
#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif
using namespace yuzen::util;

namespace
{
	void Refuse(const char *what)
	{
		Global::FatalError<<"Numerics.ManufacturedSolution: "<<what
		                  <<" would enter the equations but not the manufactured sources, so the"
		                  <<" comparison would not mean anything."<<nl;
		Global::FatalError.Abort();
	}

	//  The generated h_t + div Q and concentration residual
	using std::sin; using std::cos; using std::pow; using std::exp;
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "ManufacturedSourceTerms.inc"

	//  The n-point Gauss-Legendre rule on [0, 1]: nodes and weights
	struct GaussLegendre
	{
		std::vector<double> x, w;
		explicit GaussLegendre(int n) : x(n), w(n)
		{
			for(int i = 0; i < n; i++)
			{
				double z = std::cos(M_PI * (i + 0.75) / (n + 0.5)), dp = 1.0;
				for(int it = 0; it < 100; it++)
				{
					double p0 = 1.0, p1 = z;
					for(int k = 2; k <= n; k++)
					{
						const double p2 = ((2.0*k - 1.0) * z * p1 - (k - 1.0) * p0) / k;
						p0 = p1;	p1 = p2;
					}
					dp = n * (z * p1 - p0) / (z*z - 1.0);
					const double dz = p1 / dp;
					z -= dz;
					if(std::fabs(dz) < 1e-16) break;
				}
				x[i] = 0.5 * (1.0 - z);
				w[i] = 1.0 / ((1.0 - z*z) * dp * dp);
			}
		}
	};
	//  The depth integral of the flux over the whole film, and the integrals
	//  from the substrate up to each height, panel by panel between them.
	//  The integrands are smooth (exponentials of cosines in s), so these
	//  are exact to rounding.
	const GaussLegendre &GLDepth() { static const GaussLegendre g(24); return g; }
	const GaussLegendre &GLPanel() { static const GaussLegendre g(6);  return g; }
}

void C_ManufacturedSource::Build(const Input &inp, const NonDimensionalNumbers &nd, bool twoD, bool concZ)
{
	const Input::I_Numerics::I_ManufacturedSolution &m = inp.Numerics.ManufacturedSolution;
	_enabled = m.Enable;
	if(!_enabled) return;
	_twoD = twoD;
	_concZ = concZ;

	const Input::I_Physics &ph = inp.Physics;
	if(ph.Source.Enable) Refuse("Physics.Source (the nozzle)");
	if(!concZ && !Math::IsZero(m.ConcStratAmplitude))
		Refuse("ConcStratAmplitude, a concentration varying through the film, in the well-mixed model,");
	if(concZ)
	{
		if(ph.Evaporation.Enable)        Refuse("Physics.Evaporation in ConcZ");
		if(ph.DisjoiningPressure.Enable) Refuse("Physics.DisjoiningPressure in ConcZ, which does not have it,");
	}

	//  The substrate is the manufactured one; a shape in [Substrate] would
	//  be added to it.
	const Input::I_Substrate &sb = inp.Substrate;
	bool flat = sb.FlatPlane.empty() && sb.MexicanHat.empty() && sb.SmoothConvex.empty();
	for(size_t i = 0; i < sb.RectBump.size(); i++)     flat = flat && Math::IsZero(sb.RectBump[i].Height)    && Math::IsZero(sb.RectBump[i].Base);
	for(size_t i = 0; i < sb.EllipseBump.size(); i++)  flat = flat && Math::IsZero(sb.EllipseBump[i].Height) && Math::IsZero(sb.EllipseBump[i].Base);
	for(size_t i = 0; i < sb.PeriodicWave.size(); i++) flat = flat && Math::IsZero(sb.PeriodicWave[i].Amplitude) && Math::IsZero(sb.PeriodicWave[i].Base);
	if(!flat) Refuse("A shape in [Substrate] (give it as SubstrateAmplitude instead)");

	const double twopi = 2.0 * Math::pi();
	MmsParams &p = _p;
	p.M  = m.Mean;  p.A = m.Amplitude;  p.w = m.Omega;
	p.kx = twopi / m.Wavelength[0];
	p.ky = twoD ? twopi / m.Wavelength[1] : 0.0;
	//  ConcZ with neither ConcAmplitude nor ConcStratAmplitude keeps the
	//  concentration it was given (Initial.Concentration).
	p.CM = (concZ && Math::IsZero(m.ConcAmplitude) && Math::IsZero(m.ConcStratAmplitude))
	       ? inp.Initial.Concentration.Value : m.ConcMean;
	p.CA = m.ConcAmplitude;  p.cw = m.ConcOmega;  p.cph = m.ConcPhase;
	p.ckx = twopi / m.ConcWavelength[0];
	p.cky = twoD ? twopi / m.ConcWavelength[1] : 0.0;
	p.CB  = concZ ? m.ConcStratAmplitude : 0.0;
	p.bw  = m.ConcStratOmega;
	p.bkx = twopi / m.ConcStratWavelength[0];
	p.bky = twoD ? twopi / m.ConcStratWavelength[1] : 0.0;
	p.SB = m.SubstrateMean;  p.SA = m.SubstrateAmplitude;
	p.skx = twopi / m.SubstrateWavelength[0];
	p.sky = twoD ? twopi / m.SubstrateWavelength[1] : 0.0;

	p.rho  = inp.Fluid.Density.Value;
	p.sig0 = inp.Fluid.SurfaceTension.BaseValue;
	p.sigc = inp.Fluid.SurfaceTension.Coef;
	p.c0s  = inp.Fluid.SurfaceTension.BaseConc;
	p.mu0  = inp.Fluid.Viscosity.BaseValue;
	p.mue  = inp.Fluid.Viscosity.Exponent;
	p.c0m  = inp.Fluid.Viscosity.BaseConc;
	p.D0   = inp.Fluid.Diffusivity.BaseValue;
	p.De   = inp.Fluid.Diffusivity.Exponent;
	p.c0D  = inp.Fluid.Diffusivity.BaseConc;

	p.g   = ph.VerticalGravity.Enable   ? ph.VerticalGravity.Acceleration       : 0.0;
	p.gx  = ph.HorizontalGravity.Enable ? ph.HorizontalGravity.Acceleration[0] : 0.0;
	p.gy  = (ph.HorizontalGravity.Enable && twoD) ? ph.HorizontalGravity.Acceleration[1] : 0.0;
	p.Om2 = ph.Centrifugal.Enable ? ph.Centrifugal.RotationRate * ph.Centrifugal.RotationRate : 0.0;
	p.ax  = ph.Centrifugal.RotationAxis.size() > 0 ? ph.Centrifugal.RotationAxis[0] : 0.0;
	p.ay  = ph.Centrifugal.RotationAxis.size() > 1 ? ph.Centrifugal.RotationAxis[1] : 0.0;

	_L0 = nd.ExtentScale();
	_H0 = nd.HeightScale();
	_T0 = nd.TimeScale();

	//  The disjoining pressure in pascals: the bracket of the dimensionless
	//  flux carries Dp grad Pi next to grad lap H, whose scale is
	//  sigma0 h0 / L0^3, so P = Dp sigma0 h0 / L0^2.
	p.P  = ph.DisjoiningPressure.Enable ? nd.Dp() * p.sig0 * _H0 / (_L0 * _L0) : 0.0;
	p.hs = inp.Initial.SurfaceShape.PrecursorRatio * _H0;
	p.pn = ph.DisjoiningPressure.ExponentN;
	p.pm = ph.DisjoiningPressure.ExponentM;

	p.cap  = ph.LaplacePressure.Enable ? 1.0 : 0.0;
	p.ms   = ph.SolutoCapillarity.Enable ? 1.0 : 0.0;
	//  ConcZ always convects and diffuses its concentration; the two
	//  switches are the well-mixed model's.
	p.conv = (concZ || ph.ConcField.EnableConvection) ? 1.0 : 0.0;
	p.diff = (concZ || ph.ConcField.EnableDiffusion) ? 1.0 : 0.0;
	//  With slip the residuals of the slip mobilities (mms_generate.py),
	//  where the model has them too (NonDimensionalNumbers: a positive length)
	p.lam  = nd.Mob().Slip() ? ph.Slip.Length : 0.0;

	_evap   = ph.Evaporation.Enable;
	_E0     = ph.Evaporation.Coef;
	_raoult = (ph.Evaporation.ConcFunc != "PowerN");
	_negpow   = (ph.Evaporation.SpatialFunc == "negPow");
	_Mr     = inp.Fluid.MolarMass.Resin / inp.Fluid.MolarMass.Solvent;
	_powN   = ph.Evaporation.PowerN;
	_concR  = ph.Evaporation.ConcR;
	_negN   = ph.Evaporation.NegPowExponent;
	_rmax   = ph.Evaporation.NegPowRmax;
	_hsEvap = p.hs;
	if(_evap && _negpow)
	{
		//  The profile is blended into a constant over the last five cells
		//  before NegPowRmax, which depends on the grid; keep the domain well
		//  inside it.
		const double rx = std::max(std::fabs(inp.Grid.Xmin), std::fabs(inp.Grid.Xmax));
		const double ry = twoD ? std::max(std::fabs(inp.Grid.Ymin), std::fabs(inp.Grid.Ymax)) : 0.0;
		if(std::sqrt(rx*rx + ry*ry) > 0.8 * _rmax) Refuse("A negPow profile reaching its blended rim");
	}
}

double C_ManufacturedSource::H(double x, double y, double t) const
{
	return _p.M + _p.A * cos(_p.kx * x - _p.w * t) * cos(_p.ky * Y(y));
}

double C_ManufacturedSource::Sub(double x, double y) const
{
	return _p.SB + _p.SA * cos(_p.skx * x) * cos(_p.sky * Y(y));
}

double C_ManufacturedSource::C(double x, double y, double t) const
{
	return _p.CM + _p.CA * cos(_p.ckx * x - _p.cw * t + _p.cph) * cos(_p.cky * Y(y));
}

//------------------------------------------------------------------
//  E = Coef f_conc(c) f_space(x, y) f_h(h), written from theory.md 3.7
//  rather than from Evaporation.cpp, so that the two are independent.
//------------------------------------------------------------------
double C_ManufacturedSource::Evap(double h, double c, double x, double y) const
{
	if(!_evap) return 0.0;
	double fc;
	if(_raoult)
	{
		const double chi = (c / _Mr) / (c / _Mr + (1.0 - c));	// mole fraction of resin
		fc = 1.0 - chi;
	}
	else fc = pow(std::max(1.0 - c - _concR, 0.0), _powN);

	double fs = 1.0;
	if(_negpow)
	{
		const double r = std::sqrt(x*x + Y(y)*Y(y));
		const double k = 0.5 * (1.0 - _negN) * (2.0 - _negN);
		fs = std::min(k * pow(1.0 - r / (_rmax * 1.001), -_negN), 10.0);
	}

	double fh = 1.0;
	if(_hsEvap > 0) fh = std::min((h - _hsEvap) / (19.0 * _hsEvap), 1.0);
	return _E0 * fc * fs * fh;
}

double C_ManufacturedSource::CZ(double x, double y, double zeta, double t) const
{
	return C(x, y, t) + _p.CB * cos(_p.bkx * x - _p.bw * t) * cos(_p.bky * Y(y)) * cos(M_PI * zeta);
}

double C_ManufacturedSource::CZFluxG(double x, double y, double zeta, double t) const
{
	if(!_enabled) return 0.0;
	const double yy = Y(y);
	const double h  = H(x, yy, t);
	const double ph = _p.kx * x - _p.w * t;
	const double hx = -_p.A * _p.kx * sin(ph) * cos(_p.ky * yy);
	const double hy = -_p.A * _p.ky * cos(ph) * sin(_p.ky * yy);
	const double sx = -_p.SA * _p.skx * sin(_p.skx * x) * cos(_p.sky * yy);
	const double sy = -_p.SA * _p.sky * cos(_p.skx * x) * sin(_p.sky * yy);
	const double c  = CZ(x, yy, zeta, t);
	const double D  = _p.D0 * exp(_p.De * (c - _p.c0D));
	const double cx = MmsCZcx(_p, x, yy, zeta, t);
	const double cy = MmsCZcy(_p, x, yy, zeta, t);
	const double cz = MmsCZcz(_p, x, yy, zeta, t);
	const double hA = zeta * hx + sx, hAy = zeta * hy + sy;
	return _p.diff * D * (cz / h - hA * (cx - hA / h * cz) - hAy * (cy - hAy / h * cz));
}

double C_ManufacturedSource::S(double x, double y, double t) const
{
	if(!_enabled) return 0.0;
	const double yy = Y(y);
	if(_concZ)
	{
		//  h_t + int_0^1 (1 - s) div(h K(s)) ds (mms_generate.py, build_cz)
		const GaussLegendre &g = GLDepth();
		double q = 0.0;
		for(size_t i = 0; i < g.x.size(); i++)
			q += g.w[i] * (1.0 - g.x[i]) * MmsCZDv(_p, x, yy, g.x[i], t);
		return MmsCZHt(_p, x, yy, t) + q;
	}
	if(_p.lam > 0) return MmsHeightResidualSlip(_p, x, yy, t) + Evap(H(x, yy, t), C(x, yy, t), x, yy);
	return MmsHeightResidual(_p, x, yy, t) + Evap(H(x, yy, t), C(x, yy, t), x, yy);
}
double C_ManufacturedSource::SC(double x, double y, double t) const
{
	if(!_enabled) return 0.0;
	const double yy = Y(y);
	const double h = H(x, yy, t), c = C(x, yy, t);
	if(_p.lam > 0) return MmsConcResidualSlip(_p, x, yy, t) - Evap(h, c, x, yy) * c / h;
	return MmsConcResidual(_p, x, yy, t) - Evap(h, c, x, yy) * c / h;
}

void C_ManufacturedSource::SCZColumn(double x, double y, double t, const double *zeta, int n, double *out) const
{
	if(!_enabled) { for(int k = 0; k < n; k++) out[k] = 0.0; return; }
	const double yy = Y(y);
	const double h = H(x, yy, t), ht = MmsCZHt(_p, x, yy, t);
	const GaussLegendre &g = GLPanel();
	//  From the substrate up: u = int_0^zeta Kx, v = int_0^zeta Ky, and
	//  int_0^zeta (zeta - s) Dv(s) ds = zeta F0 - F1 with F0 = int Dv and
	//  F1 = int s Dv, so that Omega = -(zeta F0 - F1 + zeta h_t) / h.
	double u = 0, v = 0, F0 = 0, F1 = 0, z0 = 0;
	for(int k = 0; k < n; k++)
	{
		const double z1 = zeta[k], dz = z1 - z0;
		for(size_t i = 0; i < g.x.size(); i++)
		{
			const double s = z0 + dz * g.x[i], w = dz * g.w[i];
			const double dv = MmsCZDv(_p, x, yy, s, t);
			u  += w * MmsCZKx(_p, x, yy, s, t);
			v  += w * MmsCZKy(_p, x, yy, s, t);
			F0 += w * dv;
			F1 += w * s * dv;
		}
		z0 = z1;
		const double Om = -(z1 * F0 - F1 + z1 * ht) / h;
		out[k] = MmsCZLocal(_p, x, yy, z1, t) + u * MmsCZcx(_p, x, yy, z1, t)
		       + v * MmsCZcy(_p, x, yy, z1, t) + (Om + z1 * ht / h) * MmsCZcz(_p, x, yy, z1, t);
	}
}
