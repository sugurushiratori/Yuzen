/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  NonDimensionalNumbers.cpp : computation of the scales and the non-dimensional numbers
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#include "YuzenUtil.h"
using namespace yuzen::util;
#include "NonDimensionalNumbers.h"
#include <cmath>
#include <sstream>

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif


NonDimensionalNumbers::NonDimensionalNumbers()
{
	_T0 = 1;
	_U0 = 1;
	_L0 = 1;
	_h0 = 1;
	_Df = 1;
	_Ev = 1;

	_As = 1;
	_Ma = 1;
	_Bo = 1;
	_Pgx = 1;
	_Pgy = 1;
	_Sp = 1;
	_Dp = 1;
	_HStar = 0;
	_Pe = 1;
	_Mr = 1;
}

void NonDimensionalNumbers::Build(Input& inp, double H0)
{
	//------------------------------
	//  Take the capillary length as the in-plane length scale.
	//------------------------------

	double sig0	= inp.Fluid.SurfaceTension.BaseValue;
	double rho	= inp.Fluid.Density.Value;
	double GG = 9.80665;

	_L0 = sqrt(sig0 / (rho*GG));

	if(inp.LengthScaleType.compare("custom") == 0)
	{
		yuzen::util::Global::Info << "Length scale 'custom' is selected." << endl;
		_L0 = inp.LengthScale;
	}

	_h0 = H0;

	double eps	= _h0 / _L0;
	double mu0	= inp.Fluid.Viscosity.BaseValue;
	double sigC	= inp.Fluid.SurfaceTension.Coef;
	double Gz	= inp.Physics.VerticalGravity.Acceleration;
	double Gx	= inp.Physics.HorizontalGravity.Acceleration[0];
	double Gy	= inp.Physics.HorizontalGravity.Acceleration[1];
	double omega= inp.Physics.Centrifugal.RotationRate;
	double D0	= inp.Fluid.Diffusivity.BaseValue;
	double E0	= inp.Physics.Evaporation.Coef;


	_T0 = mu0*pow(_L0, 4.0) / (sig0 * pow(_h0, 3.0));

	_U0 = _L0 / _T0;

	_As = pow(eps, -2.0);
	_Ma = sigC / sig0;
	_Bo = rho * Gz * pow(_L0, 2.0) / sig0;

	_Pgx= rho * Gx * pow(_L0, 2.0) / sig0 / eps;
	_Pgy= rho * Gy * pow(_L0, 2.0) / sig0 / eps;

	_Sp = rho *pow(omega, 2.0) * pow(_L0, 4.0) / (sig0 * _h0);
	//------------------------------------------------------------------
	//  Disjoining pressure. The coefficient B of
	//  Pi = B[(H*/H)^n - (H*/H)^m] is fixed so that the equilibrium contact
	//  angle comes out as the theta_e given in the input.
	//
	//  At equilibrium, the difference in surface energy equals the integral of
	//  the disjoining pressure:
	//      sigma(1-cos(theta_e)) = B h* (n-m) / ((n-1)(m-1))
	//  Under the lubrication approximation theta_e is small, so
	//  1-cos(theta_e) ~ theta_e^2/2 and
	//      B = (n-1)(m-1) / (2 h*(n-m)) sigma theta_e^2
	//  Non-dimensionalizing, that is divided by the pressure scale
	//  sigma*h0/L^2, which gives
	//      Dp = (L/h0)^2 (n-1)(m-1) / (2 H*(n-m)) theta_e^2
	//  (L/h0)^2 is the square of the inverse aspect ratio, i.e. the same
	//  quantity as _As.
	//------------------------------------------------------------------
	_HStar = inp.Initial.SurfaceShape.PrecursorRatio;
	_Dp = 0;
	if(inp.Physics.DisjoiningPressure.Enable)
	{
		const double nn = inp.Physics.DisjoiningPressure.ExponentN;
		const double mm = inp.Physics.DisjoiningPressure.ExponentM;
		const double th = inp.Physics.DisjoiningPressure.ContactAngle * Math::pi() / 180.0;

		if(Math::IsZero(_HStar))
		{
			Global::FatalError<<"Physics.DisjoiningPressure requires Initial.SurfaceShape.PrecursorRatio > 0."<<nl;
			Global::FatalError.Abort();
		}

		_Dp = _As * (nn-1.0)*(mm-1.0) / (2.0*_HStar*(nn-mm)) * th*th;
	}
	_Pe = E0 * _h0 / D0;
	_Mr = inp.Fluid.MolarMass.Resin / inp.Fluid.MolarMass.Solvent;

	_Ev = E0 * _T0 / _h0;
	_Df = D0 * _T0 / pow(_L0, 2.0);


	if( ! inp.Physics.VerticalGravity.Enable) _Bo = 0;
	if( ! inp.Physics.HorizontalGravity.Enable) _Pgx = _Pgy = 0;
	if( ! inp.Physics.Centrifugal.Enable) _Sp = 0;
	if( ! inp.Physics.SolutoCapillarity.Enable) _Ma = 0;
	//  Navier slip: Lambda = lambda / h0 (0, and the kernels' no-slip code
	//  paths, unless the block enables it with a positive length)
	if(inp.Physics.Slip.Enable) _mob.Set(inp.Physics.Slip.Length / _h0);

	Input iev = inp;
	iev.Physics.Evaporation.Coef = 1.0; iev.Fluid.MolarMass.Resin = _Mr; iev.Fluid.MolarMass.Solvent = 1.0;
	if(inp.Physics.Evaporation.Enable)	_evap.Build(iev, ExtentScale());
	_surface_tension.Build(inp.Fluid.SurfaceTension.BaseConc, 1.0, _Ma);
	_viscosity.Build(inp.Fluid.Viscosity.BaseConc, 1.0, inp.Fluid.Viscosity.Exponent);
	_diffusivity.Build(inp.Fluid.Diffusivity.BaseConc, 1.0, inp.Fluid.Diffusivity.Exponent);

}


std::string NonDimensionalNumbers::ToString() const
{
	std::ostringstream oss;
	oss<<"######### Scaling factors ####################"<<std::endl;
	oss<<std::scientific;
	oss<<" Length scale [m]     = "<< _L0 << std::endl;
	oss<<" Height scale [m]     = "<< _h0 << std::endl;
	oss<<" Time   scale [s]     = "<< _T0 << std::endl;
	oss<<" Velocity scale [m/s] = "<< _U0 << std::endl;
	oss<<std::endl;
	oss<<"######### Nondimensional numbers ####################"<<std::endl;
	oss<<" As = "<< _As << std::endl;
	oss<<" Ma = "<< _Ma << std::endl;
	oss<<" Bo = "<< _Bo << std::endl;
	oss<<" Pgx= "<< _Pgx << std::endl;
	oss<<" Pgy= "<< _Pgy << std::endl;
	oss<<" Sp = "<< _Sp << std::endl;
	oss<<" Dp = "<< _Dp << std::endl;
	oss<<" Pe = "<< _Pe << std::endl;
	oss<<" Df = "<< _Df << std::endl;
	oss<<" Ev = "<< _Ev << std::endl;
	oss<<" Mr = "<< _Mr << std::endl;
	//  Only with slip, so that the listing of a run without it is as it was
	if(_mob.Slip()) oss<<" Lambda = "<< _mob.Lambda() << "  (slip length / h0)" << std::endl;
	oss<<"#####################################################"<<std::endl;
	return oss.str();
}
