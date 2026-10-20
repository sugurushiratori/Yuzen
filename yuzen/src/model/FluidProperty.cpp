/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  FluidProperty.cpp : implementation of the concentration dependence of the fluid properties
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#include "FluidProperty.h"
#include <cmath>

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif


C_Fluid::C_Fluid(std::string name)
{
	_name = name;
}

C_Fluid::C_Fluid(const C_Fluid& f)
{
	_density			= f._density;
	_surface_tension	= f._surface_tension;
	_viscosity			= f._viscosity;
	_diffusivity		= f._diffusivity;
	_M_resin			= f._M_resin;
	_M_solvent			= f._M_solvent;
	_name				= f._name;

}

C_Fluid& C_Fluid::operator=(const C_Fluid& f)
{
	_density			= f._density;
	_surface_tension	= f._surface_tension;
	_viscosity			= f._viscosity;
	_diffusivity		= f._diffusivity;
	_M_resin			= f._M_resin;
	_M_solvent			= f._M_solvent;
	_name				= f._name;
	return *this;
}

void C_Fluid::Build(const yuzen::Input::I_Fluid &f)
{
	_name = f.Name;
	_density.Build(f.Density.Value);
	_surface_tension.Build(f.SurfaceTension.BaseConc, f.SurfaceTension.BaseValue, f.SurfaceTension.Coef);
	_viscosity.Build(f.Viscosity.BaseConc, f.Viscosity.BaseValue, f.Viscosity.Exponent);
	_diffusivity.Build(f.Diffusivity.BaseConc, f.Diffusivity.BaseValue, f.Diffusivity.Exponent);
	_M_resin.Build(f.MolarMass.Resin);
	_M_solvent.Build(f.MolarMass.Solvent);
}
