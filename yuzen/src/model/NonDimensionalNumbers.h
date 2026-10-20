/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  NonDimensionalNumbers.h : length/height/time scales and the non-dimensional numbers (As, Ma, Bo, Sp, Dp, Pe, ...)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _NONDIMENSIONAL_NUMBERS_INCLUDED_
#define _NONDIMENSIONAL_NUMBERS_INCLUDED_

#include <string>
#include <cmath>
using std::string;

#include "YuzenUtil.h"

#include "Input.h"
#include "FluidProperty.h"
#include "Evaporation.h"
#include "Mobility.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class Field_WM;
class Field_CZ;


//--------------------------------------------------------------------------
//  Scales and non-dimensional numbers.
//
//  The input is written in dimensional units (m, s, Pa s, ...), but the whole
//  computation is carried out in non-dimensional form. The conversion factors
//  and the non-dimensional numbers that appear in the governing equations are
//  all built here in one place.
//
//  How the scales are taken:
//    length   L0 = sqrt(sigma/rho g)   (the capillary length; "custom" can be
//                                       selected in the input instead)
//    height   h0 = Initial.SurfaceShape.TotalHeight
//    time     T0 = mu L0^4 / (sigma h0^3)
//    velocity U0 = L0 / T0
//  The lubrication approximation rests on the aspect ratio eps = h0/L0 being
//  small, and As = eps^-2.
//
//  When a term is disabled (Enable = false in the input), the corresponding
//  non-dimensional number is left at 0. The physical terms then simply pick up
//  a zero coefficient; no branching is needed.
//--------------------------------------------------------------------------
class NonDimensionalNumbers
{
public:
	/**
	 * @brief Set every scale and group to 1 (HStar to 0) as placeholders.
	 *
	 * These are not "term switched off" values: a group only becomes 0 when
	 * Build() finds the term disabled in the input. Build() must be called
	 * before any of these numbers means anything.
	 */
	NonDimensionalNumbers();

	/**
	 * @brief Fix every scale and non-dimensional group from the dimensional input.
	 * @param inp  input that has already been read, in SI units (m, s, Pa s, ...)
	 * @param H0   film height scale [m], normally Initial.SurfaceShape.TotalHeight
	 *
	 * Also constructs the concentration dependence of the fluid properties
	 * (viscosity, surface tension, diffusivity) and the evaporation model.
	 */
	void Build(Input& inp, double H0);

	//--------------------------------------------------------------
	//  Non-dimensional numbers. A disabled term leaves its number at 0.
	//--------------------------------------------------------------
	// The non-dimensional groups that appear in the governing equations. A term
	// switched off in the input leaves its group at 0, so the physical terms
	// only pick up a zero coefficient and need no branching.
	//   As     inverse squared aspect ratio (L0/h0)^2; strength of capillary pressure
	//   Ma     Marangoni number; strength of the concentration dependence of surface tension
	//   Bo     Bond number; gravity normal to the substrate
	//   Pgx    in-plane gravity in x, non-zero when the substrate is tilted
	//   Pgy    in-plane gravity in y
	//   Sp     spin number; centrifugal force in spin coating
	//   Dp     disjoining-pressure coefficient, fixed by the equilibrium contact angle
	//   HStar  precursor film thickness h*, the length in Pi(h) = (h*/h)^n - (h*/h)^m
	//   Pe     Peclet number; ratio of evaporation to diffusion
	//   Mr     molar mass ratio, resin over solvent
	//   Df     diffusion number
	//   Ev     evaporation number
	double As() const {return _As;}		//inverse squared aspect ratio (L0/h0)^2; strength of the capillary pressure
	double Ma() const {return _Ma;}		//Marangoni number; strength of the concentration dependence of the surface tension
	double Bo() const {return _Bo;}		//Bond number; gravity normal to the substrate
	double Pgx() const {return _Pgx;}	//in-plane gravity (x), for a tilted substrate
	double Pgy() const {return _Pgy;}	//in-plane gravity (y)
	double Sp() const {return _Sp;}		//Spin number; centrifugal force (spin coating)
	double Dp() const {return _Dp;}		//coefficient of the disjoining pressure; fixed by the equilibrium contact angle theta_e
	/** @brief The mobilities of the film, with the slip length Lambda = lambda/h0 (Physics.Slip; Mobility.h). */
	const Mobility& Mob() const {return _mob;}
	/** @brief Lambda = lambda / h0, the dimensionless slip length (0: no slip). */
	double Lambda() const {return _mob.Lambda();}
	double HStar() const {return _HStar;}	//precursor film thickness (non-dimensional). Initial.SurfaceShape.PrecursorRatio
	double Pe() const {return _Pe;}		//Peclet number; ratio of evaporation to diffusion
	double Mr() const {return _Mr;}		//molar mass ratio (resin/solvent)
	double Df() const {return _Df;}		//diffusion number
	double Ev() const {return _Ev;}		//evaporation number

	//--------------------------------------------------------------
	//  Concentration-dependent fluid properties (all non-dimensional). Each equals 1 at the reference concentration.
	//--------------------------------------------------------------
	// Concentration-dependent fluid properties, all NON-DIMENSIONAL: each is
	// normalized to 1 at the reference concentration, so these are ratios to the
	// reference value, not Pa s or m^2/s. C_Fluid returns the dimensional ones.
	// Viscosity and diffusivity vary exponentially with concentration, so drying
	// makes the viscosity climb steeply; surface tension varies linearly, and its
	// slope is what Ma is built from.
	double Viscosity(double conc) const  {return _viscosity(conc);};
	double Diffusivity(double conc) const {return _diffusivity(conc);}
	double SurfaceTension(double conc) const {return _surface_tension(conc);};

	// Evaporative flux at a point, non-dimensional. It depends on both the film
	// height and the concentration, so there is one overload per field type.
	// Zero everywhere when evaporation is disabled in the input.
	double EvapRate(Field_WM1D &ff, Xedge xe) const {return _evap(ff, xe);}
	double EvapRate(Field_WM2D &ff, Zface zf) const {return _evap(ff, zf);}
	double EvapRate(Field_CZ &ff, Zface zf) const {return _evap(ff, zf);}

	//--------------------------------------------------------------
	//  Conversion factors between dimensional and non-dimensional quantities. Used to restore dimensions on output.
	//--------------------------------------------------------------
	// Conversion factors between dimensional and non-dimensional quantities;
	// multiply a computed value by these to put dimensions back on it for output.
	//   L0 = sqrt(sigma / rho g), the capillary length, unless the input picks a custom one
	//   h0 = the film height scale
	//   T0 = mu L0^4 / (sigma h0^3),  U0 = L0 / T0
	double TimeScale() const {return _T0;}		// [s]
	double VelocityScale() const {return _U0;}	// [m/s]
	double HeightScale() const {return _h0;}	// [m]
	double ExtentScale() const {return _L0;}	// [m] in-plane direction

	/**
	 * @brief Format the scales and non-dimensional groups for the log.
	 * @return The "Scaling factors" and "Nondimensional numbers" listing that
	 *         the model writes to the Info stream, i.e. to standard output.
	 *         Redirect the run to keep it (yuzen case.toml > run.log).
	 */
	std::string ToString() const;

private:
	ConcLinearProperty _surface_tension;
	ConcExponentProperty _viscosity;
	ConcExponentProperty _diffusivity;

	C_Evaporation _evap;

	double _T0;
	double _U0;
	double _L0;
	double _h0;

	double _As;
	double _Ma;
	double _Bo;
	double _Pgx;
	double _Pgy;
	double _Sp;
	double _Dp;
	double _HStar;
	Mobility _mob;
	double _Pe;
	double _Mr;

	double _Df;
	double _Ev;

	NonDimensionalNumbers& operator=(const NonDimensionalNumbers&);
	NonDimensionalNumbers(const NonDimensionalNumbers&);

};

#ifdef _DEF_YUZEN_
	}	// END: namespace yuzen
#endif

#endif //_NONDIMENSIONAL_NUMBERS_INCLUDED_
