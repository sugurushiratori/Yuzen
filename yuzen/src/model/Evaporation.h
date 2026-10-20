/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Evaporation.h : evaporation rate models (Raoult / uniform / negative power)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _EVAPORATION_INCLUDED_
#define _EVAPORATION_INCLUDED_

#include "YuzenUtil.h"

#include "Input.h"
#include "MyMath.h"
#include "LubFields.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif


//--------------------------------------------------------------------------
//  Evaporation rate of the solvent.
//
//  The evaporation rate is written as a product of a concentration-dependent
//  part and a position-dependent part.
//
//      E(x, c) = Coef * f_conc(c) * f_space(x)
//
//  It enters the film height equation as -E, and the concentration equation as
//  Ev * E / h, which is the effect of the solute becoming more concentrated.
//  Because it is **divided by h**, the source term of the concentration
//  equation grows by orders of magnitude wherever the film is thin, as it is in
//  the precursor film. That is one of the reasons why the coupling becomes
//  stiff in computations that include a contact line.
//--------------------------------------------------------------------------
class C_Evaporation
{
public:
	//  Form of the concentration dependence.
	//    tc_Raoult : Raoult's law; proportional to the mole fraction of the solvent
	//    tc_PowerN : power law
	enum tConcFunc
	{
		tc_Raoult,
		tc_PowerN
	};
	//  Form of the in-plane distribution.
	//    ts_Uniform : uniform
	//    ts_NegPow  : a radial negative power of the distance to the rim of
	//                 the pool, k (1 - r/rmax')^(-NegPowExponent), blended
	//                 into a constant over 5 cells before NegPowRmax
	enum tSpatialFunc
	{
		ts_Uniform,
		ts_NegPow
	};

public:
	/** @brief Construct with evaporation disabled, so operator() returns 0 until Build() enables it. */
	C_Evaporation();

	/**
	 * @brief Fix the coefficient and the functional forms from the [Physics.Evaporation] input block.
	 * @param inp    input that has already been read
	 * @param scale  in-plane length scale L0 [m], used to non-dimensionalize the
	 *               argument of the spatial distribution
	 */
	void Build(const Input &inp, double scale);

	/**
	 * @brief Evaporative flux at a point, non-dimensional.
	 * @return E = Coef * f_conc(c) * f_space(x), POSITIVE for evaporation, i.e.
	 *         for the film losing solvent. It enters the film-height equation as
	 *         -E and the concentration equation as Ev * E / h.
	 *
	 * The flux depends on both the film height and the concentration, so there is
	 * one overload per field type; they all evaluate the same expression. When a
	 * precursor film is set (PrecursorRatio non-zero) a height-dependent
	 * correction factor multiplies the result as well, which suppresses
	 * evaporation from the precursor film. The Field_CZ overload uses the
	 * concentration at the free surface rather than a depth-averaged one. All
	 * three return 0 when evaporation is disabled in the input.
	 */
	double operator()(Field_WM1D &f, Xedge xe) const;
	double operator()(Field_WM2D &f, Zface zf) const;
	double operator()(Field_CZ &f, Zface zf) const;

	/** @brief Whether evaporation is enabled in the input; if not, operator() returns 0 everywhere. */
	bool Enabled() const {return _enabled;}

	/**
	 * @brief An evaporation rate limited to the solvent there is.
	 * @param evp  the rate (non-dimensional, as it enters h_t = -evp)
	 * @param h    the film thickness the rate is taken at
	 * @param c    the solute fraction there (the depth mean)
	 * @param dt   the time step
	 * @return evp, but at most the solvent thickness h (1 - c) over dt: a
	 *         step cannot take out more solvent than the film holds. A rate
	 *         that is not positive (no evaporation) is returned as it is.
	 *
	 * c is a mass fraction bounded by 1: at c = 1 the film is all solute
	 * and nothing is left to evaporate. Raoult's law goes to zero there by
	 * itself, but explicitly in c: a step longer than the time the film
	 * takes to dry would take it past c = 1 (and PowerN with an exponent of
	 * 0, or a ConcR < 0, does not go to zero at all). Each model takes its
	 * rate through this at every state it evaluates it at, which keeps c
	 * at or below 1 without cutting any solute away.
	 */
	static double SolventLimit(double evp, double h, double c, double dt)
	{
		if(!(evp > 0.0)) return evp;
		const double cc = (c < 0.0) ? 0.0 : ((c > 1.0) ? 1.0 : c);
		const double solvent = (h > 0.0) ? h*(1.0 - cc) : 0.0;
		if(evp*dt > solvent) evp = solvent/dt;
		return evp;
	}

	// Should a correction function that depends on the film thickness be introduced?

private:

private:
	double _M_resin;
	double _M_solvent;
	double _E0;
	double _powerN;
	double _NegPowExponent;
	double _NegPowRmax;
	C_Grid2D _G;

	double _h_ast;

	double _ConcR;

	bool _enabled;
	bool _sw_Hast_modify;

	tConcFunc _t_ConcFunc;
	tSpatialFunc _t_SpatialFunc;

	double CalChi(const double cc) const;	//Raoult's law
	double PowerN(const double cc) const;
	double NegPowFunc(const double cc, double rr) const;

	double ModifyingCoef(double hh) const;
	double Evaluate(double hh, double cc, double xx, double yy) const;

private:
	C_Evaporation(const C_Evaporation&);
	void operator=(const C_Evaporation&);
};

#ifdef _DEF_YUZEN_
	}	// END: namespace yuzen
#endif

#endif //_EVAPORATION_INCLUDED_
