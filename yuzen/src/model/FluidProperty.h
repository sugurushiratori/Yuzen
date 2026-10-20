/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  FluidProperty.h : concentration dependence of viscosity, surface tension and diffusivity
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _FLUID_PROPERTY_INCLUDED_
#define _FLUID_PROPERTY_INCLUDED_

#include <string>
#include <cmath>
#include <cassert>
using std::string;

#include "YuzenUtil.h"

#include "Input.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//--------------------------------------------------------------------------
//  Common base of the fluid properties.
//
//  The properties differ only in how they depend on the concentration, so
//  there is one derived class per kind of dependence and they are all read
//  through operator(). The caller (a physical term) does not need to know
//  which kind of dependence it is using.
//    ConstantProperty      does not depend on concentration  (density, molar mass)
//    ConcLinearProperty    linear about a reference concentration       (surface tension)
//    ConcExponentProperty  exponential about a reference concentration  (viscosity, diffusivity)
//
//  All of them return baseValue at the reference concentration baseConc. As
//  drying proceeds and the concentration rises, the viscosity grows
//  exponentially; that is the strongest nonlinearity in this solver.
//--------------------------------------------------------------------------
class C_Property
{
public:
	C_Property() : _constructed(false) {}
	C_Property(const C_Property& c) : _constructed(c._constructed) {}
	C_Property& operator=(const C_Property& c) {_constructed = c._constructed; return *this;}
	virtual ~C_Property(){}
protected:
	bool _constructed;
};

//  A property that does not depend on concentration.
class ConstantProperty : public C_Property
{
public:
	ConstantProperty() : C_Property(), _value(0) {}
	/**
	 * @brief Construct with a value but WITHOUT marking the property as constructed.
	 *
	 * The _constructed flag stays false, so operator() still asserts. Call
	 * Build() instead when the property is to be read back.
	 */
	explicit ConstantProperty(double val) : C_Property(), _value(val) {}
	/** @brief Copy-construct, carrying the value and the constructed flag over. */
	explicit ConstantProperty(const ConstantProperty& c) : C_Property(c), _value(c._value) {}
	~ConstantProperty(){}

	ConstantProperty& operator=(const ConstantProperty& c){_value = c._value; return *this;}

	/** @brief Set the value and mark the property as constructed. */
	void Build(double v)
	{
		_value = v;
		_constructed = true;
	}

	/** @brief The stored value. Asserts that Build() has been called. */
	double operator()() const
	{
		assert(_constructed);
		return _value;
	}

private:
	double _value;
};

//  A property that depends linearly on concentration.
//  value = baseValue + coef * (c - baseConc)
//  Used for the surface tension. The Marangoni number is built from coef.
class ConcLinearProperty : public C_Property
{
public:
	ConcLinearProperty() : C_Property(), _baseConc(0.5), _baseValue(0), _coef(0) {}
	ConcLinearProperty(double baseConc, double baseVal, double coef)
		: C_Property(), _baseConc(baseConc), _baseValue(baseVal), _coef(coef) {}
	ConcLinearProperty(const ConcLinearProperty& c)
		: C_Property(c), _baseConc(c._baseConc), _baseValue(c._baseValue), _coef(c._coef) {}

	~ConcLinearProperty(){}

	ConcLinearProperty& operator=(const ConcLinearProperty& c)
	{
		_constructed = c._constructed;
		_baseConc = c._baseConc;
		_baseValue = c._baseValue;
		_coef = c._coef;
		return *this;
	}

	/**
	 * @brief Set the reference point and the slope of the linear dependence.
	 * @param baseConc  reference concentration, where the property equals baseVal
	 * @param baseVal   value at the reference concentration
	 * @param coef      d(value)/d(concentration); for surface tension this is what Ma is built from
	 */
	void Build(double baseConc, double baseVal, double coef)
	{
		_baseConc = baseConc;
		_baseValue = baseVal;
		_coef = coef;
		_constructed = true;
	}

	/** @brief Evaluate baseValue + coef * (conc - baseConc). Asserts that Build() has been called. */
	double operator()(double conc) const
	{
		assert(_constructed);
		return _baseValue + _coef*(conc - _baseConc);
	}
private:
	double _baseConc;
	double _baseValue;
	double _coef;
};

//  A property that depends exponentially on concentration.
//  value = baseValue * exp(exponent * (c - baseConc))
//  Used for the viscosity and the diffusivity. The exponent of the viscosity
//  can reach about 10, so a concentration change of only 0.1 multiplies the
//  viscosity by e.
class ConcExponentProperty : public C_Property
{
public:
	ConcExponentProperty() : C_Property(),  _baseConc(0.5), _baseValue(0), _exponent(0) {}
	ConcExponentProperty(double baseConc, double baseVal, double exponent)
		: C_Property(), _baseConc(baseConc), _baseValue(baseVal), _exponent(exponent) {}
	ConcExponentProperty(const ConcExponentProperty& c)
		: C_Property(c), _baseConc(c._baseConc), _baseValue(c._baseValue), _exponent(c._exponent) {}

	~ConcExponentProperty(){}

	ConcExponentProperty& operator=(const ConcExponentProperty& c)
	{
		_constructed = c._constructed;
		_baseConc = c._baseConc;
		_baseValue = c._baseValue;
		_exponent = c._exponent;
		return *this;
	}

	/**
	 * @brief Set the reference point and the exponent of the exponential dependence.
	 * @param baseConc  reference concentration, where the property equals baseVal
	 * @param baseVal   value at the reference concentration
	 * @param exponent  rate of the exponential; for viscosity it can reach about 10,
	 *                  so a concentration change of 0.1 multiplies the viscosity by e
	 */
	void Build(double baseConc, double baseVal, double exponent)
	{
		_baseConc = baseConc;
		_baseValue = baseVal;
		_exponent = exponent;
		_constructed = true;
	}

	/** @brief Evaluate baseValue * exp(exponent * (conc - baseConc)). Asserts that Build() has been called. */
	double operator()(double conc) const
	{
		assert(_constructed);
		return _baseValue * std::exp(_exponent*(conc - _baseConc));
	}

private:
	double _baseConc;
	double _baseValue;
	double _exponent;
};




/////////////
//--------------------------------------------------------------------------
//  The full set of properties of a fluid. It holds the input [Fluid] block as
//  it stands.
//
//  What this class returns is the **dimensional** value. The non-dimensional
//  properties used in the computation are held separately by
//  NonDimensionalNumbers (normalized so that they equal 1 at the reference
//  concentration).
//--------------------------------------------------------------------------
class C_Fluid
{
public:
	C_Fluid() : _name("unknown") {}
	/**
	 * @brief Construct with a name only; every property is left unbuilt.
	 * @param name  name of the fluid, as it appears in the log
	 *
	 * Build() still has to be called before any property can be read.
	 */
	explicit C_Fluid(std::string name);
	/** @brief Copy-construct, carrying the name and every property over. */
	explicit C_Fluid(const C_Fluid& f);
	C_Fluid& operator=(const C_Fluid& c);
	~C_Fluid(){}

	/** @brief Name of the fluid, used in the log. Defaults to "unknown" until Build() runs. */
	const string& Name() const {return _name;}

	/**
	 * @brief Assemble the properties from the [Fluid] input block.
	 * @param f  the [Fluid] block, in SI units
	 */
	void Build(const Input::I_Fluid &f);

	// Property lookups. Everything here is DIMENSIONAL, in the SI units shown.
	// The non-dimensional counterparts used in the computation live on
	// NonDimensionalNumbers, normalized to 1 at the reference concentration.
	// Density and the molar masses do not depend on concentration; surface
	// tension varies linearly with it, viscosity and diffusivity exponentially.
	double Density() const                   {return _density();}		// [kg/m³]

	double SurfaceTension(double conc) const {return _surface_tension(conc);}	// [N/m]
	double Diffusivity(double conc) const    {return _diffusivity(conc);}		// [m²/s]
	double Viscosity(double conc) const      {return _viscosity(conc);}			// [Pa·s]
	double M_solvent() const                 {return _M_solvent();}	// molar mass of the solvent
	double M_resin() const                   {return _M_resin();}	// molar mass of the resin

private:
	string _name;
	ConstantProperty _density;
	ConcLinearProperty _surface_tension;
	ConcExponentProperty _viscosity;
	ConcExponentProperty _diffusivity;
	ConstantProperty _M_solvent, _M_resin;

};

#ifdef _DEF_YUZEN_
	}	// END: namespace yuzen
#endif

#endif //_FLUID_PROPERTY_INCLUDED_
