/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  FieldInfo.h : a visitor that takes the aggregates of a field (volume, mean
*                concentration, maxima and minima)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _FIELD_INFO_INCLUDED_
#define _FIELD_INFO_INCLUDED_


#include <list>

#include "Grids.h"
#include "Input.h"
#include "FieldVisitor.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class Field_WM1D;
class Field_WM2D;
class Field_CZ;

//--------------------------------------------------------------------------
//  A visitor that takes the aggregates of a field.
//
//  Calling field.Accept(info) sweeps the field once and fills in every value
//  below. They are the values printed to History.dat and the ones used for the
//  conservation check (VolErr / ResinVolErr). Only the interior of the
//  computational domain is aggregated; the ghost layers are not included.
//--------------------------------------------------------------------------
class FieldInfo : public FieldVisitor
{
public:
	/** @brief Builds an empty summary, with every aggregate set to zero. */
	FieldInfo();
	/**
	 * @brief Copies an already computed summary. Declared explicit, so a
	 *        FieldInfo is never copied into a by-value parameter by accident.
	 */
	explicit FieldInfo(const FieldInfo&);
	FieldInfo& operator=(const FieldInfo&);

	/**
	 * @brief Sweeps the field once and fills in every aggregate below.
	 *
	 * Called through field.Accept(info); there is one overload per kind of
	 * field. Each call recomputes the aggregates from scratch, so the results
	 * describe the last field visited and not an accumulation over several.
	 */
	//  Sweeps the field and aggregates it. One overload per kind of field.
	virtual void Visit(Field_WM1D &);
	virtual void Visit(Field_WM2D &);
	virtual void Visit(Field_CZ &);

	//--------------------------------------------------------------
	// Aggregates of the last field visited. All non-dimensional, and only
	// meaningful once Visit has run; before that they are zero.
	//--------------------------------------------------------------
	//  Aggregates (all non-dimensional). Read them after calling Visit.
	//--------------------------------------------------------------
	/** @brief Volume of liquid in the domain, the integral of h. */
	double FluidVol() const {return _fluidVolume;}	// volume of liquid (int h)
	/** @brief Volume of resin in the domain, the integral of h*c. */
	double ResinVol() const {return _resinVolume;}	// volume of resin (int h c)
	/** @brief Concentration averaged over the volume, ResinVol() / FluidVol(). */
	double AveConc() const {return _averageConc;}	// volume-averaged concentration
	/** @brief Largest and smallest film thickness h over the domain. */
	double MaxH() const {return _maxH;}				// film thickness h
	double MinH() const {return _minH;}
	/** @brief Largest and smallest concentration c over the domain. */
	double MaxC() const {return _maxC;}				// concentration c
	double MinC() const {return _minC;}
	/** @brief Largest and smallest free-surface height h+s over the domain. */
	double MaxHS() const {return _maxHS;}			// free-surface height h+s
	double MinHS() const {return _minHS;}
	//  (What crosses the open edges is not taken from a field: the model
	//  sums the face fluxes of each step as the equations took them,
	//  HeightEquationBase::EdgeOutflowVolume.)

private:
	double _fluidVolume;
	double _resinVolume;
	double _averageConc;
	double _maxH, _minH;
	double _maxC, _minC;
	double _maxHS, _minHS;
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_FIELD_INFO_INCLUDED_
