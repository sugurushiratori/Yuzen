/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  WellMixedTraits.h : the part of the list of "what differs between
*                      dimensions" that does depend on the physics-term
*                      classes
*
*  Everything to do with the grid lives in GridTraits.h; this file merely
*  derives from that and adds the typedefs of the physics-term and equation
*  classes.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _WELL_MIXED_TRAITS_INCLUDED_
#define _WELL_MIXED_TRAITS_INCLUDED_

#include "GridTraits.h"

#include "VerticalGravityKernel.h"
#include "DisjoiningPressureKernel.h"
#include "SolutoCapillarityKernel.h"
#include "CentrifugalKernel.h"
#include "ParallelGravityKernel.h"

#include "LaplacePressure_1D.h"
#include "Source_1D.h"

#include "LaplacePressure_2D.h"
#include "Source.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class WellMixed_ConcEquation_1D;
class WellMixed_ConcEquation_2D;
class WellMixed_HeightEquation_1D;
class WellMixed_HeightEquation_2D;

/**
 * @brief 1-D traits: the grid traits, plus the concrete physics-term and
 *        equation classes used by the 1-D well-mixed model.
 */
struct WM1D_Traits : public WM1D_Grid
{
	typedef PhysicsTermBase<WM1D_Grid>  Physics;
	typedef WellMixed_ConcEquation_1D   ConcEq;
	typedef VectorField1D<>            VelField;

	typedef LaplacePressure_WM1D                LaplacePressure;
	typedef VerticalGravityKernel<WM1D_Grid>    VerticalGravity;
	typedef ParallelGravityKernel<WM1D_Grid>    ParallelGravity;
	typedef CentrifugalKernel<WM1D_Grid>        Centrifugal;
	typedef SolutoCapillarityKernel<WM1D_Grid>  SolutoCapillarity;
	typedef DisjoiningPressureKernel<WM1D_Grid> DisjoiningPressure;
	typedef Source_WM1D                         Source;

	typedef WellMixed_HeightEquation_1D HeightEq;
	typedef Field_WM1D_Out              FieldOut;
};

/**
 * @brief 2-D counterpart of WM1D_Traits; the members mirror it one for one.
 */
struct WM2D_Traits : public WM2D_Grid
{
	typedef PhysicsTermBase<WM2D_Grid>  Physics;
	typedef WellMixed_ConcEquation_2D   ConcEq;
	typedef VectorField2D<>            VelField;

	typedef LaplacePressure_WM2D                LaplacePressure;
	typedef VerticalGravityKernel<WM2D_Grid>    VerticalGravity;
	typedef ParallelGravityKernel<WM2D_Grid>    ParallelGravity;
	typedef CentrifugalKernel<WM2D_Grid>        Centrifugal;
	typedef SolutoCapillarityKernel<WM2D_Grid>  SolutoCapillarity;
	typedef DisjoiningPressureKernel<WM2D_Grid> DisjoiningPressure;
	typedef Source_WM2D                         Source;

	typedef WellMixed_HeightEquation_2D HeightEq;
	typedef Field_WM2D_Out              FieldOut;
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_WELL_MIXED_TRAITS_INCLUDED_
