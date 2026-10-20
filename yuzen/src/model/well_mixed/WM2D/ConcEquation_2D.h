/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ConcEquation_2D.h : the concentration equation of the 2-D well-mixed model
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _WM_CONC_EQUATION_2D_INCLUDED_
#define _WM_CONC_EQUATION_2D_INCLUDED_

#include "ConcEquationBase.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//--------------------------------------------------------------------------
//  The concentration equation of the 2-D well-mixed model: the
//  finite-volume balance of the solute h c, solved for the change of c
//  over the step with the face volume fluxes of the film-height equation
//  (ConcEquationBase.h, ConcStep.h). Everything is in the base; the grid
//  traits WM2D_Grid supply the dimension.
//--------------------------------------------------------------------------
class WellMixed_ConcEquation_2D : public WellMixed_ConcEquationBase<WM2D_Grid>
{
public:
	WellMixed_ConcEquation_2D() : WellMixed_ConcEquationBase<WM2D_Grid>() {}

private:
	WellMixed_ConcEquation_2D(const WellMixed_ConcEquation_2D&);
	void operator=(const WellMixed_ConcEquation_2D&);
};


#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_WM_CONC_EQUATION_2D_INCLUDED_
