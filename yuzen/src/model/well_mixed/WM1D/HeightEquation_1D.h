/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  HeightEquation_1D.h : the common part has moved to the template in HeightEquationBase.h
*
*  All that remains here is the boundary condition of the velocity field.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _WM_HEIGHT_EQUATION_1D_INCLUDED_
#define _WM_HEIGHT_EQUATION_1D_INCLUDED_

#include "HeightEquationBase.h"
#include "Evaporation.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class WellMixed_HeightEquation_1D : public WellMixed_HeightEquationBase<WM1D_Traits>
{
public:
	WellMixed_HeightEquation_1D() : WellMixed_HeightEquationBase<WM1D_Traits>() {}

	/**
	 * @brief Fill the ghost layers of a velocity field at the x boundaries.
	 *
	 * Called once per physics term and once more for the mixed velocity.
	 * How the interior value is extended outwards depends on the boundary
	 * type: symmetry and pinning force zero, periodic copies from the
	 * opposite end, outlet and inlet extrapolate with zero gradient.
	 * @param vv  velocity field whose ghost nodes are overwritten
	 */
	void SetVelocityBC(VectorField1D<> &vv);

private:
	WellMixed_HeightEquation_1D(const WellMixed_HeightEquation_1D&);
	void operator=(const WellMixed_HeightEquation_1D&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_WM_HEIGHT_EQUATION_1D_INCLUDED_
