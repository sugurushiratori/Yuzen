/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  HeightEquation_2D.h : the common part has moved to the template in HeightEquationBase.h
*
*  All that remains here is the boundary condition of the velocity field
*  and the number of grid points in Y.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _WM_HEIGHT_EQUATION_2D_INCLUDED_
#define _WM_HEIGHT_EQUATION_2D_INCLUDED_

#include "HeightEquationBase.h"
#include "Evaporation.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class WellMixed_HeightEquation_2D : public WellMixed_HeightEquationBase<WM2D_Traits>
{
public:
	WellMixed_HeightEquation_2D() : WellMixed_HeightEquationBase<WM2D_Traits>(), _Ny(1) {}

	/**
	 * @brief Fill the ghost layers of a velocity field at the four boundaries.
	 *
	 * Called once per physics term and once more for the mixed velocity.
	 * The x and y edges are swept separately, and how the interior value is
	 * extended outwards depends on the boundary type: symmetry and pinning
	 * force the normal component to zero while copying the tangential one,
	 * periodic copies from the opposite end, outlet and inlet extrapolate
	 * with zero gradient.
	 * @param vv  velocity field whose ghost nodes are overwritten
	 */
	void SetVelocityBC(VectorField2D<> &vv);

protected:
	//  The 2-D case also needs the number of grid points in y (the base
	//  class holds only _Nx)
	virtual void BuildGridExtent() { _Ny = G->ny(); }

	int _Ny;

private:
	WellMixed_HeightEquation_2D(const WellMixed_HeightEquation_2D&);
	void operator=(const WellMixed_HeightEquation_2D&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_WM_HEIGHT_EQUATION_2D_INCLUDED_
