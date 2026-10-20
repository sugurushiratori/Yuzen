/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LaplacePressure_1D.h : the common part has moved to LaplacePressureKernel<T>
*
*  All that remains here is SetLinearCoef, CalVelocity and the pinned-face
*  test. (The 2-D biharmonic operator has a cross term with no 1-D
*  counterpart, and the velocity computation writes out a third-order
*  difference directly.)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _LAPLACE_PRESSURE_1D_INCLUDED_
#define _LAPLACE_PRESSURE_1D_INCLUDED_

#include "LaplacePressureKernel.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//  The skeleton lives in LaplacePressureKernel. All that remains here is
//  what is specific to this dimension and this term: the edges of the grid,
//  the range over which the velocity is written and the points it is
//  written to, and the coefficients and the driving field.
class LaplacePressure_WM1D : public LaplacePressureKernel<WM1D_Grid>
{
public:
	LaplacePressure_WM1D() : LaplacePressureKernel<WM1D_Grid>() {}
	virtual ~LaplacePressure_WM1D(){}

	/**
	 * @brief Compute the capillary contribution to the film velocity.
	 *
	 * The third difference of the free surface h + s is written out directly
	 * and divided by the face-averaged film height; the result goes into the
	 * x component of U_LaplacePressure on every interior Xm face.
	 * @param field  h and s are read, U_LaplacePressure is written back
	 */
	virtual void CalVelocity(Field_WM1D &field);

protected:
	/**
	 * @brief Record the edges of the grid; in 1-D only the x direction
	 *        has any.
	 */
	virtual void SetExtent() { _Nx = _G->nx(); }
	/**
	 * @brief Build the finite-difference coefficients of the biharmonic
	 *        operator (1/dx^4) and the coefficients CV used for the velocity.
	 *
	 * Unlike the 2-D case, no cross term (1/dx^2dy^2) appears.
	 */
	virtual void SetLinearCoef();
	/**
	 * @brief Whether the face touches a pinned boundary. In 1-D only the two
	 *        ends in x are examined.
	 */
	virtual bool PinnedFace(Xedge &s, tEdge tedge) const;
	/**
	 * @brief The third difference across face d, from differences of the
	 *        values (see LaplacePressureKernel::FaceGradient).
	 */
	virtual double FaceGradient(const double *v, int d) const;

private:
	LaplacePressure_WM1D(const LaplacePressure_WM1D&);
	void operator=(const LaplacePressure_WM1D&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_LAPLACE_PRESSURE_1D_INCLUDED_
