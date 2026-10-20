/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LaplacePressure_2D.h : the common part has moved to LaplacePressureKernel<T>
*
*  All that remains here is SetLinearCoef, CalVelocity and the pinned-face
*  test. (The 2-D biharmonic operator has a cross term with no 1-D
*  counterpart, and the velocity computation writes out a third-order
*  difference directly.)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _LAPLACE_PRESSURE_2D_INCLUDED_
#define _LAPLACE_PRESSURE_2D_INCLUDED_

#include "LaplacePressureKernel.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//  The skeleton lives in LaplacePressureKernel. All that remains here is
//  what is specific to this dimension and this term: the edges of the grid,
//  the range over which the velocity is written and the points it is
//  written to, and the coefficients and the driving field.
class LaplacePressure_WM2D : public LaplacePressureKernel<WM2D_Grid>
{
public:
	LaplacePressure_WM2D() : LaplacePressureKernel<WM2D_Grid>() {}
	virtual ~LaplacePressure_WM2D(){}

	/**
	 * @brief Compute the capillary contribution to the film velocity.
	 *
	 * The third differences of the free surface h + s are written out
	 * directly in both directions, including the cross terms that carry the
	 * dx/dy anisotropy, and divided by the face-averaged film height. The
	 * results go into the x and y components of U_LaplacePressure on the Xm
	 * and Ym faces. On a pseudo 1-D grid (Ny == 1) the y part is dropped.
	 * @param field  h and s are read, U_LaplacePressure is written back
	 */
	virtual void CalVelocity(Field_WM2D &field);

protected:
	/** @brief Record the edges of the grid, in both x and y. */
	virtual void SetExtent() { _Nx = _G->nx(); _Ny = _G->ny(); }
	/**
	 * @brief Build the finite-difference coefficients of the biharmonic
	 *        operator.
	 *
	 * Besides the 1/dx^4 and 1/dy^4 terms, this includes the cross term
	 * 1/dx^2dy^2, which has no 1-D counterpart.
	 */
	virtual void SetLinearCoef();
	/**
	 * @brief Whether the face touches a pinned boundary. The two ends in y
	 *        are examined as well as those in x.
	 */
	virtual bool PinnedFace(Zface &s, tEdge tedge) const;
	/**
	 * @brief The third difference across face d, from differences of the
	 *        values (see LaplacePressureKernel::FaceGradient).
	 */
	virtual double FaceGradient(const double *v, int d) const;

private:
	LaplacePressure_WM2D(const LaplacePressure_WM2D&);
	void operator=(const LaplacePressure_WM2D&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_LAPLACE_PRESSURE_2D_INCLUDED_
