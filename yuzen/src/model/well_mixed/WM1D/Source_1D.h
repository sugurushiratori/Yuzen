/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Source_1D.h : the 1-D version of the liquid source (the nozzle). New.
*
*  The dimension-independent part (reading the input, the nozzle position,
*  the ramp-up of the flow rate) lives in SourceBase; all this file does is
*  walk the grid and compute each cell's covered fraction.
*
*  The nozzle position is the 2-D expression with y dropped, that is
*      xn = ArmX + ArmLength * cos(theta)
*  With ArmLength = 0 it is fixed at ArmX. The input block is shared with
*  the 2-D version.
*
*  Where 2-D takes the covered fraction to be "the fraction of the cell's
*  area that falls within the nozzle's circle", 1-D takes it to be "the
*  fraction of the cell's width that falls within |x - xn| <= Radius". The
*  integral simply loses one dimension; the idea is the same.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _SOURCE_1D_INCLUDED_
#define _SOURCE_1D_INCLUDED_

#include "SourceBase.h"
#include "FiniteDifferences.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class Source_WM1D : public SourceBase
{
public:
	Source_WM1D() : SourceBase(), _G(0) {}

	/**
	 * @brief Read the nozzle input, non-dimensionalize it and prepare the coverage field.
	 * @param inp_phys  the Physics block of the input file; the nozzle settings are read from here
	 * @param g  1-D grid on which the per-cell coverage is measured
	 * @param nonDim    the set of non-dimensional numbers, used to non-dimensionalize the nozzle settings
	 * @param ts        time-step object, from which the nozzle motion in time is taken
	 */
	void Build(Input::I_Physics &inp_phys, C_Grid1D *g, NonDimensionalNumbers *nonDim, C_TimeSteps *ts);

	/**
	 * @brief Liquid supplied to one cell at the current time.
	 * @return flow rate now times that cell's coverage; zero while the source is disabled
	 */
	double operator()(const Xedge& xe) const;

	/**
	 * @brief Recompute the coverage of every cell for the present nozzle position.
	 *
	 * Call it whenever the nozzle may have moved; the height equation does so
	 * once per Newton iteration. Skipped entirely while the source is disabled.
	 */
	void updateFraction();

protected:
	C_Grid1D *_G;
	ScalarField1D<> _fraction;

	// The fraction of the cell's width that falls within the interval the
	// nozzle covers
	double calFraction(const Xedge& xe, double xn);

private:
	Source_WM1D(const Source_WM1D&);
	void operator=(const Source_WM1D&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_SOURCE_1D_INCLUDED_
