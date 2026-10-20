/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Source.h : the nozzle source of the 2-D well-mixed model
*
*  The dimension-independent part is in SourceBase; this is the part that
*  walks the grid and computes each cell's covered fraction.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _SOURCE_INCLUDED_
#define _SOURCE_INCLUDED_

#include "SourceBase.h"
#include "FiniteDifferences.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif


class Source_WM2D : public SourceBase
{
public:
	Source_WM2D() : SourceBase(), _G(0) {}

	/**
	 * @brief Read the nozzle input, non-dimensionalize it and prepare the coverage field.
	 * @param inp_phys  the Physics block of the input file; the nozzle settings are read from here
	 * @param g  2-D grid on which the per-cell coverage is measured
	 * @param nonDim    the set of non-dimensional numbers, used to non-dimensionalize the nozzle settings
	 * @param ts        time-step object, from which the nozzle motion in time is taken
	 */
	void Build(Input::I_Physics &inp_phys, C_Grid2D *g, NonDimensionalNumbers *nonDim, C_TimeSteps *ts);

	/**
	 * @brief Liquid supplied to one cell at the current time.
	 * @return flow rate now times that cell's coverage; zero while the source is disabled
	 */
	double operator()(const Zface& zf) const;

	/**
	 * @brief Recompute the coverage of every cell for the present nozzle position.
	 *
	 * Call it whenever the nozzle may have moved; the height equation does so
	 * once per Newton iteration. Cells further than the nozzle radius plus two
	 * mesh widths away are set to zero without measuring the overlap.
	 */
	void updateFraction();

protected:
	C_Grid2D *_G;
	ScalarField2D<> _fraction;

	// How much of the cell is covered by the nozzle's circle (an area ratio)
	double calFraction(const Zface& zf, double xn, double yn);

private:
	Source_WM2D(const Source_WM2D&);
	void operator=(const Source_WM2D&);
};



#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_SOURCE_INCLUDED_
