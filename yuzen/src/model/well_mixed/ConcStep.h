/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ConcStep.h : what the film-height equation hands the concentration
*               equation of the well-mixed model for one solve
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _WM_CONC_STEP_INCLUDED_
#define _WM_CONC_STEP_INCLUDED_

#include <vector>
#include "Vector.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//--------------------------------------------------------------------------
//  The film side of one solve of the well-mixed concentration equation.
//
//  The equation is the balance of the solute h c, cell by cell,
//
//      h1 c1 - h0 c0 = -dt sum_f [theta G_f(c1) + (1 - theta) G_f(c0)]
//                      + dt S c_S
//                      - dt (theta r1 c1 + (1 - theta) r0 c0),
//
//  G_f = F_f c_f - K_f (c_hi - c_lo) the solute flux through face f, F_f
//  the volume flux of the film-height equation through that face, S the
//  liquid the nozzle (and a manufactured source) supplies and r the
//  sideways outflow of the 1-D centrifugal term. c_S is the concentration
//  the supplied liquid arrives at: that of the film where it lands,
//  theta c1 + (1 - theta) c0 (Physics.Source.ConcType = "local", the
//  default), or the nozzle's Physics.Source.Conc ("supplied", srcSupplied;
//  the liquid is then mixed through the thickness of the column it lands
//  on). Evaporation takes solvent only and does not appear: it is in h1.
//  When the step is completed
//  (commit) the fluxes are those of the film-height equation's own step,
//
//      h1 = h0 - dt sum_f [theta F_f^1 + (1 - theta) F_f^0]
//              - dt (theta r1 + (1 - theta) r0) - dt (E - S)
//
//  cell by cell (HeightEquationBase::ConsistentFluxes), and the two are one
//  finite-volume step: the solute changes by the fluxes through the edges
//  and the sources only, and a uniform c stays uniform (unless liquid of
//  another concentration is supplied). In the Newton
//  passes before, h1 is the film of the pass and F^1 its fluxes, which
//  give that film only once the iteration has converged.
//
//  Every array is indexed by row (ghosts included). F0 and F1 point to
//  T::nAxis arrays each; F[a][i] is the flux through the face between row
//  i and its neighbor on the high side along axis a.
//--------------------------------------------------------------------------
struct WMConcStep
{
	const Vector<double> *h0;			// the film at the start of the step
	const Vector<double> *h1;			// the film at its end
	const std::vector<double> *F0;		// the face volume fluxes at the start (T::nAxis arrays)
	const std::vector<double> *F1;		// and at the film of the end
	const Vector<double> *r0, *r1;		// the sideways outflow, start and end
	const Vector<double> *src;			// the liquid supplied over the step, a rate
	const Vector<double> *evp;			// the evaporation over the step, a rate
	bool srcSupplied;					// the supplied liquid arrives at srcConc, not at the local c
	double srcConc;						// its concentration then (Physics.Source.Conc)
	bool commit;						// h1 is the film of the fluxes: the step is committed

	WMConcStep() : h0(0), h1(0), F0(0), F1(0), r0(0), r1(0), src(0), evp(0),
	               srcSupplied(false), srcConc(0), commit(false) {}
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_WM_CONC_STEP_INCLUDED_
