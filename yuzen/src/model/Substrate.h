/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Substrate.h : generation of the substrate shape
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _DEF_YUZEN_SUBSTRATE_INCLUDED_
#define _DEF_YUZEN_SUBSTRATE_INCLUDED_

#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "grid/YuzenGrid.h"
#include "Input.h"
//#include "Grids.h"


#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

using namespace yuzen::grid;

class C_Grid1D;
class C_Grid2D;

/////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------
//  Generation of the substrate shape (the relief of the bare surface).
//
//  The film thickness h is measured from the substrate, so when the substrate
//  has relief the free surface sits at h + s. What the capillary pressure
//  responds to is that h + s, so the substrate shape feeds into the physical
//  terms (LaplacePressureKernel looks at s).
//--------------------------------------------------------------------------
class C_Substrate
{
public:
	/**
	 * @brief Write the substrate profile s over a grid.
	 * @param g        grid to evaluate the profile on
	 * @param profile  output field, overwritten
	 * @param subs     [Substrate] input block, describing the relief
	 * @param scale    in-plane length scale L0 [m], used to non-dimensionalize the
	 *                 in-plane positions in the input
	 *
	 * The free surface sits at h + s, and that is what the capillary pressure
	 * responds to, so the substrate shape feeds into the physical terms.
	 * One overload per grid dimension; they do the same thing.
	 */
	static void MakeSubstrateProfile(const C_Grid1D &g, ScalarField1D<> &profile, const Input::I_Substrate &subs, double scale);
	static void MakeSubstrateProfile(const C_Grid2D &g, ScalarField2D<> &profile, const Input::I_Substrate &subs, double scale);
};


#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_DEF_YUZEN_SUBSTRATE_INCLUDED_
