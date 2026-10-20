/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  InitialCondition.h : setting of the initial film thickness distribution
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _INITIAL_CONDITION_INCLUDED_
#define _INITIAL_CONDITION_INCLUDED_


#include "LubFields.h"


#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//--------------------------------------------------------------------------
//  Write the initial film thickness distribution and the substrate shape into
//  a field.
//
//  The profile is built by superposing the shapes listed in the input under
//  [Initial.SurfaceShape] (rectangular bump, elliptic bump, plane, and so on).
//  How they are combined is decided by the CreateType of each shape
//  ("stack" adds it, "cut" carves it away).
//
//  What is written here is the **dimensional** value. After this call the model
//  divides it by TotalHeight to make it non-dimensional.
//
//  The precursor film (PrecursorRatio) is laid down everywhere as a floor.
//  When the disjoining pressure is used it is the h* of that term, so it
//  cannot be zero.
//--------------------------------------------------------------------------
class C_InitialSurface
{
public:
	/**
	 * @brief Write the initial film height and the substrate shape into a field.
	 * @param field  field to fill; its height and substrate arrays are overwritten
	 * @param inpH   [Initial.SurfaceShape] block, listing the shapes to superpose
	 * @param inpS   [Substrate] block, describing the substrate relief
	 * @param scale  in-plane length scale L0 [m]
	 *
	 * The height written here is DIMENSIONAL; the caller divides it by
	 * TotalHeight afterwards to make it non-dimensional. The precursor film
	 * (PrecursorRatio) is laid down everywhere as a floor, and it doubles as the
	 * h* of the disjoining pressure, so it cannot be zero when that term is on.
	 * One overload per field type; they do the same thing.
	 */
	static void Set(Field_WM1D& field, Input::I_Initial::I_SurfaceShape &inpH, Input::I_Substrate &inpS, double scale);
	static void Set(Field_WM2D& field, Input::I_Initial::I_SurfaceShape &inpH, Input::I_Substrate &inpS, double scale);
	static void Set(Field_CZ& field, Input::I_Initial::I_SurfaceShape &inpH, Input::I_Substrate &inpS, double scale);
private:
	static void Set(const C_Grid1D &g, ScalarField1D<> &height, ScalarField1D<> &substrate,Input::I_Initial::I_SurfaceShape &inpH, Input::I_Substrate &inpS, double scale);
	static void Set(const C_Grid2D &g, ScalarField2D<> &height, ScalarField2D<> &substrate,Input::I_Initial::I_SurfaceShape &inpH, Input::I_Substrate &inpS, double scale);
};







#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_INITIAL_CONDITION_INCLUDED_
