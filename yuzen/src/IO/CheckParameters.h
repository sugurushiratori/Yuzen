/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CheckParameters.h : consistency check of the input values
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _DEF_CHECK_PARAMETERS_
#define _DEF_CHECK_PARAMETERS_

#include "YuzenUtil.h"
#include "Input.h"


#ifdef _DEF_YUZEN_
namespace yuzen{
#endif


//--------------------------------------------------------------------------
//  Consistency check of the input values.
//
//  The range of each individual value is checked by Input::Read as it reads.
//  What is looked at here are the **combinations between values**: the kind of
//  contradiction where each value is valid on its own but they cannot be used
//  together -- and the values with which the equations mean nothing (a
//  viscosity, density, surface tension or height scale of zero or below, a
//  time step or an output interval of zero, a shape that divides by a zero
//  width), which the manual's key tables call positive or non-zero and
//  tools/doc_keys.py therefore does not hold Input.cpp to.
//--------------------------------------------------------------------------
class CheckParameters
{
public:
	CheckParameters(){}

	/**
	 * @brief Checks the input parameters against each other.
	 *
	 * Per-value ranges have already been checked while reading; what is looked
	 * for here are combinations that are each valid on their own but cannot be
	 * used together.
	 *
	 * @param inp  the parsed input; it is not modified
	 * @return true if the settings are consistent. On false a warning naming
	 *         the conflicting parameters has been printed, and the caller
	 *         (C_Yuzen) stops the run there.
	 */
	//  Prints a warning and returns false if there is a contradiction. The
	//  caller (C_Yuzen) stops the run there.
	static bool Check(Input& inp);

private:
	CheckParameters(const CheckParameters&);
	void operator=(const CheckParameters&);
};


#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif


#endif //_DEF_CHECK_PARAMETERS_
