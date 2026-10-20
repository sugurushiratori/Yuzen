/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  SolverModel.cpp : factory that builds a model from the input settings
*
*  This file is the only one that knows the concrete models. The top-level
*  driver (Yuzen.cpp) sees nothing but ISolverModel. Adding a new model means
*  adding one line here.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "SolverModel.h"

#include "WellMixed_Model_1D.h"
#include "WellMixed_Model_2D.h"
#include "ConcZ_Model_2D.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

ISolverModel* MakeSolverModel(const std::string &concModel, const std::string &dimension)
{
	const bool x  = (dimension == "xh");
	const bool xy = (dimension == "xyh");

	//  Models that take the concentration to be uniform across the film thickness.
	if(concModel == "WM")
	{
		if(x)  return new WellMixed_Model_1D();
		if(xy) return new WellMixed_Model_2D();
		return 0;
	}

	//  Models that solve the concentration across the film thickness as well.
	//  "xh"  : the film is one-dimensional, plus a concentration profile across the thickness
	//  "xyh" : the film is two-dimensional, plus a concentration profile across the thickness
	//          (pseudo three-dimensional)
	//  The same class handles both. The dimension of the film follows from the
	//  grid's ny(); when ny()==1 the y-direction terms never appear.
	if(concModel == "CZ")
	{
		if(x || xy) return new ConcZ_Model_2D();
		return 0;
	}

	//  When no model is given, or the string is not one we know, decide on the dimension alone.
	if(x) return new WellMixed_Model_1D();
	return new WellMixed_Model_2D();
}

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif
