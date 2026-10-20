/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  SolverModel.h : the seam at which the model being solved is swapped
*
*  There is a single job here - assemble the height equation and the
*  concentration equation and march them in time - but three implementations
*  of it, differing in the dimension and in how the concentration is treated.
*  Which one is used is decided by the input:
*    Physics.ConcField.Model  = "WM" | "CZ"
*    Grid.Dimension           = "xh" | "xyh"
*  The top-level driver (C_Yuzen) knows nothing but this interface.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _SOLVER_MODEL_INCLUDED_
#define _SOLVER_MODEL_INCLUDED_

#include <string>

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class Input;

//--------------------------------------------------------------------------
//  Common interface of the models that can be solved.
//--------------------------------------------------------------------------
class ISolverModel
{
public:
	virtual ~ISolverModel(){}

	/** @brief Name of the model, as written to the log. */
	virtual const char* Name() const = 0;

	/**
	 * @brief Assemble the grid, fields, equations, boundary conditions and output from the input.
	 * @param inp  input that has already been read
	 */
	virtual void Build(Input &inp) = 0;

	/** @brief Advance in time until a termination condition is met. */
	virtual void Run() = 0;

	/**
	 * @brief Whether the run ended because the solution diverged.
	 *
	 * Checked by the top-level driver after Run() so that a diverged
	 * calculation leaves a non-zero exit status instead of looking like a
	 * normal finish. The default is false, so a model that does not detect
	 * divergence keeps its previous behavior.
	 */
	virtual bool Diverged() const { return false; }
};

/**
 * @brief Create the solver model selected by the input (factory).
 * @param concModel  Physics.ConcField.Model ("WM" / "CZ")
 * @param dimension  Grid.Dimension ("xh" / "xyh")
 * @return The newly created model; the caller owns it and must delete it.
 *         Returns 0 if the combination is not a valid one.
 */
ISolverModel* MakeSolverModel(const std::string &concModel, const std::string &dimension);

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_SOLVER_MODEL_INCLUDED_
