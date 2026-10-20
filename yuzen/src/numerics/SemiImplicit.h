/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  SemiImplicit.h : semi-implicit time integration (the nonlinear coefficients are frozen at the previous step)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _SEMI_IMPLICIT_INCLUDED_
#define _SEMI_IMPLICIT_INCLUDED_

#include "ImplicitSolver.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif


//--------------------------------------------------------------------------
//  Semi-implicit time integration.
//
//  The nonlinear coefficients (the mobility and the material properties) are
//  frozen at the values of the previous step, and the next values are
//  obtained by solving the system once. There is no Newton iteration, so a
//  step is cheap; but where the nonlinearity is strong, accuracy is only
//  obtained by making the time step small.
//
//  The concentration equation is solved this way. On the film-height
//  equation it is used by the ConcZ model alone: the height in the
//  well-mixed model uses C_NewtonKantorovich.
//--------------------------------------------------------------------------
class C_SemiImplicitEquation : public C_ImplicitSolver
{
public:
	/** @brief Construct with the defaults inherited from C_ImplicitSolver. */
	C_SemiImplicitEquation() : C_ImplicitSolver() {}
	virtual ~C_SemiImplicitEquation(){}

	/** @brief Set up the working storage and the matrix structure. */
	virtual void Initialize();

	/**
	 * @brief Advance one time step, doing the trials and the step adjustment
	 *        together.
	 * @param iter   [out] number of iterations used
	 * @param resid  [out] residual reached
	 * @return 0 on success; non-zero when the step could not be advanced.
	 */
	virtual int TimeAdvance(int &iter, double &resid);

	/**
	 * @brief Solve once only, adjusting nothing and committing nothing.
	 * @param iter   [out] number of iterations used
	 * @param resid  [out] residual reached
	 * @return true when the solve succeeded: false when the linear solver
	 *         did not converge (without step control: at the one step
	 *         there is) or AfterTimeStepTrial() failed.
	 *
	 * This is the entry point called on every iteration from the outer
	 * Newton loop of the film-height equation. Neither the time step nor the
	 * advance to the next step is touched here; use CommitStep() for that.
	 */
	bool TimeStepTrial(int &iter, double &resid);

	/**
	 * @brief Hooks called before and after each trial, for a derived class
	 *        to refresh its fields.
	 * @return false to abort the trial. AfterTimeStepTrial() runs after a
	 *         solve the step test accepts (every solve, without step
	 *         control); false from it fails the trial: the step is cut and
	 *         solved again under step control, and without it
	 *         TimeStepTrial() returns false.
	 */
	virtual bool EveryTimeStepTrial(){return true;}
	virtual bool AfterTimeStepTrial(){return true;}

protected:
	virtual double CheckTimeStepCriterion() {return 0.1;}

	virtual int Solve(int &iter, double &resid);

	//  The pieces of Solve(), for a derived class that stages its own step.
	virtual void MakeMatrix_IA();
	virtual void MakeRHS();
	bool SolveExplicit();

private:

	virtual void SetFlagBC() = 0;
	virtual void SetMatrixA() = 0;

	C_SemiImplicitEquation(const C_SemiImplicitEquation& c);
	void operator=(const C_SemiImplicitEquation& c);

};


#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_SEMI_IMPLICIT_INCLUDED_
