/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  NewtonKantorovich.h : nonlinear time integration by Newton-Kantorovich iteration
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _DEF_NEWTON_KANTOROVICH_INCLUDED_
#define _DEF_NEWTON_KANTOROVICH_INCLUDED_

#include <cstdlib>

#include "ImplicitSolver.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif


//--------------------------------------------------------------------------
//  Nonlinear time integration by Newton-Kantorovich iteration.
//
//  In the film-height equation the mobility h^3/mu(c) depends on the
//  unknowns, so a Newton iteration inside each step resolves the
//  nonlinearity as well. One iteration consists of
//
//    1. EveryNKLoop()      update the fields from the current estimate and
//                          solve the coupled side (the concentration)
//    2. SetMatrixA()       assemble the matrix A
//    3. SetMatrixJacobi()  assemble the Jacobian (once every
//                          _jacobiInterval iterations)
//    4. the linear solve   obtain the correction _correctVec
//    5. add it to _phi_guess
//
//  The convergence test is taken on the **larger** of the film-height
//  correction and the correction of the coupled side. That is what makes the
//  iteration converge to the solution of the height equation and the
//  concentration equation solved fully implicitly as one system
//  (docs/PHASE6.md). Both an absolute test (_EpsilonNK) and a test relative
//  to the first correction of the step (_EpsilonRelNK) are applied. Why the
//  relative test is needed is explained in the comments in the .cpp file.
//--------------------------------------------------------------------------
class C_NewtonKantorovich : public C_ImplicitSolver
{
public:
	/**
	 * @brief Construct with the default iteration limits and tolerances.
	 *
	 * The environment variable YUZEN_DEBUG_RESIDUAL is read here; when it is
	 * set, the breakdown of the convergence test is written to standard
	 * error during the run.
	 */
	C_NewtonKantorovich() : C_ImplicitSolver() {
		_bcg_resid = 0;
		_bcg_iter = 0;
		_bcg_iter_step = 0;
		_MaxIterNK = 30;
		_EpsilonNK = 1.0e-10;
		_EpsilonRelNK = 1.0e-6;
		_jacobiInterval = 1;
		_nkIter = 0;
		_MinIterNK = 1;
		_coupledCorrection = 0;
		_stagnationCuts = 0;
		_warnedStagnation = false;
		_dbgResid = (std::getenv("YUZEN_DEBUG_RESIDUAL") != 0);
	}
	virtual ~C_NewtonKantorovich(){}

	/** @brief Set up the working storage for the Newton iteration. */
	virtual void Initialize();
	/**
	 * @brief Advance one time step, iterating Newton-Kantorovich until the
	 *        correction meets the convergence test.
	 * @param iter      [out] number of iterations used
	 * @param resid_eq  [out] residual reached
	 * @return 0 on success; non-zero when the step could not be advanced.
	 */
	virtual int TimeAdvance(int &iter, double &resid_eq);

	/**
	 * @brief Upper limit on the outer (Newton) iteration.
	 *
	 * Reaching it is reported as "did not converge", so with an adaptive
	 * time step the step then shrinks.
	 */
	void SetMaxIterNK(int nn) {_MaxIterNK = nn;}

	/**
	 * @brief Minimum number of outer iterations.
	 * @param nn  the minimum, clamped up to at least 1.
	 *
	 * With the default of 1 the loop leaves on the convergence test alone.
	 * Because the correction of the coupled side is part of that test, 1 is
	 * normally the right value; this is kept for experiments that need a
	 * fixed iteration count.
	 */
	void SetMinIterNK(int nn) {_MinIterNK = (nn < 1 ? 1 : nn);}
	/**
	 * @brief Absolute convergence tolerance: converged once the correction
	 *        falls below this.
	 *
	 * The correction has a floor coming from cancellation, and that floor is
	 * proportional to the time step. Testing on the absolute value alone
	 * therefore couples convergence to the step size ("shrink the step and
	 * it will get there"), so this is normally used together with the
	 * relative test below (docs/PHASE7.md 8.5).
	 */
	void SetEpsilonNK(double eps){_EpsilonNK = eps;}
	/**
	 * @brief Relative convergence tolerance, measured against the **first
	 *        correction** of the step.
	 * @param eps  the ratio; zero or negative disables the relative test and
	 *             leaves only the absolute one.
	 *
	 * This test is independent of the scale of the solution, of the time
	 * step and of the grid.
	 */
	void SetEpsilonRelNK(double eps){_EpsilonRelNK = eps;}

	/**
	 * @brief Modified Newton: rebuild the Jacobian once every n iterations
	 *        instead of every iteration.
	 * @param n  the interval, clamped up to at least 1. The default of 1
	 *           rebuilds every iteration, which is the unapproximated method.
	 *
	 * Measured on the 1D leveling_1d_reg case, one Newton iteration splits as
	 * SetMatrixA 19% / SetMatrixJacobi 43% / IA+RHS 2% / linear solve 37%,
	 * so assembling the Jacobian is the largest single item. The Jacobian
	 * only sets how fast the iteration converges, not what it converges to,
	 * so rebuilding it less often leaves the solution unchanged, though the
	 * iteration count can rise.
	 */
	void SetJacobiInterval(int n) { _jacobiInterval = (n < 1 ? 1 : n); }

	//  Iteration count and residual of the most recent linear solve, as
	//  reported in CalcInfo.dat.
	int BCG_Iter() const {return _bcg_iter;}
	double BCG_Resid() const {return _bcg_resid;}

	//  A counter kept for performance measurements. It takes no part in the
	//  numerics at all.
	//  Where BCG_Iter() is the count in the last Newton iteration, this is
	//  the total of the linear iterations over every Newton iteration of one
	//  time step. It is used in scaling measurements, to confirm that the
	//  amount of work is the same.
	int BCG_IterOfStep() const {return _bcg_iter_step;}

	//  Steps that failed because the correction stagnated far from
	//  convergence (each one cut and solved again, or the run stopped at a
	//  fixed step). JobInfo.txt reports it.
	long long StagnationCuts() const {return _stagnationCuts;}

	//  The last Newton correction. The linear solver starts the next one from
	//  it (solveWithGuess), so a restart file carries it.
	const Vector<double>& LastCorrection() const {return _correctVec;}
	void SetLastCorrection(const Vector<double> &v) {_correctVec = v;}

protected:
	CRS_Matrix _crsJacobi;
	Vector<double> _correctVec;
	Vector<double> _phi_guess;
	Vector<double> _Ax_guess;

	double _bcg_resid;
	int _bcg_iter;
	int _bcg_iter_step;		// total linear iterations within one time step (for measurement)
	double _EpsilonNK;
	double _EpsilonRelNK;	// tolerance relative to the first correction
	int _MaxIterNK;
	int _MinIterNK;
	int _jacobiInterval;	// interval at which the Jacobian is rebuilt (1 = every iteration)
	int _nkIter;			// number of the Newton iteration within the time step
	//  opt ew (util/Opt.h): the right-hand side of the last pass and its linear
	//  tolerance (Eisenstat-Walker); opt extrap (util/Opt.h): how far the start of
	//  the iteration was moved from the last step
	double _ewRhsPrev = 0, _ewEta = 0;
	double _extrapMove = 0;
	bool _extrapBlocked = false;	// this attempt starts from phi^n (a retry)

	//------------------------------------------------------------------
	//  The correction of the other field this equation is coupled to.
	//  It is updated inside EveryNKLoop, before the film height is solved,
	//  and it enters the convergence test. In a model with no coupling it
	//  stays zero, so the test is then decided by the film height alone.
	//------------------------------------------------------------------
	//  The correction of the coupled side (the concentration). EveryNKLoop
	//  writes it on every iteration, and the convergence test is taken on
	//  the larger of this and the film-height correction.
	double _coupledCorrection;

	//  Steps failed on a stagnation far from convergence, and whether the
	//  first of them under step control has been reported
	long long _stagnationCuts;
	bool _warnedStagnation;

	//  Write the breakdown of the convergence test to standard error
	//  (environment variable YUZEN_DEBUG_RESIDUAL)
	bool _dbgResid;

	virtual double CheckTimeStepCriterion() {return 0.1;}

	virtual int Solve(int &iter, double &resid);

private:
	virtual void MakeMatrix_IA();
	virtual void MakeRHS();

	bool SolveExplicit();

	virtual void SetFlagBC() = 0;
	bool NKLoop(int &iter, double &resid);
	virtual void SetMatrixA() = 0;
	virtual void SetMatrixJacobi() = 0;
	virtual void AllocAdd();

	//--------------------------------------------------------------
	//  Hooks specific to Newton-Kantorovich.
	//  This is where the film-height equation (HeightEquationBase) couples
	//  to the concentration.
	//--------------------------------------------------------------
	//  Called before the iteration begins, on every attempt at the step
	//  (retries included), and before A is assembled for the explicit half:
	//  a derived class that assembles A from a working field must put that
	//  field back to the start of the step here, and the coupled side takes
	//  its reference values.
	virtual void BeforeNKLoop(){}

	//  Called at the top of every iteration. It updates the fields from the
	//  current estimate, solves the coupled side (the concentration equation)
	//  again, and puts that side's correction into _coupledCorrection.
	//  Returning false abandons the iteration at that point.
	virtual bool EveryNKLoop(){return true;}

	//  Called once after the time step has been committed, so that the
	//  coupled side can commit its own step.
	virtual bool AfterNKLoop(){return true;}

private:
	C_NewtonKantorovich(const C_NewtonKantorovich& c);
	void operator=(const C_NewtonKantorovich& c);

};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_DEF_NEWTON_KANTOROVICH_INCLUDED_
