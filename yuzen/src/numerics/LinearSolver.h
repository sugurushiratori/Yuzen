/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LinearSolver.h : the front door to the linear solver
*
*  It holds the settings and delegates the actual solving to an
*  ILinearSolverBackend (the Strategy pattern; see LinearSolverBackend.h).
*  Which implementation is used is decided by the input file through
*    Numerics.*.LinearSolver.Backend = "eigen" | "amgcl" | "auto"
*  and the default in the input file is "auto".
*
*  Changing the backend changes how the rounding enters, so results do not
*  agree bit for bit across backends; they agree only to within the
*  convergence tolerance. **When comparing against a reference run, use the
*  same backend for both.**
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _LINEAR_SOLVER_INCLUDED_
#define _LINEAR_SOLVER_INCLUDED_

#include <string>

#include "MyMath.h"
#include "Vector.h"
#include "LinearSolverBackend.h"

class CRS_Matrix;

//--------------------------------------------------------------------------
//  The front door to the linear solver.
//
//  All it holds is the settings (preconditioner, convergence test, iteration
//  limit and so on) and a reference to the implementation that was chosen.
//  The settings are passed to the implementation on every solve, so changing
//  one after Build() still takes effect.
//--------------------------------------------------------------------------
class LinearSolver
{
public:
	//  The enumerations the strings in the input file map to. They are
	//  defined in LinearSolverBackend.h.
	typedef LinearSolverSettings::tPrecondition tPrecondition;

public:
	/**
	 * @brief Construct with the default settings and the "eigen" backend.
	 *
	 * The backend is replaced later by SetBackend(), which Configure() calls
	 * with the name taken from the input file.
	 */
	LinearSolver();
	/** @brief Destroy the solver, deleting the backend it owns. */
	~LinearSolver();

	/**
	 * @brief Allocate the working storage.
	 * @param nRow     number of rows
	 * @param nonZero  number of non-zero entries (whole matrix)
	 * @param nonZL    number of non-zero entries (lower triangle)
	 * @param nonZU    number of non-zero entries (upper triangle)
	 */
	void Build(int nRow, int nonZero, int nonZL, int nonZU);

	/**
	 * @brief Solve crs * vec_x = rhs.
	 * @param crs    coefficient matrix
	 * @param rhs    right-hand side
	 * @param vec_x  solution. It is **also used as the initial guess**: the
	 *               iteration starts from whatever is already in it.
	 * @param iter   [out] number of iterations
	 * @param resid  [out] residual
	 * @return 0 = converged, 1 = hit the iteration limit,
	 *         2 = diverged or broke down numerically
	 *
	 * If a NaN left over from a previous failure is still in vec_x, the
	 * solve keeps breaking down. The caller must throw vec_x away after a
	 * failure, as C_NewtonKantorovich::Solve does.
	 */
	int Solve(CRS_Matrix &crs, Vector<double> &rhs, Vector<double> &vec_x, int &iter, double &resid);

	/**
	 * @brief Apply every setting from the LinearSolver block of the input
	 *        file at once.
	 * @param p  the LinearSolver block of the input file
	 * @param inPlane  number of in-plane dimensions of the film (1 or 2); the
	 *                 "auto" backend resolves differently for each
	 *
	 * I_LinearSolver is generated as a separate type on the HeightEquation
	 * side and on the ConcEquation side, so this takes it as a template
	 * parameter.
	 */
	template<class T_LinearSolverParam>
	void Configure(const T_LinearSolverParam &p, int inPlane = 2)
	{
		SetBackend(p.Backend, inPlane);
		SetPrecondition(p.PreconditionType);
		_set.scaling            = p.RowScaling;
		_set.epsilon            = p.Epsilon;
		_set.maxIteration       = p.MaxIterations;
	}

	//------------------------------------------------------------------
	//  Set the options from the strings in the input file. A string that is
	//  not recognized falls back to the default.
	//------------------------------------------------------------------
	/**
	 * @brief Choose the backend by name.
	 * @param s        "eigen" / "amgcl" / "auto"
	 * @param inPlane  in-plane dimensions of the film (1 or 2); only "auto"
	 *                 uses it, to pick the backend that is faster for that
	 *                 shape. See MakeLinearSolverBackend().
	 */
	void SetBackend(const std::string &s, int inPlane = 2);
	/** @brief Mark this as the solver of a concentration equation (opt concpc, util/Opt.h). */
	void SetConcRole(bool c) { _set.conc = c; }
	/** @brief Choose the preconditioner: "None" / "Jacobi" / "ILU". */
	void SetPrecondition(const std::string &s);

	/** @brief Name of the backend currently selected. */
	std::string BackendName() const;

	//------------------------------------------------------------------
	//  Individual settings, named and used exactly as in MyBiCGStab.
	//
	//  Each has a const and a non-const form; the non-const form returns a
	//  reference, so it can be assigned to directly (one can write
	//  solver.Epsilon() = 1e-12). Normally everything is set at once from the
	//  input file through Configure() above, so these are touched directly
	//  only in tests and experiments.
	//------------------------------------------------------------------
	/**
	 * @brief Individual settings, named and used exactly as in MyBiCGStab.
	 *
	 * Each has a const and a non-const form; the non-const form returns a
	 * reference, so it can be assigned to directly
	 * (solver.Epsilon() = 1e-12). The settings are passed to the backend on
	 * every solve, so changing one after Build() still takes effect.
	 */
	const tPrecondition& PreconditionType() const {return _set.precondition;}
	      tPrecondition& PreconditionType()       {return _set.precondition;}
	const bool& Scaling() const {return _set.scaling;}
	      bool& Scaling()       {return _set.scaling;}
	const double& Epsilon() const {return _set.epsilon;}
	      double& Epsilon()       {return _set.epsilon;}
	const int& MaxIteration() const {return _set.maxIteration;}
	      int& MaxIteration()       {return _set.maxIteration;}

private:
	LinearSolverSettings  _set;			// the complete set of settings
	ILinearSolverBackend *_backend;		// the thing that does the solving (defaults to "eigen")

	//  The arguments of Build() are kept, so that the working storage can be
	//  allocated again if the implementation is reselected after Build() has
	//  run. The model side calls Build() first and Configure() afterwards.
	bool _constructed;
	int  _nRow, _nonZero, _nonZL, _nonZU;

private:
	LinearSolver(const LinearSolver&);
	LinearSolver& operator=(const LinearSolver&);
};

#endif //_LINEAR_SOLVER_INCLUDED_
