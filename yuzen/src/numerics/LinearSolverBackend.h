/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LinearSolverBackend.h : the seam at which the linear-solver implementation
*                          is swapped out
*
*  This is the Strategy pattern. There is one job, solving a linear system,
*  and two implementations of it, chosen by the input file through
*    Numerics.*.LinearSolver.Backend = "eigen" | "amgcl" | "auto"
*  The caller (LinearSolver) knows nothing but this interface.
*
*  The implementation classes are confined to LinearSolverBackend.cpp and can
*  only be created from outside through MakeLinearSolverBackend(). The AMGCL
*  headers are confined there as well, so whether AMGCL is present never
*  surfaces in this header.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _LINEAR_SOLVER_BACKEND_INCLUDED_
#define _LINEAR_SOLVER_BACKEND_INCLUDED_

#include <string>

#include "Vector.h"

class CRS_Matrix;

//--------------------------------------------------------------------------
//  The complete set of solver settings.
//
//  The settings are passed in on every solve. The implementation does not
//  hold on to them, because the caller sometimes rewrites MaxIteration()
//  after Build() has already run.
//--------------------------------------------------------------------------
struct LinearSolverSettings
{
	enum tPrecondition { tp_None, tp_Jacobi, tp_ILU };

	tPrecondition precondition;			// kind of preconditioner
	bool          scaling;				// whether diagonal scaling is applied
	int           maxIteration;			// upper limit on the iteration
	double        epsilon;				// convergence tolerance
	bool          conc;					// the solver of a concentration equation (opt concpc, util/Opt.h)

	LinearSolverSettings()
	 : precondition(tp_ILU), scaling(true), maxIteration(100), epsilon(1e-12), conc(false)
	{}
};

//--------------------------------------------------------------------------
//  The interface every solver implementation shares.
//--------------------------------------------------------------------------
class ILinearSolverBackend
{
public:
	virtual ~ILinearSolverBackend(){}

	/** @brief The name written in the input file ("eigen" / "amgcl"). */
	virtual const char* Name() const = 0;

	/**
	 * @brief Allocate the working storage.
	 * @param nRow     number of rows
	 * @param nonZero  number of non-zero entries (whole matrix)
	 * @param nonZL    number of non-zero entries (lower triangle)
	 * @param nonZU    number of non-zero entries (upper triangle)
	 */
	virtual void Build(int nRow, int nonZero, int nonZL, int nonZU) = 0;

	/**
	 * @brief Solve crs * x = rhs.
	 * @param s      solver settings
	 * @param crs    coefficient matrix
	 * @param rhs    right-hand side
	 * @param x      solution. It is **also used as the initial guess**.
	 * @param iter   [out] number of iterations
	 * @param resid  [out] residual
	 * @return 0 = converged, 1 = hit the iteration limit,
	 *         2 = diverged or broke down numerically
	 *
	 * If a NaN left over from a previous failure is still in x, the solve
	 * keeps breaking down, so the caller must throw x away after a failure.
	 */
	virtual int Solve(const LinearSolverSettings &s, CRS_Matrix &crs,
	                  Vector<double> &rhs, Vector<double> &x,
	                  int &iter, double &resid) = 0;
};

/**
 * @brief Create a solver backend from its name (Factory).
 * @param name    "eigen" / "amgcl" / "auto"
 * @param inPlane number of in-plane dimensions of the film: 1 for a 1-D film
 *                (Grid.Dimension = "xh"), 2 for a 2-D one ("xyh"). Only
 *                "auto" looks at it. The default of 2 keeps the previous
 *                choice for callers that do not pass it.
 * @return The backend that was created. The caller owns it and must delete it.
 *
 * "auto" picks by the dimension of the film (measured on 112 cores, up to
 * one socket of threads):
 *
 *   1-D  ->  "eigen".  AMGCL has almost nothing to coarsen on a
 *            tridiagonal-like system: about 100 linear iterations a solve
 *            where Eigen spends 12, and 5.7 times slower at Nx = 8000. It
 *            pulls ahead only at ~1e5 cells on many threads.
 *   2-D  ->  "amgcl".  Level at 128^2; from 256^2 up 1.1 to 1.5 times faster
 *            on a socket of threads, the count a run takes by default
 *            (eigen is 1.2 to 1.6 times faster on one thread).
 *
 * The numbers are in the manual, in the chapter on parallelism ("Which
 * backend to ask for"). An unrecognized
 * name gives "eigen". Asking for "amgcl" in a build made without AMGCL warns
 * and returns "eigen".
 */
ILinearSolverBackend* MakeLinearSolverBackend(const std::string &name,
                                              int inPlane = 2);

#endif //_LINEAR_SOLVER_BACKEND_INCLUDED_
