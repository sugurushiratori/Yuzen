/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_linear_solver.cpp : regression tests for the linear solver
*
*  The problem is 1-D Poisson, whose analytic solution is known, and the
*  tests check the solver against it: that it converges to that solution,
*  that the error falls at second order under mesh refinement, and that the
*  answer -- down to the last bit, the iteration count and the residual --
*  does not depend on how many threads are used.
*
*  They run through ILinearSolverBackend rather than against one
*  implementation, so they apply to whichever backend is named here.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "doctest/doctest.h"
#include "CRS_Matrix.h"
#include "LinearSolverBackend.h"

#include <memory>

#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace {

// u'' = 2 on [-1, 1], u(-1) = u(1) = 0   ->   u(x) = x^2 - 1
//
// Second-order central differences at the interior points; at both ends,
// the average with the ghost point is set to 0.
struct Poisson1D
{
	int n;
	double dx;
	CRS_Matrix mat;
	Vector<double> rhs;

	explicit Poisson1D(int n_) : n(n_), dx(2.0/double(n_-2)), rhs(n_)
	{
		mat.Build(n);
		mat.Diag(0) = 0.5;  mat(0,1) = 0.5;  rhs[0] = 0.0;
		for(int i = 1; i < n-1; ++i)
		{
			mat(i, i-1) =  1.0/(dx*dx);
			mat.Diag(i) = -2.0/(dx*dx);
			mat(i, i+1) =  1.0/(dx*dx);
			rhs[i] = 2.0;
		}
		mat.Diag(n-1) = 0.5; mat(n-1, n-2) = 0.5; rhs[n-1] = 0.0;
	}

	double x(int i) const { return -1.0 - 0.5*dx + double(i) * dx; }
	double exact(int i) const { const double xx = x(i); return xx*xx - 1.0; }
};

//  The backend under test. "eigen" is the one that is bit-identical whatever
//  the thread count, which is what the last test below pins down; "amgcl" is
//  not, and is not built into every configuration.
const char *kBackend = "eigen";

int solve(Poisson1D &p, Vector<double> &x,
          LinearSolverSettings::tPrecondition pre, int &iter, double &resid)
{
	std::unique_ptr<ILinearSolverBackend> bcg(MakeLinearSolverBackend(kBackend));
	REQUIRE(bcg);
	bcg->Build(p.mat.RowSize(), p.mat.NonZero(),
	           p.mat.NonZ_LU().first, p.mat.NonZ_LU().second);

	LinearSolverSettings s;
	s.epsilon      = 1e-13;
	s.maxIteration = 10 * p.n;
	s.precondition = pre;
	s.scaling      = true;

	x.resize(p.n);
	x = 0.0;
	return bcg->Solve(s, p.mat, p.rhs, x, iter, resid);
}

} // namespace

TEST_CASE("Linear solver: matches the analytic solution u = x^2 - 1 of 1-D Poisson")
{
	const int n = 202;			// 200 interior points + the two ends

	LinearSolverSettings::tPrecondition pres[3] = {
		LinearSolverSettings::tp_None, LinearSolverSettings::tp_Jacobi,
		LinearSolverSettings::tp_ILU };
	const char* names[3] = {"None", "Jacobi", "ILU"};

	for(int k = 0; k < 3; ++k)
	{
		INFO("precondition = " << names[k]);
		Poisson1D p(n);
		Vector<double> x;
		int iter = 0; double resid = 0;
		const int ch = solve(p, x, pres[k], iter, resid);

		CHECK(ch == 0);				// it converged
		CHECK(resid < 1e-11);

		// It must agree with the analytic solution at the interior points.
		// The discretization error is O(dx^2). Since dx = 2/200 = 0.01 we
		// take 1e-4 as the upper bound (the measured value is 2.5e-5).
		// The order itself is checked by the convergence test below.
		double err = 0.0;
		for(int i = 1; i < n-1; ++i) err = std::max(err, std::fabs(x[i] - p.exact(i)));
		CAPTURE(err);
		CHECK(err < 1e-4);
	}
}

TEST_CASE("Linear solver: the error falls at second order under mesh refinement")
{
	double prev_err = -1.0;
	double prev_dx  = 0.0;
	for(int n = 52; n <= 402; n = (n-2)*2 + 2)
	{
		Poisson1D p(n);
		Vector<double> x;
		int iter = 0; double resid = 0;
		REQUIRE(solve(p, x, LinearSolverSettings::tp_ILU, iter, resid) == 0);

		double err = 0.0;
		for(int i = 1; i < n-1; ++i) err = std::max(err, std::fabs(x[i] - p.exact(i)));

		if(prev_err > 0.0)
		{
			// Halving dx should reduce the error to about 1/4.
			// The rounding of the convergence test rides on top of that,
			// so we test loosely for 1/2 or less.
			CAPTURE(prev_dx); CAPTURE(p.dx); CAPTURE(prev_err); CAPTURE(err);
			CHECK(err < prev_err * 0.5);
		}
		prev_err = err;
		prev_dx  = p.dx;
	}
}

TEST_CASE("Linear solver: solving the same problem twice gives exactly the same answer")
{
	// The ILU construction and its forward/backward substitution have no
	//          data races: a race would make the result change with the
	//          thread schedule.
	const int n = 202;
	Poisson1D p1(n), p2(n);
	Vector<double> x1, x2;
	int i1 = 0, i2 = 0; double r1 = 0, r2 = 0;

	REQUIRE(solve(p1, x1, LinearSolverSettings::tp_ILU, i1, r1) == 0);
	REQUIRE(solve(p2, x2, LinearSolverSettings::tp_ILU, i2, r2) == 0);

	CHECK(i1 == i2);
	CHECK(r1 == r2);
	for(int i = 0; i < n; ++i) { CAPTURE(i); CHECK(x1[i] == x2[i]); }
}

#if defined(_OPENMP)
TEST_CASE("Linear solver: the answer is bit-identical whatever the thread count")
{
	const int n = 402;
	const int saved = omp_get_max_threads();

	omp_set_num_threads(1);
	Poisson1D p1(n);
	Vector<double> ref;
	int iref = 0; double rref = 0;
	REQUIRE(solve(p1, ref, LinearSolverSettings::tp_ILU, iref, rref) == 0);

	for(int th = 2; th <= 8; th *= 2)
	{
		omp_set_num_threads(th);
		CAPTURE(th);
		Poisson1D p(n);
		Vector<double> x;
		int iter = 0; double resid = 0;
		REQUIRE(solve(p, x, LinearSolverSettings::tp_ILU, iter, resid) == 0);

		CHECK(iter == iref);		// even the iteration count must match
		CHECK(resid == rref);
		for(int i = 0; i < n; ++i) { CAPTURE(i); CHECK(x[i] == ref[i]); }
	}
	omp_set_num_threads(saved);
}
#endif

TEST_CASE("Linear solver: a zero right-hand side returns the trivial solution")
{
	const int n = 52;
	Poisson1D p(n);
	p.rhs = 0.0;
	Vector<double> x;
	int iter = 0; double resid = 0;
	CHECK(solve(p, x, LinearSolverSettings::tp_ILU, iter, resid) == 0);
	//  Only that the answer is the trivial one. How many iterations it takes
	//  to get there is up to the backend: the solver this file used to test
	//  short-circuited a zero right-hand side and returned in 0 iterations,
	//  which is an optimization, not a property of the problem.
	for(int i = 0; i < n; ++i) CHECK(x[i] == doctest::Approx(0.0));
}

TEST_CASE("Linear solver: Solve restores the matrix and the right-hand side it was given")
{
	// The diagonal scaling in the preconditioner destroys the matrix and
	// the right-hand side, so Solve() restores them from a backup.
	// This pins down that contract.
	const int n = 52;
	Poisson1D p(n);

	CRS_Matrix mat_before(p.mat);
	Vector<double> rhs_before(p.rhs);

	Vector<double> x;
	int iter = 0; double resid = 0;
	REQUIRE(solve(p, x, LinearSolverSettings::tp_ILU, iter, resid) == 0);

	for(int i = 0; i < n; ++i)
	{
		CAPTURE(i);
		CHECK(p.rhs[i] == rhs_before[i]);
		CHECK(p.mat.Diag(i) == mat_before.Diag(i));
		CHECK(p.mat.NonZero(i) == mat_before.NonZero(i));
	}
}
