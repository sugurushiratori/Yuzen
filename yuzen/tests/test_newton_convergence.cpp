/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_newton_convergence.cpp : pins down how the outer Newton-Kantorovich
*                                iteration exits
*
*  The thickness equation and the concentration equation are coupled by
*  block Gauss-Seidel, and only once that has converged do we have "the
*  solution of the two solved fully implicitly together". A convergence
*  test that looked only at the correction to the thickness would let the
*  loop exit while the concentration is still moving, and a fixed minimum
*  number of iterations would only paper over that.
*
*  Here we do not solve the real equations. We feed the per-iteration
*  corrections from a script and examine only how the loop exits. What we
*  check is:
*    - it does not exit while the correction on the coupled side is large
*    - it exits immediately once both are small
*    - NewtonMinIterations can force a minimum iteration count
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "doctest/doctest.h"

#include "NewtonKantorovich.h"
#include "TimeSteps.h"

#include <vector>

using namespace yuzen;

namespace {

//------------------------------------------------------------------
//  A Newton solver that does nothing but return the per-iteration
//  corrections from a script. It performs neither a linear solve nor any
//  matrix assembly.
//------------------------------------------------------------------
class ScriptedNK : public C_NewtonKantorovich
{
public:
	ScriptedNK(const std::vector<double> &dH, const std::vector<double> &dC)
		: _dH(dH), _dC(dC), _call(0) {}

	int Calls() const { return _call; }

private:
	std::vector<double> _dH;	// thickness correction returned at each iteration
	std::vector<double> _dC;	// coupled-side correction returned at each iteration
	int _call;

	double Pick(const std::vector<double> &v) const
	{
		if(v.empty()) return 0.0;
		const size_t k = (size_t)_call;
		return (k < v.size()) ? v[k] : v.back();
	}

	virtual void SetFlagBC() {}
	virtual void SetMatrixA() {}
	virtual void SetMatrixJacobi() {}

	//  Pretend the coupled side (the concentration equation in the real
	//  solver) has been solved, and report the correction
	virtual bool EveryNKLoop()
	{
		_coupledCorrection = Pick(_dC);
		return true;
	}

	//  Instead of a linear solve, put in the scripted correction as it is
	virtual int Solve(int &iter, double &resid)
	{
		const double d = Pick(_dH);
		for(int i = 0; i < _correctVec.size(); ++i) _correctVec[i] = d;
		iter = 1;
		resid = 0.0;
		++_call;
		return 0;
	}
};


//  Returns the iteration count
int RunOneStep(ScriptedNK &nk, int minIter)
{
	C_TimeSteps ts;
	ts.Build(0.1, 0.0, 1.0, false);

	nk.Build(8, ts, C_ImplicitSolver::CrankNicholson, false);
	nk.SetEpsilonNK(1e-10);
	nk.SetMaxIterNK(30);
	nk.SetMinIterNK(minIter);
	nk.Initialize();

	int iter = 0;
	double resid = 0;
	nk.TimeAdvance(iter, resid);
	return iter;
}

const double TINY = 1e-20;		// a correction counted as converged
const double BIG  = 1.0;		// a correction that is still moving

} // namespace


TEST_CASE("Newton: it does not exit while the coupled side is still moving")
{
	//  The thickness is converged throughout. The coupled side only
	//  settles at the third iteration.
	std::vector<double> dH(1, TINY);
	std::vector<double> dC;
	dC.push_back(BIG);
	dC.push_back(BIG);
	dC.push_back(TINY);

	ScriptedNK nk(dH, dC);
	CHECK(RunOneStep(nk, 1) == 3);
}


TEST_CASE("Newton: it exits in one iteration when both corrections are small")
{
	std::vector<double> dH(1, TINY);
	std::vector<double> dC(1, TINY);

	ScriptedNK nk(dH, dC);
	CHECK(RunOneStep(nk, 1) == 1);
}


TEST_CASE("Newton: with no coupling the thickness alone decides the exit")
{
	//  In a model that does not report a coupled-side correction (ConcZ,
	//  for example), _coupledCorrection stays 0 and the test is decided by
	//  the thickness alone.
	std::vector<double> dH;
	dH.push_back(BIG);
	dH.push_back(BIG);
	dH.push_back(TINY);
	std::vector<double> dC;			// empty = always 0

	ScriptedNK nk(dH, dC);
	CHECK(RunOneStep(nk, 1) == 3);
}


TEST_CASE("Newton: NewtonMinIterations sets the minimum iteration count")
{
	//  Even when it converges immediately, the loop runs up to the
	//  requested count, which comes from the input (NewtonMinIterations), not from a
	//  number fixed in the expression.
	std::vector<double> dH(1, TINY);
	std::vector<double> dC(1, TINY);

	{ ScriptedNK nk(dH, dC); CHECK(RunOneStep(nk, 6) == 6); }
	{ ScriptedNK nk(dH, dC); CHECK(RunOneStep(nk, 3) == 3); }

	//  Values of 0 or below are rounded up to 1
	{ ScriptedNK nk(dH, dC); CHECK(RunOneStep(nk, 0) == 1); }
}


//==================================================================
//  Behavior when the iteration hits the floor set by the accuracy of the
//  linear solver
//
//  Newton converges quadratically, so the correction drops sharply, but
//  once it reaches the accuracy of the linear solve (about 1e-10 at this
//  problem size) it stops decreasing and merely bounces around at random.
//
//  The original test was only "converged once it falls below the absolute
//  value NewtonEpsilon", and the real inputs used NewtonEpsilon = 1e-12. Since the
//  floor is 1e-10, that value is never reached: the loop runs all the way
//  to NewtonMaxIterations and reports "failed to converge". The caller sees that and
//  shrinks the time step. In other words, the solver spun its wheels even
//  though it had converged well enough, and destroyed the time step on top
//  of that. This alone made computations with a contact line 100 times
//  slower.
//==================================================================
namespace {

//  Quadratic convergence -> scatter at the 1e-10 floor -> converged (if
//  retried)
std::vector<double> FloorScript()
{
	std::vector<double> v;
	v.push_back(2.108252e-03);
	v.push_back(3.576681e-05);
	v.push_back(7.291922e-07);
	v.push_back(1.669462e-08);
	v.push_back(3.868198e-10);
	//  The floor starts here. It no longer decreases, it just bounces
	v.push_back(9.020539e-11);
	v.push_back(8.684096e-11);
	v.push_back(1.567753e-10);
	v.push_back(9.732429e-11);
	v.push_back(3.120448e-11);
	v.push_back(4.350973e-11);
	v.push_back(1.190828e-10);
	v.push_back(2.507970e-11);
	v.push_back(1.002455e-10);
	v.push_back(4.546102e-11);
	v.push_back(7.090800e-11);
	v.push_back(7.484886e-11);
	v.push_back(9.904685e-11);
	v.push_back(9.818438e-11);
	v.push_back(1.619510e-10);
	return v;
}

//  Advance one step with time-step adjustment enabled, and return the
//  iteration count and the final dt
void RunAdj(ScriptedNK &nk, double epsAbs, double epsRel, int maxIter,
            int &iter, double &dt_end)
{
	C_TimeSteps ts;
	ts.Build(0.1, 0.0, 1.0, true);
	ts.SetMinDT(1e-8);
	ts.SetMaxDT(1.0);

	nk.Build(8, ts, C_ImplicitSolver::CrankNicholson, true);
	nk.SetEpsilonNK(epsAbs);
	nk.SetEpsilonRelNK(epsRel);
	nk.SetMaxIterNK(maxIter);
	nk.SetMinIterNK(1);
	nk.SetAdjustTimeChangeFactor(0.5);
	nk.Initialize();

	double resid = 0;
	iter = 0;
	nk.TimeAdvance(iter, resid);
	dt_end = ts.DT().now;
}

} // namespace


TEST_CASE("Newton: the relative test exits without chasing an unreachable absolute")
{
	//  NewtonEpsilon = 1e-12 is below the floor (1e-10), so the absolute test
	//  can never be met. The relative test (1e-6 times the first
	//  correction) converges at the fifth iteration.
	//    2.108e-03 * 1e-6 = 2.11e-09 > 3.868e-10  -> converged at the
	//    fifth iteration
	std::vector<double> dC;
	ScriptedNK nk(FloorScript(), dC);

	int iter = 0; double dt = 0;
	RunAdj(nk, 1e-12, 1e-6, 50, iter, dt);

	CHECK(iter == 5);
	//  It converged, so the time step is not shrunk
	CHECK(dt == doctest::Approx(0.1));
}


TEST_CASE("Newton: with the relative test off, the correction floor exits by stagnation")
{
	//  NewtonEpsilonRel = 0 disables the relative test. Even so, the loop is
	//  cut off once "the smallest value so far is not improved upon"
	//  persists. It does not run all the way to NewtonMaxIterations.
	std::vector<double> dC;
	ScriptedNK nk(FloorScript(), dC);

	int iter = 0; double dt = 0;
	RunAdj(nk, 1e-12, 0.0, 50, iter, dt);

	CHECK(iter < 25);			// it does not run to 50
	CHECK(iter > 5);			// it does look a few times after reaching the floor
	//  It has dropped by more than 4 orders of magnitude from the first
	//  correction, so this counts as "solved". The step is not shrunk.
	CHECK(dt == doctest::Approx(0.1));
}


//==================================================================
//  A stagnation far from convergence fails the step.
//
//  Kept, such a step is not a solution of the step (a thick nozzle film
//  kept the residual as a volume error of 2.2e-8), so under step control it
//  is cut and solved again, counted, like a linear solve that did not
//  converge. (At a fixed step it stops the run, FailedFixedStep, which ends
//  the process: tools/features.py checks that.)
//==================================================================
namespace {

//  n iterations of a correction that never drops, then one that converges
std::vector<double> StallThenConverge(int trials, int perTrial, double v)
{
	std::vector<double> s;
	for(int t = 0; t < trials; ++t)
		for(int k = 0; k < perTrial; ++k) s.push_back(v);
	s.push_back(1e-3);
	s.push_back(TINY);
	return s;
}

} // namespace


TEST_CASE("Newton: a constant correction stops on stagnation, within NewtonMaxIterations, and the step is cut")
{
	//  The correction is constant, so the stagnation test (no improvement
	//  for 6 iterations in a row) stops the first trial at its 7th
	//  iteration, NewtonMaxIterations = 7 not reached. It has not dropped at
	//  all from the first one: the step fails, is cut to half and converges.
	std::vector<double> dC;
	ScriptedNK nk(StallThenConverge(1, 7, BIG), dC);

	int iter = 0; double dt = 0;
	RunAdj(nk, 1e-10, 1e-6, 7, iter, dt);

	CHECK(nk.Calls() == 9);			// 7 in the failed trial, 2 in the second
	CHECK(iter == 2);
	CHECK(dt == doctest::Approx(0.05));
	CHECK(nk.StagnationCuts() == 1);
}


TEST_CASE("Newton: a correction that never drops is not reported as converged")
{
	//  Two trials stagnate at the first correction (as on the contact-line
	//  cases, where the correction stopped at 0.13 to 150 times the first):
	//  both are cut, and the third, at a quarter of the step, converges.
	std::vector<double> dC;
	ScriptedNK nk(StallThenConverge(2, 7, 1.0), dC);

	int iter = 0; double dt = 0;
	RunAdj(nk, 1e-12, 1e-6, 50, iter, dt);

	CHECK(nk.Calls() == 16);
	CHECK(dt == doctest::Approx(0.025));
	CHECK(nk.StagnationCuts() == 2);
}


TEST_CASE("Newton: a stagnation within 1000 eps |phi| is accepted")
{
	//  A film of 1e5 (eps |phi| = 2.2e-11) whose first correction, 1e-9, is
	//  already close to its floor: the correction stops at 1e-10, neither
	//  below NewtonEpsilon (1e-12) nor four decades below the first, but
	//  within 1000 eps |phi|, which the film cannot be updated more finely
	//  than by much. That is the floor, not a failure: the step is accepted
	//  and not cut.
	std::vector<double> dH;
	dH.push_back(1e-9);
	for(int k = 0; k < 10; ++k) dH.push_back(1e-10);
	std::vector<double> dC;
	ScriptedNK nk(dH, dC);

	C_TimeSteps ts;
	ts.Build(0.1, 0.0, 1.0, true);
	ts.SetMinDT(1e-8);
	ts.SetMaxDT(1.0);
	nk.Build(8, ts, C_ImplicitSolver::CrankNicholson, true);
	Vector<double> phi(8);
	phi = 1e5;
	nk.SetX(phi);
	nk.SetEpsilonNK(1e-12);
	nk.SetEpsilonRelNK(1e-6);
	nk.SetMaxIterNK(50);
	nk.SetMinIterNK(1);
	nk.SetAdjustTimeChangeFactor(0.5);
	nk.Initialize();

	int iter = 0; double resid = 0;
	nk.TimeAdvance(iter, resid);
	CHECK(iter == 8);			// the drop to 1e-10, then six that do not improve on it
	CHECK(ts.DT().now == doctest::Approx(0.1));
	CHECK(nk.StagnationCuts() == 0);

	//  The same stagnation on a film of 1e-8 (the corrections added to a
	//  film that starts at 0) is far from its floor, and the step is cut.
	std::vector<double> dH2(dH);
	dH2.push_back(1e-3);
	dH2.push_back(TINY);
	ScriptedNK nk2(dH2, dC);
	int iter2 = 0; double dt2 = 0;
	RunAdj(nk2, 1e-12, 1e-6, 50, iter2, dt2);
	CHECK(dt2 == doctest::Approx(0.05));
	CHECK(nk2.StagnationCuts() == 1);
}
