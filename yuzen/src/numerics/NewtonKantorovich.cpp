/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  NewtonKantorovich.cpp : implementation of the Newton-Kantorovich iteration
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "Opt.h"
#include <iostream>
#include "Profile.h"
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <limits>

#include "YuzenUtil.h"
using namespace yuzen::util;
using namespace Global;

#include "NewtonKantorovich.h"
#include "TimeSteps.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif



void C_NewtonKantorovich::AllocAdd()
{
	// Allocate the additional storage
//	_crsJacobi.Build(_nRows);
	_phi_guess.resize(_nRows);
	_Ax_guess.resize(_nRows);
	_correctVec.resize(_nRows);
}

void C_NewtonKantorovich::Initialize()
{
	assert(_constructed);

	_crsJacobi = _crsA;
	_crsIA = _crsA;
	ApplyA(_Ax_old, _phi_old);

	_phi_old2 = _phi_old;
	_phi_guess = _phi_old;
	_phi_new =  _phi_guess;
	_Ax_guess = _Ax_old;

	_initialized = true;
}

bool C_NewtonKantorovich::SolveExplicit()
{
	// for Explicit
	double delta_t = _timeStep->DT().now;

	for(int ix = 0; ix < _Ax_old.size(); ix++)
	{
		if(_crsA.IsBC(ix))
			_phi_guess[ix] = _phi_old[ix] -  _Ax_old[ix];
		else
			_phi_guess[ix] = _phi_old[ix] - delta_t * _Ax_old[ix];
	}
	return true;
}

int C_NewtonKantorovich::Solve(int &iter, double &resid)
{
	// Solve the system [ _crsIA * _correctVec = _Rhs ] once
	assert(_initialized);

	{ YUZEN_PROF("h.assemble"); SetMatrixA(); }	// assemble the matrix

	// The modified Newton method. With _jacobiInterval == 1 (the default) the
	// Jacobian is rebuilt every time. When it is not rebuilt, _crsIA from the
	// previous iteration can be reused as well, because MakeMatrix_IA depends
	// on nothing but _crsJacobi and dt.
	if(_nkIter % _jacobiInterval == 0)
	{
		YUZEN_PROF("h.jacobian");
		SetMatrixJacobi();		// assemble the matrix of the correction
		MakeMatrix_IA();			// assemble the matrix for the implicit time advance
	}
	++_nkIter;

	{
		YUZEN_PROF("h.rhs");
		ApplyA(_Ax_guess, _phi_guess);	// evaluate the explicit part
		MakeRHS();					// evaluate the right-hand side
	}

	//--------------------------------------------------------------
	//  opt ew (util/Opt.h): the linear tolerance of this pass from the progress of
	//  the Newton iteration (Eisenstat and Walker 1996, choice 2):
	//      eta_k = 0.9 (|F_k| / |F_k-1|)^2,
	//  safeguarded by 0.9 eta_k-1^2 while that is above 0.1, between the
	//  tolerance of the input (Epsilon) and kEtaMax; the first pass of a
	//  step at kEtaMax. F is the right-hand side of the pass, the residual
	//  of the step at the current guess.
	//--------------------------------------------------------------
	const double epsInput = _bcg.Epsilon();
	if(yuzen::opt::On(yuzen::opt::EW))
	{
		const double kEtaMax = 1.0e-4;
		const double rn = _Rhs.AbsMax();
		double eta = kEtaMax;
		if(_nkIter > 1 && _ewRhsPrev > 0)
		{
			const double ratio = rn / _ewRhsPrev;
			eta = 0.9 * ratio * ratio;
			const double safe = 0.9 * _ewEta * _ewEta;
			if(safe > 0.1 && safe > eta) eta = safe;
		}
		if(eta > kEtaMax) eta = kEtaMax;
		if(eta < epsInput) eta = epsInput;
		_ewEta = eta;
		_ewRhsPrev = rn;
		_bcg.Epsilon() = eta;
	}

	// Solve the system
	int ch;
	{ YUZEN_PROF("h.linsolve"); ch = _bcg.Solve(_crsIA, _Rhs,_correctVec, iter, resid); }
	_bcg.Epsilon() = epsInput;
	//  The boundary rows exactly, not to the tolerance of the solve: a
	//  mirror ghost a little off its mirror makes the capillary flux of an
	//  edge face, a third difference over it, a leak of volume
	//  (CRS_Matrix::SatisfyBCRows)
	if(ch == 0) _crsIA.SatisfyBCRows(_Rhs, _correctVec);
	if(ch != 0)
	{
		Global::Warning<<"Newton-Kantorovich: BiCGStab did not converge: ret = "<< ch << ", iter = "<<iter<< ", resid = "<< std::scientific <<resid<< nl;

		//--------------------------------------------------------------
		//  Throw away the solution vector of the failed solve.
		//
		//  The linear solver also uses _correctVec as its initial guess
		//  (solveWithGuess, in the case of Eigen). If the solve diverged and
		//  left a NaN in it, handing that to the next call makes every
		//  subsequent solve fail with a NaN, however far the time step is
		//  reduced: once a single step could not be solved, the time step
		//  would collapse to its lower limit and the run abort.
		//--------------------------------------------------------------
		_correctVec = 0;
	}

	return ch;
}

bool C_NewtonKantorovich::NKLoop(int &iter, double &resid_eq)
{
	// The Newton-Kantorovich iteration
	assert(_initialized);

	_nkIter = 0;				// at the top of a time step the Jacobian is always rebuilt
	_bcg_iter_step = 0;			// total linear iterations (a measurement counter; it does not affect the numerics)
	_coupledCorrection = 0;		// reset the correction of the coupled side
	_failCause.clear();			// what made this trial fail, if it does (FailedFixedStep)
	_phi_guess = _phi_old;		// start the estimate from the old step (inside EveryNKLoop the fields are updated from phi_guess)

	//--------------------------------------------------------------
	//  opt extrap (util/Opt.h): start from the film extrapolated linearly from the
	//  last two steps instead, phi^n + (dt / dt_old)(phi^n - phi^n-1), unless
	//  that leaves a cell that had liquid without any (then from phi^n as
	//  before). The relative convergence tests below stay relative to the
	//  first correction, now the error of the extrapolation: taken relative
	//  to the move of the extrapolation instead, the iteration stopped a
	//  pass earlier and the volume of the thick nozzle film drifted to
	//  VolErr 5e-11 (tools/analytic.py nozzle_balance; 7e-14 so).
	//--------------------------------------------------------------
	_extrapMove = 0;
	if(yuzen::opt::On(yuzen::opt::Extrap) && !_extrapBlocked && _phi_old2.size() == _phi_old.size())
	{
		const double dtNow = _timeStep->DT().now, dtOld = _timeStep->DT().old;
		const double r = (dtOld > 0) ? dtNow / dtOld : 0.0;
		bool ok = (r > 0);
		double move = 0;
		for(int i = 0; i < _phi_old.size() && ok; ++i)
		{
			const double g = _phi_old[i] + r * (_phi_old[i] - _phi_old2[i]);
			if(!_crsA.IsBC(i) && _phi_old[i] > 0 && !(g > 0)) ok = false;
			const double m = std::fabs(g - _phi_old[i]);
			if(m > move) move = m;
		}
		if(ok && move > 0)
		{
			for(int i = 0; i < _phi_old.size(); ++i)
				_phi_guess[i] = _phi_old[i] + r * (_phi_old[i] - _phi_old2[i]);
			_extrapMove = move;
		}
	}

	BeforeNKLoop();				// tell the coupled side that a time step is starting
	SetMatrixA();				// assemble the matrix
	ApplyA(_Ax_old, _phi_old);	// evaluate the explicit part

	if(Math::IsZero(_theta)) return SolveExplicit();	// explicit scheme: no system is solved

	//####### Loop ########
	double resid_eq_old;
	bool NK_done = false;
	double alpha_NK = 1.0;
	int ch;

	int outer_iter;

	//--------------------------------------------------------------
	//  ** Do not use the absolute value alone as the convergence test. **
	//
	//  The Newton iteration converges quadratically, so the correction falls
	//  like
	//      2.1e-03 -> 3.6e-05 -> 7.3e-07 -> 1.7e-08 -> 3.9e-10 -> 9.0e-11
	//  From there on, however, it hits a floor set by the truncation error
	//  and does nothing but wander erratically between 2e-11 and 1.6e-10.
	//  The height of that floor is roughly eps_mach * dt / dt_explicit, so it
	//  rises as the time step is made finer.
	//
	//  Testing on the absolute value _EpsilonNK alone (typically 1e-12) means
	//  that under conditions where the floor sits at 1e-10 the test is never
	//  met: the loop runs all the way to NewtonMaxIterations and reports "did not
	//  converge". The caller sees that and shrinks the time step, so the
	//  solution spins uselessly and loses its step size even though it had
	//  converged perfectly well. On a computation that resolves the contact
	//  line this made the run 100 times slower.
	//
	//  Two further tests were therefore added.
	//
	//  (1) A relative convergence test. Once the correction has fallen to
	//      _EpsilonRelNK times the first correction of that step it counts as
	//      converged (the default is 1e-6). This means "the nonlinearity of
	//      this time step has been resolved to six digits", so it depends
	//      neither on the scale of the solution nor on the step size.
	//
	//  (2) A stagnation test, taken on the number of iterations that fail to
	//      improve on the smallest value seen so far. Testing on the ratio
	//      |resid - resid_old| instead would not fire when the values scatter
	//      on the floor, where successive differences run 50 to 300 %.
	//--------------------------------------------------------------
	double resid_first = 0;		// the first correction of this step
	double resid_best  = 1e40;	// the smallest value seen so far
	int stall_count = 0;		// how many times the smallest value was not improved on
	const int stall_limit = 5;

	resid_eq_old = 1e40;
	for(outer_iter = 0; outer_iter < _MaxIterNK; outer_iter++)
	{
		bool everyOk;
		{ YUZEN_PROF("h.couple"); everyOk = EveryNKLoop(); }
		if(everyOk == false) 	{	NK_done = false; break;	}

		ch = Solve(_bcg_iter, _bcg_resid);
		_bcg_iter_step += _bcg_iter;

		//--------------------------------------------------------------
		//  The convergence test is taken on the larger of the film-height
		//  correction and the correction of the coupled side.
		//
		//  Looking at the film-height correction alone would let the loop
		//  leave while the concentration was still moving. Measuring the
		//  coupled correction as well means the loop leaves early under
		//  conditions where the coupling does not matter, and keeps iterating
		//  as long as necessary under conditions where it does.
		//
		//  Once it has converged, the height equation and the concentration
		//  equation are satisfied at the same time. In other words, what this
		//  iteration converges to is exactly the solution of the two solved
		//  together, fully implicitly.
		//--------------------------------------------------------------
		resid_eq = _correctVec.AbsMax();

		//  Output for when one wants to see whether it is the film height or
		//  the coupling that is holding convergence back. It is written to
		//  standard error only when YUZEN_DEBUG_RESIDUAL is set.
		//    [resid] iteration number, height correction, concentration correction
		if(_dbgResid)
			std::fprintf(stderr, "[resid] %d %.6e %.6e\n",
			             outer_iter, resid_eq, _coupledCorrection);

		//--------------------------------------------------------------
		//  Output for seeing where the correction levels off. It is written
		//  to standard error only when YUZEN_DEBUG_FLOOR is set.
		//
		//  The Newton correction does not fall all the way to machine
		//  precision. When the right-hand side
		//      Rhs = (phi_old - phi_guess) - dt[theta A phi_guess + (1-theta) A phi_old + B]
		//  is formed, A is a fourth-order difference, so its coefficients are
		//  of order 1/dx^4: terms of size |A phi| ~ 1e15 cancel against each
		//  other to produce an answer of size 1e5. The absolute error left
		//  behind is eps_mach * ||A|| * |phi|, and it is multiplied by dt as
		//  it enters Rhs. The floor on the correction is therefore about
		//
		//      floor ~ eps_mach * dt / dt_explicit,   dt_explicit ~ 3 mu dx^4 / H^3
		//
		//  that is, proportional to how far the step exceeds the stability
		//  limit of the explicit scheme. The accuracy of the linear solver
		//  (bcg_resid) has nothing to do with it; it measures 1e-25.
		//--------------------------------------------------------------
		static const bool dbgFloor = (std::getenv("YUZEN_DEBUG_FLOOR") != 0);
		if(dbgFloor)
		{
			double axmax = 0, rhsmax = 0, phimax = 0;
			const double dt_now = _timeStep->DT().now;
			for(int i = 0; i < _Rhs.size(); ++i)
			{
				if(fabs(_Ax_guess[i])  > axmax ) axmax  = fabs(_Ax_guess[i]);
				if(fabs(_Rhs[i])       > rhsmax) rhsmax = fabs(_Rhs[i]);
				if(fabs(_phi_guess[i]) > phimax) phimax = fabs(_phi_guess[i]);
			}
			std::fprintf(stderr,
			  "[floor] it=%d dt=%.3e |dphi|=%.4e |Rhs|=%.4e |dt*A.phi|=%.4e "
			  "|phi|=%.4e bcg_it=%d bcg_resid=%.3e\n",
			  outer_iter, dt_now, resid_eq, rhsmax, dt_now*axmax, phimax,
			  _bcg_iter, _bcg_resid);
		}

		if(_coupledCorrection > resid_eq) resid_eq = _coupledCorrection;

		if(ch != 0) // && resid_eq > max_resid_eq)
		{
			// Stop at once when the solve did not converge and the correction
			// is too large
			Info<<"*";
			NK_done = false;
			break;
		}
		if(outer_iter == 0) resid_first = resid_eq;

		//  Did this iteration improve on the smallest value so far?
		if(resid_eq < resid_best * 0.9) { resid_best = resid_eq; stall_count = 0; }
		else                            { stall_count++; }

		if(stall_count > stall_limit)
		{
			//----------------------------------------------------------
			//  The correction has stopped decreasing. Either it has hit the
			//  floor set by the accuracy of the linear solver, or something
			//  non-smooth is at work (the upwind direction switching, or
			//  evaporation and the concentration hitting their limits).
			//  Either way it will not decrease by iterating further, so the
			//  loop stops here.
			//
			//  If it has fallen far enough, the step counts as solved and
			//  the loop leaves quietly. Far enough is any of
			//    - the absolute test,
			//    - four decades below the first correction: the
			//      nonlinearity of the step resolved to four digits,
			//    - within 1000 eps |phi|: the film cannot be updated more
			//      finely than eps |phi| (ch-convergence, "Machine
			//      precision is the limit"; the floors measured there are
			//      about that), for a step whose first correction was
			//      already near it.
			//  Otherwise the iteration has stopped far from a solution: at
			//  the default tolerances the shipped contact-line cases
			//  (wet_static*) stopped at 0.13 to 160 times the first
			//  correction and the thick nozzle film with liquid supplied at
			//  its own concentration at 8e-4 to 1.9 times; none between
			//  1e-6 (the relative test) and 8e-4. Kept, such a step is not
			//  a solution of the step (a volume error of 2e-8 on the nozzle
			//  film), so it fails like a linear solve that did not
			//  converge: cut and solved again under step control (the first
			//  one warned, all counted), and it stops the run at a fixed
			//  step (FailedFixedStep).
			//  The four decades are not tied to NewtonEpsilonRel on
			//  purpose: with it at 1e-11 (NewtonEpsilon 1e-16) the same
			//  nozzle film levels off at 1e-11 to 1e-4 of the first
			//  correction in two steps of three (the limited upwinding at
			//  Courant numbers up to 150 is not smooth), and failing those
			//  would chase the relative test as an absolute-only test
			//  chases its threshold.
			//----------------------------------------------------------
			const double floorPhi = 1000.0 * std::numeric_limits<double>::epsilon() * _phi_guess.AbsMax();
			const bool solved = (resid_eq < _EpsilonNK)
			                 || (resid_first > 0 && resid_eq < resid_first * 1.0e-4)
			                 || (resid_eq <= floorPhi);
			if(_dbgResid)
				std::fprintf(stderr, "[stall] %d %.6e %.6e %.6e %d\n",
				             outer_iter, resid_eq, resid_first, floorPhi, (int)solved);
			if(!solved)
			{
				std::ostringstream why;
				why<<"the Newton correction stagnated far from convergence (iteration "<<outer_iter+1
				   <<", correction "<<std::scientific<<std::setprecision(6)<<resid_eq
				   <<" against "<<resid_first<<" in the first)";
				_failCause = why.str();
				_stagnationCuts++;
				//  (not from an attempt from the extrapolated film: TimeAdvance
				//  makes it again from phi^n)
				if(_adjustTimeStep && !_warnedStagnation && _extrapMove == 0)
				{
					Global::Warning<<"Newton-Kantorovich: "<<_failCause<<"; the step is cut and solved again"
					                 " (JobInfo.txt counts the steps cut for this; later ones are not reported here)."<<nl;
					_warnedStagnation = true;
				}
			}
			NK_done = solved;
			break;
		}

		if(resid_eq > resid_eq_old) alpha_NK = 0.5;
		else alpha_NK = 1.0;

		_phi_guess +=  _correctVec * alpha_NK;

		//  Convergence test: absolute (_EpsilonNK) or relative (_EpsilonRelNK)
		if(outer_iter >= _MinIterNK-1)
		{
			if(resid_eq < _EpsilonNK) { NK_done = true; break; }
			if(_EpsilonRelNK > 0 && resid_first > 0
			   && resid_eq < resid_first * _EpsilonRelNK) { NK_done = true; break; }
		}
		resid_eq_old = resid_eq;
	}
	iter = outer_iter+1;
	return NK_done;
}

int C_NewtonKantorovich::TimeAdvance(int &iter, double &resid)
{
	_phi_old2 = _phi_old;
	_phi_old = _phi_new;
	//  NKLoop assembles A at the start of the step and forms A h^n again
	//  before anything reads it; with opt assembly (util/Opt.h) this product with
	//  the A of the last step is not formed
	if(!yuzen::opt::On(yuzen::opt::Assembly)) ApplyA(_Ax_old, _phi_old);

	bool timeStepDone = false;
	double criterion;
	bool ch;
	iter = 0;
	do
	{
		//  An attempt that started from the extrapolated film and failed is
		//  made once more from phi^n, at the same step, before the step is
		//  cut or the run stopped: a source that turns on (a nozzle) is not
		//  in the last two steps, and from their extrapolation the Newton
		//  correction of the first step of tools/analytic.py nozzle_area
		//  stagnated at a fixed step. The first attempt leaves no count or
		//  warning behind.
		const long long cuts0 = _stagnationCuts;
		const bool warned0 = _warnedStagnation;
		ch = NKLoop(iter,resid);
		if(!ch && _extrapMove > 0)
		{
			_stagnationCuts   = cuts0;
			_warnedStagnation = warned0;
			_extrapBlocked = true;
			ch = NKLoop(iter,resid);
			_extrapBlocked = false;
		}

		if(_adjustTimeStep)
		{
			if(ch) criterion = CheckTimeStepCriterion();
			else criterion = -1;
			timeStepDone = TimeStepAdjustment(criterion);
		}
		else
		{
			//  With the step fixed there is no shorter one to retry at. A
			//  step whose iteration failed (the linear solve, no
			//  convergence in NewtonMaxIterations, or a correction that
			//  stagnated far from convergence) stops the run: taken as it stands --
			//  the guess of the failed iteration, often the film it started
			//  from -- the run would end normally (tools/proptest.py seed
			//  3094).
			if(!ch) FailedFixedStep("Newton-Kantorovich");
			timeStepDone = true;
		}

	}while(!timeStepDone);

	_phi_new = _phi_guess;

	bool afterOk;
	{ YUZEN_PROF("h.after"); afterOk = AfterNKLoop(); }
	if(afterOk == false) Info<<"*";

	return 0;
}

void C_NewtonKantorovich::MakeMatrix_IA()
{
	double delta_t = _timeStep->DT().now;
	int ix, jy;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(jy) if(yuzen::opt::Par(_crsJacobi.RowSize()))
#endif
	for(ix = 0; ix < _crsJacobi.RowSize(); ix++)
	{
		if(_crsJacobi.IsBC(ix))
		{
			_crsIA.Diag(ix) = _crsJacobi.Diag(ix);
			for(jy = 1; jy < _crsJacobi.NonZero(ix); jy++)
				_crsIA.Element(ix,jy) =  _crsJacobi.Element(ix,jy);
		}
		else
		{
			_crsIA.Diag(ix) = 1.0 + _theta * delta_t * _crsJacobi.Diag(ix);
			for(jy = 1; jy < _crsJacobi.NonZero(ix); jy++)
				_crsIA.Element(ix,jy) =  _theta * delta_t * _crsJacobi.Element(ix,jy);
		}
	}
}

void C_NewtonKantorovich::MakeRHS()
{
	double delta_t = _timeStep->DT().now;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_crsA.RowSize()))
#endif
	for(int ix = 0; ix < _crsA.RowSize(); ix++)
	{
		if(_crsA.IsBC(ix)) continue;
		_Rhs[ix] = _phi_old[ix] -_phi_guess[ix] - (delta_t * _theta * _Ax_guess[ix]) - (delta_t*(1.0 - _theta)*_Ax_old[ix])  - delta_t * _VecB[ix];
	}
}
