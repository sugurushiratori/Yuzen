/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  SemiImplicit.cpp : implementation of the semi-implicit time integration
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "Opt.h"
#include "YuzenUtil.h"
#include "Profile.h"
using namespace yuzen::util;

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "SemiImplicit.h"
#include "TimeSteps.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif


void C_SemiImplicitEquation::Initialize()
{
	assert(_constructed);

	SetFlagBC();
	SetMatrixA();
	_crsA.MatVec(_Ax_old, _phi_old);
	_crsIA = _crsA;		//give it the same matrix structure
	MakeMatrix_IA();
	MakeRHS();

	_phi_old2 = _phi_old;
	_phi_new =  _phi_old;

	_initialized = true;
}

bool C_SemiImplicitEquation::SolveExplicit()
{
	// for Explicit
	double delta_t = _timeStep->DT().now;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_Ax_old.size()))
#endif
	for(int ix = 0; ix < _Ax_old.size(); ix++)
	{
		if( _crsA.IsBC(ix))
			 _phi_new[ix] = _phi_old[ix] -  _Ax_old[ix];
		else
			_phi_new[ix] = _phi_old[ix] - delta_t * (_Ax_old[ix] + _VecB[ix]);
	}
	return 0;
}


int C_SemiImplicitEquation::Solve(int &iter, double &resid)
{
	// Solve the system [ _crsIA * _phi_new = _Rhs ] once

	assert(_initialized);

	{ YUZEN_PROF("h.assemble"); SetMatrixA(); }	// assemble the matrix
	if(Math::IsZero(_theta))			// explicit scheme: no system is solved
		return SolveExplicit();
	{
		YUZEN_PROF("h.rhs");
		MakeMatrix_IA();					// assemble the matrix for the implicit time advance
		MakeRHS();						// evaluate the right-hand side
	}

	// Solve the system
	int ch;
	{ YUZEN_PROF("h.linsolve"); ch = _bcg.Solve(_crsIA, _Rhs, _phi_new, iter, resid); }
	if(ch != 0)
	{
		Global::Warning<<"Semi-Implicit: BiCGStab did not converge: ret = "<< ch << ", iter = "<<iter<< ", resid = "<< std::scientific <<resid<< nl;

		//  _phi_new is also the initial guess of the next solve. A solve that
		//  broke down leaves NaN in it, and every solve after it -- each retry
		//  of the step at a smaller dt -- would fail at once with a NaN, the
		//  step be cut down to MinDt and the run abort (tools/proptest.py).
		//  Start the next one from the unknowns at the start of the step
		//  instead, as the Newton side resets its guess
		//  (C_NewtonKantorovich::Solve).
		bool finite = true;
		for(int i = 0; i < _phi_new.size() && finite; i++) finite = std::isfinite(_phi_new[i]);
		if(!finite) _phi_new = _phi_old;
	}

	return ch;
}

bool C_SemiImplicitEquation::TimeStepTrial(int &iter, double &resid)
{
	// Solve one step's worth without advancing the time
	assert(_initialized);

	bool timeStepDone = false;
	bool failed = false;
	double criterion;
	int ch;

	do
	{
		_failCause.clear();
		if(EveryTimeStepTrial() == false) { AfterTimeStepTrial(); break; }

		ch = Solve(iter, resid);
		//  The rest of the step (AfterTimeStepTrial: in ConcZ the
		//  concentration, solved after the film) belongs to the trial: when
		//  it fails, the step has failed as surely as when the solve above
		//  did, and a concentration solve that did not converge is not taken
		//  as it stands.
		if(_adjustTimeStep)
		{
			if(ch == 0) criterion = CheckTimeStepCriterion();
			else criterion = -1;
			//  (for a step the test accepts; a rejected one is solved again,
			//  the concentration with it)
			if(criterion > 0 && !AfterTimeStepTrial()) criterion = -1;
			timeStepDone = TimeStepAdjustment(criterion);
		}
		else
		{
			//  Without step control there is no shorter step to retry at. A
			//  failed solve is reported as a failure, not as a success with
			//  its answer taken: the equation that owns the
			//  step then stops the run, naming the cause (FailedFixedStep),
			//  and an equation solved within the step of another (the
			//  well-mixed concentration, in every Newton pass of the film)
			//  fails that step, which the step control, when it is on, cuts
			//  and retries.
			timeStepDone = true;
			failed = (ch != 0) || !AfterTimeStepTrial();
		}

	}while(!timeStepDone);

	return timeStepDone && !failed;

}

int C_SemiImplicitEquation::TimeAdvance(int &iter, double &resid)
{
	_phi_old2 = _phi_old;
	_phi_old = _phi_new;

	_crsA.MatVec(_Ax_old, _phi_old);	// evaluate the explicit part

	bool ch = TimeStepTrial(iter, resid);
	return (ch ? 0 : -1);
}

void C_SemiImplicitEquation::MakeMatrix_IA()
{
	double delta_t = _timeStep->DT().now;
	int ix, jy;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(jy) if(yuzen::opt::Par(_crsA.RowSize()))
#endif
	for(ix = 0; ix < _crsA.RowSize(); ix++)
	{
		if(_crsA.IsBC(ix))
		{
			_crsIA.Diag(ix) = _crsA.Diag(ix);
			for(jy = 1; jy < _crsA.NonZero(ix); jy++)
				_crsIA.Element(ix,jy) =  _crsA.Element(ix,jy);
		}
		else
		{
			_crsIA.Diag(ix) = 1.0 + _theta * delta_t * _crsA.Diag(ix);
			for(jy = 1; jy < _crsA.NonZero(ix); jy++)
				_crsIA.Element(ix,jy) =  _theta * delta_t * _crsA.Element(ix,jy);
		}
	}
}

void C_SemiImplicitEquation::MakeRHS()
{
	double delta_t = _timeStep->DT().now;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_crsA.RowSize()))
#endif
	for(int ix = 0; ix < _crsA.RowSize(); ix++)
	{
		if(_crsA.IsBC(ix)) continue;
		_Rhs[ix] = _phi_old[ix] - delta_t*(1.0 - _theta)*_Ax_old[ix] - delta_t*_VecB[ix];
	}
}
