/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ImplicitSolver.cpp : implementation of the implicit time-integration base
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "Opt.h"
#include "YuzenUtil.h"
using namespace yuzen::util;

#include "ImplicitSolver.h"
#include "TimeSteps.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

C_ImplicitSolver::C_ImplicitSolver()
	:
	_nRows(1),
	_timeStep(0),
	_constructed(false),
	_initialized(false),
	_adjustTimeStep(false),
	_adjustTimeKeepSteps(50),
	_adjustTimeCount(0),
	_adjustTimeChangeFactor(0.5),
	_adjustTimeCriterion(1e-3),
	_adjustTimeScaleFloor(0.0),
	_implicit_type(CrankNicholson)
{
	CheckTheta();
}


void C_ImplicitSolver::Build(const int n, C_TimeSteps &ts, tImplicitType t, bool asjustTimeStep)
{
	_nRows = n;
	_adjustTimeStep = asjustTimeStep;
	_implicit_type = t;
	_timeStep = &ts;
	CheckTheta();
	Alloc();

}

void C_ImplicitSolver::CheckTheta()
{
	switch(_implicit_type)
	{
	case ForwardEuler:		_theta = 0.0; break;
	case CrankNicholson:	_theta = 0.5; break;
	case BackwardEuler:		_theta = 1.0; break;
	}
}

void C_ImplicitSolver::Alloc()
{
	_phi_old2.resize(_nRows);
	_phi_old.resize(_nRows);
	_phi_new.resize(_nRows);
	_Ax_old.resize(_nRows);
	_Rhs.resize(_nRows);
	_VecB.resize(_nRows);
	_crsA.Build(_nRows);
	_crsIA.Build(_nRows);

	SetFlagBC();
	SetMatrixA();				//establish the matrix structure

	_crsIA = _crsA;

	_bcg.Build(_crsA.RowSize(), _crsA.NonZero(), _crsA.NonZ_LU().first, _crsA.NonZ_LU().second);

	AllocAdd();

	_constructed = true;
}


void C_ImplicitSolver::AllocAdd()
{
	// Overridden by derived classes
}


void C_ImplicitSolver::SetX(Vector<double> &f)
{
	assert(_constructed);

	_phi_old = f;
	_phi_old2 = _phi_old;
	_phi_new = _phi_old;
}
void C_ImplicitSolver::SetX(ScalarField1D<> &f)
{
	assert(_constructed);

	int ix;
	for(Xedge xe = f.Grid().First<Xedge>(); xe != f.Grid().Last<Xedge>(); ++xe)
	{
		ix = f.Grid().ix(xe);
		_phi_old[ix] = f(xe);
	}
	_phi_old2 = _phi_old;
	_phi_new = _phi_old;
}
void C_ImplicitSolver::SetX(ScalarField2D<> &f)
{
	assert(_constructed);

	int ix, jy;
	for(Zface zf = f.Grid().First<Zface>(); zf != f.Grid().Last<Zface>(); ++zf)
	{
		ix = f.Grid().ix(zf);
		jy = f.Grid().iy(zf);
		_phi_old[jy*f.Grid().nx() + ix] = f(zf);
	}
	_phi_old2 = _phi_old;
	_phi_new = _phi_old;
}
void C_ImplicitSolver::SetX(ScalarField3D<> &f)
{
	assert(_constructed);

	int ix, jy, kz;
	int nx, ny;
	nx = f.Grid().nx();
	ny = f.Grid().ny();
	for(Cell cl = f.Grid().First<Cell>(); cl != f.Grid().Last<Cell>(); ++cl)
	{
		ix = f.Grid().ix(cl);
		jy = f.Grid().iy(cl);
		kz = f.Grid().iz(cl);
		_phi_old[ kz*(nx*ny) + jy*nx + ix] = f(cl);
	}
	_phi_old2 = _phi_old;
	_phi_new = _phi_old;
}

bool C_ImplicitSolver::TimeStepAdjustment(double criterion)
{
	double dt_max = _timeStep->DT().range_max;
	double dt_min = _timeStep->DT().range_min;
	double dt_next = _timeStep->DT().now;

	bool done = false;

	// The criterion has to exceed this before dt is allowed to grow
	const double allow_increase = 0.90;

	if(criterion > 0.0)
	{
		done = true;
		if(criterion > allow_increase && dt_next < dt_max)
		{
			_adjustTimeCount++;
			if(_adjustTimeCount >= _adjustTimeKeepSteps)
			{
				_adjustTimeCount = 0;
				dt_next /= _adjustTimeChangeFactor;
				if(dt_next > dt_max) dt_next = dt_max;
				//  Only asked for: the step just accepted stays DT().now
				//  until the clock has advanced by it (C_TimeSteps::
				//  operator++). Changed at once, DT().now would make the
				//  clock add the lengthened step, which has not been solved,
				//  and the end of the step (AfterNKLoop / AfterTimeStepTrial,
				//  the sums of History) would see it too.
				if(_timeStep->RequestNextDT(dt_next))
					Global::Info<<"dt+ ("<<std::scientific<<dt_next<<")"<<nl;
			}
		}
	}
	else
	{
		dt_next *= _adjustTimeChangeFactor;
		if(dt_next < dt_min)
		{
			Global::FatalError<<"TimeStepAdjustment(). Under flow: the requested step fell below Numerics.TimeStepControl.MinDt."<<nl;
			Global::FatalError.Abort();
		}

		if(_timeStep->ChangeDT(dt_next))
			Global::Info<<"dt- ("<<std::scientific<<dt_next<<")"<<nl;
	}
	return done;
}

void C_ImplicitSolver::FailedFixedStep(const char *who) const
{
	if(_failCause.empty())
		Global::FatalError<<who<<": the step did not converge at the fixed step (Warnings.txt says why), ";
	else
		Global::FatalError<<who<<": the step did not converge at the fixed step: "<<_failCause<<", ";
	Global::FatalError<<"and with Numerics.TimeStepControl.Enable = false it cannot be retried at a shorter one. "
	                    "Shorten Numerics.InitialDt or enable the time-step control."<<nl;
	Global::FatalError.Abort();
}

//////////////////////////////////////////////////////////////////
//  The "rate of change" used for the steady-state test
//
//  ** Do not base the test on the change per step. **
//  max|phi_new - phi_old| / max|phi_new| becomes small merely because the
//  time step is small, so a run whose step has collapsed is mistaken for a
//  "steady state" and stops without a word. On a computation that resolves
//  the contact line, with dt ~ 1e-11, it would stop after 12 steps.
//
//  What is returned instead is that quantity divided by the time step: the
//  relative rate of change per unit of (dimensionless) time. Being
//  independent of the step size, it gives the threshold a physical meaning,
//  namely "this is no longer moving". The threshold is set by
//  Numerics.SteadyStateCriterion.
//////////////////////////////////////////////////////////////////
double C_ImplicitSolver::ChangeRateOfVariable()
{
	if(_timeStep->TimeStep() < 5) return 1;

	//  The step that took _phi_old to _phi_new. This is called once the
	//  clock has advanced (CheckTermination, CalcInfo), when DT().now is
	//  already the step after it.
	const double dt = _timeStep->LastStepDT();
	if(Math::IsZero(dt)) return 1;

	Vector<double> vec(_phi_old.size());
	vec = _phi_new - _phi_old;
	const double scale = _phi_new.AbsMax();
	if(Math::IsZero(scale)) return 0;

	return vec.AbsMax() / (scale * dt);
}


//////////////////////////////////////////////////////////////////
//  The time-step test based on the local truncation error
//
//  The local truncation error is estimated from a second-order difference in
//  time, and what is returned is the margin against the tolerance th. A
//  positive value accepts the step; a negative one shrinks it.
//
//  The error is divided by the local value with a floor on the weight.
//
//  The weight w_i = |h_old_i| alone, which makes the error a purely local
//  relative quantity, breaks down on computations with a thin precursor
//  film. The precursor (H* = 0.02) is 1/50 of the drop itself (H ~ 1), so a
//  slight wobble inside the precursor carries 50 times the weight of one in
//  the body of the drop. The internal structure of the contact line affects
//  neither the spreading rate nor the contact angle, and yet it would decide
//  the time step: refining the grid by a factor of two (dx 0.6 -> 0.3 um)
//  drops dt by a factor of 56, and neither raising TimeStepControl.Criterion
//  by 20x nor raising InitialDt by 10x moves it.
//
//  The weight used is
//
//      w_i = max( |h_old_i| ,  floor * max|h_old| )
//
//  where floor is Numerics.TimeStepControl.ScaleFloor.
//    floor = 0    : the default (a purely local relative error)
//    floor = 1    : normalized by the global maximum alone (close to an
//                   absolute error)
//    0 < floor < 1: the body of the drop is still judged relatively, while
//                   only the overestimation in the thin regions is damped
//
//  It means "wherever the film is thinner than floor times the maximum, do
//  not judge it any more strictly than that", which doubles as the guard
//  against dividing by zero in a dry region (h -> 0).
//////////////////////////////////////////////////////////////////
double C_ImplicitSolver::LocalTruncationCriterion(const Vector<double> &phi_new) const
{
	const double dt0 = _timeStep->DT().old;
	const double dt  = _timeStep->DT().now;

	double th = _adjustTimeCriterion;
	if(Math::IsZero(th)) th = 1.0e-2;

	// When dt0 is zero, as on the first step, no test is possible. Report no
	// constraint.
	if(Math::IsZero(dt0)) return 1.0;

	//  The floor on the weight, as a fraction of the maximum of _phi_old.
	const double wfloor = _adjustTimeScaleFloor * _phi_old.AbsMax();

	Vector<double> Ek(phi_new.size());
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(phi_new.size()))
#endif
	for(int ii = 0; ii < phi_new.size(); ++ii)
	{
		const double hm = _phi_old2[ii];
		const double ho = _phi_old[ii];
		const double hp = phi_new[ii];

		double w = fabs(ho);
		if(w < wfloor) w = wfloor;

		const double denom = (dt0 + dt) * w;
		// w == 0, which means a completely dry cell with the floor at zero
		// as well, is left out of the error estimate.
		Ek[ii] = Math::IsZero(denom)
		       ? 0.0
		       : 2.0 * (dt/dt0) * (dt0*hp + dt*hm - (dt0+dt)*ho) / denom;
	}
	return (1.0 - Ek.AbsMax()/th);
}
