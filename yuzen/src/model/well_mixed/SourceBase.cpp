/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  SourceBase.cpp : the dimension-independent part of the source. For the
*                   details, see the top of SourceBase.h
*
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "SourceBase.h"
#include <cmath>

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

using namespace yuzen::util;

//------------------------------------------------------------------
//  Read the input and non-dimensionalize the lengths, the times and the
//  flow rate, each by its own scale.
//  If the non-dimensionalization does not match the grid's, the nozzle
//  position ends up displaced on the grid.
//------------------------------------------------------------------
void SourceBase::BuildParams(Input::I_Physics &inp_phys, NonDimensionalNumbers *nonDim, C_TimeSteps *ts)
{
	_NonDim   = nonDim;
	_timeStep = ts;

	_enabled      = inp_phys.Source.Enable;
	_Radius       = inp_phys.Source.Radius / _NonDim->ExtentScale();
	_flowRate     = inp_phys.Source.InletRate / (_NonDim->HeightScale() / _NonDim->TimeScale());
	_Conc         = inp_phys.Source.Conc;
	_concSupplied = (inp_phys.Source.ConcType == "supplied");
	_timeConstant = inp_phys.Source.TimeConstant / _NonDim->TimeScale();

	_armLength    = inp_phys.Source.ArmLength / _NonDim->ExtentScale();
	_armX         = inp_phys.Source.ArmCenter[0] / _NonDim->ExtentScale();
	_armY         = inp_phys.Source.ArmCenter[1] / _NonDim->ExtentScale();
	_rotationRate = inp_phys.Source.RotationRate * _NonDim->TimeScale();
	_startAngle   = inp_phys.Source.StartAngle;
	//  Times are divided by the time scale, like TimeConstant (multiplied by
	//  it, a RampTime of 1 s would mean T0^2 / (1 s) instead of 1 s).
	_slopeTime    = inp_phys.Source.RampTime / _NonDim->TimeScale();
	_mainTime     = inp_phys.Source.PlateauTime / _NonDim->TimeScale();
	_numSubRegion = inp_phys.Source.SubRegions;
}

//------------------------------------------------------------------
//  The nozzle position at the tip of the rotating arm.
//  The rotation rate ramps up over RampTime, stays constant for PlateauTime,
//  and then falls over another RampTime.
//------------------------------------------------------------------
double SourceBase::EvalTime() const
{
	return _useEvalTime ? _evalTime : _timeStep->NowTime();
}

void SourceBase::NozzlePos(double &xn, double &yn) const
{
	double tt = EvalTime();
	double t0 = _timeStep->StartTime();

	//  The angle is the integral of the rotation rate, omega f(t), where f
	//  rises linearly from 0 to 1 over Ts, holds 1 for Tm and falls back to 0
	//  over another Ts. omega f(t) (t - t0), the rate of the moment times
	//  the whole time elapsed, would turn the arm twice as far as it should
	//  during the ramp-up, and during the ramp-down, where f falls while t
	//  keeps growing, swing the arm backwards at up to omega (Ts + Tm) / Ts
	//  and leave it at StartAngle once the rotation has stopped (with Ts = 0,
	//  jump straight back at t = Tm).
	const double Ts = _slopeTime;
	const double Tm = _mainTime;
	const double t  = Math::Max(0.0, tt - t0);
	double swept;		// the integral of f from 0 to t
	if(t <= Ts)                swept = (Ts > 0) ? 0.5 * t * t / Ts : 0.0;
	else if(t <= Ts + Tm)      swept = 0.5 * Ts + (t - Ts);
	else if(t <= 2.0 * Ts + Tm)
	{
		const double tau = t - Ts - Tm;
		swept = 0.5 * Ts + Tm + tau - 0.5 * tau * tau / Ts;
	}
	else                       swept = Ts + Tm;

	// position of nozzle
	double theta = _startAngle + _rotationRate * swept;
	xn = _armX + _armLength * cos(theta);
	yn = _armY + _armLength * sin(theta);
}

//------------------------------------------------------------------
//  The flow rate at the current time. It ramps up as an exponential
//  relaxation with time constant tau.
//------------------------------------------------------------------
double SourceBase::FlowRateNow() const
{
	double tt = EvalTime();

	// actual flow rate
	double tau = _timeConstant;
	double coeff_Tau = 1.0;

	if(!Math::IsLessEpsilon(tau))
	{
		coeff_Tau = 1.0 - exp(- tt / tau);
	}

	return _flowRate * coeff_Tau;
}
