/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_source.cpp : the motion of the nozzle on its rotating arm
*
*  The arm turns at omega f(t), with f rising linearly from 0 to 1 over
*  RampTime Ts, holding 1 for PlateauTime Tm and falling back to 0 over another
*  Ts. Its angle is therefore StartAngle + omega * (integral of f), which is
*  compared here with a numerical integral of f, and with the angle the arm
*  must stop at, StartAngle + omega (Ts + Tm).
*
*  An angle of omega f(t) t would give twice the turn during the ramp-up, a
*  backward swing during the ramp-down, and the arm back at StartAngle once
*  the rotation has stopped. None of the regression cases turns the arm, so
*  nothing else would notice.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "doctest/doctest.h"

#include <cmath>

#include "SourceBase.h"
#include "TimeSteps.h"

using namespace yuzen;

namespace {

//  Exposes the nozzle position and sets the arm up without an input file.
class ArmProbe : public SourceBase
{
public:
	ArmProbe(C_TimeSteps *ts, double omega, double theta0, double Ts, double Tm)
	{
		_timeStep = ts;
		_armLength = 1.0; _armX = 0.0; _armY = 0.0;
		_rotationRate = omega; _startAngle = theta0;
		_slopeTime = Ts; _mainTime = Tm;
	}
	double Angle(double t)
	{
		SetEvalTime(t);
		double x, y;
		NozzlePos(x, y);
		return std::atan2(y, x);
	}
};

double RateFactor(double t, double Ts, double Tm)
{
	if(t < 0) return 0;
	if(t < Ts) return t / Ts;
	if(t < Ts + Tm) return 1.0;
	if(t < 2*Ts + Tm) return 1.0 - (t - Ts - Tm) / Ts;
	return 0.0;
}

//  The angle wrapped into (-pi, pi], for comparing with atan2
double Wrap(double a)
{
	const double pi = std::acos(-1.0);
	a = std::fmod(a + pi, 2*pi);
	if(a < 0) a += 2*pi;
	return a - pi;
}

}	// namespace

TEST_CASE("Source: the arm angle is the integral of the rotation rate")
{
	C_TimeSteps ts;
	ts.Build(1e-3, 0.0, 10.0, false);

	const double omega = 0.7, theta0 = 0.3, Ts = 0.8, Tm = 1.5;
	ArmProbe arm(&ts, omega, theta0, Ts, Tm);

	//  Midpoint rule on a fine grid, compared at a set of times that cover
	//  all four phases
	const int n = 200000;
	const double T = 2*Ts + Tm + 1.0, h = T / n;
	double swept = 0;
	int next = 0;
	const double checks[] = {0.2, 0.8, 1.3, 2.3, 2.9, 3.1, 3.5, T};
	for(int i = 0; i < n && next < 8; ++i)
	{
		swept += RateFactor((i + 0.5) * h, Ts, Tm) * h;
		if((i + 1) * h >= checks[next] - 1e-12)
		{
			CAPTURE(checks[next]);
			CHECK(arm.Angle((i + 1) * h) == doctest::Approx(Wrap(theta0 + omega * swept)).epsilon(1e-6));
			++next;
		}
	}
	CHECK(next == 8);

	//  Once the rotation has stopped the arm stays where it stopped
	const double stop = Wrap(theta0 + omega * (Ts + Tm));
	CHECK(arm.Angle(2*Ts + Tm + 0.5) == doctest::Approx(stop).epsilon(1e-12));
	CHECK(arm.Angle(2*Ts + Tm + 5.0) == doctest::Approx(stop).epsilon(1e-12));
}

TEST_CASE("Source: with no ramp the arm turns at omega and stops at PlateauTime")
{
	C_TimeSteps ts;
	ts.Build(1e-3, 0.0, 10.0, false);

	const double omega = 0.5, theta0 = -0.2, Tm = 2.0;
	ArmProbe arm(&ts, omega, theta0, 0.0, Tm);

	CHECK(arm.Angle(0.0) == doctest::Approx(theta0).epsilon(1e-12));
	CHECK(arm.Angle(1.0) == doctest::Approx(Wrap(theta0 + omega * 1.0)).epsilon(1e-12));
	CHECK(arm.Angle(Tm)  == doctest::Approx(Wrap(theta0 + omega * Tm)).epsilon(1e-12));
	CHECK(arm.Angle(Tm + 3.0) == doctest::Approx(Wrap(theta0 + omega * Tm)).epsilon(1e-12));
}
