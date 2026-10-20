/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_timesteps.cpp : tests for the time-step management
*
*  test_timesteps.cpp  --  regression tests for numerics/TimeSteps.h
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "doctest/doctest.h"
#include "TimeSteps.h"
#include <vector>
#include <algorithm>

using namespace yuzen;

TEST_CASE("C_TimeSteps: basic time advance")
{
	C_TimeSteps ts;
	ts.Build(0.1, 0.0, 1.0, false);

	CHECK(ts.StartTime() == doctest::Approx(0.0));
	CHECK(ts.EndTime()   == doctest::Approx(1.0));
	CHECK(ts.NowTime()   == doctest::Approx(0.0));
	CHECK(ts.DT().now    == doctest::Approx(0.1));
	CHECK(ts.TimeStep()  == 0u);
	CHECK(ts.IsRunning() == true);

	int steps = 0;
	while(ts.IsRunning() && steps < 1000) { ++ts; ++steps; }

	// IsRunning() tests _nowTime <= _endTime with "less than or equal", so
	// the step at t = EndTime is executed as well: 11 steps from 0.0 to
	// 1.0.
	CHECK(steps == 11);
	CHECK(ts.TimeStep() == 11u);
	CHECK(ts.NowTime() == doctest::Approx(1.1));
	CHECK(ts.IsRunning() == false);
}

TEST_CASE("C_TimeSteps: EndTimeStep is past StartTimeStep")
{
	C_TimeSteps ts;
	ts.Build(0.1, 0.0, 1.0, false);
	CHECK(ts.StartTimeStep() == 0u);
	CHECK(ts.EndTimeStep() > ts.StartTimeStep());
}

TEST_CASE("C_TimeSteps: output intervals")
{
	C_TimeSteps ts;
	ts.Build(0.1, 0.0, 1.0, false);
	ts.Add_OutputItem("field", 0.5);

	int outs = 0, steps = 0;
	while(ts.IsRunning() && steps < 1000)
	{
		if(ts.TimeForOut("field")) { ++outs; ts.UpdateItem("field"); }
		++ts; ++steps;
	}
	// Output is expected around t = 0, 0.5 and 1.0. We allow 2 to 4 times
	// to absorb differences in the implementation.
	CHECK(outs >= 2);
	CHECK(outs <= 4);
}

TEST_CASE("C_TimeSteps: ChangeDT rejects a value outside the range")
{
	C_TimeSteps ts;
	ts.Build(0.1, 0.0, 1.0, true);
	ts.SetMinDT(1e-6);
	ts.SetMaxDT(1.0);

	CHECK(ts.ChangeDT(0.05) == true);
	CHECK(ts.DT().now == doctest::Approx(0.05));

	CHECK(ts.ChangeDT(1e-9) == false);	// below the lower limit
	CHECK(ts.ChangeDT(10.0) == false);	// above the upper limit
	CHECK(ts.DT().now == doctest::Approx(0.05));	// unchanged once rejected
}

//  The times one output item is written at, stepping by dt from the clock's
//  present time to its end; the schedule as Update() leaves it.
static std::vector<double> OutputTimes(C_TimeSteps &ts, const char *key)
{
	std::vector<double> t;
	int steps = 0;
	while(ts.IsRunning() && steps++ < 100000)
	{
		if(ts.TimeForOut(key)) { t.push_back(ts.NowTime()); ts.UpdateItem(key); }
		++ts;
	}
	return t;
}

TEST_CASE("C_OutputIntervals: a restart without a state file keeps to the schedule")
{
	//  Straight through: every 0.1 with dt = 0.01.
	C_TimeSteps u;
	u.Build(0.01, 0.0, 1.005, false);
	u.Add_OutputItem("field", 0.1);
	const std::vector<double> tu = OutputTimes(u, "field");
	REQUIRE(tu.size() == 11);

	//  Restarted at 0.23 from the field alone. A schedule left at zero would
	//  make the restart write every half interval until it had caught up.
	C_TimeSteps r;
	r.Build(0.01, 0.0, 1.005, false);
	r.Add_OutputItem("field", 0.1);
	r.Set(0.23);
	r.OutputItem("field").RestartAt(0.23);
	const std::vector<double> tr = OutputTimes(r, "field");
	REQUIRE(tr.size() == 9);			// 0.23 once more, then 0.3 ... 1.0
	CHECK(tr[0] == doctest::Approx(0.23));
	for(size_t k = 1; k < tr.size(); k++)
		CHECK(tr[k] == doctest::Approx(0.2 + 0.1 * (double)k));

	//  Restarted at an output time, given as a .pvd gives it (six digits,
	//  just short of it): not written a second time one step later.
	C_TimeSteps p;
	p.Build(0.01, 0.0, 1.005, false);
	p.Add_OutputItem("field", 0.1);
	p.Set(0.2999999);
	p.OutputItem("field").RestartAt(0.2999999);
	const std::vector<double> tp = OutputTimes(p, "field");
	REQUIRE(tp.size() == 8);
	CHECK(tp[1] == doctest::Approx(0.4));
}

TEST_CASE("C_TimeSteps: Restore puts back the clock, and MaxSteps counts this run")
{
	C_TimeSteps ts;
	ts.Build(0.01, 0.0, 10.0, true);
	ts.SetMaxSteps(5);
	ts.Restore(0.5, 50, 0.02, 0.015, 0.015);
	CHECK(ts.NowTime() == 0.5);
	CHECK(ts.TimeStep() == 50u);
	CHECK(ts.DT().now == 0.02);
	CHECK(ts.DT().old == 0.015);
	int steps = 0;
	while(ts.IsRunning() && steps < 100) { ++ts; ++steps; }
	CHECK(steps == 5);
}

TEST_CASE("C_TimeSteps: an output off the schedule at the end leaves the schedule alone")
{
	C_TimeSteps ts;
	ts.Build(0.01, 0.0, 1.0, false);
	ts.Add_OutputItem("field", 0.1);
	for(int k = 0; k < 23; k++)
	{
		if(ts.TimeForOut("field")) ts.UpdateItem("field");
		++ts;
	}
	const double next = ts.OutputItem("field").NextScheduled();
	const double prev = ts.OutputItem("field").PreviousRealOut();
	ts.UpdateItemAtEnd("field");		// t = 0.23, not due
	CHECK(ts.WrittenAt("field"));
	CHECK(ts.OutputItem("field").NextScheduled() == next);
	CHECK(ts.OutputItem("field").PreviousRealOut() == prev);
	CHECK(ts.OutputItem("field").Count() == 4u);
}

//  The order of the advance: the clock adds the step that was solved, and
//  only then takes up the step the control asked for. Changing DT().now as
//  soon as a step is accepted would make the clock add the lengthened step,
//  which has not been solved, and run ahead.
TEST_CASE("C_TimeSteps: the clock adds the step solved, then takes up the next")
{
	C_TimeSteps ts;
	ts.Build(0.25, 0.0, 100.0, true);
	ts.SetMinDT(1e-6);
	ts.SetMaxDT(10.0);

	//  Step 1, accepted, and the step after it lengthened.
	CHECK(ts.RequestNextDT(0.5) == true);
	CHECK(ts.DT().now == 0.25);			// still the step being solved
	CHECK(ts.PendingDT() == 0.5);
	++ts;
	CHECK(ts.NowTime() == 0.25);			// not 0.5
	CHECK(ts.DT().now == 0.5);
	CHECK(ts.DT().old == 0.25);
	CHECK(ts.LastStepDT() == 0.25);
	CHECK(ts.PendingDT() == 0.0);

	//  Step 2, rejected once and solved again at 0.375 (the clock does not
	//  move), then accepted. The retry follows step 1: DT().old stays the
	//  step step 1 was taken with, not the attempt that was rejected.
	CHECK(ts.ChangeDT(0.375) == true);
	CHECK(ts.NowTime() == 0.25);
	CHECK(ts.LastStepDT() == 0.25);
	CHECK(ts.DT().old == 0.25);
	++ts;
	CHECK(ts.NowTime() == 0.625);
	CHECK(ts.DT().now == 0.375);			// nothing asked for: kept
	CHECK(ts.LastStepDT() == 0.375);
	CHECK(ts.DT().old == 0.375);

	//  Step 3, unchanged: DT().old is the step before, not the value before
	//  the last change.
	++ts;
	CHECK(ts.NowTime() == 1.0);
	CHECK(ts.LastStepDT() == 0.375);
	CHECK(ts.DT().old == 0.375);

	//  A request outside the limits is refused and changes nothing.
	CHECK(ts.RequestNextDT(20.0) == false);
	CHECK(ts.PendingDT() == 0.0);
	CHECK(ts.DT().max_in_run == 0.5);
	CHECK(ts.DT().min_in_run == 0.25);
}

TEST_CASE("C_TimeSteps: with a varying step the time is the sum of the steps solved")
{
	C_TimeSteps ts;
	ts.Build(1e-3, 0.0, 1.0, true);
	ts.SetMinDT(1e-9);
	ts.SetMaxDT(0.05);
	double sum = 0;
	int k = 0;
	while(ts.IsRunning() && k < 100000)
	{
		//  A pattern like that of the control: a cut now and then (a retry
		//  of the same step), a lengthening every third accepted step.
		if(k % 7 == 3) ts.ChangeDT(ts.DT().now * 0.75);
		const double solved = ts.DT().now;
		if(k % 3 == 2) ts.RequestNextDT(std::min(solved / 0.75, ts.MaxDT()));
		++ts;
		sum += solved;
		CHECK(ts.LastStepDT() == solved);
		k++;
	}
	CHECK(ts.NowTime() == sum);			// the same additions, in the same order
	CHECK(ts.NowTime() > 1.0);
	CHECK(ts.NowTime() < 1.0 + ts.MaxDT());
}

TEST_CASE("C_TimeSteps: Restore drops a step asked for and not taken up")
{
	C_TimeSteps ts;
	ts.Build(0.01, 0.0, 10.0, true);
	ts.RequestNextDT(0.02);
	ts.Restore(0.5, 50, 0.03, 0.015, 0.015);
	++ts;
	CHECK(ts.NowTime() == doctest::Approx(0.53));
	CHECK(ts.DT().now == 0.03);
}

TEST_CASE("C_TimeSteps: listed output times, landed on (ExactTimes)")
{
	C_TimeSteps ts;
	ts.Build(0.3, 0.0, 2.0, false);
	ts.SetMaxDT(10.0);
	ts.SetMinDT(1e-12);
	const std::vector<double> listed = {0.05, 0.1, 0.7, 1.0};
	ts.Add_OutputList("field", listed, true);
	CHECK(ts.OutputItem("field").IsList());
	CHECK(ts.OutputItem("field").Interval() == -4.0);

	std::vector<double> written, steps;
	int n = 0;
	while(ts.IsRunning() && n < 100)
	{
		if(ts.TimeForOut("field")) { written.push_back(ts.NowTime()); ts.UpdateItem("field"); }
		ts.ApplyListedTimes();
		steps.push_back(ts.DT().now);
		++ts; ++n;
	}
	//  The start, then each listed time to rounding, once
	REQUIRE(written.size() == 5u);
	CHECK(written[0] == 0.0);
	for(size_t i = 0; i < listed.size(); i++)
		CHECK(written[i + 1] == doctest::Approx(listed[i]).epsilon(1e-14));
	//  Shortened to 0.05, 0.05, then back to 0.3 (not left at 0.05): 0.4,
	//  0.7 (cut), 1.0 (0.3), then 0.3 on
	REQUIRE(steps.size() >= 6u);
	CHECK(steps[0] == doctest::Approx(0.05));
	CHECK(steps[1] == doctest::Approx(0.05));
	CHECK(steps[2] == doctest::Approx(0.3));
	CHECK(steps[3] == doctest::Approx(0.3));
	CHECK(steps[4] == doctest::Approx(0.3));
	CHECK(steps[5] == doctest::Approx(0.3));
	//  The step solved is what the clock and LastStepDT() saw
	CHECK(ts.NowTime() == doctest::Approx(2.2));
}

TEST_CASE("C_TimeSteps: listed output times, at the closest step")
{
	C_TimeSteps ts;
	ts.Build(0.3, 0.0, 2.0, false);
	ts.SetMaxDT(10.0);
	ts.SetMinDT(1e-12);
	//  0.05 and 0.06 are closer together than any step: written once
	ts.Add_OutputList("field", {0.05, 0.06, 1.0}, false);
	std::vector<double> written, steps;
	int n = 0;
	while(ts.IsRunning() && n < 200)
	{
		if(ts.TimeForOut("field")) { written.push_back(ts.NowTime()); ts.UpdateItem("field"); }
		ts.ApplyListedTimes();
		steps.push_back(ts.DT().now);
		++ts; ++n;
	}
	//  Held to the gap between the listed times around it: 0.05 first, then
	//  0.01, then 0.94, then the step of the run again (0.3) past the last
	CHECK(steps[0] == doctest::Approx(0.05));
	CHECK(written.size() == 4u);		// 0, ~0.05, ~0.06, ~1.0
	CHECK(written.back() == doctest::Approx(1.0).epsilon(0.5));
	CHECK(steps.back() == doctest::Approx(0.3));
}

TEST_CASE("C_OutputIntervals: a listed schedule restored from a state file")
{
	C_OutputIntervals it;
	it.BuildList("field", {0.1, 0.2, 0.4}, true, 0.0);
	it.Restore(0.15, 0.2, 3, 0.15);			// was waiting for 0.2
	CHECK(it.NextScheduled() == 0.2);
	CHECK(it.ListCap(0.15) == doctest::Approx(0.05));
	it.RestoreListAfter(0.25, 3, 0.25);		// no state: the first after 0.25
	CHECK(it.NextScheduled() == 0.4);
}
