/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  TimeSteps.h : management of the time, the time step and the output schedule
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _DEF_TIMESTEPS_INCLUDED_
#define _DEF_TIMESTEPS_INCLUDED_

#include <map>
#include <string>
#include <list>
#include <limits>
#include <cstdint>
#include <cmath>
#include <vector>

#include "Vector.h"
#include "MyMath.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif


//--------------------------------------------------------------------------
//  Manages the output schedule for one kind of output.
//
//  This is a small state machine that realizes "write every so many seconds"
//  on top of a time step that varies. The step does not necessarily land
//  exactly on a scheduled time, so the output is written at **the step
//  closest to the next scheduled time** - neither after overshooting it nor
//  before reaching it. On top of that, at least half an interval is left
//  since the last output, so that a step finer than the output interval does
//  not produce a burst of outputs.
//
//  Field, Probe, History and CalcInfo each hold one of these.
//--------------------------------------------------------------------------
class C_OutputIntervals
{

public:
	C_OutputIntervals() :
		 _key("undefined"),
		  _outputInterval(1e40),
		  _previousRealOut(-1e40),
		  _nextScheduled(0),
		  _outputCount(0),
		  _realign(false),
		  _lastWritten(-1e40),
		  _idx(0),
		  _exact(false)
		 {}

	C_OutputIntervals(const C_OutputIntervals& c);

	virtual ~C_OutputIntervals(){}
	explicit C_OutputIntervals(const std::string& key) :
		 _key(key),
		  _outputInterval(1e40),
		  _previousRealOut(-1e40),
		  _nextScheduled(0),
		  _outputCount(0),
		  _realign(false),
		  _lastWritten(-1e40),
		  _idx(0),
		  _exact(false)
		 {}

	/**
	 * @brief Set the kind of output, the interval between outputs and the
	 *        time the schedule starts from.
	 * @param key             name of the kind, such as "Field"
	 * @param outputInterval  interval between outputs (dimensionless time)
	 * @param startTime       time the first scheduled output is due
	 */
	void Build(const std::string& key, double outputInterval, double startTime)
	{
		_key = key;
		_outputInterval  = outputInterval;
		_previousRealOut = -1e40;
		_nextScheduled   = startTime;
		_realign         = false;
		_lastWritten     = -1e40;
		_times.clear();
		_idx   = 0;
		_exact = false;
	}

	/**
	 * @brief Set the kind of output to a list of times instead of an
	 *        interval (Output.FieldData.Times or LogSchedule).
	 * @param times      the times, ascending, all after startTime
	 *                   (dimensionless)
	 * @param exact      land a step on each of them (ExactTimes) rather
	 *                   than write at the step closest to it
	 * @param startTime  the start of the run: the state there is written
	 *                   first, as with an interval
	 *
	 * Interval() then returns minus the number of times, which no interval
	 * can be: a restart state file tells the two kinds apart by it.
	 */
	void BuildList(const std::string& key, const std::vector<double>& times, bool exact, double startTime)
	{
		Build(key, -(double)times.size(), startTime);
		_times = times;
		_exact = exact;
		_idx   = 0;
	}

	/** @brief Whether the outputs follow a list of times. */
	bool IsList() const {return !_times.empty();}
	/** @brief Whether a step is landed on each listed time. */
	bool Exact() const {return _exact;}

	/**
	 * @brief The largest step the list allows from `now`: with ExactTimes,
	 *        what is left to the next listed time; otherwise the gap
	 *        between the listed times the next output falls in (as an
	 *        interval caps the step). A huge value past the last time.
	 */
	double ListCap(double now) const
	{
		if(_times.empty() || _idx >= _times.size()) return 1e300;
		if(_exact) return _times[_idx] - now;
		const double prev = (_idx == 0) ? _listStart : _times[_idx - 1];
		return _times[_idx] - prev;
	}
	/** @brief The next listed time not yet written (huge past the last one). */
	double NextListed() const {return (_idx < _times.size()) ? _times[_idx] : 1e300;}

	/** @brief Name of the kind of output, such as "Field". */
	const std::string& Key() const {return _key;}

	/**
	 * @brief How many outputs have been written so far.
	 * @return The count, which becomes the sequence number in the file name.
	 */
	unsigned int Count()    const {return _outputCount;}
	/** @brief Carry the file sequence number over on restart. */
	void SetCount(int count) {_outputCount = count;}

	/** @brief Interval between outputs. */
	double Interval() const {return _outputInterval;}

	/**
	 * @brief Whether this step is the one to write an output at.
	 * @param nowTime  the current time
	 * @param deltaT   the time step about to be taken
	 * @return true when both conditions hold: at least half an interval has
	 *         passed since the last output (which stops bursts of outputs),
	 *         and now is closer to the scheduled time than the next step
	 *         would be.
	 *
	 * Calling this does not change any state. Call Update() once the output
	 * has actually been written.
	 */
	bool TimeForOut(double nowTime, double deltaT) const
	{
		if(!_times.empty())
		{
			//  A list: the start, then each listed time. With ExactTimes a
			//  step lands on it (to rounding, or just past it when the
			//  remainder was below the smallest step); otherwise the step
			//  closest to it, as with an interval but without the
			//  half-interval guard, since listed times may be closer
			//  together than a step.
			if(_nextScheduled >= 1e300) return false;
			const double tol = 1e-9 * Math::abs(_nextScheduled);
			if(nowTime >= _nextScheduled - tol) return true;
			if(_exact) return false;
			return Math::abs(nowTime - _nextScheduled) < Math::abs(nowTime + deltaT - _nextScheduled);
		}
		if( (( nowTime - _previousRealOut) > 0.5*_outputInterval) &&
			( Math::abs( nowTime - _nextScheduled) < Math::abs( nowTime + deltaT - _nextScheduled) )) return true;
		else return false;
	}

	/**
	 * @brief Whether the state at this exact time has already been written.
	 *
	 * Update() records the time it was called with, so the comparison is an
	 * exact one between two copies of the same value. The final output at
	 * the end of a run asks this before writing, so that a run whose last
	 * scheduled output happens to land on the end time does not get the
	 * same state twice.
	 */
	bool WrittenAt(double nowTime) const {return _lastWritten == nowTime;}

	/**
	 * @brief Record that an output was written and move the schedule on to
	 *        the next due time.
	 */
	void Update(double nowTime)
	{
		_previousRealOut = nowTime;
		_lastWritten     = nowTime;
		if(!_times.empty())
		{
			//  On to the first listed time after this one: a step that
			//  passed several (closer together than the step) writes once.
			//  The time written for is passed in any case: without
			//  ExactTimes the step closest to it can fall short of it.
			if(!_exact && _idx < _times.size() && _times[_idx] <= _nextScheduled) _idx++;
			AdvanceListPast(nowTime);
			_realign = false;
			_outputCount++;
			return;
		}
		if(_realign)
		{
			//  The first multiple of the interval after the restart time
			//  (which a .pvd gives to 12 digits: a time just short of a
			//  multiple is that multiple), with nothing held against it.
			_nextScheduled = (std::floor(nowTime / _outputInterval + 0.01) + 1.0) * _outputInterval;
			_previousRealOut = -1e40;
			_realign = false;
		}
		else _nextScheduled += _outputInterval;
		_outputCount++;
	}

	/**
	 * @brief Record an output written off the schedule (the one at the end
	 *        of a run, at a time the schedule did not ask for): counted, but
	 *        the schedule is left as it was, so that a run restarted from
	 *        that output writes at the times an uninterrupted run would.
	 */
	void UpdateOffSchedule(double nowTime) { _lastWritten = nowTime; _outputCount++; }

	/**
	 * @brief Start the schedule at a restart that has no state file to take
	 *        it from: the state restarted from is written at once, and the
	 *        outputs after it fall on the multiples of
	 *        the interval again, from the first one after the restart time
	 *        (where an uninterrupted run writes them, to within a step).
	 *
	 * Left where Build() put it, at time zero, the schedule would make every
	 * step after a restart at t whose distance from the last output exceeds
	 * half an interval an output, until it had caught up with t: twice as
	 * many files as asked for, off the times an uninterrupted run writes at.
	 */
	void RestartAt(double now)
	{
		_previousRealOut = -1e40;
		_nextScheduled   = now;
		_realign         = true;
		_lastWritten     = -1e40;
		//  A list goes on from the first listed time after the restart
		//  (Update, after the state restarted from is written)
	}

	//  The schedule as a restart file keeps it, and putting it back.
	double PreviousRealOut() const {return _previousRealOut;}
	double LastWritten()     const {return _lastWritten;}
	double NextScheduled()   const {return _nextScheduled;}
	void Restore(double previousRealOut, double nextScheduled, unsigned int count, double lastWritten)
	{
		_lastWritten     = lastWritten;
		_previousRealOut = previousRealOut;
		_nextScheduled   = nextScheduled;
		_outputCount     = count;
		_realign         = false;
		if(!_times.empty())
		{
			//  The listed time the run being continued was waiting for
			_idx = 0;
			while(_idx < _times.size() && _times[_idx] < nextScheduled * (1.0 - 1e-12)) _idx++;
			_nextScheduled = NextListed();
		}
	}
	/**
	 * @brief A list continued from a state that did not keep it: the
	 *        outputs go on from the first listed time after `now`.
	 */
	void RestoreListAfter(double now, unsigned int count, double lastWritten)
	{
		_lastWritten     = lastWritten;
		_previousRealOut = now;
		_outputCount     = count;
		_realign         = false;
		AdvanceListPast(now);
	}

private:
	void AdvanceListPast(double nowTime)
	{
		const double tol = 1e-9 * Math::abs(nowTime);
		if(_idx == 0) _listStart = nowTime;
		while(_idx < _times.size() && _times[_idx] <= nowTime + tol) _idx++;
		_nextScheduled = NextListed();
	}

private:
	std::string _key;
	double _outputInterval;
	double _previousRealOut;
	double _nextScheduled;
	unsigned int _outputCount;
	bool _realign;		// see RestartAt
	double _lastWritten;	// the time of the last output, on the schedule or off it
	std::vector<double> _times;	// the listed times (empty: an interval)
	size_t _idx;			// the first listed time not yet written
	bool _exact;			// land a step on each listed time
	double _listStart = 0;		// the time of the first output (the start)

private:
	void operator=(const C_OutputIntervals& c);
};

//--------------------------------------------------------------------------
//  Management of the time, the time step and the output schedule.
//
//  Every time and every step held here is **dimensionless**. Seconds given
//  in the input file are divided by the time scale T0 on the model side
//  before they are passed in.
//
//  The step is changed through ChangeDT() and is capped by the upper and
//  lower limits (SetMaxDT / SetMinDT). How the step is actually decided is a
//  judgment made on the ImplicitSolver side; this class only changes it as
//  it is told.
//--------------------------------------------------------------------------
class C_TimeSteps
{
public:
	//  The complete set of time steps. Both now and old are kept because the
	//  estimate of the local truncation error uses a two-level difference
	//  (ImplicitSolver::LocalTruncationCriterion).
	//
	//  now is the step being solved, from the first attempt at it to the
	//  advance of the clock (++), which adds exactly that step to the time;
	//  after ++ it is the step the next one will be attempted with. old is
	//  the step the last accepted step was solved with (= LastStepDT()):
	//  the distance between the two states the truncation estimate
	//  differences against -- not the value now had before its last change
	//  (after a step left unchanged, the step before the change; on the
	//  retry of a rejected attempt, the rejected step), over which the
	//  estimate would difference a step that has not been taken, leaving a
	//  first-order term (proportional to the rate of change times the
	//  mismatch) in what should be a second difference and holding the step
	//  down.
	struct t_deltaTime
	{
		double now;			// the step being solved (the next one, after ++)
		double old;			// the step the last accepted step was solved with
		double initial;		// initial  time interval
		double range_max;
		double range_min;
		double max_in_run;
		double min_in_run;
	};

	C_TimeSteps() : _constructed(false), _sw_adaptDT(false), _dT_changed(false), _pendingDT(0), _lastStepDT(0) {}
private:
	C_TimeSteps(const C_TimeSteps& cts);
	C_TimeSteps& operator=(const C_TimeSteps& cts);
public:
	C_TimeSteps(double deltaT, double startTime, double endTime, bool adaption) : _dT_changed(false), _pendingDT(0), _lastStepDT(0)
	{
		Build(deltaT, startTime, endTime, adaption);
	}

	virtual ~C_TimeSteps()
	{
		_outItemList.clear();
		_constructed=false;
	}

	/**
	 * @brief Set the initial step and the start and end times.
	 * @param deltaT     initial time step (dimensionless)
	 * @param startTime  start time (dimensionless)
	 * @param endTime    end time (dimensionless)
	 * @param adaption   whether the step is allowed to vary
	 */
	void Build(double deltaT, double startTime, double endTime, bool adaption);

	/** @brief Whether the time step is allowed to vary. */
	const bool& EnableAdaption() const {return _sw_adaptDT;}
	      bool& EnableAdaption()       {return _sw_adaptDT;}

	/** @brief The whole set of time steps: now, old, initial and the limits. */
	const t_deltaTime& DT() const {return _deltaT;}

	/**
	 * @brief Upper and lower limits on the time step.
	 *
	 * The upper limit is normally MaxDtFactor * dx^4 and the lower limit is
	 * MinDt from the input file. Going below the lower limit makes
	 * TimeStepAdjustment abort the run.
	 */
	void SetMaxDT(double v){ _deltaT.range_max = v;}
	void SetMinDT(double v){ _deltaT.range_min = v;}
	double MaxDT() const {return _deltaT.range_max;}
	double MinDT() const {return _deltaT.range_min;}

	//  Times, all dimensionless.
	double StartTime() const {return _startTime;}
	double EndTime()   const {return _endTime;}
	double NowTime()   const {return _nowTime;}

	//  Step numbers.
	std::uint64_t StartTimeStep() const {return _startTimeStep;}
	std::uint64_t EndTimeStep()   const {return _endTimeStep;}
	std::uint64_t TimeStep()      const {return _timeStep;}

	/**
	 * @brief Stop after this many steps, whatever the end time says.
	 * @param n  maximum number of steps to take, counted from the start step
	 *
	 * Corresponds to the input key Numerics.MaxSteps. Build() leaves the
	 * limit at the largest representable value, so a run is bounded by the
	 * end time alone unless this is called.
	 */
	void SetMaxSteps(std::uint64_t n)
	{
		//  IsRunning() tests _timeStep <= _endTimeStep, so stopping after
		//  exactly n steps means setting the limit to n-1.
		_endTimeStep = (n == 0) ? _startTimeStep : _startTimeStep + n - 1;
	}

	/**
	 * @brief Register a kind of output.
	 * @param key       "Field" / "Probe" / "GeneralHistory" / "CalcInfo"
	 * @param interval  interval between outputs of that kind
	 */
	void Add_OutputItem(const std::string& key, double interval);
	/**
	 * @brief Register a kind of output written at listed times.
	 * @param times  ascending, after the start time (dimensionless)
	 * @param exact  land a step on each (ApplyListedTimes)
	 */
	void Add_OutputList(const std::string& key, const std::vector<double>& times, bool exact);
	/**
	 * @brief Look up a registered kind of output.
	 * @return The schedule for that kind. Passing a name that was never
	 *         registered aborts the run.
	 */
	C_OutputIntervals&  OutputItem(const std::string& key);
	/** @brief Whether that kind of output is due at the current time. */
	bool TimeForOut(const std::string& key) {return OutputItem(key).TimeForOut(_nowTime, _deltaT.now);}
	/**
	 * @brief Record that an output of that kind was written.
	 *
	 * Must be called whenever TimeForOut() returned true.
	 */
	void UpdateItem(const std::string& key) {return OutputItem(key).Update(_nowTime);}
	/**
	 * @brief Record an output written at the end of a run: on the schedule
	 *        when the schedule had it due (TimeForOut), off it otherwise.
	 */
	void UpdateItemAtEnd(const std::string& key)
	{
		if(TimeForOut(key)) UpdateItem(key);
		else OutputItem(key).UpdateOffSchedule(_nowTime);
	}

	/** @brief Whether this kind of output already holds the present time. */
	bool WrittenAt(const std::string& key) {return OutputItem(key).WrittenAt(_nowTime);}

	/** @brief Carry the file sequence number of one kind over on restart. */
	void SetCount(const std::string& key, int count){OutputItem(key).SetCount(count);}

	/** @brief Whether neither the end time nor the end step has been reached yet. */
	bool IsRunning() const;

	/** @brief Set the current time directly, used when restarting. */
	void Set(double now){_nowTime = now;}

	/**
	 * @brief Put back the clock of a run that is being continued exactly
	 *        (from a field*.state): the time, the step count, the time step,
	 *        the one before it and the step the last step was solved with.
	 *
	 * Numerics.MaxSteps keeps counting the steps of this run: the limit is
	 * moved on by the steps already taken.
	 */
	void Restore(double now, std::uint64_t step, double dtNow, double dtOld, double lastStep)
	{
		_nowTime = now;
		if(_endTimeStep != std::numeric_limits<std::uint64_t>::max())
			_endTimeStep += step - _startTimeStep;
		_startTimeStep = step;
		_timeStep = step;
		_deltaT.now = dtNow;
		_deltaT.old = dtOld;
		_pendingDT  = 0;
		_lastStepDT = lastStep;
	}

	/** @brief Set the start time. */
	void SetStartTime(double st){_startTime = st;}

	/**
	 * @brief Advance one step: add the step that was solved (DT().now) to
	 *        the time, increment the step number, and only then take up the
	 *        step asked for by RequestNextDT() for the next one.
	 *
	 * LastStepDT() becomes the step just taken. Does nothing until Build()
	 * has been called.
	 *
	 * A lengthened step takes effect here, after the clock has advanced by
	 * the step solved. Lengthened as soon as a step is accepted, before the
	 * clock advances, it would make the clock add a step that has not been
	 * solved: the reported time would run ahead (by 0.37 % at the end of
	 * leveling_1d_reg), and with it every time-dependent term of the steps
	 * after.
	 */
	void operator++();
	void operator++(int){ operator++();}

	/**
	 * @brief Change the step of the attempt being made: the first attempt
	 *        (the initial step capped) or a retry after a rejected attempt.
	 * @param dt_new  the requested step; outside the upper and lower limits
	 *                it is refused
	 * @return true when the step actually changed, so that the caller can
	 *         skip logging when it did not.
	 *
	 * The clock has not moved, and DT().old stays the step last taken: the
	 * retry follows that step, not the attempt that was rejected.
	 */
	bool ChangeDT(double dt_new);

	/**
	 * @brief Ask for the step after this one (the step being solved has
	 *        been accepted and may be lengthened): taken up by ++, after
	 *        the clock has advanced by the step that was solved.
	 * @return false when the step is outside the limits (nothing changes).
	 */
	bool RequestNextDT(double dt_next);

	/**
	 * @brief Before a step: the listed output times (BuildList) cap it.
	 *
	 * Without ExactTimes the step is held to the gap between the listed
	 * times the next output falls in, as an interval holds it to the
	 * interval. With ExactTimes a step that would pass the next listed time
	 * is shortened to land on it, and the step taken up after it is the one
	 * the run was on before (operator++): a listed time is not a reason to
	 * keep the steps short. A remainder below 1e-9 of the time is rounding
	 * and is not stepped to.
	 */
	void ApplyListedTimes();

	/** @brief The step asked for by RequestNextDT() and not yet taken up (0: none). */
	double PendingDT() const {return _pendingDT;}

	/**
	 * @brief The step the last accepted step was solved with: what took the
	 *        solution from its previous state to the present one, once the
	 *        clock has advanced (0 before the first step).
	 */
	double LastStepDT() const {return _lastStepDT;}

protected:
	bool _constructed;
	bool _sw_adaptDT;
	bool _dT_changed;

	t_deltaTime _deltaT{};
	double _pendingDT;	// the next step asked for (RequestNextDT), 0 when none
	double _cutFrom = 0;	// the step a listed time shortened (ApplyListedTimes), 0: none
	double _cutTo   = 0;	// what it was shortened to
	double _lastStepDT;	// the step the last accepted step was solved with

	std::uint64_t _startTimeStep = 0;
	std::uint64_t _endTimeStep = 0;
	std::uint64_t _timeStep = 0;

	double _nowTime = 0.0;
	double _startTime = 0.0;
	double _endTime = 0.0;

	typedef std::list<C_OutputIntervals> tOutputItemList;
	tOutputItemList _outItemList;

};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_DEF_TIMESTEPS_INCLUDED_
