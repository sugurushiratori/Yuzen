/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  TimeSteps.cpp : implementation of the time, time-step and output-schedule management
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <fstream>
#include <cstdlib>
#include <limits>

#include "YuzenUtil.h"
using namespace yuzen::util;

#include "TimeSteps.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

C_OutputIntervals::C_OutputIntervals(const C_OutputIntervals& c)
{
	this->_key = c._key;
	this->_nextScheduled = c._nextScheduled;
	this->_outputCount = c._outputCount;
	this->_outputInterval = c._outputInterval;
	this->_previousRealOut = c._previousRealOut;
	this->_realign = c._realign;
	this->_lastWritten = c._lastWritten;
	this->_times = c._times;
	this->_idx = c._idx;
	this->_exact = c._exact;
	this->_listStart = c._listStart;
}

void C_TimeSteps::Build(double deltaT, double startTime, double endTime, bool adaption)
{
	_sw_adaptDT = adaption;

	_deltaT.now     = deltaT;
	_deltaT.old     = deltaT;
	_deltaT.initial = deltaT;
	_deltaT.range_max = 1e10;
	_deltaT.range_min = 1e-15;
	_deltaT.max_in_run = deltaT;
	_deltaT.min_in_run = deltaT;

	_startTime = startTime;
	_endTime   = endTime;

	_startTimeStep = 0;
	_endTimeStep   = std::numeric_limits<std::uint64_t>::max();

	_nowTime  = _startTime;
	_timeStep = _startTimeStep;
	_pendingDT = 0;
	_lastStepDT = 0;

	_constructed = true;
}

void C_TimeSteps::Add_OutputItem(const std::string& key, double interval)
{
	C_OutputIntervals item;
	item.Build(key, interval, _startTime);
	_outItemList.push_back(item);
}

void C_TimeSteps::Add_OutputList(const std::string& key, const std::vector<double>& times, bool exact)
{
	C_OutputIntervals item;
	item.BuildList(key, times, exact, _startTime);
	_outItemList.push_back(item);
}

bool C_TimeSteps::IsRunning() const
{
	if( _constructed
		&&(_startTime     <= _nowTime  && _nowTime  <= _endTime)
		&&(_startTimeStep <= _timeStep && _timeStep <= _endTimeStep) ) return true;

	else return false;
}

C_OutputIntervals&  C_TimeSteps::OutputItem(const std::string& key)
{
	tOutputItemList::iterator it = _outItemList.begin();
	while(it != _outItemList.end())
	{
		if( it->Key() == key) break;
		++it;
	}
	assert(it != _outItemList.end());

	return *it;
}

////////////////////////////////////////////////////////

bool C_TimeSteps::ChangeDT(double dt_new)
{
	if(dt_new < _deltaT.range_min || _deltaT.range_max < dt_new) return false;

	//  DT().old is left alone: it is the step last taken (see the header).
	_deltaT.now = dt_new;
	if(_deltaT.max_in_run < dt_new) _deltaT.max_in_run = dt_new;
	if(_deltaT.min_in_run > dt_new) _deltaT.min_in_run = dt_new;

	return true;
}

bool C_TimeSteps::RequestNextDT(double dt_next)
{
	if(dt_next < _deltaT.range_min || _deltaT.range_max < dt_next) return false;
	_pendingDT = dt_next;
	return true;
}

void C_TimeSteps::operator++()
{
	if(!_constructed) return;

	//  The clock first, by the step that was solved ...
	_nowTime += _deltaT.now;
	_timeStep++;
	_lastStepDT = _deltaT.now;
	_deltaT.old = _deltaT.now;

	//  ... then the step for the next one.
	if(_pendingDT > 0)
	{
		_deltaT.now = _pendingDT;
		if(_deltaT.max_in_run < _pendingDT) _deltaT.max_in_run = _pendingDT;
		if(_deltaT.min_in_run > _pendingDT) _deltaT.min_in_run = _pendingDT;
		_pendingDT = 0;
	}

	//  A step shortened to land on a listed output time (ApplyListedTimes):
	//  the next one goes back to the step the run was on, unless the step
	//  control asked for less than the step just solved, or the step was
	//  retried shorter still (both say the step is too long, not the list).
	//  The step solved stays what it was in DT().old and LastStepDT(): the
	//  error estimate of the next step differences over the real distance.
	if(_cutFrom > 0)
	{
		const double solved = _lastStepDT;
		if(solved >= _cutTo * (1.0 - 1e-12) && _deltaT.now >= solved * (1.0 - 1e-12))
		{
			double back = _cutFrom;
			if(back > _deltaT.range_max) back = _deltaT.range_max;
			if(back > _deltaT.now) _deltaT.now = back;
		}
		_cutFrom = 0;
		_cutTo   = 0;
	}
}

void C_TimeSteps::ApplyListedTimes()
{
	for(tOutputItemList::iterator it = _outItemList.begin(); it != _outItemList.end(); ++it)
	{
		if(!it->IsList()) continue;
		const double cap = it->ListCap(_nowTime);
		if(cap >= 1e299) continue;
		if(it->Exact())
		{
			if(cap <= 1e-9 * Math::abs(it->NextListed())) continue;	// there already, to rounding
			if(_deltaT.now > cap)
			{
				if(_cutFrom == 0) _cutFrom = _deltaT.now;
				_deltaT.now = cap;
				_cutTo = cap;
			}
		}
		else
		{
			//  A fixed step held down by a short gap goes back up as the
			//  gaps widen (an adaptive one is grown by the step control)
			if(!_sw_adaptDT)
			{
				double want = _deltaT.initial;
				if(want > _deltaT.range_max) want = _deltaT.range_max;
				if(_deltaT.now < want) _deltaT.now = (want < cap) ? want : cap;
			}
			if(_deltaT.now > cap)
			{
				_deltaT.now = cap;
				if(_deltaT.min_in_run > cap) _deltaT.min_in_run = cap;
			}
		}
	}
}
