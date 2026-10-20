/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CourantCheck.h : the Courant number of a limited convection scheme,
*                   watched over the run
*
*  The limited schemes (Numerics.Convection Koren, VanLeer) keep what they
*  carry bounded, with Crank-Nicolson, only for a Courant number up to one:
*  beyond it the explicit half of the step is no longer a positive
*  combination of the values it starts from (tools/advection.py: Koren
*  overshoots a step in c by 16 % at 1.5). The check only reports; it does
*  not shorten the step or change any result.
*
*  The model hands it, after every step, the largest rate R (1/time) over
*  the cells of the state the step ended on; the Courant number of the step
*  is dt max(R_start, R_end), R_start the rate the step before ended on.
*  The first step above one is reported in Warnings.txt as it happens, and
*  at the end of the run the number of such steps and the largest Courant
*  number; the largest is written to the run log and JobInfo.txt as well.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _COURANT_CHECK_INCLUDED_
#define _COURANT_CHECK_INCLUDED_

#include "YuzenUtil.h"
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class C_CourantCheck
{
public:
	C_CourantCheck() : _on(false), _rate(0), _max(0), _over(0), _steps(0) {}

	/**
	 * @brief Watch one equation.
	 * @param key     the input key of the scheme (for the messages),
	 *                e.g. "Numerics.Convection.Concentration"
	 * @param scheme  its value as given ("Koren", "VanLeer")
	 * @param what    what is carried, e.g. "the concentration"
	 * @param on      false: the scheme is not a limited one, or nothing is
	 *                carried; nothing is watched or reported
	 */
	void Build(const std::string &key, const std::string &scheme, const std::string &what, bool on)
	{ _key = key; _scheme = scheme; _what = what; _on = on; }

	bool On() const { return _on; }

	/** @brief The rate of the state the first step starts from. */
	void Start(double rate) { _rate = rate; }

	/**
	 * @brief The step just taken.
	 * @param rate  the largest rate of the state it ended on (dimensionless)
	 * @param dt    the step (dimensionless)
	 * @param step  the step's number in this run, from 1
	 * @param time  the time the step ended at, in seconds
	 */
	void Step(double rate, double dt, long long step, double time)
	{
		if(!_on) return;
		const double c = dt * std::max(_rate, rate);
		_rate = rate;
		++_steps;
		if(c > _max) _max = c;
		if(!(c > 1.0)) return;
		if(_over++ == 0)
		{
			std::ostringstream m;
			m<<std::setprecision(3)<<"Courant number "<<c<<" at step "<<step<<" (t = "<<time
			 <<" s) with "<<_key<<" = \""<<_scheme<<"\": a limited scheme keeps "<<_what
			 <<" free of over- and undershoots only up to 1. Reported once; the count of such steps is given"
			 <<" at the end of the run.";
			yuzen::util::Global::Warning<<m.str()<<nl;
		}
	}

	/**
	 * @brief At the end of the run: the count, and the largest number.
	 * @param tag  the equation, as JobInfo.txt names it ("conc", "height")
	 */
	void Report(const std::string &tag) const
	{
		if(!_on) return;
		std::ostringstream v;
		v<<std::setprecision(3)<<_max;
		if(_over > 0)
			yuzen::util::Global::Warning<<"Courant number above 1 in "<<(int)_over<<" of "<<(int)_steps
			               <<" steps (largest "<<v.str()<<") with "<<_key<<" = \""<<_scheme<<"\"."<<nl;
		yuzen::util::Global::Info<<"Largest Courant number of "<<_key<<" = \""<<_scheme<<"\": "<<v.str()<<nl;
		yuzen::util::Global::JobInfo.Add("MaxCourant_" + tag, _max, "Largest Courant number <" + tag + ">");
	}

private:
	std::string _key, _scheme, _what;
	bool _on;
	double _rate;		// the rate the last step ended on
	double _max;		// the largest Courant number so far
	long long _over;	// the steps above 1
	long long _steps;
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_COURANT_CHECK_INCLUDED_
