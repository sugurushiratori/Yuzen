/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  RestartState.h : the state a run carries from one step to the next that
*                   is not in the field files (field*.state)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _DEF_RESTART_STATE_INCLUDED_
#define _DEF_RESTART_STATE_INCLUDED_

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class C_TimeSteps;

//--------------------------------------------------------------------------
//  What a field file does not hold, but the next step depends on.
//
//  A field file holds h, s and c. A run also carries, from one step to the
//  next: the time and the step count, the time step and the one before it
//  (the adaptive control and its truncation-error estimate need both, and the
//  film of the step before), the counters of that control, the output
//  schedule, the running totals of History.dat (evaporated and outflowed
//  volume, the initial volumes the errors are taken against, the outflow rate
//  of the last step for the trapezoidal sum), the count of the steady-state
//  test, the explicit half of the concentration step (well-mixed) and the
//  concentration of the step before (ConcZ, whose corrector extrapolates c
//  from the last two steps). A restart from the field alone starts all of
//  that afresh, which changes the answer and resets the balances.
//
//  "dt" holds the step about to be taken and, twice, the step the last step
//  was solved with (the dt0 of the truncation estimate, and the step the
//  steady-state rate divides by). A file that lacks the third takes the
//  first in its place.
//
//  Every field output writes this file beside the field files, with the same
//  number: field00012.state next to field00012.vts / .dat. A restart looks
//  for it next to Restart.FileName and, when it is there and was written at
//  Restart.RestartTime, continues from it exactly.
//
//  Layout: text lines "key v1 v2 ..." (scalars, written with 17 significant
//  digits, which a double survives exactly), "text key word", lines
//  "vector name n" announcing
//  the arrays, a line "binary", and then the arrays as raw little-endian
//  doubles in the order announced. Unknown keys are ignored, so later
//  versions can add to it without breaking older files.
//--------------------------------------------------------------------------
class C_RestartState
{
public:
	C_RestartState() {}

	void Clear() { _scalars.clear(); _vectors.clear(); _order.clear(); _texts.clear(); }

	void Set(const std::string &key, double v) { _scalars[key] = std::vector<double>(1, v); }
	void Set(const std::string &key, const std::vector<double> &v) { _scalars[key] = v; }
	/** @brief An array (written in binary). */
	void SetVector(const std::string &key, const std::vector<double> &v);
	template<class V> void SetVectorOf(const std::string &key, const V &v)
	{
		std::vector<double> w((size_t)v.size());
		for(size_t i = 0; i < w.size(); i++) w[i] = v[(int)i];
		SetVector(key, w);
	}

	/** @brief A word (no blanks), such as the model name. */
	void SetText(const std::string &key, const std::string &v) { _texts[key] = v; }
	std::string Text(const std::string &key) const
	{ std::map<std::string, std::string>::const_iterator it = _texts.find(key); return it == _texts.end() ? "" : it->second; }

	bool   Has(const std::string &key) const { return _scalars.count(key) != 0; }
	bool   HasVector(const std::string &key) const { return _vectors.count(key) != 0; }
	/** @brief The k-th value of a key (dflt if absent). */
	double Get(const std::string &key, size_t k = 0, double dflt = 0) const;
	const std::vector<double>& Vector(const std::string &key) const;
	/** @brief Copy an array into v, if it is there and of v's size. */
	template<class V> bool GetVectorInto(const std::string &key, V &v) const
	{
		if(!HasVector(key)) return false;
		const std::vector<double> &w = Vector(key);
		if((size_t)v.size() != w.size()) return false;
		for(size_t i = 0; i < w.size(); i++) v[(int)i] = w[i];
		return true;
	}

	//--------------------------------------------------------------
	//  The fields themselves, every point of the storage (ghosts included),
	//  as the run held them. A field file holds them in the units of the
	//  output (micrometers when Output.Dimensioned) and, for ASCII ConcZ,
	//  without the outer ghost layer, so reading it back gives the state to
	//  rounding, not to the bit; the copy here makes the continuation exact.
	//--------------------------------------------------------------
	template<class F> void SetField(const std::string &key, const F &f)
	{
		const int n = f.Grid().imax();
		std::vector<double> w((size_t)n);
		for(int i = 1; i <= n; i++) w[(size_t)(i-1)] = f[i];
		SetVector(key, w);
	}
	/**
	 * @brief Replace a field read from a field file by the copy kept here,
	 *        if the two agree where the field file has it.
	 * @param key     the name the copy was kept under (SetField)
	 * @param f       the field read from the field file
	 * @param inside  per linear index (1..imax), whether the point is one
	 *                the field file holds (the real cells)
	 * @param tol     largest difference, relative to the largest value, that
	 *                is still the same state written to fewer digits
	 * @param diff    the difference found
	 * @return 1 replaced; 0 no copy here (or not of this size), f untouched;
	 *         -1 the field file holds something else (it was edited), f
	 *         untouched
	 */
	template<class F> int RefineField(const std::string &key, F &f, const std::vector<char> &inside,
	                                  double tol, double &diff) const
	{
		diff = 0;
		const int n = f.Grid().imax();
		if(!HasVector(key) || Vector(key).size() != (size_t)n || inside.size() != (size_t)n + 1) return 0;
		const std::vector<double> &w = Vector(key);
		double scale = 0;
		for(int i = 1; i <= n; i++)
		{
			if(!inside[(size_t)i]) continue;
			const double d = w[(size_t)(i-1)] - f[i];
			diff  = std::max(diff, d < 0 ? -d : d);
			scale = std::max(scale, f[i] < 0 ? -f[i] : f[i]);
		}
		if(scale > 0) diff /= scale;
		if(!(diff <= tol)) return -1;
		for(int i = 1; i <= n; i++) f[i] = w[(size_t)(i-1)];
		return 1;
	}

	bool Write(const std::string &path) const;
	/** @brief false when the file is missing or not a state file. */
	bool Read(const std::string &path);

	/**
	 * @brief The state file that goes with a field file: the extension
	 *        (.vts or .dat) replaced by .state, or .state appended.
	 */
	static std::string PathFor(const std::string &fieldFile);

	//--------------------------------------------------------------
	//  The part every model has: the clock and the output schedule.
	//--------------------------------------------------------------
	/** @brief Record the time, step, dt and the schedule of each output item. */
	void PutClock(C_TimeSteps &ts, const std::vector<std::string> &items);
	/**
	 * @brief Put them back.
	 * @param ts          the clock to set
	 * @param items       the output items whose schedules to restore
	 * @param restoreDt   take dt and the dt before it from the file (adaptive
	 *                    time step); otherwise the step of the input stands
	 * @param fieldItem   the item whose count is the file number, which
	 *                    Restart.FileIndex sets instead
	 * @param fieldCount  that count (Restart.FileIndex)
	 * An item whose interval was changed in the input since starts its
	 * schedule afresh at the restart time.
	 */
	void TakeClock(C_TimeSteps &ts, const std::vector<std::string> &items, bool restoreDt,
	               const std::string &fieldItem, int fieldCount) const;

	/**
	 * @brief Find and check the state file of a restart.
	 * @param fieldFile  Restart.FileName
	 * @param model      the model name the file must have been written by
	 * @param tRestart   Restart.RestartTime made dimensionless
	 * @param message    what was decided, for the log
	 * @return true when the file is there, of this model, and was written at
	 *         tRestart (to the 12 digits a .pvd gives it, or to six): the
	 *         run then continues from it. Otherwise only the field is taken,
	 *         as a restart always did.
	 */
	bool Load(const std::string &fieldFile, const std::string &model, double tRestart, std::string &message);
	/**
	 * @brief Whether the values read can be taken as they stand: every value
	 *        finite, the counters whole numbers in range, dt positive.
	 * @param bad  set to what is wrong when false is returned
	 */
	bool Sane(std::string &bad) const;
	/** @brief Whether the last Load() found a state file at all. */
	bool Found() const { return _found; }

private:
	std::map<std::string, std::vector<double> > _scalars;
	std::map<std::string, std::vector<double> > _vectors;
	std::vector<std::string> _order;
	std::map<std::string, std::string> _texts;
	bool _found = false;
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif
