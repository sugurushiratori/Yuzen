/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  RestartState.cpp : reading and writing field*.state (see RestartState.h)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "RestartState.h"
#include "YuzenUtil.h"
#include "TimeSteps.h"

#include <algorithm>
#include <cerrno>
#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

static const char *kMagic = "# Yuzen restart state";

void C_RestartState::SetVector(const std::string &key, const std::vector<double> &v)
{
	if(!_vectors.count(key)) _order.push_back(key);
	_vectors[key] = v;
}

double C_RestartState::Get(const std::string &key, size_t k, double dflt) const
{
	std::map<std::string, std::vector<double> >::const_iterator it = _scalars.find(key);
	if(it == _scalars.end() || it->second.size() <= k) return dflt;
	return it->second[k];
}

const std::vector<double>& C_RestartState::Vector(const std::string &key) const
{
	static const std::vector<double> empty;
	std::map<std::string, std::vector<double> >::const_iterator it = _vectors.find(key);
	return it == _vectors.end() ? empty : it->second;
}

bool C_RestartState::Write(const std::string &path) const
{
	//  Written to a temporary name and renamed, so that a run killed while
	//  writing never leaves a state file that looks whole and is not.
	const std::string tmp = path + ".tmp";
	errno = 0;
	{
		std::ofstream f(tmp.c_str(), std::ios::binary);
		if(!f.good()) return false;
		f << kMagic << " 1\n";
		char buf[64];
		for(std::map<std::string, std::vector<double> >::const_iterator it = _scalars.begin(); it != _scalars.end(); ++it)
		{
			f << it->first;
			for(size_t k = 0; k < it->second.size(); k++)
			{
				std::snprintf(buf, sizeof(buf), " %.17g", it->second[k]);
				f << buf;
			}
			f << "\n";
		}
		for(size_t k = 0; k < _order.size(); k++)
			f << "vector " << _order[k] << " " << _vectors.find(_order[k])->second.size() << "\n";
		for(std::map<std::string, std::string>::const_iterator it = _texts.begin(); it != _texts.end(); ++it)
			f << "text " << it->first << " " << it->second << "\n";
		f << "binary\n";
		for(size_t k = 0; k < _order.size(); k++)
		{
			const std::vector<double> &v = _vectors.find(_order[k])->second;
			if(!v.empty()) f.write((const char*)v.data(), (std::streamsize)(v.size() * sizeof(double)));
		}
		//  Closed before it is tested: the end of the buffer is written by
		//  close(), and on a full disk that is where the write fails. (Tested
		//  before closing, a state file cut short was renamed into place.)
		f.close();
		if(f.fail()) { std::remove(tmp.c_str()); return false; }
	}
	std::string why;
	if(yuzen::util::System::ReplaceFile(tmp, path, why)) return true;
	std::remove(tmp.c_str());
	return false;
}

bool C_RestartState::Read(const std::string &path)
{
	Clear();
	std::ifstream f(path.c_str(), std::ios::binary);
	if(!f.good()) return false;
	std::string line;
	if(!std::getline(f, line) || line.compare(0, std::string(kMagic).size(), kMagic) != 0) return false;

	std::vector<std::pair<std::string, size_t> > arrays;
	bool binary = false;
	while(std::getline(f, line))
	{
		if(line == "binary") { binary = true; break; }
		std::istringstream is(line);
		std::string key;
		if(!(is >> key)) continue;
		if(key == "text")
		{
			std::string name, val;
			is >> name >> val;
			_texts[name] = val;
			continue;
		}
		if(key == "vector")
		{
			std::string name; size_t n = 0;
			if(!(is >> name >> n)) return false;
			arrays.push_back(std::make_pair(name, n));
			continue;
		}
		std::vector<double> v;
		std::string tok;
		while(is >> tok)
		{
			//  A word that is not a number is a damaged file, not a zero
			char *end = 0;
			v.push_back(std::strtod(tok.c_str(), &end));
			if(end == tok.c_str() || *end != '\0') { Clear(); return false; }
		}
		_scalars[key] = v;
	}
	if(!binary) return false;
	//  The lengths of the arrays come from the file: one longer than what is
	//  left of it is a damaged file, and is refused rather than allocated as
	//  it stands (a length near 2^64 would throw std::length_error out of the
	//  run).
	const std::streampos here = f.tellg();
	f.seekg(0, std::ios::end);
	const std::streampos fileEnd = f.tellg();
	f.seekg(here);
	if(here < 0 || fileEnd < here) { Clear(); return false; }
	size_t left = (size_t)(fileEnd - here);
	for(size_t k = 0; k < arrays.size(); k++)
	{
		if(arrays[k].second > left / sizeof(double)) { Clear(); return false; }
		left -= arrays[k].second * sizeof(double);
	}
	for(size_t k = 0; k < arrays.size(); k++)
	{
		std::vector<double> v(arrays[k].second);
		if(!v.empty()) f.read((char*)v.data(), (std::streamsize)(v.size() * sizeof(double)));
		if(!f) { Clear(); return false; }
		SetVector(arrays[k].first, v);
	}
	return true;
}

std::string C_RestartState::PathFor(const std::string &fieldFile)
{
	const size_t slash = fieldFile.find_last_of("/\\");
	const size_t dot = fieldFile.find_last_of('.');
	if(dot != std::string::npos && (slash == std::string::npos || dot > slash))
		return fieldFile.substr(0, dot) + ".state";
	return fieldFile + ".state";
}

void C_RestartState::PutClock(C_TimeSteps &ts, const std::vector<std::string> &items)
{
	Set("time", ts.NowTime());
	Set("step", (double)ts.TimeStep());
	std::vector<double> dt(3);
	dt[0] = ts.DT().now; dt[1] = ts.DT().old; dt[2] = ts.LastStepDT();
	Set("dt", dt);
	for(size_t k = 0; k < items.size(); k++)
	{
		C_OutputIntervals &it = ts.OutputItem(items[k]);
		std::vector<double> v(5);
		v[0] = it.PreviousRealOut(); v[1] = it.NextScheduled(); v[2] = it.Count(); v[3] = it.Interval();
		v[4] = it.LastWritten();
		Set("out." + items[k], v);
	}
}

void C_RestartState::TakeClock(C_TimeSteps &ts, const std::vector<std::string> &items, bool restoreDt,
                               const std::string &fieldItem, int fieldCount) const
{
	const double now = Get("time");
	const double dtNow = restoreDt ? Get("dt", 0) : ts.DT().now;
	//  The step the last step was solved with; a file that does not hold
	//  it takes the step about to be taken in its place.
	ts.Restore(now, (std::uint64_t)Get("step"),
	           dtNow,
	           restoreDt ? Get("dt", 1) : ts.DT().now,
	           restoreDt ? Get("dt", 2, dtNow) : ts.DT().now);
	for(size_t k = 0; k < items.size(); k++)
	{
		C_OutputIntervals &it = ts.OutputItem(items[k]);
		const std::string key = "out." + items[k];
		const unsigned int count = (items[k] == fieldItem) ? (unsigned int)fieldCount
		                                                   : (unsigned int)Get(key, 2);
		const double last = Get(key, 4, Get(key, 0));
		if(Has(key) && Get(key, 3) == it.Interval())
			it.Restore(Get(key, 0), Get(key, 1), count, last);
		else if(it.IsList())
			it.RestoreListAfter(now, count, last);
		else
			it.Restore(now, now + it.Interval(), count, last);
	}
}

bool C_RestartState::Load(const std::string &fieldFile, const std::string &model, double tRestart, std::string &message)
{
	const std::string path = PathFor(fieldFile);
	_found = false;
	if(!Read(path))
	{
		//  A file that is there and cannot be read (cut short by a full
		//  disk, say) is reported as damaged, and with a warning, rather
		//  than as missing
		const bool there = std::ifstream(path.c_str()).good();
		_found = there;
		message = (there ? path + " is damaged and is not used" : "no state file " + path)
		          + ": the run starts from the field of " + fieldFile
		          + " alone. The time step, the counters of its control, the running totals of"
		          " History.dat and the steady-state count start afresh, so a run continued this"
		          " way is not, to the bit, the run it would have been without the interruption.";
		return false;
	}
	_found = true;
	char buf[512];
	const double t = Get("time"), dt = Get("dt");
	const double tol = 1e-5 * std::max(std::fabs(t), std::fabs(dt));
	if(Text("model") != model)
	{
		message = path + " was written by another model: only the field is taken.";
		Clear();
		return false;
	}
	if(!(std::fabs(tRestart - t) <= tol))
	{
		std::snprintf(buf, sizeof(buf), "%s was written at t = %.9g (dimensionless), RestartTime"
		              " gives %.9g: only the field is taken, as a fresh start from that state.",
		              path.c_str(), t, tRestart);
		message = buf;
		Clear();
		return false;
	}
	//  A state file that was damaged (edited, or written by something else)
	//  is not continued from: its counters are cast to integers by the
	//  models, where a NaN or a value out of range is undefined behavior,
	//  and its values seed the run.
	std::string bad;
	if(!Sane(bad))
	{
		message = path + " is damaged (" + bad + "): only the field is taken, as a fresh start"
		          " from that state.";
		Clear();
		return false;
	}
	message = "continuing exactly from " + path;
	return true;
}

bool C_RestartState::Sane(std::string &bad) const
{
	char buf[256];
	for(std::map<std::string, std::vector<double> >::const_iterator it = _scalars.begin(); it != _scalars.end(); ++it)
		for(size_t k = 0; k < it->second.size(); k++)
			if(!std::isfinite(it->second[k]))
			{
				std::snprintf(buf, sizeof(buf), "%.60s is not a finite number", it->first.c_str());
				bad = buf;
				return false;
			}
	for(std::map<std::string, std::vector<double> >::const_iterator it = _vectors.begin(); it != _vectors.end(); ++it)
		for(size_t k = 0; k < it->second.size(); k++)
			if(!std::isfinite(it->second[k]))
			{
				std::snprintf(buf, sizeof(buf), "%.60s holds a value that is not a finite number", it->first.c_str());
				bad = buf;
				return false;
			}
	//  What is cast to an integer: a whole number that fits, not negative.
	struct Count { const char *key; size_t k; double max; };
	std::vector<Count> counts;
	counts.push_back(Count{"step", 0, 9.0e18});
	counts.push_back(Count{"steadyCount", 0, 2147483647.0});
	counts.push_back(Count{"adjustCount.height", 0, 2147483647.0});
	for(std::map<std::string, std::vector<double> >::const_iterator it = _scalars.begin(); it != _scalars.end(); ++it)
		if(it->first.compare(0, 4, "out.") == 0) counts.push_back(Count{it->first.c_str(), 2, 4294967295.0});
	for(size_t c = 0; c < counts.size(); c++)
	{
		const double v = Get(counts[c].key, counts[c].k, 0);
		if(!(v >= 0 && v <= counts[c].max && v == std::floor(v)))
		{
			std::snprintf(buf, sizeof(buf), "%.60s = %g is not a count", counts[c].key, v);
			bad = buf;
			return false;
		}
	}
	//  The time steps the run would go on with
	if(Has("dt") && !(Get("dt", 0) > 0 && Get("dt", 1, 1) > 0))
	{
		bad = "dt is not positive";
		return false;
	}
	return true;
}
