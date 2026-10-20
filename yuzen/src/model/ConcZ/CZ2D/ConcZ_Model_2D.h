/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ConcZ_Model_2D.h : driver of the ConcZ model (a 1-D or 2-D film plus a concentration profile across the thickness)
*
*  ConcZ_Model.h
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _CONCZ_MODEL_2D_INCLUDED_
#define _CONCZ_MODEL_2D_INCLUDED_

#include "YuzenUtil.h"
using namespace yuzen::util;
#include <memory>
#include <vector>
#include "CZ_HeightEquation_2D.h"
#include "CZ_ConcEquation_2D.h"

#include "SolverModel.h"

#include "Output.h"
#include "FieldInfo.h"
#include "Probe.h"
#include "RestartState.h"
#include "CourantCheck.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif


/////////////////////////////////////////////////////////////////////////////
class ConcZ_Model_2D : public ISolverModel
{
public:
	ConcZ_Model_2D();
	virtual ~ConcZ_Model_2D(){}

	/**
	 * @brief Name of this model as it appears in the log.
	 *
	 * The film is 1-D ("ConcZ-2D": x and the thickness direction) or 2-D
	 * ("ConcZ-3D": x, y and the thickness direction), decided by the grid.
	 */
	virtual const char* Name() const { return _name; }
	/** @brief True when the run was stopped by the divergence check. */
	virtual bool Diverged() const { return _diverged; }
	/**
	 * @brief Build the whole calculation from the input.
	 *
	 * Grid, field, properties, non-dimensional numbers, both equations,
	 * boundary condition, probes and output are set up here, and the initial
	 * condition is applied. Both film dimensions are supported: Ny is used as
	 * it stands when Grid.Dimension is "xyh", and with "xh" the film is
	 * one-dimensional, so an Ny other than 1 is reset to 1 with a warning.
	 */
	virtual void Build(Input& inp);
	/**
	 * @brief Run the time-step loop to the end.
	 *
	 * Stops when the end time is reached, when the film has dried out or when
	 * the solution has gone steady.
	 */
	virtual void Run();

private:
	Input* _inputParam;
	Input::I_Grid _gridNonDim;

	Field_CZ _field;

	C_Fluid _fluid;
	C_TimeSteps _timeStep;
	NonDimensionalNumbers _NonDim;


	BoundaryCondition2D _boundary;

	bool _swConcField;

	CZ_HeightEquation_2D _heightEq;
	CZ_ConcEquation_2D _concEq;
	//  The Courant number of a limited concentration scheme (CourantCheck.h)
	C_CourantCheck _courantConc;

	int _iter_heightEq;
	double _resid_heightEq;

	//  Counters for performance measurement. They play no part in the numerics.
	//  The Elapsed Time of JobInfo includes building the grid, the initial field
	//  and the initial output, so the time-step loop is timed separately here.
	//  tools/scaling/ reads these.
	long long _nSteps;			// number of time steps carried out
	long long _sumNKIter;		// total number of Newton iterations
	long long _sumBCGIter;		// total number of linear iterations
	long long _sumCoupling;		// total number of coupling passes taken again
	//  The solute moved into cells below c = 0 over the run, the steps
	//  that moved any, and the solute made up there
	//  (CZ_HeightEquation_2D::RedistributeDeficit); its warning is given once
	double    _deficitMoved, _deficitAdded;
	long long _deficitSteps;
	bool      _warnedAdd;
	double    _loopSeconds;		// wall-clock time of the time-step loop [s]
	double    _sumStepDt;		// sum of the steps solved in the loop (dimensionless)
	double    _loopStartTime;	// time the loop started from (dimensionless)

	OutputVTK _outputVTK;
	OutputNeutral _outputNeutral;
	OutputTimeline _outputHistory, _outputCalcInfo;
	// One per probe. OutputTimeline cannot be copied, so they are held in a
	// vector of unique_ptr, which also releases them.
	std::vector<std::unique_ptr<OutputTimeline> > _probeOut;
	C_PVD_File _pvdFile;

	ProbeList _probeList;
	FieldInfo _fieldInfo, _fieldInfoIni;

	FilePath itemField;
	FilePath itemProbe;
	FilePath itemHistory;
	FilePath itemCalc;

	//  Find, on the 3-D grid, the cell nearest each probe position given in the
	//  input, and create the output destination for it.
	void BuildProbes();
	//  Write out whichever items (field, probes, timeline, calculation info) have
	//  reached their output time.
	/**
	 * @brief Write whichever of the periodic outputs are due.
	 * @param atEnd  true for the single call made after the time-step loop,
	 *               where every output not already holding the present time
	 *               is due whatever the schedule says.
	 */
	void FreqOutput(bool atEnd = false);

	/** @brief Whether one kind of output is due now. */
	bool OutputDue(const std::string& key, bool atEnd)
	{ return atEnd ? !_timeStep.WrittenAt(key) : _timeStep.TimeForOut(key); }
	/** @brief Record that one kind of output was written. */
	void MarkWritten(const std::string& key, bool atEnd)
	{ if(atEnd) _timeStep.UpdateItemAtEnd(key); else _timeStep.UpdateItem(key); }

	//  Exact continuation (field*.state, see RestartState.h)
	void WriteState(const std::string &path);
	//  Replace the field read from the field file by the exact copy in the
	//  state file, if the two are the same state; true when it was.
	bool RefineFromState();
	std::vector<std::string> OutputItems() const
	{ std::vector<std::string> v; v.push_back(itemField.c_str()); v.push_back(itemProbe.c_str());
	  v.push_back(itemHistory.c_str()); v.push_back(itemCalc.c_str()); return v; }
	C_RestartState _rstate;
	bool _continuing;

	//  Counter for the steady-state test. It counts how many steps in a row have
	//  been quiet, so it has to be kept across steps.
	int _steadyCount;

	//  The balance of History.dat: the volume and the solute at the start
	//  of the whole run, and the running totals of what evaporated and of
	//  what the 1-D centrifugal term flung sideways (liquid and solute)
	double _vol0, _resin0;
	double _evaporatedVolume, _outflowedVolume, _outflowedResin;
	ScalarField3D<> _cStart;	// c at the start of the step (1-D with the centrifugal term)

	//  The name written to the log. It follows from the grid's ny() (set in Build()).
	const char *_name;

	//------------------------------------------------------------------
	//  Detection of divergence.
	//  _stepStatus is the return value of the last TimeAdvance (0 means success),
	//  so that a step that has not converged is not simply carried on from.
	//  Once _diverged is set, Run() breaks out of its loop and the exit code
	//  reports it.
	//------------------------------------------------------------------
	int  _stepStatus;
	bool _diverged;

	//  How far the concentration and the film height may stray outside their
	//  physical range before the run is taken to have diverged. A healthy
	//  calculation can stray slightly through discretization error, so the
	//  allowance is set loosely, to catch only excursions of an order of
	//  magnitude.
	static const double kConcSlack;
	static const double kHeightSlack;
	//  Termination test. It looks at three things in turn: the time, drying out,
	//  and the solution going steady.
	bool CheckTermination();
	/** @brief Whether the field _fieldInfo holds has met the dry-out test. */
	bool DriedOut() const;
	//  If restart is enabled, read the field, the time and the file number back in.
	void FillGhostColumns();
	void BuildRestart();

	ConcZ_Model_2D(const ConcZ_Model_2D&);
	void operator=(const ConcZ_Model_2D&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_CONCZ_MODEL_2D_INCLUDED_
