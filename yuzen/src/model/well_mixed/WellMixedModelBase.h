/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  WellMixedModelBase.h : the dimension-independent part of the model driver
*
*  The 1-D and 2-D well-mixed models (WellMixed_Model_1D.cpp, _2D.cpp)
*  share everything here. What differs between them is
*
*    - the types (grid, site, field, boundary condition, equations), which
*      the traits of WellMixedTraits.h supply, and
*    - the small differences that follow from whether there is a Y
*      direction (the arguments of the boundary condition, the smallest
*      mesh spacing, the volume of a cell, locating the probes and their
*      column headings, the informational output), left to the derived
*      classes as virtual functions.
*
*  The dry-out test follows the input (EndConcType / EndConcValue) in both,
*  and the time-step loop accumulates its sums with the step the equations
*  were solved with (StepDt).
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _WELL_MIXED_MODEL_BASE_INCLUDED_
#define _WELL_MIXED_MODEL_BASE_INCLUDED_

#include "YuzenUtil.h"
#include "Profile.h"
#include <algorithm>
#include <iomanip>
#include <sstream>
#include "WellMixedTraits.h"
#include "HeightEquationBase.h"
#include "BoundaryCondition.h"
#include "FluidProperty.h"
#include "NonDimensionalNumbers.h"
#include "TimeSteps.h"
#include "Output.h"
#include "Probe.h"
#include "FieldInfo.h"
#include "InitialCondition.h"
#include "Restart.h"
#include "RestartState.h"
#include "CourantCheck.h"
#include "EnergyBudget.h"
#include "Input.h"

#include <cmath>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace yuzen::util;
using namespace yuzen::util::Global;

#include "SolverModel.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//--------------------------------------------------------------------------
//  The driver of the well-mixed model.
//
//  It looks after one computation from beginning to end.
//    Build()  assemble the grid, the fields, the fluid properties, the
//             non-dimensional numbers, the equations and the output, and
//             put in the initial condition (or the restart)
//    Run()    step through time, writing output and testing for termination
//
//  There are three ways a run ends: the specified time is reached, the film
//  dries out, or the solution goes steady. The steady test asks "how much
//  would the shape change over all the remaining time, at the rate it is
//  changing now" (see the comments in CheckTermination).
//--------------------------------------------------------------------------
template<class T>
class WellMixed_ModelBase : public ISolverModel
{
public:
	typedef typename T::Grid     Grid;
	typedef typename T::Site     Site;
	typedef typename T::Field    Field;
	typedef typename T::FieldOut FieldOut;
	typedef typename T::Boundary Boundary;
	typedef typename T::HeightEq HeightEq;
	typedef typename T::ConcEq   ConcEq;

	/** @brief Zero the accumulators; Build() constructs everything else. */
	WellMixed_ModelBase();
	virtual ~WellMixed_ModelBase(){}

	/**
	 * @brief Assemble a whole run from the input.
	 *
	 * Builds the fluid properties, the non-dimensional numbers, the grid and
	 * fields, the boundary condition, the time steps, the concentration and
	 * height equations, the probes and the output writers, and applies the
	 * initial condition.
	 */
	virtual void Build(Input& inp);
	/**
	 * @brief Run the time-step loop to the end.
	 *
	 * Reads a restart file if one was requested, then steps until the target
	 * time is reached, the film dries out, or the solution goes steady.
	 */
	virtual void Run();
	/** @brief True when the run was stopped by the divergence check. */
	virtual bool Diverged() const { return _diverged; }

protected:
	//--------------------------------------------------------------
	//  What differs with the dimension (implemented by the derived class)
	//--------------------------------------------------------------
	/**
	 * @brief Build the boundary condition.
	 *
	 * The order of the arguments differs with the dimension, which is why
	 * this lives here.
	 */
	virtual void   BuildBoundaryCondition(Input &inp, Grid *G) = 0;
	/**
	 * @brief Zero the initial velocity at one point. The number of
	 *        components differs with the dimension.
	 */
	virtual void   ZeroInitialVelocity(Grid *G, const Site &s) = 0;
	/** @brief The representative length used to cap the time step. */
	virtual double MinMeshSize(Grid *G) const = 0;
	/** @brief The volume of one cell (dx in 1-D, dx*dy in 2-D). */
	virtual double CellMeasure(Grid *G) const = 0;
	/** @brief Log the extent of the domain and the number of grid points. */
	virtual void   InfoGridExtent() const = 0;
	/**
	 * @brief Find the grid point nearest to each probe point given in the
	 *        input.
	 */
	virtual void   FindProbes(Grid *G) = 0;
	/**
	 * @brief Add the column headings of the velocity components of the probe
	 *        output (Ux only in 1-D, Ux and Uy in 2-D).
	 */
	virtual void   AddVelocityLabels(std::vector<std::string> &labels, bool dimensioned) const = 0;
	/**
	 * @brief Write the coordinates of the probe into the header of its
	 *        output file.
	 */
	virtual void   WriteProbeHeader(std::ostringstream &oss, ProbeList::iterator it) const = 0;
	/**
	 * @brief Read an ASCII restart file.
	 *
	 * Only this one differs between the dimensions, and only in the name of
	 * the function it calls.
	 */
	virtual bool   ReadRestartAscii(Input::I_Restart &inp, double coef) = 0;

	Input* _inputParam;
	Input::I_Grid _gridNonDim;

	Field _field, _field_guess;
	FieldOut _field_out;

	C_Fluid _fluid;
	C_TimeSteps _timeStep;
	NonDimensionalNumbers _NonDim;

	OutputVTK _outputVTK;
	OutputNeutral _outputNeutral;
	OutputTimeline _outputHistory, _outputCalcInfo;
	//  Output.Energy (EnergyBudget.h): Energy.dat, written with History.dat.
	//  F and D of the state the run holds; the dissipated energy int D dt
	//  (trapezoidal over the steps, as Crank-Nicolson weights them); the
	//  balance (F^{n+1} - F^n)/dt + (D^n + D^{n+1})/2 of the last step, and
	//  whether there was one. Of the step, not of the interval between two
	//  rows, so that a run stopped and continued, which writes a row where
	//  it stopped, gives every row the value the uninterrupted run gives it.
	OutputTimeline _outputEnergy;
	bool _energyOn = false;
	double _enF = 0, _enD = 0, _enDiss = 0, _enBal = 0;
	bool _enStep = false;
	PotentialParams EnergyParams() const;
	void EnergyEval();
	// One per probe. OutputTimeline cannot be copied, so they are held as
	// a vector of unique_ptr, which also releases them.
	std::vector<std::unique_ptr<OutputTimeline> > _probeOut;
	C_PVD_File _pvdFile;

	ProbeList _probeList;
	FieldInfo _fieldInfo, _fieldInfoIni;

	Boundary _boundary;

	double _evaporatedVolume;
	double _suppliedVolume;		// what the nozzle (Physics.Source) has added
	double _suppliedResin;		// and the solute that came with it
	double _outflowedVolume;
	double _outflowedResin;

	FilePath itemField;
	FilePath itemProbe;
	FilePath itemHistory;
	FilePath itemCalc;

	int _NK_iter;
	double _NK_resid;

	//  Counters used for performance measurement. They take no part in the
	//  numerics at all.
	//  They exist so that the time-step loop alone can be measured. The
	//  Elapsed Time in JobInfo covers building the grid, the initial field
	//  and the initial output as well, so on a larger grid the setup stops
	//  being negligible and the figure is useless as a measure of parallel
	//  scaling.
	long long _nSteps;			// the number of time steps performed
	long long _sumNKIter;		// the total number of Newton iterations
	long long _sumBCGIter;		// the total number of linear iterations
	double    _loopSeconds;		// the wall-clock time of the time-step loop [s]
	double    _sumStepDt;		// the sum of the steps solved in the loop (dimensionless)
	double    _loopStartTime;	// the time the loop started from (dimensionless)

	HeightEq _heightEq;
	//  The Courant numbers of the limited schemes (CourantCheck.h)
	C_CourantCheck _courantConc, _courantHeight;
	//  The solute moved from cells above c = 1 to cells with room over the
	//  run, the steps that moved any, and the solute that had to be cut;
	//  the same for the cells below c = 0, filled from cells holding
	//  solute, and what had to be made up there
	//  (HeightEquationBase::RedistributeExcess); the warning of a cut, or of
	//  solute made up, is given once
	double    _excessMoved, _excessTrimmed;
	long long _excessSteps;
	bool      _warnedTrim;
	double    _deficitMoved, _deficitAdded;
	long long _deficitSteps;
	bool      _warnedAdd;
	/** @brief Report the solute moved and cut at c = 1, and moved and made up at c = 0 (end of the run). */
	void ReportExcess();
	ConcEq   _concEq;
	bool _swConcField;

	/**
	 * @brief The evaporation at this instant, as a volume per unit time.
	 *
	 * It is the time-step side that accumulates it.
	 */
	double CalEvaporatedVolume();
	/**
	 * @brief Locate the probe points and prepare a time-series output for
	 *        each of them.
	 */
	void BuildProbes();
	/**
	 * @brief Write out only those items whose output time has arrived
	 *        (the field, the probes, the general time series and the
	 *        computation information).
	 */
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

	//------------------------------------------------------------------
	//  Exact continuation (field*.state, see RestartState.h)
	//------------------------------------------------------------------
	/** @brief Write the state a restart from the field just written needs. */
	void WriteState(const std::string &path);
	/**
	 * @brief Replace the field read from the field file by the exact copy in
	 *        the state file, if the two are the same state.
	 * @return true when h, s and c were all replaced
	 */
	bool RefineFromState();
	/** @brief The output items, in the order they are kept in a state file. */
	std::vector<std::string> OutputItems() const
	{ std::vector<std::string> v; v.push_back(itemField.c_str()); v.push_back(itemProbe.c_str());
	  v.push_back(itemHistory.c_str()); v.push_back(itemCalc.c_str()); return v; }
	C_RestartState _rstate;		// read by BuildRestart, applied there and in Run()
	bool _continuing;			// the restart continues from a state file
	//  The volumes the conservation errors of History.dat are taken against:
	//  those of the initial state, carried over a restart by the state file
	double _vol0, _resin0;

	//  The counter of the steady test. It counts how many steps in a row
	//  have been quiet, so it has to be kept across steps.
	int _steadyCount;
	//  The run was stopped because the field broke down (see CheckTermination)
	bool _diverged;
	//  How far outside [0, 1] the concentration, and below zero the
	//  dimensionless film thickness, may stray before the run is taken to
	//  have diverged (the values of ConcZ)
	static constexpr double kConcSlack   = 0.1;
	static constexpr double kHeightSlack = 1.0e-3;
	/**
	 * @brief Test for termination: the specified time is reached, the film
	 *        has dried out, or the solution has gone steady.
	 */
	bool CheckTermination();
	/** @brief Whether the field _fieldInfo holds has met the dry-out test. */
	bool DriedOut() const;
	/**
	 * @brief If a restart was requested, read the field, the time and the
	 *        file number back in.
	 */
	void BuildRestart();

	/**
	 * @brief Check that the film reaches zero at every Pinning edge.
	 *
	 * Pinning holds the film thickness at zero on its edge -- the condition
	 * is h = 0 at the face between the ghost layer and the first real cell --
	 * so the film it starts from has to go to zero there as well. If it does
	 * not, the first step meets a jump of the whole film thickness inside one
	 * cell; Newton stagnates, the step is cut until it underflows MinDt, and
	 * the run stops with nothing in its messages that points at the wall.
	 *
	 * The film is extrapolated linearly from the first two real cells to that
	 * face. More than a tenth of the thickest film there stops the run, with
	 * the numbers; otherwise the value is reported in the run log, so it can be
	 * seen on every run that pins an edge.
	 */
	void CheckPinnedEdges();

private:
	WellMixed_ModelBase(const WellMixed_ModelBase&);
	void operator=(const WellMixed_ModelBase&);
};


//==================================================================
//  Implementation
//==================================================================

template<class T>
WellMixed_ModelBase<T>::WellMixed_ModelBase()
	: _inputParam(0), _evaporatedVolume(0), _suppliedVolume(0), _suppliedResin(0), _outflowedVolume(0), _outflowedResin(0),
	  _NK_iter(0), _NK_resid(0),
	  _nSteps(0), _sumNKIter(0), _sumBCGIter(0), _loopSeconds(0), _sumStepDt(0), _loopStartTime(0),
	  _excessMoved(0), _excessTrimmed(0), _excessSteps(0), _warnedTrim(false),
	  _deficitMoved(0), _deficitAdded(0), _deficitSteps(0), _warnedAdd(false),
	  _swConcField(false), _continuing(false), _vol0(0), _resin0(0), _steadyCount(0),
	  _diverged(false)
{}

template<class T>
void WellMixed_ModelBase<T>::Build(Input& inp)
{
	//----------------------------------------------------
	// Copy the input parameters
	//----------------------------------------------------
	_inputParam = &inp;
	_swConcField = (inp.Physics.ConcField.Model.compare("WM") == 0);

	//----------------------------------------------------
	// Build the fluid properties
	//----------------------------------------------------
	_fluid.Build(_inputParam->Fluid);

	//----------------------------------------------------
	// Build the non-dimensional numbers
	//----------------------------------------------------
	_NonDim.Build(inp, inp.Initial.SurfaceShape.TotalHeight);

	//----------------------------------------------------
	// Non-dimensionalize the grid
	//----------------------------------------------------
	_gridNonDim = _inputParam->Grid;
	_gridNonDim.Xmin /= _NonDim.ExtentScale();
	_gridNonDim.Xmax /= _NonDim.ExtentScale();
	_gridNonDim.Ymin /= _NonDim.ExtentScale();
	_gridNonDim.Ymax /= _NonDim.ExtentScale();

	//----------------------------------------------------
	// Build the field (this also sets the substrate shape internally)
	//----------------------------------------------------
	_field.Build(_gridNonDim, inp.Substrate, _NonDim.ExtentScale());
	Grid *G = _field.GetGridPtr();

	_field_out.Build(_field, _NonDim.ExtentScale(), _NonDim.HeightScale(), _NonDim.VelocityScale(), inp.Output.Dimensioned);

	//----------------------------------------------------
	// Build the boundary condition
	//----------------------------------------------------
	BuildBoundaryCondition(inp, G);


	//----------------------------------------------------
	// Initialize the surface shape
	//----------------------------------------------------
	C_InitialSurface::Set(_field, inp.Initial.SurfaceShape, _inputParam->Substrate, _NonDim.ExtentScale());

	//----------------------------------------------------
	// Scaling and initial settings
	//----------------------------------------------------
	for(Site st = G->template First<Site>(); st != G->template Last<Site>(); ++st)
	{
		_field.HH(st) /= _inputParam->Initial.SurfaceShape.TotalHeight;
		_field.SS(st) /= _inputParam->Initial.SurfaceShape.TotalHeight;
		_field.CC(st) = _inputParam->Initial.Concentration.Value;
		ZeroInitialVelocity(G, st);
	}

	//  With a manufactured solution the film, the substrate and the
	//  concentration start from it, ghosts included (verification only)
	{
		C_ManufacturedSource mms;
		mms.Build(inp, _NonDim, T::nAxis == 2);
		if(mms.Enabled())
		{
			for(Site st = G->template First<Site>(); st != G->template Last<Site>(); ++st)
			{
				double xx, yy;
				T::SiteXY(G, st, xx, yy);
				_field.HH(st) = mms.Hnd(xx, yy, 0.0);
				_field.SS(st) = mms.Subnd(xx, yy);
				_field.CC(st) = mms.Cnd(xx, yy, 0.0);
			}
		}
	}

	//----------------------------------------------------
	//  Make the initial shape and the substrate shape consistent with the
	//  boundary condition
	//----------------------------------------------------
	_boundary.ApplyToField(_field);

	//----------------------------------------------------
	//  Build the time steps
	//----------------------------------------------------
	double deltaT = _inputParam->Numerics.InitialDt;
	double startTime = 0;
	double endTime = inp.Numerics.EndTime;

	double T0 = _NonDim.TimeScale();
	deltaT		/= T0;
	startTime	/= T0;
	endTime		/= T0;

	_timeStep.Build(deltaT, startTime, endTime, inp.Numerics.TimeStepControl.Enable);

	//  Numerics.MaxSteps: stop after the given number of steps. The default
	//  (two billion) is effectively unlimited, so by default the run ends on
	//  EndTime alone.
	if(inp.Numerics.MaxSteps > 0)
		_timeStep.SetMaxSteps((std::uint64_t)inp.Numerics.MaxSteps);

	//----------------------------------------------------
	// Set the largest and smallest time step
	//----------------------------------------------------
	double max_dt, min_dt;
	double factor = inp.Numerics.TimeStepControl.MaxDtFactor;
	double dxMin = MinMeshSize(G);
	max_dt = factor * pow(dxMin, 4.0);
	min_dt = inp.Numerics.TimeStepControl.MinDt;
	_timeStep.SetMaxDT(max_dt);
	_timeStep.SetMinDT(min_dt);

	_field_guess = _field;

	//----------------------------------------------------
	// Build the equation of the concentration field
	//----------------------------------------------------
	if(_swConcField)
	{
		_concEq.Build(inp, &_fluid, &_timeStep, &_field_guess, &_NonDim, &_boundary);
		_concEq.SetX(_field_guess.CC());
		_concEq.Initialize();
	}

	_field = _field_guess;

	//----------------------------------------------------
	// Build the equation of the height function
	//----------------------------------------------------
	_heightEq.Build(inp, &_fluid, &_timeStep, &_concEq, &_field, &_field_guess, &_NonDim, &_boundary);
	_heightEq.SetX(_field.HH());
	_heightEq.Initialize();


	//----------------------------------------------------
	// Build the probes
	//----------------------------------------------------
	BuildProbes();

	//----------------------------------------------------
	// Set up the output data items
	//----------------------------------------------------
	itemField	= "Field";
	itemProbe	= "Probe";
	itemHistory	= "GeneralHistory";
	itemCalc	= "CalcInfo";

	if(_inputParam->Output.FieldData.Listed())
	{
		//  Output.FieldData.Times or LogSchedule (dimensionless here)
		std::vector<double> listed(_inputParam->Output.FieldData.Times);
		for(size_t k = 0; k < listed.size(); k++) listed[k] /= T0;
		_timeStep.Add_OutputList(itemField.c_str(), listed, _inputParam->Output.FieldData.ExactTimes);
	}
	else
	_timeStep.Add_OutputItem(itemField.c_str(), _inputParam->Output.FieldData.TimeInterval / T0);
	_timeStep.Add_OutputItem(itemProbe.c_str(),   _inputParam->Output.ProbeList.TimeInterval / T0);
	_timeStep.Add_OutputItem(itemHistory.c_str(), _inputParam->Output.GeneralHistory.TimeInterval / T0);
	_timeStep.Add_OutputItem(itemCalc.c_str(),    _inputParam->Output.CalcInfo.TimeInterval / T0);

	//----------------------------------------------------
	// Build the output classes
	//----------------------------------------------------
	std::vector<std::string> histLabels, calcLabels;
	//  History.dat carries the raw scales: volumes in m^3 and heights in m,
	//  not the mm/um of the field files. Say so in the header, since that is
	//  exactly the difference somebody comparing two files will trip over.
	const bool dim = _inputParam->Output.Dimensioned;
	//  In 1-D, per meter of the width of the strip
	const char *uV = dim ? (_inputParam->Grid.Dimension == "xh" ? " [m^3/m]" : " [m^3]") : " [-]";
	const char *uH = dim ? " [m]"   : " [-]";
	histLabels.push_back(std::string("Volume") + uV);
	histLabels.push_back(std::string("ResinVol") + uV);
	histLabels.push_back("AveC [-]");
	histLabels.push_back(std::string("EvaporatedVol") + uV);
	histLabels.push_back(std::string("OutflowedVol") + uV);
	histLabels.push_back(std::string("OutflowedResin") + uV);
	histLabels.push_back("VolErr [-]");	histLabels.push_back("ResinVolErr [-]");
	histLabels.push_back(std::string("maxHS") + uH);
	histLabels.push_back(std::string("minHS") + uH);
	histLabels.push_back(std::string("maxH") + uH);
	histLabels.push_back(std::string("minH") + uH);
	histLabels.push_back("maxC [-]");	histLabels.push_back("minC [-]");
	//  Last, so that the columns before it keep their places
	histLabels.push_back(std::string("SuppliedVol") + uV);
	histLabels.push_back(std::string("SuppliedResin") + uV);

	calcLabels.push_back("NK_iter");	calcLabels.push_back("NK_resid");
	calcLabels.push_back("BCG_iter");	calcLabels.push_back("BCG_resid");
	calcLabels.push_back("deltaT");
	//  Write out the very quantity the steady test uses. It gives a guide
	//  when choosing Numerics.SteadyStateCriterion.
	calcLabels.push_back("ChangeRate");	// the relative rate of change of the film height [1/s]

	_outputVTK.Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), "field", 0);
	_outputNeutral.Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), "field", 0);
	_outputHistory.Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), "History", (int)histLabels.size(), histLabels);
	_outputCalcInfo.Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), "CalcInfo", (int)calcLabels.size(), calcLabels);
	_energyOn = _inputParam->Output.Energy.Enable;
	if(_energyOn)
	{
		//  The energy on the scale sigma0 h0^2 (per meter of width in 1-D),
		//  the power on that over T0
		const bool oneD = _inputParam->Grid.Dimension == "xh";
		const char *uE = dim ? (oneD ? " [J/m]" : " [J]") : " [-]";
		const char *uP = dim ? (oneD ? " [W/m]" : " [W]") : " [-]";
		std::vector<std::string> enLabels;
		enLabels.push_back(std::string("FreeEnergy") + uE);
		enLabels.push_back(std::string("Dissipation") + uP);
		enLabels.push_back(std::string("Dissipated") + uE);
		enLabels.push_back(std::string("Balance") + uP);
		_outputEnergy.Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), "Energy", (int)enLabels.size(), enLabels);
	}

	double st_time = _timeStep.StartTime();
	double freq = _timeStep.OutputItem(itemField.c_str()).Interval();
	int fst = 0;
	if(_inputParam->Output.Dimensioned)
	{
		st_time *= _NonDim.TimeScale();
		freq    *= _NonDim.TimeScale();
	}
	_pvdFile.Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), _inputParam->JobName, "field", st_time, freq, fst);


	//----------------------------------------------------
	//  Record in JobInfo.txt how the system was solved.
	//
	//  Without this there is no way to tell from the results what
	//  Backend = "auto" resolved to, or whether "amgcl" fell back to eigen
	//  in a build without AMGCL. A regression comparison has to be made
	//  with the backends matched, so this is written down.
	//----------------------------------------------------
	JobInfo<<"\n###### Solver Info ######\n";
	JobInfo.Add("Model",           std::string(Name()),        "Model");
	JobInfo.Add("BackendHeight",   _heightEq.BackendName(),    "Linear solver <height>");
	if(_swConcField)
		JobInfo.Add("BackendConc", _concEq.BackendName(),      "Linear solver <conc>");
	//  The face interpolation of the convective terms (Numerics.Convection)
	JobInfo.Add("ConvectionHeight", _heightEq.ConvectionKey(), "Convection scheme <height>");
	if(_swConcField)
		JobInfo.Add("ConvectionConc", C_Interpolation_base4::Key(_concEq.Scheme().Scheme()), "Convection scheme <conc>");

	//----------------------------------------------------
	// Informational output
	//----------------------------------------------------
	Info()<<std::scientific;
	Info<<"###### Solver #####"<<nl;
	Info<<" Linear solver <height> = "<<_heightEq.BackendName()<<nl;
	if(_swConcField)
		Info<<" Linear solver <conc>   = "<<_concEq.BackendName()<<nl;
	Info<<" Convection <height>    = "<<_heightEq.ConvectionKey()<<nl;
	if(_swConcField)
		Info<<" Convection <conc>      = "<<C_Interpolation_base4::Key(_concEq.Scheme().Scheme())<<nl;
	Info<<"######## Domain & Grids ########"<<nl;
	InfoGridExtent();
	Info<<_NonDim.ToString()<<nl;
	Info<<"###### Non-dimensional time info #####"<<nl;
	Info<<" deltaT    = "<<deltaT <<nl;
	Info<<" startTime = "<<startTime <<nl;
	Info<<" endTime   = "<<endTime <<nl;
	Info<<" max dt    = "<<_timeStep.MaxDT()<<nl;
	Info<<" min dt    = "<<_timeStep.MinDT()<<nl;
	Info<<"--- output intervals ----"<<nl;
	if(_inputParam->Output.FieldData.Listed())
		Info<<" FieldData = "<<_inputParam->Output.FieldData.Times.size()<<" listed times, "
		    <<_inputParam->Output.FieldData.Times.front() / T0<<" to "<<_inputParam->Output.FieldData.Times.back() / T0
		    <<(_inputParam->Output.FieldData.ExactTimes ? " (exact)" : "")<<nl;
	else
		Info<<" FieldData = "<<_inputParam->Output.FieldData.TimeInterval / T0 <<nl;
	Info<<" Probe     = "<<_inputParam->Output.ProbeList.TimeInterval / T0 <<nl;
	Info<<" History   = "<<_inputParam->Output.GeneralHistory.TimeInterval / T0 <<nl;
	Info<<" CalcInfo  = "<<_inputParam->Output.CalcInfo.TimeInterval / T0 <<nl;
	Info<<"######################################"<<nl;


}

template<class T>
void WellMixed_ModelBase<T>::BuildProbes()
{
	//-------------------------
	//  Locate the probe points
	//-------------------------
	Grid *G = _field.GetGridPtr();
	FindProbes(G);



	std::vector<std::string> probeLabels;
	ostringstream oss;

	if(_inputParam->Output.Dimensioned)
	{
		probeLabels.push_back("HS [um]");	probeLabels.push_back("Conc [-]");
		AddVelocityLabels(probeLabels, true);
		probeLabels.push_back("HH [um]");	probeLabels.push_back("SS [um]");
	}
	else
	{
		probeLabels.push_back("HS");	probeLabels.push_back("Conc");
		AddVelocityLabels(probeLabels, false);
		probeLabels.push_back("HH");	probeLabels.push_back("SS");
	}
	char head[2000];
	unsigned int ii = 0;

	_probeOut.clear();
	for(size_t ii = 0; ii < _probeList.size(); ii++)
		_probeOut.push_back(std::unique_ptr<OutputTimeline>(new OutputTimeline()));
	for(ProbeList::iterator it = _probeList.begin(); it != _probeList.end(); ++it)
	{
		sprintf_s(head, "%s%02u", "Probe", ii);
		oss.clear();
		oss.str("");

		WriteProbeHeader(oss, it);

		_probeOut[ii]->Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), head, (int)probeLabels.size(), probeLabels, oss.str());
		ii++;
	}
}

template<class T>
void WellMixed_ModelBase<T>::FreqOutput(bool atEnd)
{
	double volScale = 1;
	double heightScale = 1;
	double timeScale   = 1;
	double velocityScale = 1;
	if(_inputParam->Output.Dimensioned)
	{
		timeScale = _NonDim.TimeScale();
		//  A 1-D film is a strip 1 m wide in y: its volumes are per meter
		//  of width (the sum of h dx), whatever the length scale or
		//  Ymin/Ymax
		volScale = _NonDim.ExtentScale() * _NonDim.HeightScale();
		if(_inputParam->Grid.Dimension != "xh") volScale *= _NonDim.ExtentScale();
		heightScale =  _NonDim.HeightScale();
		velocityScale =  _NonDim.VelocityScale();
	}
	//  Probes follow the field files and report thicknesses in um, which is
	//  what their header has always claimed. History.dat keeps the raw meter
	//  scale (heightScale above), because its volumes are in m^3 and mixing
	//  the two in one file would be worse than the difference between files.
	const double probeHeightScale = heightScale * (_inputParam->Output.Dimensioned ? 1e6 : 1.0);

	char wfile_head[2000];

	// Write the field file
	std::string statePath;		// the state file to write once every output is done
	if(OutputDue(itemField.c_str(), atEnd))
	{
		statePath = C_RestartState::PathFor(_inputParam->Output.FieldData.NeutralFile
		                                    ? _outputNeutral.FileName() : _outputVTK.FileName());
		Info<<"Field file " << _timeStep.OutputItem(itemField.c_str()).Count()  << ", Time = "<<_timeStep.NowTime() << nl;

		_field_out.ConvertFromStg(_field);
		if(_inputParam->Output.FieldData.WritePotential)
		{
			//  The driving potential (LubFields.h, PotentialParams); dimensioned,
			//  on the pressure scale sigma0 h0 / L0^2 [Pa]
			PotentialParams pp;
			pp.capillary = _inputParam->Physics.LaplacePressure.Enable;
			pp.dp    = _NonDim.Dp();
			pp.hstar = _NonDim.HStar();
			pp.n     = _inputParam->Physics.DisjoiningPressure.ExponentN;
			pp.m     = _inputParam->Physics.DisjoiningPressure.ExponentM;
			pp.bo    = _NonDim.Bo();
			pp.nd    = &_NonDim;
			pp.scale = _inputParam->Fluid.SurfaceTension.BaseValue * _NonDim.HeightScale()
			         / (_NonDim.ExtentScale() * _NonDim.ExtentScale());
			_field_out.ConvertPotential(_field, pp);
		}

		// The ASCII file
		//  A file that cannot be written stops the run here, before the
		//  state file that would go with it (StopOnOutputError).
		if(_inputParam->Output.FieldData.NeutralFile)
			if(!_outputNeutral.Write(_field_out)) StopOnOutputError();

		//The VTK file
		if(_inputParam->Output.FieldData.VTKFile )
		{
			if(!_outputVTK.Write(_field_out)) StopOnOutputError();
			if(!_pvdFile.Add(_timeStep.NowTime() * timeScale)) StopOnOutputError();
		}
		MarkWritten(itemField.c_str(), atEnd);
	}


	//Write the time series at the probe points
	if(OutputDue(itemProbe.c_str(), atEnd))
	{
		ostringstream oss;
		int p = 0;
		for(std::list<PointProbe>::iterator it = _probeList.begin(); it != _probeList.end(); ++it)
		{
			_field.Accept(*it);
			oss.clear();
			oss.str("");
			oss<<std::scientific<<_timeStep.NowTime() * timeScale;
			oss<<"\t"<<it->HS() * probeHeightScale <<"\t"<<it->CC()
				<<"\t"<<it->Ux() * velocityScale
				<<"\t"<<it->Uy() * velocityScale
				<<"\t"<<it->HH() * probeHeightScale <<"\t"<<it->SS() * probeHeightScale;
			if(!_probeOut[p]->Write(oss)) StopOnOutputError();
			p++;
		}
		MarkWritten(itemProbe.c_str(), atEnd);
	}

	//Write the general time series
	if(OutputDue(itemHistory.c_str(), atEnd))
	{
		ostringstream oss;
		double vf, vf0, vr, vr0, evaporated, outVol, outResin, volErr, resinErr;
		_field.Accept(_fieldInfo);
		vf = _fieldInfo.FluidVol();	vf0 = _fieldInfoIni.FluidVol();
		vr = _fieldInfo.ResinVol(); vr0 = _fieldInfoIni.ResinVol();
		evaporated = _evaporatedVolume;
		outVol     = _outflowedVolume;
		outResin   = _outflowedResin;

		vf0 = _vol0; vr0 = _resin0;
		//  What the nozzle supplied, and the solute that came with it, are
		//  taken off
		volErr   = ( vf + outVol + evaporated - _suppliedVolume)/vf0 -1.0;
		resinErr = ( vr + outResin - _suppliedResin)/vr0 -1.0;

		oss<<std::scientific<<_timeStep.NowTime() * timeScale;
		oss<<"\t"<< vf * volScale<< "\t" << vr * volScale << "\t"<< _fieldInfo.AveConc();
		oss<<"\t"<< evaporated  * volScale<< "\t" << outVol * volScale << "\t" << outResin * volScale;
		oss<<"\t"<< volErr << "\t" << resinErr;
		oss<<"\t"<<_fieldInfo.MaxHS() * heightScale<<"\t"<<_fieldInfo.MinHS() * heightScale;
		oss<<"\t"<<_fieldInfo.MaxH() * heightScale<<"\t"<<_fieldInfo.MinH() * heightScale;
		oss<<"\t"<<_fieldInfo.MaxC()<<"\t"<<_fieldInfo.MinC();
		oss<<"\t"<< _suppliedVolume * volScale << "\t" << _suppliedResin * volScale;
		if(!_outputHistory.Write(oss)) StopOnOutputError();
		if(_energyOn)
		{
			//  The scales of the energy and the power (Build); the balance
			//  is that of the interval since the last row, 0 on the first
			double eS = 1, pS = 1;
			if(_inputParam->Output.Dimensioned)
			{
				eS = _inputParam->Fluid.SurfaceTension.BaseValue * _NonDim.HeightScale() * _NonDim.HeightScale();
				if(_inputParam->Grid.Dimension == "xh") eS /= _NonDim.ExtentScale();
				pS = eS / _NonDim.TimeScale();
			}
			ostringstream oe;
			oe<<std::scientific<<std::setprecision(17)<<_timeStep.NowTime() * timeScale;
			oe<<"\t"<<_enF * eS<<"\t"<<_enD * pS<<"\t"<<_enDiss * eS<<"\t"<<(_enStep ? _enBal : 0.0) * pS;
			if(!_outputEnergy.Write(oe)) StopOnOutputError();
		}
		MarkWritten(itemHistory.c_str(), atEnd);
	}

	//Write the computation information
	if(OutputDue(itemCalc.c_str(), atEnd))
	{
		ostringstream oss;
		_field.Accept(_fieldInfo);
		sprintf_s(wfile_head, "%s%s", JobInfo.ResultDir().c_str(), "/CalcInfo.dat");

		//  The time and the time step are written on the same scale as in
		//  History.dat ([s]).
		oss<<std::scientific<<_timeStep.NowTime() * timeScale;
		oss<<"\t"<<_NK_iter << "\t" <<_NK_resid
			<<"\t"<< _heightEq.BCG_Iter() << "\t" << _heightEq.BCG_Resid()
			<<"\t" <<_timeStep.DT().now * timeScale
			<<"\t" <<_heightEq.ChangeRateOfVariable() / timeScale;
		if(!_outputCalcInfo.Write(oss)) StopOnOutputError();
		MarkWritten(itemCalc.c_str(), atEnd);
	}

	//  Last, so that it holds the schedule as every item has just left it.
	if(!statePath.empty() && (_inputParam->Output.FieldData.NeutralFile || _inputParam->Output.FieldData.VTKFile))
		WriteState(statePath);
}

template<class T>
bool WellMixed_ModelBase<T>::RefineFromState()
{
	Grid *G = _field.GetGridPtr();
	std::vector<char> inside((size_t)G->imax() + 1, 0);
	for(int ii = G->imin(); ii <= G->imax(); ii++)
		inside[(size_t)ii] = (G->CheckGhost(G->template Of<Site>(ii)) == Grid::Inside);
	//  The field file gives the state to the digits it was written with
	//  (18 significant ones, and one rounding in the conversion from the
	//  output units); more than this is a different state.
	const double tol = 1e-12;
	double dh = 0, ds = 0, dc = 0;
	const int rh = _rstate.RefineField("field.h", _field.HH(), inside, tol, dh);
	const int rs = _rstate.RefineField("field.s", _field.SS(), inside, tol, ds);
	const int rc = _rstate.RefineField("field.c", _field.CC(), inside, tol, dc);
	if(rh == 1 && rs == 1 && rc == 1)
	{
		Info<<"Restart: the field is taken to the bit from the state file (the field file"
		    <<" agrees with it to "<<std::max(dh, std::max(ds, dc))<<")."<<nl;
		return true;
	}
	if(rh < 0 || rs < 0 || rc < 0)
		Global::Warning<<"Restart: "<<_inputParam->Restart.FileName<<" is not the state saved beside it"
		               <<" (h, s, c differ by "<<dh<<", "<<ds<<", "<<dc<<" of their largest values):"
		               <<" the field is taken from the field file, the rest from the state file."<<nl;
	return false;
}

//------------------------------------------------------------------
//  Output.Energy (EnergyBudget.h)
//------------------------------------------------------------------
template<class T>
PotentialParams WellMixed_ModelBase<T>::EnergyParams() const
{
	PotentialParams pp;
	pp.capillary = _inputParam->Physics.LaplacePressure.Enable;
	pp.dp    = _NonDim.Dp();
	pp.hstar = _NonDim.HStar();
	pp.n     = _inputParam->Physics.DisjoiningPressure.ExponentN;
	pp.m     = _inputParam->Physics.DisjoiningPressure.ExponentM;
	pp.bo    = _NonDim.Bo();
	pp.nd    = &_NonDim;
	return pp;
}

template<class T>
void WellMixed_ModelBase<T>::EnergyEval()
{
	const PotentialParams pp = EnergyParams();
	const bool px = _inputParam->Boundary.XminBC.Type == "Periodic";
	const bool py = _inputParam->Boundary.YminBC.Type == "Periodic";
	_enF = energy_budget::FreeEnergy(_field, pp, px, py);
	std::vector<double> phi;
	energy_budget::Potential(_field, pp, phi);
	_enD = _heightEq.FluxPower(_field, phi) * CellMeasure(_field.GetGridPtr());
}

template<class T>
void WellMixed_ModelBase<T>::WriteState(const std::string &path)
{
	C_RestartState st;
	std::string model = Name();
	std::replace(model.begin(), model.end(), ' ', '_');
	st.SetText("model", model);
	st.PutClock(_timeStep, OutputItems());
	st.Set("steadyCount", _steadyCount);
	std::vector<double> v(2);
	v[0] = _vol0; v[1] = _resin0;							st.Set("volume0", v);
	st.Set("evaporated", _evaporatedVolume);
	st.Set("supplied", _suppliedVolume);
	st.Set("suppliedResin", _suppliedResin);
	v[0] = _outflowedVolume; v[1] = _outflowedResin;		st.Set("outflowed", v);
	st.Set("adjustCount.height", _heightEq.AdjustTimeCount());
	if(_energyOn)
	{
		std::vector<double> e(3);
		e[0] = _enDiss; e[1] = _enBal; e[2] = _enStep ? 1 : 0;
		st.Set("energy", e);
	}
	st.SetVectorOf("height.old", _heightEq.GetXOld());
	st.SetVectorOf("height.correction", _heightEq.LastCorrection());
	st.SetField("field.h", _field.HH());
	st.SetField("field.s", _field.SS());
	st.SetField("field.c", _field.CC());
	if(!st.Write(path))
	{
		Global::FatalError<<"Could not write the restart state "<<path<<": "<<System::LastIOError()<<"."<<nl;
		StopOnOutputError();
	}
}

//------------------------------------------------------------------
//  The film at each pinned edge, extrapolated to the wall face, for the
//  check below. One overload per field type; each appends
//  (edge name, worst wall thickness, thickest film) for every Pinning edge.
//------------------------------------------------------------------
namespace pinning_check {

struct Edge { const char *name; double wall; double where; };

//  Linear extrapolation from the first two real cells (centers half a cell
//  and one and a half cells in) to the wall face.
inline double AtWall(double h1, double h2) { return h1 - 0.5 * (h2 - h1); }

inline void Collect(const Field_WM1D &f, const BoundaryCondition1D &bc, std::vector<Edge> &out)
{
	const C_Grid1D &G = f.GetGrid();
	const int i0 = G.IminReal(), i1 = G.ImaxReal() - 1;
	if(bc.TypeXmin() == BoundaryCondition::Pinning)
	{
		Edge e = {"Xmin", AtWall(f.HH(G.At<Xedge>(i0)), f.HH(G.At<Xedge>(i0 + 1))), G.X(G.At<Xedge>(i0))};
		out.push_back(e);
	}
	if(bc.TypeXmax() == BoundaryCondition::Pinning)
	{
		Edge e = {"Xmax", AtWall(f.HH(G.At<Xedge>(i1)), f.HH(G.At<Xedge>(i1 - 1))), G.X(G.At<Xedge>(i1))};
		out.push_back(e);
	}
}

inline void Collect(const Field_WM2D &f, const BoundaryCondition2D &bc, std::vector<Edge> &out)
{
	const C_Grid2D &G = f.GetGrid();
	const int i0 = G.IminReal(), i1 = G.ImaxReal() - 1;
	const int j0 = G.JminReal(), j1 = G.JmaxReal() - 1;
	//  Along each pinned edge the worst row is the one reported.
	struct Side { BoundaryCondition::tBoundaryCondition t; const char *name; bool alongY; int a, b; };
	const Side sides[4] = {
		{bc.TypeXmin(), "Xmin", true,  i0, i0 + 1},
		{bc.TypeXmax(), "Xmax", true,  i1, i1 - 1},
		{bc.TypeYmin(), "Ymin", false, j0, j0 + 1},
		{bc.TypeYmax(), "Ymax", false, j1, j1 - 1}};
	for(int k = 0; k < 4; ++k)
	{
		if(sides[k].t != BoundaryCondition::Pinning) continue;
		Edge worst = {sides[k].name, -1e300, 0.0};
		const int lo = sides[k].alongY ? j0 : i0, hi = sides[k].alongY ? j1 : i1;
		for(int m = lo; m <= hi; ++m)
		{
			const Zface z1 = sides[k].alongY ? G.At<Zface>(sides[k].a, m) : G.At<Zface>(m, sides[k].a);
			const Zface z2 = sides[k].alongY ? G.At<Zface>(sides[k].b, m) : G.At<Zface>(m, sides[k].b);
			const double w = AtWall(f.HH(z1), f.HH(z2));
			if(w > worst.wall)
			{
				worst.wall = w;
				worst.where = sides[k].alongY ? G.Y(z1) : G.X(z1);
			}
		}
		out.push_back(worst);
	}
}

inline bool TwoD(const Field_WM1D &) { return false; }
inline bool TwoD(const Field_WM2D &) { return true; }

inline double MaxH(const Field_WM1D &f)
{
	const C_Grid1D &G = f.GetGrid();
	double m = 0;
	for(int i = G.IminReal(); i < G.ImaxReal(); ++i) m = std::max(m, f.HH(G.At<Xedge>(i)));
	return m;
}

inline double MaxH(const Field_WM2D &f)
{
	const C_Grid2D &G = f.GetGrid();
	double m = 0;
	for(int j = G.JminReal(); j < G.JmaxReal(); ++j)
		for(int i = G.IminReal(); i < G.ImaxReal(); ++i) m = std::max(m, f.HH(G.At<Zface>(i, j)));
	return m;
}

} // namespace pinning_check

//------------------------------------------------------------------
//  The rates of the Courant check (CourantCheck.h): the largest, over the
//  real cells, of the sum over the axes of the larger speed of the cell's
//  two faces over the spacing.
//    conc    the velocity the concentration is advected by (the mixed
//            velocity).
//    height  the characteristic speed dQ/dh of the fluxes the scheme
//            carries, the upwinded ones: each is its velocity times h_f,
//            so 3 |u| for the in-plane gravity and the centrifugal force
//            (Q ~ h^3) and 2 |u| for the Marangoni flux (Q ~ h^2); added,
//            as each is upwinded on its own. The capillary and disjoining
//            fluxes, which are not upwinded, are left out.
//------------------------------------------------------------------
namespace courant_check {

template<class V, class S>
inline double HeightSpeedX(const V &u, S s)
{
	return 3.0 * std::fabs(u.U_HorizontalGravity().x(s)) + 3.0 * std::fabs(u.U_Centrifugal().x(s))
	     + 2.0 * std::fabs(u.U_SolutoCapillarity().x(s));
}

template<class V, class S>
inline double HeightSpeedY(const V &u, S s)
{
	return 3.0 * std::fabs(u.U_HorizontalGravity().y(s)) + 3.0 * std::fabs(u.U_Centrifugal().y(s))
	     + 2.0 * std::fabs(u.U_SolutoCapillarity().y(s));
}

inline void Rates(const Field_WM1D &f, double &conc, double &height)
{
	const C_Grid1D &G = f.GetGrid();
	const Velocity_WM1D &u = f.UU();
	const double dx = G.dx();
	conc = height = 0;
	for(int i = G.IminReal(); i < G.ImaxReal(); ++i)
	{
		const Xedge xe = G.At<Xedge>(i);
		const Node m = G.HalfXm(xe), p = G.HalfXp(xe);
		conc   = std::max(conc,   std::max(std::fabs(u.Ux(m)), std::fabs(u.Ux(p))) / dx);
		height = std::max(height, std::max(HeightSpeedX(u, m), HeightSpeedX(u, p)) / dx);
	}
}

inline void Rates(const Field_WM2D &f, double &conc, double &height)
{
	const C_Grid2D &G = f.GetGrid();
	const Velocity_WM2D &u = f.UU();
	const double dx = G.dx(), dy = G.dy();
	conc = height = 0;
	for(int j = G.JminReal(); j < G.JmaxReal(); ++j)
		for(int i = G.IminReal(); i < G.ImaxReal(); ++i)
		{
			const Zface z = G.At<Zface>(i, j);
			const Yedge xm = G.HalfXm(z), xp = G.HalfXp(z);
			const Xedge ym = G.HalfYm(z), yp = G.HalfYp(z);
			conc = std::max(conc, std::max(std::fabs(u.Ux(xm)), std::fabs(u.Ux(xp))) / dx
			                    + std::max(std::fabs(u.Uy(ym)), std::fabs(u.Uy(yp))) / dy);
			height = std::max(height, std::max(HeightSpeedX(u, xm), HeightSpeedX(u, xp)) / dx
			                        + std::max(HeightSpeedY(u, ym), HeightSpeedY(u, yp))  / dy);
		}
}

} // namespace courant_check

template<class T>
void WellMixed_ModelBase<T>::CheckPinnedEdges()
{
	std::vector<pinning_check::Edge> edges;
	pinning_check::Collect(_field, _boundary, edges);
	if(edges.empty()) return;
	//  The check stops an initial condition the input builds. A restart
	//  continues from the film of a field file, most often one a run wrote,
	//  whose edge has evolved under the run's own physics (pushed against
	//  the wall and dried, the film there can pass a tenth of the thickest
	//  legitimately), and the run straight through goes on there
	//  (tools/proptest.py seed 1206). On a restart the same line is
	//  reported, as a warning when it is past the limit, and the run goes on.
	const bool restart = _inputParam->Restart.Enable;

	const double H0 = _NonDim.HeightScale() * 1e6;		// to um
	const double L0 = _NonDim.ExtentScale() * 1e3;		// to mm
	const double hmax = pinning_check::MaxH(_field);
	const double limit = 0.1;							// of the thickest film

	bool bad = false;
	for(size_t k = 0; k < edges.size(); ++k)
	{
		const double frac = (hmax > 0) ? edges[k].wall / hmax : 0.0;
		std::ostringstream msg;
		msg<<std::fixed<<std::setprecision(3)
		   <<"Pinning "<<edges[k].name<<": the "<<(restart ? "restarted" : "initial")
		   <<" film reaches "<<edges[k].wall * H0
		   <<" um at the wall ("<<std::setprecision(1)<<100.0 * frac<<" % of its "
		   <<std::setprecision(3)<<hmax * H0<<" um maximum";
		if(pinning_check::TwoD(_field))
			msg<<", worst at "<<(edges[k].name[0] == 'X' ? "y" : "x")<<" = "<<edges[k].where * L0<<" mm";
		msg<<")";
		if(frac > limit && restart)
			Global::Warning<<msg.str()<<": more than "<<(int)(100 * limit)<<" %, taken as it is"
			               <<" (the check stops an initial condition, not a restart)."<<nl;
		else if(frac > limit)
		{
			bad = true;
			Global::FatalError<<msg.str()<<"."<<nl;
		}
		else
			Info<<msg.str()<<" -- ok."<<nl;
	}
	if(bad)
	{
		Global::FatalError
			<<"Pinning holds the film thickness at zero on its edge, so the initial film"
			<<" has to go to zero there too: more than "<<(int)(100 * limit)<<" % of the"
			<<" film thickness left at the wall is a step inside one cell, which the"
			<<" first time steps cannot get across. Shape the film down to the wall ---"
			<<" Initial.SurfaceShape.EdgeBump with Location on that edge does this --- and"
			<<" keep Initial.SurfaceShape.PrecursorRatio small, since the film is never"
			<<" thinner than that anywhere else."<<nl;
		Global::FatalError.Abort();
	}
}

template<class T>
void WellMixed_ModelBase<T>::BuildRestart()
{
	//----------------------------------------------------
	// Read the restart data
	//----------------------------------------------------
	if(_inputParam->Restart.Enable)
	{
		Info<<"######## Restart enabled. Reading....   ########"<<nl;

		double coef = 1.0;
		if(_inputParam->Output.Dimensioned) coef = 1e-6 / _NonDim.HeightScale();

		if(_inputParam->Restart.ASCII)
		{
			if( ReadRestartAscii(_inputParam->Restart, coef) == false )Global::FatalError.Abort();
		}
		else
		{
			if( C_Restart::Read(_inputParam->Restart, _field, coef) == false )Global::FatalError.Abort();
		}


		double tt = _inputParam->Restart.RestartTime;
		double T0 = _NonDim.TimeScale();
		tt /= T0;
		_timeStep.Set(tt);
		_timeStep.SetCount(itemField.c_str(), _inputParam->Restart.FileIndex);
		_outputVTK.SetFNum(_inputParam->Restart.FileIndex);
		_outputNeutral.SetFNum(_inputParam->Restart.FileIndex);
		_pvdFile.SetFNum(_inputParam->Restart.FileIndex);

		//  What the field file does not hold (see RestartState.h)
		std::string model = Name(), msg;
		std::replace(model.begin(), model.end(), ' ', '_');
		_continuing = _rstate.Load(_inputParam->Restart.FileName.c_str(), model, tt, msg);
		bool exact = false;
		if(_continuing)
		{
			Info<<"Restart: "<<msg<<"."<<nl;
			exact = RefineFromState();
		}
		else
		{
			//  A field file without a state file (for example one made by
			//  hand to start from) is expected; one with a state file
			//  that does not fit is worth a warning.
			//  on the terminal as well: the run is then not the exact continuation asked for
			if(_rstate.Found()) { Global::Warning<<"Restart: "<<msg<<nl; Info<<"Restart: "<<msg<<nl; }
			else                Info<<"Restart: "<<msg<<nl;
			const std::vector<std::string> items = OutputItems();
			for(size_t k = 0; k < items.size(); k++) _timeStep.OutputItem(items[k]).RestartAt(tt);
		}

		//  The state file holds the ghosts as the run had them; otherwise
		//  they come from the boundary conditions.
		if(!exact) _boundary.ApplyToField(_field);
		if(_swConcField) _concEq.SetX(_field.CC());
		_heightEq.SetX(_field.HH());

		if(_continuing)
		{
			_rstate.TakeClock(_timeStep, OutputItems(), _inputParam->Numerics.TimeStepControl.Enable,
			                  itemField.c_str(), _inputParam->Restart.FileIndex);
			_steadyCount = (int)_rstate.Get("steadyCount");
			_evaporatedVolume = _rstate.Get("evaporated");
			_suppliedVolume   = _rstate.Get("supplied", 0, 0.0);
			_suppliedResin    = _rstate.Get("suppliedResin", 0, 0.0);
			_outflowedVolume  = _rstate.Get("outflowed", 0);
			_outflowedResin   = _rstate.Get("outflowed", 1);
			if(_energyOn)
			{
				//  A state written without Output.Energy starts the
				//  dissipated energy afresh, from the restart
				_enDiss = _rstate.Get("energy", 0, 0.0);
				_enBal  = _rstate.Get("energy", 1, 0.0);
				_enStep = _rstate.Get("energy", 2, 0.0) != 0;
			}
			_heightEq.SetAdjustTimeCount((int)_rstate.Get("adjustCount.height"));
			Vector<double> hold(_heightEq.GetXOld().size());
			if(_rstate.GetVectorInto("height.old", hold)) _heightEq.SetXOld(hold);
			if(_rstate.GetVectorInto("height.correction", hold)) _heightEq.SetLastCorrection(hold);
		}

		_field_guess = _field;

		Info<<"ASCII:       "<<(_inputParam->Restart.ASCII ? "true" : "false")<<nl;
		Info<<"File Name:   "<<_inputParam->Restart.FileName<<nl;
		Info<<"RestartTime: "<<_inputParam->Restart.RestartTime<<nl;
		Info<<"FileIndex:   "<<_inputParam->Restart.FileIndex<<nl;
		Info<<"######## Restart setting applied.   ########"<<nl;
	}
}

template<class T>
void WellMixed_ModelBase<T>::Run()
{
	// Cap the time step at the smallest of the output intervals
	double max_dt = _timeStep.MaxDT();
	double val;
	double T0 = _NonDim.TimeScale();
	//  Listed field outputs cap each step by themselves (ApplyListedTimes)
	if(!_inputParam->Output.FieldData.Listed())
	{ val = _inputParam->Output.FieldData.TimeInterval / T0;	if(val < max_dt) max_dt = val; }
	val = _inputParam->Output.ProbeList.TimeInterval / T0;		if(val < max_dt) max_dt = val;
	val = _inputParam->Output.GeneralHistory.TimeInterval / T0;	if(val < max_dt) max_dt = val;
	val = _inputParam->Output.CalcInfo.TimeInterval / T0;		if(val < max_dt) max_dt = val;

	_timeStep.SetMaxDT(max_dt);

	if(_timeStep.MaxDT() < _timeStep.DT().initial)
	{
		// The case where dt is larger than dtmax
		Global::Info<<"initial DT is larger than Max-DT. DT is set to Max-DT ("<<_timeStep.MaxDT()<<")."<<nl;
		_timeStep.ChangeDT(_timeStep.MaxDT());
	}

	BuildRestart();
	CheckPinnedEdges();

	//Compute the initial velocity field
	Global::Info << " calVelocity"<<nl;
	_heightEq.CalVelocity();
	Global::Info << " set field"<<nl;
	_field = _field_guess;
	//  (The concentration equation forms the explicit half of every step
	//  from the state the step starts from, HeightEquationBase::BeforeNKLoop;
	//  a run continued from a state needs nothing more of it.)

	//Probe the initial values
	Global::Info << " Probes "<<nl;
	_field.Accept(_fieldInfo);
	for(ProbeList::iterator it = _probeList.begin(); it != _probeList.end(); ++it)
		_field.Accept(*it);

	_fieldInfoIni = _fieldInfo;
	if(_energyOn) EnergyEval();
	_vol0   = _fieldInfoIni.FluidVol();
	_resin0 = _fieldInfoIni.ResinVol();
	if(_continuing)
	{
		//  The balances go on against the initial state of the whole run
		_vol0   = _rstate.Get("volume0", 0);
		_resin0 = _rstate.Get("volume0", 1);
		//  History, CalcInfo, the probes and the .pvd go on in the files the
		//  run being continued wrote, when the restart writes into the same
		//  directory: their rows up to the restart time are kept.
		const double ts = _inputParam->Output.Dimensioned ? _NonDim.TimeScale() : 1.0;
		const double tmax = _timeStep.NowTime() * ts * (1.0 + 2e-6);
		if(!_outputHistory.Resume(tmax)) StopOnOutputError();
		if(!_outputCalcInfo.Resume(tmax)) StopOnOutputError();
		if(_energyOn && !_outputEnergy.Resume(tmax)) StopOnOutputError();
		for(size_t k = 0; k < _probeOut.size(); k++) if(!_probeOut[k]->Resume(tmax)) StopOnOutputError();
		if(!_pvdFile.Resume(tmax)) StopOnOutputError();
		_rstate.Clear();
	}

	//  The Courant numbers of the limited schemes: watched and reported,
	//  never acted on (CourantCheck.h)
	_courantConc.Build("Numerics.Convection.Concentration", _inputParam->Numerics.Convection.Concentration,
	                   "the concentration", _swConcField && _inputParam->Physics.ConcField.EnableConvection
	                                        && _concEq.Scheme().IsLimited());
	_courantHeight.Build("Numerics.Convection.Height", _inputParam->Numerics.Convection.Height,
	                     "the film thickness", _heightEq.LimitedUpwind());
	{
		double rc, rh;
		courant_check::Rates(_field, rc, rh);
		_courantConc.Start(rc);
		_courantHeight.Start(rh);
	}

	//The time-step loop
	Global::Info << "---- Time step loop ----"<<nl;
	const double loopStart = Perf::NowSeconds();		// for measurement only; unrelated to the numerics
	_loopStartTime = _timeStep.NowTime();
	_sumStepDt = 0;
	//  A run continued from the step at which the steady-state test was met
	//  -- it was stopped on EndTime or Numerics.MaxSteps at that same step --
	//  takes no further step, as the uninterrupted run took none.
	const bool steadyAlready = _continuing && _inputParam->Numerics.SteadyStateCriterion > 0
	                        && _steadyCount > _inputParam->Numerics.SteadyStateSteps;
	if(steadyAlready) Info<<"\n############## Steady state (when restarted). #############"<<nl;
	//  The same for the dry-out test: a run stopped on Numerics.MaxSteps at
	//  the step that dried out, and continued, must not take a step the
	//  uninterrupted run, which stopped there, did not take
	//  (tools/proptest.py). The field is the one that run ended on.
	const bool dryAlready = _continuing && !steadyAlready && DriedOut();
	if(dryAlready) Info<<"\n############## Dry out (when restarted). #############"<<nl;
	while(!steadyAlready && !dryAlready && _timeStep.IsRunning())
	{
		{ YUZEN_PROF("output"); FreqOutput(); }
		_timeStep.ApplyListedTimes();

		Info<<".";
		{ YUZEN_PROF("step"); _heightEq.TimeAdvance(_NK_iter, _NK_resid); }
		++_nSteps;
		_sumNKIter  += _NK_iter;
		_sumBCGIter += _heightEq.BCG_IterOfStep();

		//  Accumulate the evaporated, the supplied and the total outflowed
		//  volume. The time step used is the one the step was solved with
		//  (StepDt, the dt the source term was assembled with; DT().now
		//  until ++), not (t+dt)-t, which does not equal dt once rounding is
		//  taken into account.
		const double dtStep = _heightEq.StepDt();
		_sumStepDt += dtStep;
		_field.Accept(_fieldInfo);
		//  What crossed the open edges (FreshInlet, Outlet) in the step, as
		//  the step's own face fluxes put it through them: the volume flux of
		//  the film-height equation and the solute flux of the concentration
		//  equation through every edge face, at the step's Crank-Nicolson
		//  weighting, so that the balances of History.dat close to rounding
		//  on an open domain as on a closed one. Only the step committed
		//  counts (a rejected attempt is never completed).
		{
			double edgeResin = 0;
			_outflowedVolume += _heightEq.EdgeOutflowVolume(CellMeasure(_field.GetGridPtr()), &edgeResin);
			_outflowedResin  += edgeResin;
		}
		_evaporatedVolume += CalEvaporatedVolume() * dtStep;
		_suppliedVolume   += _heightEq.SuppliedVolume(CellMeasure(_field.GetGridPtr()));
		_suppliedResin    += _heightEq.SuppliedResin(CellMeasure(_field.GetGridPtr()));
		//  The solute the commit moved from above c = 1, and what it cut;
		//  a cut beyond rounding is reported as it happens, once
		{
			const double cm = CellMeasure(_field.GetGridPtr());
			const double mv = _heightEq.ExcessMoved(cm), cut = _heightEq.ExcessTrimmed(cm);
			if(mv > 0.0) { _excessMoved += mv; ++_excessSteps; }
			_excessTrimmed += cut;
			if(cut > 1e-12 * std::fabs(_resin0) && !_warnedTrim)
			{
				std::ostringstream m;
				m<<std::setprecision(3)<<"The concentration passed 1 at step "<<_nSteps<<" (t = "
				 <<(_timeStep.NowTime() + dtStep) * T0<<" s) with more solute than the film has room for: "
				 <<cut / _resin0<<" of the initial solute was cut at c = 1 and is missing from ResinVolErr."
				 <<" Reported once; the total is given at the end of the run.";
				Global::Warning<<m.str()<<nl;
				_warnedTrim = true;
			}
			//  The same below c = 0
			const double dm = _heightEq.DeficitMoved(cm), add = _heightEq.DeficitAdded(cm);
			if(dm > 0.0) { _deficitMoved += dm; ++_deficitSteps; }
			_deficitAdded += add;
			if(add > 1e-12 * std::fabs(_resin0) && !_warnedAdd)
			{
				std::ostringstream m;
				m<<std::setprecision(3)<<"The concentration fell below 0 at step "<<_nSteps<<" (t = "
				 <<(_timeStep.NowTime() + dtStep) * T0<<" s) with less solute than the film needs to fill it: "
				 <<add / _resin0<<" of the initial solute was made up at c = 0 and is an error of ResinVolErr."
				 <<" Reported once; the total is given at the end of the run.";
				Global::Warning<<m.str()<<nl;
				_warnedAdd = true;
			}
		}
		//  In 1-D, the liquid (and the solute) the centrifugal term flung
		//  sideways into the y that is not solved for: it has left the
		//  strip, through its long sides, and is counted as outflowed
		{
			double sideResin = 0;
			_outflowedVolume += _heightEq.SideOutflowVolume(CellMeasure(_field.GetGridPtr()), &sideResin);
			_outflowedResin  += sideResin;
		}
		//  The free energy and the dissipation of the state the step
		//  reached, and the energy dissipated over it at the step's
		//  Crank-Nicolson weighting
		if(_energyOn)
		{
			const double f0 = _enF, d0 = _enD;
			EnergyEval();
			_enDiss += 0.5 * (d0 + _enD) * dtStep;
			_enBal = (_enF - f0) / dtStep + 0.5 * (d0 + _enD);
			_enStep = true;
		}

		++_timeStep;

		if(_courantConc.On() || _courantHeight.On())
		{
			double rc, rh;
			courant_check::Rates(_field, rc, rh);
			const double tEnd = _timeStep.NowTime() * T0;
			_courantConc.Step(rc, dtStep, _nSteps, tEnd);
			_courantHeight.Step(rh, dtStep, _nSteps, tEnd);
		}

		//Check the termination conditions
		if(CheckTermination()) break;

		//  Ctrl-C (SIGINT) or SIGTERM: stop after this step, and write the
		//  state reached below as at the end of a run, so that it can be
		//  continued from (util/YuzenUtil.h, Interrupt)
		if(Interrupt::Requested())
		{
			Info<<"\n############## Stopped by "<<Interrupt::Name()<<". #############"<<nl;
			break;
		}
	}

	_loopSeconds = Perf::NowSeconds() - loopStart;

	//  The loop writes the state it is about to advance, so without this the
	//  final state -- the one at EndTime, which is the answer the run was
	//  asked for -- is kept only when the output schedule happens to land on
	//  it exactly. Anything already written at this time is not written twice.
	//  A diverged run is the exception: what it holds is not a solution.
	if(!_diverged) FreqOutput(true);

	if(!_pvdFile.Close()) StopOnOutputError();

	//  Record the measured figures for the time-step loop alone in
	//  JobInfo.txt. The scaling measurements (tools/scaling/) read these.
	JobInfo<<nl<<nl<<"###### Performance ######"<<nl;
	JobInfo.Add("TimeLoopSeconds", _loopSeconds,        "Time loop [s]");
	JobInfo.Add("TimeStepCount",   (int)_nSteps,        "Time steps");
	JobInfo.Add("NKIterTotal",     (int)_sumNKIter,     "Newton iterations");
	JobInfo.Add("BCGIterTotal",    (int)_sumBCGIter,    "Linear iterations");
	JobInfo.Add("NewtonStagnationCuts", (int)_heightEq.StagnationCuts(), "Steps failed: the Newton correction stagnated");
	_courantConc.Report("conc");
	_courantHeight.Report("height");
	ReportExcess();
	//  The sum of the steps the equations were solved with, and the times
	//  the clock went from and to, on the time scale of the outputs and to
	//  every digit: the clock is the sum of the steps solved, to rounding
	//  (tools/features.py, adaptive_clock).
	{
		const double ts = _inputParam->Output.Dimensioned ? _NonDim.TimeScale() : 1.0;
		auto full = [](double v) { std::ostringstream o; o<<std::setprecision(17)<<v; return o.str(); };
		JobInfo.Add("StepDtSum",     full(_sumStepDt * ts),          "Sum of the steps solved");
		JobInfo.Add("LoopStartTime", full(_loopStartTime * ts),      "Time at the start of the loop");
		JobInfo.Add("LoopEndTime",   full(_timeStep.NowTime() * ts), "Time at the end of the loop");
	}

	if(_diverged)
		Info<<"\n############## Program terminated after divergence. #############"<<nl;
	else if(Interrupt::Requested())
		Info<<"\n############## Program stopped by "<<Interrupt::Name()<<"; its output is complete. #############"<<nl;
	else
		Info<<"\n############## Program normally terminated. #############"<<nl;

}

//------------------------------------------------------------------
//  The solute moved from above c = 1, and cut there, and moved into cells
//  below c = 0, and made up there, over the run: to the run log and
//  JobInfo.txt, as fractions of the initial solute; a cut or a make-up
//  beyond rounding to Warnings.txt as well. An unbounded convection scheme
//  over- or undershooting at a sharp front moves solute; nothing else does
//  beyond rounding (the evaporation is limited to the solvent there is).
//------------------------------------------------------------------
template<class T>
void WellMixed_ModelBase<T>::ReportExcess()
{
	if(!_swConcField) return;
	const double r0 = (_resin0 != 0.0) ? std::fabs(_resin0) : 1.0;
	const double moved = _excessMoved / r0, cut = _excessTrimmed / r0;
	if(_excessSteps > 0 || _excessTrimmed > 0)
	{
		std::ostringstream v;
		v<<std::setprecision(3)<<moved<<" of the initial in "<<_excessSteps<<" steps (cut "<<cut<<")";
		Info<<"Solute moved from c > 1 to the cells around: "<<v.str()<<nl;
	}
	if(_excessTrimmed > 1e-12 * r0)
	{
		std::ostringstream m;
		m<<std::setprecision(3)<<"Solute cut at c = 1 over the run: "<<cut
		 <<" of the initial (missing from ResinVolErr).";
		Global::Warning<<m.str()<<nl;
	}
	JobInfo.Add("ConcExcessMoved", moved, "Solute moved from c > 1 [-]");
	JobInfo.Add("ConcExcessCut",   cut,   "Solute cut at c = 1 [-]");

	const double dmoved = _deficitMoved / r0, added = _deficitAdded / r0;
	if(_deficitSteps > 0 || _deficitAdded > 0)
	{
		std::ostringstream v;
		v<<std::setprecision(3)<<dmoved<<" of the initial in "<<_deficitSteps<<" steps (made up "<<added<<")";
		Info<<"Solute moved to c < 0 from the cells around: "<<v.str()<<nl;
	}
	if(_deficitAdded > 1e-12 * r0)
	{
		std::ostringstream m;
		m<<std::setprecision(3)<<"Solute made up at c = 0 over the run: "<<added
		 <<" of the initial (an error of ResinVolErr).";
		Global::Warning<<m.str()<<nl;
	}
	JobInfo.Add("ConcDeficitMoved", dmoved, "Solute moved to c < 0 [-]");
	JobInfo.Add("ConcDeficitAdded", added,  "Solute made up at c = 0 [-]");
}

template<class T>
double WellMixed_ModelBase<T>::CalEvaporatedVolume()
{
	if (! _inputParam->Physics.Evaporation.Enable) return 0;

	//  The rate the film-height equation took out over the step just
	//  solved (its Crank-Nicolson average, capped as it was capped), so
	//  that VolErr balances to rounding with evaporation too. The uncapped
	//  rate of the end of the step would not: the caps -- at half the film,
	//  at the solvent there is -- would make the difference a balance error.
	return _heightEq.EvaporationRate(CellMeasure(_field.GetGridPtr()));
}

//------------------------------------------------------------------
//  The dry-out test, on the field _fieldInfo was last taken of. The
//  threshold and what is measured follow EndConcType and EndConcValue.
//------------------------------------------------------------------
template<class T>
bool WellMixed_ModelBase<T>::DriedOut() const
{
	double EndConcValue = _inputParam->Physics.Evaporation.EndConcValue;
	double checkVal = 0.0;
	if(_inputParam->Physics.Evaporation.EndConcType == "mean")
	{
		checkVal = _fieldInfo.AveConc();
	}
	else if(_inputParam->Physics.Evaporation.EndConcType == "max")
	{
		checkVal = _fieldInfo.MaxC();
	}
	else
	{
		checkVal = _fieldInfo.MinC();
	}
	return checkVal > EndConcValue;
}

template<class T>
bool WellMixed_ModelBase<T>::CheckTermination()
{
	//  The count of the steady-state test is taken first, before any test
	//  that can end the run: the step that ends it (on EndTime or on
	//  Numerics.MaxSteps) is counted like every other, so that a run
	//  continued from the state written after that step counts on from
	//  where an uninterrupted run would be. (Taken after the time test, the
	//  last step of each interrupted part would go uncounted.)
	const double steady_eps  = _inputParam->Numerics.SteadyStateCriterion;
	const int    steady_need = _inputParam->Numerics.SteadyStateSteps;
	if(steady_eps > 0)
	{
		double tspan = _timeStep.EndTime() - _timeStep.StartTime();
		if(!(tspan > 0)) tspan = 1.0;		// no scaling when EndTime is unspecified
		if(_heightEq.ChangeRateOfVariable() * tspan < steady_eps) _steadyCount++;
		else _steadyCount = 0;
	}

	//------------------------------------------------------------------
	//  Detection of divergence, before every other test: once the field
	//  holds NaN the comparisons of the dry-out and steady-state tests lose
	//  their meaning, and the last step of a run is tested as well as any
	//  other.
	//
	//  Without it, a step whose Newton iteration does not converge at a
	//  fixed dt leaves a film that goes negative and then NaN, and the run
	//  goes on to EndTime and ends with exit code 0 and NaN in History.dat
	//  and the field files (tools/proptest.py).
	//  The thresholds are those of ConcZ: loose, so that what a healthy run
	//  strays by through discretization error passes, and an excursion of an
	//  order of magnitude does not. _fieldInfo was taken of this step's field
	//  just before the call.
	//------------------------------------------------------------------
	{
		const double vol  = _fieldInfo.FluidVol();
		const double maxC = _fieldInfo.MaxC();
		//  The lowest c of the start of the step as well: a field file read
		//  with a cell far below 0 must stop the run at its first step, which
		//  may mix the cell back into range by itself (or fill it from the
		//  cells around, RedistributeExcess)
		const double minC = _swConcField ? std::min(_fieldInfo.MinC(), _heightEq.MinStartConc())
		                                 : _fieldInfo.MinC();
		const double minH = _fieldInfo.MinH();
		const char *why = 0;
		if(!std::isfinite(vol) || !std::isfinite(maxC) || !std::isfinite(minC) || !std::isfinite(minH))
			why = "the field is not a number";
		else if(_swConcField && (maxC > 1.0 + kConcSlack || minC < -kConcSlack))
			why = "the concentration left [0,1] by more than the allowed slack";
		else if(minH < -kHeightSlack)
			why = "the film thickness went negative";
		else if(vol < 0.0)
			why = "the fluid volume went negative";
		if(why)
		{
			Info<<"\n############## Diverged: "<<why<<". #############"<<nl;
			Global::Warning<<Name()<<": diverged at step "<<_nSteps<<": "<<why
			               <<" (FluidVol = "<<vol<<", C = ["<<minC<<", "<<maxC
			               <<"], minH = "<<minH<<")."<<nl;
			_diverged = true;
			return true;
		}
	}

	if(_timeStep.IsRunning() == false)
	{
		//  So that the reason for stopping is not mistaken, report
		//  separately whether the run ended on time or was cut short by
		//  Numerics.MaxSteps.
		if(_timeStep.TimeStep() > _timeStep.EndTimeStep())
			Info<<"\n############## Step count reached Numerics.MaxSteps. #############"<<nl;
		else
			Info<<"\n############## Time reached to the EndTime. #############"<<nl;
		return true;
	}


	// Concentration check
	if(DriedOut())
	{
		Info<<"\n############## Dry out. #############"<<nl;
		return true;
	}

	//------------------------------------------------------------------
	//  The steady-state test.
	//
	//  ** Do not base this test on the change per time step. **
	//  max|dh|/max|h| gets small simply because the time step is small, so
	//  in a computation with dt = 1e-11 the solution is judged "steady"
	//  right after the start and the run ends silently.
	//
	//  What is looked at here is
	//
	//      (relative rate of change per unit time) x (total time span of
	//                                                the computation)
	//    = "how much would the shape change over all the remaining time, at
	//       the rate it is changing now"
	//
	//  which is a non-dimensional quantity. It depends neither on the time
	//  step nor on the scaling. The default of 1e-3 means "it would move
	//  only 0.1 % even if the run went all the way to EndTime". The test
	//  can be switched off with Numerics.SteadyStateCriterion <= 0.
	//------------------------------------------------------------------
	if(steady_eps > 0)
	{
		if(_steadyCount > steady_need)
		{
			Info<<"\n############## Steady state. #############"<<nl;
			return true;
		}
	}
	return false;
}

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_WELL_MIXED_MODEL_BASE_INCLUDED_
