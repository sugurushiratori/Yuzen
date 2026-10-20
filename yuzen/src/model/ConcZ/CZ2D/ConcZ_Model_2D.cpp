/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ConcZ_Model_2D.cpp : implementation of the driver of the ConcZ model
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "Profile.h"
#include "ConcZ_Model_2D.h"
#include "Restart.h"
#include "ManufacturedSource.h"
#include <iomanip>
#include <sstream>

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

using namespace yuzen::util;
using namespace Global;



//  Allowance of the divergence test. The physical range is [0,1] for the
//  concentration and positive for the film thickness, and the run is taken to
//  have diverged only when it leaves that range by an order of magnitude.
const double ConcZ_Model_2D::kConcSlack   = 0.1;
const double ConcZ_Model_2D::kHeightSlack = 1.0e-3;

ConcZ_Model_2D::ConcZ_Model_2D()
{
	_stepStatus = 0;
	_diverged = false;
	_name = "ConcZ-2D";
	_iter_heightEq = 0;
	_resid_heightEq = 0;
	_steadyCount = 0;
	_nSteps = 0;
	_sumNKIter = 0;
	_sumBCGIter = 0;
	_sumCoupling = 0;
	_deficitMoved = _deficitAdded = 0;
	_deficitSteps = 0;
	_warnedAdd = false;
	_loopSeconds = 0;
	_sumStepDt = 0;
	_loopStartTime = 0;
	_continuing = false;
	_vol0 = _resin0 = 0;
	_evaporatedVolume = _outflowedVolume = _outflowedResin = 0;
}

void ConcZ_Model_2D::Build(Input& inp)
{
	//----------------------------------------------------
	//Copy the input parameters.
	//----------------------------------------------------
	_inputParam = &inp;
	_swConcField = true;

	//  With Grid.Dimension = "xh" the film is one-dimensional, so Ny has to be
	//  1. With "xyh" it is two-dimensional (pseudo three-dimensional: the film
	//  extends in x-y and the concentration also varies across the thickness),
	//  and Ny is used as it stands.
	if(inp.Grid.Dimension == "xh" && inp.Grid.Ny != 1)
	{
		Global::Warning<<"Grid.Dimension = \"xh\" requires Ny = 1. Ny=1 is forced."<<nl;
		inp.Grid.Ny = 1;
	}

	//----------------------------------------------------
	//Build the fluid properties.
	//----------------------------------------------------
	_fluid.Build(_inputParam->Fluid);

	//----------------------------------------------------
	//Build the non-dimensional numbers.
	//----------------------------------------------------
	_NonDim.Build(inp, inp.Initial.SurfaceShape.TotalHeight);

	//----------------------------------------------------
	//Make the grid non-dimensional.
	//----------------------------------------------------
	_gridNonDim = _inputParam->Grid;
	_gridNonDim.Xmin /= _NonDim.ExtentScale();
	_gridNonDim.Xmax /= _NonDim.ExtentScale();
	_gridNonDim.Ymin /= _NonDim.ExtentScale();
	_gridNonDim.Ymax /= _NonDim.ExtentScale();

	//----------------------------------------------------
	// Build the field (the substrate shape is set inside it as well).
	//----------------------------------------------------
	_field.Build(_gridNonDim, inp.Substrate, _NonDim.ExtentScale());
	C_Grid2D *G = _field.GetGridPtr2();
	int nx = G->nx();

	//  Name it according to whether the film is one- or two-dimensional
	//  (the class behind it is the same).
	_name = (G->ny() != 1) ? "ConcZ-3D" : "ConcZ-2D";

	//----------------------------------------------------
	//Build the boundary conditions.
	//----------------------------------------------------
	BoundaryCondition::tBoundaryCondition tXmin, tXmax, tYmin, tYmax;
	tXmin = BoundaryCondition::FromStr(inp.Boundary.XminBC.Type);
	tXmax = BoundaryCondition::FromStr(inp.Boundary.XmaxBC.Type);
	tYmin = BoundaryCondition::FromStr(inp.Boundary.YminBC.Type);
	tYmax = BoundaryCondition::FromStr(inp.Boundary.YmaxBC.Type);
	_boundary.Build(nx, G->ny(), tXmin, tXmax, tYmin, tYmax);

	//----------------------------------------------------
	// Initialize the surface shape.
	//----------------------------------------------------
	C_InitialSurface::Set(_field, _inputParam->Initial.SurfaceShape, _inputParam->Substrate, _NonDim.ExtentScale());

	//----------------------------------------------------
	//Scaling and initial setting
	//----------------------------------------------------
	for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf)
	{
		_field.HH(zf) /= _inputParam->Initial.SurfaceShape.TotalHeight;
		_field.SS(zf) /= _inputParam->Initial.SurfaceShape.TotalHeight;
	}
	//  With a manufactured solution the film and the substrate start from it
	//  (verification only)
	{
		C_ManufacturedSource mms;
		const bool twoD = (G->ny() != 1);
		mms.Build(*_inputParam, _NonDim, twoD, true);
		if(mms.Enabled())
			for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf)
			{
				const double xx = G->X(zf), yy = twoD ? G->Y(zf) : 0.0;
				_field.HH(zf) = mms.Hnd(xx, yy, 0.0);
				_field.SS(zf) = mms.Subnd(xx, yy);
			}
	}
	//  The ghost columns of the substrate and of the film, the way the
	//  well-mixed model fills them (BoundaryCondition2D::ApplyToField):
	//  mirrored at a wall, taken from the other end across a periodic pair.
	//  They held the shape functions evaluated outside the domain: across a
	//  periodic pair that is right only for a shape that repeats with the
	//  domain, and at a wall it is the shape's own continuation rather than
	//  its mirror image.
	FillGhostColumns();

	C_Grid3D *G3 = _field.GetGridPtr3();
	for(Cell cl = G3->First<Cell>(); cl != G3->Last<Cell>(); ++cl){
		_field.CC(cl) = _inputParam->Initial.Concentration.Value;
		_field.Ux(G3->HalfXm(cl)) = 0;
		_field.Uy(G3->HalfYm(cl)) = 0;
		_field.Uz(G3->HalfZm(cl)) = 0;
	}
	//  A manufactured concentration, at the height zeta of each level (the
	//  ghost levels included: its part that varies through the film is even
	//  about the substrate and the surface, as the no-flux conditions there
	//  make the ghosts)
	{
		C_ManufacturedSource mms;
		const bool twoD = (G->ny() != 1);
		mms.Build(*_inputParam, _NonDim, twoD, true);
		const Input::I_Numerics::I_ManufacturedSolution &m = _inputParam->Numerics.ManufacturedSolution;
		if(mms.Enabled() && !(Math::IsZero(m.ConcAmplitude) && Math::IsZero(m.ConcStratAmplitude)))
			for(int kz = 0; kz < G3->nz(); kz++)
			for(int jy = 0; jy < G3->ny(); jy++)
			for(int ix = 0; ix < G3->nx(); ix++)
			{
				const Zface zf = G->At<Zface>(ix, jy);
				const double zeta = G3->dz() * (-1.5 + (double)kz);
				_field.CC(G3->At<Cell>(ix, jy, kz)) = mms.CZnd(G->X(zf), twoD ? G->Y(zf) : 0.0, zeta, 0.0);
			}
	}

	//----------------------------------------------------
	//Build the time steps.
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
	//  (two billion) is unlimited in practice, so by default the run ends on
	//  EndTime alone.
	if(inp.Numerics.MaxSteps > 0)
		_timeStep.SetMaxSteps((std::uint64_t)inp.Numerics.MaxSteps);

	//----------------------------------------------------
	//Build the equation of the concentration field.
	//----------------------------------------------------
	_concEq.Build(inp, &_fluid, &_timeStep, &_field, &_NonDim, &_boundary);
	_concEq.SetX(_field.CC());
	_concEq.Initialize();

	//----------------------------------------------------
	//Build the equation of the height function.
	//----------------------------------------------------
	_heightEq.Build(inp, &_fluid, &_timeStep, &_concEq, &_field, &_NonDim, &_boundary);
	_heightEq.SetX(_field.HH());
	_heightEq.Initialize();



	//----------------------------------------------------
	//Set up the output data.
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
	_timeStep.Add_OutputItem(itemProbe.c_str(), _inputParam->Output.ProbeList.TimeInterval / T0);
	_timeStep.Add_OutputItem(itemHistory.c_str(), _inputParam->Output.GeneralHistory.TimeInterval / T0);
	_timeStep.Add_OutputItem(itemCalc.c_str(), _inputParam->Output.CalcInfo.TimeInterval / T0);


	//----------------------------------------------------
	// Build the output classes.
	//----------------------------------------------------
	std::vector<std::string> histLabels, calcLabels;
	//  *** The labels must match the columns that are actually written. ***
	//  The 14 labels of the well-mixed side, in the same places: FreqOutput
	//  writes EvaporatedVol / OutflowedVol / OutflowedResin / VolErr /
	//  ResinVolErr after the others. A label without its column shifts every
	//  column after it (tools/proptest.py checks that the counts agree).
	//  As on the well-mixed side: History.dat keeps the raw scales (m^3, m),
	//  unlike the field files, so the header says which it is.
	const bool dim = _inputParam->Output.Dimensioned;
	//  In 1-D, per meter of the width of the strip
	const char *uV = dim ? (_inputParam->Grid.Dimension == "xh" ? " [m^3/m]" : " [m^3]") : " [-]";
	const char *uH = dim ? " [m]"   : " [-]";
	histLabels.push_back(std::string("Volume") + uV);
	histLabels.push_back(std::string("ResinVol") + uV);
	histLabels.push_back("AveC [-]");
	histLabels.push_back(std::string("maxHS") + uH);
	histLabels.push_back(std::string("minHS") + uH);
	histLabels.push_back(std::string("maxH") + uH);
	histLabels.push_back(std::string("minH") + uH);
	histLabels.push_back("maxC [-]");	histLabels.push_back("minC [-]");
	//  The balance, as on the well-mixed side. ConcZ has no open edge, so
	//  OutflowedVol is the liquid the 1-D centrifugal term flings sideways,
	//  into the y that is not solved for (zero in 2-D), and OutflowedResin
	//  the solute it carries.
	histLabels.push_back(std::string("EvaporatedVol") + uV);
	histLabels.push_back(std::string("OutflowedVol") + uV);
	histLabels.push_back(std::string("OutflowedResin") + uV);
	histLabels.push_back("VolErr [-]");	histLabels.push_back("ResinVolErr [-]");

	//  CalcInfo.dat: four labels after the time, for the four columns
	//  written, not the six of the well-mixed side (NK_iter and NK_resid
	//  have no column here). The film-height equation of ConcZ takes no
	//  Newton iteration: what is written is the linear iteration count and
	//  residual of the step, then dt and the rate.
	calcLabels.push_back("BCG_iter");	calcLabels.push_back("BCG_resid");
	calcLabels.push_back("deltaT");
	calcLabels.push_back("ChangeRate");	// relative rate of change of the film thickness [1/s]

	_outputVTK.Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), "field", 0);
	_outputNeutral.Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), "field", 0);
	_outputHistory.Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), "History", (int)histLabels.size(), histLabels);
	_outputCalcInfo.Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), "CalcInfo", (int)calcLabels.size(), calcLabels);

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
	//  Record how it was solved in JobInfo.txt (as WellMixedModelBase does).
	//----------------------------------------------------
	JobInfo<<"\n###### Solver Info ######\n";
	JobInfo.Add("Model",         std::string(Name()),     "Model");
	JobInfo.Add("BackendHeight", _heightEq.BackendName(), "Linear solver <height>");
	JobInfo.Add("BackendConc",   _concEq.BackendName(),   "Linear solver <conc>");
	JobInfo.Add("ConvectionConc", C_Interpolation_base4::Key(_concEq.Scheme().Scheme()), "Convection scheme <conc>");

	//----------------------------------------------------
	// Information output
	//----------------------------------------------------
	Info()<<std::scientific;
	Info<<"###### Solver #####"<<nl;
	Info<<" Linear solver <height> = "<<_heightEq.BackendName()<<nl;
	Info<<" Linear solver <conc>   = "<<_concEq.BackendName()<<nl;
	Info<<"######## Domain & Grids ########"<<nl;
	Info<<" (Nx, Ny) = ("<<_inputParam->Grid.Nx<<" , "<<_inputParam->Grid.Ny <<")" <<nl;
	Info<<"xrange["<<_inputParam->Grid.Xmin<<" , "<<_inputParam->Grid.Xmax <<"]" <<nl;
	Info<<"yrange["<<_inputParam->Grid.Ymin<<" , "<<_inputParam->Grid.Ymax <<"]" <<nl;
	Info<<_NonDim.ToString()<<nl;
	Info<<"###### Non-dimensional time info #####"<<nl;
	Info<<" deltaT    = "<<deltaT <<nl;
	Info<<" startTime = "<<startTime <<nl;
	Info<<" endTime   = "<<endTime <<nl;
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


	//----------------------------------------------------
	//Build the probes.
	//----------------------------------------------------
	BuildProbes();

}

void ConcZ_Model_2D::BuildProbes()
{
	//-------------------------
	//  Find the probe points.
	//-------------------------
	double px, py, pz;
	double xx, yy, zz;
	C_Grid3D *G = _field.GetGridPtr3();
	C_Grid2D *G2p = _field.GetGridPtr2();
	for(unsigned int i = 0; i < _inputParam->Output.ProbeList.Probe.size(); i++)
	{
		//  px and py are meters, like the well-mixed models', made
		//  dimensionless before they are compared with the grid; pz is the
		//  height above the substrate in meters, turned into the zeta of the
		//  grid with the film thickness the probe's column starts with. The
		//  nearest cell within one mesh width is taken, not the first found.
		px = _inputParam->Output.ProbeList.Probe[i].px / _NonDim.ExtentScale();
		py = _inputParam->Output.ProbeList.Probe[i].py / _NonDim.ExtentScale();
		pz = _inputParam->Output.ProbeList.Probe[i].pz / _NonDim.HeightScale();
		const bool twoDp = (G2p->ny() != 1);
		int bx = -1, by = -1; double dbest = 1e300;
		for(Zface zf = G2p->First<Zface>(); zf != G2p->Last<Zface>(); ++zf)
		{
			const double ex = fabs(G2p->X(zf) - px), ey = twoDp ? fabs(G2p->Y(zf) - py) : 0.0;
			if(ex < G2p->dx() && (!twoDp || ey < G2p->dy()))
			{
				const double d = (ex/G2p->dx())*(ex/G2p->dx()) + (twoDp ? (ey/G2p->dy())*(ey/G2p->dy()) : 0.0);
				if(d < dbest) { dbest = d; bx = G2p->ix(zf); by = G2p->iy(zf); }
			}
		}
		if(bx < 0) continue;
		const double hcol = _field.HH(G2p->At<Zface>(bx, by));
		const double zeta = (hcol > 0) ? pz / hcol : 0.0;
		int bz = -1; double dz = 1e300;
		for(Cell cl = G->First<Cell>(); cl != G->Last<Cell>(); ++cl)
		{
			if(G->ix(cl) != bx || G->iy(cl) != by || G->CheckGhost(cl) != C_Grid3D::Inside) continue;
			zz = G->Z(cl);
			if(fabs(zz - zeta) < dz) { dz = fabs(zz - zeta); bz = G->iz(cl); }
		}
		(void)xx; (void)yy;
		if(bz >= 0) _probeList.push_back(PointProbe(bx, by, bz));
	}

	std::vector<std::string> probeLabels;
	//  Units in the labels, and the same ones the well-mixed probes use.
	if(_inputParam->Output.Dimensioned)
	{
		probeLabels.push_back("HS [um]");	probeLabels.push_back("Conc [-]");
		probeLabels.push_back("Ux [m/s]");	probeLabels.push_back("Uy [m/s]");
		probeLabels.push_back("HH [um]");	probeLabels.push_back("SS [um]");
	}
	else
	{
		probeLabels.push_back("HS [-]");	probeLabels.push_back("Conc [-]");
		probeLabels.push_back("Ux [-]");	probeLabels.push_back("Uy [-]");
		probeLabels.push_back("HH [-]");	probeLabels.push_back("SS [-]");
	}
	char head[2000];
	ostringstream oss;
	unsigned int ii = 0;

	_probeOut.clear();
	for(size_t k = 0; k < _probeList.size(); k++)
		_probeOut.push_back(std::unique_ptr<OutputTimeline>(new OutputTimeline()));
	for(ProbeList::iterator it = _probeList.begin(); it != _probeList.end(); ++it)
	{
		sprintf_s(head, "%s%02u", "Probe", ii);
		oss.clear();
		oss.str("");
		//  The position as the well-mixed probes give it (mm when the output
		//  is dimensioned), and the zeta of the cell through the film
		{
			const double sc = _inputParam->Output.Dimensioned ? _NonDim.ExtentScale() * 1e3 : 1.0;
			oss<<(_inputParam->Output.Dimensioned ? "# (x, y) [mm] = (" : "# (x, y) [-] = (")
				<< _field.GetGrid3().XOf<Cell>(it->Ix(),it->Jy(), it->Kz()) * sc <<", "
				<< _field.GetGrid3().YOf<Cell>(it->Ix(),it->Jy(), it->Kz()) * sc <<"), zeta = "
				<< _field.GetGrid3().ZOf<Cell>(it->Ix(),it->Jy(), it->Kz()) <<" ";
		}

		_probeOut[ii]->Build(_timeStep, _NonDim, _inputParam->Output.Dimensioned, JobInfo.ResultDir(), head, (int)probeLabels.size(), probeLabels, oss.str());
		ii++;
	}

}

void ConcZ_Model_2D::FreqOutput(bool atEnd)
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
	//  As on the well-mixed side: probes report thicknesses in um like the
	//  field files, History.dat keeps the raw meter scale beside its m^3.
	const double probeHeightScale = heightScale * (_inputParam->Output.Dimensioned ? 1e6 : 1.0);

	//  The velocity of the state that is written. What the field holds after
	//  a step is the velocity of that step's midpoint (the corrector's, which
	//  the concentration step advects with), half a step behind the film and
	//  the concentration written beside it. Recomputing it here changes
	//  nothing in the solution: every step computes its own from its start.
	if(OutputDue(itemField.c_str(), atEnd) || OutputDue(itemProbe.c_str(), atEnd)
	   || OutputDue(itemHistory.c_str(), atEnd) || OutputDue(itemCalc.c_str(), atEnd))
		this->_heightEq.CalVelocity();

	//Output of the field file
	std::string statePath;		// the state file to write once every output is done
	if(OutputDue(itemField.c_str(), atEnd))
	{
		statePath = C_RestartState::PathFor(_inputParam->Output.FieldData.NeutralFile
		                                    ? _outputNeutral.FileName() : _outputVTK.FileName());
		Info<<"Field file " << _timeStep.OutputItem(itemField.c_str()).Count()  << ", Time = "<<_timeStep.NowTime() << nl;

		// ASCII file
		//  A file that cannot be written stops the run here, before the
		//  state file that would go with it (StopOnOutputError).
		if( _inputParam->Output.FieldData.NeutralFile)
			if(!_outputNeutral.Write(_field)) StopOnOutputError();

		// VTK file
		if(_inputParam->Output.FieldData.VTKFile)
		{
			if(!_outputVTK.Write(_field)) StopOnOutputError();
			if(!_pvdFile.Add(_timeStep.NowTime() * timeScale)) StopOnOutputError();
		}

		MarkWritten(itemField.c_str(), atEnd);
	}

	//Output of the time-series data at the probe points
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

	//Output of the general time-series data
	if(OutputDue(itemHistory.c_str(), atEnd))
	{
		ostringstream oss;
		double vf, vr;
		_field.Accept(_fieldInfo);

		vf = _fieldInfo.FluidVol();
		vr = _fieldInfo.ResinVol();

		oss<<std::scientific<<_timeStep.NowTime() * timeScale;
		oss<<"\t"<< vf * volScale<< "\t" << vr * volScale << "\t"<< _fieldInfo.AveConc();
		oss<<"\t"<<_fieldInfo.MaxHS() * heightScale<<"\t"<<_fieldInfo.MinHS() * heightScale;
		oss<<"\t"<<_fieldInfo.MaxH() * heightScale<<"\t"<<_fieldInfo.MinH() * heightScale;
		oss<<"\t"<<_fieldInfo.MaxC()<<"\t"<<_fieldInfo.MinC();
		const double volErr   = (vf + _outflowedVolume + _evaporatedVolume) / _vol0 - 1.0;
		const double resinErr = (vr + _outflowedResin) / _resin0 - 1.0;
		oss<<"\t"<<_evaporatedVolume * volScale<<"\t"<<_outflowedVolume * volScale
		   <<"\t"<<_outflowedResin * volScale;
		oss<<"\t"<<volErr<<"\t"<<resinErr;
		if(!_outputHistory.Write(oss)) StopOnOutputError();
		MarkWritten(itemHistory.c_str(), atEnd);
	}

	//Output of the data about the calculation
	if(OutputDue(itemCalc.c_str(), atEnd))
	{
		ostringstream oss;
		_field.Accept(_fieldInfo);

		//  The time and the time step are written on the same scale as in History.dat ([s]).
		oss<<std::scientific<<_timeStep.NowTime() * timeScale;
		oss<< "\t" << _iter_heightEq<< "\t" << _resid_heightEq
			<< "\t" <<_timeStep.DT().now * timeScale
			<< "\t" <<_heightEq.ChangeRateOfVariable() / timeScale;
		if(!_outputCalcInfo.Write(oss)) StopOnOutputError();
		MarkWritten(itemCalc.c_str(), atEnd);
	}

	//  Last, so that it holds the schedule as every item has just left it.
	if(!statePath.empty() && (_inputParam->Output.FieldData.NeutralFile || _inputParam->Output.FieldData.VTKFile))
		WriteState(statePath);
}

bool ConcZ_Model_2D::RefineFromState()
{
	const C_Grid2D &G2 = _field.GetGrid2();
	const C_Grid3D &G3 = _field.GetGrid3();
	std::vector<char> in2((size_t)G2.imax() + 1, 0), in3((size_t)G3.imax() + 1, 0);
	for(int ii = G2.imin(); ii <= G2.imax(); ii++)
		in2[(size_t)ii] = (G2.CheckGhost(G2.Of<Zface>(ii)) == C_Grid2D::Inside);
	for(int ii = G3.imin(); ii <= G3.imax(); ii++)
		in3[(size_t)ii] = (G3.CheckGhost(G3.Of<Cell>(ii)) == C_Grid3D::Inside);
	//  The field files give the state to the digits they were written with
	//  (all of them in a .vts, 18 significant ones in a field*.dat, and one
	//  rounding in the conversion from the output units); more than this is
	//  a different state.
	const double tol = 1e-12;
	double dh = 0, ds = 0, dc = 0;
	const int rh = _rstate.RefineField("field.h", _field.HH(), in2, tol, dh);
	const int rs = _rstate.RefineField("field.s", _field.SS(), in2, tol, ds);
	const int rc = _rstate.RefineField("field.c", _field.CC(), in3, tol, dc);
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

void ConcZ_Model_2D::WriteState(const std::string &path)
{
	C_RestartState st;
	st.SetText("model", Name());
	st.PutClock(_timeStep, OutputItems());
	st.Set("steadyCount", _steadyCount);
	st.Set("adjustCount.height", _heightEq.AdjustTimeCount());
	st.SetVectorOf("height.old", _heightEq.GetXOld());
	std::vector<double> v(2);
	v[0] = _vol0; v[1] = _resin0;							st.Set("volume0", v);
	st.Set("evaporated", _evaporatedVolume);
	v[0] = _outflowedVolume; v[1] = _outflowedResin;		st.Set("outflowed", v);
	//  The corrector extrapolates c to the midpoint of the next step from
	//  the starts of this step and the one before.
	if(_heightEq.HasPreviousConc())
	{
		st.Set("conc.previousDt", _heightEq.PreviousDt());
		st.SetVectorOf("conc.previous", _heightEq.PreviousConc());
	}
	st.SetField("field.h", _field.HH());
	st.SetField("field.s", _field.SS());
	st.SetField("field.c", _field.CC());
	if(!st.Write(path))
	{
		Global::FatalError<<"Could not write the restart state "<<path<<": "<<System::LastIOError()<<"."<<nl;
		StopOnOutputError();
	}
}
void ConcZ_Model_2D::BuildRestart()
{

	//----------------------------------------------------
	// Reading in for a restart
	//----------------------------------------------------
	if(_inputParam->Restart.Enable)
	{
		Info<<"######## Restart enabled. Reading....   ########"<<nl;

		double coef = 1.0;
		if(_inputParam->Output.Dimensioned) coef = 1e-6 / _NonDim.HeightScale();
		//  Restart.ASCII picks the reader: a field*.dat handed to the binary
		//  one would fail with nothing on screen.
		const bool ok = _inputParam->Restart.ASCII
			? C_Restart::ReadAscii_CZ(_inputParam->Restart, _field, coef)
			: C_Restart::Read(_inputParam->Restart, _field, coef);
		if(!ok) Global::FatalError.Abort();
		//  A field*.dat holds the film at the wall faces, from which the
		//  reader recovers the ghost column next to the wall, but not the
		//  outer one: that held the initial film, and the first step of the
		//  height equation, whose stencil reaches both, started from it (the
		//  cell at the wall then stayed 1e-8 off for the rest of the run).
		//  Both columns are what the boundary condition makes them, a mirror
		//  image or the other end of a periodic pair, as the run had them.
		if(_inputParam->Restart.ASCII) FillGhostColumns();

		//  RestartTime is in seconds, like every time in the input, and the
		//  time held by _timeStep is dimensionless, as in the well-mixed
		//  model. (Taken as it stands, a restart at 0.007 s would start at
		//  0.007 T0 -- 218.7 s for concz_1d_reg -- past EndTime, and the run
		//  would stop before taking a step.)
		_timeStep.Set(_inputParam->Restart.RestartTime / _NonDim.TimeScale());
		_timeStep.SetCount(itemField.c_str(), _inputParam->Restart.FileIndex);
		_outputVTK.SetFNum(_inputParam->Restart.FileIndex);
		_outputNeutral.SetFNum(_inputParam->Restart.FileIndex);
		_pvdFile.SetFNum(_inputParam->Restart.FileIndex);

		//  What the field file does not hold (see RestartState.h)
		std::string msg;
		_continuing = _rstate.Load(_inputParam->Restart.FileName.c_str(), Name(),
		                           _inputParam->Restart.RestartTime / _NonDim.TimeScale(), msg);
		if(_continuing)
		{
			Info<<"Restart: "<<msg<<"."<<nl;
			RefineFromState();
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
			for(size_t k = 0; k < items.size(); k++)
				_timeStep.OutputItem(items[k]).RestartAt(_timeStep.NowTime());
		}

		if(_swConcField) _concEq.SetX(_field.CC());
		_heightEq.SetX(_field.HH());
		//  Every step starts from the copy of the field the height equation
		//  keeps; it still held the initial condition.
		_heightEq.TakeFieldAsStart();

		if(_continuing)
		{
			_rstate.TakeClock(_timeStep, OutputItems(), _inputParam->Numerics.TimeStepControl.Enable,
			                  itemField.c_str(), _inputParam->Restart.FileIndex);
			_steadyCount = (int)_rstate.Get("steadyCount");
			_heightEq.SetAdjustTimeCount((int)_rstate.Get("adjustCount.height"));
			Vector<double> hold(_heightEq.GetXOld().size());
			if(_rstate.GetVectorInto("height.old", hold)) _heightEq.SetXOld(hold);
			Vector<double> cprev(_concEq.GetX().size());
			if(_rstate.Has("conc.previousDt") && _rstate.GetVectorInto("conc.previous", cprev))
				_heightEq.SetPreviousConc(cprev, _rstate.Get("conc.previousDt"));
			//  The balance goes on against the initial state of the whole run
			{
				Field_CZ &f = _field;
				FieldInfo fi;
				f.Accept(fi);
				_vol0   = _rstate.Has("volume0") ? _rstate.Get("volume0", 0) : fi.FluidVol();
				_resin0 = _rstate.Has("volume0") ? _rstate.Get("volume0", 1) : fi.ResinVol();
				_evaporatedVolume = _rstate.Has("evaporated") ? _rstate.Get("evaporated") : 0.0;
				_outflowedVolume  = _rstate.Has("outflowed") ? _rstate.Get("outflowed", 0) : 0.0;
				_outflowedResin   = _rstate.Has("outflowed") ? _rstate.Get("outflowed", 1) : 0.0;
			}
			//  The probes, History, CalcInfo and the .pvd go on in the files
			//  of the run being continued, if this one writes into them.
			const double ts = _inputParam->Output.Dimensioned ? _NonDim.TimeScale() : 1.0;
			const double tmax = _timeStep.NowTime() * ts * (1.0 + 2e-6);
			if(!_outputHistory.Resume(tmax)) StopOnOutputError();
			if(!_outputCalcInfo.Resume(tmax)) StopOnOutputError();
			for(size_t k = 0; k < _probeOut.size(); k++) if(!_probeOut[k]->Resume(tmax)) StopOnOutputError();
			if(!_pvdFile.Resume(tmax)) StopOnOutputError();
			_rstate.Clear();
		}

		Info<<"ASCII:       "<<(_inputParam->Restart.ASCII ? "true" : "false")<<nl;
		Info<<"File Name:   "<<_inputParam->Restart.FileName<<nl;
		Info<<"RestartTime: "<<_inputParam->Restart.RestartTime<<nl;
		Info<<"FileIndex:   "<<_inputParam->Restart.FileIndex<<nl;
		Info<<"######## Restart setting applied.   ########"<<nl;
	}
}

void ConcZ_Model_2D::Run()
{
	// Match the maximum time step to the smallest output interval.
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
		// The case where dt is larger than dtmax.
		Global::Info<<"initial DT is larger than Max-DT. DT is set to Max-DT ("<<_timeStep.MaxDT()<<")."<<nl;
		_timeStep.ChangeDT(_timeStep.MaxDT());
	}

	BuildRestart();

	//Compute the initial velocity field.
	this->_heightEq.CalVelocity();

	//Probe the initial values.
	_field.Accept(_fieldInfo);
	for(ProbeList::iterator it = _probeList.begin(); it != _probeList.end(); ++it)
		_field.Accept(*it);

	_fieldInfoIni = _fieldInfo;
	//  (A run continued from a state file took them from it, BuildRestart.)
	if(!_continuing)
	{
		_vol0   = _fieldInfoIni.FluidVol();
		_resin0 = _fieldInfoIni.ResinVol();
		_evaporatedVolume = _outflowedVolume = _outflowedResin = 0;
	}

	//  The Courant number of a limited concentration scheme: watched and
	//  reported, never acted on (CourantCheck.h). Its rate is that of the
	//  volume fluxes of the step, so the first step has none before it.
	_courantConc.Build("Numerics.Convection.Concentration", _inputParam->Numerics.Convection.Concentration,
	                   "the concentration", _concEq.Scheme().IsLimited());

	//Time-step loop
	const double loopStart = Perf::NowSeconds();		// for measurement only; it plays no part in the numerics
	_loopStartTime = _timeStep.NowTime();
	_sumStepDt = 0;
	//  A run continued from the step at which the steady-state test was met
	//  -- it was stopped on EndTime or Numerics.MaxSteps at that same step --
	//  takes no further step, as the uninterrupted run took none.
	const bool steadyAlready = _continuing && _inputParam->Numerics.SteadyStateCriterion > 0
	                        && _steadyCount > _inputParam->Numerics.SteadyStateSteps;
	if(steadyAlready) Info<<"\n############## Steady state (when restarted). #############"<<nl;
	//  The same for the dry-out test (see WellMixedModelBase.h). _fieldInfo
	//  was taken of the field just above.
	const bool dryAlready = _continuing && !steadyAlready && DriedOut();
	if(dryAlready) Info<<"\n############## Dry out (when restarted). #############"<<nl;
	while(!steadyAlready && !dryAlready && _timeStep.IsRunning())
	{
		{ YUZEN_PROF("output"); FreqOutput(); }
		_timeStep.ApplyListedTimes();

		Info<<".";
		//  The concentration at the start of the step, for the solute the
		//  sideways outflow of 1-D carries (ColumnSinks)
		if(_heightEq.HasSideOutflow()) _cStart = _field.CC();
		{ YUZEN_PROF("step"); _stepStatus = _heightEq.TimeAdvance(_iter_heightEq, _resid_heightEq); }
		++_nSteps;
		_sumStepDt += _heightEq.StepDt();
		//  What left the columns other than through a face, as the step took it
		{
			C_Grid3D *G3 = _field.GetGridPtr3();
			double ev = 0, sideR = 0;
			//  The cell measure of FieldInfo: dx dy, dx alone in 1-D
			const double cm = G3->dx() * (_field.GetGridPtr2()->ny() > 1 ? G3->dy() : 1.0);
			_outflowedVolume  += _heightEq.ColumnSinks(cm, _field,
			                        _heightEq.HasSideOutflow() ? _cStart : _field.CC(), ev, &sideR);
			_outflowedResin   += sideR;
			_evaporatedVolume += ev;
			//  The solute moved into cells below c = 0, and made up there;
			//  a make-up beyond rounding is reported as it happens, once
			const double dm = _heightEq.DeficitMoved(cm), add = _heightEq.DeficitAdded(cm);
			if(dm > 0.0) { _deficitMoved += dm; ++_deficitSteps; }
			_deficitAdded += add;
			if(add > 1e-12 * std::fabs(_resin0) && !_warnedAdd)
			{
				std::ostringstream m;
				m<<std::setprecision(3)<<"The concentration fell below 0 at step "<<_nSteps
				 <<" with less solute than the film needs to fill it: "<<add / _resin0
				 <<" of the initial solute was made up at c = 0 and is an error of ResinVolErr."
				 <<" Reported once; the total is given at the end of the run.";
				Global::Warning<<m.str()<<nl;
				_warnedAdd = true;
			}
		}
		//  The film-height equation of ConcZ is semi-implicit, so there is no Newton
		//  iteration. The iter that TimeAdvance returns is the number of linear
		//  iterations as it stands.
		_sumBCGIter += _iter_heightEq;
		_sumCoupling += _heightEq.CouplingPasses();

		_timeStep++;

		if(_courantConc.On())
			_courantConc.Step(_concEq.CourantRate(), _heightEq.StepDt(), _nSteps, _timeStep.NowTime() * T0);

		//Check the termination conditions.
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
	//  state at EndTime is kept only when the output schedule happens to land
	//  on it exactly. A diverged run is the exception: what it holds is not a
	//  solution and writing it would put nonsense in the result files.
	if(!_diverged) FreqOutput(true);

	//  Record in JobInfo.txt what was measured for the time-step loop alone.
	JobInfo<<nl<<nl<<"###### Performance ######"<<nl;
	JobInfo.Add("TimeLoopSeconds", _loopSeconds,     "Time loop [s]");
	JobInfo.Add("TimeStepCount",   (int)_nSteps,     "Time steps");
	JobInfo.Add("NKIterTotal",     (int)_sumNKIter,  "Newton iterations");	// always 0, since the scheme is semi-implicit
	JobInfo.Add("BCGIterTotal",    (int)_sumBCGIter, "Linear iterations");
	JobInfo.Add("CouplingPassTotal", (int)_sumCoupling, "Film-concentration passes taken again");
	JobInfo.Add("CouplingStepCuts", (int)_heightEq.CouplingCuts(), "Steps failed: the passes did not settle");
	//  The solute moved into cells below c = 0, and made up there, as
	//  fractions of the initial solute (an unbounded convection scheme
	//  undershooting at a sharp front)
	{
		const double r0 = (_resin0 != 0.0) ? std::fabs(_resin0) : 1.0;
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
	_courantConc.Report("conc");
	//  The sum of the steps solved and the times the clock went from and to,
	//  on the time scale of the outputs: they agree to rounding
	//  (tools/features.py, adaptive_clock).
	{
		const double ts = _inputParam->Output.Dimensioned ? _NonDim.TimeScale() : 1.0;
		auto full = [](double v) { std::ostringstream o; o<<std::setprecision(17)<<v; return o.str(); };
		JobInfo.Add("StepDtSum",     full(_sumStepDt * ts),          "Sum of the steps solved");
		JobInfo.Add("LoopStartTime", full(_loopStartTime * ts),      "Time at the start of the loop");
		JobInfo.Add("LoopEndTime",   full(_timeStep.NowTime() * ts), "Time at the end of the loop");
	}

	//  If the loop was left because of divergence, this must not claim a
	//  "normal termination". The exit code is separated too, so that a batch run
	//  or a test can tell that it failed.
	if(_diverged)
		Info<<"\n############## Program terminated after divergence. #############"<<nl;
	else if(Interrupt::Requested())
		Info<<"\n############## Program stopped by "<<Interrupt::Name()<<"; its output is complete. #############"<<nl;
	else
		Info<<"\n############## Program normally terminated. #############"<<nl;
}

bool ConcZ_Model_2D::DriedOut() const
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

bool ConcZ_Model_2D::CheckTermination()
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
	//  Detection of divergence.
	//
	//  *** This test must come before the dry-out test and the steady-state test, ***
	//  *** and before the time test. ***
	//  After the time test, a field that broke down in the last step of a run
	//  would end it as "Step count reached" with exit code 0, and the final
	//  output would write the broken field (tools/proptest.py).
	//  The other two compare a value against a threshold, so once the field has
	//  been polluted with NaN the comparison itself loses its meaning and a
	//  diverged calculation terminates normally as a "Dry out" or a
	//  "Steady state" (at 240x41, dt=1e-3 a run really did end with exit code 0
	//   while ResinVol = 5.6e+21 and minC = -1.9e+27).
	//
	//  Two things are looked at.
	//    (1) Whether the time integration itself returned a failure (the
	//        return value of TimeAdvance).
	//    (2) Whether the field is finite and within a physically possible range.
	//        The concentration should be in [0,1] and the film thickness should
	//        be positive, so straying far outside that is taken as divergence.
	//        The thresholds are loose (a negative concentration is allowed down
	//        to -0.1, and so on) because a healthy calculation can stray slightly
	//        through discretization error; what is to be caught here is an
	//        excursion of an order of magnitude.
	//------------------------------------------------------------------
	if(_stepStatus != 0)
	{
		Info<<"\n############## Diverged: the time step did not converge. #############"<<nl;
		Global::Warning<<"ConcZ: time step "<<_nSteps<<" failed (ret = "<<_stepStatus
		               <<", iter = "<<_iter_heightEq<<", resid = "<<_resid_heightEq<<")."<<nl;
		_diverged = true;
		return true;
	}

	{
		//  *** The aggregation must be taken again here. ***
		//  The only thing that updates _fieldInfo is FreqOutput(), and that runs
		//  only at an output time. Reading it without taking it again would mean
		//  looking at the stale values from the last output, and a field that broke
		//  in between would go unnoticed (at 60x11, dt=1e-2, FluidVol swinging from
		//   1e-8 to negative was slipping through). The dry-out test below uses the
		//  same values, so this makes that one look at the current field too.
		_field.Accept(_fieldInfo);

		const double vol  = _fieldInfo.FluidVol();
		const double maxC = _fieldInfo.MaxC();
		//  The lowest c of the start of the step as well: a field file read
		//  with a cell far below 0 must stop the run at its first step, whose
		//  fill (RedistributeDeficit) may have brought the cell back to 0
		const double minC = std::min(_fieldInfo.MinC(), _heightEq.MinStartConc());
		const double minH = _fieldInfo.MinH();
		const char *why = 0;

		//  A NaN is not equal to itself, which is what this form picks up on.
		if(!(vol == vol) || !(maxC == maxC) || !(minC == minC) || !(minH == minH))
			why = "the field is not a number";
		else if(maxC > 1.0 + kConcSlack || minC < -kConcSlack)
			why = "the concentration left [0,1] by more than the allowed slack";
		else if(minH < -kHeightSlack)
			why = "the film thickness went negative";
		//  There is also a way of breaking in which the total volume of liquid goes
		//  negative. At 60x11, dt=1e-2 the concentration stayed within [0,1] while
		//  FluidVol went from 1e-8 to -3.8e-10, slipping through all three tests
		//  above. Decreasing through evaporation is normal, so only a swing into
		//  negative values is looked at.
		else if(vol < 0.0)
			why = "the fluid volume went negative";

		if(why)
		{
			Info<<"\n############## Diverged: "<<why<<". #############"<<nl;
			Global::Warning<<"ConcZ: diverged at step "<<_nSteps<<": "<<why
			               <<" (FluidVol = "<<vol<<", C = ["<<minC<<", "<<maxC
			               <<"], minH = "<<minH<<")."<<nl;
			_diverged = true;
			return true;
		}
	}

	if(_timeStep.IsRunning() == false)
	{
		//  Report separately whether the run ended on the time or was cut short by
		//  Numerics.MaxSteps, so that the reason is not mistaken.
		if(_timeStep.TimeStep() > _timeStep.EndTimeStep())
			Info<<"\n############## Step count reached Numerics.MaxSteps. #############"<<nl;
		else
			Info<<"\n############## Time reached to the EndTime. #############"<<nl;
		return true;
	}

	//  Dry-out test. The threshold and what is measured follow EndConcType and
	//  EndConcValue in the input.
	if(DriedOut())
	{
		Info<<"\n############## Dry out. #############"<<nl;
		return true;
	}

	//  The steady-state test asks how far the solution would move over all the
	//  remaining time at its present rate (see the comment in
	//  WellMixedModelBase::CheckTermination).
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



//  Fill the two ghost layers of h and s on every edge from the real cells:
//  Periodic takes the cells at the other end, anything else (the CZ model
//  accepts Symmetry and Periodic only) mirrors about the wall face. x first
//  on every row, then y on every column, so the corners come out mirrored
//  (or wrapped) in both.
void ConcZ_Model_2D::FillGhostColumns()
{
	C_Grid2D *G = _field.GetGridPtr2();
	const int nx = G->nx(), ny = G->ny();
	const bool perX = (_boundary.TypeXmin() == BoundaryCondition::Periodic);
	const bool perY = (_boundary.TypeYmin() == BoundaryCondition::Periodic);

	auto copy = [&](int ix, int jy, int fx, int fy)
	{
		_field.HH(G->At<Zface>(ix, jy)) = _field.HH(G->At<Zface>(fx, fy));
		_field.SS(G->At<Zface>(ix, jy)) = _field.SS(G->At<Zface>(fx, fy));
	};
	for(int jy = 0; jy < ny; jy++)
	{
		if(perX) { copy(0, jy, nx-4, jy); copy(1, jy, nx-3, jy); copy(nx-2, jy, 2, jy); copy(nx-1, jy, 3, jy); }
		else     { copy(0, jy, 3,    jy); copy(1, jy, 2,    jy); copy(nx-2, jy, nx-3, jy); copy(nx-1, jy, nx-4, jy); }
	}
	if(ny == 1) return;
	for(int ix = 0; ix < nx; ix++)
	{
		if(perY) { copy(ix, 0, ix, ny-4); copy(ix, 1, ix, ny-3); copy(ix, ny-2, ix, 2); copy(ix, ny-1, ix, 3); }
		else     { copy(ix, 0, ix, 3);    copy(ix, 1, ix, 2);    copy(ix, ny-2, ix, ny-3); copy(ix, ny-1, ix, ny-4); }
	}
}
