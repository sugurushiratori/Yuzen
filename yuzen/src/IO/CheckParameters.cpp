/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CheckParameters.cpp : checks that the combination of input values is not
*                        self-contradictory
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;

#include "YuzenUtil.h"
using namespace yuzen::util;
#include "CheckParameters.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

enum tPosi{	Xmin, Xmax, Ymin, Ymax, Other};

tPosi GetPosi(std::string str)
{
	if(str.compare("Xmin") == 0) return Xmin;
	if(str.compare("Xmax") == 0) return Xmax;
	if(str.compare("Ymin") == 0) return Ymin;
	if(str.compare("Ymax") == 0) return Ymax;
	return Other;
}

void Message(std::string param1, std::string param2)
{
	Sout<<"Input mismatch: "<<param1<<" and "<<param2<<nl;
}

//------------------------------------------------------------------
//  Values with which the equations mean nothing.
//
//  Taken as given, each of these would divide by zero or set up a problem
//  without a meaning: a viscosity of zero or below makes the time scale
//  zero or negative and the run end after no step with exit code 0; a time
//  step, an output interval or a step-control factor of zero stops time
//  from advancing; a negative surface tension makes the capillary length
//  the root of a negative number; a film of zero thickness everywhere has
//  no liquid in it, and History.dat fills with NaN (tools/proptest.py and
//  yuzen/tests/fuzz feed such edge values). They are refused here, with
//  the key and the value. What is refused is what the
//  manual's key tables already called positive, non-zero or non-negative,
//  and nothing that has a physical meaning however unusual (a negative
//  gravity, a condensing film, an exponent of any sign).
//------------------------------------------------------------------
namespace {

struct Refusals
{
	bool ok = true;
	void Refuse(const std::string &key, double value, const char *why)
	{
		Sout<<"Input mismatch: "<<key<<" = "<<value<<": "<<why<<nl;
		ok = false;
	}
	void Positive(const std::string &key, double v, const char *what = "must be positive")
	{
		if(!(v > 0)) Refuse(key, v, what);
	}
	void NonNegative(const std::string &key, double v)
	{
		if(!(v >= 0)) Refuse(key, v, "must not be negative");
	}
	void NonZero(const std::string &key, double v)
	{
		if(v == 0) Refuse(key, v, "must not be zero: the shape divides by it");
	}
};

//  The shapes of Initial.SurfaceShape and of Substrate have the same members
//  but are different types.
template<class S> void CheckShapes(const S &s, const std::string &where, Refusals &r)
{
	for(size_t i = 0; i < s.RectBump.size(); i++)
		for(int k = 0; k < 2; k++)
			r.Positive(where + ".RectBump[" + std::to_string(i) + "].Width[" + std::to_string(k) + "]",
			           s.RectBump[i].Width[k]);
	for(size_t i = 0; i < s.EllipseBump.size(); i++)
		for(int k = 0; k < 2; k++)
			r.NonZero(where + ".EllipseBump[" + std::to_string(i) + "].Width[" + std::to_string(k) + "]",
			          s.EllipseBump[i].Width[k]);
	for(size_t i = 0; i < s.PeriodicWave.size(); i++)
		for(int k = 0; k < 2; k++)
			r.NonZero(where + ".PeriodicWave[" + std::to_string(i) + "].Wavelength[" + std::to_string(k) + "]",
			          s.PeriodicWave[i].Wavelength[k]);
	for(size_t i = 0; i < s.SmoothConvex.size(); i++)
		for(int k = 0; k < 2; k++)
			r.NonZero(where + ".SmoothConvex[" + std::to_string(i) + "].Width[" + std::to_string(k) + "]",
			          s.SmoothConvex[i].Width[k]);
}

//  c -> base exp(e (c - c0)) (viscosity, diffusivity) at both ends of [0, 1].
double Exponential(double base, double e, double c0, double c) { return base*std::exp(e*(c - c0)); }

}  // namespace

static bool CheckMeaningful(Input &inp)
{
	Refusals r;
	const Input::I_Fluid &f = inp.Fluid;

	r.Positive("Fluid.Density.Value", f.Density.Value);
	//  k (1 - r/R')^(-n) with k = (1 - n)(2 - n)/2 has a mean of 1 over the
	//  pool only for n < 1; at n >= 1 the rate is zero or negative inside it.
	if(inp.Physics.Evaporation.Enable && inp.Physics.Evaporation.SpatialFunc == "negPow"
	   && !(inp.Physics.Evaporation.NegPowExponent < 1.0))
		r.Refuse("Physics.Evaporation.NegPowExponent", inp.Physics.Evaporation.NegPowExponent,
		         "must be less than 1");
	r.Positive("Fluid.SurfaceTension.BaseValue", f.SurfaceTension.BaseValue);
	//  The linear sigma(c) may go below zero somewhere in [0, 1] (the stress
	//  cases take Coef = 0.1 at a base of 0.03): it is an approximation about
	//  BaseConc, and only the concentrations a run reaches matter.
	r.Positive("Fluid.Viscosity.BaseValue", f.Viscosity.BaseValue);
	for(int c = 0; c <= 1; c++)
		if(f.Viscosity.BaseValue > 0
		   && !std::isfinite(Exponential(f.Viscosity.BaseValue, f.Viscosity.Exponent, f.Viscosity.BaseConc, c)))
			r.Refuse("Fluid.Viscosity.Exponent", f.Viscosity.Exponent,
			         "makes the viscosity overflow within c in [0, 1]");
	//  D0 = 0 is "no diffusion" to the well-mixed model; ConcZ divides by it
	//  (the Peclet number of its evaporation condition).
	if(inp.Physics.ConcField.Model.compare("CZ") == 0)
		r.Positive("Fluid.Diffusivity.BaseValue", f.Diffusivity.BaseValue,
		           "must be positive with Physics.ConcField.Model = \"CZ\"");
	else
		r.NonNegative("Fluid.Diffusivity.BaseValue", f.Diffusivity.BaseValue);
	for(int c = 0; c <= 1; c++)
		if(!std::isfinite(Exponential(f.Diffusivity.BaseValue, f.Diffusivity.Exponent, f.Diffusivity.BaseConc, c)))
			r.Refuse("Fluid.Diffusivity.Exponent", f.Diffusivity.Exponent,
			         "makes the diffusivity overflow within c in [0, 1]");
	r.Positive("Fluid.MolarMass.Solvent", f.MolarMass.Solvent);
	r.Positive("Fluid.MolarMass.Resin", f.MolarMass.Resin);

	//  The scales
	r.Positive("Initial.SurfaceShape.TotalHeight", inp.Initial.SurfaceShape.TotalHeight,
	           "must be positive: it is the height scale h0");
	const double hs = inp.Initial.SurfaceShape.PrecursorRatio;
	if(!(hs >= 0 && hs < 1))
		r.Refuse("Initial.SurfaceShape.PrecursorRatio", hs,
		         "must be in [0, 1): it is the precursor film as a fraction of TotalHeight");
	if(inp.LengthScaleType.compare("custom") == 0)
		r.Positive("LengthScale", inp.LengthScale, "must be positive with LengthScaleType = \"custom\"");

	//  The domain
	if(!std::isfinite(inp.Grid.Xmax - inp.Grid.Xmin))
		r.Refuse("Grid.Xmax - Grid.Xmin", inp.Grid.Xmax - inp.Grid.Xmin, "must be finite");
	if(inp.Grid.Dimension.compare("xyh") == 0)
	{
		if(!(inp.Grid.Ymax > inp.Grid.Ymin))
			r.Refuse("Grid.Ymax", inp.Grid.Ymax, "must exceed Grid.Ymin in a 2-D film");
		else if(!std::isfinite(inp.Grid.Ymax - inp.Grid.Ymin))
			r.Refuse("Grid.Ymax - Grid.Ymin", inp.Grid.Ymax - inp.Grid.Ymin, "must be finite");
	}

	CheckShapes(inp.Initial.SurfaceShape, "Initial.SurfaceShape", r);
	CheckShapes(inp.Substrate, "Substrate", r);
	for(size_t i = 0; i < inp.Initial.SurfaceShape.EdgeBump.size(); i++)
	{
		const std::string e = "Initial.SurfaceShape.EdgeBump[" + std::to_string(i) + "].";
		r.Positive(e + "Alpha", inp.Initial.SurfaceShape.EdgeBump[i].alpha);
		r.Positive(e + "InnerLength", inp.Initial.SurfaceShape.EdgeBump[i].Lin);
	}

	//  Time: a step, an interval or a factor of zero stops it.
	const Input::I_Numerics &n = inp.Numerics;
	r.Positive("Numerics.InitialDt", n.InitialDt);
	r.NonNegative("Numerics.EndTime", n.EndTime);
	r.Positive("Numerics.TimeStepControl.MaxDtFactor", n.TimeStepControl.MaxDtFactor);
	r.Positive("Numerics.TimeStepControl.MinDt", n.TimeStepControl.MinDt);
	r.Positive("Numerics.TimeStepControl.Criterion", n.TimeStepControl.Criterion);
	r.Positive("Numerics.HeightEquation.NewtonEpsilon", n.HeightEquation.NewtonEpsilon);
	if(n.HeightEquation.NewtonMaxIterations < 1)
		r.Refuse("Numerics.HeightEquation.NewtonMaxIterations", n.HeightEquation.NewtonMaxIterations, "must be at least 1");
	r.Positive("Numerics.HeightEquation.LinearSolver.Epsilon", n.HeightEquation.LinearSolver.Epsilon);
	r.Positive("Numerics.ConcEquation.LinearSolver.Epsilon", n.ConcEquation.LinearSolver.Epsilon);
	r.Positive("Output.FieldData.TimeInterval", inp.Output.FieldData.TimeInterval);
	r.Positive("Output.GeneralHistory.TimeInterval", inp.Output.GeneralHistory.TimeInterval);
	r.Positive("Output.CalcInfo.TimeInterval", inp.Output.CalcInfo.TimeInterval);
	r.Positive("Output.ProbeList.TimeInterval", inp.Output.ProbeList.TimeInterval);

	//  Field output at listed times (Times, or a LogSchedule resolved into
	//  them by the reader): instead of TimeInterval, ascending, positive;
	//  those after EndTime are dropped, and at least one must remain
	{
		Input::I_Output::I_FieldData &fd = inp.Output.FieldData;
		const std::string key = fd.LogSchedule.Present ? "Output.FieldData.LogSchedule" : "Output.FieldData.Times";
		if(fd.Listed())
		{
			if(fd.HasTimeInterval)
				r.Refuse("Output.FieldData.TimeInterval", fd.TimeInterval,
				         "cannot be given with Output.FieldData.Times or LogSchedule: give one of them");
			for(size_t i = 0; i < fd.Times.size(); i++)
			{
				if(!(fd.Times[i] > 0))
				{
					r.Refuse(key + "[" + std::to_string(i) + "]", fd.Times[i], "must be positive (the state at 0 is always written)");
					break;
				}
				if(i > 0 && !(fd.Times[i] > fd.Times[i-1]))
				{
					r.Refuse(key + "[" + std::to_string(i) + "]", fd.Times[i], "must be after the time before it: the times must be in ascending order");
					break;
				}
			}
			if(r.ok)
			{
				std::vector<double> kept;
				for(size_t i = 0; i < fd.Times.size(); i++)
					if(fd.Times[i] <= n.EndTime) kept.push_back(fd.Times[i]);
				if(kept.empty())
					r.Refuse(key + "[0]", fd.Times[0], "every listed time is after Numerics.EndTime: no field would be written on the list");
				else
					fd.Times.swap(kept);
			}
		}
		if(fd.WritePotential && inp.Physics.ConcField.Model == "CZ")
			r.Refuse("Output.FieldData.WritePotential", 1,
			         "is not implemented for Physics.ConcField.Model = \"CZ\"");
		if(inp.Output.Energy.Enable && inp.Physics.ConcField.Model == "CZ")
			r.Refuse("Output.Energy.Enable", 1,
			         "is not implemented for Physics.ConcField.Model = \"CZ\"");
		if(!fd.Listed() && fd.ExactTimes)
			r.Refuse("Output.FieldData.ExactTimes", 1, "applies to Output.FieldData.Times only, which is not given");
	}

	//  Navier slip (Mobility.h). CZ solves the velocity profile across the
	//  film, which would need the uniform slip velocity added to it: not
	//  implemented, so refused by name rather than ignored.
	r.NonNegative("Physics.Slip.Length", inp.Physics.Slip.Length);
	if(inp.Physics.Slip.Enable && inp.Physics.ConcField.Model.compare("CZ") == 0)
		r.Refuse("Physics.Slip.Enable", 1,
		         "Physics.Slip is not implemented for Physics.ConcField.Model = \"CZ\"");
	if(inp.Physics.Slip.Enable && inp.Physics.Slip.Length == 0)
		Global::Warning<<"Physics.Slip.Enable = true with Physics.Slip.Length = 0: no slip, the run is the same as without the block."<<nl;
	if(inp.Physics.Source.Enable)
	{
		r.NonNegative("Physics.Source.TimeConstant", inp.Physics.Source.TimeConstant);
		//  The liquid supplied at Conc is mixed into the film by it: a mass
		//  fraction, as the film's own concentration is
		const double cs = inp.Physics.Source.Conc;
		if(inp.Physics.Source.ConcType == "supplied" && !(cs >= 0 && cs <= 1))
			r.Refuse("Physics.Source.Conc", cs,
			         "must be in [0, 1] with Physics.Source.ConcType = \"supplied\": it is the mass fraction of the supplied liquid");
	}
	if(inp.Restart.Enable)
	{
		r.NonNegative("Restart.RestartTime", inp.Restart.RestartTime);
		if(inp.Restart.FileIndex < 0)
			r.Refuse("Restart.FileIndex", inp.Restart.FileIndex, "must not be negative");
	}
	return r.ok;
}

bool CheckParameters::Check(Input& inp)
{
	bool ch = true;
	// Dimension
	if(inp.Grid.Dimension.compare("xyh") == 0)
	{
		if(inp.Grid.Ny < 2) Message("Grid.Dimension", "Grid.Ny"), ch = false;
	}

	if(inp.Physics.ConcField.Model.compare("CZ") == 0)
	{
		if(inp.Grid.Nz < 2) Message("Physics.ConcField.Model", "Grid.Nz"), ch = false;
	}

	//  The edges that exist: x always, y only for a 2-D film.
	const bool twoD = (inp.Grid.Dimension.compare("xyh") == 0);
	struct EdgeBC { const char *key; const std::string *type; };
	const EdgeBC edges[4] = {
		{"Boundary.XminBC.Type", &inp.Boundary.XminBC.Type},
		{"Boundary.XmaxBC.Type", &inp.Boundary.XmaxBC.Type},
		{"Boundary.YminBC.Type", &inp.Boundary.YminBC.Type},
		{"Boundary.YmaxBC.Type", &inp.Boundary.YmaxBC.Type}};
	const int nEdges = twoD ? 4 : 2;

	//  Periodic pairs one edge with the one opposite; on one edge alone it has
	//  nothing to pair with. It is checked here for both axes and both
	//  dimensions: a one-sided Periodic that got through would run to the end
	//  and give a wrong answer without a word.
	for(int k = 0; k < nEdges; k += 2)
	{
		const bool a = edges[k].type->compare("Periodic") == 0;
		const bool b = edges[k + 1].type->compare("Periodic") == 0;
		if(a != b)
		{
			Sout<<"Input mismatch: "<<edges[k].key<<" and "<<edges[k + 1].key
			    <<": Periodic must be set on both edges of a pair or on neither."<<nl;
			ch = false;
		}
	}

	//  The centrifugal force pushes the film along the arm x - x0 (and y -
	//  y0 in 2-D), which grows across the domain and jumps by its length
	//  across a periodic seam: a rotating film is not periodic, and the seam
	//  takes in liquid on one side that it does not give out on the other
	//  (VolErr drifts). Refused, except under a manufactured solution, whose
	//  source is formed for the arm of each cell (tools/mms.py, wm_all_*,
	//  cz_all_*: a periodic domain with every term; the seam is then a
	//  consistent discretization of the manufactured problem). A strip of a
	//  film wider than the domain is an Outlet at each end (tools/analytic.py
	//  spin_1d).
	if(inp.Physics.Centrifugal.Enable && inp.Physics.Centrifugal.RotationRate != 0
	   && !inp.Numerics.ManufacturedSolution.Enable)
	{
		for(int k = 0; k < nEdges; ++k)
		{
			if(edges[k].type->compare("Periodic") != 0) continue;
			Sout<<"Input mismatch: Physics.Centrifugal.Enable = true and "<<edges[k].key
			    <<" = \"Periodic\": the centrifugal force grows with the distance from the"
			    <<" rotation axis, so a rotating film is not periodic (across the seam the"
			    <<" arm jumps by the length of the domain, and liquid is created or lost"
			    <<" there). Use Symmetry or Pinning for a wall, Outlet for a film that"
			    <<" goes on beyond the edge."<<nl;
			ch = false;
			break;
		}
	}

	//  The CZ model builds its boundary rows for two kinds only: tied to the
	//  opposite side when periodic, mirrored otherwise. Pinning, FreshInlet and
	//  Outlet were therefore taken as Symmetry without a word. Refuse them.
	if(inp.Physics.ConcField.Model.compare("CZ") == 0)
	{
		for(int k = 0; k < nEdges; ++k)
		{
			const std::string &t = *edges[k].type;
			if(t.compare("Symmetry") != 0 && t.compare("Periodic") != 0)
			{
				Sout<<"Input mismatch: Physics.ConcField.Model = \"CZ\" and "<<edges[k].key
				    <<" = \""<<t<<"\": the CZ model implements only Symmetry and Periodic"
				    <<" boundaries, and would treat this one as Symmetry."<<nl;
				ch = false;
			}
		}
	}

	//  The film-height equation of CZ has no upwinded flux (its body-force
	//  and Marangoni fluxes are built from the face mobility, centered), so a
	//  scheme for it would be read by nobody.
	if(inp.Physics.ConcField.Model.compare("CZ") == 0 && inp.Numerics.Convection.HeightGiven)
	{
		Sout<<"Input mismatch: Physics.ConcField.Model = \"CZ\" and Numerics.Convection.Height = \""
		    <<inp.Numerics.Convection.Height<<"\": the CZ film-height equation has no upwinded"
		    <<" convective term for the scheme to act on. Remove the key (Numerics.Convection."
		    <<"Concentration does apply to CZ)."<<nl;
		ch = false;
	}

	if(inp.Initial.SurfaceShape.EdgeBump.size() !=0)
	{
		std::string str;
		for(size_t ii = 0; ii < inp.Initial.SurfaceShape.EdgeBump.size(); ii++)
		{
			// str is set here, outside the switch: a statement placed before
			//          the first case of a switch is never executed. The index
			//          goes through std::to_string, since "...EdgeBump[" + ii
			//          would be pointer arithmetic into the string literal, not
			//          concatenation.
			str = "Initial.SurfaceShape.EdgeBump[" + std::to_string(ii) + "].Location";
			switch(GetPosi(inp.Initial.SurfaceShape.EdgeBump[ii].Location))
			{
				case Xmin:
					if(inp.Boundary.XminBC.Type.compare("Pinning") != 0) Message(str, "Boundary.XminBC.Type"), ch = false;
					break;
				case Xmax:
					if(inp.Boundary.XmaxBC.Type.compare("Pinning") != 0) Message(str, "Boundary.XmaxBC.Type"), ch = false;
					break;
				case Ymin:
					if(inp.Boundary.YminBC.Type.compare("Pinning") != 0) Message(str, "Boundary.YminBC.Type"), ch = false;
					break;
				case Ymax:
					if(inp.Boundary.YmaxBC.Type.compare("Pinning") != 0) Message(str, "Boundary.YmaxBC.Type"), ch = false;
					break;
				default:
					break;
			}
		}
	}

	if(!CheckMeaningful(inp)) ch = false;

	return ch;
}
