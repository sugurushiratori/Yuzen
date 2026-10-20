/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  fuzz_restart.cpp : libFuzzer harness of the restart readers. A restart
*  file is an input like any other: it may be truncated, of another grid, of
*  another model, edited by hand or not a restart file at all.
*
*  The first byte picks the reader, the rest is the file:
*    0  field*.dat of WM-1D      (C_Restart::ReadAscii_WM1D)
*    1  field*.dat of WM-2D      (C_Restart::ReadAscii_WM2D)
*    2  .vts of WM-1D            (C_Restart::Read, Field_WM1D)
*    3  .vts of WM-2D            (C_Restart::Read, Field_WM2D)
*    4  .vts of CZ, 1-D film     (C_Restart::Read, Field_CZ)
*    5  .vts of CZ, 2-D film
*    6  field*.dat of CZ, 1-D    (C_Restart::ReadAscii_CZ)
*    7  field*.dat of CZ, 2-D
*    8  field*.state             (C_RestartState::Load, then what the models
*                                 take from it: TakeClock and the counters)
*  The fields are built once, from small grids, as the models build them.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "fuzz_common.h"

#include "Input.h"
#include "LubFields.h"
#include "Restart.h"
#include "RestartState.h"
#include "TimeSteps.h"

#include <cmath>
#include <fstream>
#include <vector>

using namespace yuzen;

namespace {

const char *kTemplate = R"(
SolverName = "Yuzen-1.0"
JobName = "fuzz"
[Grid]
Dimension = "@DIM@"
Xmin = 0
Xmax = 1e-3
Ymin = 0
Ymax = 1e-3
Nx = @NX@
Ny = @NY@
Nz = @NZ@
[Initial.Concentration]
Type = "uniform"
Value = 0.5
[Initial.SurfaceShape]
TotalHeight = 10e-6
[Boundary.XminBC]
[Boundary.XmaxBC]
[Fluid]
Name = "f"
[Fluid.Density]
Value = 1000.0
[Fluid.SurfaceTension]
BaseValue = 30.0e-3
[Fluid.Viscosity]
BaseValue = 100.0e-3
[Fluid.Diffusivity]
BaseValue = 1.0e-9
[Fluid.MolarMass]
[Numerics.ConcEquation.LinearSolver]
[Numerics.HeightEquation.LinearSolver]
[Numerics.TimeStepControl]
[Output.CalcInfo]
[Output.GeneralHistory]
[Physics.ConcField]
Model = "@MODEL@"
[Numerics]
InitialDt = 1e-3
EndTime = 0.02
[Output.FieldData]
)";

std::string Subst(std::string s, const std::string &key, const std::string &val)
{
	for(size_t p = s.find(key); p != std::string::npos; p = s.find(key, p + val.size()))
		s.replace(p, key.size(), val);
	return s;
}

Input *MakeInput(const char *dim, int nx, int ny, int nz, const char *model)
{
	std::string t = kTemplate;
	t = Subst(t, "@DIM@", dim);
	t = Subst(t, "@NX@", std::to_string(nx));
	t = Subst(t, "@NY@", std::to_string(ny));
	t = Subst(t, "@NZ@", std::to_string(nz));
	t = Subst(t, "@MODEL@", model);
	const std::string path = yzfuzz::TempPath("yz_fuzz_restart_tpl", ".toml");
	yzfuzz::WriteFile(path, (const uint8_t*)t.data(), t.size());
	Input *inp = new Input();
	if(!YZFUZZ_TRAP_EXIT(inp->Read(path)))
	{
		std::fprintf(stderr, "fuzz_restart: the template input was refused\n");
		std::abort();
	}
	return inp;
}

struct Fields
{
	Input *in1, *in2, *inz1, *inz2;
	Field_WM1D wm1;
	Field_WM2D wm2;
	Field_CZ   cz1, cz2;
};

Fields &F()
{
	static Fields *f = 0;
	if(!f)
	{
		f = new Fields();
		f->in1  = MakeInput("xh",  8, 1, 1, "WM");
		f->in2  = MakeInput("xyh", 6, 5, 1, "WM");
		f->inz1 = MakeInput("xh",  6, 1, 4, "CZ");
		f->inz2 = MakeInput("xyh", 5, 4, 3, "CZ");
		f->wm1.Build(f->in1->Grid, f->in1->Substrate, 1.0);
		f->wm2.Build(f->in2->Grid, f->in2->Substrate, 1.0);
		f->cz1.Build(f->inz1->Grid, f->inz1->Substrate, 1.0);
		f->cz2.Build(f->inz2->Grid, f->inz2->Substrate, 1.0);
	}
	return *f;
}

//  What the well-mixed model does with a state file it has accepted
//  (WellMixedModelBase.h, BuildRestart and Run), short of the solver objects.
void TakeState(const C_RestartState &st)
{
	C_TimeSteps ts(1e-3, 0, 1, true);
	std::vector<std::string> items;
	items.push_back("Field"); items.push_back("Probe");
	items.push_back("History"); items.push_back("Calc");
	for(size_t k = 0; k < items.size(); k++) ts.Add_OutputItem(items[k], 0.1);
	st.TakeClock(ts, items, true, "Field", 3);
	volatile int a = (int)st.Get("steadyCount");
	volatile int b = (int)st.Get("adjustCount.height");
	volatile double c = st.Get("evaporated") + st.Get("outflowed", 0) + st.Get("outflowed", 1)
	                  + st.Get("volume0", 0) + st.Get("volume0", 1) + st.Get("conc.previousDt");
	(void)a; (void)b; (void)c;
	Vector<double> v(10);
	(void)st.GetVectorInto("height.old", v);
	double diff = 0;
	Field_WM1D &f = F().wm1;
	(void)st.RefineField("field.h", f.HH(), std::vector<char>((size_t)f.GetGrid().imax() + 1, 1), 1e-6, diff);
}

}  // namespace

extern "C" int LLVMFuzzerInitialize(int *, char ***)
{
	(void)F();
	return 0;
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
	if(size < 1) return 0;
	const int which = data[0] % 9;
	static const std::string dat   = yzfuzz::TempPath("yz_fuzz_restart", ".dat");
	static const std::string vts   = yzfuzz::TempPath("yz_fuzz_restart", ".vts");
	static const std::string state = C_RestartState::PathFor(dat);
	const std::string &path = (which == 8) ? state : (which >= 2 && which <= 5) ? vts : dat;
	if(!yzfuzz::WriteFile(path, data + 1, size - 1)) return 0;

	Fields &f = F();
	Input::I_Restart r;
	r.Enable = true; r.ASCII = false; r.FileName = path; r.RestartTime = 0; r.FileIndex = 0;

	const bool returned = YZFUZZ_TRAP_EXIT(
		switch(which)
		{
		case 0: (void)C_Restart::ReadAscii_WM1D(r, f.wm1, 1.0); break;
		case 1: (void)C_Restart::ReadAscii_WM2D(r, f.wm2, 1.0); break;
		case 2: (void)C_Restart::Read(r, f.wm1, 1.0); break;
		case 3: (void)C_Restart::Read(r, f.wm2, 1.0); break;
		case 4: (void)C_Restart::Read(r, f.cz1, 1.0); break;
		case 5: (void)C_Restart::Read(r, f.cz2, 1.0); break;
		case 6: (void)C_Restart::ReadAscii_CZ(r, f.cz1, 1.0); break;
		case 7: (void)C_Restart::ReadAscii_CZ(r, f.cz2, 1.0); break;
		case 8:
		{
			C_RestartState st;
			std::string msg;
			//  The model names are those of the models; the fuzzer finds them
			//  in the dictionary.
			if(st.Load(dat, "Well-Mixed-1D", 0.0, msg)) TakeState(st);
			break;
		}
		});
	if(!returned && yzfuzz::LastExitCode() == 0) std::abort();
	return 0;
}
