/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CZ_HeightEquation_2D.cpp : implementation of the film-height equation of the ConcZ model
*
*  HeightEquation.cpp
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "Opt.h"
#include <cstdio>
#include "Profile.h"
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include "TimeSteps.h"
#include "CZ_HeightEquation_2D.h"


#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

using namespace yuzen::util;
using namespace Global;


//##########################################################
void CZ_HeightEquation_2D::Build(Input &inp, C_Fluid *prop, C_TimeSteps *ts, CZ_ConcEquation_2D *concEq,
	Field_CZ *ff, NonDimensionalNumbers *nonDim, BoundaryCondition2D *bc)
{
	//Set up the pointers.
	_timeStep = ts;
	_field = ff;

	_prop = prop;
	_boundary = bc;
	_NonDim = nonDim;

	_concEq = concEq;
	swConcField = true;


	//Set up everything to do with the grid.
	_G2 = ff->GetGridPtr2();
	_G3 = ff->GetGridPtr3();
	_Nx = _G3->nx();
	_Ny = _G3->ny();
	_Nz = _G3->nz();

	//Switches of the physical models
	_sw_LaplacePressure		= inp.Physics.LaplacePressure.Enable;
	_sw_SolutoCapillarity	= inp.Physics.SolutoCapillarity.Enable;
	_sw_Centrifugal			= inp.Physics.Centrifugal.Enable;
	_sw_VerticalGravity		= inp.Physics.VerticalGravity.Enable;
	_sw_HorizontalGravity	= inp.Physics.HorizontalGravity.Enable;

	//  Center of rotation of the centrifugal force. It is held in
	//  non-dimensional form (treated the same way as in
	//  CentrifugalKernel::Build on the well-mixed side).
	_rotCenter[0] = _rotCenter[1] = 0.0;
	for(size_t a = 0; a < inp.Physics.Centrifugal.RotationAxis.size() && a < 2; ++a)
		_rotCenter[a] = inp.Physics.Centrifugal.RotationAxis[a] / nonDim->ExtentScale();

	_dtdh.resize(_Nx*_Ny);
	{
		const size_t n2 = (size_t)(_Nx*_Ny);
		_qCapX.assign(n2, 0.0); _qExpX.assign(n2, 0.0);
		_qCapY.assign(n2, 0.0); _qExpY.assign(n2, 0.0);
		_qX.assign(n2, 0.0);    _qY.assign(n2, 0.0);
		_colSrc.assign(n2, 0.0); _sideSink.assign(n2, 0.0); _evapSink.assign(n2, 0.0);
	}
	//  The concentration equation takes its volume fluxes from these
	concEq->LinkFilmFluxes(&_qX, &_qY, &_sideSink, &_ppSide, &_colSrc);
	_mms.Build(inp, *nonDim, _Ny != 1, true);

	_fdm.Build(C_Interpolation_base4::Central_2nd);

	_VecB = 0;

	//  The start of the step, which SetMatrixA reads (the evaporation
	//  limit), from the first assembly on
	_field_tmp = *_field;

	//Build the equation solver.
	//This function calls SetMatrix, so it must not run until every model has been built.
	C_ImplicitSolver::Build(ff->HH(), *ts, C_ImplicitSolver::CrankNicholson, inp.Numerics.TimeStepControl.Enable);


	//------------------------
	//  Everything to do with the linear solver
	//------------------------
	//  Set the preconditioner, the scaling and the convergence criterion
	//  together. What decides the shape of the linear system is the in-plane
	//  dimension of the film. ConcZ also has a grid in the film-thickness
	//  direction, but that has no bearing on the "auto" choice.
	this->_bcg.Configure(inp.Numerics.HeightEquation.LinearSolver, (_Ny != 1) ? 2 : 1);

	//  Below this the midpoint of the coupling passes counts as settled
	//  whatever it moved by in the first pass: ten times what the
	//  concentration solve may leave of its residual (Numerics.ConcEquation.
	//  LinearSolver.Epsilon), and no less than 1e-10 (c is a mass fraction).
	_couplingFloor = std::max(1e-10, 10.0*inp.Numerics.ConcEquation.LinearSolver.Epsilon);

	if(inp.Numerics.TimeStepControl.Enable)
	{
		this->SetAdjustTimeChangeFactor(inp.Numerics.TimeStepControl.ChangeFactor);
		this->SetAdjustTimeCriterion(inp.Numerics.TimeStepControl.Criterion);
		this->SetAdjustTimeScaleFloor(inp.Numerics.TimeStepControl.ScaleFloor);
		this->SetAdjustTimeKeepSteps(inp.Numerics.TimeStepControl.KeepSteps);
	}

	// Copy the field into the temporary storage.
	_field_tmp = *_field;

}

bool CZ_HeightEquation_2D::EveryTimeStepTrial()
{
	//Restore the field before the time-step trial.
	_field->operator=(_field_tmp);
	return true;
}

//  After a trial has been accepted: write the film height back into the field,
//  hand its dh/dt to the concentration equation and then advance the
//  concentration by the same step. The concentration is kept at or below 1
//  without taking solute away (RedistributeExcess), and the field is put
//  aside in the temporary storage ready for the next trial.
//
//  *** Two things here keep a concentration stratified through the film
//  from limiting the step as dx^2. ***
//
//  (1) The velocity the solute is carried with is formed here, from the film
//  halfway to the end of the step and the concentration the corrector took,
//  not the one the corrector's SetMatrixA() leaves, formed from the film
//  halfway to the *prediction*. On the short waves the capillary pressure
//  is stiff: within the step the film takes up the shape at which its
//  capillary flux cancels the Marangoni flux. The predicted film has done
//  so for the Marangoni stress of c at the start of the step, the end state
//  for that of the extrapolated c, and a velocity with its capillary part
//  from the first and its Marangoni part from the second carries, besides
//  the flux the film moved with, the Marangoni flux of the difference
//  c_mid - c_n. The zeta velocity at the free surface is then not zero but
//  minus the divergence of that over h, and it moves a concentration
//  stratified through the film across its own gradient:
//      c_t ~ D_M lap(c_mid - c_n),  D_M = h (dsigma/dc) (dc/dzeta) / (2 mu),
//  with c_mid - c_n = (c_n - c_(n-1))/2. A mode of the shortest wave then
//  goes as (g - 1)(g + 2 D_M dt/dx^2) = 0 from step to step, which is
//  stable only for dt < dx^2 / (2 D_M); in cz_strat_1d (tools/mms.py) that
//  is 4e-4 s at 128 cells and 1e-4 s at 256 (D_M = 5e-9 m^2/s).
//
//  (2) With the Marangoni term on, the film and the concentration are then
//  solved again, one stage each, with every coefficient at the midpoint the
//  pass before found, until that midpoint settles: the Marangoni stress, the
//  surface tension and the viscosity are then those of the midpoint of the
//  step itself instead of an extrapolation, and the step is the implicit
//  midpoint rule for the coupled pair. Without it a second limit remains
//  (the Marangoni flow formed from an extrapolated concentration; over
//  0.08 s of cz_strat_1d, 2e-3 s at 64 cells, 1e-3 s at 128 and 5e-4 s at
//  256 go unstable). The passes stop once the midpoint moves by less
//  than a tenth of what the extrapolation missed it by, so what is left over
//  is O(dt^2) and the order is kept.
//
//  The passes are a fixed-point iteration, and the coupling is stiff on the
//  short waves: plain, it converges only while dt is below a limit that goes
//  as dx^2 again (larger than that of (1), but reached in an evaporating
//  film: 480 cells over 1 mm, 5e-4 s). Each pass therefore takes the
//  midpoint Aitken-relaxed, cm + w (c_mid(cm) - cm), with w from the change
//  of the residual between passes (w = 1 on the first), which damps the
//  stiff modes: in that film every step tried converged, up to 5e-3 s (100
//  to 200 times the step at which the scheme without (1) and (2) diverged),
//  in three or four passes.
//
//  Linearized about a flat film with c0(zeta) through it, a pass maps the
//  midpoint c of a wave exp(ikx) through its surface value alone (rank one),
//  with the gain lambda_k of tools/cz_coupling.py, and
//      1 - lambda_k = det(I - dt/2 A_k) / [(1 + a_k) det(B_k)]
//  (A_k the linearized system, 1 + a_k and B_k the Crank-Nicolson factors
//  of the film and of the diffusion through it, both positive): lambda_k >
//  1 only where a wave grows faster than 2/dt, and lambda_k < 0 where
//  sigma_c dc/dzeta < 0. The relaxation tends to 1/(1 - lambda) of the wave
//  that holds the residual, and so removes one wave of any gain, but not
//  waves with gains on both sides of 1 together; a viscosity that depends
//  on c couples every level as well, which is measured, not derived. Passes
//  that do not settle fail the step (cut under step control, the run
//  stopped at a fixed step), as a solve that does not converge does
//  (ch-verification, "The coupling passes of CZ"). The midpoint counts as settled, too, once it
//  moves by less than _couplingFloor, what the concentration solve leaves.
//  YUZEN_DEBUG_COUPLING prints every pass.
namespace {
//  DT().now set to the step's dt for as long as this lives, then restored
//  (ChangeDT keeps old and now, so the pair leaves the two as they were).
struct StepDtScope
{
	C_TimeSteps *ts; double next; bool swapped;
	StepDtScope(C_TimeSteps *t, double stepDt) : ts(t), next(t->DT().now), swapped(false)
	{
		if(stepDt > 0.0 && next != stepDt) swapped = ts->ChangeDT(stepDt);
	}
	~StepDtScope() { if(swapped) ts->ChangeDT(next); }
};
}

bool CZ_HeightEquation_2D::AfterTimeStepTrial()
{
	//  Everything below belongs to the step that was solved, with its dt.
	//  When the time-step control has just lengthened dt, DT().now holds the
	//  next step's (TimeStepAdjustment runs before this); taken by the
	//  conservative update, dh/dt and the concentration, the film would lose
	//  the evaporation of the longer step, more than the solvent the limit of
	//  the solved one allows, and c would pass 1 (tools/analytic.py dried_cz
	//  with the time-step control on: 1.013).
	StepDtScope stepDt(_timeStep, _solvedDt);
	ConservativeUpdate();
	_field->Update_H(_phi_new);
	_couplingPasses = 0;

	if(swConcField)
	{
		int ch = 0, iter = 0; double resid = 0;
		//  c at the start of the step, and at its midpoint: for the first pass
		//  extrapolated from the last two starts, as the corrector took it.
		Vector<double> c0 = _concEq->GetX();
		Vector<double> cm;
		if(!MidpointConc(cm)) cm = c0;
		//  What a failed step has to leave as it was, to be solved again
		//  from it at a shorter step
		_concEq->KeepStepStart();
		const Vector<double> ccPrevKeep = _ccPrev;
		const double dtPrevKeep = _dtPrev;
		const bool ccPrevSetKeep = _ccPrevSet;
		_ccPrev = c0;					// for the next step's extrapolation
		_dtPrev = _timeStep->DT().now;
		_ccPrevSet = true;

		static const bool dbgCoupling = (std::getenv("YUZEN_DEBUG_COUPLING") != 0);
		const int maxPass = (_sw_SolutoCapillarity) ? 1 + kMaxCouplingPasses : 1;
		double e0 = 0, e = 0;
		double omega = 1.0;				// Aitken relaxation of the passes
		Vector<double> r, rPrev;		// residual c_mid(cm) - cm of this pass and the last
		int pass;
		for(pass = 0; pass < maxPass; pass++)
		{
			Vector<double> mid = _phi_old; mid += _phi_new; mid *= 0.5;
			if(pass > 0)
			{
				//  The film again from the start of the step, with every
				//  coefficient at the midpoint the last pass found.
				_field->Update_H(mid);
				_field->Update_C(cm);
				int it2 = 0; double r2 = 0;
				ch = SolveStage(it2, r2);
				_field->Update_C(c0);
				if(ch != 0) break;
				ConservativeUpdate();
				mid = _phi_old; mid += _phi_new; mid *= 0.5;
			}
			CalDtDh();
			_concEq->SetDtDh(GetDtDh());
			//  The velocity from this film: (1) above
			_field->Update_H(mid);
			_field->Update_C(cm);
			CalVelocity();
			_field->Update_C(c0);
			//  With the passes to come, the first one needs no predictor of
			//  its own: it is taken in one stage at the extrapolated midpoint.
			{
				YUZEN_PROF("c.solve");
				if(pass > 0)          ch = _concEq->SolveAtMidpoint(cm, iter, resid);
				else if(maxPass > 1)  ch = _concEq->AdvanceAtMidpoint(cm, iter, resid);
				else                  ch = _concEq->TimeAdvance(iter, resid);
			}
			if(ch != 0) break;
			if(maxPass == 1) break;
			_couplingPasses = pass;

			//  The midpoint of the concentration this pass came to, and how far
			//  it is from the one it was taken with.
			r = c0; r += _concEq->GetX(); r *= 0.5; r -= cm;
			e = r.AbsMax();
			if(pass == 0) e0 = e;
			//  Converged once the midpoint has settled to a tenth of what the
			//  extrapolation missed it by (O(dt^2), so the order is kept
			//  whatever the number of passes).
			const bool settled = !(e > kCouplingTol*e0) || e <= _couplingFloor;
			if(pass > 0 && !settled)
			{
				double num = 0, den = 0;
				for(int i = 0; i < r.size(); i++)
				{
					const double d = r[i] - rPrev[i];
					num += rPrev[i]*d;
					den += d*d;
				}
				if(den > 0) omega = -omega*num/den;
			}
			//  YUZEN_DEBUG_COUPLING: each pass to standard error,
			//    [couple] time, pass, |midpoint moved|, that over the first, relaxation
			if(dbgCoupling)
				std::fprintf(stderr, "[couple] %.6e %d %.6e %.6e %.6e\n", _timeStep->NowTime(),
				             pass, e, (e0 > 0) ? e/e0 : 0.0, settled ? 0.0 : omega);
			if(settled) break;
			for(int i = 0; i < cm.size(); i++) cm[i] += omega*r[i];
			rPrev = r;
		}
		//  Passes that did not settle leave a step that is not the implicit
		//  midpoint step, and where they diverge one far from it (the
		//  coupling's gain is past 1 only where a wave grows faster than 2/dt;
		//  ch-verification, "The coupling passes of CZ"): the step has failed,
		//  as when a solve does not converge.
		if(ch == 0 && maxPass > 1 && pass == maxPass)
		{
			std::ostringstream why;
			why<<"the film and the concentration did not settle in "<<kMaxCouplingPasses
			   <<" coupling passes (the midpoint moved by "<<e<<", against "<<e0<<" in the first)";
			_failCause = why.str();
			_couplingCuts++;
			if(_adjustTimeStep && !_warnedCoupling)
			{
				Global::Warning<<"ConcZ: "<<_failCause<<"; the step is cut and solved again"
				               " (JobInfo.txt counts the steps cut for this; later ones are not reported here)."<<nl;
				_warnedCoupling = true;
			}
			ch = -1;
		}
		_field->Update_H(_phi_new);
		if(ch != 0)
		{
			//  A linear solve that did not converge, of the concentration or
			//  of the film in a coupling pass, or passes that did not settle:
			//  the step has failed. It is cut and solved again under step
			//  control, and stops the run at a fixed step (TimeAdvance).
			//  Everything this step moved on is put back for that.
			if(_failCause.empty())
				_failCause = "the concentration solve did not converge (Numerics.ConcEquation.LinearSolver)";
			_concEq->BackToStepStart();
			_ccPrev = ccPrevKeep;
			_dtPrev = dtPrevKeep;
			_ccPrevSet = ccPrevSetKeep;
			_couplingPasses = 0;
			return false;
		}
		//  c at or below 1, on the equation's own state (the field follows
		//  it), without taking solute away (RedistributeExcess). A cap of the
		//  field alone would lose the solute above it from the field while the
		//  equation went on with it (tools/proptest.py seed 1118).
		{
			double cmaxLeft = 0.0;
			RedistributeExcess(_concEq->GetX(), cmaxLeft);
			if(cmaxLeft > 1.0 + 1e-12 && !_warnedExcess)
			{
				Global::Warning<<"ConcZ: a column holds more solute than its film (largest c "<<cmaxLeft
				               <<"); c is left above 1 there rather than the solute cut."<<nl;
				_warnedExcess = true;
			}
			//  ... and at or above 0, without making solute
			//  (the lowest c the step started from first, for the
			//  divergence test: MinStartConc)
			{
				const bool has_y = (_Ny != 1);
				const int jlo = has_y ? 2 : 0, jhi = has_y ? _Ny-2 : 1;
				_minStart = 0.0;
				for(int ix = 2; ix < _Nx-2; ix++)
				for(int jy = jlo; jy < jhi; jy++)
				for(int kz = 2; kz < _Nz-2; kz++)
				{
					const double c0 = _field_tmp.CC(_G3->At<Cell>(ix, jy, kz));
					if(c0 < _minStart || !(c0 == c0)) _minStart = c0;
				}
			}
			RedistributeDeficit(_concEq->GetX(), _deficitMoved, _deficitAdded);
		}
		_field->Update_C(_concEq->GetX());
		//  The ghosts of the end state from its own conditions
		_concEq->RefreshGhosts();
		_field->Update_C(_concEq->GetX());
	}

	// Copy the updated field into the temporary storage.
	_field_tmp = (*_field);

	return true;
}

double CZ_HeightEquation_2D::ColumnMeanConc(Field_CZ &f, int ix, int jy) const
{
	double s = 0.0;
	for(int kz = 2; kz < _Nz-2; kz++) s += f.CC(_G3->At<Cell>(ix, jy, kz));
	return s / (double)(_Nz-4);
}

//  The cells of a column have the same volume, h dzeta, so moving c from
//  one cell of a column to another moves solute and keeps the column's sum.
//  With the evaporation limited to the solvent the column holds (SetMatrixA)
//  the mean c of a column stays at or below 1; a cell can still pass 1 where
//  the solvent is drawn out of the top cell faster than diffusion brings it
//  up within the step. What it holds above 1 then goes to the cell below:
//  the solvent that evaporated beyond what the top cell had came from the
//  cells beneath it. Where nothing passes 1 nothing is touched.
double CZ_HeightEquation_2D::RedistributeExcess(Vector<double> &c, double &cmaxLeft) const
{
	const bool has_y = (_Ny != 1);
	const int jlo = has_y ? 2 : 0, jhi = has_y ? _Ny-2 : 1;
	const int kb = 2, kt = _Nz-3;
	double moved = 0.0;
	cmaxLeft = 0.0;
	for(int ix = 2; ix < _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		bool any = false;
		for(int kz = kb; kz <= kt; kz++) if(c[IL(ix, jy, kz)] > 1.0) { any = true; break; }
		if(!any) continue;
		//  Down, from the surface to the substrate
		for(int kz = kt; kz > kb; kz--)
		{
			double &ck = c[IL(ix, jy, kz)];
			if(ck > 1.0) { const double e = ck - 1.0; ck = 1.0; c[IL(ix, jy, kz-1)] += e; moved += e; }
		}
		//  and what reached the substrate back up into the cells with room
		double &cb = c[IL(ix, jy, kb)];
		if(cb > 1.0)
		{
			double e = cb - 1.0;
			cb = 1.0;
			for(int kz = kb+1; kz <= kt && e > 0.0; kz++)
			{
				double &ck = c[IL(ix, jy, kz)];
				const double take = std::min(e, 1.0 - ck);
				if(take > 0.0) { ck += take; e -= take; moved += take; }
			}
			//  A column with more solute than film: spread evenly
			if(e > 0.0)
				for(int kz = kb; kz <= kt; kz++) c[IL(ix, jy, kz)] += e / (double)(kt - kb + 1);
		}
		for(int kz = kb; kz <= kt; kz++) cmaxLeft = std::max(cmaxLeft, c[IL(ix, jy, kz)]);
	}
	return moved;
}

//  c at or above 0 without making solute (see the declaration). The
//  cells of a column have the same volume, h dzeta, so within a column c
//  is moved as it stands; between columns the solute of a column is h
//  times the mean c of its cells. Where no real cell is below 0 nothing is
//  touched. The ghosts are formed from the real cells afterwards
//  (RefreshGhosts).
void CZ_HeightEquation_2D::RedistributeDeficit(Vector<double> &c, double &moved, double &added) const
{
	moved = added = 0.0;
	const bool has_y = (_Ny != 1);
	const int jlo = has_y ? 2 : 0, jhi = has_y ? _Ny-2 : 1;
	const int kb = 2, kt = _Nz-3;
	const double nLay = (double)(kt - kb + 1);
	auto H = [this](int ix, int jy) { return _field->HH(_G2->At<Zface>(ix, jy)); };

	bool any = false;
	for(int ix = 2; ix < _Nx-2 && !any; ix++)
	for(int jy = jlo; jy < jhi && !any; jy++)
	for(int kz = kb; kz <= kt; kz++) if(c[IL(ix, jy, kz)] < 0.0) { any = true; break; }
	if(!any) return;

	//  Within each column: its cells below 0 from its cells above, in
	//  proportion to what they hold; a column with less than no solute in
	//  all is left with its mean in every cell
	bool columns = false;
	for(int ix = 2; ix < _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		double def = 0.0, pos = 0.0;
		for(int kz = kb; kz <= kt; kz++)
		{
			const double ck = c[IL(ix, jy, kz)];
			if(ck < 0.0) def -= ck; else pos += ck;
		}
		if(def == 0.0) continue;
		moved += def * H(ix, jy) / nLay;
		if(pos >= def)
		{
			const double f = def / pos;
			for(int kz = kb; kz <= kt; kz++)
			{
				double &ck = c[IL(ix, jy, kz)];
				ck = (ck < 0.0) ? 0.0 : ck - f * ck;
			}
		}
		else
		{
			const double m = (pos - def) / nLay;
			for(int kz = kb; kz <= kt; kz++) c[IL(ix, jy, kz)] = m;
			columns = true;
		}
	}
	if(!columns) return;

	//  Between columns: a column with less than no solute, -h m, takes it
	//  from the columns next to it in proportion to their solute h mean(c)
	//  (each scaled as a whole, so no cell of theirs goes below 0); the same
	//  pass structure as the well-mixed RedistributeExcess
	auto colMean = [&](int ix, int jy) {
		double s = 0.0;
		for(int kz = kb; kz <= kt; kz++) s += c[IL(ix, jy, kz)];
		return s / nLay;
	};
	auto scale = [&](int ix, int jy, double g) {
		for(int kz = kb; kz <= kt; kz++) c[IL(ix, jy, kz)] *= g;
	};
	auto real = [&](int ix, int jy) { return ix >= 2 && ix < _Nx-2 && jy >= jlo && jy < jhi; };
	const int nb = has_y ? 4 : 2;
	const int dxs[4] = {-1, 1, 0, 0}, dys[4] = {0, 0, -1, 1};
	std::vector<double> room((size_t)(_Nx * _Ny), 0.0), take((size_t)(_Nx * _Ny), 0.0);
	bool left = true;
	for(int pass = 0; pass < 4 && left; ++pass)
	{
		for(int ix = 2; ix < _Nx-2; ix++)
		for(int jy = jlo; jy < jhi; jy++)
		{
			const double m = colMean(ix, jy), h = H(ix, jy);
			room[(size_t)IL(ix, jy)] = (m > 0.0 && h > 0.0) ? h * m : 0.0;
			take[(size_t)IL(ix, jy)] = 0.0;
		}
		left = false;
		for(int ix = 2; ix < _Nx-2; ix++)
		for(int jy = jlo; jy < jhi; jy++)
		{
			const double m = colMean(ix, jy), h = H(ix, jy);
			if(!(m < 0.0)) continue;
			if(!(h > 0.0)) { added -= h * m; scale(ix, jy, 0.0); continue; }
			const double e = -h * m;
			double r = 0.0;
			for(int q = 0; q < nb; q++)
				if(real(ix + dxs[q], jy + dys[q])) r += room[(size_t)IL(ix + dxs[q], jy + dys[q])];
			if(!(r > 0.0)) { left = true; continue; }
			const double f = (e < r) ? e / r : 1.0;
			for(int q = 0; q < nb; q++)
				if(real(ix + dxs[q], jy + dys[q])) take[(size_t)IL(ix + dxs[q], jy + dys[q])] += f;
			const double mNew = (f < 1.0) ? 0.0 : -(e - r) / h;
			for(int kz = kb; kz <= kt; kz++) c[IL(ix, jy, kz)] = mNew;
			if(f >= 1.0) left = true;
		}
		for(int ix = 2; ix < _Nx-2; ix++)
		for(int jy = jlo; jy < jhi; jy++)
		{
			const double t = take[(size_t)IL(ix, jy)];
			if(t > 0.0)
			{
				scale(ix, jy, 1.0 - t);
				if(colMean(ix, jy) < 0.0) left = true;
			}
		}
	}
	if(!left) return;

	//  What is left: from every column holding solute, in proportion to it
	double lack = 0.0, total = 0.0;
	for(int ix = 2; ix < _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		const double m = colMean(ix, jy), h = H(ix, jy);
		if(m < 0.0) { lack -= h * m; scale(ix, jy, 0.0); }
		else if(h > 0.0) total += h * m;
	}
	const double f = (total > 0.0) ? std::min(1.0, lack / total) : 0.0;
	if(f > 0.0)
		for(int ix = 2; ix < _Nx-2; ix++)
		for(int jy = jlo; jy < jhi; jy++)
			if(H(ix, jy) > 0.0 && colMean(ix, jy) > 0.0) scale(ix, jy, 1.0 - f);
	if(total < lack) added += lack - total;
}

void CZ_HeightEquation_2D::SetFlagBC()
{
	for(Node nd = _G2->First<Node>(); nd != _G2->Last<Node>(); ++nd)
	{
		if(_G2->CheckGhost(nd) != C_Grid2D::Inside)
			this->SetAsBoundary(_G2->IL(nd));
	}
}
//-----------------------------------------------------------------------
//  Mobility at a face (positivity preserving scheme)
//
//  Eq. (23) of Diez & Kondic, J. Comput. Phys. 183 (2002) 274-306,
//
//      D = 2 h1^2 h2^2 / (h1 + h2)
//
//  times the mean of 1/mu at the face (below). It is the same formula
//  as in model/well_mixed/LaplacePressureKernel.h on the well-mixed side.
//
//  *** Do not simplify this into an arithmetic mean. ***
//  In this form the mobility goes to 0 as O(h^2) when either side thins, so
//  the flux vanishes at a face on the dry side and h cannot go negative. With
//  an arithmetic mean a mobility of h^3/2 would remain even when one side is
//  completely dry, and liquid would be drawn out of the dry region.
//  (The property is pinned down by tests/test_positivity_mobility.cpp.)
//
//  The viscosity at the face. The concentration, and with it the viscosity,
//  varies through the film, and the flux is the depth integral of the
//  velocity profile CalVelocity() builds with that viscosity:
//
//      Q = h^3 b int_0^1 (1 - zeta)^2 / mu dzeta + h^2 tau int_0^1 (1 - zeta) / mu dzeta
//
//  for a driving gradient b and a surface stress tau. So the capillary,
//  gravity and centrifugal mobility takes the mean of 1/mu weighted by
//  3 (1 - zeta)^2, and the Marangoni flux the mean weighted by 2 (1 - zeta)
//  (FaceInvViscosity). With a uniform viscosity both are 1/mu and the flux
//  is h^3 b/(3 mu) + h^2 tau/(2 mu); only the h factors carry the
//  positivity preserving form, so that property is untouched.
//
//  The viscosity of the depth-averaged concentration would not do: it is
//  not the depth integral of the profile when the concentration varies
//  through the film, so the film-height equation would move a different
//  flux from the one the velocity profile -- and through continuity the
//  vertical velocity -- carries, the zeta velocity at the free surface would
//  not vanish and the solute column would not follow its own film.
//-----------------------------------------------------------------------
//  The mean of 1/mu over the film at a face, weighted by (1 - zeta)^power:
//  power 2 for the mobility, 1 for the Marangoni flux. At each level the
//  viscosity is that of the concentration of the two cells either side, as
//  in CalVelocity(); the weights are normalized by their own sum, so a
//  uniform viscosity gives 1/mu exactly. Only the real cells of the column,
//  2 .. _Nz-3 (_Nz counts the two ghost layers at each end). When Ny=1 a
//  face in the Y direction is never asked for.
double CZ_HeightEquation_2D::FaceInvViscosity(int ix, int jy, tFaceDir dir, int power) const
{
	const int ax = (dir == FaceX) ? 1 : 0;
	const int ay = (dir == FaceY) ? 1 : 0;
	const double dzeta = _G3->dz();

	double num = 0, den = 0;
	for(int kz = 2; kz < _Nz-2; kz++)
	{
		const double c = 0.5*(_field->CC(_G3->At<Cell>(ix-ax, jy-ay, kz)) + _field->CC(_G3->At<Cell>(ix, jy, kz)));
		const double r = 1.0 - dzeta*(-1.5 + (double)kz);
		const double w = (power == 2) ? r*r : r;
		num += w / _NonDim->Viscosity(c);
		den += w;
	}
	return num / den;
}

//  The surface concentration at a face, for the surface tension of the
//  capillary pressure: the mean of the two columns' surface values, each the
//  mean of the top real cell (_Nz-3) and the ghost above it (_Nz-2), as in
//  CalVelocity() and MarangoniFlux(), not the depth average: the velocity
//  profile, and the capillary pressure at a free surface, take the surface
//  tension of the surface.
double CZ_HeightEquation_2D::FaceSurfConc(int ix, int jy, tFaceDir dir) const
{
	const int ax = (dir == FaceX) ? 1 : 0;
	const int ay = (dir == FaceY) ? 1 : 0;
	return 0.25*(_field->CC(_G3->At<Cell>(ix-ax, jy-ay, _Nz-2)) + _field->CC(_G3->At<Cell>(ix-ax, jy-ay, _Nz-3))
	           + _field->CC(_G3->At<Cell>(ix   , jy   , _Nz-2)) + _field->CC(_G3->At<Cell>(ix   , jy   , _Nz-3)));
}

//  Whether the face lies on a wall (a Symmetry edge; BoundaryCondition::
//  IsWall): the x face ix = 2 or _Nx-2, the y face jy = 2 or _Ny-2.
//  No liquid crosses it, through any term.
//
//  The ghosts mirror the film and the substrate at a wall, so the capillary,
//  the vertical gravity and the Marangoni flux of its face vanish by
//  themselves (to rounding, and to the residual of the solve that sets the
//  ghost rows); but the in-plane gravity and the centrifugal flux are the
//  face mobility, which the mirrored ghost leaves at the film's own, times a
//  force that does not depend on h, and they would take liquid through the
//  wall (tools/proptest.py seed 1015: VolErr 2e-5). The face of a wall therefore
//  has no mobility and no Marangoni flux (FaceMobility, MarangoniFlux): in
//  the matrix, on the right-hand side and in the face fluxes of the
//  conservative update alike; CalVelocity() gives it no velocity, so the
//  layer fluxes of the concentration equation through it are zero too (the
//  well-mixed counterpart: UpwindKernelBase::WallFace).
bool CZ_HeightEquation_2D::WallFace(int ix, int jy, tFaceDir dir) const
{
	if(dir == FaceX)
		return (ix == 2      && BoundaryCondition::IsWall(_boundary->TypeXmin()))
		    || (ix == _Nx-2  && BoundaryCondition::IsWall(_boundary->TypeXmax()));
	return (jy == 2      && BoundaryCondition::IsWall(_boundary->TypeYmin()))
	    || (jy == _Ny-2  && BoundaryCondition::IsWall(_boundary->TypeYmax()));
}

double CZ_HeightEquation_2D::FaceMobility(int ix, int jy, tFaceDir dir) const
{
	if(WallFace(ix, jy, dir)) return 0.0;
	return OpenFaceMobility(ix, jy, dir);
}

double CZ_HeightEquation_2D::OpenFaceMobility(int ix, int jy, tFaceDir dir) const
{
	const int ax = (dir == FaceX) ? 1 : 0;
	const int ay = (dir == FaceY) ? 1 : 0;

	const double h1 = _field->HH(_G2->At<Zface>(ix-ax, jy-ay));
	const double h2 = _field->HH(_G2->At<Zface>(ix   , jy   ));

	const double denom = h1 + h2;
	if(denom <= 0.0) return 0.0;		//  no liquid on either side

	//  Eq. (23). For uniform h it comes out as h^3.
	const double mob = 2.0*(h1*h1)*(h2*h2) / denom;

	//  The flow rate is h^3/(3 mu), so the 1/3 is needed. On the well-mixed side
	//  it is the _Magnitude = 1/3 of LaplacePressureKernel.
	//  (That integrating the velocity profile in z produces the factor 1/3 has
	//   been checked by integrating u = k[(h+s)S1 - S2].)
	return mob / 3.0 * FaceInvViscosity(ix, jy, dir, 2);
}

//  Marangoni flux at a face.
//
//  It is a shear flow driven by the gradient of the surface tension, so the
//  flux is h^2/(2 mu) dsigma/dx (with 1/mu the mean weighted by 2 (1 - zeta),
//  FaceInvViscosity). Unlike the mobility this goes as h^2, so for
//  positivity it has to be written in a form that vanishes when either side
//  thins. h1 h2 is the positivity preserving counterpart of h^2: it is 0
//  whenever either of the two is.
double CZ_HeightEquation_2D::MarangoniFlux(int ix, int jy, tFaceDir dir) const
{
	if(WallFace(ix, jy, dir)) return 0.0;
	const int ax = (dir == FaceX) ? 1 : 0;
	const int ay = (dir == FaceY) ? 1 : 0;

	const double h1 = _field->HH(_G2->At<Zface>(ix-ax, jy-ay));
	const double h2 = _field->HH(_G2->At<Zface>(ix   , jy   ));
	if(h1 + h2 <= 0.0) return 0.0;

	//  Gradient of the surface concentration. When there is a profile across the
	//  thickness the surface value is milder than the average, and it is that
	//  difference that acts here.
	const double ds = (dir == FaceX) ? _G3->dx() : _G3->dy();
	//  The surface value is the mean of the top real cell (_Nz-3) and the
	//  ghost above it (_Nz-2), as in CalVelocity; _Nz-2 alone is the ghost,
	//  half a cell above the surface.
	const double cs1 = 0.5*(_field->CC(_G3->At<Cell>(ix-ax, jy-ay, _Nz-2)) + _field->CC(_G3->At<Cell>(ix-ax, jy-ay, _Nz-3)));
	const double cs2 = 0.5*(_field->CC(_G3->At<Cell>(ix   , jy   , _Nz-2)) + _field->CC(_G3->At<Cell>(ix   , jy   , _Nz-3)));
	const double dcds = (cs2 - cs1)/ds;

	//  The positivity preserving combination that stands in for h^2
	const double h2eff = h1*h2;

	return 0.5*h2eff * _NonDim->As()*_NonDim->Ma()*dcds
	       * FaceInvViscosity(ix, jy, dir, 1);
}

//  Assemble the coefficient matrix from the conservative form of the film-height
//  equation. The capillary pressure is linear in h, so it goes into the matrix;
//  the substrate relief, gravity, the centrifugal force, Marangoni and
//  evaporation go onto the right-hand side. The velocity field is needed by the
//  convection of the concentration equation; it is formed where that is
//  taken, in AfterTimeStepTrial(), from the film of the step, and not here.
void CZ_HeightEquation_2D::SetMatrixA()
{
	_solvedDt = _timeStep->DT().now;
	_crsA.ClearData();
	_VecB = 0;

	double Evap;
	double sPo,spo,smo,sMo;
	int ix, jy;
	double dx = _G3->dx();
	Zface zf;

	//-------------------------------------------------------------------
	//  Conservative form. h_t = -dQ/dx - E, with the flux at a face
	//
	//      Q = D_face * [ sigma d^3(h+s)/dx^3 - Bo d(h+s)/dx + Pgx + Sp x ]
	//          + Ma_face                                   (Marangoni)
	//
	//  D_face is the positivity preserving mobility returned by FaceMobility().
	//  The capillary-pressure term is linear in h, so it goes into the matrix and
	//  the rest onto the right-hand side.
	//
	//  Writing it in conservative form buys three things: the positivity
	//  preserving scheme can be placed on the face mobility, the volume is
	//  conserved exactly, and the capillary term is solved implicitly, so no
	//  explicit stability condition is imposed on the time step.
	//-------------------------------------------------------------------
	const double rdx  = 1.0/dx;
	const double rdx2 = rdx*rdx;
	const double rdx4 = rdx2*rdx2;

	//  Whether this is one-dimensional is read off ny() == 1. In that case the
	//  grid carries no ghost layer in y (SetSite in Grids.h makes the same test
	//  and writes only 5 points).
	//
	//  In one dimension neither the y faces nor the cross terms ever appear, and
	//  the body of the if(has_y) below is never executed. The matrix holds only
	//  the columns that were assigned to (CRS_Matrix::Insert does a push_back), so
	//  a one-dimensional row stays at 5 points, and the number of non-zeros, the
	//  work arrays of the ILU and the cost of one BiCGSTAB iteration all come out
	//  the same as in an implementation written for one dimension only.
	const bool has_y = (_Ny != 1);
	const double dy   = has_y ? _G3->dy() : 1.0;
	const double rdy  = 1.0/dy;
	const double rdy2 = rdy*rdy;
	const double rdy4 = rdy2*rdy2;
	const double rdx2dy2 = rdx2*rdy2;

	const int jlo = has_y ? 2       : 0;
	const int jhi = has_y ? _Ny-2   : 1;

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(ix,jy,Evap,sPo,spo,smo,sMo,zf)
#endif
	for(ix = 2; ix < _Nx-2; ix++)
	for(jy = jlo; jy < jhi; jy++)
	{
		zf = _G2->At<Zface>(ix,jy);

		//  Mobility of the left and right faces. On the boundary side only one side
		//  exists, so it falls to 0.
		const double DL = FaceMobility(ix  , jy, FaceX);	// face ix-1|ix
		const double DR = FaceMobility(ix+1, jy, FaceX);	// face ix|ix+1

		//  Capillary pressure: -d/dx [ D sigma d^3(h+s)/dx^3 ]
		//  The third difference at a face uses 4 points, which makes a 5-point
		//  stencil per row. The coefficients follow from the weights D_L and D_R.
		//  The surface tension of each face, as in y. (With one value, that
		//  of the left face, for both, the two rows that share a face would
		//  see different fluxes through it, and with a surface tension that
		//  depends on the concentration the film would not converge.)
		const double sigL = _NonDim->SurfaceTension(FaceSurfConc(ix  , jy, FaceX));
		const double sigR = _NonDim->SurfaceTension(FaceSurfConc(ix+1, jy, FaceX));
		const double aL = _sw_LaplacePressure ? DL*sigL*rdx4 : 0.0;
		const double aR = _sw_LaplacePressure ? DR*sigR*rdx4 : 0.0;

		//  *** Sign convention ***
		//  What the solver solves is (I + theta dt A) h^{n+1} = h^n
		//  - (1-theta) dt A h^n - dt B, which is in effect h_t = -A h - B. The matrix
		//  A therefore takes the coefficients of h_t **with their signs reversed**.
		//
		//  Expanding h_t = -d/dx[D d^3h/dx^3] over 5 points gives the coefficients
		//      i-2: -D_L,  i-1: 3D_L+D_R,  i: -3D_L-3D_R,
		//      i+1: D_L+3D_R,  i+2: -D_R       (all of them /dx^4)
		//  so it is those with their signs reversed that go into A.
		_crsA(IL(ix, jy), IL(ix-2, jy)) = +aL;
		_crsA(IL(ix, jy), IL(ix-1, jy)) = -3.0*aL -aR;
		_crsA.Diag(IL(ix, jy))          = +3.0*aL +3.0*aR;
		_crsA(IL(ix, jy), IL(ix+1, jy)) = -aL -3.0*aR;
		_crsA(IL(ix, jy), IL(ix+2, jy)) = +aR;

		//---------------------------------------------------------------
		//  The y direction and the cross terms. They are not executed in one
		//  dimension, so the matrix stays at 5 points there.
		//
		//  In two dimensions the third derivative that enters the flux at a face
		//  differs from the one-dimensional case. An x face carries
		//      d^3h/dx^3 + d^3h/(dx dy^2)
		//  and a y face carries
		//      d^3h/dy^3 + d^3h/(dy dx^2)
		//  The second term in each is the x-y cross term, and dropping it would mean
		//  no longer discretizing grad^2 grad^2. C(0)..C(3) of
		//  LaplacePressure_2D.cpp on the well-mixed 2D side are the same thing.
		//---------------------------------------------------------------
		if(has_y)
		{
			const double DB = FaceMobility(ix, jy  , FaceY);	// face jy-1|jy
			const double DT = FaceMobility(ix, jy+1, FaceY);	// face jy|jy+1
			const double sigB = _NonDim->SurfaceTension(FaceSurfConc(ix, jy  , FaceY));
			const double sigT = _NonDim->SurfaceTension(FaceSurfConc(ix, jy+1, FaceY));

			const double bB = _sw_LaplacePressure ? DB*sigB*rdy4 : 0.0;
			const double bT = _sw_LaplacePressure ? DT*sigT*rdy4 : 0.0;

			//  The 5 points of d^3h/dy^3 (same form as in x)
			_crsA(IL(ix, jy), IL(ix, jy-2)) += +bB;
			_crsA(IL(ix, jy), IL(ix, jy-1)) += -3.0*bB -bT;
			_crsA.Diag(IL(ix, jy))          += +3.0*bB +3.0*bT;
			_crsA(IL(ix, jy), IL(ix, jy+1)) += -bB -3.0*bT;
			_crsA(IL(ix, jy), IL(ix, jy+2)) += +bT;

			//  The cross terms: d^3h/(dx dy^2) at an x face and d^3h/(dy dx^2) at a y
			//  face. Both carry a weight of 1/(dx^2 dy^2) and include the four diagonal
			//  points.
			const double cL = _sw_LaplacePressure ? DL*sigL*rdx2dy2 : 0.0;
			const double cR = _sw_LaplacePressure ? DR*sigR*rdx2dy2 : 0.0;
			const double cB = _sw_LaplacePressure ? DB*sigB*rdx2dy2 : 0.0;
			const double cT = _sw_LaplacePressure ? DT*sigT*rdx2dy2 : 0.0;

			//  From the x faces: -d/dx[ D d^3h/(dx dy^2) ]
			_crsA(IL(ix, jy), IL(ix-1, jy-1)) += +cL;
			_crsA(IL(ix, jy), IL(ix-1, jy  )) += -2.0*cL;
			_crsA(IL(ix, jy), IL(ix-1, jy+1)) += +cL;
			_crsA(IL(ix, jy), IL(ix+1, jy-1)) += +cR;
			_crsA(IL(ix, jy), IL(ix+1, jy  )) += -2.0*cR;
			_crsA(IL(ix, jy), IL(ix+1, jy+1)) += +cR;
			_crsA(IL(ix, jy), IL(ix, jy-1)) += -cL -cR;
			_crsA.Diag(IL(ix, jy))          += +2.0*cL +2.0*cR;
			_crsA(IL(ix, jy), IL(ix, jy+1)) += -cL -cR;

			//  From the y faces: -d/dy[ D d^3h/(dy dx^2) ]
			_crsA(IL(ix, jy), IL(ix-1, jy-1)) += +cB;
			_crsA(IL(ix, jy), IL(ix  , jy-1)) += -2.0*cB;
			_crsA(IL(ix, jy), IL(ix+1, jy-1)) += +cB;
			_crsA(IL(ix, jy), IL(ix-1, jy+1)) += +cT;
			_crsA(IL(ix, jy), IL(ix  , jy+1)) += -2.0*cT;
			_crsA(IL(ix, jy), IL(ix+1, jy+1)) += +cT;
			_crsA(IL(ix, jy), IL(ix-1, jy)) += -cB -cT;
			_crsA.Diag(IL(ix, jy))          += +2.0*cB +2.0*cT;
			_crsA(IL(ix, jy), IL(ix+1, jy)) += -cB -cT;
		}

		//---------------------------------------------------------------
		//  Right-hand side: substrate relief, gravity, centrifugal force, Marangoni
		//  and evaporation. None of them is linear in h (or they do not contain h at
		//  all), so they are placed explicitly.
		//---------------------------------------------------------------
		sPo = _field->SS(_G2->At<Zface>(ix+2,jy));
		spo = _field->SS(_G2->At<Zface>(ix+1,jy));
		smo = _field->SS(_G2->At<Zface>(ix-1,jy));
		sMo = _field->SS(_G2->At<Zface>(ix-2,jy));
		const double soo = _field->SS(zf);

		//  Build up the contribution to h_t first, and reverse its sign at the end
		//  when putting it into B (the same convention as above: h_t = -A h - B).
		double dhdt = 0.0;
		double side = 0.0;		// the sideways outflow of 1-D (centrifugal)

		//  Capillary pressure through the substrate relief: -d/dx[ D sigma d^3 s/dx^3 ]
		if(_sw_LaplacePressure)
		{
			const double d3sL = (        spo -3.0*soo +3.0*smo -sMo)*rdx*rdx2;
			const double d3sR = (sPo -3.0*spo +3.0*soo -smo        )*rdx*rdx2;
			dhdt -= (DR*sigR*d3sR - DL*sigL*d3sL)*rdx;
		}

		//  Vertical gravity: the flux is -D Bo d(h+s)/dx, so
		//  h_t = +d/dx[D Bo d(h+s)/dx]. It contains h but is evaluated explicitly,
		//  at the coefficients of the stage.
		if(_sw_VerticalGravity)
		{
			const double hsP = _field->HH(_G2->At<Zface>(ix+1,jy)) + spo;
			const double hsO = _field->HH(zf) + soo;
			const double hsM = _field->HH(_G2->At<Zface>(ix-1,jy)) + smo;
			const double gR = DR*_NonDim->Bo()*(hsP - hsO)*rdx;
			const double gL = DL*_NonDim->Bo()*(hsO - hsM)*rdx;
			dhdt += (gR - gL)*rdx;
		}

		//  In-plane gravity and the centrifugal force: the driving force does not
		//  depend on h, so the divergence of the flux is taken as it stands.
		if(_sw_HorizontalGravity)
		{
			const double fR = DR*_NonDim->Pgx();
			const double fL = DL*_NonDim->Pgx();
			dhdt -= (fR - fL)*rdx;
		}
		if(_sw_Centrifugal)
		{
			const double xR = _G2->XOf<Yedge>(ix+1, jy) - _rotCenter[0];
			const double xL = _G2->XOf<Yedge>(ix  , jy) - _rotCenter[0];
			dhdt -= (DR*_NonDim->Sp()*xR - DL*_NonDim->Sp()*xL)*rdx;
			//  In 1-D the film is flung sideways as well, into the y that is
			//  not solved for: d/dy of the flux D Sp (y - y0) is D Sp. The
			//  well-mixed model counts this outflow too, so the two models
			//  agree on a spinning strip.
			//  (The mobility of the column itself: that of its two faces
			//  as the film has them, a wall or not.)
			if(_Ny == 1)
			{
				side = 0.5*(OpenFaceMobility(ix, jy, FaceX) + OpenFaceMobility(ix+1, jy, FaceX))*_NonDim->Sp();
				dhdt -= side;
			}
		}

		//  Marangoni: the flux amounts to h^2/(2 mu) * dsigma/dx. It is set by the
		//  gradient of the surface concentration, so the concentration profile across
		//  the thickness acts here.
		if(_sw_SolutoCapillarity)
		{
			dhdt -= (MarangoniFlux(ix+1, jy, FaceX) - MarangoniFlux(ix, jy, FaceX))*rdx;
		}

		//---------------------------------------------------------------
		//  The y direction of the right-hand side. Not executed in one dimension.
		//---------------------------------------------------------------
		if(has_y)
		{
			const double DB = FaceMobility(ix, jy  , FaceY);
			const double DT = FaceMobility(ix, jy+1, FaceY);
			const double sigB = _NonDim->SurfaceTension(FaceSurfConc(ix, jy  , FaceY));
			const double sigT = _NonDim->SurfaceTension(FaceSurfConc(ix, jy+1, FaceY));

			//  Capillary pressure through the substrate relief. Both the third
			//  derivative in y and the cross term d^3 s/(dy dx^2) on the y faces are
			//  needed.
			if(_sw_LaplacePressure)
			{
				const double soP = _field->SS(_G2->At<Zface>(ix,jy+2));
				const double sop = _field->SS(_G2->At<Zface>(ix,jy+1));
				const double som = _field->SS(_G2->At<Zface>(ix,jy-1));
				const double soM = _field->SS(_G2->At<Zface>(ix,jy-2));
				const double d3sB = (        sop -3.0*soo +3.0*som -soM)*rdy*rdy2;
				const double d3sT = (soP -3.0*sop +3.0*soo -som        )*rdy*rdy2;
				dhdt -= (DT*sigT*d3sT - DB*sigB*d3sB)*rdy;

				//  The cross terms: d^3 s/(dx dy^2) at an x face and d^3 s/(dy dx^2) at a
				//  y face.
				const double smm = _field->SS(_G2->At<Zface>(ix-1,jy-1));
				const double smp = _field->SS(_G2->At<Zface>(ix-1,jy+1));
				const double spm = _field->SS(_G2->At<Zface>(ix+1,jy-1));
				const double spp = _field->SS(_G2->At<Zface>(ix+1,jy+1));

				//  x difference of d^2 s/dy^2 on faces L and R
				const double xL_ = ((som -2.0*soo +sop) - (smm -2.0*smo +smp))*rdx*rdy2;
				const double xR_ = ((spm -2.0*spo +spp) - (som -2.0*soo +sop))*rdx*rdy2;
				dhdt -= (DR*sigR*xR_ - DL*sigL*xL_)*rdx;

				//  y difference of d^2 s/dx^2 on faces B and T
				const double yB_ = ((smo -2.0*soo +spo) - (smm -2.0*som +spm))*rdy*rdx2;
				const double yT_ = ((smp -2.0*sop +spp) - (smo -2.0*soo +spo))*rdy*rdx2;
				dhdt -= (DT*sigT*yT_ - DB*sigB*yB_)*rdy;
			}

			if(_sw_VerticalGravity)
			{
				const double sop = _field->SS(_G2->At<Zface>(ix,jy+1));
				const double som = _field->SS(_G2->At<Zface>(ix,jy-1));
				const double hsT = _field->HH(_G2->At<Zface>(ix,jy+1)) + sop;
				const double hsO = _field->HH(zf) + soo;
				const double hsB = _field->HH(_G2->At<Zface>(ix,jy-1)) + som;
				const double gT = DT*_NonDim->Bo()*(hsT - hsO)*rdy;
				const double gB = DB*_NonDim->Bo()*(hsO - hsB)*rdy;
				dhdt += (gT - gB)*rdy;
			}

			if(_sw_HorizontalGravity)
				dhdt -= (DT*_NonDim->Pgy() - DB*_NonDim->Pgy())*rdy;

			if(_sw_Centrifugal)
			{
				const double yT = _G2->YOf<Xedge>(ix, jy+1) - _rotCenter[1];
				const double yB = _G2->YOf<Xedge>(ix, jy  ) - _rotCenter[1];
				dhdt -= (DT*_NonDim->Sp()*yT - DB*_NonDim->Sp()*yB)*rdy;
			}

			if(_sw_SolutoCapillarity)
				dhdt -= (MarangoniFlux(ix, jy+1, FaceY) - MarangoniFlux(ix, jy, FaceY))*rdy;
		}

		//  Evaporation removes film, so it enters h_t with a negative sign.
		//  No more than the solvent the column holds, at the state the rate
		//  is taken at and at the start of the step (the step removes
		//  Evap dt from the column, and with its solute conserved the
		//  column is then left with h (1 - mean c) >= 0: its mean c stays
		//  at or below 1; see C_Evaporation::SolventLimit)
		Evap = _NonDim->Ev() * _NonDim->EvapRate(*_field, zf);
		if(Evap > 0.0)
		{
			const double dtE = _timeStep->DT().now;
			Evap = C_Evaporation::SolventLimit(Evap, _field->HH(zf), ColumnMeanConc(*_field, ix, jy), dtE);
			Evap = C_Evaporation::SolventLimit(Evap, _field_tmp.HH(zf), ColumnMeanConc(_field_tmp, ix, jy), dtE);
		}
		dhdt -= Evap;

		//  The manufactured source, at the middle of the step
		double mmsSrc = 0.0;
		if(_mms.Enabled())
		{
			mmsSrc = _mms.Snd(_G2->X(zf), has_y ? _G2->Y(zf) : 0.0,
			                  _timeStep->NowTime() + 0.5*_timeStep->DT().now);
			dhdt += mmsSrc;
		}

		//  What of h_t is not the divergence of a face flux (ConservativeUpdate)
		_colSrc[IL(ix, jy)]   = mmsSrc - Evap - side;
		_sideSink[IL(ix, jy)] = side;
		_evapSink[IL(ix, jy)] = Evap;

		_VecB(IL(ix, jy)) = -dhdt;
	}
	StoreFaceFluxParts();
	SetBoundaryCondition();
}

//-----------------------
//  Setting of the boundary conditions
//-----------------------
//  Write the ghost-layer rows directly. If periodic, they are tied to the real
//  cell on the opposite side; otherwise the real cell on the inside is
//  mirrored.
//  In one dimension (ny()==1) there are no y ghost layers, so only the two
//  layers in x are written.
void CZ_HeightEquation_2D::SetBoundaryCondition()
{
	const bool has_y = (_Ny != 1);
	const int jlo = has_y ? 0     : 0;
	const int jhi = has_y ? _Ny   : 1;

	const bool perX = (_boundary->TypeXmin()==BoundaryCondition::Periodic
	                && _boundary->TypeXmax()==BoundaryCondition::Periodic);

	//  The x ghost layers. In two dimensions they are placed on every y row.
	for(int jy = jlo; jy < jhi; jy++)
	{
		if(perX)
		{
			_crsA.Diag(IL(0, jy)) = 1.0;	_crsA(IL(0,jy), IL(_Nx-4,jy)) = -1.0;
			_crsA.Diag(IL(1, jy)) = 1.0;	_crsA(IL(1,jy), IL(_Nx-3,jy)) = -1.0;

			_crsA.Diag(IL(_Nx-2, jy)) = 1.0;	_crsA(IL(_Nx-2, jy), IL(2,jy)) = -1.0;
			_crsA.Diag(IL(_Nx-1, jy)) = 1.0;	_crsA(IL(_Nx-1, jy), IL(3,jy)) = -1.0;
		}
		else
		{
			_crsA.Diag(IL(0, jy)) = 1.0;	_crsA(IL(0,jy), IL(3,jy)) = -1.0;
			_crsA.Diag(IL(1, jy)) = 1.0;	_crsA(IL(1,jy), IL(2,jy)) = -1.0;

			_crsA.Diag(IL(_Nx-1, jy)) = 1.0;	_crsA(IL(_Nx-1, jy), IL(_Nx-4,jy)) = -1.0;
			_crsA.Diag(IL(_Nx-2, jy)) = 1.0;	_crsA(IL(_Nx-2, jy), IL(_Nx-3,jy)) = -1.0;
		}
	}

	if(!has_y) return;

	//  The y ghost layers. The corners have already been written on the x side,
	//  so only the interior x values are looped over.
	const bool perY = (_boundary->TypeYmin()==BoundaryCondition::Periodic
	                && _boundary->TypeYmax()==BoundaryCondition::Periodic);

	for(int ix = 2; ix < _Nx-2; ix++)
	{
		if(perY)
		{
			_crsA.Diag(IL(ix, 0)) = 1.0;	_crsA(IL(ix,0), IL(ix,_Ny-4)) = -1.0;
			_crsA.Diag(IL(ix, 1)) = 1.0;	_crsA(IL(ix,1), IL(ix,_Ny-3)) = -1.0;

			_crsA.Diag(IL(ix,_Ny-2)) = 1.0;	_crsA(IL(ix,_Ny-2), IL(ix,2)) = -1.0;
			_crsA.Diag(IL(ix,_Ny-1)) = 1.0;	_crsA(IL(ix,_Ny-1), IL(ix,3)) = -1.0;
		}
		else
		{
			_crsA.Diag(IL(ix, 0)) = 1.0;	_crsA(IL(ix,0), IL(ix,3)) = -1.0;
			_crsA.Diag(IL(ix, 1)) = 1.0;	_crsA(IL(ix,1), IL(ix,2)) = -1.0;

			_crsA.Diag(IL(ix,_Ny-1)) = 1.0;	_crsA(IL(ix,_Ny-1), IL(ix,_Ny-4)) = -1.0;
			_crsA.Diag(IL(ix,_Ny-2)) = 1.0;	_crsA(IL(ix,_Ny-2), IL(ix,_Ny-3)) = -1.0;
		}
	}
}
//-----------------------------------------------------------------------
//  The face fluxes of the film-height equation, and the step written
//  through them.
//
//  The concentration equation is in conservative form on the moving zeta
//  grid (CZ_ConcEquation_2D.cpp): its cells exchange liquid through their
//  faces, and the volume flux through the zeta faces is what makes the
//  volume of every cell change by what the film height says. For that to
//  hold exactly (the geometric conservation law: a uniform concentration
//  stays uniform), the layer fluxes of a column have to add up to the flux
//  this equation moved the film with -- face by face -- and the film has to
//  have moved by the divergence of those fluxes. So:
//
//  StoreFaceFluxParts() keeps, from the assembly of the stage that is
//  solved, the coefficient of the capillary term (linear in h, implicit)
//  and the explicit rest of the flux of every face; ConservativeUpdate()
//  forms the flux with the Crank-Nicolson h of the solve and writes the end
//  of the step as h_old + dt (-div q + the column sources). That is the
//  solution of the linear system up to the solver's residual, which is
//  thereby taken out of the volume balance as well.
//-----------------------------------------------------------------------
void CZ_HeightEquation_2D::StoreFaceFluxParts()
{
	const bool has_y = (_Ny != 1);
	const double dx = _G3->dx();
	const double dy = has_y ? _G3->dy() : 1.0;
	const double rdx = 1.0/dx, rdx2 = rdx*rdx;
	const double rdy = 1.0/dy, rdy2 = rdy*rdy;
	const int jlo = has_y ? 2      : 0;
	const int jhi = has_y ? _Ny-2  : 1;
	auto S = [&](int i, int j) { return _field->SS(_G2->At<Zface>(i, j)); };
	auto H = [&](int i, int j) { return _field->HH(_G2->At<Zface>(i, j)); };

	//  x faces: face ix lies between ix-1 and ix, ix = 2 .. _Nx-2
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for
#endif
	for(int ix = 2; ix <= _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		const double D = FaceMobility(ix, jy, FaceX);
		double cap = 0.0, expl = 0.0;
		if(_sw_LaplacePressure)
		{
			cap = D * _NonDim->SurfaceTension(FaceSurfConc(ix, jy, FaceX));
			double d3s = (S(ix+1,jy) -3.0*S(ix,jy) +3.0*S(ix-1,jy) -S(ix-2,jy))*rdx*rdx2;
			if(has_y)
				d3s += ((S(ix,jy+1) -2.0*S(ix,jy) +S(ix,jy-1)) - (S(ix-1,jy+1) -2.0*S(ix-1,jy) +S(ix-1,jy-1)))*rdx*rdy2;
			expl += cap*d3s;
		}
		if(_sw_VerticalGravity)
			expl -= D*_NonDim->Bo()*((H(ix,jy) + S(ix,jy)) - (H(ix-1,jy) + S(ix-1,jy)))*rdx;
		if(_sw_HorizontalGravity)
			expl += D*_NonDim->Pgx();
		if(_sw_Centrifugal)
			expl += D*_NonDim->Sp()*(_G2->XOf<Yedge>(ix, jy) - _rotCenter[0]);
		if(_sw_SolutoCapillarity)
			expl += MarangoniFlux(ix, jy, FaceX);
		_qCapX[IL(ix, jy)] = cap;
		_qExpX[IL(ix, jy)] = expl;
	}
	if(!has_y) return;

	//  y faces: face jy lies between jy-1 and jy, jy = 2 .. _Ny-2
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for
#endif
	for(int ix = 2; ix < _Nx-2; ix++)
	for(int jy = 2; jy <= _Ny-2; jy++)
	{
		const double D = FaceMobility(ix, jy, FaceY);
		double cap = 0.0, expl = 0.0;
		if(_sw_LaplacePressure)
		{
			cap = D * _NonDim->SurfaceTension(FaceSurfConc(ix, jy, FaceY));
			double d3s = (S(ix,jy+1) -3.0*S(ix,jy) +3.0*S(ix,jy-1) -S(ix,jy-2))*rdy*rdy2;
			d3s += ((S(ix-1,jy) -2.0*S(ix,jy) +S(ix+1,jy)) - (S(ix-1,jy-1) -2.0*S(ix,jy-1) +S(ix+1,jy-1)))*rdy*rdx2;
			expl += cap*d3s;
		}
		if(_sw_VerticalGravity)
			expl -= D*_NonDim->Bo()*((H(ix,jy) + S(ix,jy)) - (H(ix,jy-1) + S(ix,jy-1)))*rdy;
		if(_sw_HorizontalGravity)
			expl += D*_NonDim->Pgy();
		if(_sw_Centrifugal)
			expl += D*_NonDim->Sp()*(_G2->YOf<Xedge>(ix, jy) - _rotCenter[1]);
		if(_sw_SolutoCapillarity)
			expl += MarangoniFlux(ix, jy, FaceY);
		_qCapY[IL(ix, jy)] = cap;
		_qExpY[IL(ix, jy)] = expl;
	}
}

void CZ_HeightEquation_2D::ConservativeUpdate()
{
	const bool has_y = (_Ny != 1);
	const double dx = _G3->dx();
	const double dy = has_y ? _G3->dy() : 1.0;
	const double rdx = 1.0/dx, rdx2 = rdx*rdx;
	const double rdy = 1.0/dy, rdy2 = rdy*rdy;
	const int jlo = has_y ? 2      : 0;
	const int jhi = has_y ? _Ny-2  : 1;
	//  The step being solved: DT().now stays that until the clock advances,
	//  not the lengthened next step after a step on which the time-step
	//  control lengthened dt (the film would be moved on by a step the solve
	//  was not made with).
	const double dt = _timeStep->DT().now;
	_stepDt = dt;

	//  The ghosts of the solve hold the copies only to the solver's residual;
	//  made exact first, so that the two sides of a periodic seam see the
	//  same flux
	CopyGhosts(_phi_new);

	//  The h the capillary term of the solve was applied to
	Vector<double> hc = _phi_new;
	hc *= _theta;
	{
		Vector<double> t = _phi_old;
		t *= (1.0 - _theta);
		hc += t;
	}
	auto Hc = [&](int i, int j) { return hc[IL(i, j)]; };

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for
#endif
	for(int ix = 2; ix <= _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		double d3h = (Hc(ix+1,jy) -3.0*Hc(ix,jy) +3.0*Hc(ix-1,jy) -Hc(ix-2,jy))*rdx*rdx2;
		if(has_y)
			d3h += ((Hc(ix,jy+1) -2.0*Hc(ix,jy) +Hc(ix,jy-1)) - (Hc(ix-1,jy+1) -2.0*Hc(ix-1,jy) +Hc(ix-1,jy-1)))*rdx*rdy2;
		_qX[IL(ix, jy)] = _qCapX[IL(ix, jy)]*d3h + _qExpX[IL(ix, jy)];
	}
	if(has_y)
	{
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for
#endif
		for(int ix = 2; ix < _Nx-2; ix++)
		for(int jy = 2; jy <= _Ny-2; jy++)
		{
			double d3h = (Hc(ix,jy+1) -3.0*Hc(ix,jy) +3.0*Hc(ix,jy-1) -Hc(ix,jy-2))*rdy*rdy2;
			d3h += ((Hc(ix-1,jy) -2.0*Hc(ix,jy) +Hc(ix+1,jy)) - (Hc(ix-1,jy-1) -2.0*Hc(ix,jy-1) +Hc(ix+1,jy-1)))*rdy*rdx2;
			_qY[IL(ix, jy)] = _qCapY[IL(ix, jy)]*d3h + _qExpY[IL(ix, jy)];
		}
	}

	//  The end of the step through the fluxes
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for
#endif
	for(int ix = 2; ix < _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		double div = (_qX[IL(ix+1, jy)] - _qX[IL(ix, jy)])*rdx;
		if(has_y) div += (_qY[IL(ix, jy+1)] - _qY[IL(ix, jy)])*rdy;
		_phi_new[IL(ix, jy)] = _phi_old[IL(ix, jy)] + dt*(_colSrc[IL(ix, jy)] - div);
	}
	//  and its ghosts
	CopyGhosts(_phi_new);
}

double CZ_HeightEquation_2D::ColumnSinks(double cellMeasure, const Field_CZ &f, const ScalarField3D<> &cStart,
                                         double &evaporated, double *sideResin) const
{
	const bool has_y = (_Ny != 1);
	const int jlo = has_y ? 2      : 0;
	const int jhi = has_y ? _Ny-2  : 1;
	const double dzeta = _G3->dz();
	double ev = 0, side = 0, sideR = 0;
	for(int ix = 2; ix < _Nx-2; ix++)
	for(int jy = jlo; jy < jhi; jy++)
	{
		ev += _evapSink[IL(ix, jy)];
		const double sd = _sideSink[IL(ix, jy)];
		if(sd == 0.0) continue;
		side += sd;
		if(!sideResin) continue;
		//  The weights of the layers as CZ_ConcEquation_2D::BuildVolumeFluxes
		//  takes them (the velocity profile of the column; uniform without it)
		double wsum = 0.0;
		const bool pp = !_ppSide.empty();
		if(pp)
			for(int kz = 2; kz < _Nz-2; kz++)
				wsum += 0.5*(_ppSide[_G3->At<Xface>(ix, jy, kz).Index()]
				           + _ppSide[_G3->At<Xface>(ix+1, jy, kz).Index()]) * dzeta;
		double cm = 0.0;
		for(int kz = 2; kz < _Nz-2; kz++)
		{
			double o = 1.0;
			if(pp && wsum > 0.0)
				o = 0.5*(_ppSide[_G3->At<Xface>(ix, jy, kz).Index()]
				       + _ppSide[_G3->At<Xface>(ix+1, jy, kz).Index()]) / wsum;
			const Cell cl = _G3->At<Cell>(ix, jy, kz);
			cm += o * 0.5*(cStart(cl) + f.CC(cl)) * dzeta;
		}
		sideR += sd * cm;
	}
	evaporated = ev * cellMeasure * _stepDt;
	if(sideResin) *sideResin = sideR * cellMeasure * _stepDt;
	return side * cellMeasure * _stepDt;
}

//  The ghosts of h from the boundary rows (copies across a wall or a seam;
//  a corner reads an edge, so the rows are swept twice)
void CZ_HeightEquation_2D::CopyGhosts(Vector<double> &h) const
{
	for(int sweep = 0; sweep < 2; sweep++)
	for(int r = 0; r < _crsA.RowSize(); r++)
	{
		if(!_crsA.IsBC(r)) continue;
		double sum = 0.0;
		for(int j = 1; j < _crsA.NonZero(r); j++)
			sum -= _crsA.Element(r, j) * h[_crsA.Jcol(r, j)];
		h[r] = sum / _crsA.Diag(r);
	}
}

//  Rate of change of the film height over this step. It is used for the
//  advection that accounts for the zeta grid of the concentration equation
//  stretching with the film.
void CZ_HeightEquation_2D::CalDtDh()
{
	_dtdh  = _phi_new;
	_dtdh -= _phi_old;
	_dtdh /= _timeStep->DT().now;
}

//-----------------------
//  Computation of the velocity field
//
//  The viscosity varies across the thickness through the concentration, so the
//  in-plane velocity cannot be written analytically. 1/mu and z/mu are
//  integrated across the thickness by Simpson's rule (sum1, sum2) and
//  multiplied by the driving forces of whichever physical terms are enabled
//  (capillary pressure, Marangoni, vertical gravity, horizontal gravity,
//  centrifugal force) to form u(x,z). The vertical velocity w then follows
//  from integrating the continuity equation across the thickness. Finally the
//  ghost layers are filled according to the kind of boundary.
//-----------------------
void CZ_HeightEquation_2D::CalVelocity()
{
	YUZEN_PROF("velocity");
	double hpo, hoo, hmo, hMo, hL;
	double spo, soo, smo, sMo, sL;
	double sigma, dcdx_suf, dx3hs, dxhs;
	double zT, zM, dx, dz, dz_tmp, uu;
	double muB, muT, muM;
	double fB1, fM1, fT1, fB2, fM2, fT2;
	double Cpoo,Cooo,Cmoo,CMoo,Cpom,Coom,Cmom,CMom;
	double xL; //yL;

	double sum1, sum2;
	Xface xf; Yface yf; Zface zf;
	int ix, jy, kz;
	dx = _G3->dx();
	//  dx is a grid constant, so pow() need only be called once. Placed inside
	//  the ix loop it would be evaluated Nx times per call. The value is the
	//  same, so the result is unchanged bit for bit.
	const double dx3 = pow(dx,3.0);

	//  In 1-D with the centrifugal term the velocity profile PP is kept at
	//  every x face, for the sideways divergence Sp PP of the vertical
	//  velocity below.
	const bool sideOut = (_Ny == 1 && _sw_Centrifugal);
	if(sideOut && _ppSide.size() != (size_t)(_G3->imax() + 1)) _ppSide.assign(_G3->imax() + 1, 0.0);

	//Cal Ux

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(hpo,hoo,hmo,hMo,hL,spo,soo,smo,sMo,sL,sigma,dcdx_suf,dx3hs,dxhs,zT,zM,dz,dz_tmp,uu, \
								muB,muT,muM,fB1,fM1,fT1,fB2,fM2,fT2,Cpoo,Cooo,Cmoo,CMoo,Cpom,Coom,Cmom,CMom, \
								xL,sum1,sum2,xf,yf,zf,ix,jy,kz)
#endif
	//  Every row of the in-plane grid, the ghost rows included, so that the
	//  velocity is there wherever the concentration equation and the output
	//  look for it. (jy = 0 alone is right in 1-D, where that is the only
	//  row, but in 2-D it is a ghost row, and every real row would be left
	//  with Ux = 0 -- no x advection of the solute at all.)
	for(ix = 2; ix < _Nx-1; ix++)
	for(jy = 0; jy < _Ny; jy++)
	{

		hpo  = _field->HH(_G2->At<Zface>(ix+1, jy));	hoo = _field->HH(_G2->At<Zface>(ix  , jy));
		hmo  = _field->HH(_G2->At<Zface>(ix-1, jy));	hMo = _field->HH(_G2->At<Zface>(ix-2, jy));
		spo  = _field->SS(_G2->At<Zface>(ix+1, jy));	soo = _field->SS(_G2->At<Zface>(ix  , jy));
		smo  = _field->SS(_G2->At<Zface>(ix-1, jy));	sMo = _field->SS(_G2->At<Zface>(ix-2, jy));
//		hL = (-1.0*hMo +9.0*hmo +9.0*hoo -1.0*hpo)/16.0;
//		sL = (-1.0*sMo +9.0*smo +9.0*soo -1.0*spo)/16.0;
		hL = 0.5*(hmo+hoo);
		sL = 0.5*(smo+soo);

		//  From the center of rotation, as in the film-height equation and in
		//  y. (Taken from the coordinate itself, an axis off x = 0 would move
		//  the film, but not the velocity the solute is carried with.)
		xL = _G2->XOf<Yedge>(ix, jy) - _rotCenter[0];
		//yL = _G2->YOf<Zface>(ix, jy);

		Cpoo = _field->CC(_G3->At<Cell>(ix+1, jy, _Nz-2));
		Cooo = _field->CC(_G3->At<Cell>(ix  , jy, _Nz-2));
		Cmoo = _field->CC(_G3->At<Cell>(ix-1, jy, _Nz-2));
		CMoo = _field->CC(_G3->At<Cell>(ix-2, jy, _Nz-2));
		Cpom = _field->CC(_G3->At<Cell>(ix+1, jy, _Nz-3));
		Coom = _field->CC(_G3->At<Cell>(ix  , jy, _Nz-3));
		Cmom = _field->CC(_G3->At<Cell>(ix-1, jy, _Nz-3));
		CMom = _field->CC(_G3->At<Cell>(ix-2, jy, _Nz-3));


		dcdx_suf = 0.5*((-Cpoo +27.0*Cooo -27.0*Cmoo +CMoo)+(-Cpom +27.0*Coom -27.0*Cmom +CMom))/(24.0*dx);
		dx3hs = ( (hpo+spo)- 3.0*(hoo+soo)+ 3.0*(hmo+smo)-(hMo+sMo))/dx3;
		dxhs  = (-(hpo+spo)+27.0*(hoo+soo)-27.0*(hmo+smo)+(hMo+sMo))/(24.0*dx);
		sigma = _NonDim->SurfaceTension(0.25*(Cooo+Coom+Cmoo+Cmom));

		// Integrate by Simpson's rule.
		Cooo = _field->CC(_G3->At<Cell>(ix  , jy, 2));
		Cmoo = _field->CC(_G3->At<Cell>(ix-1, jy, 2));
		Coom = _field->CC(_G3->At<Cell>(ix  , jy, 1));
		Cmom = _field->CC(_G3->At<Cell>(ix-1, jy, 1));
		muB = _NonDim->Viscosity(0.25*(Cooo+Coom+Cmoo+Cmom));
		fB1 = 1.0/muB;
		fB2 =  sL/muB;
		dz = _G3->dz() * hL;
		sum1 = sum2 = 0;

		//--------------------------------------------------------------
		//  The coefficients of the driving forces do not depend on kz, so they are
		//  formed once per face, outside the kz loop, and ((hL+sL)*sum1 - sum2),
		//  which four of the terms share, once for all of them. Keep the order in
		//  which the products are associated: the values are then the same bit for
		//  bit as when each is formed per kz and per term.
		//--------------------------------------------------------------
		//  In 2-D the capillary pressure gradient at an x face is
		//  d/dx lap(h+s) = d^3/dx^3 + d^3/(dx dy^2): the cross term, as the
		//  film-height equation has it; without it the solute would be carried
		//  by a capillary flow the film does not have (in 2-D it is found from
		//  the rows either side; the outermost ghost rows, which nothing reads,
		//  go without it).
		double dxyyhs = 0.0;
		if(_Ny != 1 && jy >= 1 && jy <= _Ny-2)
		{
			const double dy2 = _G3->dy()*_G3->dy();
			const double fOp = _field->HH(_G2->At<Zface>(ix  , jy+1)) + _field->SS(_G2->At<Zface>(ix  , jy+1));
			const double fOm = _field->HH(_G2->At<Zface>(ix  , jy-1)) + _field->SS(_G2->At<Zface>(ix  , jy-1));
			const double fMp = _field->HH(_G2->At<Zface>(ix-1, jy+1)) + _field->SS(_G2->At<Zface>(ix-1, jy+1));
			const double fMm = _field->HH(_G2->At<Zface>(ix-1, jy-1)) + _field->SS(_G2->At<Zface>(ix-1, jy-1));
			dxyyhs = ((fOp - 2.0*(hoo+soo) + fOm) - (fMp - 2.0*(hmo+smo) + fMm)) / (dx*dy2);
		}
		const double kLap  = _sw_LaplacePressure   ? sigma*(dx3hs + dxyyhs)             : 0.0;
		const double kMa   = _sw_SolutoCapillarity ? _NonDim->As()*_NonDim->Ma()*dcdx_suf : 0.0;
		const double kBo   = _sw_VerticalGravity   ? -_NonDim->Bo()*dxhs                : 0.0;
		const double kPg   = _sw_HorizontalGravity ? _NonDim->Pgx()                     : 0.0;
		const double kSp   = _sw_Centrifugal       ? _NonDim->Sp() * xL                 : 0.0;
		const double hLsL  = hL + sL;

		//  In the film-thickness direction the linear indices are spaced nx*ny apart.
		//  kz runs from 2 to _Nz-1 and never wraps, so the index need only be formed
		//  once and then advanced by that stride.
		const int strideZ = _G3->nx() * _G3->ny();
		int iCoo = _G3->template At<Cell>(ix  , jy, 2).Index();
		int iCmo = _G3->template At<Cell>(ix-1, jy, 2).Index();
		int iXf  = _G3->template At<Xface>(ix, jy, 2).Index();
		for(kz = 2; kz < _Nz; kz++)
		{
			zT = sL + dz*(-1.5+(double)kz);
			//  The midpoint of the Simpson interval below zT. The first
			//  interval is the half cell from the substrate to the first
			//  cell center, so its midpoint is a quarter cell up, not the
			//  substrate itself: taken at the substrate, the profile would be
			//  off by h^2 dzeta^2/(12 mu) at every height even with a
			//  uniform viscosity, for which Simpson's rule is exact.
			zM = (kz == 2) ? sL + 0.25*dz : sL + dz*(-2.0+(double)kz);
			//  Carry the values of the level below down as they are: the lower two of
			//  the four points hold the same values as the upper two of the previous
			//  level, so only the upper two are read. kz==2 is skipped because the
			//  preparatory block has already put the correct four points in place.
			if(kz > 2)
			{
				Coom = Cooo;
				Cmom = Cmoo;
				iCoo += strideZ;	iCmo += strideZ;
				Cooo = _field->CC(Cell(iCoo));
				Cmoo = _field->CC(Cell(iCmo));
			}
			//  (At kz == 2 the quarter-cell point: halfway between the
			//  substrate value and the first cell center.)
			muM = (kz == 2) ? _NonDim->Viscosity(0.375*(Cooo+Cmoo) + 0.125*(Coom+Cmom))
			                : _NonDim->Viscosity(0.25*(Cooo+Coom+Cmoo+Cmom));
			muT = _NonDim->Viscosity(0.50*(Cooo+Cmoo));

			fT1 = 1.0/muT;
			fT2 = zT /muT;
			fM1 = 1.0/muM;
			fM2 = zM /muM;
			dz_tmp = dz;
			if(kz == 2) dz_tmp = 0.5*dz;
			sum1 += (fB1 +fT1 +4.0*fM1)*dz_tmp/6.0;
			sum2 += (fB2 +fT2 +4.0*fM2)*dz_tmp/6.0;

			fB1 = fT1;
			fB2 = fT2;

			//  The capillary pressure, vertical gravity, in-plane gravity and the
			//  centrifugal force all use this same combination.
			const double PP = hLsL*sum1 - sum2;
			uu = 0;
			if( _sw_LaplacePressure)   uu += kLap*PP;
			if( _sw_SolutoCapillarity) uu += kMa*sum1;
			if( _sw_VerticalGravity)   uu += kBo*PP;
			if( _sw_HorizontalGravity) uu += kPg*PP;
			if( _sw_Centrifugal)       uu += kSp*PP;
			if(sideOut) _ppSide[iXf] = PP;

			xf = Xface(iXf);
			_field->Ux(xf) = uu;
			iXf += strideZ;
		}
	}

	//-------------------------------------------------------------------
	//  Velocity in the y direction. Not executed in one dimension (ny()==1).
	//
	//  It is structured exactly like the x loop above; only the orientation of
	//  the faces differs. Nothing in the order of operations on the x side has
	//  been touched, so the one-dimensional results are unchanged bit for bit.
	//-------------------------------------------------------------------
	if(_Ny != 1)
	{
		const double dy  = _G3->dy();
		const double dy3 = pow(dy,3.0);
		const int strideZ = _G3->nx() * _G3->ny();

		double hop, hoo_, hom, hoM, hB;
		double sop, soo_, som, soM, sB;
		double sigmaY, dcdy_suf, dy3hs, dyhs;
		double zT_, zM_, dz_, dzt, vv;
		double mB, mT, mM;
		double gB1, gM1, gT1, gB2, gM2, gT2;
		double s1, s2, yB;
		double Copo,Cooo_,Como,CoMo,Copm,Coom_,Comm,CoMm;
		int jyy, kzz;

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(hop,hoo_,hom,hoM,hB,sop,soo_,som,soM,sB,sigmaY,dcdy_suf,dy3hs,dyhs, \
								zT_,zM_,dz_,dzt,vv,mB,mT,mM,gB1,gM1,gT1,gB2,gM2,gT2,s1,s2,yB, \
								Copo,Cooo_,Como,CoMo,Copm,Coom_,Comm,CoMm,jyy,kzz)
#endif
		//  Every column, the ghost columns included, as for Ux.
		for(int ixx = 0; ixx < _Nx; ixx++)
		for(jyy = 2; jyy < _Ny-1; jyy++)
		{
			hop  = _field->HH(_G2->At<Zface>(ixx, jyy+1));	hoo_ = _field->HH(_G2->At<Zface>(ixx, jyy  ));
			hom  = _field->HH(_G2->At<Zface>(ixx, jyy-1));	hoM  = _field->HH(_G2->At<Zface>(ixx, jyy-2));
			sop  = _field->SS(_G2->At<Zface>(ixx, jyy+1));	soo_ = _field->SS(_G2->At<Zface>(ixx, jyy  ));
			som  = _field->SS(_G2->At<Zface>(ixx, jyy-1));	soM  = _field->SS(_G2->At<Zface>(ixx, jyy-2));

			hB = 0.5*(hom+hoo_);
			sB = 0.5*(som+soo_);

			yB = _G2->YOf<Xedge>(ixx, jyy) - _rotCenter[1];

			Copo = _field->CC(_G3->At<Cell>(ixx, jyy+1, _Nz-2));
			Cooo_= _field->CC(_G3->At<Cell>(ixx, jyy  , _Nz-2));
			Como = _field->CC(_G3->At<Cell>(ixx, jyy-1, _Nz-2));
			CoMo = _field->CC(_G3->At<Cell>(ixx, jyy-2, _Nz-2));
			Copm = _field->CC(_G3->At<Cell>(ixx, jyy+1, _Nz-3));
			Coom_= _field->CC(_G3->At<Cell>(ixx, jyy  , _Nz-3));
			Comm = _field->CC(_G3->At<Cell>(ixx, jyy-1, _Nz-3));
			CoMm = _field->CC(_G3->At<Cell>(ixx, jyy-2, _Nz-3));

			dcdy_suf = 0.5*((-Copo +27.0*Cooo_ -27.0*Como +CoMo)+(-Copm +27.0*Coom_ -27.0*Comm +CoMm))/(24.0*dy);
			dy3hs = ( (hop+sop)- 3.0*(hoo_+soo_)+ 3.0*(hom+som)-(hoM+soM))/dy3;
			dyhs  = (-(hop+sop)+27.0*(hoo_+soo_)-27.0*(hom+som)+(hoM+soM))/(24.0*dy);
			sigmaY = _NonDim->SurfaceTension(0.25*(Cooo_+Coom_+Como+Comm));

			Cooo_= _field->CC(_G3->At<Cell>(ixx, jyy  , 2));
			Como = _field->CC(_G3->At<Cell>(ixx, jyy-1, 2));
			Coom_= _field->CC(_G3->At<Cell>(ixx, jyy  , 1));
			Comm = _field->CC(_G3->At<Cell>(ixx, jyy-1, 1));
			mB = _NonDim->Viscosity(0.25*(Cooo_+Coom_+Como+Comm));
			gB1 = 1.0/mB;
			gB2 =  sB/mB;
			dz_ = _G3->dz() * hB;
			s1 = s2 = 0;

			//  The cross term d^3/(dy dx^2), as in x
			double dyxxhs = 0.0;
			if(ixx >= 1 && ixx <= _Nx-2)
			{
				const double dx2 = _G3->dx()*_G3->dx();
				const double fpO = _field->HH(_G2->At<Zface>(ixx+1, jyy  )) + _field->SS(_G2->At<Zface>(ixx+1, jyy  ));
				const double fmO = _field->HH(_G2->At<Zface>(ixx-1, jyy  )) + _field->SS(_G2->At<Zface>(ixx-1, jyy  ));
				const double fpM = _field->HH(_G2->At<Zface>(ixx+1, jyy-1)) + _field->SS(_G2->At<Zface>(ixx+1, jyy-1));
				const double fmM = _field->HH(_G2->At<Zface>(ixx-1, jyy-1)) + _field->SS(_G2->At<Zface>(ixx-1, jyy-1));
				dyxxhs = ((fpO - 2.0*(hoo_+soo_) + fmO) - (fpM - 2.0*(hom+som) + fmM)) / (dy*dx2);
			}
			const double kLapY = _sw_LaplacePressure   ? sigmaY*(dy3hs + dyxxhs)              : 0.0;
			const double kMaY  = _sw_SolutoCapillarity ? _NonDim->As()*_NonDim->Ma()*dcdy_suf : 0.0;
			const double kBoY  = _sw_VerticalGravity   ? -_NonDim->Bo()*dyhs                  : 0.0;
			const double kPgY  = _sw_HorizontalGravity ? _NonDim->Pgy()                       : 0.0;
			const double kSpY  = _sw_Centrifugal       ? _NonDim->Sp() * yB                   : 0.0;
			const double hBsB  = hB + sB;

			int iCoo2 = _G3->template At<Cell>(ixx, jyy  , 2).Index();
			int iCom2 = _G3->template At<Cell>(ixx, jyy-1, 2).Index();
			int iYf2  = _G3->template At<Yface>(ixx, jyy, 2).Index();

			for(kzz = 2; kzz < _Nz; kzz++)
			{
				zT_ = sB + dz_*(-1.5+(double)kzz);
				zM_ = (kzz == 2) ? sB + 0.25*dz_ : sB + dz_*(-2.0+(double)kzz);	// as in x
				if(kzz > 2)
				{
					Coom_ = Cooo_;
					Comm  = Como;
					iCoo2 += strideZ;	iCom2 += strideZ;
					Cooo_ = _field->CC(Cell(iCoo2));
					Como  = _field->CC(Cell(iCom2));
				}
				mM = (kzz == 2) ? _NonDim->Viscosity(0.375*(Cooo_+Como) + 0.125*(Coom_+Comm))
				                : _NonDim->Viscosity(0.25*(Cooo_+Coom_+Como+Comm));
				mT = _NonDim->Viscosity(0.50*(Cooo_+Como));

				gT1 = 1.0/mT;
				gT2 = zT_/mT;
				gM1 = 1.0/mM;
				gM2 = zM_/mM;
				dzt = dz_;
				if(kzz == 2) dzt = 0.5*dz_;
				s1 += (gB1 +gT1 +4.0*gM1)*dzt/6.0;
				s2 += (gB2 +gT2 +4.0*gM2)*dzt/6.0;

				gB1 = gT1;
				gB2 = gT2;

				const double QQ = hBsB*s1 - s2;
				vv = 0;
				if( _sw_LaplacePressure)   vv += kLapY*QQ;
				if( _sw_SolutoCapillarity) vv += kMaY*s1;
				if( _sw_VerticalGravity)   vv += kBoY*QQ;
				if( _sw_HorizontalGravity) vv += kPgY*QQ;
				if( _sw_Centrifugal)       vv += kSpY*QQ;

				_field->Uy(Yface(iYf2)) = vv;
				iYf2 += strideZ;
			}
		}
	}

	//------------------------------------------------------------------
	//  Ghost values of the velocity. The in-plane components have been
	//  computed on every face from the wall face to the wall face, in every
	//  row and column; what is left is the two layers below the substrate,
	//  the faces outside the walls, and Uz.
	//
	//  Below the substrate (kz = 0, 1): no slip, so the value mirrors the
	//  first level with its sign reversed and the mean across the substrate
	//  is zero. This is taken on every face, the wall and seam faces
	//  included, so that the mean at the bottom of the first cell next to a
	//  periodic seam is not taken against a ghost that nobody has set.
	//------------------------------------------------------------------
	const bool perX = (_boundary->TypeXmin()==BoundaryCondition::Periodic);
	const bool perY = (_Ny != 1) && (_boundary->TypeYmin()==BoundaryCondition::Periodic);

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(ix,jy) if(yuzen::opt::Par(_Nx*_Ny))
#endif
	for(ix = 2; ix < _Nx-1; ix++)
	for(jy = 0; jy < _Ny; jy++)
	{
		_field->Ux(_G3->At<Xface>(ix, jy, 0)) =  0;
		_field->Ux(_G3->At<Xface>(ix, jy, 1)) =  - _field->Ux(_G3->At<Xface>(ix, jy, 2));
	}

	//  Outside the x walls. Periodic: the ghost faces from the other end.
	//  The two seam faces each keep the velocity computed on them, which is
	//  the flux the film-height equation takes through them from that side:
	//  for every periodic force the two are the same, and for one that is
	//  not periodic -- the centrifugal term, whose (x - x0) jumps by the
	//  length of the domain across the seam -- copying the far one made the
	//  continuity of the first column disagree with its film height, a
	//  first-order error in the solute. Otherwise no flux through the wall
	//  face and the ghost faces mirrored with the sign reversed.
	for(kz = 0; kz < _Nz; kz++)
	for(jy = 0; jy < _Ny; jy++)
	{
		if(perX)
		{
			_field->Ux(_G3->At<Xface>(1, jy, kz)) = _field->Ux(_G3->At<Xface>(_Nx-3, jy, kz));
			_field->Ux(_G3->At<Xface>(0, jy, kz)) = _field->Ux(_G3->At<Xface>(_Nx-4, jy, kz));
			_field->Ux(_G3->At<Xface>(_Nx-1, jy, kz)) = _field->Ux(_G3->At<Xface>(3, jy, kz));
		}
		else
		{
			//  No velocity through the wall face, as no flux of the
			//  film-height equation goes through it (WallFace). With the film
			//  and the substrate mirrored at the wall, the capillary and
			//  Marangoni parts of it vanish there by themselves; the in-plane
			//  gravity and the centrifugal part do not, and are set to zero
			//  here as the film-height equation sets their flux to zero there.
			//  The two have to agree: a wall face with a velocity and
			//  no film flux, or the other way round, makes the continuity
			//  equation of the wall cell disagree with its film height.
			if(WallFace(2, jy, FaceX))     _field->Ux(_G3->At<Xface>(2, jy, kz)) = 0.0;
			if(WallFace(_Nx-2, jy, FaceX)) _field->Ux(_G3->At<Xface>(_Nx-2, jy, kz)) = 0.0;
			_field->Ux(_G3->At<Xface>(1, jy, kz)) = - _field->Ux(_G3->At<Xface>(3, jy, kz));
			_field->Ux(_G3->At<Xface>(0, jy, kz)) = - _field->Ux(_G3->At<Xface>(4, jy, kz));
			_field->Ux(_G3->At<Xface>(_Nx-1, jy, kz)) = - _field->Ux(_G3->At<Xface>(_Nx-3, jy, kz));
		}
	}

	// The same for Uy (two dimensions only)
	if(_Ny != 1)
	{
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(ix,jy) if(yuzen::opt::Par(_Nx*_Ny))
#endif
		for(ix = 0; ix < _Nx; ix++)
		for(jy = 2; jy < _Ny-1; jy++)
		{
			_field->Uy(_G3->At<Yface>(ix, jy, 0)) =  0;
			_field->Uy(_G3->At<Yface>(ix, jy, 1)) =  - _field->Uy(_G3->At<Yface>(ix, jy, 2));
		}

		for(kz = 0; kz < _Nz; kz++)
		for(ix = 0; ix < _Nx; ix++)
		{
			if(perY)
			{
				//  (The seam faces keep their own velocity, as in x.)
				_field->Uy(_G3->At<Yface>(ix, 1, kz)) = _field->Uy(_G3->At<Yface>(ix, _Ny-3, kz));
				_field->Uy(_G3->At<Yface>(ix, 0, kz)) = _field->Uy(_G3->At<Yface>(ix, _Ny-4, kz));
				_field->Uy(_G3->At<Yface>(ix, _Ny-1, kz)) = _field->Uy(_G3->At<Yface>(ix, 3, kz));
			}
			else
			{
				//  No velocity through the wall face, as in x.
				if(WallFace(ix, 2, FaceY))     _field->Uy(_G3->At<Yface>(ix, 2, kz)) = 0.0;
				if(WallFace(ix, _Ny-2, FaceY)) _field->Uy(_G3->At<Yface>(ix, _Ny-2, kz)) = 0.0;
				_field->Uy(_G3->At<Yface>(ix, 1, kz)) = - _field->Uy(_G3->At<Yface>(ix, 3, kz));
				_field->Uy(_G3->At<Yface>(ix, 0, kz)) = - _field->Uy(_G3->At<Yface>(ix, 4, kz));
				_field->Uy(_G3->At<Yface>(ix, _Ny-1, kz)) = - _field->Uy(_G3->At<Yface>(ix, _Ny-3, kz));
			}
		}
	}

	// Cal Uz
	//  The vertical velocity w from the continuity equation du/dx + dv/dy +
	//  dw/dz = 0, the derivatives at fixed z, integrated up from w = 0 on the
	//  substrate, in every real column. The velocities are held at fixed
	//  zeta, and d/dx at fixed z is d/dx at fixed zeta less A d/dzeta, A =
	//  (zeta h_x + s_x)/h, so
	//
	//      w(zeta) = -h int_0^zeta div u dzeta' + (zeta h_x + s_x) u(zeta)
	//                - h_x int_0^zeta u dzeta'   (+ the same in y),
	//
	//  div u taken at fixed zeta. It is written for the output; the
	//  concentration equation forms the flux through its zeta faces itself,
	//  from the layer fluxes and the h_t of the step (the discrete geometric
	//  conservation law, CZ_ConcEquation_2D.cpp), which is the discrete
	//  counterpart of -[div(h int_0^zeta u) + zeta h_t]. Two things matter:
	//  - The integral through cell kz is stored on the face at its top
	//    (kz+1), not at its bottom (face kz), which would make w at every
	//    height the one a cell higher up (first order in dzeta).
	//  - Without the last two terms what is stored is not w but
	//    -h int div u at fixed zeta; with A u subtracted from it the zeta
	//    velocity would not vanish at the free surface of a sloping film,
	//    but be ((zeta h_x + s_x) u - h_x int u)/h, and the solute would be
	//    advected through the surface of its own film.
	//  In two dimensions the y divergence is needed too (dropping it makes
	//  the mass fail to balance).
	double dudx, uR, uL, wB, wT;
	dx = _G3->dx();
	const double dyU  = (_Ny != 1) ? _G3->dy() : 1.0;
	const double dzeta = _G3->dz();
	const int juz_lo = (_Ny != 1) ? 2     : 0;
	const int juz_hi = (_Ny != 1) ? _Ny-2 : 1;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(jy,kz,hoo,dudx,dz,uR,uL,wB,wT)
#endif
	for(ix = 2; ix < _Nx-2; ix++)
	for(jy = juz_lo; jy < juz_hi; jy++)
	{
		hoo = _field->HH(_G2->At<Zface>(ix, jy));
		dz = dzeta * hoo;
		//  The slopes at the column center, as the concentration equation
		//  takes them for its metric
		const double hx0 = (_field->HH(_G2->At<Zface>(ix-2, jy)) - 8.0*_field->HH(_G2->At<Zface>(ix-1, jy))
		                  + 8.0*_field->HH(_G2->At<Zface>(ix+1, jy)) - _field->HH(_G2->At<Zface>(ix+2, jy)))/(12.0*dx);
		const double sx0 = (_field->SS(_G2->At<Zface>(ix-2, jy)) - 8.0*_field->SS(_G2->At<Zface>(ix-1, jy))
		                  + 8.0*_field->SS(_G2->At<Zface>(ix+1, jy)) - _field->SS(_G2->At<Zface>(ix+2, jy)))/(12.0*dx);
		double hy0 = 0, sy0 = 0;
		if(_Ny != 1)
		{
			hy0 = (_field->HH(_G2->At<Zface>(ix, jy-2)) - 8.0*_field->HH(_G2->At<Zface>(ix, jy-1))
			     + 8.0*_field->HH(_G2->At<Zface>(ix, jy+1)) - _field->HH(_G2->At<Zface>(ix, jy+2)))/(12.0*dyU);
			sy0 = (_field->SS(_G2->At<Zface>(ix, jy-2)) - 8.0*_field->SS(_G2->At<Zface>(ix, jy-1))
			     + 8.0*_field->SS(_G2->At<Zface>(ix, jy+1)) - _field->SS(_G2->At<Zface>(ix, jy+2)))/(12.0*dyU);
		}
		wB = 0;
		double ui = 0, vi = 0;		// int_0^zeta u dzeta, int_0^zeta v dzeta
		//  Cell kz, between the faces kz (below) and kz+1 (above); the
		//  top one is the ghost above the free surface.
		for(kz = 2; kz <= _Nz-2; kz++)
		{
			uL = _field->Ux(_G3->At<Xface>(ix  , jy, kz));
			uR = _field->Ux(_G3->At<Xface>(ix+1, jy, kz));
			dudx = (uR-uL)/dx;
			ui += dzeta * 0.5*(uL + uR);
			//  u and v on the face above, as the concentration equation takes them
			const double uf = 0.25*(uL + uR + _field->Ux(_G3->At<Xface>(ix, jy, kz+1)) + _field->Ux(_G3->At<Xface>(ix+1, jy, kz+1)));
			double vf = 0;
			if(_Ny != 1)
			{
				const double vS = _field->Uy(_G3->At<Yface>(ix, jy  , kz));
				const double vN = _field->Uy(_G3->At<Yface>(ix, jy+1, kz));
				dudx += (vN-vS)/dyU;
				vi += dzeta * 0.5*(vS + vN);
				vf = 0.25*(vS + vN + _field->Uy(_G3->At<Yface>(ix, jy, kz+1)) + _field->Uy(_G3->At<Yface>(ix, jy+1, kz+1)));
			}
			else if(_sw_Centrifugal)
			{
				//  The sideways outflow of 1-D, dv/dy = Sp PP (see the height
				//  equation), at the cell from its two faces
				dudx += _NonDim->Sp() * 0.5*(_ppSide[_G3->At<Xface>(ix, jy, kz).Index()]
				                            + _ppSide[_G3->At<Xface>(ix+1, jy, kz).Index()]);
			}
			wT = wB - dz*dudx;
			wB = wT;
			const double zf = dzeta * (double)(kz - 1);
			_field->Uz(_G3->At<Zface>(ix, jy, kz+1)) = wT + (zf*hx0 + sx0)*uf - hx0*ui
			                                              + (zf*hy0 + sy0)*vf - hy0*vi;
		}
	}

	// B.C. for Uz: zero on and below the substrate, in every real column.
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(ix,jy) if(yuzen::opt::Par(_Nx*_Ny))
#endif
	for(ix = 2; ix < _Nx-2; ix++)
	for(jy = juz_lo; jy < juz_hi; jy++)
	{
		_field->Uz(_G3->At<Zface>(ix, jy, 0)) = 0;
		_field->Uz(_G3->At<Zface>(ix, jy, 1)) = 0;
		_field->Uz(_G3->At<Zface>(ix, jy, 2)) = 0;
	}
	//  The ghost columns: Uz is even about a wall and periodic across a seam.
	//  In x on every row of real y, then in y on every column (corners
	//  included, from the x ghosts just set).
	for(kz = 0; kz < _Nz; kz++)
	for(jy = juz_lo; jy < juz_hi; jy++)
	{
		if(perX)
		{
			_field->Uz(_G3->At<Zface>(0, jy, kz))		= _field->Uz(_G3->At<Zface>(_Nx-4, jy, kz));
			_field->Uz(_G3->At<Zface>(1, jy, kz))		= _field->Uz(_G3->At<Zface>(_Nx-3, jy, kz));
			_field->Uz(_G3->At<Zface>(_Nx-1, jy, kz))	= _field->Uz(_G3->At<Zface>(3, jy, kz));
			_field->Uz(_G3->At<Zface>(_Nx-2, jy, kz))	= _field->Uz(_G3->At<Zface>(2, jy, kz));
		}
		else
		{
			_field->Uz(_G3->At<Zface>(0, jy, kz))		= _field->Uz(_G3->At<Zface>(3, jy, kz));
			_field->Uz(_G3->At<Zface>(1, jy, kz))		= _field->Uz(_G3->At<Zface>(2, jy, kz));
			_field->Uz(_G3->At<Zface>(_Nx-1, jy, kz))	= _field->Uz(_G3->At<Zface>(_Nx-4, jy, kz));
			_field->Uz(_G3->At<Zface>(_Nx-2, jy, kz))	= _field->Uz(_G3->At<Zface>(_Nx-3, jy, kz));
		}
	}
	if(_Ny != 1)
	for(kz = 0; kz < _Nz; kz++)
	for(ix = 0; ix < _Nx; ix++)
	{
		if(perY)
		{
			_field->Uz(_G3->At<Zface>(ix, 0, kz))		= _field->Uz(_G3->At<Zface>(ix, _Ny-4, kz));
			_field->Uz(_G3->At<Zface>(ix, 1, kz))		= _field->Uz(_G3->At<Zface>(ix, _Ny-3, kz));
			_field->Uz(_G3->At<Zface>(ix, _Ny-1, kz))	= _field->Uz(_G3->At<Zface>(ix, 3, kz));
			_field->Uz(_G3->At<Zface>(ix, _Ny-2, kz))	= _field->Uz(_G3->At<Zface>(ix, 2, kz));
		}
		else
		{
			_field->Uz(_G3->At<Zface>(ix, 0, kz))		= _field->Uz(_G3->At<Zface>(ix, 3, kz));
			_field->Uz(_G3->At<Zface>(ix, 1, kz))		= _field->Uz(_G3->At<Zface>(ix, 2, kz));
			_field->Uz(_G3->At<Zface>(ix, _Ny-1, kz))	= _field->Uz(_G3->At<Zface>(ix, _Ny-4, kz));
			_field->Uz(_G3->At<Zface>(ix, _Ny-2, kz))	= _field->Uz(_G3->At<Zface>(ix, _Ny-3, kz));
		}
	}
}

//-----------------------------------------------------------------------
//  One time step of the film height, second order in dt.
//
//  The film-height equation is linearized: its matrix and its explicit
//  terms (gravity, the centrifugal and Marangoni fluxes, the capillary
//  pressure of the substrate relief, evaporation) are formed from the state
//  the step is taken from. Solved once like that, the step is first order in
//  dt, whatever the Crank-Nicolson weighting. So it is solved twice: once
//  with the coefficients of the start of the step, which gives a predicted
//  end state, and again from the same start with every coefficient and every
//  explicit term formed at the midpoint between the two. That is the
//  implicit midpoint rule for the whole right-hand side, and it is second
//  order. The velocity the concentration equation is carried with is formed
//  afterwards, in AfterTimeStepTrial, from the film halfway to the end
//  state (not to the prediction), so that it carries the flux the film
//  moved with.
//
//  The explicit part A phi_old is formed from the matrix of the stage being
//  solved, not from the matrix of the step before, which would put the
//  operator a step behind.
//-----------------------------------------------------------------------
int CZ_HeightEquation_2D::SolveStage(int &iter, double &resid)
{
	{ YUZEN_PROF("h.assemble"); SetMatrixA(); }
	_crsA.MatVec(_Ax_old, _phi_old);
	if(Math::IsZero(_theta)) return SolveExplicit();
	{ YUZEN_PROF("h.rhs"); MakeMatrix_IA(); MakeRHS(); }
	int ch;
	{ YUZEN_PROF("h.linsolve"); ch = _bcg.Solve(_crsIA, _Rhs, _phi_new, iter, resid); }
	if(ch != 0)
	{
		Global::Warning<<"Semi-Implicit: BiCGStab did not converge: ret = "<< ch << ", iter = "<<iter<< ", resid = "<< std::scientific <<resid<< nl;
		_failCause = "the film-height solve did not converge (Numerics.HeightEquation.LinearSolver)";
	}
	return ch;
}

int CZ_HeightEquation_2D::TimeAdvance(int &iter, double &resid)
{
	const int r = C_SemiImplicitEquation::TimeAdvance(iter, resid);
	//  With the step fixed there is no shorter one to retry at, and a step
	//  whose film or concentration solve has failed is not taken as it
	//  stands: the run stops.
	if(r != 0 && !_adjustTimeStep) FailedFixedStep("ConcZ");
	return r;
}

//  The concentration at the midpoint of this step, extrapolated half a step
//  on from the starts of this step and the one before. False (cm untouched)
//  at the first step, or without a concentration field: the concentration
//  of the field is then used as it stands.
bool CZ_HeightEquation_2D::MidpointConc(Vector<double> &cm) const
{
	if(!(swConcField && _ccPrevSet)) return false;
	const double r = 0.5 * _timeStep->DT().now / _dtPrev;
	cm = _concEq->GetX();
	for(int i = 0; i < cm.size(); i++) cm[i] += r * (cm[i] - _ccPrev[i]);
	return true;
}

int CZ_HeightEquation_2D::Solve(int &iter, double &resid)
{
	assert(_initialized);
	int ch = SolveStage(iter, resid);
	if(ch != 0) return ch;

	//  The midpoint state: h halfway between the start and the prediction,
	//  and the concentration extrapolated half a step on from the last two
	//  (the concentration of this step is solved after the film, so there is
	//  no prediction of it yet; at the first step it is left where it is).
	Vector<double> mid = _phi_old;
	mid += _phi_new;
	mid *= 0.5;
	_field->Update_H(mid);
	Vector<double> cm;
	const bool extrapolate = MidpointConc(cm);
	if(extrapolate) _field->Update_C(cm);

	int iter2 = 0;
	double resid2 = 0;
	ch = SolveStage(iter2, resid2);
	iter += iter2;
	resid = resid2;

	//  Back to the start of the step; AfterTimeStepTrial writes the end.
	_field->Update_H(_phi_old);
	if(extrapolate)
	{
		Vector<double> c0 = _concEq->GetX();
		_field->Update_C(c0);
	}
	return ch;
}

double CZ_HeightEquation_2D::CheckTimeStepCriterion()
{
	return LocalTruncationCriterion(_phi_new);
}
