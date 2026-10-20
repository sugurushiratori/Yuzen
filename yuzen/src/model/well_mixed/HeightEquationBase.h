/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  HeightEquationBase.h : the dimension-independent part of the height
*                         equation driver
*
*  The 1-D and 2-D film-height equations (WM1D/HeightEquation_1D.cpp,
*  WM2D/HeightEquation_2D.cpp) differ in two respects only:
*
*    (a) the types of the grid, the site, the field and the boundary
*        condition, and the map from a site to its row number, which the
*        traits of WellMixedTraits.h supply
*    (b) the body of SetVelocityBC, the boundary condition of the velocity
*        field, which is left pure virtual because the expression itself
*        differs with the dimension
*
*  Everything else is kept here as a single copy, the source, the
*  evaporation and the concentration clamps included.
*
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _HEIGHT_EQUATION_BASE_INCLUDED_
#define _HEIGHT_EQUATION_BASE_INCLUDED_

#include "Opt.h"
#include "ManufacturedSource.h"
#include "Profile.h"
#include "YuzenUtil.h"
#include "NewtonKantorovich.h"
#include "WellMixedTraits.h"
#include "BoundaryCondition.h"
#include "NonDimensionalNumbers.h"
#include "FluidProperty.h"
#include "ConcStep.h"

#include <algorithm>
#include <vector>

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class C_Fluid;

template<class T>
class WellMixed_HeightEquationBase : public C_NewtonKantorovich
{
public:
	typedef typename T::Grid     Grid;
	typedef typename T::Site     Site;
	typedef typename T::Field    Field;
	typedef typename T::Boundary Boundary;
	typedef typename T::ConcEq   ConcEq;
	typedef typename T::VelField VelField;

	/** @brief Clear every pointer; Build() supplies the real objects. */
	WellMixed_HeightEquationBase()
		: C_NewtonKantorovich(),
		  G(0), _Nx(1), _field(0), _field_guess(0), _prop(0), _NonDim(0), _boundary(0),
		  _swConcField(false), _concEq(0),
		  _bcAppliedToA(false), _bcAppliedToJacobi(false)
	{ _stepDt = 0; _sideOld = _sideOldResin = _sideStep = _sideStepResin = 0; _excessMoved = _excessTrimmed = 0;
	  _openEdges = false; _edgeOld = _edgeStep = _edgeStepResin = 0;
	  _deficitMoved = _deficitAdded = 0; _minStart = 0; }

	virtual ~WellMixed_HeightEquationBase(){}

	//  Build every physics term and line the active ones up in _physics.
	//  concEq receives the concentration equation (for the coupling).
	//  fieldGuess is the intermediate field of the Newton iteration, and it
	//  is the one the physics terms read.
	/**
	 * @brief Build every physics term, keep the active ones, and set up the
	 *        implicit solver and its convergence criteria.
	 * @param inp     the whole input file; the blocks this equation needs are read from it
	 * @param prop    fluid properties (viscosity and the concentration dependence)
	 * @param ts      time-step object supplying dt and the current time
	 * @param concEq     the concentration equation, kept for the coupling
	 * @param ff         the committed field
	 * @param fieldGuess the working field of the Newton iteration; this is the
	 *                   one the physics terms read
	 * @param nonDim  non-dimensional numbers of the problem
	 * @param bc      boundary condition applied to the rows of this equation
	 */
	void Build(Input &inp, C_Fluid *prop, C_TimeSteps *ts, ConcEq *concEq,
	           Field *ff, Field *fieldGuess, NonDimensionalNumbers *nonDim, Boundary *bc);

	//  Call SetMatrixA / SetJacobiMatrix on each active physics term in
	//  turn. The order in which the contributions are added is the order of
	//  _physics, and it fixes the order of the floating-point additions, so
	//  changing that order changes the results bit for bit.
	/**
	 * @brief Rebuild the matrix A for the current guess, then fill B with the
	 *        evaporation and source terms.
	 *
	 * The boundary-condition rows are written once and preserved on every
	 * later call. Evaporation is capped so that no step removes more than
	 * half the local film thickness, nor more than the solvent it holds.
	 */
	virtual void SetMatrixA();
	/**
	 * @brief Copy A into the Jacobian and let each active term add its own
	 *        derivative on top.
	 */
	virtual void SetMatrixJacobi();

	//--------------------------------------------------------------
	//  Coupling with the concentration. These are the hooks of the Newton
	//  iteration (see NewtonKantorovich.h).
	//--------------------------------------------------------------
	/**
	 * @brief Record the concentration the step starts from, as the reference
	 *        against which later iterations are measured.
	 */
	virtual void BeforeNKLoop();
	//  Every iteration: update the field with the guessed film height ->
	//  solve the concentration equation once -> put the resulting correction
	//  into _coupledCorrection. That correction enters the convergence test,
	//  so what the iteration converges to is the solution of the height and
	//  the concentration solved fully implicitly as one system.
	/**
	 * @brief One outer iteration: update the field with the current height
	 *        guess, solve the concentration equation once, and report how far
	 *        the concentration moved.
	 * @return false if the concentration solve failed, which aborts the step
	 */
	virtual bool EveryNKLoop();
	/**
	 * @brief Commit the step: fix the concentration side as well, without
	 *        solving it again, and copy the working field back.
	 */
	virtual bool AfterNKLoop();

	/**
	 * @brief The quantity the adaptive time step is judged on.
	 * @return the local truncation criterion evaluated on the height guess
	 */
	virtual double CheckTimeStepCriterion() { return LocalTruncationCriterion(_phi_guess); }

	/** @brief Mark the ghost rows as boundary-condition rows. */
	virtual void SetFlagBC();

	/**
	 * @brief Have every active term write the velocity it produces, then
	 *        combine them into the mixed velocity.
	 *
	 * The per-term boundary condition is applied to each velocity field and
	 * to the mixed one.
	 */
	void CalVelocity();

protected:
	/**
	 * @brief Apply the boundary condition of a velocity field.
	 *
	 * The expression itself differs with the dimension, so this is
	 * implemented by the derived class.
	 */
	virtual void SetVelocityBC(VelField &vv) = 0;

	/**
	 * @brief Dimension-specific initialization performed once the grid has
	 *        been received (the 2-D case picks up _Ny here).
	 */
	virtual void BuildGridExtent() {}

	/** @brief Map a site to its row number. */
	int Index(const Site &s) const { return T::Index(G, s); }

	/**
	 * @brief Cap the concentration of the working field at 1.
	 *
	 * Common to both dimensions, in EveryNKLoop and AfterNKLoop. It acts on
	 * the field the coefficients (viscosity, Marangoni, evaporation) are
	 * taken from, not on the state of the concentration equation, so it
	 * takes no solute away; the committed state is brought to c <= 1 by
	 * RedistributeExcess.
	 */
	void ClampConc();

	/**
	 * @brief Bring the concentration of the step's end state into
	 *        0 <= c <= 1 without taking solute away or making any.
	 *
	 * The solute a real cell holds above c = 1, h (c - 1), goes to its
	 * face neighbors that have room, h (1 - c), in proportion to that room;
	 * what they cannot take stays and goes on in the next pass, from the
	 * cells it then sits in. After a few passes what is left is spread over
	 * every real cell with room, in proportion to it. Only a film holding
	 * more solute than liquid in all is left with any: that is cut, and
	 * counted (ExcessTrimmed). A cell below c = 0 is filled the same way,
	 * the other way round: its deficit, -h c, is taken from its face
	 * neighbors in proportion to the solute h c they hold, then from the
	 * whole film; only a film with less than no solute in all keeps a
	 * remainder, which is made up and counted (DeficitAdded). The sum of
	 * h c over the real cells is kept to rounding; where no real cell is
	 * outside [0, 1] nothing is touched. Ghost rows outside [0, 1] are set
	 * to the bound (they hold no solute of the film).
	 *
	 * @param cx  the state of the concentration equation, row by row
	 */
	void RedistributeExcess(Vector<double> &cx);
	/**
	 * @brief The passes of RedistributeExcess, on the real cells: above
	 *        c = 1 (lower false) or below c = 0 (lower true).
	 */
	void MoveExcess(Vector<double> &cx, bool lower);
public:
	/**
	 * @brief The solute, h (c - 1) summed over the real cells above 1, that
	 *        the step last committed moved to cells with room
	 *        (RedistributeExcess), times the cell measure.
	 */
	double ExcessMoved(double cellMeasure) const { return _excessMoved * cellMeasure; }
	/**
	 * @brief The solute that step could not place and cut (a film with more
	 *        solute than liquid in all, or a cell above 1 without film),
	 *        times the cell measure. It is missing from the solute balance
	 *        (ResinVolErr).
	 */
	double ExcessTrimmed(double cellMeasure) const { return _excessTrimmed * cellMeasure; }
	/**
	 * @brief The solute, -h c summed over the real cells below 0, that the
	 *        step last committed took from cells holding solute to fill
	 *        them (RedistributeExcess), times the cell measure.
	 */
	double DeficitMoved(double cellMeasure) const { return _deficitMoved * cellMeasure; }
	/**
	 * @brief The solute that step had to make up at c = 0 (a film with less
	 *        than no solute in all, or a cell below 0 without film), times
	 *        the cell measure. It is an error of the solute balance
	 *        (ResinVolErr).
	 */
	double DeficitAdded(double cellMeasure) const { return _deficitAdded * cellMeasure; }
	/**
	 * @brief The lowest concentration of a real cell at the start of the
	 *        step last solved (0 if none was below). The divergence test
	 *        looks at it as well as at the end of the step: a field file
	 *        read with a cell far below 0 must stop the run at its first
	 *        step, and that step may mix the cell back into range by itself
	 *        (or RedistributeExcess fill it). The concentration a solve
	 *        gives before the fill is not looked at: at a sharp front an
	 *        unbounded scheme undershoots by more than the slack within one
	 *        step, and the fill is what keeps that from growing.
	 */
	double MinStartConc() const { return _minStart; }
	/**
	 * @brief Whether a real cell of the film h has no liquid in it (h <= 0,
	 *        or not a number). The balance of h c means nothing there: the
	 *        concentration is not solved with such a film, and the film
	 *        itself is left to the divergence test of the model.
	 */
	bool FilmNotPositive(const Vector<double> &h) const
	{
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			const int i = Index(s);
			if( ! _crsA.IsBC(i) && !(h[i] > 0.0) ) return true;
		}
		return false;
	}
	/**
	 * @brief The power of the face fluxes against a potential (Output.Energy):
	 *        sum over the interior rows of phi times the rate the fluxes give
	 *        the row, -(the difference of its face fluxes) - its sideways
	 *        outflow, in the state f; per unit cell measure. phi is indexed
	 *        by the grid's linear index (G->imin() .. G->imax()). With phi the
	 *        potential, -dF/dH, this is the dissipation -dF/dt of the scheme
	 *        itself (EnergyBudget.h).
	 */
	double FluxPower(const Field &f, const std::vector<double> &phi)
	{
		std::vector<double> F[T::nAxis];
		Vector<double> r;
		FaceVolumeFluxes(f, F, r);
		double p = 0;
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			const int i = Index(s);
			if( _crsA.IsBC(i) ) continue;
			double div = 0;
			for(int a = 0; a < T::nAxis; ++a)
				div += F[a][(size_t)i] - F[a][(size_t)Index(T::Adj(G, s, T::Dir(2*a)))];
			p += phi[(size_t)ii] * (-div - r[i]);
		}
		return p;
	}
	/**
	 * @brief The volume the nozzle supplied in the step last solved, as the
	 *        film-height equation took it: the source of every interior
	 *        row, at the time of the middle of that step, times the cell
	 *        measure and that step's dt.
	 */
	double SuppliedVolume(double cellMeasure) const
	{
		if(!_Source.Enabled()) return 0;
		double q = 0;
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			if( ! _crsA.IsBC(Index(s)) ) q += _Source(s);
		}
		return q * cellMeasure * _stepDt;
	}
	/**
	 * @brief The solute that came with the volume the nozzle supplied in the
	 *        step last solved (SuppliedVolume): each row's supply times the
	 *        concentration the liquid arrived at.
	 *
	 * With Physics.Source.ConcType = "local" the supplied liquid arrives at
	 * the concentration it lands on: the balance of h c the concentration
	 * equation solves takes the supply S of the film-height equation with
	 * S (c_0 + c_1)/2 of solute (ConcStep.h), c_1 the concentration that
	 * solve gave (before RedistributeExcess, which moves solute between cells
	 * only). With "supplied" it arrives at Physics.Source.Conc, and the
	 * balance takes S Conc, exactly what is counted here. Without the
	 * concentration equation c does not change, and the liquid is counted at
	 * the c it lands on whatever ConcType says.
	 */
	double SuppliedResin(double cellMeasure) const
	{
		if(!_Source.Enabled()) return 0;
		if(_Source.ConcSupplied() && _swConcField)
		{
			const double cs = _Source.Conc();
			double q = 0;
			for(int ii = G->imin(); ii <= G->imax(); ii++)
			{
				Site s = G->template Of<Site>(ii);
				if( ! _crsA.IsBC(Index(s)) ) q += _Source(s) * cs;
			}
			return q * cellMeasure * _stepDt;
		}
		const bool start = (_cStepStart.size() == _nRows);
		const bool end   = (_cStepEnd.size() == _nRows);
		double q = 0;
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			const int i = Index(s);
			if( _crsA.IsBC(i) ) continue;
			const double c1 = end ? _cStepEnd[i] : _field->CC(s);
			const double c0 = start ? _cStepStart[i] : c1;
			q += _Source(s) * 0.5 * (c0 + c1);
		}
		return q * cellMeasure * _stepDt;
	}
	/**
	 * @brief The evaporation the film-height equation took out in the step
	 *        last solved, as a rate (volume per time): the evaporation of
	 *        every interior row as B was assembled with it, times the cell
	 *        measure.
	 */
	double EvaporationRate(double cellMeasure) const
	{
		if(_evapB.size() != _nRows) return 0;
		double q = 0;
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			if( ! _crsA.IsBC(Index(s)) ) q += _evapB[Index(s)];
		}
		return q * cellMeasure;
	}
	/**
	 * @brief The liquid the 1-D centrifugal term flung sideways, into the y
	 *        that is not solved for, in the step last solved, as the
	 *        film-height equation took it out: its theta-weighted average of
	 *        the rate at the start of the step and at the film the step
	 *        ended on, times the cell measure and that step's dt (StepDt).
	 *        *resin receives the solute that liquid carried. 0 in 2-D and
	 *        without the centrifugal term.
	 */
	double SideOutflowVolume(double cellMeasure, double *resin) const
	{
		if(resin) *resin = _sideStepResin * cellMeasure * _stepDt;
		return _sideStep * cellMeasure * _stepDt;
	}
	/**
	 * @brief The liquid that left through the open edges (FreshInlet,
	 *        Outlet) in the step last solved, as both equations took it: the
	 *        theta-weighted volume flux of the film-height equation through
	 *        every edge face of an open side, outward (the face fluxes of
	 *        the film at the start of the step and of the film it ended on),
	 *        times the cell measure and that step's dt. *resin receives the
	 *        solute the concentration equation put through those faces in
	 *        the same step (WellMixed_ConcEquationBase::EdgeOutflowRate).
	 *        Negative for what came in. 0 without an open edge.
	 */
	double EdgeOutflowVolume(double cellMeasure, double *resin) const
	{
		if(resin) *resin = _edgeStepResin * cellMeasure * _stepDt;
		return _edgeStep * cellMeasure * _stepDt;
	}
	/** @brief The dt of the step last solved (see _stepDt). */
	double StepDt() const { return _stepDt; }
	/** @brief The input key of the face interpolation of the upwinded fluxes. */
	std::string ConvectionKey() const { return C_Interpolation_base4::Key(_ParallelGravity.Scheme().Scheme()); }
	/** @brief Whether a limited scheme carries an upwinded flux that is switched on. */
	bool LimitedUpwind() const
	{
		return _ParallelGravity.Scheme().IsLimited()
		    && (_ParallelGravity.Enabled() || _Centrifugal.Enabled() || _SolutoCapillarity.Enabled());
	}
protected:

	//------------------------------------------------------------------
	//  Return how far the concentration moved during one outer iteration,
	//  and update the reference value to the current one. Used by the
	//  convergence test.
	//  The measurement is taken on the values after the cap (ClampConc) has
	//  been applied. Those are the values that enter the coefficients of the
	//  physics terms, and they are also what decides whether the coupling
	//  has settled.
	//------------------------------------------------------------------
	double ConcCorrection();

	Vector<double> _concPrev;	// the concentration at the previous outer iteration
	//  The concentration iterate as the equation solved it, before ClampConc
	//  caps the field: what a limited scheme freezes its limiter at
	//  (ConcEq::LinkLimiterState). Equal to the field wherever c <= 1.
	Vector<double> _cIterate;
	/** @brief _cIterate from a field, row by row (ghosts included). */
	void SetIterate(const Field &f);

	// Physical models
	typename T::LaplacePressure   _LaplacePressure;
	typename T::VerticalGravity   _VerticalGravity;
	typename T::ParallelGravity   _ParallelGravity;
	typename T::Centrifugal       _Centrifugal;
	typename T::SolutoCapillarity _SolutoCapillarity;
	typename T::DisjoiningPressure _DisjoiningPressure;
	typename T::Source            _Source;

	//  The part of A of the terms that do not evaluate A h themselves
	//  (PhysicsTermBase::ConservativeForm): ApplyA multiplies it out and
	//  adds the others' conservative A h
	CRS_Matrix _crsArest;
	virtual void ApplyA(Vector<double> &y, const Vector<double> &x) const;

	//  The source of the method of manufactured solutions (verification only)
	C_ManufacturedSource _mms;
	//  The evaporation at the start of the step, the explicit half of its
	//  Crank-Nicolson average (see SetMatrixA)
	Vector<double> _evapOld;
	//  The dt the source term B was last assembled with: that of the step
	//  last solved, which is DT().now until the clock advances (++) and
	//  DT().old after it.
	double _stepDt;
	//  The sideways outflow of the 1-D centrifugal term (SideOutflowVolume):
	//  its rate (and that of the solute) at the start of the step, the
	//  explicit half A(h^n) h^n takes, and the theta-weighted rate of the
	//  step last solved
	double _sideOld, _sideOldResin, _sideStep, _sideStepResin;
	//  Whether an edge is open (FreshInlet or Outlet); the volume flux out
	//  through the open edges at the start of the step, the explicit half
	//  of the step (EdgeFluxRate); and the theta-weighted rates of the step
	//  last solved, of the liquid and of the solute (EdgeOutflowVolume)
	bool _openEdges;
	double _edgeOld, _edgeStep, _edgeStepResin;
	/**
	 * @brief The volume flux out through the edge faces of the open sides
	 *        in the state f, per unit cell measure: the sum over the terms of
	 *        their FaceFlux, as FaceVolumeFluxes forms the flux of every
	 *        face, so that it is the flux the step's fluxes have there
	 *        (ConsistentFluxes does not touch an edge face).
	 */
	double EdgeFluxRate(const Field &f) const;
	//  The solute (per unit cell measure) the last commit moved from cells
	//  above c = 1 to cells with room, and what it had to cut; and what it
	//  moved into cells below c = 0, and had to make up (RedistributeExcess)
	double _excessMoved, _excessTrimmed;
	double _deficitMoved, _deficitAdded;
	//  The lowest c of a real cell at the start of the step (MinStartConc)
	double _minStart;
	//  The concentration at the start of the step (with the nozzle on), and
	//  the one the committed solve gave, for the solute the nozzle supplies
	//  (SuppliedResin)
	Vector<double> _cStepStart, _cStepEnd;
	//  The evaporation B was last assembled with (SetMatrixA): that of the
	//  solve that gave the film of the step, what EvaporationRate reports
	//  and what the concentration equation completes the step with
	Vector<double> _evapB;

	//------------------------------------------------------------------
	//  The film side of the concentration equation (ConcStep.h)
	//------------------------------------------------------------------
	//  The volume flux of every face, once per face, from the physics terms'
	//  FaceFlux in the state f: F[a][i] the face between row i and its
	//  neighbor on the high side along axis a (the low cell's own face for
	//  a real cell, the real cell's for a ghost below it); r the sideways
	//  outflow of every real row
	void FaceVolumeFluxes(const Field &f, std::vector<double> *F, Vector<double> &r);
	//  The evaporation (Crank-Nicolson average, as B takes it) and the
	//  supply (nozzle and manufactured source, at the middle of the step) of
	//  every real row, with the film f at the end
	void StepSources(const Field &f);
	//  The face fluxes of the end of the step made those of the step the
	//  film took (see the implementation): F1 += phi / theta, phi running
	//  sums through the inner faces whose divergence is what the fluxes
	//  miss of h0 - h1 (the residual of the step at the film it ended on,
	//  less its mean, which is the volume error of the step and stays as
	//  that)
	void ConsistentFluxes(double dt);
	std::vector<double> _F0[T::nAxis], _F1[T::nAxis];
	Vector<double> _r0, _r1, _srcStep, _evpStep;
	WMConcStep _cstep;
	//  The evaporation rate at a site, capped so that no step removes more
	//  than half of the local film, nor more than the solvent there is
	double EvapAt(const Field &f, Site &s, double dt) const;

	//------------------------------------------------------------------
	//  The list of pointers used to iterate over the physics terms as a
	//  group. The objects themselves remain the members above (ownership is
	//  not transferred).
	//  ** Do not change the order of this list. ** Adding into a CRS_Matrix
	//  with += creates the columns in insertion order, so changing the order
	//  changes the order in which the entries are stored, which changes the
	//  rounding of the summations in ScalingDiag / Residual.
	//------------------------------------------------------------------
	std::vector<typename T::Physics*> _physics;
	/**
	 * @brief Build the list of pointers above.
	 *
	 * The order of the list is itself meaningful, so a new term is appended
	 * at the end.
	 */
	void RegisterPhysics();

	Grid* G;
	int   _Nx;

	Field *_field;
	Field *_field_guess;

	C_Fluid *_prop;
	NonDimensionalNumbers *_NonDim;
	Boundary *_boundary;

	bool    _swConcField;
	ConcEq *_concEq;

	//  Whether the boundary condition has already been written into A. These
	//  are held per instance; making them function-local statics would stop
	//  several cases from being run within one process.
	bool _bcAppliedToA;
	bool _bcAppliedToJacobi;

private:
	WellMixed_HeightEquationBase(const WellMixed_HeightEquationBase&);
	void operator=(const WellMixed_HeightEquationBase&);
};


//==================================================================
//  Implementation
//==================================================================

template<class T>
void WellMixed_HeightEquationBase<T>::RegisterPhysics()
{
	_physics.clear();
	_physics.push_back(&_LaplacePressure);
	_physics.push_back(&_VerticalGravity);
	_physics.push_back(&_ParallelGravity);
	_physics.push_back(&_SolutoCapillarity);
	_physics.push_back(&_Centrifugal);
	//  Appended at the end. When the term is disabled its SetMatrixA is not
	//  called at all, so the results of a computation that does not use the
	//  disjoining pressure do not change by a single bit.
	_physics.push_back(&_DisjoiningPressure);
}


template<class T>
void WellMixed_HeightEquationBase<T>::Build(Input &inp, C_Fluid *prop, C_TimeSteps *ts, ConcEq *concEq,
		Field *ff, Field *fieldGuess, NonDimensionalNumbers *nonDim, Boundary *bc)
{
	//Set up the pointers
	_timeStep = ts;
	_field = ff;
	_field_guess = fieldGuess;
	_prop = prop;
	_boundary = bc;
	_NonDim = nonDim;
	_concEq = concEq;
	_swConcField = (inp.Physics.ConcField.Model.compare("WM") == 0);
	if(_swConcField && _concEq) { _concEq->LinkStep(&_cstep); _concEq->LinkLimiterState(&_cIterate); }

	//Set up everything related to the grid
	G = ff->GetGridPtr();
	_Nx = G->nx();
	BuildGridExtent();

	_VecB = 0;

	//Build the physical models
	_LaplacePressure.Build(inp.Physics, G, _NonDim, _boundary);
	_VerticalGravity.Build(inp.Physics, G, _NonDim, _boundary);
	_ParallelGravity.Build(inp.Physics, G, _NonDim, _boundary);
	_Centrifugal.Build(inp.Physics, G, _NonDim, _boundary);
	_SolutoCapillarity.Build(inp.Physics, G, _NonDim, _boundary);
	_DisjoiningPressure.Build(inp.Physics, G, _NonDim, _boundary);
	//  The face interpolation of the upwinded fluxes (Numerics.Convection.Height)
	{
		const C_Interpolation_base4::t_Scheme hs =
			C_Interpolation_base4::FromKeyOr(inp.Numerics.Convection.Height, C_Interpolation_base4::Koren);
		_ParallelGravity.SetConvectionScheme(hs);
		_Centrifugal.SetConvectionScheme(hs);
		_SolutoCapillarity.SetConvectionScheme(hs);
	}
	_Source.Build(inp.Physics, G, _NonDim, _timeStep);
	_mms.Build(inp, *_NonDim, T::nAxis == 2);

	RegisterPhysics();

	_openEdges = false;
	for(int d = 0; d < T::nDir; ++d)
		_openEdges = _openEdges || BoundaryCondition::IsOpen(_boundary->TypeOfSide(d));

	//Build the equation solver.
	//This function calls SetMatrix internally, so it must come after every
	//model has been built.
	C_ImplicitSolver::Build(ff->HH(), *ts, C_ImplicitSolver::CrankNicholson, inp.Numerics.TimeStepControl.Enable);

	//------------------------
	//  Everything related to the linear solver
	//------------------------
	//  Set the preconditioner, the scaling and the convergence criteria in
	//  one go. What fixes the shape of the linear system is the in-plane
	//  dimension: T::nAxis is 1 in one dimension and 2 in two. The "auto"
	//  setting chooses the backend from this.
	this->_bcg.Configure(inp.Numerics.HeightEquation.LinearSolver, (int)T::nAxis);

	// Set the convergence criteria of the Newton-Kantorovich iteration
	this->SetMaxIterNK(inp.Numerics.HeightEquation.NewtonMaxIterations);
	this->SetMinIterNK(inp.Numerics.HeightEquation.NewtonMinIterations);
	this->SetEpsilonNK(inp.Numerics.HeightEquation.NewtonEpsilon);
	this->SetEpsilonRelNK(inp.Numerics.HeightEquation.NewtonEpsilonRel);
	this->SetJacobiInterval(inp.Numerics.HeightEquation.JacobianUpdateInterval);

	if(inp.Numerics.TimeStepControl.Enable)
	{
		this->SetAdjustTimeChangeFactor(inp.Numerics.TimeStepControl.ChangeFactor);
		this->SetAdjustTimeCriterion(inp.Numerics.TimeStepControl.Criterion);
		this->SetAdjustTimeScaleFloor(inp.Numerics.TimeStepControl.ScaleFloor);
		this->SetAdjustTimeKeepSteps(inp.Numerics.TimeStepControl.KeepSteps);
	}
}


template<class T>
void WellMixed_HeightEquationBase<T>::SetFlagBC()
{
	//Set the flags that mark the boundary-condition rows
	for(Site s = G->template First<Site>(); s != G->template Last<Site>(); ++s)
	{
		if(G->CheckGhost(s) != Grid::Inside)
			this->SetAsBoundary(Index(s));
	}
}


template<class T>
void WellMixed_HeightEquationBase<T>::SetMatrixA()
{
	int jj;

	if(!_bcAppliedToA)
	{
		//Set the boundary-condition part of the matrix (once, at initialization)
		_boundary->Apply(_crsA, _Rhs, BoundaryCondition::Height, *_field);
		_bcAppliedToA = true;
	}

	//--------------------
	// Clear the matrix (the boundary-condition part is kept)
	//--------------------
	//  A term that evaluates A h itself (ConservativeForm) goes into A
	//  alone; the others go into _crsArest, the part of A that ApplyA
	//  multiplies out, and are added to A after it.
	const bool fused = yuzen::opt::On(yuzen::opt::Regions);
	if(fused)
	{
		//  opt regions (util/Opt.h): both cleared in one parallel loop
		if(_crsArest.RowSize() != _crsA.RowSize()) _crsArest = _crsA;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(jj) if(yuzen::opt::Par(_crsA.RowSize()))
#endif
		for(int irow = 0; irow < _crsA.RowSize(); irow++)
		{
			if( _crsA.IsBC(irow) ) continue;
			for(jj = 0; jj < _crsA.NonZero(irow); jj++) _crsA.Element(irow, jj) = 0;
			for(jj = 0; jj < _crsArest.NonZero(irow); jj++) _crsArest.Element(irow, jj) = 0;
		}
	}
	else
	{
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(jj) if(yuzen::opt::Par(_crsA.RowSize()))
#endif
	for(int irow = 0; irow < _crsA.RowSize(); irow++)
	{
		if( ! _crsA.IsBC(irow) )
		{
			for(jj = 0; jj < _crsA.NonZero(irow); jj++)
				_crsA.Element(irow, jj) = 0;
		}
	}

	if(_crsArest.RowSize() != _crsA.RowSize()) _crsArest = _crsA;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(jj) if(yuzen::opt::Par(_crsArest.RowSize()))
#endif
	for(int irow = 0; irow < _crsArest.RowSize(); irow++)
	{
		if( ! _crsArest.IsBC(irow) )
		{
			for(jj = 0; jj < _crsArest.NonZero(irow); jj++)
				_crsArest.Element(irow, jj) = 0;
		}
	}
	}

	//--------------------------------------
	// Add the contribution of each physical model that is in use
	//--------------------------------------
	for(size_t k = 0; k < _physics.size(); ++k)
		if(_physics[k]->Enabled())
			_physics[k]->SetMatrixA(_physics[k]->ConservativeForm() ? _crsA : _crsArest, *_field_guess);
	//  Row by row: each entry of A takes one addition from its row of
	//  _crsArest, so the rows may go in any order (in parallel with
	//  opt regions, util/Opt.h)
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(jj) if(yuzen::opt::ParNew(_crsArest.RowSize()))
#endif
	for(int irow = 0; irow < _crsArest.RowSize(); irow++)
	{
		if( _crsArest.IsBC(irow) ) continue;
		for(jj = 0; jj < _crsArest.NonZero(irow); jj++)
		{
			const double v = _crsArest.Element(irow, jj);
			if(jj == 0) _crsA.Diag(irow) += v;
			else        _crsA(irow, _crsArest.Jcol(irow, jj)) += v;
		}
	}

	//--------------------
	// Set the source term
	//--------------------
	//  B is taken as a Crank-Nicolson average like the terms in A: the
	//  evaporation at the start of the step (_evapOld, set in BeforeNKLoop)
	//  and at the current guess of its end, weighted 1 - theta and theta.
	//  The start of the step alone would make a drying film first order in
	//  dt however well the rest is integrated. The liquid supply
	//  and the manufactured source are evaluated at the middle of the step.
	Site s;
	double evp, src, xx, yy;
	const double dt = _timeStep->DT().now;
	const double tmid = _timeStep->NowTime() + 0.5 * dt;
	_stepDt = dt;
	if(_evapOld.size() != _nRows) { _evapOld.resize(_nRows); _evapOld = 0; }
	if(_evapB.size() != _nRows) { _evapB.resize(_nRows); _evapB = 0; }
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(s,evp,src,xx,yy) if(yuzen::opt::Par((G->imax() - G->imin() + 1)))
#endif
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		s = G->template Of<Site>(ii);
		if( ! _crsA.IsBC(Index(s)) )
		{
			evp = (1.0 - _theta) * _evapOld[Index(s)] + _theta * EvapAt(*_field_guess, s, dt);
			_evapB[Index(s)] = evp;
			src = _Source(s);
			if(_mms.Enabled())
			{
				T::SiteXY(G, s, xx, yy);
				src += _mms.Snd(xx, yy, tmid);
			}
			_VecB(Index(s)) = evp -src;
		}
	}
}

//------------------------------------------------------------------
//  A x for the explicit parts of the Newton iteration, x the film height A
//  was last assembled at: the rows of the boundary conditions and the terms
//  that have no conservative form as the product with their entries, the
//  others (the capillary pressure and the vertical gravity) face by face.
//  The volume a step adds is then the supply, the evaporation and the
//  sideways outflow to rounding of the fluxes; the product of A with h
//  would be off by eps |A| |h| per cell, which on a thick film over small
//  cells is not small: VolErr -5.4e-5 in nozzle_1d_rotating, the film
//  1.2 mm thick over cells 3.3 um wide (tools/analytic.py nozzle_balance).
//  The upwinded terms stay products: their entries are of order D/dx, not
//  D/dx^4, and what they leave is at rounding of the film.
//------------------------------------------------------------------
template<class T>
void WellMixed_HeightEquationBase<T>::ApplyA(Vector<double> &y, const Vector<double> &x) const
{
	bool any = false;
	for(size_t k = 0; k < _physics.size(); ++k)
		if(_physics[k]->Enabled() && _physics[k]->ConservativeForm()) any = true;
	if(!any || _crsArest.RowSize() != _crsA.RowSize()) { _crsA.MatVec(y, x); return; }

	_crsArest.MatVec(y, x);
	//  The boundary rows as A has them
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::ParNew(_crsA.RowSize()))
#endif
	for(int irow = 0; irow < _crsA.RowSize(); irow++)
	{
		if( ! _crsA.IsBC(irow) ) continue;
		double v = 0;
		for(int jj = 0; jj < _crsA.NonZero(irow); jj++)
			v += _crsA.Element(irow, jj) * x[_crsA.Jcol(irow, jj)];
		y[irow] = v;
	}
	for(size_t k = 0; k < _physics.size(); ++k)
		if(_physics[k]->Enabled() && _physics[k]->ConservativeForm())
			_physics[k]->AddConservative(y, x, _crsA, *_field_guess);
}

template<class T>
double WellMixed_HeightEquationBase<T>::EvapAt(const Field &f, Site &s, double dt) const
{
	double evp = _NonDim->Ev() * _NonDim->EvapRate(const_cast<Field&>(f), s);
	if(evp <= 0) return evp;
	//  Never evaporate more than half of the film thickness in one step. A
	//  guess of the Newton iteration can pass through h < 0; the cap is then
	//  0, not 0.5 h / dt, which would be a source of liquid (and one that
	//  would act even with the evaporation switched off).
	const double hh = f.HH(s);
	if(evp * dt > 0.5 * hh) evp = (hh > 0) ? 0.5*hh/dt : 0.0;
	//  Nor more than the solvent there is (C_Evaporation::SolventLimit):
	//  each end of the Crank-Nicolson step taking at most h (1 - c) / dt,
	//  the solvent left at the end is at least a third of that at the
	//  start, so c stays below 1 and a film dries toward it without
	//  passing it, and nothing beyond rounding is left above 1 for
	//  RedistributeExcess to move (tools/proptest.py seed 1206).
	return C_Evaporation::SolventLimit(evp, hh, f.CC(s), dt);
}

//------------------------------------------------------------------
//  The volume flux of every face, once (see the declaration). Each slot is
//  written by one real cell only, so the loop runs in parallel and the
//  result does not depend on the number of threads.
//------------------------------------------------------------------
template<class T>
void WellMixed_HeightEquationBase<T>::FaceVolumeFluxes(const Field &f, std::vector<double> *F, Vector<double> &r)
{
	for(int a = 0; a < T::nAxis; ++a) F[a].assign((size_t)_nRows, 0.0);
	if(r.size() != _nRows) r.resize(_nRows);
	r = 0;
	std::vector<Site> svec(T::nStencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel firstprivate(svec)
	{
#pragma omp for
#endif
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		Site s = G->template Of<Site>(ii);
		const int i = Index(s);
		if( _crsA.IsBC(i) ) continue;
		G->SetSite(svec, s);
		for(int a = 0; a < T::nAxis; ++a)
		{
			const int il = Index(T::Adj(G, s, T::Dir(2*a)));
			const bool ghostBelow = _crsA.IsBC(il);
			double hi = 0, lo = 0;
			for(size_t k = 0; k < _physics.size(); ++k)
			{
				if(!_physics[k]->Enabled()) continue;
				hi += _physics[k]->FaceFlux(s, 2*a + 1, f, svec);
				if(ghostBelow) lo += _physics[k]->FaceFlux(s, 2*a, f, svec);
			}
			F[a][(size_t)i] = hi;
			if(ghostBelow) F[a][(size_t)il] = lo;
		}
		double rr = 0;
		for(size_t k = 0; k < _physics.size(); ++k)
			if(_physics[k]->Enabled()) rr += _physics[k]->SideRate(s, f);
		r[i] = rr;
	}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}

//------------------------------------------------------------------
//  The volume flux out through the open edges (see the declaration). The
//  terms are summed per face in the order FaceVolumeFluxes sums them, and
//  the faces in order over the cells, so the value is that of the face
//  fluxes of the step to the bit and does not depend on the threads.
//------------------------------------------------------------------
template<class T>
double WellMixed_HeightEquationBase<T>::EdgeFluxRate(const Field &f) const
{
	if(!_openEdges) return 0;
	std::vector<Site> svec(T::nStencil);
	double q = 0;
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		Site s = G->template Of<Site>(ii);
		const int i = Index(s);
		if( _crsA.IsBC(i) ) continue;
		bool any = false;
		for(int d = 0; d < T::nDir && !any; ++d)
			any = _crsA.IsBC(Index(T::Adj(G, s, T::Dir(d))));
		if(!any) continue;
		G->SetSite(svec, s);
		for(int d = 0; d < T::nDir; ++d)
		{
			if( ! _crsA.IsBC(Index(T::Adj(G, s, T::Dir(d)))) ) continue;
			if( ! BoundaryCondition::IsOpen(_boundary->TypeOfSide(d)) ) continue;
			double F = 0;
			for(size_t k = 0; k < _physics.size(); ++k)
				if(_physics[k]->Enabled()) F += _physics[k]->FaceFlux(s, d, f, svec);
			q += (d % 2 == 1) ? F : -F;
		}
	}
	return q;
}

//------------------------------------------------------------------
//  The evaporation and the supply of the step, row by row, as SetMatrixA
//  puts them into B (evp - src) for the film f
//------------------------------------------------------------------
template<class T>
void WellMixed_HeightEquationBase<T>::StepSources(const Field &f)
{
	const double dt = _timeStep->DT().now;
	const double tmid = _timeStep->NowTime() + 0.5 * dt;
	if(_evapOld.size() != _nRows) { _evapOld.resize(_nRows); _evapOld = 0; }
	if(_evpStep.size() != _nRows) _evpStep.resize(_nRows);
	if(_srcStep.size() != _nRows) _srcStep.resize(_nRows);
	_evpStep = 0;
	_srcStep = 0;
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		Site s = G->template Of<Site>(ii);
		const int i = Index(s);
		if( _crsA.IsBC(i) ) continue;
		_evpStep[i] = (1.0 - _theta) * _evapOld[i] + _theta * EvapAt(f, s, dt);
		double src = _Source(s);
		if(_mms.Enabled())
		{
			double xx, yy;
			T::SiteXY(G, s, xx, yy);
			src += _mms.Snd(xx, yy, tmid);
		}
		_srcStep[i] = src;
	}
}

//------------------------------------------------------------------
//  The fluxes of the step the film took. The Newton iteration ends on a
//  film h1 whose step is h0 - h1 = dt [theta div F + (1 - theta) div F0 +
//  sources] only to the residual of the step there: the fluxes F at h1 are
//  not quite the ones the iteration solved with (it linearized them about
//  the film before its last correction), and on a film over cells far
//  smaller than it is thick the residual is at the rounding of fluxes of
//  enormous size (eps dt / dt_explicit: tools/analytic.py nozzle_balance,
//  the film 1.2 mm over 3.3 um cells), so the film of the fluxes, h0 - dt
//  (...), is not one to keep (it went negative there). The film stays as
//  the iteration left it and the fluxes take up the residual instead: rho,
//  per real cell, is what the divergence of the fluxes misses of the step
//  (its mean over the cells, the volume error of the step, at rounding, is
//  left alone), and phi, added to F1 / theta, is a flux through the inner
//  faces whose divergence is rho: along y the sum of each row of cells,
//  spread evenly over the row, then along x within each row what is left
//  of it, both as running sums (in 1-D the running sum of rho alone). It
//  is exact to the rounding of the sums, needs no solve, and touches no
//  edge face: the fluxes through the edges, and so the balance of what
//  crosses them, are those of the film, and nothing goes through a wall.
//  phi is of the size of the residual.
//------------------------------------------------------------------
template<class T>
void WellMixed_HeightEquationBase<T>::ConsistentFluxes(double dt)
{
	std::vector<double> rho((size_t)_nRows, 0.0);
	double sum = 0;
	int n = 0;
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		Site s = G->template Of<Site>(ii);
		const int i = Index(s);
		if( _crsA.IsBC(i) ) continue;
		double div1 = 0, div0 = 0;
		for(int a = 0; a < T::nAxis; ++a)
		{
			const size_t il = (size_t)Index(T::Adj(G, s, T::Dir(2*a)));
			div1 += _F1[a][(size_t)i] - _F1[a][il];
			div0 += _F0[a][(size_t)i] - _F0[a][il];
		}
		const double r = (_phi_old[i] - _phi_guess[i]) / dt
		               - ( _theta * (div1 + _r1[i]) + (1.0 - _theta) * (div0 + _r0[i])
		                 + (_evpStep[i] - _srcStep[i]) );
		rho[(size_t)i] = r;
		sum += r;
		++n;
	}
	if(n == 0) return;
	const double mean = sum / n;

	int nx, ny;
	T::SetExtent(G, nx, ny);
	const double w = 1.0 / _theta;
	if(T::nAxis == 1)
	{
		//  The real cells are 2 .. nx-3; the face between ix and ix+1 is
		//  F[0][ix], the last one (nx-3) an edge face
		double acc = 0;
		for(int ix = 2; ix < nx - 3; ix++)
		{
			acc += rho[(size_t)ix] - mean;
			_F1[0][(size_t)ix] += w * acc;
		}
		return;
	}
	const int m = nx - 4;		// real cells in a row
	std::vector<double> rowSum((size_t)ny, 0.0);
	for(int jy = 2; jy <= ny - 3; jy++)
	{
		double r = 0;
		for(int ix = 2; ix <= nx - 3; ix++) r += rho[(size_t)(jy*nx + ix)] - mean;
		rowSum[(size_t)jy] = r;
	}
	//  Along y: the face between rows jy and jy+1 carries the sum of the
	//  rows up to jy, a share 1/m of it through every column
	std::vector<double> &Fy = _F1[T::nAxis - 1];
	double accY = 0;
	for(int jy = 2; jy < ny - 3; jy++)
	{
		accY += rowSum[(size_t)jy] / m;
		for(int ix = 2; ix <= nx - 3; ix++) Fy[(size_t)(jy*nx + ix)] += w * accY;
	}
	//  Along x, within each row: what the y faces did not take
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(ny))
#endif
	for(int jy = 2; jy <= ny - 3; jy++)
	{
		const double share = rowSum[(size_t)jy] / m;
		double acc = 0;
		for(int ix = 2; ix < nx - 3; ix++)
		{
			acc += (rho[(size_t)(jy*nx + ix)] - mean) - share;
			_F1[0][(size_t)(jy*nx + ix)] += w * acc;
		}
	}
}


template<class T>
void WellMixed_HeightEquationBase<T>::SetMatrixJacobi()
{
	int jj;

	if(!_bcAppliedToJacobi)
	{
		_crsJacobi = _crsA;		//Copy the structure of crsA

		//Set the boundary-condition part of the matrix (once, at initialization)
		_boundary->Apply(_crsJacobi, _Rhs, BoundaryCondition::HeightJacobi, *_field);

		_bcAppliedToJacobi = true;
	}

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(jj) if(yuzen::opt::Par(_crsA.RowSize()))
#endif
	for(int irow = 0; irow <_crsA.RowSize(); irow++)		// Note that G.imin() is a management index!!
	{
		if( ! _crsA.IsBC(irow) )
		{
			for(jj = 0; jj < _crsA.NonZero(irow); jj++)
				_crsJacobi.Element(irow, jj) = _crsA.Element(irow, jj);
		}
	}

	// Add the contribution of each physical model that is in use
	for(size_t k = 0; k < _physics.size(); ++k)
		if(_physics[k]->Enabled()) _physics[k]->SetJacobiMatrix(_crsJacobi, *_field_guess);
}


//------------------------------------------------------------------
//  Cap the concentration of the working field at 1 (for the coefficients)
//------------------------------------------------------------------
template<class T>
void WellMixed_HeightEquationBase<T>::ClampConc()
{
	Site s;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(s) if(yuzen::opt::Par((G->imax() - G->imin() + 1)))
#endif
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		s = G->template Of<Site>(ii);
		if(_field_guess->CC(s) > 1.0) _field_guess->CC(s) = 1.0;
	}
}


//------------------------------------------------------------------
//  0 <= c <= 1 without taking solute away or making any (see the
//  declaration). The film is the one the step ended on (_field_guess, its
//  h already updated), the one the solute balance (FieldInfo) sums h c
//  over.
//
//  Each pass is independent of the order of the cells: the room of every
//  cell is taken before the pass moves anything, the cells outside the
//  bound give (or take) and are not touched otherwise, and what a cell
//  receives (or gives) is applied after every cell has been seen. A pass
//  that leaves a cell outside (its neighbors full, or empty, or asked by
//  two at once beyond what they have) is followed by another, from where
//  the solute then is, up to four; what is left after them is spread over
//  the film.
//------------------------------------------------------------------
template<class T>
void WellMixed_HeightEquationBase<T>::RedistributeExcess(Vector<double> &cx)
{
	_excessMoved = _excessTrimmed = 0;
	_deficitMoved = _deficitAdded = 0;
	bool over = false, under = false;
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		Site s = G->template Of<Site>(ii);
		const int i = Index(s);
		if(_crsA.IsBC(i)) continue;
		if(cx[i] > 1.0) over = true;
		if(cx[i] < 0.0) under = true;
	}
	//  Above 1 first: what it moves goes only to cells below 1, and filling
	//  a cell below 0 then takes only from cells above 0, so neither undoes
	//  the other
	if(over)  MoveExcess(cx, false);
	if(under) MoveExcess(cx, true);
	if(over || under)
	{
		//  The ghosts from the real cells as they now are, by the boundary
		//  rows of the solve. Left as solved, a Periodic or Symmetry ghost
		//  would hold the value from before the move, and the explicit half
		//  of the next step would carry a flux through that edge that its
		//  partner does not (1e-8 of the solute in a front crossing a
		//  periodic seam)
		_concEq->GetMatrixIA().SatisfyBCRows(_concEq->GetRHS(), cx);
	}
	//  Ghost rows still outside [0, 1] (an Outlet extrapolating a real cell
	//  at a bound) are set to the bound: they hold no solute of the film,
	//  and the boundary condition of the next solve sets them again
	for(int i = 0; i < cx.size(); i++)
	{
		if(cx[i] > 1.0) cx[i] = 1.0;
		if(under && cx[i] < 0.0) cx[i] = 0.0;
	}
}

//  The passes. Above 1 (lower false) a cell holds e = h (c - 1) too much
//  and a cell below 1 has room h (1 - c); below 0 (lower true) a cell
//  lacks e = -h c and a cell above 0 has h c to give. In either case the
//  cell outside goes to the bound and its neighbors give up (or take) the
//  fraction e / (their room) of their room.
template<class T>
void WellMixed_HeightEquationBase<T>::MoveExcess(Vector<double> &cx, bool lower)
{
	bool any = true;
	const double bound = lower ? 0.0 : 1.0;
	auto outside = [lower](double c) { return lower ? (c < 0.0) : (c > 1.0); };
	//  What a cell holds beyond the bound, and the room a cell has, as
	//  solute (both >= 0)
	auto beyond = [lower](double h, double c) { return lower ? -h * c : h * (c - 1.0); };
	auto roomOf = [lower](double h, double c) { return lower ? h * c : h * (1.0 - c); };
	double &moved = lower ? _deficitMoved : _excessMoved;
	double &cut   = lower ? _deficitAdded : _excessTrimmed;

	std::vector<double> room((size_t)cx.size(), 0.0), gain((size_t)cx.size(), 0.0);
	const int nPass = 4;
	for(int pass = 0; pass < nPass && any; ++pass)
	{
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			const int i = Index(s);
			const double h = _field_guess->HH(s);
			const bool inside = lower ? (cx[i] > 0.0) : (cx[i] < 1.0);
			room[i] = (!_crsA.IsBC(i) && h > 0.0 && inside) ? roomOf(h, cx[i]) : 0.0;
			gain[i] = 0.0;
		}
		any = false;
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			const int i = Index(s);
			if(_crsA.IsBC(i) || !outside(cx[i])) continue;
			const double h = _field_guess->HH(s);
			//  A cell without film holds no solute to move: the bound there
			//  changes the sum by what is counted
			if(!(h > 0.0)) { cut += beyond(h, cx[i]); cx[i] = bound; continue; }
			const double e = beyond(h, cx[i]);
			if(pass == 0) moved += e;
			double r = 0.0;
			for(int d = 0; d < T::nDir; ++d) r += room[Index(T::Adj(G, s, T::Dir(d)))];
			if(!(r > 0.0)) { any = true; continue; }
			const double f = (e < r) ? e / r : 1.0;
			for(int d = 0; d < T::nDir; ++d)
			{
				const int j = Index(T::Adj(G, s, T::Dir(d)));
				gain[j] += f * room[j];
			}
			//  What the neighbors could not take (or give) stays, for the
			//  next pass
			if(f < 1.0) cx[i] = bound;
			else        { cx[i] = lower ? -(e - r) / h : 1.0 + (e - r) / h; any = true; }
		}
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			const int i = Index(s);
			if(gain[i] > 0.0)
			{
				const double dc = gain[i] / _field_guess->HH(s);
				cx[i] += lower ? -dc : dc;
				if(outside(cx[i])) any = true;
			}
		}
	}
	if(!any) return;

	//  Left after the passes: spread over every real cell with room, in
	//  proportion to it (c + f (1 - c) stays at or below 1, and c (1 - f) at
	//  or above 0, for f <= 1)
	double left = 0.0, total = 0.0;
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		Site s = G->template Of<Site>(ii);
		const int i = Index(s);
		if(_crsA.IsBC(i)) continue;
		const double h = _field_guess->HH(s);
		if(outside(cx[i])) { left += beyond(h, cx[i]); cx[i] = bound; }
		else if(h > 0.0)   total += roomOf(h, cx[i]);
	}
	const double f = (total > 0.0) ? std::min(1.0, left / total) : 0.0;
	if(f > 0.0)
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			const int i = Index(s);
			if(_crsA.IsBC(i) || !(_field_guess->HH(s) > 0.0)) continue;
			if(lower)
			{
				if(!(cx[i] > 0.0)) continue;
				cx[i] -= f * cx[i];
				if(cx[i] < 0.0) cx[i] = 0.0;
			}
			else
			{
				if(!(cx[i] < 1.0)) continue;
				cx[i] += f * (1.0 - cx[i]);
				if(cx[i] > 1.0) cx[i] = 1.0;
			}
		}
	//  More solute than liquid in all (above 1), or less than none (below
	//  0): the rest is cut, or made up
	if(total < left) cut += left - total;
}


template<class T>
void WellMixed_HeightEquationBase<T>::SetIterate(const Field &f)
{
	if(_cIterate.size() != _nRows) _cIterate.resize(_nRows);
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		Site s = G->template Of<Site>(ii);
		_cIterate[Index(s)] = f.CC(s);
	}
}


//------------------------------------------------------------------
//  Measure how far the concentration moved during the outer iteration
//------------------------------------------------------------------
template<class T>
double WellMixed_HeightEquationBase<T>::ConcCorrection()
{
	double dmax = 0;
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		Site s = G->template Of<Site>(ii);
		const int i = Index(s);
		const double cc = _field_guess->CC(s);
		const double d  = std::fabs(cc - _concPrev[i]);
		if(d > dmax) dmax = d;
		_concPrev[i] = cc;
	}
	return dmax;
}


//------------------------------------------------------------------
//  Called at the start of every attempt at a time step: the first, and each
//  retry after the step was rejected (a Newton iteration that failed, or a
//  converged one whose truncation criterion was too large).
//------------------------------------------------------------------
template<class T>
void WellMixed_HeightEquationBase<T>::BeforeNKLoop()
{
	//  Start the working field from the state the step starts from. NKLoop
	//  builds A, and with it the explicit half A(h^n) h^n, from the working
	//  field before the first iteration writes into it. After a rejected
	//  attempt that field holds the attempt's last iterate: a converged but
	//  rejected answer, or a diverged one with values of 1e76 and NaN in it.
	//  An explicit half built from it would make the retry either take that
	//  wrong half into its answer or fail at once, and the step be cut again
	//  for nothing. After an accepted step the two fields are
	//  already equal (AfterNKLoop copies one into the other), so this
	//  changes nothing there.
	*_field_guess = *_field;
	//  Nothing moved above c = 1 yet in this step (RedistributeExcess)
	_excessMoved = _excessTrimmed = 0;

	//  The evaporation at the start of the step, for the explicit half of B
	const double dt = _timeStep->DT().now;
	if(_evapOld.size() != _nRows) _evapOld.resize(_nRows);
	_evapOld = 0;
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		Site s = G->template Of<Site>(ii);
		if( ! _crsA.IsBC(Index(s)) ) _evapOld[Index(s)] = EvapAt(*_field, s, dt);
	}
	_Source.SetEvalTime(_timeStep->NowTime() + 0.5 * dt);
	if(_Source.Enabled())
	{
		if(_cStepStart.size() != _nRows) _cStepStart.resize(_nRows);
		_cStepStart = 0;
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			_cStepStart[Index(s)] = _field->CC(s);
		}
	}
	//  The sideways outflow at the start of the step, as A(h^n) h^n has it
	_sideOld = _Centrifugal.SideOutflowRate(_crsA, *_field, &_sideOldResin);
	//  And what flows out through the open edges there
	_edgeOld = EdgeFluxRate(*_field);

	if(!_swConcField) return;

	//  The limiter of the first pass from the state the step starts from
	SetIterate(*_field);

	//  Set the reference of the iteration to the concentration at the start
	//  of this step
	if(_concPrev.size() != _nRows) _concPrev.resize(_nRows);
	_minStart = 0;
	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		Site s = G->template Of<Site>(ii);
		const double c = _field->CC(s);
		_concPrev[Index(s)] = c;
		if( ! _crsA.IsBC(Index(s)) && (c < _minStart || !(c == c)) ) _minStart = c;
	}

	//  The film side of the concentration equation: the face fluxes and the
	//  sideways outflow at the start of the step (the explicit half), the
	//  film there, and the solves of the step from d = 0
	FaceVolumeFluxes(*_field, _F0, _r0);
	_cstep.h0  = &_phi_old;
	_cstep.F0  = _F0;
	_cstep.F1  = _F1;
	_cstep.r0  = &_r0;
	_cstep.r1  = &_r1;
	_cstep.src = &_srcStep;
	_cstep.evp = &_evpStep;
	//  The concentration the supplied liquid arrives at (ConcStep.h)
	_cstep.srcSupplied = _Source.Enabled() && _Source.ConcSupplied();
	_cstep.srcConc     = _Source.Conc();
	_concEq->ResetGuess();
}


template<class T>
bool WellMixed_HeightEquationBase<T>::EveryNKLoop()
{
	//Called once per time step
	_field_guess->Update_H(_phi_guess);

	_Source.updateFraction();

	bool ch;
	int iter; double resid;
	if(_swConcField)
	{
		//  The concentration with the fluxes of this pass's film, and that
		//  film as the film of the end (its fluxes give it exactly only
		//  once the iteration has converged; the step is completed with the
		//  fluxes made those of the step the film took, AfterNKLoop). The
		//  velocities are not needed for it: they are formed for the state
		//  the step ends on (AfterNKLoop).
		//  A pass whose film has a cell without liquid leaves the
		//  concentration as it is (the balance of h c is not one there); if
		//  the film the iteration ends on is so too, the divergence test
		//  stops the run on the film
		if(FilmNotPositive(_phi_guess))
		{
			_coupledCorrection = 0;
			return true;
		}
		FaceVolumeFluxes(*_field_guess, _F1, _r1);
		StepSources(*_field_guess);
		_cstep.h1 = &_phi_guess;
		_cstep.commit = false;
		ch = _concEq->TimeStepTrial(iter, resid);
		//  A concentration solve that did not converge fails this Newton
		//  pass and with it the step: cut and retried under step control,
		//  the run stopped at a fixed step (C_NewtonKantorovich::TimeAdvance),
		//  never taken as it stands.
		if(!ch)
		{
			_failCause = "the concentration solve did not converge (Numerics.ConcEquation.LinearSolver)";
			return false;
		}

		_field_guess->Update_C(_concEq->GetX());
		//  The limiter of the next pass from the iterate as solved
		_cIterate = _concEq->GetX();

		ClampConc();

		//  Hand how far the concentration moved during this iteration to
		//  the convergence test
		_coupledCorrection = ConcCorrection();
	}
	return true;
}


template<class T>
bool WellMixed_HeightEquationBase<T>::AfterNKLoop()
{
	_field_guess->Update_H(_phi_guess);

	if(!_swConcField)
	{
		CalVelocity();
		//  The sideways outflow of the step, as the solve took it: the
		//  implicit half at the film the step ended on
		double rr = 0;
		const double r = _Centrifugal.SideOutflowRate(_crsA, *_field_guess, &rr);
		_sideStep      = _theta * r  + (1.0 - _theta) * _sideOld;
		_sideStepResin = _theta * rr + (1.0 - _theta) * _sideOldResin;
		//  The same of the open edges (no solute without the concentration)
		_edgeStep      = _theta * EdgeFluxRate(*_field_guess) + (1.0 - _theta) * _edgeOld;
		_edgeStepResin = 0;
		_field->operator=(*_field_guess);
		return true;
	}

	//--------------------------------------------------------------
	//  The step with the concentration. The Newton passes each solved the
	//  concentration with the film of that pass, so the iteration has
	//  converged on the two together. What the step keeps is one
	//  finite-volume step of both (ConcStep.h): the film the iteration
	//  ended on; the face volume fluxes and the sideways outflow at that
	//  film, with the concentration of its last pass, made the fluxes of
	//  the step the film took (ConsistentFluxes: they miss it by the
	//  residual of the step there); the evaporation and the supply as the
	//  last pass took them; and the concentration solved once more with
	//  exactly those. The solute of every cell then changes by the solute
	//  fluxes through its faces and its sources only, and a uniform
	//  concentration stays uniform.
	//
	//  Do not call the concentration's TimeAdvance() here: it would solve
	//  one more step on top of this one.
	//
	//  All of it is over the step that was solved (_stepDt). When the
	//  time-step control has just lengthened dt, DT().now already holds the
	//  next step's (TimeStepAdjustment runs before AfterNKLoop). DT().now is
	//  set back for this and restored after it (ChangeDT keeps old and now,
	//  so the pair leaves them as they were).
	//--------------------------------------------------------------
	const double dtNext = _timeStep->DT().now;
	const bool swapDt = (dtNext != _stepDt) && _timeStep->ChangeDT(_stepDt);

	FaceVolumeFluxes(*_field_guess, _F1, _r1);
	//  What the film-height equation's fluxes take out through the open
	//  edges: those of the edge faces of _F1 and _F0, which ConsistentFluxes
	//  leaves as they are
	_edgeStep = _theta * EdgeFluxRate(*_field_guess) + (1.0 - _theta) * _edgeOld;
	StepSources(*_field_guess);
	//  The evaporation as B took it in the last pass, the one the film's
	//  step holds (and the one the volume balance counts, EvaporationRate)
	if(_evapB.size() == _nRows) _evpStep = _evapB;
	ConsistentFluxes(_stepDt);
	_cstep.h1 = &_phi_guess;
	_cstep.commit = true;

	//  A limited scheme freezes its limiter at the last iterate; the solve
	//  is repeated with it frozen at its own answer until that settles
	//  (a few passes), so that the face values of the step are the limited
	//  face values of its answer
	const Vector<double> lastPass = _concEq->GetX();
	const bool limited = _concEq->ConvectionOn() && _concEq->Scheme().IsLimited();
	//  A film with a cell without liquid: the last pass stands, and the
	//  divergence test of the model looks at the film
	const bool noFilm = FilmNotPositive(_phi_guess);
	bool solved = true;
	for(int k = 0; k < (limited ? 8 : 1) && !noFilm; ++k)
	{
		int it2; double r2;
		if(!_concEq->TimeStepTrial(it2, r2)) { solved = false; break; }
		if(!limited) break;
		const Vector<double> &x = _concEq->GetX();
		double dmax = 0;
		for(int i = 0; i < x.size(); i++) dmax = std::max(dmax, std::fabs(x[i] - _cIterate[i]));
		_cIterate = x;
		if(dmax <= 1e-12) break;
	}
	if(swapDt) _timeStep->ChangeDT(dtNext);
	if(!solved)
	{
		//  The last pass stands
		_concEq->GetX() = lastPass;
		Global::Warning<<"The concentration solve that completes the step did not converge;"
		               <<" the step keeps the last Newton pass, and its solute balance is off by that."<<nl;
	}

	//  The sideways outflow of the step, as both equations took it (the
	//  solute at the concentration the solve gave)
	{
		const Vector<double> &x = _concEq->GetX();
		double r = 0, rr = 0;
		for(int ii = G->imin(); ii <= G->imax(); ii++)
		{
			Site s = G->template Of<Site>(ii);
			const int i = Index(s);
			if( _crsA.IsBC(i) ) continue;
			r  += _r1[i];
			rr += _r1[i] * x[i];
		}
		_sideStep      = _theta * r  + (1.0 - _theta) * _sideOld;
		_sideStepResin = _theta * rr + (1.0 - _theta) * _sideOldResin;
	}
	//  The solute through the open edges, as the solve that completed the
	//  step put it through them (none if the last pass stands)
	_edgeStepResin = (solved && !noFilm) ? _concEq->EdgeOutflowRate() : 0.0;
	_cStepEnd = _concEq->GetX();

	//  0 <= c <= 1 on the equation's own state, before the step is kept, so
	//  that the next step starts from the field that is written. With the
	//  evaporation limited to the solvent there is (EvapAt) c leaves [0, 1]
	//  only by rounding, and by what an unbounded advection scheme (KK,
	//  UTOPIA, Central4, QUICK; a limited one beyond Courant 1) over- or
	//  undershoots by at a sharp front. A cap would cut the solute beyond
	//  the bound (0.7 % of it as a front from c = 0.99 to 0.3 crossed a
	//  periodic domain of 100 cells once with UTOPIA; tools/analytic.py
	//  conc_overshoot); RedistributeExcess moves it to and from the cells
	//  around instead.
	RedistributeExcess(_concEq->GetX());
	_concEq->CommitStep();
	_field_guess->Update_C(_concEq->GetX());
	_cIterate = _concEq->GetX();

	ClampConc();
	//  The velocities of the state kept (the output, the outflow of the
	//  balance)
	CalVelocity();

	_field->operator=(*_field_guess);

	return true;
}


template<class T>
void WellMixed_HeightEquationBase<T>::CalVelocity()
{
	YUZEN_PROF("velocity");
	// Each term writes into a velocity field of its own, so the order does
	// not affect the results
	for(size_t k = 0; k < _physics.size(); ++k)
		if(_physics[k]->Enabled()) _physics[k]->CalVelocity(*_field_guess);

	SetVelocityBC(_field_guess->UU().U_LaplacePressure());
	SetVelocityBC(_field_guess->UU().U_SolutoCapillarity());
	SetVelocityBC(_field_guess->UU().U_VerticalGravity());
	SetVelocityBC(_field_guess->UU().U_HorizontalGravity());
	SetVelocityBC(_field_guess->UU().U_Centrifugal());
	SetVelocityBC(_field_guess->UU().U_DisjoiningPressure());

	_field_guess->UU().Cal_MixedVelocity();
	SetVelocityBC(_field_guess->UU().UU());
}

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_HEIGHT_EQUATION_BASE_INCLUDED_
