/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CZ_HeightEquation_2D.h : film-height equation of the ConcZ model (semi-implicit)
*
*  CZ_HeightEquation.h
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _CZ_HEIGHT_EQUATION_2D_INCLUDED_
#define _CZ_HEIGHT_EQUATION_2D_INCLUDED_


#include "ManufacturedSource.h"
#include "LubFields.h"
#include "NonDimensionalNumbers.h"
#include "InitialCondition.h"
#include "BoundaryCondition.h"
#include "FiniteDifferences.h"
#include "CZ_ConcEquation_2D.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

/////////////////////////////////////////
/**
 * @brief Film-height equation of the ConcZ model.
 *
 * The film is advanced in conservative form, \f$h_t = -\partial_x Q - E\f$,
 * with the flux at a face built from the positivity preserving mobility of
 * Diez and Kondic (2002), Eq.~(23). Writing it this way keeps the volume
 * exactly constant, lets the mobility switch the flux off at the edge of a
 * dry region, and leaves the capillary term implicit, so the time step
 * carries no explicit stability limit from it.
 *
 * Because the viscosity varies through the thickness with the local
 * concentration, the in-plane velocity is not analytic; CalVelocity()
 * integrates it through the film, and the concentration equation advects
 * with it. The flux of the film-height equation is the depth integral of
 * that profile: the surface tension is that of the surface concentration,
 * the Marangoni stress is driven by its gradient, and the viscosity enters
 * as the depth means of 1/mu that the integral weights it with.
 *
 * The two equations are coupled step by step. After each accepted trial the
 * concentration equation is advanced by the same step, and the film's rate of
 * change is passed to it so that its stretching zeta grid follows the film.
 */
class CZ_HeightEquation_2D : public C_SemiImplicitEquation
{
public:
	CZ_HeightEquation_2D() : _dtPrev(1.0), _ccPrevSet(false), _stepDt(0), _couplingPasses(0), _couplingCuts(0), _couplingFloor(1e-10), _warnedExcess(false), _deficitMoved(0), _deficitAdded(0), _minStart(0), _warnedCoupling(false), _solvedDt(0.0) {}

	//------------------------------------------------------------------
	//  Orientation of a face. A face always lies between cell (ix,jy) and its
	//  neighbor on the negative side: between it and (ix-1,jy) for FaceX,
	//  between it and (ix,jy-1) for FaceY.
	//  When Ny = 1 no FaceY face ever appears, so the y-direction contribution
	//  simply drops out.
	//------------------------------------------------------------------
	enum tFaceDir { FaceX = 0, FaceY = 1 };

	/**
	 * @brief Assemble the equation and take hold of the concentration equation.
	 * @param inp     the whole input file; the blocks this equation needs are read from it
	 * @param prop    fluid properties (viscosity and the concentration dependence)
	 * @param ts      time-step object supplying dt and the current time
	 * @param concEq  concentration equation advanced together with the film
	 * @param ff      ConcZ field: h and s on the 2-D grid, c and u on the 3-D grid
	 * @param nonDim  non-dimensional numbers of the problem
	 * @param bc      boundary condition applied to the rows of this equation
	 */
	void Build(Input &inp, C_Fluid *prop, C_TimeSteps *ts, CZ_ConcEquation_2D *concEq,
		Field_CZ *ff, NonDimensionalNumbers *nonDim, BoundaryCondition2D *bc);

	/** @brief dh/dt of the step just computed, as handed to the concentration equation. */
	const Vector<double>& GetDtDh() const {return _dtdh;}
	/** @brief Form dh/dt for this step as the difference of the film height over dt. */
	void CalDtDh();

	/**
	 * @brief Time-step limit from the local truncation error of the film height.
	 * @return the criterion the time-step controller compares against
	 */
	virtual double CheckTimeStepCriterion();

	/**
	 * @brief Compute the velocity field through the thickness of the film.
	 *
	 * The in-plane velocity is integrated in z by Simpson's rule, with the
	 * viscosity taken from the local concentration, and receives whichever of
	 * the capillary, Marangoni, gravity and centrifugal terms are enabled.
	 * The vertical velocity then follows from continuity, and the ghost layers
	 * are filled according to the boundary type.
	 */
	void CalVelocity();

	/**
	 * @brief Mobility at the face on the negative side of cell (ix, jy), in the
	 *        positivity preserving form of Diez & Kondic (2002), Eq.~(23).
	 *
	 * @param ix   x index of the cell the face belongs to
	 * @param jy   y index of that cell (only 0 is used when Ny = 1)
	 * @param dir  which face: FaceX is the one shared with (ix-1, jy), FaceY
	 *             the one shared with (ix, jy-1)
	 * @return    2 h1^2 h2^2 / (h1 + h2) / 3 times the mean of 1/mu at the face
	 *
	 * The film-height equation is written in conservative form, h_t = -dQ/dx,
	 * and this is the coefficient that multiplies the driving force to give
	 * the flux Q at a face. The form matters: it vanishes as O(h^2) when
	 * either side thins, so the flux switches itself off at the edge of a dry
	 * region and the scheme cannot drive h negative. The manual's
	 * discretization chapter derives it; tests/test_positivity_mobility.cpp
	 * pins the property down.
	 *
	 * The viscosity enters as the mean of 1/mu through the film weighted by
	 * 3 (1 - zeta)^2, which makes the flux h^3 b int (1 - zeta)^2/mu dzeta the
	 * depth integral of the velocity profile of CalVelocity().
	 */
	double FaceMobility(int ix, int jy, tFaceDir dir) const;
	/** @brief FaceMobility() as the film has it, a wall or not (the sideways outflow of 1-D). */
	double OpenFaceMobility(int ix, int jy, tFaceDir dir) const;
	/**
	 * @brief Whether the face is that of a wall (a Symmetry edge): no flux
	 *        and no velocity of any term crosses it.
	 */
	bool WallFace(int ix, int jy, tFaceDir dir) const;
	/** @brief The mean c of the real cells of column (ix, jy): its solute over h. */
	double ColumnMeanConc(Field_CZ &f, int ix, int jy) const;
	/**
	 * @brief Keep c at or below 1 without taking solute away: in each
	 *        column, what a cell holds above c = 1 is passed on to the cell
	 *        below it, and what reaches the substrate back up to the cells
	 *        that have room. Returns how much solute (in c times cells) was
	 *        moved, and the largest c left above 1 (a column whose mean is
	 *        above 1, which the evaporation limit does not let happen).
	 */
	double RedistributeExcess(Vector<double> &c, double &cmaxLeft) const;
	/**
	 * @brief Keep c at or above 0 without making solute (an unbounded
	 *        advection scheme undershoots at a sharp front). A column whose
	 *        cells hold solute in all fills its cells below 0 from them, in
	 *        proportion to what each holds; a column with less than none
	 *        takes what it lacks from the columns next to it that hold
	 *        solute, in proportion to it (their profiles scaled), in up to
	 *        four passes, and then from every such column. What a film with
	 *        less than no solute in all still lacks is made up.
	 * @param c      [in,out] the concentration of every layer, made >= 0
	 * @param moved [out] the solute moved into cells below 0, -c h/n
	 *               summed (n the layers), per unit cell measure
	 * @param added  [out] the solute made up
	 */
	void RedistributeDeficit(Vector<double> &c, double &moved, double &added) const;

	/**
	 * @brief Mean of 1/mu through the film at a face, weighted by
	 *        (1 - zeta)^power (2: the mobility, 1: the Marangoni flux).
	 */
	double FaceInvViscosity(int ix, int jy, tFaceDir dir, int power) const;
	/** @brief Surface concentration at the face, for the surface tension. */
	double FaceSurfConc(int ix, int jy, tFaceDir dir) const;
	/**
	 * @brief Marangoni flux at a face.
	 *
	 * Driven by the gradient of the surface concentration: h^2/2 times the
	 * mean of 1/mu through the film weighted by 2 (1 - zeta), the depth
	 * integral of the shear flow the surface stress drives.
	 */
	double MarangoniFlux(int ix, int jy, tFaceDir dir) const;

	/**
	 * @brief Passes of the film and the concentration taken again, in the
	 *        last step, with the concentration of the midpoint the pass
	 *        before found (0: the first pass, with the extrapolated one,
	 *        was kept). Only with the Marangoni term on.
	 */
	int CouplingPasses() const {return _couplingPasses;}
	/** @brief Steps that failed because the coupling passes did not settle. */
	long long CouplingCuts() const {return _couplingCuts;}
	/**
	 * @brief The solute the last step moved into cells below c = 0, and the
	 *        solute it had to make up there (RedistributeDeficit), times
	 *        the cell measure.
	 */
	double DeficitMoved(double cellMeasure) const {return _deficitMoved * cellMeasure;}
	double DeficitAdded(double cellMeasure) const {return _deficitAdded * cellMeasure;}
	/**
	 * @brief The lowest concentration of a real cell at the start of the
	 *        last step (0 if none was below). The divergence test looks at
	 *        it as well as at the end of the step: a field file read with a
	 *        cell far below 0 must stop the run at its first step, which
	 *        RedistributeDeficit may have filled.
	 */
	double MinStartConc() const {return _minStart;}

	/** @brief Mark the ghost rows so they are treated as boundary rows. */
	virtual void SetFlagBC();

	/** @brief Before each trial step: restore the field saved at the end of the last one. */
	virtual bool EveryTimeStepTrial();
	/**
	 * @brief After an accepted trial: commit the film and advance the concentration.
	 *
	 * The new film height is written to the field, dh/dt is passed to the
	 * concentration equation, which is then advanced by the same step with
	 * the velocity of the film halfway through it. With the Marangoni term
	 * on, the film and the concentration are then solved again, one stage
	 * each, with the concentration of the midpoint the pass before found,
	 * until that midpoint settles (the coupling through the surface
	 * concentration is then implicit). The concentration is kept at or
	 * below 1 without taking solute away (RedistributeExcess), and the
	 * field is saved for the next trial.
	 * @return false if a linear solve of the film or of the concentration
	 *         did not converge: the step has failed, and the concentration
	 *         equation and the extrapolation are left as they were before
	 *         it, so that it can be solved again at a shorter step
	 */
	virtual bool AfterTimeStepTrial();

	/**
	 * @brief One time step (C_SemiImplicitEquation::TimeAdvance); with the
	 *        time-step control off, a step that fails stops the run with
	 *        its cause (FailedFixedStep).
	 */
	virtual int TimeAdvance(int &iter, double &resid);

protected:
	//  One step in two stages (see the .cpp): an Euler-like predictor with the
	//  coefficients of the start of the step, then the same system with every
	//  coefficient and every explicit term taken at the midpoint.
	virtual int Solve(int &iter, double &resid);
	int SolveStage(int &iter, double &resid);
	//  c at the midpoint of the step, extrapolated from the last two starts
	bool MidpointConc(Vector<double> &cm) const;

public:
	/**
	 * @brief Take the field as it now stands as the start of the next step.
	 *
	 * Each trial starts from a copy of the field saved at the end of the
	 * step before (EveryTimeStepTrial). A restart writes the field it reads
	 * into the field, not into that copy, so this has to follow it: without
	 * it the first step after a restart began from the initial condition
	 * the run was built with.
	 */
	void TakeFieldAsStart() { _field_tmp = *_field; }

	//  The concentration at the start of the previous step and that step's
	//  dt, which the corrector extrapolates from; a restart file carries them.
	bool HasPreviousConc() const { return _ccPrevSet; }
	const Vector<double>& PreviousConc() const { return _ccPrev; }
	double PreviousDt() const { return _dtPrev; }
	void SetPreviousConc(const Vector<double> &c, double dt) { _ccPrev = c; _dtPrev = dt; _ccPrevSet = true; }

	/** @brief The dt the end of the last step was formed with (ConservativeUpdate). */
	double StepDt() const { return _stepDt; }

	/**
	 * @brief What left the columns in the step last solved other than
	 *        through a face, as ConservativeUpdate took it out: the volume
	 *        that evaporated (evaporated) and the liquid the 1-D centrifugal
	 *        term flung sideways, into the y that is not solved for (the
	 *        value returned), each times the cell measure and StepDt. With
	 *        sideResin, the solute that liquid carried: per column, the
	 *        concentration of its layers weighted as the concentration
	 *        equation spreads the outflow through them, the mean of that at
	 *        the start of the step (cStart) and at its end (the field f).
	 */
	double ColumnSinks(double cellMeasure, const Field_CZ &f, const ScalarField3D<> &cStart,
	                   double &evaporated, double *sideResin) const;
	/** @brief Whether the film is flung sideways (the 1-D centrifugal term). */
	bool HasSideOutflow() const { return _sw_Centrifugal && _Ny == 1; }

	// Grid index to matrix row: (i, j) on the 2-D grid, (i, j, k) on the 3-D one
	inline int IL(const int i, const int j) const {return j*_Nx + i;}
	inline int IL(const int i, const int j, const int k) const {return k*(_Nx*_Ny) +j*_Nx + i;}

private:
	//  The source of the method of manufactured solutions (verification only)
	C_ManufacturedSource _mms;
	//  The concentration at the start of the previous step and that step's
	//  dt, to extrapolate the concentration to the midpoint of this one
	Vector<double> _ccPrev;
	double _dtPrev;
	bool   _ccPrevSet;
	//  The dt ConservativeUpdate formed the end of the last step with
	double _stepDt;
	//  Passes of the coupling taken again in the last step
	int    _couplingPasses;
	//  Steps whose passes did not settle (each one failed: cut, or the run
	//  stopped at a fixed step), and the midpoint change below which the
	//  passes count as settled in any case
	long long _couplingCuts;
	double _couplingFloor;
	//  Whether the warning of RedistributeExcess has been given
	bool   _warnedExcess;
	//  What the last step moved into cells below 0, and made up
	double _deficitMoved, _deficitAdded;
	double _minStart;			// MinStartConc
	//  Whether a step cut for passes that did not settle has been reported
	bool   _warnedCoupling;
	//  The dt the film-height equation was last assembled and solved with
	//  (SetMatrixA). DT().now may already hold the next step's when
	//  AfterTimeStepTrial runs: TimeStepAdjustment lengthens it as soon as
	//  the step is accepted.
	double _solvedDt;
	//  At most this many of them, until the midpoint concentration moves by
	//  no more than kCouplingTol of what the extrapolation missed it by
	//  (Aitken-relaxed; see AfterTimeStepTrial)
	static const int    kMaxCouplingPasses = 8;
	static constexpr double kCouplingTol = 0.1;
	//  The velocity profile at the x faces, for the sideways outflow of 1-D
	std::vector<double> _ppSide;

	//  The flux of each face of the stage solved last (StoreFaceFluxParts):
	//  the coefficient of the capillary term, applied to h, and the rest.
	//  x face ix lies between ix-1 and ix, y face jy between jy-1 and jy;
	//  both indexed IL(ix, jy).
	std::vector<double> _qCapX, _qExpX, _qCapY, _qExpY;
	//  The flux with the h of the solve (ConservativeUpdate), which the
	//  concentration equation's layer fluxes add up to
	std::vector<double> _qX, _qY;
	//  The part of h_t of each column that is not a flux divergence
	//  (evaporation, the manufactured source, the sideways outflow of 1-D),
	//  and that sideways outflow on its own
	std::vector<double> _colSrc, _sideSink;
	//  ... and the evaporation of each column on its own (the balance of
	//  History.dat, ColumnSinks)
	std::vector<double> _evapSink;
	void StoreFaceFluxParts();
	//  Write the end of the step as h_old + dt (-div q + _colSrc), with the
	//  flux q of the Crank-Nicolson h of the solve
	void ConservativeUpdate();
	void CopyGhosts(Vector<double> &h) const;

	// Pointer to grid
	C_Grid2D* _G2;
	C_Grid3D* _G3;
	int _Nx, _Ny, _Nz;

	// Pointers to field data
	Field_CZ *_field;
	Field_CZ _field_tmp;

	// Pointer to fluid property
	C_Fluid *_prop;

	// Nondimensional numbers
	NonDimensionalNumbers *_NonDim;

	C_Evaporation _evap;

	// ConcEquation
	bool swConcField;
	CZ_ConcEquation_2D *_concEq;


	// Boundary Condition
	BoundaryCondition2D *_boundary;

	// Physics
	bool _sw_LaplacePressure;
	bool _sw_SolutoCapillarity;
	bool _sw_Centrifugal;
	bool _sw_VerticalGravity;
	bool _sw_HorizontalGravity;

	//  Center of rotation of the centrifugal force (non-dimensional). As in the
	//  CentrifugalKernel on the well-mixed side, the driving force is the face
	//  coordinate minus this.
	double _rotCenter[2];

	C_FiniteDifferences _fdm;

	Vector<double> _dtdh;


	//  Assemble the kinematic condition of the free surface into the coefficient
	//  matrix: a fourth-order central difference of the surface velocity, and on
	//  the right-hand side (vertical velocity at the surface - evaporation -
	//  contribution of the substrate relief).
	virtual void SetMatrixA();
	//  Overwrite the rows at both ends in x directly (periodic or symmetric).
	void SetBoundaryCondition();

private:
	CZ_HeightEquation_2D(const CZ_HeightEquation_2D&);
	void operator=(const CZ_HeightEquation_2D&);

};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_CZ_HEIGHT_EQUATION_2D_INCLUDED_
