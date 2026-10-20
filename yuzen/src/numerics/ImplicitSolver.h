/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ImplicitSolver.h : base of implicit time integration (theta method, time-step control, steady-state test)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _IMPLICIT_SOLVER_INCLUDED_
#define _IMPLICIT_SOLVER_INCLUDED_

#include <string>

#include "grid/YuzenGrid.h"

#include "MyMath.h"
#include "CRS_Matrix.h"
#include "LinearSolver.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

using namespace yuzen::grid;

class C_TimeSteps;

//--------------------------------------------------------------------------
//  Base of every equation that is advanced in time implicitly.
//
//  The form that is solved is
//
//      d(phi)/dt = -A phi - B
//
//  where A is the part proportional to the unknowns (a matrix) and B is
//  everything else (a vector). Time is discretized with the theta method,
//  where theta = 0 / 0.5 / 1 correspond to forward Euler, Crank-Nicolson and
//  backward Euler respectively. **Only the terms put into A receive the theta
//  weighting**, so a term put into B is effectively backward Euler on its
//  own. Which of the two a physical term goes into therefore bears directly
//  on the accuracy in time.
//
//  Derived classes:
//    C_SemiImplicitEquation  freezes the nonlinear coefficients at the values
//                            of the previous step and solves once
//    C_NewtonKantorovich      resolves the nonlinearity as well, by Newton
//                            iteration
//
//  This class itself holds the matrices and the solution vectors, the
//  time-step control (the test on the local truncation error), and the
//  steady-state test.
//--------------------------------------------------------------------------
class C_ImplicitSolver
{
public:
	//  The weighting in time. theta = 0 / 0.5 / 1.
	typedef enum
	{
		ForwardEuler,
		CrankNicholson,
		BackwardEuler
	}tImplicitType;

	C_ImplicitSolver();
	virtual ~C_ImplicitSolver(){}

	/**
	 * @brief Name of the linear-solver backend this equation ended up using.
	 *
	 * This is the resolved name, not the string from the input file: with
	 * Backend = "auto" it reports whichever implementation was chosen, and
	 * after a fall-back it reports the one that actually ran. It is written
	 * to JobInfo.txt so that a result directory records how it was solved.
	 */
	std::string BackendName() const { return _bcg.BackendName(); }

	/**
	 * @brief Put the matrix and the initial values into a state from which
	 *        the computation can start.
	 *
	 * Call this after Build() and SetX().
	 */
	virtual void Initialize() = 0;

	/**
	 * @brief Advance one time step.
	 * @param iter   [out] number of iterations used along the way
	 * @param resid  [out] residual reached along the way
	 * @return 0 on success; a non-zero value reports a failure to advance.
	 *
	 * Both outputs are the values printed to CalcInfo.dat. When the time
	 * step is adaptive, the trials and the adjustment of the step happen
	 * inside this call.
	 */
	virtual int TimeAdvance(int &iter, double &resid) = 0;

	/**
	 * @brief Allocate the storage and build the matrix structure for a given
	 *        number of unknowns.
	 * @param n               number of unknowns
	 * @param ts              the time-step object shared by the equations
	 * @param t               weighting in time (theta method)
	 * @param asjustTimeStep  set this and the equation becomes the one that
	 *                        decides the time step. Only the film-height
	 *                        equation passes true; the concentration
	 *                        equation follows the height.
	 */
	void Build(const int n, C_TimeSteps &ts, tImplicitType t = CrankNicholson, bool asjustTimeStep = false);

	/** @brief Build() taking the number of unknowns from the length of a vector. */
	void Build(const Vector<double>  &f, C_TimeSteps &ts, tImplicitType t = CrankNicholson, bool asjustTimeStep = false)
	{
		Build(f.size(), ts, t, asjustTimeStep);
	}

	/**
	 * @brief Build() taking the number of unknowns from the grid of a field.
	 *
	 * The ghost points are counted as unknowns as well.
	 */
	template<class GRID, class SPOT>
	void Build(const Field<GRID, SPOT, double> &f,  C_TimeSteps &ts, tImplicitType t = CrankNicholson, bool asjustTimeStep = false)
	{
		const GRID& g = f.Grid();
		int n = g.imax() - g.imin() +1;
		Build(n, ts, t, asjustTimeStep);
	}

	//--------------------------------------------------------------
	//  Time-step control. These correspond to Numerics.TimeStepControl in
	//  the input file.
	//--------------------------------------------------------------
	/**
	 * @brief How many steps the criterion has to keep being met before the
	 *        time step is allowed to grow.
	 */
	void SetAdjustTimeKeepSteps(int s){_adjustTimeKeepSteps = s;}
	/**
	 * @brief Factor the time step is multiplied by when it shrinks.
	 * @param v  the factor; growing the step divides by it instead.
	 */
	void SetAdjustTimeChangeFactor(double v){_adjustTimeChangeFactor = v;}
	/** @brief Local truncation error that is tolerated. */
	void SetAdjustTimeCriterion(double v){_adjustTimeCriterion = v;}
	/**
	 * @brief Floor on the weight used when the error is measured, as a
	 *        fraction of the overall maximum.
	 * @param v  the fraction, clamped into [0, 1].
	 *
	 * See the comment on LocalTruncationCriterion() for the details.
	 */
	void SetAdjustTimeScaleFloor(double v){_adjustTimeScaleFloor = (v < 0 ? 0 : (v > 1 ? 1 : v));}

	/**
	 * @brief Relative rate of change per unit time, used to decide whether
	 *        the solution has reached a steady state.
	 * @return The rate, divided by the time step and so independent of it.
	 *         See the comments in the implementation.
	 */
	double ChangeRateOfVariable();

	/** @brief Which weighting in time is in use (theta = 0 / 0.5 / 1). */
	tImplicitType GetType() const {return _implicit_type;}

	//--------------------------------------------------------------
	//  Moving the unknowns in and out: from a field into the internal
	//  vectors, or the other way round. SetX is used to set the initial
	//  values and to feed the state back in when restarting.
	//--------------------------------------------------------------
	/**
	 * @brief Copy a field, or a plain vector, into the internal unknowns.
	 *
	 * Used to set the initial values and to feed the state back in on
	 * restart. Overloaded per field type so that a mismatch does not compile.
	 */
	void SetX(Vector<double> &f);
	void SetX(ScalarField1D<> &f);
	void SetX(ScalarField2D<> &f);
	void SetX(ScalarField3D<> &f);

	//--------------------------------------------------------------
	//  What a restart file has to carry beyond the unknowns (see
	//  C_RestartState): the unknowns at the start of the last step (the
	//  truncation-error estimate of the adaptive step needs them once the
	//  next step has begun) and the counter of quiet steps before dt may
	//  grow.
	//--------------------------------------------------------------
	const Vector<double>& GetXOld() const {return _phi_old;}
	void SetXOld(const Vector<double> &v) {_phi_old = v;}
	int  AdjustTimeCount() const {return _adjustTimeCount;}
	void SetAdjustTimeCount(int n) {_adjustTimeCount = n;}

	/** @brief Read the current unknowns back out, by element or as a whole. */
	       double   GetX(int i) const {return _phi_new[i];}
	Vector<double>& GetX()      {return _phi_new;}
    void  GetX(Vector<double> &xx){ xx = _phi_new;}

	//------------------------------------------------------------------
	//  Commit one time step without solving again.
	//
	//  Use this after building up a solution through repeated calls to
	//  TimeStepTrial(), to commit that result as the answer for this step.
	//  TimeAdvance() cannot be used when the trials are being driven from
	//  outside, because it solves one more step on top of committing.
	//------------------------------------------------------------------
	/**
	 * @brief Accept the current solution as the answer for this time step,
	 *        without solving again.
	 *
	 * Use this after driving the solution with repeated calls to
	 * TimeStepTrial(). TimeAdvance() cannot be used for that, because it
	 * solves one more step on top of committing.
	 */
	void CommitStep()
	{
		_phi_old2 = _phi_old;
		_phi_old  = _phi_new;
		_crsA.MatVec(_Ax_old, _phi_old);	// the explicit part of the next step
	}

	//  Read-only views of the assembled system, for tests and diagnostics.
	//  GetMatrixA is the operator A of the governing equation, GetMatrixIA
	//  the matrix actually handed to the linear solver.
	const CRS_Matrix& GetMatrixA()  const {return _crsA;}
	const CRS_Matrix& GetMatrixIA() const {return _crsIA;}
	const Vector<double>& GetRHS() const {return _Rhs;}



protected:
	int _nRows;

	C_TimeSteps* _timeStep;

	double _theta;		// how implicit the scheme is, in [0:1]

	bool _constructed;
	bool _initialized;
	bool _adjustTimeStep;
	int _adjustTimeKeepSteps;
	//  Counter for deciding whether dt may be increased. The height equation
	//  and the concentration equation are separate instances, so each holds
	//  its own; a function-local static would let the two interfere.
	int _adjustTimeCount;
	double _adjustTimeChangeFactor;
	double _adjustTimeCriterion;
	double _adjustTimeScaleFloor;


	Vector<double> _phi_old2, _phi_old, _phi_new;	//the unknowns
	Vector<double> _Ax_old;	 //matrix * unknowns

	Vector<double> _Rhs;	//right-hand side of the system that is finally solved
	Vector<double> _VecB;	//the part of the governing equation that does not depend on the unknowns
	CRS_Matrix _crsA;		//coefficient matrix of the part of the governing equation that does depend on them
	CRS_Matrix _crsIA;		//coefficient matrix of the system that is finally solved
	LinearSolver _bcg;		//the linear solver (its implementation is chosen by Backend in the input file)

protected:
	virtual void MakeMatrix_IA() = 0;			//assemble the matrix for the time advance
	virtual void SetMatrixA() = 0;				//assemble the coefficient matrix of the governing equation
	virtual void MakeRHS() = 0;					//assemble the right-hand side vector
	//  A x, for the explicit parts A(phi_old) phi_old and A(phi_guess)
	//  phi_guess of the Newton iteration (x the unknowns A was last
	//  assembled at). The film-height equation of the well-mixed model
	//  evaluates the stiff terms face by face instead (conservative to
	//  rounding of the fluxes); everything else is the product with A.
	virtual void ApplyA(Vector<double> &y, const Vector<double> &x) const { _crsA.MatVec(y, x); }

	virtual double CheckTimeStepCriterion() = 0;
	virtual bool TimeStepAdjustment(double criterion);
	//  A step that failed with the time-step control off: abort the run,
	//  since there is no shorter step to retry it at. The message names
	//  _failCause when the failing part set it.
	void FailedFixedStep(const char *who) const;
	//  What made the last trial fail, for that message: set by the part
	//  that failed (the concentration solve of a coupled step, say), cleared
	//  at the start of every trial. Empty: see Warnings.txt.
	std::string _failCause;

	//------------------------------------------------------------------
	//  The time-step test based on the local truncation error. How the
	//  weighting of the error is decided is explained in the comments in
	//  the .cpp file.
	//
	//  Dividing by zero in a dry region (ho -> 0) or on the first step
	//  (dt0 = 0) would send a NaN through Ek.AbsMax() and wreck the whole
	//  time-step control, so both cases are guarded against.
	//------------------------------------------------------------------
	double LocalTruncationCriterion(const Vector<double> &phi_new) const;

	bool IsBoundary(const int i_row) const {return _crsA.IsBC(i_row);}
	void SetAsBoundary(const int i_row)
	{
		_crsA.SetAsBC(i_row);
		_crsIA.SetAsBC(i_row);
	}

	virtual int Solve(int &iter, double &resid) = 0;			//solve the system

private:
	tImplicitType _implicit_type;

	void Alloc();
	void CheckTheta();
	virtual void SetFlagBC() = 0;				//mark the rows that correspond to boundary conditions
	virtual void AllocAdd();

private:
	C_ImplicitSolver(const C_ImplicitSolver& c);
	C_ImplicitSolver& operator=(const C_ImplicitSolver& c);

};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_IMPLICIT_SOLVER_INCLUDED_
