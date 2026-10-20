/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  SourceBase.h : the dimension-independent part of the liquid source
*                 (the nozzle)
*
*  What 1-D and 2-D have in common is kept here:
*    - reading the input and non-dimensionalizing it
*    - the position of the nozzle (the tip of the rotating arm)
*    - the ramp-up of the flow rate (an exponential relaxation with time
*      constant tau)
*  What differs between the dimensions,
*    - the part that walks the grid and computes each cell's covered
*      fraction _fraction
*  lives in the derived classes (Source_WM1D / Source_WM2D).
*
*  The nozzle position in 1-D is the 2-D expression with y dropped, that is
*      xn = ArmX + ArmLength * cos(theta)
*  With ArmLength = 0 the nozzle is fixed at ArmX, so a "nozzle that does
*  not move" can be written with the same input.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _SOURCE_BASE_INCLUDED_
#define _SOURCE_BASE_INCLUDED_

#include "YuzenUtil.h"
#include "NonDimensionalNumbers.h"
#include "TimeSteps.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class SourceBase
{
public:
	/** @brief Zero every nozzle parameter; BuildParams() reads the real ones. */
	SourceBase()
	 : _enabled(false), _NonDim(0), _timeStep(0),
	   _Radius(0), _flowRate(0), _Conc(0), _concSupplied(false), _timeConstant(0),
	   _armLength(0), _armX(0), _armY(0), _rotationRate(0), _startAngle(0),
	   _slopeTime(0), _mainTime(0), _numSubRegion(1), _evalTime(-1), _useEvalTime(false) {}
	virtual ~SourceBase(){}

	/** @brief Whether a source is switched on in the input. */
	bool   Enabled() const { return _enabled; }

	/**
	 * @brief Concentration of the supplied liquid (Physics.Source.Conc).
	 *
	 * Enters the concentration equation only when ConcSupplied() holds; with
	 * the default ConcType = "local" the liquid arrives at the concentration
	 * of the film where it lands, and this value is not used.
	 */
	double Conc() const { return _Conc; }
	/**
	 * @brief Whether the supplied liquid arrives at Conc()
	 *        (Physics.Source.ConcType = "supplied") rather than at the local
	 *        concentration ("local"): the solute balance of the well-mixed
	 *        concentration equation then takes S Conc() of solute with the
	 *        supply S (ConcStep.h).
	 */
	bool ConcSupplied() const { return _concSupplied; }

	/**
	 * @brief The time at which the nozzle position and the flow rate are
	 *        taken (dimensionless). The film-height equation sets the middle
	 *        of the step, so that the supply enters its Crank-Nicolson update
	 *        at second order; until it is set, the current time is used.
	 */
	void SetEvalTime(double t) { _evalTime = t; _useEvalTime = true; }

protected:
	double EvalTime() const;
	/**
	 * @brief Read the input and non-dimensionalize it. Called from the
	 *        derived class' Build.
	 */
	void BuildParams(Input::I_Physics &inp_phys, NonDimensionalNumbers *nonDim, C_TimeSteps *ts);

	/**
	 * @brief The nozzle position at the current time (the tip of the
	 *        rotating arm).
	 */
	void NozzlePos(double &xn, double &yn) const;

	/**
	 * @brief The flow rate at the current time, with the ramp-up of time
	 *        constant tau applied.
	 */
	double FlowRateNow() const;

	bool _enabled;
	NonDimensionalNumbers* _NonDim;
	C_TimeSteps* _timeStep;

	//  The nozzle parameters, read from the input and held in
	//  non-dimensional form.
	//    _Radius        the radius of the nozzle opening
	//    _flowRate      the steady-state supply flow rate
	//    _Conc          the concentration of the supplied liquid
	//    _concSupplied  whether the liquid arrives at _Conc (ConcType =
	//                   "supplied") or at the local concentration ("local")
	//    _timeConstant  the time constant tau of the flow-rate ramp-up
	//    _armLength     the length of the rotating arm (0 means a nozzle
	//                   that does not move)
	//    _armX, _armY   the center of the rotating arm
	//    _rotationRate  the rotation rate of the arm
	//    _startAngle    the initial angle of the arm
	//    _slopeTime     the time the rotation rate takes to ramp up and to
	//                   ramp down
	//    _mainTime      the time the rotation rate stays constant
	//    _numSubRegion  the number of sub-divisions within a cell used when
	//                   measuring its covered fraction
	double _Radius;
	double _flowRate;
	double _Conc;
	bool   _concSupplied;
	double _timeConstant;

	double _armLength;
	double _armX, _armY;
	double _rotationRate;
	double _startAngle;
	double _slopeTime;
	double _mainTime;

	int _numSubRegion;

	double _evalTime;		// see SetEvalTime
	bool   _useEvalTime;

private:
	SourceBase(const SourceBase&);
	void operator=(const SourceBase&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_SOURCE_BASE_INCLUDED_
