/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Probe.h : a probe that records the time series at a given point
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _PROBE_INCLUDED_
#define _PROBE_INCLUDED_

#include "Grids.h"
#include "Input.h"
#include "FieldVisitor.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class Field_WM;
class Field_CZ;

//--------------------------------------------------------------------------
// Records the values at one grid point so they can be written as a time
// series. One instance per probe listed in the input file.
//--------------------------------------------------------------------------
class PointProbe : public FieldVisitor
{
public:
	/**
	 * @brief Builds a probe watching one grid point.
	 * @param ix  x grid index of the point. The input file gives coordinates;
	 *            they are rounded to the nearest grid point before they get here.
	 * @param jy  y grid index; ignored in the 1-D case
	 * @param kz  z grid index; ignored in the 1-D and 2-D cases
	 *
	 * The recorded values start at zero, so reading them before the first
	 * Visit gives zero rather than an undefined value.
	 */
	// The value members (_hh through _uz) are all initialized here, so that
	// reading them before Visit() is defined.
	PointProbe(int ix, int jy = 0, int kz = 0)
		: _ix(ix), _jy(jy), _kz(kz),
		  _hh(0), _ss(0), _cc(0), _ux(0), _uy(0), _uz(0)
	{
	}
	/**
	 * @brief Picks the values at the watched point out of the field and keeps
	 *        them. Called through field.Accept(probe); one overload per kind
	 *        of field. Each call replaces the values from the previous one.
	 */
	//  Picks the values at that grid point out of the field and keeps them
	virtual void Visit(Field_WM1D &);
	virtual void Visit(Field_WM2D &);
	virtual void Visit(Field_CZ &);

	// Grid indices of the watched point, fixed at construction.
	//  Indices of the grid point being watched. The input file gives a
	//  coordinate, which is rounded to the nearest grid point at Build time.
	int Ix() const {return _ix;}
	int Jy() const {return _jy;}
	int Kz() const {return _kz;}

	//--------------------------------------------------------------
	// Values picked up at the watched point. All non-dimensional, and only
	// meaningful once Visit has run; before that they are zero.
	//--------------------------------------------------------------
	//  Values picked up (non-dimensional). Read them after calling Visit.
	//--------------------------------------------------------------
	/** @brief Film thickness h at the watched point. */
	double HH() const {return _hh;}			// film thickness
	/** @brief Substrate shape s at the watched point. */
	double SS() const {return _ss;}			// substrate shape
	/** @brief Height of the free surface, h + s. */
	double HS() const {return _hh + _ss;}	// height of the free surface
	/** @brief Concentration c at the watched point. */
	double CC() const {return _cc;}			// concentration
	/**
	 * @brief Components of the combined advection velocity at the watched point.
	 *
	 * This is the mixed velocity the concentration equation advects with: the
	 * velocities of the individual physical terms added together, not the
	 * velocity of any one term.
	 */
	double Ux() const {return _ux;}			// combined (mixed) velocity
	double Uy() const {return _uy;}
	double Uz() const {return _uz;}



private:
	int _ix, _jy, _kz;
	double _hh;
	double _ss;
	double _cc;
	double _ux;
	double _uy;
	double _uz;

	// The default copy constructor is allowed
	//PointProbe(const PointProbe&);

	void operator=(const PointProbe&);
};

/** The probes of one run, in the order they were listed in the input file. */
typedef std::list<PointProbe> ProbeList;

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_PROBE_INCLUDED_
