/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Source_1D.cpp : the 1-D version of the source. For the details, see the
*                  top of Source_1D.h
*
*  It is laid out the same way as the 2-D version (WM2D/Source.cpp).
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "Source_1D.h"
#include <cmath>

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

using namespace yuzen::util;
using namespace Global;

void Source_WM1D::Build(Input::I_Physics &inp_phys, C_Grid1D *g, NonDimensionalNumbers *nonDim, C_TimeSteps *ts)
{
	_G = g;

	_fraction.Build(*_G);
	// updateFraction() is skipped while the source is disabled, so the
	// field is explicitly zeroed here
	for(Xedge xe = _G->First<Xedge>(); xe != _G->Last<Xedge>(); ++xe) _fraction(xe) = 0.0;

	BuildParams(inp_phys, nonDim, ts);

	updateFraction();

	if(_enabled)
		Sout << "number of sub-regions for calculation of source term fraction: " << _numSubRegion << nl;
}

void Source_WM1D::updateFraction()
{
	// While the source is disabled operator() always returns 0, so the
	// computation itself is skipped
	if(!_enabled) return;

	double xn, yn;
	NozzlePos(xn, yn);	// the y component is unused in 1-D

	double dx = _G->dx();

	for(Xedge xe = _G->First<Xedge>(); xe != _G->Last<Xedge>(); ++xe)
	{
		double xc = _G->X(xe);
		double rr = fabs(xc - xn);

		if(rr < _Radius + 2.0 * dx)
		{
			_fraction(xe) = calFraction(xe, xn);
		}
		else
		{
			_fraction(xe) = 0.0;
		}
	}
}

//------------------------------------------------------------------
//  Divide the cell's width into _numSubRegion equal parts and count the
//  fraction that falls within the interval the nozzle covers. Each part is
//  tested at its center: tested at its left end, it would move the supply
//  dx/(2 SubRegions) to the right of the nozzle -- half a cell with the
//  default SubRegions = 1 -- so that a nozzle at the center of a
//  symmetric film fed it asymmetrically.
//  This is the one-dimensional version of the 2-D area ratio
//  (calFraction).
//------------------------------------------------------------------
double Source_WM1D::calFraction(const Xedge& xe, double Xn)
{
	double xc = _G->X(xe);
	double dx = _G->dx();

	double dx_sub = dx / (double)_numSubRegion;
	double xL = xc - 0.5 * dx;

	double alpha = 0.0;
	for(int ix = 0; ix < _numSubRegion; ix++)
	{
		double xx = xL + ((double)ix + 0.5) * dx_sub;
		// distance from nozzle center
		double rr = fabs(xx - Xn);
		double Fs = (rr <= _Radius ? 1.0 : 0.0);

		alpha += Fs;
	}

	alpha /= (double)_numSubRegion;

	return alpha;
}

double Source_WM1D::operator()(const Xedge& xe) const
{
	if (!_enabled) return 0.0;

	double phi = FlowRateNow();

	double source = phi * _fraction(xe);

	return source;
}
