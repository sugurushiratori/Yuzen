/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Source.cpp : the implementation of the source (inlet) term
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/



#include "YuzenUtil.h"
#include "Source.h"
#include <cmath>

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

using namespace yuzen::util;
using namespace Global;

void Source_WM2D::Build(Input::I_Physics &inp_phys, C_Grid2D *g, NonDimensionalNumbers *nonDim, C_TimeSteps *ts)
{
	_G = g;

	_fraction.Build(*_G);
	// updateFraction() is skipped while the source is disabled, so the
	// field is explicitly zeroed here
	for(Zface zf = _G->First<Zface>(); zf != _G->Last<Zface>(); ++zf) _fraction(zf) = 0.0;

	BuildParams(inp_phys, nonDim, ts);

	updateFraction();

	if(_enabled)
		Sout << "number of sub-regions for calculation of source term fraction: " << _numSubRegion << nl;
}

void Source_WM2D::updateFraction()
{
	// While the source is disabled operator() always returns 0, so the walk
	//          over the whole grid is skipped (it costs about 8% of the run
	//          time of a two-dimensional case).
	if(!_enabled) return;

	double xn, yn;
	NozzlePos(xn, yn);

	double dx = _G->dx();
	double dy = _G->dy();

	for(Zface zf = _G->First<Zface>(); zf != _G->Last<Zface>(); ++zf)
	{
		double xc = _G->X(zf);
		double yc = _G->Y(zf);
		double rr = sqrt( pow(xc - xn, 2) + pow(yc - yn, 2));

		//  A cell the nozzle reaches has its center within R and half its
		//  diagonal of the nozzle's center; the test allows a whole diagonal.
		//  (It allowed 2 min(dx, dy), which cells more than 3.9 times longer
		//  than wide overreach: analytic.py, nozzle_area.)
		if(rr < _Radius + std::sqrt(dx * dx + dy * dy))
		{
			double alp = calFraction(zf, xn, yn);
			_fraction(zf) = alp;
		}
		else
		{
			_fraction(zf) = 0.0;
		}

	}
}

//------------------------------------------------------------------
//  Divide the cell into _numSubRegion x _numSubRegion equal parts and count
//  the fraction of those sub-regions that fall within the nozzle's circle,
//  each tested at its center (its lower-left corner would move the supply
//  half a part up and to the right of the nozzle).
//  This is the two-dimensional version of the 1-D width ratio
//  (Source_WM1D::calFraction); here it is an area ratio.
//------------------------------------------------------------------
double Source_WM2D::calFraction(const Zface& zf, double Xn, double Yn)
{
	double xc = _G->X(zf);
	double yc = _G->Y(zf);
	double dx = _G->dx();
	double dy = _G->dy();

	double dx_sub = dx / (double)_numSubRegion;
	double dy_sub = dy / (double)_numSubRegion;

	double xL = xc - 0.5 * dx;
	double yB = yc - 0.5 * dy;

	double alpha = 0.0;
	for(int ix = 0; ix < _numSubRegion; ix++)
	{
		double xx = xL + ((double)ix + 0.5) * dx_sub;
		for(int jy = 0; jy < _numSubRegion; jy++)
		{
			double yy = yB + ((double)jy + 0.5) * dy_sub;
			// distance from nozzle center
			double rr = sqrt( pow(xx - Xn, 2) + pow(yy - Yn, 2));
			double Fs = (rr <= _Radius ? 1.0 : 0.0);

			alpha += Fs;
		}
	}

	alpha /= (double)(_numSubRegion * _numSubRegion);

	return alpha;
}

double Source_WM2D::operator()(const Zface& zf) const
{
	if (!_enabled) return 0.0;

	double phi = FlowRateNow();

	double source = phi * _fraction(zf);

	return source;
}
