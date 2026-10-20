/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ShapeFunctions.cpp : implementation of the shape functions
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#include "ShapeFunctions.h"
#include "MyMath.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif


double RectProfile_base(double point, double center, double width, double slope)
{
	return Math::Sigmoid(point, center-0.5*width,slope*0.1)*(1.0-Math::Sigmoid(point, center+0.5*width,slope*0.1));
}


//-------------------
//   RectBump
//-------------------
void C_RectBump::Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double widthX, double widthY, double slopeX, double slopeY)
{
	_order	= order;
	_CreateType = type;
	_base	= base;
	_amp	= amp;
	_cntX	= cntX;		_cntY	= cntY;
	_widthX	= widthX;	_widthY	= widthY;
	_slopeX	= slopeX;	_slopeY	= slopeY;
}
double C_RectBump::operator()(double px, double py) const
{
	double val = 0;
	val = RectProfile_base(px, _cntX, _widthX, _slopeX)
		 *RectProfile_base(py, _cntY, _widthY, _slopeY);
	return _amp * val + _base;
}

//-------------------
//   EllipseBump
//-------------------
void C_EllipseBump::Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double widthX, double widthY, double slope)
{
	_order	= order;
	_CreateType = type;
	_base	= base;
	_amp	= amp;
	_cntX	= cntX;		_cntY	= cntY;
	_widthX	= widthX;	_widthY	= widthY;
	_slope	= slope;
}
double C_EllipseBump::operator()(double px, double py) const
{
	double rr = sqrt(pow((px-_cntX)/_widthX, 2.0) + pow((py-_cntY)/_widthY, 2.0));
	double val = RectProfile_base(rr, 0, 1.0, _slope);
	return _amp * val + _base;
}

//----------------------
//     FlatPlane
//----------------------
void C_FlatPlane::Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double angleX, double angleY)
{
	_order	= order;
	_CreateType = type;
	_base	= base;
	_amp	= amp;
	_cntX	= cntX;		_cntY	= cntY;
	_angleX	= angleX;	_angleY	= angleY;
}
double C_FlatPlane::operator()(double px, double py) const
{
	return (px-_cntX)*tan(_angleX)+(py-_cntY)*tan(_angleY) + _amp;
}


//---------------------------------------------
//   Mexican-hat (1-(x/L)^2) * exp(-b*(x/L)^2)
//---------------------------------------------
void C_MexicanHat::Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double Lx, double Ly, double beta)
{
	_order	= order;
	_CreateType = type;
	_base	= base;
	_amp	= amp;
	_cntX	= cntX;		_cntY	= cntY;
	_Lx		= Lx;		_Ly		= Ly;
	_beta	= beta;
}
double C_MexicanHat::operator()(double px, double py) const
{
	double val = 0;
	if(Math::IsZero(_Lx) || Math::IsZero(_Ly)) return 0;
	val = (1.0 -pow((px-_cntX)/_Lx, 2.0) -pow((py-_cntY)/_Ly, 2.0) )
		* exp(-_beta * (pow((px-_cntX)/_Lx, 2.0)+pow((py-_cntY)/_Ly, 2.0)));
	return _amp * val + _base;
}


//--------------------------------
// Smooth convex (1-(x/L)^2)^2
//--------------------------------
void C_SmoothConvex::Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double Lx, double Ly)
{
	_order	= order;
	_CreateType = type;
	_base	= base;
	_amp	= amp;
	_cntX	= cntX;		_cntY	= cntY;
	_Lx		= Lx;		_Ly		= Ly;
}
double C_SmoothConvex::operator()(double px, double py) const
{
	double val = 0;
	if(sqrt(pow((px-_cntX)/_Lx,2.0) +pow((py-_cntY)/_Ly,2.0)) > 1.0) return _base;
	val = pow(1.0 -pow((px-_cntX)/_Lx, 2.0) -pow((py-_cntY)/_Ly, 2.0), 2.0);
	return _amp * val + _base;
}

//--------------------------------
//   PeriodicWave
//--------------------------------
void C_PeriodicWave::Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double Lx, double Ly)
{
	_order	= order;
	_CreateType = type;
	_base	= base;
	_amp	= amp;
	_cntX	= cntX;		_cntY	= cntY;
	_Lx		= Lx;		_Ly		= Ly;
}
double C_PeriodicWave::operator()(double px, double py) const
{
	double val = 0;
	val =  cos(2.0*Math::pi()* (px-_cntX)/_Lx) * cos(2.0*Math::pi() * (py-_cntY)/_Ly);
	return _amp * val + _base;
}


//--------------------------------
//   EdgeBump
//--------------------------------
void C_EdgeBump::Build(int order, tCreateType type, double base, double amp,
		double alpha, double Lin, double Lout, double rs, tLocation loc)
{
	_order		= order;
	_CreateType = type;
	_base		= base;
	_amp		= amp;
	_Lin		= Lin;
	_Lout		= Lout;
	_alpha		= alpha;
	_rs			= rs;
	_location	= loc;
}
double C_EdgeBump::operator()(double px, double py) const
{
	double ff, gg, sig;
	double aa = _amp-_base;
	double eps = 0.1 * _Lin;
	double rr, pp;

	switch(_location)
	{
	default:
	case Xmin:	rr = px; pp = -(rr - _rs); break;
	case Xmax:	rr = px; pp = +(rr - _rs); break;
	case Ymin:	rr = py; pp = -(rr - _rs); break;
	case Ymax:	rr = py; pp = +(rr - _rs); break;
	}

	ff  = aa*(1.0 -pow(pp/_Lin , 2.0)) * exp(-1.0/_alpha *pow(pp/_Lin, 2.0)) + _base;
	//  The outer part falls from Height at the switch-over point to zero at
	//  the pinned edge: aa (1 - q (1 + base/aa)) + base with q = (pp/Lout)^2,
	//  written so that it holds for Height = Base too (a film flat at Base
	//  up to the meniscus), where the form with base/aa divides by zero.
	gg  = _amp*(1.0 -pow(pp/_Lout, 2.0));
	sig = Math::Sigmoid(pp, 0, eps);

	return (1.0 -sig)*ff + sig*gg;
}
