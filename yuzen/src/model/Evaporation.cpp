/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Evaporation.cpp : implementation of the evaporation rate models
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#include "Evaporation.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif




C_Evaporation::C_Evaporation()
{
	_enabled = false;
	_sw_Hast_modify = false;
	_M_resin = 1;
	_M_solvent = 1;
	_E0 = 1;
	_NegPowExponent = 0.25;
	_NegPowRmax  = 1;
	_ConcR      = 0.0;
	_powerN     = 1;
	_h_ast      = 0.0;
	_t_ConcFunc    = tc_Raoult;
	_t_SpatialFunc = ts_Uniform;
}

void C_Evaporation::Build(const Input &inp, double scale)
{
	_enabled = inp.Physics.Evaporation.Enable;
	_M_resin = inp.Fluid.MolarMass.Resin;
	_M_solvent = inp.Fluid.MolarMass.Solvent;
	_E0 = inp.Physics.Evaporation.Coef;

	_powerN = inp.Physics.Evaporation.PowerN;
	_NegPowExponent = inp.Physics.Evaporation.NegPowExponent;
	_NegPowRmax  = inp.Physics.Evaporation.NegPowRmax;

	_ConcR      = inp.Physics.Evaporation.ConcR;

	Input inp2 = inp;
	inp2.Grid.Xmin /= scale;
	inp2.Grid.Xmax /= scale;
	inp2.Grid.Ymin /= scale;
	inp2.Grid.Ymax /= scale;
	_G.Build(inp2.Grid);
	_NegPowRmax /= scale;

	_h_ast = inp.Initial.SurfaceShape.PrecursorRatio;

	if(Math::IsZero(_h_ast)) _sw_Hast_modify = false;
	else _sw_Hast_modify = true;

	if     (inp.Physics.Evaporation.ConcFunc.compare("Raoult") == 0) _t_ConcFunc = tc_Raoult;
	else if(inp.Physics.Evaporation.ConcFunc.compare("PowerN") == 0) _t_ConcFunc = tc_PowerN;
	else _t_ConcFunc = tc_Raoult;

	if     (inp.Physics.Evaporation.SpatialFunc.compare("uniform") == 0) _t_SpatialFunc = ts_Uniform;
	else if(inp.Physics.Evaporation.SpatialFunc.compare("negPow") == 0)  _t_SpatialFunc = ts_NegPow;
	else _t_SpatialFunc = ts_Uniform;
}

double C_Evaporation::Evaluate(double hh, double cc, double xx, double yy) const
{
	double f_spatial, f_conc, f_modify = 1;

	if(_sw_Hast_modify)
		f_modify = ModifyingCoef(hh);

	switch(_t_ConcFunc)
	{
	default:
	case tc_Raoult:
		f_conc = (1.0 - CalChi(cc));
		break;
	case tc_PowerN:
		f_conc = PowerN(cc);
		break;
	}

	switch(_t_SpatialFunc)
	{
	default:
	case ts_Uniform:
		f_spatial = 1.0;
		break;
	case ts_NegPow:
		//  Blend smoothly near the boundary: from rmx-eps a quarter cosine
		//  that leaves the curve with its value and (backward-difference)
		//  slope and levels off at rmx, then constant beyond. The cosine is
		//  not taken past rmx: there the rate would fall again and oscillate
		//  with a period of 4 eps over the rest of the domain (for
		//  NegPowExponent = 0.25, to 0.75 of its rim value at rmx+2 eps).
		double rr  = sqrt(xx*xx + yy*yy);
		double rmx = _NegPowRmax;
		double dr = Math::Min(_G.dx(), _G.dy());
		double eps = dr * 5.0;
		double fRE = NegPowFunc(cc, rmx-eps);
		double dU =  (fRE -NegPowFunc(cc,rmx-eps-dr))/dr;
		double aa = 2.0*eps / Math::pi() * dU;
		double f0 = NegPowFunc(cc,rr);
		double f1 = aa * cos(-(rmx-Math::Min(rr, rmx))/eps*0.5*Math::pi()) + fRE;
		f_spatial = (rr < rmx-eps ? f0 : f1);
		break;
	}

	return _E0 * f_conc * f_spatial * f_modify;
}
double C_Evaporation::operator()(Field_WM1D &ff, Xedge xe) const
{
	if(!_enabled) return 0;

	C_Grid1D *G = ff.GetGridPtr();
	double xx, cc, hh;
	xx = G->X(xe);

	double yy = 0;	//////////////

	cc = ff.CC(xe);
	hh = ff.HH(xe);
	return Evaluate(hh, cc, xx, yy);
}
double C_Evaporation::operator()(Field_WM2D &ff, Zface zf) const
{
	if(!_enabled) return 0;

	C_Grid2D *G = ff.GetGridPtr();
	double xx, yy, cc, hh;
	xx = G->X(zf);
	yy = G->Y(zf);
	cc = ff.CC(zf);
	hh = ff.HH(zf);
	return Evaluate(hh, cc, xx, yy);
}
double C_Evaporation::operator()(Field_CZ &ff, Zface zf) const
{
	if(!_enabled) return 0;

	C_Grid2D *G2 = ff.GetGridPtr2();
	C_Grid3D *G3 = ff.GetGridPtr3();
	double xx, yy, cc, hh;
	xx = G2->X(zf);
	yy = G2->Y(zf);

	int ix, jy, kz;
	//  The concentration at the free surface: the mean of the top real cell
	//  and the ghost above it. KmaxReal() is one past the last real cell
	//  (the range is half-open), so the pair is KmaxReal()-1, KmaxReal(),
	//  not KmaxReal(), KmaxReal()+1 -- the two ghosts, a whole cell above
	//  the surface.
	kz = G3->KmaxReal() - 1;
	ix = G2->ix(zf);
	jy = G2->iy(zf);
	Cell cl1, cl2;
	cl1 = G3->At<Cell>(ix, jy, kz);
	cl2 = G3->At<Cell>(ix, jy, kz+1);
	cc = 0.5*(ff.CC(cl1)+ff.CC(cl2));
	hh = ff.HH(zf);
	return Evaluate(hh, cc, xx, yy);
}


//  The mole fraction of the solute. c is a mass fraction and bounded by 1;
//  past it (a ghost value, the surface value of ConcZ extrapolated from the
//  top cell, or rounding) chi would come out above 1 and Raoult's law give
//  a negative rate: condensation onto a film that has no solvent left. c is
//  taken within [0, 1], so the rate is 0 from c = 1 on.
double C_Evaporation::CalChi(const double c) const
{
	const double cc = (c < 0.0) ? 0.0 : ((c > 1.0) ? 1.0 : c);
	return (cc/_M_resin) / (cc/_M_resin + (1.0-cc)/_M_solvent);
}

double C_Evaporation::PowerN(const double cc) const
{
	double val = 1.0 - cc - _ConcR;
	val = Math::Max(val, 0.0);
	return pow(val, _powerN);
}

//-------------------------------------------------
//  Radial negative power: k (1 - r/rmax')^(-n), k = (1 - n)(2 - n)/2
//-------------------------------------------------
double C_Evaporation::NegPowFunc(const double cc, double rr) const
{
	//  Cap the value so that it does not blow up.
	double nn, kk, rev_r, val;
	double max_val = 10;

	double rmx = _NegPowRmax;
	double EpsRmx = 1e-3;

	//  k makes the mean over the disc r < rmax' equal to 1.
	nn = _NegPowExponent;
	kk = 0.5*(1.0-nn)*(2.0-nn);
	rev_r = 1.0 - rr/(rmx *(1.0+EpsRmx));

	if( rev_r < 0.0 || Math::IsLessEpsilon(rev_r)) rev_r = 1e-10;
	val = kk * pow(rev_r,-nn);
	return Math::Min(val, max_val);
}

double C_Evaporation::ModifyingCoef(double hh) const
{
	double aa, val;
	aa = 20.0;
	if(Math::IsZero(_h_ast)) val = 1.0;
	else val = (hh - _h_ast) / ((aa-1.0)*_h_ast);
	//  Zero from the precursor down: below it the ramp would turn negative
	//  and condense solvent onto the film out of nothing.
	return Math::Max(0.0, Math::Min(val, 1.0));
}
