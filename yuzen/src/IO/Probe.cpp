/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Probe.cpp : implementation of the visitor that picks the values at the probe
*              point out of the field
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <fstream>

#include "Probe.h"
#include "LubFields.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif


void PointProbe::Visit(Field_WM1D &f)
{
	Xedge xe = f.GetGrid().At<Xedge>(_ix);
	_hh = f.HH(xe);
	_ss = f.SS(xe);
	_cc = f.CC(xe);
	_ux = f.Ux(f.GetGrid().At<Node>(_ix));
	_uy = 0;
	_uz = 0;
}


void PointProbe::Visit(Field_WM2D &f)
{
	Zface zf = f.GetGrid().At<Zface>(_ix, _jy);
	Xedge xe = f.GetGrid().At<Xedge>(_ix, _jy);
	Yedge ye = f.GetGrid().At<Yedge>(_ix, _jy);
	_hh = f.HH(zf);
	_ss = f.SS(zf);
	_cc = f.CC(zf);
	_ux = f.Ux(ye);
	_uy = f.Uy(xe);
	_uz = 0;
}

void PointProbe::Visit(Field_CZ &f)
{
	Zface z1 = f.GetGrid2().At<Zface>(_ix, _jy);
	Cell  cl = f.GetGrid3().At<Cell>(_ix, _jy, _kz);
	Xface xf = f.GetGrid3().At<Xface>(_ix, _jy, _kz);
	Yface yf = f.GetGrid3().At<Yface>(_ix, _jy, _kz);
	Zface zf = f.GetGrid3().At<Zface>(_ix, _jy, _kz);

	_hh = f.HH(z1);
	_ss = f.SS(z1);
	_cc = f.CC(cl);
	_ux = f.Ux(xf);
	_uy = f.Uy(yf);
	_uz = f.Uz(zf);
}
