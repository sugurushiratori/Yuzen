/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  FieldInfo.cpp : implementation of the visitor that sweeps a field and
*                  aggregates volume, mean concentration, maxima and minima
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <fstream>

#include "FieldInfo.h"
#include "LubFields.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

FieldInfo::FieldInfo() :
	_fluidVolume(0),
	_resinVolume(0),
	_averageConc(0),
	_maxH(0), _minH(0),
	_maxC(0), _minC(0),
	_maxHS(0), _minHS(0)
	{}

FieldInfo::FieldInfo(const FieldInfo& f)
{
	this->_fluidVolume = f._fluidVolume;
	this->_resinVolume = f._resinVolume;
	this->_averageConc = f._averageConc;
	this->_maxH  = f._maxH;  this->_minH = f._minH;
	this->_maxC  = f._maxC;  this->_minC = f._minC;
	this->_maxHS = f._maxHS; this->_minHS = f._minHS;
}
FieldInfo& FieldInfo::operator=(const FieldInfo& f)
{
	if(this != &f)
	{
		this->_fluidVolume = f._fluidVolume;
		this->_resinVolume = f._resinVolume;
		this->_averageConc = f._averageConc;
		this->_maxH  = f._maxH;  this->_minH = f._minH;
		this->_maxC  = f._maxC;  this->_minC = f._minC;
		this->_maxHS = f._maxHS; this->_minHS = f._minHS;
	}
	return *this;
}

void FieldInfo::Visit(Field_WM1D &f)
{
	double sumH = 0, aveC = 0, sumHr = 0;
	C_Grid1D *G = f.GetGridPtr();
	double dx = G->dx();
	double co, ho, so;
	_maxH = _maxHS = _maxC = -1e40;
	_minH = _minHS = _minC = 1e40;

	Xedge xe;

	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		xe = G->Of<Xedge>(ii);
		if(G->CheckGhost(xe) != G->Inside) continue;

		co = f.CC(xe);
		ho = f.HH(xe);
		so = f.SS(xe);
		//  The volume is the sum of the cells, as in 2-D: it is what the
		//  finite-volume film equation conserves. The mean of the two face
		//  values, 0.5*(hL+hR) dx, would bring the ghosts in at the edges,
		//  dx/4 (h_ghost - h_first): nothing at a Symmetry edge (h_ghost =
		//  h_first), but -h_first dx/2 at a Pinning one (h_ghost = -h_first),
		//  so VolErr of a pinned film would move with its first cell (~5e-5
		//  for an EdgeBump, tools/proptest.py) while the film lost nothing.
		sumH += ho*dx;
		aveC += co*ho*dx;
		sumHr+= ho*co*dx;

		_maxH  = Math::Max(_maxH, ho);
		_minH  = Math::Min(_minH, ho);
		_maxHS = Math::Max(_maxHS, ho+so);
		_minHS = Math::Min(_minHS, ho+so);
		_maxC  = Math::Max(_maxC, co);
		_minC  = Math::Min(_minC, co);
	}

	_fluidVolume = sumH;
	_resinVolume = sumHr;
	_averageConc = aveC / _fluidVolume;
}

void FieldInfo::Visit(Field_WM2D &f)
{
	double sumH = 0, aveC = 0, sumHr = 0;
	C_Grid2D *G = f.GetGridPtr();
	double dx = G->dx();
	double dy = G->dy();
	double cc, hh, ss, hr;
	_maxH = _maxHS = _maxC = -1e40;
	_minH = _minHS = _minC = 1e40;

	Zface zf;

	for(int ii = G->imin(); ii <= G->imax(); ii++)
	{
		zf = G->Of<Zface>(ii);
		if(G->CheckGhost(zf) != G->Inside) continue;

		cc = f.CC(zf);
		hh = f.HH(zf);
		ss = f.SS(zf);
		hr = hh*cc;
		sumH += hh*dx*dy;
		aveC += cc*hh*dx*dy;
		sumHr+= hr*dx*dy;

		_maxH  = Math::Max(_maxH, hh);
		_minH  = Math::Min(_minH, hh);
		_maxHS = Math::Max(_maxHS, hh+ss);
		_minHS = Math::Min(_minHS, hh+ss);
		_maxC  = Math::Max(_maxC, cc);
		_minC  = Math::Min(_minC, cc);
	}

	_fluidVolume = sumH;
	_resinVolume = sumHr;
	_averageConc = aveC / _fluidVolume;
}

void FieldInfo::Visit(Field_CZ &f)
{
	C_Grid2D *G2 = f.GetGridPtr2();
	C_Grid3D *G3 = f.GetGridPtr3();
	double dx = G3->dx();
	//  A 1-D film (Ny = 1) is a strip of unit width, as in the well-mixed
	//  model and in the equations of ConcZ: the grid's dy there is the
	//  Ymin..Ymax of the input, which a 1-D film does not use (it made the
	//  volumes of a 1-D ConcZ run (Ymax - Ymin) / L0 times those of the
	//  well-mixed model for the same film, 1100 times in concz_1d_reg)
	double dy = (G2->ny() > 1) ? G3->dy() : 1.0;
	double dz = G3->dz();
	double cc, hh, ss;
	_maxH = _maxHS = _maxC = -1e40;
	_minH = _minHS = _minC = 1e40;
	int ix,jy,kz;
	double sumH = 0;
	double sumC, sumHr=0;
	for(Zface zf = G2->First<Zface>(); zf != G2->Last<Zface>(); ++zf)
	{
		if(G2->CheckGhost(zf) != G2->Inside) continue;
		ix = G2->ix(zf);
		jy = G2->iy(zf);
		hh = f.HH(zf);
		ss = f.SS(zf);
		sumH += hh*dx*dy;
		_maxH  = Math::Max(_maxH, hh);
		_minH  = Math::Min(_minH, hh);
		_maxHS = Math::Max(_maxHS, hh+ss);
		_minHS = Math::Min(_minHS, hh+ss);

		sumC=0;
		for(kz = 2; kz < G3->nz()-2; kz++)
		{
			cc= f.CC(G3->At<Cell>(ix,jy,kz));
			sumC += cc *dz;
			_maxC  = Math::Max(_maxC, cc);
			_minC  = Math::Min(_minC, cc);
		}
		sumHr += sumC *hh*dx*dy;
	}

	_fluidVolume = sumH;
	_resinVolume = sumHr;
	_averageConc = sumHr / sumH;

}
