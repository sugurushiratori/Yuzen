/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  WellMixed_Model_2D.cpp : the common part has moved to the template in WellMixedModelBase.h
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "WellMixed_Model_2D.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

using namespace yuzen::util;
using namespace Global;

void WellMixed_Model_2D::BuildBoundaryCondition(Input &inp, Grid *G)
{
	BoundaryCondition::tBoundaryCondition tXmin, tXmax, tYmin, tYmax;
	tXmin = BoundaryCondition::FromStr(inp.Boundary.XminBC.Type);
	tXmax = BoundaryCondition::FromStr(inp.Boundary.XmaxBC.Type);
	tYmin = BoundaryCondition::FromStr(inp.Boundary.YminBC.Type);
	tYmax = BoundaryCondition::FromStr(inp.Boundary.YmaxBC.Type);
	_boundary.Build(G->nx(), G->ny(), tXmin, tXmax, tYmin, tYmax);
}

void WellMixed_Model_2D::InfoGridExtent() const
{
	Info<<" (Nx, Ny) = ("<<_inputParam->Grid.Nx<<" , "<<_inputParam->Grid.Ny <<")" <<nl;
	Info<<"xrange["<<_inputParam->Grid.Xmin<<" , "<<_inputParam->Grid.Xmax <<"]" <<nl;
	Info<<"yrange["<<_inputParam->Grid.Ymin<<" , "<<_inputParam->Grid.Ymax <<"]" <<nl;
}

//  Non-dimensionalize the probe position given in the input (which carries
//  dimensions) and search from the start of the grid for a cell within one
//  mesh width of it in both x and y. A probe for which no such cell is
//  found is not registered.
void WellMixed_Model_2D::FindProbes(Grid *G)
{
	double px, py;
	for(unsigned int i = 0; i < _inputParam->Output.ProbeList.Probe.size(); i++)
	{
		px = _inputParam->Output.ProbeList.Probe[i].px;
		py = _inputParam->Output.ProbeList.Probe[i].py;
		px /= _NonDim.ExtentScale();
		py /= _NonDim.ExtentScale();
		//  The nearest point within one mesh width in x and in y (not the
		//  first one found, which can be a neighbor)
		int bx = -1, by = -1; double dbest = 1e300;
		for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf)
		{
			const double ex = fabs(G->X(zf) - px), ey = fabs(G->Y(zf) - py);
			if(ex < G->dx() && ey < G->dy())
			{
				const double d = (ex/G->dx())*(ex/G->dx()) + (ey/G->dy())*(ey/G->dy());
				if(d < dbest) { dbest = d; bx = G->ix(zf); by = G->iy(zf); }
			}
		}
		if(bx >= 0) _probeList.push_back(PointProbe(bx, by));
	}
}

void WellMixed_Model_2D::WriteProbeHeader(std::ostringstream &oss, ProbeList::iterator it) const
{
	double xx = _field.GetGrid().XOf<Zface>(it->Ix(),it->Jy());
	double yy = _field.GetGrid().YOf<Zface>(it->Ix(),it->Jy());

	if(_inputParam->Output.Dimensioned)
		oss<<"# (x, y) [mm] = ("<< xx * _NonDim.ExtentScale() *1e3 <<", "<< yy *_NonDim.ExtentScale() *1e3 <<") ";
	else
		oss<<"# (x, y) [-] = ("<< xx <<", "<< yy <<") ";
}
