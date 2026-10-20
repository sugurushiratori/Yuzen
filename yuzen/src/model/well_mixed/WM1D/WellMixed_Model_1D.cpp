/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  WellMixed_Model_1D.cpp : the common part has moved to the template in WellMixedModelBase.h
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "WellMixed_Model_1D.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

using namespace yuzen::util;
using namespace Global;

void WellMixed_Model_1D::BuildBoundaryCondition(Input &inp, Grid *G)
{
	BoundaryCondition::tBoundaryCondition tXmin, tXmax;
	tXmin = BoundaryCondition::FromStr(inp.Boundary.XminBC.Type);
	tXmax = BoundaryCondition::FromStr(inp.Boundary.XmaxBC.Type);
	_boundary.Build(G->nx(), tXmin, tXmax);
}

void WellMixed_Model_1D::InfoGridExtent() const
{
	Info<<" Nx = "<<_inputParam->Grid.Nx <<nl;
	Info<<"xrange["<<_inputParam->Grid.Xmin<<" , "<<_inputParam->Grid.Xmax <<"]" <<nl;
}

//  Non-dimensionalize the probe position given in the input (which carries
//  dimensions) and search from the start of the grid for a cell within one
//  mesh width of it in x. A probe for which no such cell is found is not
//  registered.
void WellMixed_Model_1D::FindProbes(Grid *G)
{
	double px;
	for(unsigned int i = 0; i < _inputParam->Output.ProbeList.Probe.size(); i++)
	{
		px = _inputParam->Output.ProbeList.Probe[i].px;
		px /= _NonDim.ExtentScale();
		//  The nearest point within one mesh width, not the first one found
		//  within a mesh width, which can be the neighbor on the far side of
		//  the probe.
		int best = -1; double dbest = G->dx();
		for(Xedge xe = G->First<Xedge>(); xe != G->Last<Xedge>(); ++xe)
		{
			const double d = fabs(G->X(xe) - px);
			if(d < dbest) { dbest = d; best = G->ix(xe); }
		}
		if(best >= 0) _probeList.push_back(PointProbe(best));
	}
}

void WellMixed_Model_1D::WriteProbeHeader(std::ostringstream &oss, ProbeList::iterator it) const
{
	Grid *G = const_cast<Field&>(_field).GetGridPtr();
	//  Ix() is a grid coordinate (At<>), not a linear index (Of<>): with
	//  Of<> the header named the point one mesh width short of the probe.
	double xx = G->XOf<Xedge>(it->Ix());

	if(_inputParam->Output.Dimensioned)
		oss<<"# x [mm] = ("<< xx * _NonDim.ExtentScale() *1e3 <<") ";
	else
		oss<<"# x [-] = ("<< xx <<") ";
}
