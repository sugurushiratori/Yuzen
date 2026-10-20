/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  HeightEquation_1D.cpp : the common part has moved to the template in HeightEquationBase.h
*
*  All that remains here is the boundary condition of the velocity field
*  (the expression itself differs with the dimension).
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "HeightEquation_1D.h"
#include "ConcEquation_1D.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

using namespace yuzen::util;
using namespace Global;

//  Fill the ghost layers of the velocity field according to the type of the
//  boundary. Symmetry and pinning give 0, periodic copies the ghost faces
//  from just inside the opposite end (each seam face keeps its own), and
//  outlet and inlet extrapolate with zero gradient.
void WellMixed_HeightEquation_1D::SetVelocityBC(VectorField1D<> &vv)
{
	vv.x(G->At<Node>(0)) = 0;

	switch(_boundary->TypeXmin())
	{
	default:
	case BoundaryCondition::Symmetry:
	case BoundaryCondition::Pinning:
		vv.x(G->At<Node>(1))	= 0;
		vv.x(G->At<Node>(2))	= 0;
		break;
	case BoundaryCondition::Periodic:
		vv.x(G->At<Node>(1))	= vv.x(G->At<Node>(_Nx-3));
		break;
	case BoundaryCondition::Outlet:
		vv.x(G->At<Node>(1))	= vv.x(G->At<Node>(2));
		break;
	case BoundaryCondition::FreshInlet:
		vv.x(G->At<Node>(1))	= vv.x(G->At<Node>(2));
		vv.x(G->At<Node>(0))	= vv.x(G->At<Node>(1));
		break;
	}

	switch(_boundary->TypeXmax())
	{
	default:
	case BoundaryCondition::Symmetry:
	case BoundaryCondition::Pinning:
		vv.x(G->At<Node>(_Nx-2))	= 0.0;
		vv.x(G->At<Node>(_Nx-1))	= 0.0;
		break;
	case BoundaryCondition::Periodic:
		//  Node _Nx-2, the seam face as seen from the last cell, keeps the
		//  velocity computed on it -- the flux the film-height equation
		//  takes through it from that side -- and only the ghost face beyond
		//  comes from the other end, not the velocity of Node 2, the same
		//  face seen from the first cell. For every periodic force the two
		//  are the same; for the centrifugal one, whose arm x - x0 jumps by
		//  the length of the domain across the seam, they are not, and the
		//  last cell's solute would be carried by a velocity its film does
		//  not have: a first-order error in the concentration (in the
		//  manufactured solution wm_all_1d the order of c falls to 0.9 on
		//  finer grids).
		vv.x(G->At<Node>(_Nx-1))	= vv.x(G->At<Node>(3));
		break;
	case BoundaryCondition::FreshInlet:
		//  Node _Nx-2 is the inlet face itself (Node 2 at Xmin), and keeps
		//  the velocity the film-height equation gave it; only the ghost face
		//  beyond is extrapolated. (Overwritten with the velocity of the face
		//  one cell inside, the solute carried in through an inlet at Xmax
		//  would not be what the film carries.)
		vv.x(G->At<Node>(_Nx-1))	= vv.x(G->At<Node>(_Nx-2));
		break;
	case BoundaryCondition::Outlet:
		vv.x(G->At<Node>(_Nx-1))	= vv.x(G->At<Node>(_Nx-2));
		break;
	}

}
