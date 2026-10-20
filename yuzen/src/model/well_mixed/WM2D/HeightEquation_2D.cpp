/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  HeightEquation_2D.cpp : the common part has moved to the template in HeightEquationBase.h
*
*  All that remains here is the boundary condition of the velocity field
*  (the expression itself differs with the dimension).
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "Opt.h"
#include "YuzenUtil.h"
#include "HeightEquation_2D.h"
#include "ConcEquation_2D.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

using namespace yuzen::util;
using namespace Global;

//  Fill the ghost layers of the velocity field according to the type of the
//  boundary. The two ends in x and the two ends in y are walked separately.
//  Symmetry and pinning set the component normal to the boundary to 0 and
//  copy the tangential component from just inside. Periodic copies the
//  ghost faces from just inside the opposite end (the seam faces keep their
//  own), and outlet and inlet extrapolate with zero
//  gradient.
void WellMixed_HeightEquation_2D::SetVelocityBC(VectorField2D<> &vv)
{

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_Ny))
#endif
	for(int jy = 1; jy < _Ny-1 ; jy++)
	{
		vv.x(G->At<Yedge>(0, jy)) = 0;

		switch(_boundary->TypeXmin())
		{
		default:
		case BoundaryCondition::Symmetry:
		case BoundaryCondition::Pinning:
			vv.y(G->At<Xedge>(1, jy))	= vv.y(G->At<Xedge>(2, jy));
			vv.x(G->At<Yedge>(1, jy))	= 0;
			vv.x(G->At<Yedge>(2, jy))	= 0;
			break;
		case BoundaryCondition::Periodic:
			vv.y(G->At<Xedge>(1, jy))	= vv.y(G->At<Xedge>(_Nx-3, jy));
			vv.x(G->At<Yedge>(1, jy))	= vv.x(G->At<Yedge>(_Nx-3, jy));
			break;
		case BoundaryCondition::Outlet:
			vv.y(G->At<Xedge>(1, jy))	= vv.y(G->At<Xedge>(2, jy));
			vv.x(G->At<Yedge>(1, jy))	= vv.x(G->At<Yedge>(2, jy));
			break;
		case BoundaryCondition::FreshInlet:
			vv.y(G->At<Xedge>(1, jy))	= vv.y(G->At<Xedge>(2, jy));
			vv.x(G->At<Yedge>(1, jy))	= vv.x(G->At<Yedge>(2, jy));
			vv.x(G->At<Yedge>(0, jy))	= vv.x(G->At<Yedge>(1, jy));
			break;
		}

		switch(_boundary->TypeXmax())
		{
		default:
		case BoundaryCondition::Symmetry:
		case BoundaryCondition::Pinning:
			vv.y(G->At<Xedge>(_Nx-2, jy))	= vv.y(G->At<Xedge>(_Nx-3, jy));
			vv.x(G->At<Yedge>(_Nx-2, jy))	= 0.0;
			vv.x(G->At<Yedge>(_Nx-1, jy))	= 0.0;
			break;
		case BoundaryCondition::Periodic:
			//  Yedge _Nx-2, the seam face seen from the last column, keeps its
			//  own velocity (the flux the film-height equation uses there) and
			//  only the ghost face beyond comes from the other end; see the
			//  1-D SetVelocityBC. (Not the velocity of Yedge 2, which
			//  differs from it for the centrifugal force.)
			vv.y(G->At<Xedge>(_Nx-2, jy))	= vv.y(G->At<Xedge>(2, jy));
			vv.x(G->At<Yedge>(_Nx-1, jy))	= vv.x(G->At<Yedge>(3, jy));
			break;
		case BoundaryCondition::FreshInlet:
			//  Yedge _Nx-2 is the inlet face and keeps its velocity, as Yedge 2
			//  does at Xmin (see the 1-D SetVelocityBC).
			vv.y(G->At<Xedge>(_Nx-2, jy))	= vv.y(G->At<Xedge>(_Nx-3, jy));
			vv.x(G->At<Yedge>(_Nx-1, jy))	= vv.x(G->At<Yedge>(_Nx-2, jy));
			break;
		case BoundaryCondition::Outlet:
			vv.y(G->At<Xedge>(_Nx-2, jy))	= vv.y(G->At<Xedge>(_Nx-3, jy));
			vv.x(G->At<Yedge>(_Nx-1, jy))	= vv.x(G->At<Yedge>(_Nx-2, jy));
			break;
		}
	}

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par(_Nx))
#endif
	for(int ix = 1; ix < _Nx-1 ; ix++)
	{
		vv.y(G->At<Xedge>(ix, 0)) = 0;

		switch(_boundary->TypeYmin())
		{
		default:
		case BoundaryCondition::Symmetry:
		case BoundaryCondition::Pinning:
			vv.x(G->At<Yedge>(ix, 1))	= vv.x(G->At<Yedge>(ix, 2));
			vv.y(G->At<Xedge>(ix, 1))	= 0;
			vv.y(G->At<Xedge>(ix, 2))	= 0;
			break;
		case BoundaryCondition::Periodic:
			vv.x(G->At<Yedge>(ix, 1))	= vv.x(G->At<Yedge>(ix, _Ny-3));
			vv.y(G->At<Xedge>(ix, 1))	= vv.y(G->At<Xedge>(ix, _Ny-3));
			break;
		case BoundaryCondition::Outlet:
			vv.x(G->At<Yedge>(ix, 1))	= vv.x(G->At<Yedge>(ix, 2));
			vv.y(G->At<Xedge>(ix, 1))	= vv.y(G->At<Xedge>(ix, 2));
			break;
		case BoundaryCondition::FreshInlet:
			vv.x(G->At<Yedge>(ix, 1))	= vv.x(G->At<Yedge>(ix, 2));
			vv.y(G->At<Xedge>(ix, 1))	= vv.y(G->At<Xedge>(ix, 2));
			vv.y(G->At<Xedge>(ix, 0))	= vv.y(G->At<Xedge>(ix, 1));
			break;
		}

		switch(_boundary->TypeYmax())
		{
		default:
		case BoundaryCondition::Symmetry:
		case BoundaryCondition::Pinning:
			vv.x(G->At<Yedge>(ix, _Ny-2))	= vv.x(G->At<Yedge>(ix, _Ny-3));
			vv.y(G->At<Xedge>(ix, _Ny-2))	= 0;
			vv.y(G->At<Xedge>(ix, _Ny-1))	= 0;
			break;
		case BoundaryCondition::Periodic:
			//  Xedge _Ny-2 is the seam face and keeps its own velocity, as in x.
			vv.x(G->At<Yedge>(ix, _Ny-2))	= vv.x(G->At<Yedge>(ix, 2));
			vv.y(G->At<Xedge>(ix, _Ny-1))	= vv.y(G->At<Xedge>(ix, 3));
			break;
		case BoundaryCondition::FreshInlet:
			//  Xedge _Ny-2 is the inlet face and keeps its velocity (Xedge 2 at Ymin).
			vv.x(G->At<Yedge>(ix, _Ny-2))	= vv.x(G->At<Yedge>(ix, _Ny-3));
			vv.y(G->At<Xedge>(ix, _Ny-1))	= vv.y(G->At<Xedge>(ix, _Ny-2));
			break;
		case BoundaryCondition::Outlet:
			vv.x(G->At<Yedge>(ix, _Ny-2))	= vv.x(G->At<Yedge>(ix, _Ny-3));
			vv.y(G->At<Xedge>(ix, _Ny-1))	= vv.y(G->At<Xedge>(ix, _Ny-2));
			break;
		}

	}

}
