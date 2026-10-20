/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LaplacePressure_1D.cpp : the common part has moved to LaplacePressureKernel.h
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "LaplacePressure_1D.h"
#include "NonDimensionalNumbers.h"
#include "LubFields.h"
#include "MyMath.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

//  The mobility is set to zero on a face that touches a pinned boundary
bool LaplacePressure_WM1D::PinnedFace(Xedge &s, tEdge tedge) const
{
	if(_boundary->TypeXmin() == BoundaryCondition::Pinning && _G->ix(s) ==   1   && tedge==PhysicsFace::Xp) return true;
	if(_boundary->TypeXmin() == BoundaryCondition::Pinning && _G->ix(s) ==   2   && tedge==PhysicsFace::Xm) return true;
	if(_boundary->TypeXmax() == BoundaryCondition::Pinning && _G->ix(s) == _Nx-3 && tedge==PhysicsFace::Xp) return true;
	if(_boundary->TypeXmax() == BoundaryCondition::Pinning && _G->ix(s) == _Nx-2 && tedge==PhysicsFace::Xm) return true;
	return false;
}

//  The finite-difference coefficients of the biharmonic operator. A
//  four-point 1/dx^4 stencil is laid down for each of the left and right
//  faces (C(0), C(1)). CV(0) is the same with the sign flipped back and dx
//  multiplied in again, for the velocity computation. Unlike the 2-D case,
//  no cross term appears.
void LaplacePressure_WM1D::SetLinearCoef()
{
	double dx = _G->dx();
	double dx2 = dx*dx;
	double dx4 = dx2*dx2;

	typedef C_Grid1D gt;

	C(0)[gt::M] = +1.0/dx4;
	C(0)[gt::m] = 1.0*(-3.0/dx4);
	C(0)[gt::o] = 1.0*(+3.0/dx4);
	C(0)[gt::p] = -1.0/dx4;

	C(1)[gt::m] = -1.0/dx4;
	C(1)[gt::o] = 1.0*(+3.0/dx4);
	C(1)[gt::p] = 1.0*(-3.0/dx4);
	C(1)[gt::P] = +1.0/dx4;


	for(int s = 0; s < _n_stencil; s++)
	{
		CV(0)[s] = -C(0)[s] * dx;
	}

}

//  The third difference across face d (0: Xm, 1: Xp) of the values v on
//  the stencil of a row, (v_{b+1} - 3 v_b + 3 v_a - v_{a-1}) / dx^4 with a
//  and b the cells below and above the face, as differences of differences:
//  the cell on either side of the face computes it from the same four
//  values in the same order (LaplacePressureKernel::AddConservative), and
//  with mirror ghosts (Symmetry) it is exactly 0 on an edge face.
double LaplacePressure_WM1D::FaceGradient(const double *v, int d) const
{
	typedef C_Grid1D gt;
	const double dx = _G->dx();
	const double dx2 = dx*dx;
	const double dx4 = dx2*dx2;
	double a1, a, b, b1;
	if(d == 0) { a1 = v[gt::M]; a = v[gt::m]; b = v[gt::o]; b1 = v[gt::p]; }
	else       { a1 = v[gt::m]; a = v[gt::o]; b = v[gt::p]; b1 = v[gt::P]; }
	return ( ((b1 - b) - 2.0*(b - a)) + (a - a1) ) / dx4;
}

//  The velocity produced by the capillary pressure. The third difference of
//  the free surface h+s is written out as it stands, multiplied by the face
//  mobility DL, and divided by the face film height hL.
//  The two outermost points are skipped, because their stencils reach
//  outside the grid; the boundary condition fills them in.
void LaplacePressure_WM1D::CalVelocity(Field_WM1D &field)
{
	double DL;
	double hL;
	double Qhx, Qsx;
	Xedge xe;
	typedef C_Grid1D gt;
	double dx, rdx3;
	dx = _G->dx();
	rdx3 = 1.0/(dx*dx*dx);

	double hp, ho, hm, hM;
	double sp, so, sm, sM;

	int ix;
	//  Per-thread workspace. It is created outside the parallel region and
	//  handed out with firstprivate (declaring it inside would take the
	//  whole scope with it in a build without OpenMP).
	std::vector<Xedge> svec(_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(ix,xe,DL,hL,Qhx, Qsx, hp, ho, hm, hM, sp, so, sm, sM) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = _G->imin(); ii <= _G->imax(); ii++)
		{
			xe = _G->Of<Xedge>(ii);
			ix = _G->ix(xe);
			if(ix < 2 || _Nx-2 < ix) continue;

			DL = NonlinearTerm(xe, PhysicsFace::Xm, field);
			hL = 0.5*(field.HH(_G->ShiftXm(xe)) +  field.HH(xe));

			_G->SetSite(svec, xe);

			ho = field.HH(svec[gt::o]);		so = field.SS(svec[gt::o]);
			hp = field.HH(svec[gt::p]);		sp = field.SS(svec[gt::p]);
			hm = field.HH(svec[gt::m]);		sm = field.SS(svec[gt::m]);
			hM = field.HH(svec[gt::M]);		sM = field.SS(svec[gt::M]);


			Qhx = ( (hp-2.0*ho+hm)  -(ho-2.0*hm+hM) ) * rdx3;
			Qsx = ( (sp-2.0*so+sm)  -(so-2.0*sm+sM) ) * rdx3;

			Qhx *= _Magnitude * DL;
			Qsx *= _Magnitude * DL;

			field.UU().U_LaplacePressure().x(_G->HalfXm(xe)) = (double)((Qhx+Qsx)/hL);
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}
