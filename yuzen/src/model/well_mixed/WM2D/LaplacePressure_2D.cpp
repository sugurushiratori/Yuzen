/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LaplacePressure_2D.cpp : the common part has moved to LaplacePressureKernel.h
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "LaplacePressure_2D.h"
#include "NonDimensionalNumbers.h"
#include "LubFields.h"
#include "MyMath.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

//  The mobility is set to zero on a face that touches a pinned boundary
bool LaplacePressure_WM2D::PinnedFace(Zface &s, tEdge tedge) const
{
	if(_boundary->TypeXmin() == BoundaryCondition::Pinning && _G->ix(s) ==   1   && tedge==PhysicsFace::Xp) return true;
	if(_boundary->TypeXmin() == BoundaryCondition::Pinning && _G->ix(s) ==   2   && tedge==PhysicsFace::Xm) return true;
	if(_boundary->TypeXmax() == BoundaryCondition::Pinning && _G->ix(s) == _Nx-3 && tedge==PhysicsFace::Xp) return true;
	if(_boundary->TypeXmax() == BoundaryCondition::Pinning && _G->ix(s) == _Nx-2 && tedge==PhysicsFace::Xm) return true;
	if(_boundary->TypeYmin() == BoundaryCondition::Pinning && _G->iy(s) ==   1   && tedge==PhysicsFace::Yp) return true;
	if(_boundary->TypeYmin() == BoundaryCondition::Pinning && _G->iy(s) ==   2   && tedge==PhysicsFace::Ym) return true;
	if(_boundary->TypeYmax() == BoundaryCondition::Pinning && _G->iy(s) == _Ny-3 && tedge==PhysicsFace::Yp) return true;
	if(_boundary->TypeYmax() == BoundaryCondition::Pinning && _G->iy(s) == _Ny-2 && tedge==PhysicsFace::Ym) return true;
	return false;
}

//  The finite-difference coefficients of the biharmonic operator. For each
//  of the four faces (C(0)..C(3)), the 1/dx^4 and 1/dy^4 terms are laid
//  down, together with the cross term 1/dx^2dy^2, which has no 1-D
//  counterpart. CV(0) and CV(1) are the same with the sign flipped back and
//  dx and dy multiplied in again, for the velocity computation.
void LaplacePressure_WM2D::SetLinearCoef()
{
	double dx = _G->dx();
	double dy = _G->dy();
	double dx2 = dx*dx;
	double dy2 = dy*dy;
	double dx4 = dx2*dx2;
	double dy4 = dy2*dy2;
	double dx2dy2 = dx2*dy2;

	typedef C_Grid2D gt;

	C(0)[gt::mm] = +1.0/dx2dy2;
	C(0)[gt::om] = -1.0/dx2dy2;
	C(0)[gt::Mo] = +1.0/dx4;
	C(0)[gt::mo] = 1.0*(-3.0/dx4 -2.0/dx2dy2);
	C(0)[gt::oo] = 1.0*(+3.0/dx4 +2.0/dx2dy2);
	C(0)[gt::po] = -1.0/dx4;
	C(0)[gt::mp] = +1.0/dx2dy2;
	C(0)[gt::op] = -1.0/dx2dy2;

	C(1)[gt::om] = -1.0/dx2dy2;
	C(1)[gt::pm] = +1.0/dx2dy2;
	C(1)[gt::mo] = -1.0/dx4;
	C(1)[gt::oo] = 1.0*(+3.0/dx4 +2.0/dx2dy2);
	C(1)[gt::po] = 1.0*(-3.0/dx4 -2.0/dx2dy2);
	C(1)[gt::Po] = +1.0/dx4;
	C(1)[gt::op] = -1.0/dx2dy2;
	C(1)[gt::pp] = +1.0/dx2dy2;

	C(2)[gt::oM] = +1.0/dy4;
	C(2)[gt::mm] = +1.0/dx2dy2;
	C(2)[gt::om] = 1.0*(-3.0/dy4 -2.0/dx2dy2);
	C(2)[gt::pm] = +1.0/dx2dy2;
	C(2)[gt::mo] = -1.0/dx2dy2;
	C(2)[gt::oo] = 1.0*(+3.0/dy4 +2.0/dx2dy2);
	C(2)[gt::po] = -1.0/dx2dy2;
	C(2)[gt::op] = -1.0/dy4;

	C(3)[gt::om] = -1.0/dy4;
	C(3)[gt::mo] = -1.0/dx2dy2;
	C(3)[gt::oo] = 1.0*(+3.0/dy4 +2.0/dx2dy2);
	C(3)[gt::po] = -1.0/dx2dy2;
	C(3)[gt::mp] = +1.0/dx2dy2;
	C(3)[gt::op] = 1.0*(-3.0/dy4 -2.0/dx2dy2);
	C(3)[gt::pp] = +1.0/dx2dy2;
	C(3)[gt::oP] = +1.0/dy4;


	for(int s = 0; s < _n_stencil; s++)
	{
		CV(0)[s] = -C(0)[s] * dx;
		CV(1)[s] = -C(2)[s] * dy;
	}

}

//  The third difference across face d (0: Xm, 1: Xp, 2: Ym, 3: Yp) of the
//  values v on the stencil of a row: along the normal,
//  (v_{b+1} - 3 v_b + 3 v_a - v_{a-1}) / dn^4, plus the cross term, the
//  difference across the face of the second difference along it, / dx^2 dy^2
//  (a and b the cells below and above the face). Both are written as
//  differences of differences, so the cell on either side of the face
//  computes the same value from the same values in the same order
//  (LaplacePressureKernel::AddConservative), and with mirror ghosts
//  (Symmetry) it is exactly 0 on an edge face.
double LaplacePressure_WM2D::FaceGradient(const double *v, int d) const
{
	typedef C_Grid2D gt;
	const double dx = _G->dx(), dy = _G->dy();
	const double dx2 = dx*dx, dy2 = dy*dy;
	const double dn4 = (d < 2) ? dx2*dx2 : dy2*dy2;
	const double dx2dy2 = dx2*dy2;
	//  normal: a1 a | b b1;  cross: (am a ap) and (bm b bp) along the face
	double a1, a, b, b1, am, ap, bm, bp;
	switch(d)
	{
	case 0:  a1 = v[gt::Mo]; a = v[gt::mo]; b = v[gt::oo]; b1 = v[gt::po];
	         am = v[gt::mm]; ap = v[gt::mp]; bm = v[gt::om]; bp = v[gt::op]; break;
	case 1:  a1 = v[gt::mo]; a = v[gt::oo]; b = v[gt::po]; b1 = v[gt::Po];
	         am = v[gt::om]; ap = v[gt::op]; bm = v[gt::pm]; bp = v[gt::pp]; break;
	case 2:  a1 = v[gt::oM]; a = v[gt::om]; b = v[gt::oo]; b1 = v[gt::op];
	         am = v[gt::mm]; ap = v[gt::pm]; bm = v[gt::mo]; bp = v[gt::po]; break;
	default: a1 = v[gt::om]; a = v[gt::oo]; b = v[gt::op]; b1 = v[gt::oP];
	         am = v[gt::mo]; ap = v[gt::po]; bm = v[gt::mp]; bp = v[gt::pp]; break;
	}
	const double normal = ( ((b1 - b) - 2.0*(b - a)) + (a - a1) ) / dn4;
	const double cross  = ( ((bp - b) - (b - bm)) - ((ap - a) - (a - am)) ) / dx2dy2;
	return normal + cross;
}

//  The velocity produced by the capillary pressure. The third differences of
//  the free surface h+s are written out in x and in y, multiplied by the
//  face mobilities (DL, DB), and divided by the face film heights (hL, hB).
//  ratio2x and ratio2y are the mesh ratios; they give the weights of the
//  cross terms.
//  The two outermost points are skipped, because their stencils reach
//  outside the grid; the boundary condition fills them in.
void LaplacePressure_WM2D::CalVelocity(Field_WM2D &field)
{
	double DL, DB;
	double hL, hB;
	double Qhx, Qhy, Qsx, Qsy;
	Zface zf;
	typedef C_Grid2D gt;
	double dx, dy, rdx3, rdy3, ratio2x, ratio2y;
	dx = _G->dx();
	dy = _G->dy();
	rdx3 = 1.0/(dx*dx*dx);
	rdy3 = 1.0/(dy*dy*dy);
	ratio2x = (dx*dx)/(dy*dy);
	ratio2y = (dy*dy)/(dx*dx);

	double hpo, hoo, hmo, hMo, hop, hom, hpm, hoM, hmp, hmm;
	double spo, soo, smo, sMo, sop, som, spm, soM, smp, smm;

	int ix, jy;
	//  Per-thread workspace. It is created outside the parallel region and
	//  handed out with firstprivate (declaring it inside would take the
	//  whole scope with it in a build without OpenMP).
	std::vector<Zface> svec(_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(ix,jy,zf,DL,DB,hL,hB,Qhx,Qhy, Qsx, Qsy, \
				hpo, hoo, hmo, hMo, hop, hom, hpm, hoM, hmp, hmm, \
				spo, soo, smo, sMo, sop, som, spm, soM, smp, smm) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = _G->imin(); ii <= _G->imax(); ii++)
		{
			zf = _G->Of<Zface>(ii);
			ix = _G->ix(zf); jy = _G->iy(zf);
			if(ix < 2 || _Nx-2 < ix) continue;
			if(_Ny != 1 && (jy < 2 || _Ny-2 < jy)) continue;

			DL = NonlinearTerm(zf, PhysicsFace::Xm, field);
			DB = NonlinearTerm(zf, PhysicsFace::Ym, field);
			hL = 0.5*(field.HH(_G->ShiftXm(zf)) +  field.HH(zf));
			hB = 0.5*(field.HH(_G->ShiftYm(zf)) +  field.HH(zf));

			_G->SetSite(svec, zf);

			hop = field.HH(svec[gt::op]);		sop = field.SS(svec[gt::op]);
			hoo = field.HH(svec[gt::oo]);		soo = field.SS(svec[gt::oo]);
			hom = field.HH(svec[gt::om]);		som = field.SS(svec[gt::om]);
			hoM = field.HH(svec[gt::oM]);		soM = field.SS(svec[gt::oM]);
			hpo = field.HH(svec[gt::po]);		spo = field.SS(svec[gt::po]);
			hpm = field.HH(svec[gt::pm]);		spm = field.SS(svec[gt::pm]);
			hmo = field.HH(svec[gt::mo]);		smo = field.SS(svec[gt::mo]);
			hMo = field.HH(svec[gt::Mo]);		sMo = field.SS(svec[gt::Mo]);
			hmp = field.HH(svec[gt::mp]);		smp = field.SS(svec[gt::mp]);
			hmm = field.HH(svec[gt::mm]);		smm = field.SS(svec[gt::mm]);


			Qhx = ( ( (hpo-2.0*hoo+hmo) +(hop-2.0*hoo+hom)*ratio2x ) -( (hoo-2.0*hmo+hMo) +(hmp-2.0*hmo+hmm)*ratio2x ) ) * rdx3;
			Qhy = ( ( (hpo-2.0*hoo+hmo)*ratio2y +(hop-2.0*hoo+hom) ) -( (hpm-2.0*hom+hmm)*ratio2y +(hoo-2.0*hom+hoM) ) ) * rdy3;
			Qsx = ( ( (spo-2.0*soo+smo) +(sop-2.0*soo+som)*ratio2x ) -( (soo-2.0*smo+sMo) +(smp-2.0*smo+smm)*ratio2x ) ) * rdx3;
			Qsy = ( ( (spo-2.0*soo+smo)*ratio2y +(sop-2.0*soo+som) ) -( (spm-2.0*som+smm)*ratio2y +(soo-2.0*som+soM) ) ) * rdy3;

			Qhx *= _Magnitude * DL;
			Qhy *= _Magnitude * DB;
			Qsx *= _Magnitude * DL;
			Qsy *= _Magnitude * DB;

			if(_Ny == 1) {Qhy = Qsy = 0;}	// ##### TOTO provisional handling of the pseudo 1-D case

			field.UU().U_LaplacePressure().x(_G->HalfXm(zf)) = (double)((Qhx+Qsx)/hL);
			field.UU().U_LaplacePressure().y(_G->HalfYm(zf)) = (double)((Qhy+Qsy)/hB);
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}
