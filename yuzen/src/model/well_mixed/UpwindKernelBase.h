/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  UpwindKernelBase.h : the part shared by the physics terms that use the
*                       four-point upwind interpolation (soluto-capillarity
*                       and the body forces)
*
*  For each face, both of them have the same shape:
*    - they hold the four-point interpolation coefficients c1, and the
*      coefficients c2 that take effect with the sign of the upwind side
*    - they form the difference dc and the mean cc of the two cells that
*      straddle the face
*    - they add into every point of the stencil
*  The order of those additions is the tricky part, so the rule is written
*  down here, in one place only.
*
*  ** The per-row ordering of the directions in the Jacobian **
*  The order of the directions differs from row to row. The actual ordering
*  is
*
*    row = the diagonal : all directions in "stencil order"
*                         (Ym, Xm, Xp, Yp in 2-D)
*    row = any other    : the direction that has this row as its neighbor
*                         first, the rest in "axis order"
*
*  Example (1-D): row m is gLm*aL + gRm*aR, whereas row p is
*  gRp*aR + gLp*aL. A floating-point sum depends on the order it is taken
*  in, so do not disturb this rule.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _UPWIND_KERNEL_BASE_INCLUDED_
#define _UPWIND_KERNEL_BASE_INCLUDED_

#include "PhysicsTermBase.h"
#include "FiniteDifferences.h"
#include "NonDimensionalNumbers.h"
#include "LubFields.h"
#include "MyMath.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

template<class T>
class UpwindKernelBase : public PhysicsTermBase<T>
{
public:
	typedef typename T::Grid      Grid;
	typedef typename T::Site      Site;
	typedef typename T::Field     Field;
	typedef typename T::tEdge     tEdge;
	typedef typename T::tAdjacent tAdjacent;
	typedef typename T::Boundary  Boundary;

	//  The part shared by the physics terms that use upwinding (Marangoni
	//  and the in-plane body forces). These enter in the form
	//  "surface stress x film height" rather than as a gradient of the
	//  mobility, so, just as for advection, the difference has to be
	//  switched according to the direction of the flow.
	//  The interpolation scheme is held in _scheme, and its coefficients are
	//  laid out per direction in _c2.
	/** @brief Leave the interpolation scheme unset; a derived term selects it. */
	UpwindKernelBase() : PhysicsTermBase<T>() {}

	/**
	 * @brief Select the face interpolation (Numerics.Convection.Height) and
	 *        rebuild the coefficients. Called after Build().
	 */
	void SetConvectionScheme(C_Interpolation_base4::t_Scheme t) { _scheme.Build(t); this->SetLinearCoef(); }
	/** @brief The face interpolation in use. */
	const C_Interpolation_base4& Scheme() const { return _scheme; }

protected:
	//--------------------------------------------------------------
	//  A limited scheme (Koren, van Leer) has no fixed coefficients: the
	//  face value is formed from the four film heights of the face each
	//  time, A holds it as weights frozen at the current film (their sum
	//  times h is the flux itself, so A(h) h is exact), and the Jacobian
	//  adds the rest of the exact derivative of the flux. Each term keeps its
	//  fixed-coefficient code, untouched, for the other schemes.
	//--------------------------------------------------------------
	bool Limited() const { return _scheme.IsLimited(); }

	//--------------------------------------------------------------
	//  The face of a wall (Symmetry, Pinning) carries no flux of these
	//  terms. Their flux is not a gradient: it is a face value of the film
	//  (four-point, upwinded, reaching into the ghosts) times a driving
	//  force, and the ghosts of a wall do not make that face value zero (a
	//  Symmetry ghost mirrors the film, so h_f is the film at the wall), and
	//  the in-plane body forces would take liquid through a closed wall that
	//  way (tools/proptest.py). Each term returns zero for the
	//  mobility of such a face, in A, in the Jacobian and in the velocity
	//  alike, so the row of the cell next to the wall has no wall flux and
	//  no entry in the ghost columns from it.
	//--------------------------------------------------------------
	bool WallFace(Site &s, tEdge te) const
	{ return T::WallFace(this->_G, this->_boundary, this->_Nx, this->_Ny, s, te); }

	/** @brief The four film heights (LL, L, R, RR) of the face in direction d. */
	void FaceHeights(const std::vector<Site> &svec, int d, const Field &field, double f[4]) const
	{
		for(int j = 0; j < 4; ++j)
			f[j] = field.HH(const_cast<Site&>(svec[T::FaceStencil(d, j)]));
	}

	C_Interpolation_base4 _scheme;		// the interpolation scheme
	Vector<double> _c2[4];				// coefficients that take effect with the upwind sign (per direction)

	/** @brief Read and write access to _c2, one entry per direction. */
	Vector<double>&       C2(int d)       { return _c2[d]; }
	const Vector<double>& C2(int d) const { return _c2[d]; }

	//--------------------------------------------------------------
	//  Build the four-point interpolation coefficients for every
	//  direction.
	//  With perSpacing = true they are divided by the mesh spacing, which
	//  is the convention used on the body force side.
	//--------------------------------------------------------------
	void SetUpwindCoef(bool perSpacing)
	{
		for(int d = 0; d < T::nDir; ++d)
		{
			this->C(d) = 0;
			_c2[d].resize(this->_n_stencil);
			_c2[d] = 0;
		}

		static const C_Interpolation_base4::t_iShift pos[4] = {
			C_Interpolation_base4::LL, C_Interpolation_base4::L,
			C_Interpolation_base4::R,  C_Interpolation_base4::RR };

		for(int d = 0; d < T::nDir; ++d)
		{
			const double h = perSpacing ? T::Spacing(this->_G, d/2) : 1.0;
			for(int j = 0; j < 4; ++j)
			{
				const int idx = T::FaceStencil(d, j);
				this->C(d)[idx] = _scheme.c1(pos[j]) / h;
				_c2[d][idx]     = _scheme.c2(pos[j]) / h;
			}
		}
	}

	//--------------------------------------------------------------
	//  From the two cells that straddle the face in direction d, form the
	//  difference dc and the mean cc.
	//  The face lies between the two inner points (L and R) of the
	//  four-point interpolation.
	//--------------------------------------------------------------
	void FaceConc(const std::vector<Site> &svec, int d, const Field &field,
	              double &dc, double &cc) const
	{
		const int lo = T::FaceStencil(d, 1);
		const int hi = T::FaceStencil(d, 2);
		const double h = T::Spacing(this->_G, d/2);
		dc = (field.CC(const_cast<Site&>(svec[hi])) - field.CC(const_cast<Site&>(svec[lo]))) / h;
		cc = (field.CC(const_cast<Site&>(svec[hi])) + field.CC(const_cast<Site&>(svec[lo]))) / 2.0;
	}

	/**
	 * @brief Whether the stencil index idx is one of the four points of
	 *        direction d.
	 */
	static bool InFace(int d, int idx)
	{
		for(int j = 0; j < 4; ++j) if(T::FaceStencil(d, j) == idx) return true;
		return false;
	}

	//--------------------------------------------------------------
	//  Add into the Jacobian, row by row. The ordering follows the rule
	//  stated at the top of this file.
	//    gv(d, idx) : the Jacobian of direction d at stencil point idx
	//    a[d]       : the weight of direction d
	//--------------------------------------------------------------
	template<class GetG>
	void ScatterJacobiWide(CRS_Matrix &crsJ, int irow, Site &st,
	                       GetG gv, const double *a, double Coef)
	{
		for(int k = 0; k < T::nWide; ++k)
		{
			const int idx = T::WideIdx(k);
			double val = 0;

			if(T::IsDiag(idx))
			{
				// The diagonal runs in stencil order
				for(int j = 0; j < T::nDir; ++j)
				{
					const int d = T::StencilDir(j);
					val += gv(d, idx) * a[d];
				}
				crsJ.Diag(irow) += Coef*val;
			}
			else
			{
				// The direction that has this row as its neighbor comes
				// first; the rest follow in axis order
				const int first = T::DirOfAdj(idx);
				if(first >= 0 && InFace(first, idx)) val += gv(first, idx) * a[first];
				for(int d = 0; d < T::nDir; ++d)
				{
					if(d == first) continue;
					if(!InFace(d, idx)) continue;
					val += gv(d, idx) * a[d];
				}
				crsJ(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx))) += Coef*val;
			}
		}
	}
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_UPWIND_KERNEL_BASE_INCLUDED_
