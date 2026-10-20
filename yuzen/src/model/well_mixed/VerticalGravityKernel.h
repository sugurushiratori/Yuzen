/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  VerticalGravityKernel.h : the dimension-independent implementation of
*                            the vertical gravity term
*
*  The only difference between 1-D and 2-D is that "2-D adds the Y
*  contribution"; the expression itself is the same. The directions and the
*  stencil are therefore kept in the tables of the traits
*  (GridTraits.h), and the body is kept as a single copy.
*
*  ** On the order of summation - change it and the results stop agreeing
*  bit for bit **
*  There are two orderings of the directions, and which one is used depends
*  on the place.
*    - the summation of G in each row of SetMatrixA
*                                   : axis order (Xm, Xp, Ym, Yp)    -> T::Dir(d)
*    - the diagonal of SetJacobiMatrix, and the summations of Too / bSum
*                                   : stencil order (Ym, Xm, Xp, Yp) -> T::StencilDir(k)
*  A floating-point sum depends on the order it is taken in. Keep this
*  distinction exactly as it is.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _VERTICAL_GRAVITY_KERNEL_INCLUDED_
#define _VERTICAL_GRAVITY_KERNEL_INCLUDED_

#include "PhysicsTermBase.h"
#include "NonDimensionalNumbers.h"
#include "LubFields.h"
#include "MyMath.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

template<class T>
class VerticalGravityKernel : public PhysicsTermBase<T>
{
public:
	typedef typename T::Grid      Grid;
	typedef typename T::Site      Site;
	typedef typename T::Field     Field;
	typedef typename T::tEdge     tEdge;
	typedef typename T::tAdjacent tAdjacent;
	typedef typename T::Boundary  Boundary;

	VerticalGravityKernel() : PhysicsTermBase<T>() {}

	//  Vertical gravity. The coefficient is -Bo/3.
	//  The coefficient _Magnitude (sign included) and the
	//  finite-difference coefficients are decided here.
	//  For the contract, see PhysicsTermBase.h.
	/**
	 * @brief Switch the vertical gravity term on or off and build its
	 *        finite-difference coefficients.
	 * @param inp_phys  the Physics block of the input file
	 * @param g         grid the term is discretized on
	 * @param nonDim    the set of non-dimensional numbers; the Bond number is
	 *                  taken from here
	 * @param bc        boundary condition, used to fold the ghost points into the stencil
	 *
	 * The overall coefficient _Magnitude is -Bo/3, the sign included.
	 */
	virtual void Build(Input::I_Physics &inp_phys, Grid *g, NonDimensionalNumbers *nonDim, Boundary *bc)
	{
		this->_G = g;
		this->_NonDim = nonDim;
		this->_boundary = bc;

		this->_enabled = inp_phys.VerticalGravity.Enable;
		this->_Magnitude = -1.0 * this->_NonDim->Bo() / 3.0;

		T::SetExtent(this->_G, this->_Nx, this->_Ny);
		SetLinearCoef();
	}

	/**
	 * @brief Add the Bond-number term, together with the contribution of the
	 *        substrate topography, to the matrix A.
	 */
	virtual void SetMatrixA(CRS_Matrix &crsA, const Field &field);
	/**
	 * @brief Add the same term to the Jacobian.
	 *
	 * The derivative of the face mobility times the difference of h + s,
	 * less the topography share A already holds: the exact Jacobian
	 * (tests/test_jacobian.cpp).
	 */
	virtual void SetJacobiMatrix(CRS_Matrix &crsJ, const Field &field);

	/**
	 * @brief A h of this term evaluated face by face (see
	 *        PhysicsTermBase::ConservativeForm): the flux through a face
	 *        between cells a (below) and b (above) is
	 *        Coef (G_f (h_b - h_a) + (s_b - s_a) K (h_a + h_b)) / dn^2,
	 *        G_f the face mobility (Mobility::Face)/mu and K = G_f/(h_a + h_b),
	 *        mu at the concentration of the face -- what A holds, formed
	 *        once from the values about the face so that the two cells
	 *        that share it take the same flux.
	 */
	virtual bool ConservativeForm() const { return this->_enabled; }
	virtual void AddConservative(Vector<double> &y, const Vector<double> &x,
	                             const CRS_Matrix &crs, const Field &field) const;
	/** @brief The flux of AddConservative through face d of s, in the state field. */
	virtual double FaceFlux(Site &s, int d, const Field &field, const std::vector<Site> &svec);
	/**
	 * @brief Write the gravity-driven velocity into field.UU().U_VerticalGravity().
	 */
	virtual void CalVelocity(Field &field);

protected:
	/**
	 * @brief The three-point second-derivative coefficients, which are
	 *        fixed by the mesh spacing alone, together with the
	 *        coefficients used to compute the velocity.
	 */
	virtual void SetLinearCoef();
	/** @brief The mobility h³/μ(c) on a face. */
	virtual double NonlinearTerm(Site&, tEdge, const Field &field);
	/**
	 * @brief The above, differentiated with respect to one point (for the
	 *        Jacobian).
	 */
	virtual double NonlinearJacobi(Site&, tEdge, tAdjacent, const Field &field);

	/**
	 * @brief The share K = G_f / (h_o + h_a) / mu of each of the two cells
	 *        of face tedge in the topography flux G_f Delta s (see the
	 *        implementation).
	 */
	double TopographyShare(Site&, tEdge, const Field &field) const;
};


//==================================================================
//  Implementation
//==================================================================

template<class T>
void VerticalGravityKernel<T>::SetLinearCoef()
{
	for(int d = 0; d < T::nDir; ++d) this->C(d) = 0;

	for(int d = 0; d < T::nDir; ++d)
	{
		const double h  = T::Spacing(this->_G, d / 2);
		const double h2 = h * h;
		this->C(d)[T::AdjIdx(d)]   = +1.0 / h2;
		this->C(d)[T::CenterIdx()] = -1.0 / h2;
	}

	for(int a = 0; a < T::nAxis; ++a)
	{
		const double h = T::Spacing(this->_G, a);
		for(int s = 0; s < this->_n_stencil; s++)
			this->CV(a)[s] = -this->C(2*a)[s] * h;
	}
}


//------------------------------------------------------------------
//  h^3 / mu on a face: the positive face mobility of the capillary term
//  (Mobility::Face), as the disjoining pressure takes it, so that the flux
//  of a face is one mobility times the difference of the potential
//  (DisjoiningPressureKernel.h, EnergyBudget.h). 0 against the ghost -h of a
//  Pinning edge, which the mean (h_o^3 + h_a^3)/2 this term once took gave
//  by symmetry.
//------------------------------------------------------------------
template<class T>
double VerticalGravityKernel<T>::NonlinearTerm(Site &s, tEdge tedge, const Field &field)
{
	double Hoo = field.HH(s);
	double Coo = field.CC(s);
	double Haj = 0, Caj = 0;

	if(T::IsFace(tedge))
	{
		Site a = T::Adj(this->_G, s, tedge);
		Haj = field.HH(a);
		Caj = field.CC(a);
	}
	double Cfc = 0.5*(Coo + Caj);

	double val = this->_NonDim->Mob().Face(Hoo, Haj);
	double Visco = this->_NonDim->Viscosity(Cfc);
	return val / Visco;
}


//------------------------------------------------------------------
//  The above, differentiated with respect to one point: the cell itself
//  (tc = center) or the neighbor across the face.
//------------------------------------------------------------------
template<class T>
double VerticalGravityKernel<T>::NonlinearJacobi(Site &s, tEdge tedge, tAdjacent tc, const Field &field)
{
	double Coo = field.CC(s);
	Site a = T::IsFace(tedge) ? T::Adj(this->_G, s, tedge) : s;
	const double ho = field.HH(s), ha = field.HH(a);
	double Caj = field.CC(a);
	double Cfc = 0.5*(Coo + Caj);

	double Visco = this->_NonDim->Viscosity(Cfc);
	const Mobility &mob = this->_NonDim->Mob();
	return (((int)tc == T::CenterIdx()) ? mob.Face_d1(ho, ha) : mob.Face_d1(ha, ho)) / Visco;
}


//------------------------------------------------------------------
//  The substrate topography. The flux of a face is G_f (Delta h + Delta s),
//  G_f the face mobility; A takes G_f Delta s as K (h_lo + h_hi) Delta s,
//  K = G_f / (h_lo + h_hi), the same share for the two cells of the face
//  (with K, A h is that flux exactly). 0 where G_f is.
//------------------------------------------------------------------
template<class T>
double VerticalGravityKernel<T>::TopographyShare(Site &s, tEdge tedge, const Field &field) const
{
	const Site a = T::Adj(this->_G, s, tedge);
	const double ho = field.HH(s), ha = field.HH(a);
	const double g = this->_NonDim->Mob().Face(ho, ha);
	const double Visco = this->_NonDim->Viscosity(0.5 * (field.CC(s) + field.CC(a)));
	return (g == 0) ? 0.0 : g / (ho + ha) / Visco;
}


//------------------------------------------------------------------
//  A h face by face (see the declaration)
//------------------------------------------------------------------
template<class T>
void VerticalGravityKernel<T>::AddConservative(Vector<double> &y, const Vector<double> &x,
                                               const CRS_Matrix &crs, const Field &field) const
{
	if(!this->_enabled) return;
	const double Coef = this->_Magnitude;
	const Mobility &mob = this->_NonDim->Mob();
	Site st, lo, hi;
	int irow, d;
	double r, f, dn, hl, hh, mu, gf, tk;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(st,lo,hi,irow,d,r,f,dn,hl,hh,mu,gf,tk)
#endif
	for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
	{
		st = this->_G->template Of<Site>(ii);
		irow = T::Index(this->_G, st);
		if( crs.IsBC(irow) ) continue;
		r = 0;
		for(d = 0; d < T::nDir; ++d)
		{
			const Site a = T::Adj(this->_G, st, T::Dir(d));
			if(d % 2 == 0) { lo = a;  hi = st; }
			else           { lo = st; hi = a;  }
			dn = T::Spacing(this->_G, d / 2);
			hl = field.HH(lo);	hh = field.HH(hi);
			mu = this->_NonDim->Viscosity(0.5 * (field.CC(lo) + field.CC(hi)));
			gf = mob.Face(hl, hh);
			tk = (gf == 0) ? 0.0 : gf / (hl + hh) / mu;
			gf /= mu;
			const double xl = x[T::Index(this->_G, lo)], xh = x[T::Index(this->_G, hi)];
			f = Coef * ( (gf * (xh - xl) + (field.SS(hi) - field.SS(lo)) * (tk * (xl + xh))) / (dn * dn) );
			r += (d % 2 == 0) ? -f : f;
		}
		y[irow] += r;
	}
}


template<class T>
double VerticalGravityKernel<T>::FaceFlux(Site &st, int d, const Field &field, const std::vector<Site> &svec)
{
	(void)svec;
	if(!this->_enabled) return 0;
	const Site a = T::Adj(this->_G, st, T::Dir(d));
	const Site lo = (d % 2 == 0) ? a : st, hi = (d % 2 == 0) ? st : a;
	const double dn = T::Spacing(this->_G, d / 2);
	const double hl = field.HH(lo), hh = field.HH(hi);
	const double mu = this->_NonDim->Viscosity(0.5 * (field.CC(lo) + field.CC(hi)));
	double gf = this->_NonDim->Mob().Face(hl, hh);
	const double tk = (gf == 0) ? 0.0 : gf / (hl + hh) / mu;
	gf /= mu;
	return this->_Magnitude * ( (gf * (hh - hl) + (field.SS(hi) - field.SS(lo)) * (tk * (hl + hh))) / (dn * dn) );
}


//------------------------------------------------------------------
//  Adding into the matrix A
//------------------------------------------------------------------
template<class T>
void VerticalGravityKernel<T>::SetMatrixA(CRS_Matrix &crsA, const Field &field)
{
	if(!this->_enabled) return;

	const double Coef = this->_Magnitude;
	int irow, s, d, k;
	double G[4], Tface[4], Tsum[4], Too;

	Site st;
	//  Per-thread workspace. It is created outside the parallel region and
	//  handed out with firstprivate.
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,irow,s,d,k,G,Tface,Tsum,Too) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			irow = T::Index(this->_G, st);
			if( crsA.IsBC(irow) ) continue;

			for(d = 0; d < T::nDir; ++d)
			{
				G[d]     = NonlinearTerm  (st, T::Dir(d), field);
				Tface[d] = TopographyShare(st, T::Dir(d), field);
			}

			this->_G->SetSite(svec, st);

			for(d = 0; d < T::nDir; ++d)
			{
				Tsum[d] = 0;
				for(s = 0; s < this->_n_stencil; s++)
					Tsum[d] += this->C(d)[s] * field.SS(svec[s]);
			}

			//  The center's share of each face flux, K of that face (the
			//  same as the neighbor's), with the viscosity of that face. The
			//  summation runs in stencil order (Ym, Xm, Xp, Yp in 2-D).
			Too = 0;
			for(k = 0; k < T::nDir; ++k) { const int dd = T::StencilDir(k); Too += Tsum[dd] * Tface[dd]; }

			for(d = 0; d < T::nDir; ++d) Tsum[d] *= Tface[d];

			for(k = 0; k < T::nNarrow; ++k)
			{
				const int idx = T::NarrowIdx(k);
				double val = 0;
				// The summation of G runs in axis order (Xm, Xp, Ym, Yp)
				for(d = 0; d < T::nDir; ++d) val += this->C(d)[idx] * G[d];

				if(T::IsDiag(idx)) val += Too;
				else               val += Tsum[T::DirOfAdj(idx)];

				if(T::IsDiag(idx)) crsA.Diag(irow) += Coef*val;
				else               crsA(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx))) += Coef*val;
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}


//------------------------------------------------------------------
//  Adding into the Jacobian
//------------------------------------------------------------------
template<class T>
void VerticalGravityKernel<T>::SetJacobiMatrix(CRS_Matrix &crsJ, const Field &field)
{
	if(!this->_enabled) return;

	const double Coef = this->_Magnitude;
	int irow, s, d, k;
	double a[4], b[4], gAdj[4], gCen[4], kF[4], gaSum;

	Site st;
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,irow,s,d,k,a,b,gAdj,gCen,kF,gaSum) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			irow = T::Index(this->_G, st);
			if( crsJ.IsBC(irow) ) continue;

			this->_G->SetSite(svec, st);

			for(d = 0; d < T::nDir; ++d)
			{
				a[d] = b[d] = 0;
				for(s = 0; s < this->_n_stencil; s++)
				{
					a[d] += this->C(d)[s] * field.HH(svec[s]);
					b[d] += this->C(d)[s] * field.SS(svec[s]);
				}
				gAdj[d]  = NonlinearJacobi(st, T::Dir(d), (tAdjacent)T::AdjIdx(d),   field);
				gCen[d]  = NonlinearJacobi(st, T::Dir(d), (tAdjacent)T::CenterIdx(), field);
				kF[d]    = TopographyShare(st, T::Dir(d), field);
			}

			//  The flux of face d is G_d (a_d + b_d) (a, b the differences of h
			//  and s); its derivative at point j is dG_d/dh_j (a_d + b_d) plus
			//  G_d C(d)[j], of which A holds G_d C(d)[j] and K_d b_d: added
			//  here is dG_d/dh_j (a_d + b_d) - K_d b_d. The diagonal sum runs
			//  in stencil order.
			gaSum = 0;
			for(k = 0; k < T::nDir; ++k)
			{
				const int dd = T::StencilDir(k);
				gaSum += gCen[dd] * (a[dd] + b[dd]) - kF[dd] * b[dd];
			}

			for(k = 0; k < T::nNarrow; ++k)
			{
				const int idx = T::NarrowIdx(k);
				if(T::IsDiag(idx))
				{
					crsJ.Diag(irow) += Coef*gaSum;
				}
				else
				{
					const int dd = T::DirOfAdj(idx);
					crsJ(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx)))
						+= Coef*( gAdj[dd]*(a[dd] + b[dd]) - kF[dd]*b[dd] );
				}
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}


//------------------------------------------------------------------
//  The velocity field
//------------------------------------------------------------------
template<class T>
void VerticalGravityKernel<T>::CalVelocity(Field &field)
{
	int s, a;
	double Q, Gf;

	Site st;
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,s,a,Q,Gf) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			if(!T::InsideForVelocity(this->_G, this->_Nx, this->_Ny, st)) continue;

			this->_G->SetSite(svec, st);

			for(a = 0; a < T::nAxis; ++a)
			{
				Gf = NonlinearTerm(st, T::Dir(2*a), field);

				Q = 0;
				for(s = 0; s < this->_n_stencil; s++)
				{
					if(T::SkipInVel(s)) continue;
					Q += this->CV(a)[s] * (field.HH(svec[s]) + field.SS(svec[s]));
				}
				Q *= this->_Magnitude * Gf;

				Q = T::AxisMask(this->_Ny, a, Q);
				T::WriteFaceVel(this->_G, field.UU().U_VerticalGravity(), st, a,
				                Q / T::FaceHeight(this->_G, field, st, a));
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_VERTICAL_GRAVITY_KERNEL_INCLUDED_
