/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LaplacePressureKernel.h : the dimension-independent implementation of the
*                            capillary pressure term
*
*  ** On the order of summation - change it and the results stop agreeing
*  bit for bit **
*  There are two orderings of the directions, and on the diagonal of
*  SetJacobiMatrix both of them appear within a single expression.
*    gBo*aB + gLo*aL + gRo*aR + gTo*aT   <- stencil order (Ym, Xm, Xp, Yp)
*    sLo*bL + sRo*bR + sBo*bB + sTo*bT   <- axis order    (Xm, Xp, Ym, Yp)
*  Use the traits' StencilDir(k) / Dir(d) exactly as shown here.
*
*  SetLinearCoef and CalVelocity are left dimension-specific. The 2-D
*  biharmonic operator has a cross term (1/dx^2dy^2) with no 1-D
*  counterpart, and the velocity computation writes out a third-order
*  difference directly.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _LAPLACE_PRESSURE_KERNEL_INCLUDED_
#define _LAPLACE_PRESSURE_KERNEL_INCLUDED_

#include "PhysicsTermBase.h"
#include "NonDimensionalNumbers.h"
#include "BoundaryCondition.h"
#include "LubFields.h"
#include "MyMath.h"
#include "Opt.h"
#include <cmath>

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

template<class T>
class LaplacePressureKernel : public PhysicsTermBase<T>
{
public:
	typedef typename T::Grid      Grid;
	typedef typename T::Site      Site;
	typedef typename T::Field     Field;
	typedef typename T::tEdge     tEdge;
	typedef typename T::tAdjacent tAdjacent;
	typedef typename T::Boundary  Boundary;

	LaplacePressureKernel() : PhysicsTermBase<T>() {}

	//  Capillary pressure: the pressure produced by the curvature of the
	//  free surface h+s. The coefficient is fixed at 1/3. The coefficient
	//  _Magnitude (sign included) and the finite-difference coefficients are
	//  decided here. For the contract, see PhysicsTermBase.h.
	/**
	 * @brief Switch the capillary pressure term on or off and build its
	 *        finite-difference coefficients.
	 * @param inp_phys  the Physics block of the input file
	 * @param g         grid the term is discretized on
	 * @param nonDim    the set of non-dimensional numbers
	 * @param bc        boundary condition, used to fold the ghost points into the stencil
	 *
	 * The overall coefficient _Magnitude is fixed at 1/3; the term is driven
	 * by the curvature of the free surface h+s, so it carries no
	 * non-dimensional number of its own.
	 */
	virtual void Build(Input::I_Physics &inp_phys, Grid *g, NonDimensionalNumbers *nonDim, Boundary *bc)
	{
		this->_G = g;
		this->_NonDim = nonDim;
		this->_boundary = bc;

		this->_enabled = inp_phys.LaplacePressure.Enable;
		this->_Magnitude = 1.0 / 3.0;

		T::SetExtent(this->_G, this->_Nx, this->_Ny);
		this->SetLinearCoef();
	}

	/**
	 * @brief Add the capillary pressure to the matrix A.
	 *
	 * Being a fourth-order operator, this term reaches every point of the
	 * stencil, so it inserts over T::FullIdx rather than the narrow set.
	 */
	virtual void SetMatrixA(CRS_Matrix &crsA, const Field &field);
	/**
	 * @brief Add the capillary pressure to the Jacobian.
	 *
	 * Only the narrow stencil is filled here; see the note at the top of this
	 * file for the summation order the diagonal relies on.
	 */
	virtual void SetJacobiMatrix(CRS_Matrix &crsJ, const Field &field);

	/**
	 * @brief A h of this term evaluated face by face (see
	 *        PhysicsTermBase::ConservativeForm).
	 *
	 * The flux through a face is Coef D_f (G_f(h) + G_f(s)), D_f the
	 * mobility SetMatrixA last took for that face and G_f the third
	 * difference of the film (FaceGradient), written as differences of
	 * differences. The two cells that share a face compute it from the same
	 * values in the same order, so one takes out exactly what the other
	 * puts in. With the entries of A the sum over the cells came out as
	 * eps |A| |h| per cell -- on a film 1.2 mm thick over cells 3.3 um
	 * wide, the nozzle_1d_rotating case, a volume error of 2e-7 of the
	 * film. Mirror ghosts (Symmetry) make G_f of an edge face exactly 0.
	 * In exact arithmetic it is the product of A with h: the topography
	 * part, which A puts on the diagonal as (D_f/h_i) G_f(s) h_i, is
	 * D_f G_f(s) here.
	 */
	virtual bool ConservativeForm() const { return this->_enabled; }
	virtual void AddConservative(Vector<double> &y, const Vector<double> &x,
	                             const CRS_Matrix &crs, const Field &field) const;

	/**
	 * @brief The flux of AddConservative through face d of s, with the
	 *        mobility and the film of the state field.
	 */
	virtual double FaceFlux(Site &s, int d, const Field &field, const std::vector<Site> &svec);

protected:
	//  The mobility D of each row's faces as SetMatrixA last assembled A
	//  with, for AddConservative; indexed [direction][row]
	Vector<double> _dFace[T::nDir];

	//  opt assembly (util/Opt.h): the surface tension and the viscosity of each
	//  face of each row as SetMatrixA last evaluated them, for the Jacobian
	//  assembled right after it from the same field (the same values, so
	//  the same bits); _propValid says they are there.
	Vector<double> _sigFace[T::nDir], _visFace[T::nDir];
	bool _propValid = false;

	/** @brief The mobility of NonlinearTerm, with Sigma and Visco returned. */
	double NonlinearTermProps(Site &s, tEdge tedge, const Field &field, double &sig, double &vis)
	{
		sig = vis = 1.0;
		if(PinnedFace(s, tedge)) return 0;
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
		double denom = Hoo + Haj;
		if(Math::IsZero(denom))
		{
			Global::FatalError<<"C_LaplacePressure_WM::NonlinearTerm(), denom=0"<<nl;
			Global::FatalError.Abort();
		}
		//  With slip, the mobility D = h^3 + a h^2 of Mobility.h
		double val =  this->_NonDim->Mob().Slip() ? this->_NonDim->Mob().FaceDK(Hoo, Haj)
		                                          : 2.0*(Hoo*Hoo)*(Haj*Haj) / denom;
		sig = this->_NonDim->SurfaceTension(Cfc);
		vis = this->_NonDim->Viscosity(Cfc);
		return val * sig / vis;
	}

	/** @brief NonlinearJacobi with the face's Sigma and Visco given. */
	double NonlinearJacobiProps(Site &s, tEdge tedge, tAdjacent tc, const Field &field, double sig, double vis)
	{
		if(PinnedFace(s, tedge)) return 0;
		double Hoo = field.HH(s);
		Site a = T::Adj(this->_G, s, tedge);
		double Haj = field.HH(a);
		double H1, H2;
		if((int)tc == T::CenterIdx()) { H1 = Hoo; H2 = Haj; }
		else                           { H1 = Haj; H2 = Hoo; }
		double denom = pow(H1 + H2, 2.0);
		if(Math::IsZero(denom)) return NonlinearJacobi(s, tedge, tc, field);	// its warning
		double val =  this->_NonDim->Mob().Slip() ? this->_NonDim->Mob().FaceDK_d1(H1, H2)
		                                          : 2.0* (H1*H2*H2) * (H1+H2+H2) / denom;
		return val * sig / vis;
	}

	/**
	 * @brief The third difference G_f of the values v (on the stencil of
	 *        the row, indexed as T's stencil) across face d of that row,
	 *        with the sign of the outward normal of the face of the cell
	 *        on its low side: C(d) v is -G_f on a low face (Xm, Ym) and
	 *        +G_f on a high face. It is computed only from the values about
	 *        the face, in an order that does not depend on which of the two
	 *        cells asks.
	 */
	virtual double FaceGradient(const double *v, int d) const = 0;

	//--- dimension-specific ---
	/**
	 * @brief Whether the face tedge touches a pinned boundary, on which the
	 *        mobility is set to zero.
	 */
	virtual bool PinnedFace(Site &s, tEdge tedge) const = 0;

	/**
	 * @brief The mobility on a face: the coefficient of the form h³σ/μ that
	 *        multiplies the curvature.
	 */
	virtual double NonlinearTerm(Site&, tEdge, const Field &field);
	/**
	 * @brief The above, differentiated with respect to one point (for the
	 *        Jacobian).
	 */
	virtual double NonlinearJacobi(Site&, tEdge, tAdjacent, const Field &field);

	/**
	 * @brief The contribution of the substrate topography, which enters as
	 *        the mobility divided by the film height.
	 */
	double TopographyTerm(Site&, tEdge, const Field &field);
	/**
	 * @brief The above, differentiated with respect to one point.
	 */
	double TopographyJacobi(Site&, tEdge, tAdjacent, const Field &field);
};


//==================================================================
//  Implementation
//==================================================================

//------------------------------------------------------------------
//  The mobility on a face (positivity preserving scheme)
//
//  Equation (23) of Diez & Kondic, J. Comput. Phys. 183 (2002) 274-306:
//
//      D_k^(x) = 2 h_k^2 h_{k+1}^2 / (h_k + h_{k+1})
//
//  This is the general form of equations (21)(22) of the same paper,
//
//      D_i = (h_{i+1} - h_i) / (g_{i+1} - g_i),   g(h) = ∫ dh/D(h)
//
//  with D(h) = h^s and s = 3 substituted in.
//
//  ** Do not replace this with an arithmetic mean (such as
//  (h_k^3 + h_{k+1}^3)/2) **
//  The form above is the product of h_k h_{k+1} and the harmonic mean
//  2h_k h_{k+1}/(h_k+h_{k+1}), so it goes to zero as O(h^2) when either
//  side approaches zero. The flux vanishes on a face that borders a dry
//  region, which is what keeps the numerical solution from stepping into
//  h < 0. With an arithmetic mean a mobility of about h^3/24 survives even
//  when one side is zero; that draws liquid out of the dry region and
//  produces negative film heights. In a computation that contains a
//  contact line, this is the direct cause of the breakdown.
//  (The property is pinned down by tests/test_positivity_mobility.cpp.)
//------------------------------------------------------------------
template<class T>
double LaplacePressureKernel<T>::NonlinearTerm(Site &s, tEdge tedge, const Field &field)
{
	if(PinnedFace(s, tedge)) return 0;

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
	double denom = Hoo + Haj;

	if(Math::IsZero(denom))
	{
		Global::FatalError<<"C_LaplacePressure_WM::NonlinearTerm(), denom=0"<<nl;
		Global::FatalError.Abort();
	}

	//  With slip, the mobility D = h^3 + a h^2 of Mobility.h
		double val =  this->_NonDim->Mob().Slip() ? this->_NonDim->Mob().FaceDK(Hoo, Haj)
		                                          : 2.0*(Hoo*Hoo)*(Haj*Haj) / denom;
	double Sigma = this->_NonDim->SurfaceTension(Cfc);
	double Visco = this->_NonDim->Viscosity(Cfc);
	return val * Sigma /  Visco;
}


//------------------------------------------------------------------
//  The above, differentiated with respect to one point (tc).
//  When tc is the center, (H1,H2) = (this cell, the neighbor); when tc is
//  the neighbor, the two are swapped.
//------------------------------------------------------------------
template<class T>
double LaplacePressureKernel<T>::NonlinearJacobi(Site &s, tEdge tedge, tAdjacent tc, const Field &field)
{
	if(PinnedFace(s, tedge)) return 0;

	double Hoo = field.HH(s);
	double Coo = field.CC(s);
	double Haj = 0, Caj = 0, H1 = 0, H2 = 0;

	if(!T::IsFace(tedge))
	{
		WarningIn("NonlinearJacobi()");
	}
	else
	{
		Site a = T::Adj(this->_G, s, tedge);
		Haj = field.HH(a);
		Caj = field.CC(a);

		if     ((int)tc == T::CenterIdx())                { H1 = Hoo; H2 = Haj; }
		else if((int)tc == T::AdjIdx(T::DirIndex(tedge))) { H1 = Haj; H2 = Hoo; }
		else                                              { WarningIn("NonlinearJacobi()"); }
	}

	double Cfc = 0.5*(Coo + Caj);
	double denom = pow(H1 + H2, 2.0);

	if(Math::IsZero(denom))
		WarningIn("NonlinearJacobi(), denom=0"), Global::Warning<<"tc="<<(int)tc<<nl;

	double Sigma = this->_NonDim->SurfaceTension(Cfc);
	double Visco = this->_NonDim->Viscosity(Cfc);
	double val =  this->_NonDim->Mob().Slip() ? this->_NonDim->Mob().FaceDK_d1(H1, H2)
		                                          : 2.0* (H1*H2*H2) * (H1+H2+H2) / denom;
	return val * Sigma / Visco;
}


//------------------------------------------------------------------
//  The substrate topography term. It is defined as D/hoo.
//------------------------------------------------------------------
template<class T>
double LaplacePressureKernel<T>::TopographyTerm(Site &s, tEdge tedge, const Field &field)
{
	double Hoo = field.HH(s);
	return NonlinearTerm(s, tedge, field) / Hoo;
}

template<class T>
double LaplacePressureKernel<T>::TopographyJacobi(Site &s, tEdge tedge, tAdjacent tc, const Field &field)
{
	double Hoo = field.HH(s);
	double val = NonlinearJacobi(s, tedge, tc, field) / Hoo;
	double add = 0;

	if(tedge != T::DirCenter() && (int)tc == T::CenterIdx())
		add = NonlinearTerm(s, tedge, field) / (Hoo*Hoo);
	return val - add;
}


//------------------------------------------------------------------
//  Adding into the matrix A. The Laplace pressure reaches every point of
//  the stencil.
//------------------------------------------------------------------
template<class T>
void LaplacePressureKernel<T>::SetMatrixA(CRS_Matrix &crsA, const Field &field)
{
	if(!this->_enabled) return;

	const double Coef = this->_Magnitude;
	int irow, s, d, k;
	double D[4], Tf[4], topo;

	for(d = 0; d < T::nDir; ++d)
		if(_dFace[d].size() != crsA.RowSize()) { _dFace[d].resize(crsA.RowSize()); _dFace[d] = 0; }
	const bool cache = yuzen::opt::On(yuzen::opt::Assembly);
	if(cache)
		for(d = 0; d < T::nDir; ++d)
			if(_sigFace[d].size() != crsA.RowSize()) { _sigFace[d].resize(crsA.RowSize()); _visFace[d].resize(crsA.RowSize()); }

	Site st;
	//  Per-thread workspace. It is created outside the parallel region and
	//  handed out with firstprivate.
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,irow,s,d,k,D,Tf,topo) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			irow = T::Index(this->_G, st);
			if( crsA.IsBC(irow) ) continue;

			if(cache)
			{
				//  The face properties once; TopographyTerm is NonlinearTerm / h
				const double Hoo = field.HH(st);
				for(d = 0; d < T::nDir; ++d)
				{
					double sg, vs;
					D[d]  = NonlinearTermProps(st, T::Dir(d), field, sg, vs);
					Tf[d] = D[d] / Hoo;
					_dFace[d][irow] = D[d];
					_sigFace[d][irow] = sg;
					_visFace[d][irow] = vs;
				}
			}
			else
			for(d = 0; d < T::nDir; ++d)
			{
				D[d]  = NonlinearTerm (st, T::Dir(d), field);
				Tf[d] = TopographyTerm(st, T::Dir(d), field);
				_dFace[d][irow] = D[d];
			}

			this->_G->SetSite(svec, st);

			// The substrate topography term, added to the diagonal only.
			// The inner summation over the directions is in axis order.
			topo = 0;
			for(s = 0; s < this->_n_stencil; s++)
			{
				double c = 0;
				for(d = 0; d < T::nDir; ++d) c += this->C(d)[s] * Tf[d];
				topo += c * field.SS(svec[s]);
			}

			for(k = 0; k < T::nFull; ++k)
			{
				const int idx = T::FullIdx(k);
				double val = 0;
				for(d = 0; d < T::nDir; ++d) val += this->C(d)[idx] * D[d];	// axis order

				if(T::IsDiag(idx)) { val += topo; crsA.Diag(irow) += Coef*val; }
				else               crsA(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx))) += Coef*val;
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
	_propValid = cache;
}


//------------------------------------------------------------------
//  A h face by face (see the declaration)
//------------------------------------------------------------------
template<class T>
void LaplacePressureKernel<T>::AddConservative(Vector<double> &y, const Vector<double> &x,
                                               const CRS_Matrix &crs, const Field &field) const
{
	if(!this->_enabled) return;
	const double Coef = this->_Magnitude;
	Site st;
	int irow, s, d;
	double vh[T::nStencil], vs[T::nStencil], r, f;
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,irow,s,d,vh,vs,r,f) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			irow = T::Index(this->_G, st);
			if( crs.IsBC(irow) ) continue;

			this->_G->SetSite(svec, st);
			for(s = 0; s < this->_n_stencil; s++)
			{
				vh[s] = x[T::Index(this->_G, svec[s])];
				vs[s] = field.SS(svec[s]);
			}
			r = 0;
			for(d = 0; d < T::nDir; ++d)
			{
				f = Coef * (_dFace[d][irow] * (FaceGradient(vh, d) + FaceGradient(vs, d)));
				r += (d % 2 == 0) ? -f : f;
			}
			y[irow] += r;
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}


template<class T>
double LaplacePressureKernel<T>::FaceFlux(Site &st, int d, const Field &field, const std::vector<Site> &svec)
{
	if(!this->_enabled) return 0;
	double vh[T::nStencil], vs[T::nStencil];
	for(int s = 0; s < this->_n_stencil; s++)
	{
		vh[s] = field.HH(svec[s]);
		vs[s] = field.SS(svec[s]);
	}
	return this->_Magnitude * (NonlinearTerm(st, T::Dir(d), field) * (FaceGradient(vh, d) + FaceGradient(vs, d)));
}


//------------------------------------------------------------------
//  Adding into the Jacobian. This one fills only the narrow stencil.
//------------------------------------------------------------------
template<class T>
void LaplacePressureKernel<T>::SetJacobiMatrix(CRS_Matrix &crsJ, const Field &field)
{
	if(!this->_enabled) return;

	const double Coef = this->_Magnitude;
	int irow, s, d, k;
	double a[4], b[4], gAdj[4], gCen[4], sAdj[4], sCen[4], Hoo, gaSum, sbSum;
	const bool cache = _propValid && yuzen::opt::On(yuzen::opt::Assembly)
	                && _sigFace[0].size() == crsJ.RowSize();

	Site st;
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,irow,s,d,k,a,b,gAdj,gCen,sAdj,sCen,Hoo,gaSum,sbSum) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			irow = T::Index(this->_G, st);
			if( crsJ.IsBC(irow) ) continue;

			this->_G->SetSite(svec, st);
			Hoo = field.HH(st);

			for(d = 0; d < T::nDir; ++d)
			{
				a[d] = b[d] = 0;
				for(s = 0; s < this->_n_stencil; s++)
				{
					a[d] += this->C(d)[s] * field.HH(svec[s]);
					b[d] += this->C(d)[s] * field.SS(svec[s]);
				}
				const tAdjacent iAdj = (tAdjacent)T::AdjIdx(d);
				const tAdjacent iCen = (tAdjacent)T::CenterIdx();
				if(cache)
				{
					//  The face's Sigma and Visco as SetMatrixA had them;
					//  TopographyJacobi is NonlinearJacobi / h, less
					//  NonlinearTerm / h^2 at the center
					const double sg = _sigFace[d][irow], vs = _visFace[d][irow];
					gAdj[d] = NonlinearJacobiProps(st, T::Dir(d), iAdj, field, sg, vs);
					gCen[d] = NonlinearJacobiProps(st, T::Dir(d), iCen, field, sg, vs);
					sAdj[d] = gAdj[d] / Hoo;
					sCen[d] = gCen[d] / Hoo - _dFace[d][irow] / (Hoo*Hoo);
				}
				else
				{
				gAdj[d] = NonlinearJacobi (st, T::Dir(d), iAdj, field);
				gCen[d] = NonlinearJacobi (st, T::Dir(d), iCen, field);
				sAdj[d] = TopographyJacobi(st, T::Dir(d), iAdj, field);
				sCen[d] = TopographyJacobi(st, T::Dir(d), iCen, field);
				}
			}

			// The two summations on the diagonal run in different orders
			// (see the warning at the top of this file)
			gaSum = 0;
			for(k = 0; k < T::nDir; ++k) { const int dd = T::StencilDir(k); gaSum += gCen[dd]*a[dd]; }
			sbSum = 0;
			for(d = 0; d < T::nDir; ++d) sbSum += sCen[d]*b[d];

			for(k = 0; k < T::nNarrow; ++k)
			{
				const int idx = T::NarrowIdx(k);
				if(T::IsDiag(idx))
				{
					crsJ.Diag(irow) += Coef*( gaSum + sbSum*Hoo );
				}
				else
				{
					const int dd = T::DirOfAdj(idx);
					crsJ(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx)))
						+= Coef*( gAdj[dd]*a[dd] + sAdj[dd]*b[dd]*Hoo );
				}
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_LAPLACE_PRESSURE_KERNEL_INCLUDED_
