/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  SolutoCapillarityKernel.h : the dimension-independent implementation of
*                              the soluto-capillarity term
*
*  For everything around the four-point upwind interpolation, see
*  UpwindKernelBase.h.
*
*  Mind the order of summation (the rule at the top of UpwindKernelBase.h):
*    the additions into Val[] in SetMatrixA : stencil order
*                                             (Ym, Xm, Xp, Yp in 2-D)
*    the per-row ordering in SetJacobiMatrix : stencil order on the
*                                             diagonal, the neighboring
*                                             direction first elsewhere
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _SOLUTO_CAPILLARITY_KERNEL_INCLUDED_
#define _SOLUTO_CAPILLARITY_KERNEL_INCLUDED_

#include "Opt.h"
#include "UpwindKernelBase.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

template<class T>
class SolutoCapillarityKernel : public UpwindKernelBase<T>
{
public:
	typedef typename T::Grid      Grid;
	typedef typename T::Site      Site;
	typedef typename T::Field     Field;
	typedef typename T::tEdge     tEdge;
	typedef typename T::tAdjacent tAdjacent;
	typedef typename T::Boundary  Boundary;

	/** @brief Select the four-point upwind interpolation scheme. */
	SolutoCapillarityKernel() : UpwindKernelBase<T>()
	{
		this->_scheme.Build(C_Interpolation_base4::KK);
		this->_enabled = true;
	}

	//  Soluto-capillarity (Marangoni): the surface stress produced by the
	//  dependence of surface tension on concentration. The coefficient is
	//  As*Ma/2. The coefficient _Magnitude (sign included) and the
	//  finite-difference coefficients are decided here.
	//  For the contract, see PhysicsTermBase.h.
	/**
	 * @brief Switch the Marangoni term on or off and build its
	 *        finite-difference coefficients.
	 * @param inp_phys  the Physics block of the input file
	 * @param g         grid the term is discretized on
	 * @param nonDim    the set of non-dimensional numbers; As and Ma are
	 *                  taken from here
	 * @param bc        boundary condition, used to fold the ghost points into the stencil
	 *
	 * The overall coefficient _Magnitude is As*Ma/2. The stress comes from
	 * the dependence of surface tension on concentration, so this term
	 * couples the film-height equation to the concentration field.
	 */
	virtual void Build(Input::I_Physics &inp_phys, Grid *g, NonDimensionalNumbers *nonDim, Boundary *bc)
	{
		this->_G = g;
		this->_NonDim = nonDim;
		this->_boundary = bc;

		this->_enabled = inp_phys.SolutoCapillarity.Enable;
		this->_Magnitude = this->_NonDim->As() * this->_NonDim->Ma() / 2.0;

		T::SetExtent(this->_G, this->_Nx, this->_Ny);
		SetLinearCoef();
	}

	/**
	 * @brief Add the Marangoni stress to the matrix A.
	 *
	 * Insertion runs over the wide stencil, because the upwinded four-point
	 * interpolation reaches two cells to either side of each face.
	 */
	virtual void SetMatrixA(CRS_Matrix &crsA, const Field &field);
	/**
	 * @brief Add the Marangoni stress to the Jacobian.
	 *
	 * The row-by-row ordering is the one laid out at the top of
	 * UpwindKernelBase.h and is applied by ScatterJacobiWide.
	 */
	virtual void SetJacobiMatrix(CRS_Matrix &crsJ, const Field &field);
	/**
	 * @brief Write the Marangoni velocity into field.UU().U_SolutoCapillarity().
	 */
	virtual void CalVelocity(Field &field);

	/**
	 * @brief The flux through face d of st in the state field, as A holds
	 *        it: Coef S_d (the Marangoni coefficient dc/mu times the face
	 *        value of the film) times that face value, per mesh spacing.
	 */
	virtual double FaceFlux(Site &st, int d, const Field &field, const std::vector<Site> &svec)
	{
		if(!this->_enabled) return 0;
		double dc, cc;
		this->FaceConc(svec, d, field, dc, cc);
		const double h = T::Spacing(this->_G, d/2);
		double S, hs = 0;
		if(this->Limited())
		{
			double f[4], wt[4];
			this->FaceHeights(svec, d, field, f);
			const double up = Up(dc);
			const double hf = this->_scheme.FaceValue(f, up);
			this->_scheme.FaceWeights(f, up, -1, wt);
			S = this->WallFace(st, T::Dir(d)) ? 0.0 : dc / this->_NonDim->Viscosity(cc)
			    * (this->_NonDim->Mob().Slip() ? this->_NonDim->Mob().Shear(hf) : hf);
			for(int q = 0; q < 4; ++q) hs += wt[q] * f[q];
		}
		else
		{
			//  S as NonlinearTerm forms it, from the stencil given
			for(int q = 0; q < 4; ++q)
			{
				const int idx = T::FaceStencil(d, q);
				hs += (this->C(d)[idx] + Up(dc)*this->C2(d)[idx]) * field.HH(svec[idx]);
			}
			S = this->WallFace(st, T::Dir(d)) ? 0.0 : dc / this->_NonDim->Viscosity(cc)
			    * (this->_NonDim->Mob().Slip() ? this->_NonDim->Mob().Shear(hs) : hs);
		}
		return this->_Magnitude * (S * hs / h);
	}

protected:
	/**
	 * @brief Build the four-point upwind interpolation coefficients.
	 *
	 * They are not divided by the mesh spacing; that is this term's
	 * convention.
	 */
	virtual void SetLinearCoef() { this->SetUpwindCoef(false); }
	/**
	 * @brief The surface stress times the film height, on a face.
	 *
	 * The upwind side is chosen by the direction of the flux, Up(dc).
	 */
	virtual double NonlinearTerm(Site&, tEdge, const Field &field);
	/**
	 * @brief The above, differentiated with respect to one point (for the
	 *        Jacobian).
	 */
	virtual double NonlinearJacobi(Site&, tEdge, tAdjacent, const Field &field);

	/**
	 * @brief The direction of the Marangoni flux through a face, which picks
	 *        the upwind side: that of As Ma dc, not of dc alone.
	 *
	 * The flux is As Ma/2 h^2 dc/mu, so with a surface tension that falls
	 * with the concentration (Fluid.SurfaceTension.Coef < 0) the film moves
	 * against dc. Upwinding by the sign of dc alone would then take the
	 * downwind side: the third-order dissipation of the upwinded
	 * interpolation would change sign and become anti-dissipation. With
	 * Coef >= 0 this is Sign(dc) to the last bit.
	 */
	double Up(double dc) const { return Math::Sign(this->_Magnitude) * Math::Sign(dc); }

	//  A limited scheme (UpwindKernelBase.h): the flux Coef h_f^2 dc/mu is
	//  held in A as (sum_k wt_k h_k) h_f dc/mu, and the Jacobian adds
	//  (2 g_k - wt_k) h_f dc/mu, g_k = d h_f / d h_k.
	void SetMatrixA_Limited(CRS_Matrix &crsA, const Field &field);
	void SetJacobiMatrix_Limited(CRS_Matrix &crsJ, const Field &field);
};


//==================================================================
//  Implementation
//==================================================================

template<class T>
double SolutoCapillarityKernel<T>::NonlinearTerm(Site &st, tEdge tedge, const Field &field)
{
	if(!T::IsFace(tedge)) return 0;
	//  No flux through the face of a wall (UpwindKernelBase::WallFace)
	if(this->WallFace(st, tedge)) return 0;
	const int d = T::DirIndex(tedge);

	std::vector<Site> svec(this->_n_stencil);
	this->_G->SetSite(svec, st);

	double dc, cc;
	this->FaceConc(svec, d, field, dc, cc);

	double hh = 0;
	if(this->Limited())
	{
		double f[4];
		this->FaceHeights(svec, d, field, f);
		hh = this->_scheme.FaceValue(f, Up(dc));
	}
	else
	for(int s = 0; s < this->_n_stencil; s++)
		hh += (this->C(d)[s] + Up(dc) * this->C2(d)[s]) * field.HH(svec[s]);

	double Visco = this->_NonDim->Viscosity(cc);
	//  With slip the shear mobility over h is h + 2 Lambda (Mobility.h)
	if(this->_NonDim->Mob().Slip()) return dc / Visco * this->_NonDim->Mob().Shear(hh);
	return dc / Visco * hh;
}


template<class T>
double SolutoCapillarityKernel<T>::NonlinearJacobi(Site &st, tEdge tedge, tAdjacent tc, const Field &field)
{
	if(!T::IsFace(tedge)) return 0;
	//  No flux through the face of a wall (UpwindKernelBase::WallFace)
	if(this->WallFace(st, tedge)) return 0;
	const int d = T::DirIndex(tedge);

	std::vector<Site> svec(this->_n_stencil);
	this->_G->SetSite(svec, st);

	double dc, cc;
	this->FaceConc(svec, d, field, dc, cc);

	double dHrdH = this->C(d)[(int)tc] + Up(dc) * this->C2(d)[(int)tc];
	double Visco = this->_NonDim->Viscosity(cc);
	return dc / Visco * dHrdH;
}


template<class T>
void SolutoCapillarityKernel<T>::SetMatrixA(CRS_Matrix &crsA, const Field &field)
{
	if(!this->_enabled) return;
	if(this->Limited()) { SetMatrixA_Limited(crsA, field); return; }

	const double Coef = this->_Magnitude;
	int irow, s, j, k;
	double S[4], dc[4], cc;

	Site st;
	//  Per-thread workspace. It is created outside the parallel region and
	//  handed out with firstprivate.
	std::vector<Site>  svec(this->_n_stencil);
	Vector<double> Val(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,irow,s,j,k,S,dc,cc) firstprivate(svec,Val)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			irow = T::Index(this->_G, st);
			if( crsA.IsBC(irow) ) continue;

			this->_G->SetSite(svec, st);

			if(yuzen::opt::On(yuzen::opt::Assembly))
			{
				//  NonlinearTerm written out on this row's stencil: no
				//  stencil of its own to allocate, the face concentration
				//  once
				for(j = 0; j < T::nDir; ++j)
				{
					this->FaceConc(svec, j, field, dc[j], cc);
					if(!T::IsFace(T::Dir(j)) || this->WallFace(st, T::Dir(j))) { S[j] = 0; continue; }
					double hh = 0;
					for(s = 0; s < this->_n_stencil; s++)
						hh += (this->C(j)[s] + Up(dc[j]) * this->C2(j)[s]) * field.HH(svec[s]);
					double Visco = this->_NonDim->Viscosity(cc);
					S[j] = this->_NonDim->Mob().Slip() ? dc[j] / Visco * this->_NonDim->Mob().Shear(hh)
					                                   : dc[j] / Visco * hh;
				}
			}
			else
			for(j = 0; j < T::nDir; ++j)
			{
				S[j] = NonlinearTerm(st, T::Dir(j), field);
				this->FaceConc(svec, j, field, dc[j], cc);
			}

			Val = 0;
			// The additions run in stencil order. The negative-side
			// directions carry a negative sign.
			for(s = 0; s < this->_n_stencil; s++)
				for(k = 0; k < T::nDir; ++k)
				{
					const int d = T::StencilDir(k);
					const double h = T::Spacing(this->_G, d/2);
					const double sg = (d % 2 == 0) ? -1.0 : +1.0;
					Val[s] += sg*(this->C(d)[s] + Up(dc[d])*this->C2(d)[s])*S[d]/h;
				}

			for(k = 0; k < T::nWide; ++k)
			{
				const int idx = T::WideIdx(k);
				if(T::IsDiag(idx)) crsA.Diag(irow) += Coef*Val[idx];
				else               crsA(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx))) += Coef*Val[idx];
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}


namespace detail {
//  The "coefficient lookup" handed to ScatterJacobiWide
template<class K>
struct SolutoG
{
	K *self; const double *g;
	double operator()(int d, int idx) const { (void)self; return g[d*16 + idx]; }
};
}


template<class T>
void SolutoCapillarityKernel<T>::SetJacobiMatrix(CRS_Matrix &crsJ, const Field &field)
{
	if(!this->_enabled) return;
	if(this->Limited()) { SetJacobiMatrix_Limited(crsJ, field); return; }

	const double Coef = this->_Magnitude;
	int irow, s, d, k;
	double a[4], dc[4], cc;
	double g[4*16];

	Site st;
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,irow,s,d,k,a,dc,cc,g) firstprivate(svec)
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
				this->FaceConc(svec, d, field, dc[d], cc);
				const double h  = T::Spacing(this->_G, d/2);
				const double sg = (d % 2 == 0) ? -1.0 : +1.0;
				a[d] = 0;
				for(s = 0; s < this->_n_stencil; s++)
					a[d] += sg*(this->C(d)[s] + Up(dc[d])*this->C2(d)[s]) * field.HH(svec[s])/h;

				if(yuzen::opt::On(yuzen::opt::Assembly))
				{
					//  NonlinearJacobi with this face's dc and c, its
					//  viscosity once
					const bool none = !T::IsFace(T::Dir(d)) || this->WallFace(st, T::Dir(d));
					const double Visco = none ? 1.0 : this->_NonDim->Viscosity(cc);
					for(k = 0; k < 4; ++k)
					{
						const int idx = T::FaceStencil(d, k);
						double dHrdH = this->C(d)[idx] + Up(dc[d]) * this->C2(d)[idx];
						g[d*16 + idx] = none ? 0.0 : dc[d] / Visco * dHrdH;
					}
				}
				else
				for(k = 0; k < 4; ++k)
				{
					const int idx = T::FaceStencil(d, k);
					g[d*16 + idx] = NonlinearJacobi(st, T::Dir(d), (tAdjacent)idx, field);
				}
			}

			detail::SolutoG<SolutoCapillarityKernel<T> > gv; gv.self = this; gv.g = g;
			this->ScatterJacobiWide(crsJ, irow, st, gv, a, Coef);
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}


template<class T>
void SolutoCapillarityKernel<T>::CalVelocity(Field &field)
{
	int a;
	double U;

	Site st;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(st,a,U)
#endif
	for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
	{
		st = this->_G->template Of<Site>(ii);
		if(!T::InsideForVelocity(this->_G, this->_Nx, this->_Ny, st)) continue;

		for(a = 0; a < T::nAxis; ++a)
		{
			U = this->_Magnitude * NonlinearTerm(st, T::Dir(2*a), field);
			U = T::AxisMask(this->_Ny, a, U);
			T::WriteFaceVel(this->_G, field.UU().U_SolutoCapillarity(), st, a, U);
		}
	}
}

template<class T>
void SolutoCapillarityKernel<T>::SetMatrixA_Limited(CRS_Matrix &crsA, const Field &field)
{
	const double Coef = this->_Magnitude;
	Site st;
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			const int irow = T::Index(this->_G, st);
			if( crsA.IsBC(irow) ) continue;
			this->_G->SetSite(svec, st);

			double Val[16];
			for(int k = 0; k < 16; ++k) Val[k] = 0;
			for(int d = 0; d < T::nDir; ++d)
			{
				double dc, cc, f[4], wt[4];
				this->FaceConc(svec, d, field, dc, cc);
				this->FaceHeights(svec, d, field, f);
				const double up = Up(dc);
				const double hf = this->_scheme.FaceValue(f, up);
				this->_scheme.FaceWeights(f, up, -1, wt);
				const double S  = this->WallFace(st, T::Dir(d)) ? 0.0
				                : this->_NonDim->Mob().Slip()
				                  ? dc / this->_NonDim->Viscosity(cc) * this->_NonDim->Mob().Shear(hf)
				                  : dc / this->_NonDim->Viscosity(cc) * hf;
				const double h  = T::Spacing(this->_G, d/2);
				const double sg = (d % 2 == 0) ? -1.0 : +1.0;
				for(int q = 0; q < 4; ++q) Val[T::FaceStencil(d, q)] += sg*wt[q]*S/h;
			}
			for(int k = 0; k < T::nWide; ++k)
			{
				const int idx = T::WideIdx(k);
				if(T::IsDiag(idx)) crsA.Diag(irow) += Coef*Val[idx];
				else               crsA(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx))) += Coef*Val[idx];
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}
#endif
}

template<class T>
void SolutoCapillarityKernel<T>::SetJacobiMatrix_Limited(CRS_Matrix &crsJ, const Field &field)
{
	const double Coef = this->_Magnitude;
	Site st;
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			const int irow = T::Index(this->_G, st);
			if( crsJ.IsBC(irow) ) continue;
			this->_G->SetSite(svec, st);

			double a[4], g[4*16];
			for(int k = 0; k < 4*16; ++k) g[k] = 0;
			for(int d = 0; d < T::nDir; ++d)
			{
				double dc, cc, f[4], wt[4], gd[4];
				this->FaceConc(svec, d, field, dc, cc);
				this->FaceHeights(svec, d, field, f);
				const double up = Up(dc);
				const double hf = this->_scheme.FaceDerivative(f, up, gd);
				this->_scheme.FaceWeights(f, up, -1, wt);
				const double r  = this->WallFace(st, T::Dir(d)) ? 0.0
				                : dc / this->_NonDim->Viscosity(cc);
				const double h  = T::Spacing(this->_G, d/2);
				const double sg = (d % 2 == 0) ? -1.0 : +1.0;
				if(this->_NonDim->Mob().Slip())
				{
					//  The flux r (h_f + 2 Lambda) h_f: its derivative
					//  r ((2 h_f + 2 Lambda) gd - (h_f + 2 Lambda) wt) less
					//  the part A holds; the weight is then 1, not h_f
					const double L2 = 2.0*this->_NonDim->Mob().Lambda();
					a[d] = sg/h;
					for(int q = 0; q < 4; ++q)
						g[d*16 + T::FaceStencil(d, q)] = r*((2.0*hf + L2)*gd[q] - (hf + L2)*wt[q]);
					continue;
				}
				a[d] = sg*hf/h;
				for(int q = 0; q < 4; ++q) g[d*16 + T::FaceStencil(d, q)] = r*(2.0*gd[q] - wt[q]);
			}
			detail::SolutoG<SolutoCapillarityKernel<T> > gv; gv.self = this; gv.g = g;
			this->ScatterJacobiWide(crsJ, irow, st, gv, a, Coef);
		}
#if defined(_OPENMP) || defined(OMP)
	}
#endif
}

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_SOLUTO_CAPILLARITY_KERNEL_INCLUDED_
