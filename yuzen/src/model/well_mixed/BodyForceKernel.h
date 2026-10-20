/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  BodyForceKernel.h : the dimension-independent implementation of the body
*                      force terms (centrifugal force and in-plane gravity)
*
*  The centrifugal force and the in-plane gravity share one body of code,
*  for 1-D and 2-D alike; ParallelGravityKernel.h and CentrifugalKernel.h
*  supply what differs. Everything around the four-point upwind
*  interpolation lives in UpwindKernelBase.h.
*
*  Unlike the other physics terms, this one is assembled "one axis at a
*  time". For axis a, the contributions of the negative-side face
*  (direction 2a) and the positive-side face (direction 2a+1) are combined
*  into a single expression and added to the five stencil points of that
*  axis. In 2-D the five X points are added first and the five Y points
*  after, so the diagonal is written into twice. Do not change the shape or
*  the order of these additions (the rounding would change).
*
*  What a derived class decides:
*    - the driving field Drive / the driving field on a face DriveFace /
*      the weight Weight
*    - the per-axis coefficient Mag(a)
*    - whether the 1-D case adds the outflow in Y to the diagonal
*      (SideOutflow)
*    - where the velocity is written and with which sign
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _BODY_FORCE_KERNEL_INCLUDED_
#define _BODY_FORCE_KERNEL_INCLUDED_

#include "UpwindKernelBase.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

template<class T>
class BodyForceKernel : public UpwindKernelBase<T>
{
public:
	typedef typename T::Grid      Grid;
	typedef typename T::Site      Site;
	typedef typename T::Field     Field;
	typedef typename T::tEdge     tEdge;
	typedef typename T::tAdjacent tAdjacent;
	typedef typename T::Boundary  Boundary;
	typedef typename T::VelField  VelField;

	/** @brief Select the four-point upwind interpolation scheme. */
	BodyForceKernel() : UpwindKernelBase<T>()
	{
		this->_scheme.Build(C_Interpolation_base4::KK);
		this->_enabled = true;
	}

	//  For the contract, see PhysicsTermBase.h. Build is provided by each
	//  derived class (in-plane gravity / centrifugal force); what is here is
	//  only the skeleton the two have in common.
	/**
	 * @brief Add the body force to the matrix A, one axis at a time.
	 *
	 * For each axis the contributions of the negative and positive faces are
	 * combined into a single expression and scattered onto that axis' five
	 * stencil points. In 2-D the diagonal is therefore written twice, once
	 * for X and once for Y.
	 */
	virtual void SetMatrixA(CRS_Matrix &crsA, const Field &field);
	/**
	 * @brief Add the body force to the Jacobian, again one axis at a time.
	 *
	 * The per-row direction ordering is the one laid out at the top of
	 * UpwindKernelBase.h.
	 */
	virtual void SetJacobiMatrix(CRS_Matrix &crsJ, const Field &field);
	/**
	 * @brief Write the velocity this force produces into the slot the derived
	 *        class names in VelocityField().
	 */
	virtual void CalVelocity(Field &field);

	/**
	 * @brief The sideways outflow (SideOutflow) of the film in the state
	 *        given, as the diagonal of A takes it: the sum over the rows that
	 *        are not boundary rows of Mag h^3/mu(c), a rate per unit cell
	 *        measure; *resin, if given, receives the same sum weighted by c
	 *        (the solute the liquid carries away). 0 when there is none.
	 */
	double SideOutflowRate(const CRS_Matrix &crs, const Field &field, double *resin) const
	{
		double q = 0, qr = 0;
		if(this->_enabled && SideOutflow())
		{
			for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
			{
				Site st = this->_G->template Of<Site>(ii);
				if( crs.IsBC(T::Index(this->_G, st)) ) continue;
				const double hoo = field.HH(st), coo = field.CC(st);
				const double r = this->_NonDim->Mob().Slip() ? Mag(0)*(this->_NonDim->Mob().D(hoo)/this->_NonDim->Viscosity(coo))
				                                : Mag(0)*(hoo*hoo/this->_NonDim->Viscosity(coo))*hoo;
				q  += r;
				qr += r*coo;
			}
		}
		if(resin) *resin = qr;
		return q;
	}

	/**
	 * @brief The flux through face j of st in the state field, as A holds
	 *        it: Mag w_j h_f^2/mu(c_f) times the face value of the film
	 *        (per mesh spacing), the face value fixed-coefficient or
	 *        limited as the scheme is.
	 */
	virtual double FaceFlux(Site &st, int j, const Field &field, const std::vector<Site> &svec)
	{
		if(!this->_enabled) return 0;
		const tEdge te = T::Dir(j);
		const double w = Weight(st, te);
		double hs = 0, D = 0;
		if(this->Limited())
		{
			double hf, wt[4];
			LimitedFace(st, svec, j, field, hf, D, wt, 0);
			for(int q = 0; q < 4; ++q) hs += wt[q] * field.HH(svec[T::FaceStencil(j, q)]);
		}
		else
		{
			//  D as NonlinearTerm forms it (the upwind side from Drive), from
			//  the stencil given
			if(this->WallFace(st, te)) return 0;
			const double pp = Drive(st, te);
			double hh = 0;
			for(int q = 0; q < 4; ++q)
				hh += field.HH(svec[T::FaceStencil(j, q)])
				    * ( this->_scheme.c1((C_Interpolation_base4::t_iShift)q)
				      + Math::Sign(pp) * this->_scheme.c2((C_Interpolation_base4::t_iShift)q) );
			const double cc = 0.5 * (field.CC(st) + field.CC(T::Adj(this->_G, st, te)));
			D = (this->_NonDim->Mob().Slip() ? this->_NonDim->Mob().Upwind(hh) : hh*hh) / this->_NonDim->Viscosity(cc);
			const double dv = DriveFace(st, te);
			for(int q = 0; q < 4; ++q)
			{
				const int idx = T::FaceStencil(j, q);
				hs += (this->C(j)[idx] + Math::Sign(dv)*this->C2(j)[idx]) * field.HH(svec[idx]);
			}
		}
		return Mag(j/2) * (w * hs * D);
	}
	/** @brief The sideways outflow Mag h^3/mu of the cell (1-D centrifugal term), as SideOutflowRate takes it. */
	virtual double SideRate(Site &st, const Field &field) const
	{
		if(!(this->_enabled && SideOutflow())) return 0;
		const double hoo = field.HH(st), coo = field.CC(st);
		if(this->_NonDim->Mob().Slip()) return Mag(0)*(this->_NonDim->Mob().D(hoo)/this->_NonDim->Viscosity(coo));
		return Mag(0)*(hoo*hoo/this->_NonDim->Viscosity(coo))*hoo;
	}

protected:
	//--- what a derived class decides ---
	/**
	 * @brief The driving field on the face in direction d.
	 *
	 * Its sign is what selects the upwind side in NonlinearTerm and
	 * NonlinearJacobi.
	 */
	virtual double Drive(Site &s, tEdge te) const = 0;
	/**
	 * @brief The same quantity, but measured at the face coordinate itself,
	 *        as the matrix assembly and the velocity computation need it.
	 *
	 * On a uniform grid this is the same value as Drive, but the two do not
	 * agree down to the rounding, so they are kept separate.
	 */
	virtual double DriveFace(Site &s, tEdge te) const = 0;
	/**
	 * @brief The weight applied to the stencil (DriveFace for the
	 *        centrifugal force, 1 for the in-plane gravity).
	 */
	virtual double Weight(Site &s, tEdge te) const { (void)s; (void)te; return 1.0; }
	/** @brief The coefficient of the given axis. */
	virtual double Mag(int axis) const = 0;
	/**
	 * @brief Whether the outflow in Y is added to the diagonal; only the 1-D
	 *        centrifugal force does this.
	 */
	virtual bool SideOutflow() const { return false; }

	/**
	 * @brief The slot of the velocity field this term writes its velocity into.
	 *
	 * Each term is given a velocity field of its own, for visualization.
	 */
	virtual VelField& VelocityField(Field &field) const = 0;
	/**
	 * @brief The sign used when writing into the velocity field; it is
	 *        opposite for the centrifugal force and the in-plane gravity.
	 */
	virtual double VelocitySign() const = 0;

	Vector<double> _cv2[2];					// upwind coefficients for the velocity
	/** @brief Read and write access to _cv2, one entry per axis. */
	Vector<double>&       CV2(int a)       { return _cv2[a]; }
	const Vector<double>& CV2(int a) const { return _cv2[a]; }

	/**
	 * @brief Divide the four-point upwind interpolation coefficients by the
	 *        mesh spacing, and lay out the velocity coefficients CV / CV2
	 *        for each axis at the same time.
	 */
	virtual void SetLinearCoef()
	{
		this->SetUpwindCoef(true);
		for(int a = 0; a < T::nAxis; ++a)
		{
			const double h = T::Spacing(this->_G, a);
			_cv2[a].resize(this->_n_stencil);
			for(int s = 0; s < this->_n_stencil; s++)
			{
				this->CV(a)[s] = -this->C(2*a)[s]  * h;
				_cv2[a][s]     = -this->C2(2*a)[s] * h;
			}
		}
	}
	/**
	 * @brief h²/μ(c) on a face.
	 *
	 * The film height is brought onto the face by the four-point
	 * interpolation, and the upwind side is chosen by the sign of Drive.
	 */
	virtual double NonlinearTerm(Site&, tEdge, const Field &field);
	/**
	 * @brief The above, differentiated with respect to one point (for the
	 *        Jacobian).
	 */
	virtual double NonlinearJacobi(Site&, tEdge, tAdjacent, const Field &field);

	//  A limited scheme (UpwindKernelBase.h). The direction of every face is
	//  that of DriveFace, for A, the Jacobian and the velocity alike.
	/**
	 * @brief The flux of face j per unit Mag: its h_f and h_f^2/mu, the
	 *        frozen weights wt (per mesh spacing) and, if g is given, the
	 *        exact d h_f / d h (per mesh spacing), on the four points of the face.
	 */
	void LimitedFace(Site &st, const std::vector<Site> &svec, int j, const Field &field,
	                 double &hf, double &D, double wt[4], double *g);
	void SetMatrixA_Limited(CRS_Matrix &crsA, const Field &field);
	void SetJacobiMatrix_Limited(CRS_Matrix &crsJ, const Field &field);
	void CalVelocity_Limited(Field &field);
};


//==================================================================
//  Implementation
//==================================================================

//------------------------------------------------------------------
//  h^2 / mu on a face
//------------------------------------------------------------------
template<class T>
double BodyForceKernel<T>::NonlinearTerm(Site &st, tEdge tedge, const Field &field)
{
	//  The viscosity is taken at the concentration of the face, the mean
	//  of the two cells that share it, as every other term takes it. With
	//  the concentration of the cell the row belongs to, the two rows on
	//  either side of a face would see different fluxes through it: the
	//  scheme would not be conservative, and with a viscosity that depends
	//  on the concentration it would not converge to the equation at all.
	double cc = T::IsFace(tedge)
	          ? 0.5 * (field.CC(st) + field.CC(T::Adj(this->_G, st, tedge)))
	          : field.CC(st);
	double hh = 0;

	//  No flux through the face of a wall (UpwindKernelBase::WallFace)
	if(T::IsFace(tedge) && this->WallFace(st, tedge)) return 0;

	if(T::IsFace(tedge))
	{
		const int d = T::DirIndex(tedge);
		const double pp = Drive(st, tedge);

		std::vector<Site> svec(this->_n_stencil);
		this->_G->SetSite(svec, st);

		for(int j = 0; j < 4; ++j)
		{
			const int idx = T::FaceStencil(d, j);
			hh += field.HH(svec[idx]) * ( this->_scheme.c1((C_Interpolation_base4::t_iShift)j)
			                            + Math::Sign(pp) * this->_scheme.c2((C_Interpolation_base4::t_iShift)j) );
		}
	}

	double val = this->_NonDim->Mob().Slip() ? this->_NonDim->Mob().Upwind(hh) : hh*hh;
	double Visco = this->_NonDim->Viscosity(cc);
	return val / Visco;
}


//------------------------------------------------------------------
//  The above, differentiated with respect to one point (tc)
//------------------------------------------------------------------
template<class T>
double BodyForceKernel<T>::NonlinearJacobi(Site &st, tEdge tedge, tAdjacent tc, const Field &field)
{
	//  The viscosity is taken at the concentration of the face, as in
	//  NonlinearTerm (see there).
	double cc = T::IsFace(tedge)
	          ? 0.5 * (field.CC(st) + field.CC(T::Adj(this->_G, st, tedge)))
	          : field.CC(st);
	double hh = 0, dhdhk = 0;

	if(T::IsFace(tedge) && this->WallFace(st, tedge)) return 0;

	if(T::IsFace(tedge))
	{
		const int d = T::DirIndex(tedge);
		const double pp = Drive(st, tedge);

		std::vector<Site> svec(this->_n_stencil);
		this->_G->SetSite(svec, st);

		for(int j = 0; j < 4; ++j)
		{
			const int idx = T::FaceStencil(d, j);
			const double c = this->_scheme.c1((C_Interpolation_base4::t_iShift)j)
			               + Math::Sign(pp) * this->_scheme.c2((C_Interpolation_base4::t_iShift)j);
			hh += field.HH(svec[idx]) * c;
			if(idx == (int)tc) dhdhk = c;
		}
	}

	double val = this->_NonDim->Mob().Slip() ? this->_NonDim->Mob().dUpwind(hh)*dhdhk : 2.0*hh*dhdhk;
	double Visco = this->_NonDim->Viscosity(cc);
	return val / Visco;
}


//------------------------------------------------------------------
//  Adding into the matrix A (one axis at a time)
//------------------------------------------------------------------
template<class T>
void BodyForceKernel<T>::SetMatrixA(CRS_Matrix &crsA, const Field &field)
{
	if(!this->_enabled) return;
	if(this->Limited()) { SetMatrixA_Limited(crsA, field); return; }

	const bool side = SideOutflow();
	int irow, a, j;
	double D[4], dv[4], w[4], hoo, coo, visco;

	Site st;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for private(st,irow,a,j,D,dv,w,hoo,coo,visco)
#endif
	for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
	{
		st = this->_G->template Of<Site>(ii);
		irow = T::Index(this->_G, st);
		if( crsA.IsBC(irow) ) continue;

		for(j = 0; j < T::nDir; ++j)
		{
			dv[j] = DriveFace(st, T::Dir(j));
			w[j]  = Weight   (st, T::Dir(j));
			D[j]  = NonlinearTerm(st, T::Dir(j), field);
		}

		for(a = 0; a < T::nAxis; ++a)
		{
			const int dm = 2*a, dp = 2*a + 1;
			const double Coef = Mag(a);

			//  The four non-center points first, then the diagonal.
			//  Do not change this order.
			static const int order[5] = {0, 1, 3, 4, 2};
			for(int q = 0; q < 5; ++q)
			{
				const int idx = T::AxisStencil(a, order[q]);
				const double val =
					-w[dm]*(this->C(dm)[idx] + Math::Sign(dv[dm])*this->C2(dm)[idx])*D[dm]
					+w[dp]*(this->C(dp)[idx] + Math::Sign(dv[dp])*this->C2(dp)[idx])*D[dp];

				if(T::IsDiag(idx)) crsA.Diag(irow) += Coef*val;
				else               crsA(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx))) += Coef*val;
			}
		}

		if(side)
		{	//  Count the outflow in Y even in 1-D
			hoo = field.HH(st);
			coo = field.CC(st);
			visco = this->_NonDim->Viscosity(coo);
			crsA.Diag(irow) += this->_NonDim->Mob().Slip() ? Mag(0)*(this->_NonDim->Mob().DperH(hoo)/visco) : Mag(0)*(hoo*hoo/visco);
		}
	}
}


//------------------------------------------------------------------
//  Adding into the Jacobian (one axis at a time)
//------------------------------------------------------------------
template<class T>
void BodyForceKernel<T>::SetJacobiMatrix(CRS_Matrix &crsJ, const Field &field)
{
	if(!this->_enabled) return;
	if(this->Limited()) { SetJacobiMatrix_Limited(crsJ, field); return; }

	const bool side = SideOutflow();
	int irow, s, a, j;
	double dv[4], w[4], acc[4], hoo, coo, visco;
	double g[4*16];

	Site st;
	//  Per-thread workspace. It is created outside the parallel region and
	//  handed out with firstprivate.
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,irow,s,a,j,dv,w,acc,g,hoo,coo,visco) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			irow = T::Index(this->_G, st);
			if( crsJ.IsBC(irow) ) continue;

			this->_G->SetSite(svec, st);

			for(j = 0; j < T::nDir; ++j)
			{
				dv[j] = DriveFace(st, T::Dir(j));
				w[j]  = Weight   (st, T::Dir(j));

				const double sg = (j % 2 == 0) ? -1.0 : +1.0;
				acc[j] = 0;
				for(s = 0; s < this->_n_stencil; s++)
					acc[j] += sg*w[j]*(this->C(j)[s] + Math::Sign(dv[j])*this->C2(j)[s]) * field.HH(svec[s]);

				for(int q = 0; q < 4; ++q)
				{
					const int idx = T::FaceStencil(j, q);
					g[j*16 + idx] = NonlinearJacobi(st, T::Dir(j), (tAdjacent)idx, field);
				}
			}

			for(a = 0; a < T::nAxis; ++a)
			{
				const int dm = 2*a, dp = 2*a + 1;
				const double Coef = Mag(a);

				for(int q = 0; q < 5; ++q)
				{
					const int idx = T::AxisStencil(a, q);
					double val = 0;

					// The direction that has this point as its neighbor
					// comes first; the rest follow in that axis' order
					const int first = T::DirOfAdj(idx);
					const bool inM = this->InFace(dm, idx), inP = this->InFace(dp, idx);
					if(first == dm && inM) val += g[dm*16 + idx]*acc[dm];
					if(first == dp && inP) val += g[dp*16 + idx]*acc[dp];
					if(first != dm && inM) val += g[dm*16 + idx]*acc[dm];
					if(first != dp && inP) val += g[dp*16 + idx]*acc[dp];

					if(T::IsDiag(idx)) crsJ.Diag(irow) += Coef*val;
					else               crsJ(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx))) += Coef*val;
				}
			}

			if(side)
			{	//  Count the outflow in Y even in 1-D
				hoo = field.HH(st);
				coo = field.CC(st);
				visco = this->_NonDim->Viscosity(coo);
				//  A holds Mag h^2/mu here, so A h = Mag h^3/mu and the
				//  Jacobian needs 2 Mag h^2/mu on top.
				crsJ.Diag(irow) += this->_NonDim->Mob().Slip() ? Mag(0)*((this->_NonDim->Mob().dD(hoo) - this->_NonDim->Mob().DperH(hoo))/visco)
				                                 : Mag(0)*(2.0*hoo*hoo/visco);
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
void BodyForceKernel<T>::CalVelocity(Field &field)
{
	if(this->Limited()) { CalVelocity_Limited(field); return; }
	int s, a;
	double D, Q, h, dv, wv;

	Site st;
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,s,a,D,Q,h,dv,wv) firstprivate(svec)
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
				const int dm = 2*a;
				D  = NonlinearTerm(st, T::Dir(dm), field);
				dv = DriveFace(st, T::Dir(dm));
				wv = Weight   (st, T::Dir(dm));

				Q = 0; h = 0;
				for(s = 0; s < this->_n_stencil; s++)
				{
					if(T::SkipInVel(s)) continue;
					Q += wv*(this->CV(a)[s] + Math::Sign(dv)*this->CV2(a)[s]) * field.HH(svec[s]);
					h += (this->C(dm)[s] + Math::Sign(dv)*this->C2(dm)[s]) * field.HH(svec[s]);
				}
				//  C holds the interpolation coefficients divided by the mesh
				//  spacing (SetLinearCoef), so the sum above is h on the face
				//  over dx. Without this the velocity came out multiplied by
				//  the dimensionless dx -- 0.6 % of rho g h^2/(3 mu) on the
				//  grid of allphys_1d_reg -- and it is this velocity that
				//  carries the concentration.
				h *= T::Spacing(this->_G, a);
				Q *= Mag(a) * D;

				Q = T::AxisMask(this->_Ny, a, Q);
				//  A wall face has no flux, and its h may be zero (Pinning)
				if(this->WallFace(st, T::Dir(dm))) { Q = 0; h = 1; }
				T::WriteFaceVel(this->_G, VelocityField(field), st, a,
				                VelocitySign() * Q / h);
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif
}

//==================================================================
//  The limited schemes
//
//  The flux through face j is F_j = Mag w_j h_f^3 / mu(c_f), h_f the
//  limited face value. A holds it as w_j (sum_k wt_k h_k) h_f^2/mu, wt
//  the weights frozen at the film (FaceWeights, downwind form), so that
//  A(h) h = F(h) to rounding; the Jacobian adds w_j (3 g_k - wt_k)
//  h_f^2/mu, g_k = d h_f / d h_k (FaceDerivative), which makes J the
//  exact derivative of A(h) h wherever the limiter is differentiable.
//==================================================================
template<class T>
void BodyForceKernel<T>::LimitedFace(Site &st, const std::vector<Site> &svec, int j, const Field &field,
                                     double &hf, double &D, double wt[4], double *g)
{
	const tEdge te = T::Dir(j);
	const double dv = DriveFace(st, te);
	const double hs = T::Spacing(this->_G, j/2);
	double f[4];
	this->FaceHeights(svec, j, field, f);
	if(g)
	{
		hf = this->_scheme.FaceDerivative(f, dv, g);
		for(int q = 0; q < 4; ++q) g[q] /= hs;
	}
	else
		hf = this->_scheme.FaceValue(f, dv);
	this->_scheme.FaceWeights(f, dv, -1, wt);
	for(int q = 0; q < 4; ++q) wt[q] /= hs;
	const double cc = 0.5 * (field.CC(st) + field.CC(T::Adj(this->_G, st, te)));
	D = (this->_NonDim->Mob().Slip() ? this->_NonDim->Mob().Upwind(hf) : hf*hf) / this->_NonDim->Viscosity(cc);
	//  No flux through the face of a wall: D = 0 takes it out of A, the
	//  Jacobian (whose every entry of the face is a multiple of D) and the
	//  velocity
	if(this->WallFace(st, te)) D = 0;
}

template<class T>
void BodyForceKernel<T>::SetMatrixA_Limited(CRS_Matrix &crsA, const Field &field)
{
	const bool side = SideOutflow();
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

			double W[4][16], D[4], w[4];
			for(int j = 0; j < T::nDir; ++j)
			{
				for(int k = 0; k < 16; ++k) W[j][k] = 0;
				double hf, wt[4];
				LimitedFace(st, svec, j, field, hf, D[j], wt, 0);
				w[j] = Weight(st, T::Dir(j));
				for(int q = 0; q < 4; ++q) W[j][T::FaceStencil(j, q)] = wt[q];
			}
			for(int a = 0; a < T::nAxis; ++a)
			{
				const int dm = 2*a, dp = 2*a + 1;
				const double Coef = Mag(a);
				for(int q = 0; q < 5; ++q)
				{
					const int idx = T::AxisStencil(a, q);
					const double val = -w[dm]*W[dm][idx]*D[dm] + w[dp]*W[dp][idx]*D[dp];
					if(T::IsDiag(idx)) crsA.Diag(irow) += Coef*val;
					else               crsA(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx))) += Coef*val;
				}
			}
			if(side)
			{
				const double hoo = field.HH(st);
				crsA.Diag(irow) += this->_NonDim->Mob().Slip() ? Mag(0)*(this->_NonDim->Mob().DperH(hoo)/this->_NonDim->Viscosity(field.CC(st)))
				                                : Mag(0)*(hoo*hoo/this->_NonDim->Viscosity(field.CC(st)));
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}
#endif
}

template<class T>
void BodyForceKernel<T>::SetJacobiMatrix_Limited(CRS_Matrix &crsJ, const Field &field)
{
	const bool side = SideOutflow();
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

			double X[4][16];
			for(int j = 0; j < T::nDir; ++j)
			{
				for(int k = 0; k < 16; ++k) X[j][k] = 0;
				double hf, D, wt[4], g[4];
				LimitedFace(st, svec, j, field, hf, D, wt, g);
				//  d(D(h_f) h_f)/dh_f over D(h_f): 3 without slip
				const bool slip = this->_NonDim->Mob().Slip();
				const double r3 = !slip ? 3.0 : (this->_NonDim->Mob().Upwind(hf) != 0.0
				                                  ? this->_NonDim->Mob().dD(hf) / this->_NonDim->Mob().Upwind(hf) : 0.0);
				const double sg = (j % 2 == 0) ? -1.0 : +1.0;
				const double w  = Weight(st, T::Dir(j));
				for(int q = 0; q < 4; ++q)
					X[j][T::FaceStencil(j, q)] = slip ? sg*w*(r3*g[q] - wt[q])*D : sg*w*(3.0*g[q] - wt[q])*D;
			}
			for(int a = 0; a < T::nAxis; ++a)
			{
				const int dm = 2*a, dp = 2*a + 1;
				const double Coef = Mag(a);
				for(int q = 0; q < 5; ++q)
				{
					const int idx = T::AxisStencil(a, q);
					const double val = X[dm][idx] + X[dp][idx];
					if(T::IsDiag(idx)) crsJ.Diag(irow) += Coef*val;
					else               crsJ(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx))) += Coef*val;
				}
			}
			if(side)
			{
				const double hoo = field.HH(st);
				crsJ.Diag(irow) += this->_NonDim->Mob().Slip() ? Mag(0)*((this->_NonDim->Mob().dD(hoo) - this->_NonDim->Mob().DperH(hoo))/this->_NonDim->Viscosity(field.CC(st)))
				                                : Mag(0)*(2.0*hoo*hoo/this->_NonDim->Viscosity(field.CC(st)));
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}
#endif
}

template<class T>
void BodyForceKernel<T>::CalVelocity_Limited(Field &field)
{
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
			if(!T::InsideForVelocity(this->_G, this->_Nx, this->_Ny, st)) continue;
			this->_G->SetSite(svec, st);
			for(int a = 0; a < T::nAxis; ++a)
			{
				const int dm = 2*a;
				double hf, D, wt[4];
				LimitedFace(st, svec, dm, field, hf, D, wt, 0);
				//  The flux over h_f: Mag w h_f^2/mu, pointing the way the
				//  fixed-coefficient CalVelocity writes it
				double Q = -Weight(st, T::Dir(dm)) * Mag(a) * D;
				Q = T::AxisMask(this->_Ny, a, Q);
				T::WriteFaceVel(this->_G, VelocityField(field), st, a, VelocitySign() * Q);
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}
#endif
}

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_BODY_FORCE_KERNEL_INCLUDED_
