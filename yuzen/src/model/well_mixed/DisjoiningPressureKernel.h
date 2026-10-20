/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  DisjoiningPressureKernel.h : the disjoining pressure term. This is the
*                               term that gives the solver wettability, and
*                               hence contact lines.
*
*  As the film gets thinner, the interaction between the substrate and the
*  interface starts to take effect. The disjoining pressure expresses that
*  interaction as a pressure; here it takes a two-term power law.
*
*      Π = B [ (h* / h)^n - (h* / h)^m ]        (n > m > 1)
*      B = (n-1)(m-1) / (h* (n-m)) * σ(1-cosθe)
*        ≈ (n-1)(m-1) / (2 h* (n-m)) * σ θe²
*
*  The first term is repulsive (it tries to thicken the film) and the second
*  is attractive. Π vanishes at h = h*, and that film thickness is what
*  remains as the precursor film. The coefficient B is chosen so that the
*  equilibrium contact angle matches the θe given in the input.
*
*  Once non-dimensionalized, this term enters the film height equation in
*  the same form as the other pressure terms.
*
*      ∂H/∂T = -∂/∂X [ ... + Dp H³/(3μ̂) ∂Π̂/∂X ]
*      Π̂ = (H* / H)^n - (H* / H)^m
*      Dp = (L/h0)² (n-1)(m-1) / (2 H* (n-m)) θe²
*
*  ** Key points of the implementation **
*
*  1. The structure is exactly that of the vertical gravity term: the face
*     mobility h³/μ -- the positive one of the capillary term,
*     Mobility::Face -- multiplied by the gradient of something. For gravity
*     that something is H; for the disjoining pressure it is Π̂(H). The
*     term can therefore be written on the same skeleton as
*     VerticalGravityKernel.
*
*  2. Yuzen's framework is ∂φ/∂t = -Aφ - B, and only what is put into
*     A receives the Crank-Nicolson weighting. Π̂ is a nonlinear function
*     of H, but seeing it as Π̂ = (Π̂/H)·H lets it go into A, and that is
*     what is done here. Putting it into B (_VecB) would make this one term
*     backward Euler.
*
*  3. The Jacobian holds **both** the derivative of the mobility and the
*     derivative of the gradient. The existing gravity term keeps only the
*     derivative of the mobility (it drops the gradient side), but the
*     disjoining pressure is very stiff, Π ~ h^-n, and the gradient side
*     dominates, so it cannot be dropped.
*
*  4. The mobility is divided by the concentration-dependent viscosity
*     μ(c), just as in the other terms. This term is therefore coupled to
*     the concentration field. Near a contact line the disjoining pressure
*     is the largest term, so ∂F_h/∂c is largest there as well.
*
*  ** The order of summation **
*  This follows the same rule as the existing kernels.
*    - each row of SetMatrixA        : axis order    (Xm, Xp, Ym, Yp) -> T::Dir(d)
*    - the diagonal sums of SetJacobiMatrix : stencil order (Ym, Xm, Xp, Yp) -> T::StencilDir(k)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _DISJOINING_PRESSURE_KERNEL_INCLUDED_
#define _DISJOINING_PRESSURE_KERNEL_INCLUDED_

#include <cmath>

#include "PhysicsTermBase.h"
#include "NonDimensionalNumbers.h"
#include "LubFields.h"
#include "MyMath.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

template<class T>
class DisjoiningPressureKernel : public PhysicsTermBase<T>
{
public:
	typedef typename T::Grid      Grid;
	typedef typename T::Site      Site;
	typedef typename T::Field     Field;
	typedef typename T::tEdge     tEdge;
	typedef typename T::tAdjacent tAdjacent;
	typedef typename T::Boundary  Boundary;

	/**
	 * @brief Start from the (3, 2) pair of exponents and a precursor film of
	 *        0.01; Build() replaces all three from the input.
	 */
	DisjoiningPressureKernel() : PhysicsTermBase<T>(), _n(3.0), _m(2.0), _hstar(0.01) {}

	//  The disjoining pressure. The coefficient Dp is fixed by the
	//  equilibrium contact angle θe. The coefficient _Magnitude (sign
	//  included) and the finite-difference coefficients are decided here.
	//  For the contract, see PhysicsTermBase.h.
	/**
	 * @brief Switch the disjoining pressure term on or off, read its
	 *        exponents and precursor thickness, and build the coefficients.
	 * @param inp_phys  the Physics block of the input file; the exponents
	 *                  n and m are read from here
	 * @param g         grid the term is discretized on
	 * @param nonDim    the set of non-dimensional numbers; Dp and the
	 *                  precursor film thickness H* are taken from here
	 * @param bc        boundary condition, used to fold the ghost points into the stencil
	 *
	 * The overall coefficient _Magnitude is Dp/3. Dp is calibrated so that
	 * the equilibrium contact angle matches the one given in the input.
	 */
	virtual void Build(Input::I_Physics &inp_phys, Grid *g, NonDimensionalNumbers *nonDim, Boundary *bc)
	{
		this->_G = g;
		this->_NonDim = nonDim;
		this->_boundary = bc;

		this->_enabled = inp_phys.DisjoiningPressure.Enable;

		_n     = inp_phys.DisjoiningPressure.ExponentN;
		_m     = inp_phys.DisjoiningPressure.ExponentM;
		_hstar = nonDim->HStar();

		//  The sign follows the same arrangement as the gravity term. For
		//  the governing equation
		//    ∂H/∂T = -∂/∂X[ ... + Dp H³/(3μ̂) ∂Π̂/∂X ]
		//  the contribution to A is
		//  +Dp/3 × Σ_faces mobility × (Π̂_neighbor - Π̂_self)/dx².
		this->_Magnitude = this->_NonDim->Dp() / 3.0;

		T::SetExtent(this->_G, this->_Nx, this->_Ny);
		SetLinearCoef();
	}

	/**
	 * @brief Add the disjoining pressure to the matrix A.
	 *
	 * Pi is a nonlinear function of h, but writing it as (Pi/h)*h lets the
	 * term go into A and so receive the Crank-Nicolson weighting.
	 */
	virtual void SetMatrixA(CRS_Matrix &crsA, const Field &field);
	/**
	 * @brief Add the disjoining pressure to the Jacobian.
	 *
	 * Both the derivative of the mobility and the derivative of the gradient
	 * are included. The gradient side dominates here and dropping it, as the
	 * gravity term does, stops the Newton iteration from converging.
	 */
	virtual void SetJacobiMatrix(CRS_Matrix &crsJ, const Field &field);
	/**
	 * @brief Write the velocity this term produces into
	 *        field.UU().U_DisjoiningPressure().
	 */
	virtual void CalVelocity(Field &field);

	/**
	 * @brief The flux through face d of s in the state field: Coef times
	 *        the face mobility times the difference of Pi across the face,
	 *        what the row's A h takes of that face.
	 */
	virtual double FaceFlux(Site &st, int d, const Field &field, const std::vector<Site> &svec)
	{
		(void)svec;
		if(!this->_enabled) return 0;
		const Site a = T::Adj(this->_G, st, T::Dir(d));
		const double v = this->C(d)[T::AdjIdx(d)] * Pi(field.HH(a)) + this->C(d)[T::CenterIdx()] * Pi(field.HH(st));
		const double f = this->_Magnitude * (NonlinearTerm(st, T::Dir(d), field) * v);
		return (d % 2 == 0) ? -f : f;
	}

protected:
	double _n, _m;		// the exponents
	double _hstar;		// the precursor film thickness (non-dimensional)

	/**
	 * @brief The three-point second-derivative coefficients, which are fixed
	 *        by the mesh spacing alone. They are the same as the gravity
	 *        term's.
	 */
	virtual void SetLinearCoef();
	virtual double NonlinearTerm(Site&, tEdge, const Field &field);			// the face mobility h³/μ
	virtual double NonlinearJacobi(Site&, tEdge, tAdjacent, const Field &field);	// the above, differentiated at one point

	//--------------------------------------------------------------
	//  The disjoining pressure itself, and its derivative with respect to
	//  the film height
	//--------------------------------------------------------------
	/**
	 * @brief Π̂(H) = (H* / H)^n - (H* / H)^m.
	 * @return 0 when H <= 0
	 */
	double Pi(double h) const
	{
		if(h <= 0) return 0;
		const double r = _hstar / h;
		return std::pow(r, _n) - std::pow(r, _m);
	}

	/**
	 * @brief dPi/dH = -(1/H) [ n (H* / H)^n - m (H* / H)^m ].
	 * @return 0 when H <= 0
	 */
	double dPi(double h) const
	{
		if(h <= 0) return 0;
		const double r = _hstar / h;
		return -( _n*std::pow(r, _n) - _m*std::pow(r, _m) ) / h;
	}

	/**
	 * @brief The "disjoining pressure per unit H", which is what lets Π̂ be
	 *        written into the matrix.
	 * @return 0 when H <= 0
	 */
	double PiPerH(double h) const
	{
		if(h <= 0) return 0;
		return Pi(h) / h;
	}
};


//==================================================================
//  Implementation
//==================================================================

//  The three-point second-derivative coefficients. The same as the gravity
//  term's.
template<class T>
void DisjoiningPressureKernel<T>::SetLinearCoef()
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
//  h³/μ on a face: the positive face mobility of the capillary term
//  (Mobility::Face, Diez and Kondic), the same for the vertical gravity.
//  With the mean (h_o³ + h_a³)/2 of the two cells, as this term once took,
//  the flux of a face was not one mobility times the difference of the
//  potential sigma lap h + Dp Pi - Bo h: the dissipation of the scheme
//  was not a sum of squares, its free energy could rise near a steady
//  film (tools/analytic.py energy_saturate_mono), and a mobility of h³/2
//  survived against a dry neighbor.
//------------------------------------------------------------------
template<class T>
double DisjoiningPressureKernel<T>::NonlinearTerm(Site &s, tEdge tedge, const Field &field)
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


template<class T>
double DisjoiningPressureKernel<T>::NonlinearJacobi(Site &s, tEdge tedge, tAdjacent tc, const Field &field)
{
	double Coo = field.CC(s);
	Site a = T::IsFace(tedge) ? T::Adj(this->_G, s, tedge) : s;

	//  The derivative at the point tc: the cell itself (center) or the
	//  neighbor across the face
	const double ho = field.HH(s), ha = field.HH(a);
	double Caj = field.CC(a);
	double Cfc = 0.5*(Coo + Caj);

	double Visco = this->_NonDim->Viscosity(Cfc);
	const Mobility &mob = this->_NonDim->Mob();
	return (((int)tc == T::CenterIdx()) ? mob.Face_d1(ho, ha) : mob.Face_d1(ha, ho)) / Visco;
}


//------------------------------------------------------------------
//  Adding into the matrix A
//
//    (A h)_i = Coef Σ_d mobility_d Σ_s C(d)[s] Π̂(h_s)
//  is turned into a coefficient of h by seeing Π̂(h_s) as
//  (Π̂(h_s)/h_s) · h_s.
//------------------------------------------------------------------
template<class T>
void DisjoiningPressureKernel<T>::SetMatrixA(CRS_Matrix &crsA, const Field &field)
{
	if(!this->_enabled) return;

	const double Coef = this->_Magnitude;
	int irow, d, k;
	double G[4];

	Site st;
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,irow,d,k,G) firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = this->_G->imin(); ii <= this->_G->imax(); ii++)
		{
			st = this->_G->template Of<Site>(ii);
			irow = T::Index(this->_G, st);
			if( crsA.IsBC(irow) ) continue;

			for(d = 0; d < T::nDir; ++d) G[d] = NonlinearTerm(st, T::Dir(d), field);

			this->_G->SetSite(svec, st);

			for(k = 0; k < T::nNarrow; ++k)
			{
				const int idx = T::NarrowIdx(k);
				const double php = PiPerH(field.HH(svec[idx]));

				double val = 0;
				// The summation runs in axis order
				for(d = 0; d < T::nDir; ++d) val += this->C(d)[idx] * G[d];
				val *= php;

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
//
//    F_i = Coef Σ_d G_d A_d,      A_d = Σ_s C(d)[s] Π̂(h_s)
//    ∂F_i/∂h_j = Coef Σ_d [ (∂G_d/∂h_j) A_d + G_d C(d)[j] Π̂'(h_j) ]
//
//  The Jacobian starts as a copy of A (HeightEquationBase::SetMatrixJacobi),
//  and A already holds Coef G_d C(d)[j] Pi(h_j)/h_j, the form that lets Pi
//  be written into the matrix. What is added here is therefore the
//  derivative less that: G_d C(d)[j] (Pi'(h_j) - Pi(h_j)/h_j). The whole of
//  Pi' would count Pi/h twice; the Jacobian would be 50 % off, and the
//  Newton iteration on a contact line converge the slower for it.
//------------------------------------------------------------------
template<class T>
void DisjoiningPressureKernel<T>::SetJacobiMatrix(CRS_Matrix &crsJ, const Field &field)
{
	if(!this->_enabled) return;

	const double Coef = this->_Magnitude;
	int irow, s, d, k;
	double A[4], G[4], gAdj[4], gCen[4], gaSum, gcSum;

	Site st;
	std::vector<Site> svec(this->_n_stencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel private(st,irow,s,d,k,A,G,gAdj,gCen,gaSum,gcSum) firstprivate(svec)
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
				A[d] = 0;
				for(s = 0; s < this->_n_stencil; s++)
					A[d] += this->C(d)[s] * Pi(field.HH(svec[s]));

				G[d]    = NonlinearTerm  (st, T::Dir(d), field);
				gAdj[d] = NonlinearJacobi(st, T::Dir(d), (tAdjacent)T::AdjIdx(d),   field);
				gCen[d] = NonlinearJacobi(st, T::Dir(d), (tAdjacent)T::CenterIdx(), field);
			}

			// The two summations that enter the diagonal run in stencil order
			gaSum = 0; gcSum = 0;
			for(k = 0; k < T::nDir; ++k)
			{
				const int dd = T::StencilDir(k);
				gaSum += gCen[dd] * A[dd];
				gcSum += G[dd] * this->C(dd)[T::CenterIdx()];
			}

			for(k = 0; k < T::nNarrow; ++k)
			{
				const int idx = T::NarrowIdx(k);
				const double hj  = field.HH(svec[idx]);
				const double dpi = dPi(hj) - PiPerH(hj);

				if(T::IsDiag(idx))
				{
					crsJ.Diag(irow) += Coef*(gaSum + gcSum*dpi);
				}
				else
				{
					//  Every off-diagonal point of the narrow stencil is the
					//  neighbor of some direction, so DirOfAdj never
					//  returns a negative value here; the compiler cannot
					//  see that, so the test is written out explicitly.
					const int dd = T::DirOfAdj(idx);
					if(dd < 0) continue;
					crsJ(irow, T::Index(this->_G, T::StencilSite(this->_G, st, idx)))
						+= Coef*( gAdj[dd]*A[dd] + G[dd]*this->C(dd)[idx]*dpi );
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
void DisjoiningPressureKernel<T>::CalVelocity(Field &field)
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
					Q += this->CV(a)[s] * Pi(field.HH(svec[s]));
				}
				Q *= this->_Magnitude * Gf;

				Q = T::AxisMask(this->_Ny, a, Q);
				T::WriteFaceVel(this->_G, field.UU().U_DisjoiningPressure(), st, a,
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

#endif //_DISJOINING_PRESSURE_KERNEL_INCLUDED_
