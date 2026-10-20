/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ConcEquationBase.h : the concentration equation of the well-mixed model,
*                       for 1-D and 2-D alike (the grid traits T of
*                       GridTraits.h supply what differs)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _WM_CONC_EQUATION_BASE_INCLUDED_
#define _WM_CONC_EQUATION_BASE_INCLUDED_

#include "Opt.h"
#include <cmath>
#include "Profile.h"
#include <vector>

#include "SemiImplicit.h"
#include "FiniteDifferences.h"
#include "TimeSteps.h"
#include "LubFields.h"
#include "ManufacturedSource.h"
#include "BoundaryCondition.h"
#include "NonDimensionalNumbers.h"
#include "GridTraits.h"
#include "ConcStep.h"
#include "YuzenUtil.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class C_Fluid;

//--------------------------------------------------------------------------
//  The concentration equation of the well-mixed model.
//
//  It is the finite-volume balance of the solute h c (ConcStep.h): each
//  face carries c_f times the volume flux the film-height equation puts
//  through that face in the same step, at the same Crank-Nicolson
//  weighting, minus the diffusive flux D(c_f) h_f (c_hi - c_lo) / dx; the
//  nozzle supplies liquid at the concentration it lands on, the 1-D
//  sideways outflow takes it away at the concentration of the cell, and
//  evaporation takes solvent only. The unknown is the change d = c1 - c0
//  of the step. Subtracting c0 times the film-height step from the
//  balance gives
//
//      h1 d + theta dt [sum_f F_f^1 d_f - K_f^1 (d_hi - d_lo)
//                       + (r1 - S) d]
//        = -dt sum_f [theta F_f^1 (c0_f - c0) + (1 - theta) F_f^0 (c0_f - c0)]
//          + dt sum_f [theta K_f^1 + (1 - theta) K_f^0] (c0_hi - c0_lo)
//          + dt E c0 (+ dt h S_c, the manufactured source),
//
//  every face value c0_f - c0 formed as sum_k w_k (c0_k - c0) from the
//  weights of the face (which sum to 1): a uniform c0 gives a right-hand
//  side of exactly zero without evaporation, and d = 0. On every Newton
//  pass of the film the equation is solved with h1 the film of the pass
//  (the guess) and the fluxes there; the step is completed with the film
//  the iteration ended on and the fluxes of the step it took
//  (HeightEquationBase::AfterNKLoop, ConsistentFluxes), and the solute is
//  then formed from the face fluxes once more (Reconstruct), so that the
//  sum of h c over the cells changes by the fluxes through the edges and
//  the sources alone, to rounding, whatever the tolerance of the solve.
//
//  The face value is the scheme's (Numerics.Convection.Concentration),
//  upwinded by the sign of the face's volume flux. A limited scheme
//  freezes its limiter at the latest iterate (LinkLimiterState), in the
//  form of FaceWeights that keeps the frozen matrix an M-matrix.
//--------------------------------------------------------------------------
template<class T>
class WellMixed_ConcEquationBase : public C_SemiImplicitEquation
{
public:
	typedef typename T::Grid     Grid;
	typedef typename T::Site     Site;
	typedef typename T::Field    Field;
	typedef typename T::Boundary Boundary;

	WellMixed_ConcEquationBase()
		: C_SemiImplicitEquation(), _limiterC(0), _step(0), _G(0), _field(0), _prop(0),
		  _NonDim(0), _boundary(0), _swConvection(true), _swDiffusion(true), _edgeOut(0) {}
	virtual ~WellMixed_ConcEquationBase() {}

	/**
	 * @brief Assemble the equation from the grid, properties, field and boundary condition.
	 * @param inp     the whole input file; the blocks this equation needs are read from it
	 * @param prop    fluid properties
	 * @param ts      time-step object supplying dt and the current time
	 * @param ff      field the equation reads its coefficients from
	 * @param nonDim  non-dimensional numbers (diffusivity)
	 * @param bc      boundary condition applied to the concentration rows
	 */
	void Build(Input &inp, C_Fluid *prop, C_TimeSteps *ts, Field *ff,
	           NonDimensionalNumbers *nonDim, Boundary *bc);

	/** @brief The face interpolation of the advection. */
	const C_Interpolation_base4& Scheme() const {return _scheme;}
	/** @brief The time weighting of the equation (0.5: Crank-Nicolson). */
	double Theta() const {return _theta;}
	/**
	 * @brief The concentration a limited scheme freezes its limiter at, per
	 *        row: the latest iterate as the equation solved it (the field
	 *        holds the same values except where c > 1, which the film-height
	 *        equation caps there for the coefficients). Without it the
	 *        field is taken.
	 */
	void LinkLimiterState(const Vector<double> *c) {_limiterC = c;}
	/** @brief The film side of the solves (ConcStep.h). */
	void LinkStep(const WMConcStep *s) {_step = s;}
	/** @brief Start the solves of a step from d = 0. */
	void ResetGuess() { if(_dGuess.size() == _nRows) _dGuess = 0; }
	/** @brief Whether the advection is switched on. */
	bool ConvectionOn() const {return _swConvection;}
	/**
	 * @brief The solute that left through the open edges (FreshInlet,
	 *        Outlet) in the committed solve last made, as a rate per unit
	 *        cell measure: the theta-weighted solute flux of every edge face
	 *        of an open side, outward, as the balance of the step took it
	 *        (Reconstruct). 0 when no committed solve was made.
	 */
	double EdgeOutflowRate() const {return _edgeOut;}

protected:
	virtual void SetFlagBC();
	/** @brief The boundary rows, and the structure of the stencil of the others. */
	virtual void SetMatrixA();
	/** @brief One solve of the step, for d, then c1 = c0 + d (see the class). */
	virtual int Solve(int &iter, double &resid);
	virtual double CheckTimeStepCriterion() { return LocalTruncationCriterion(_phi_new); }

	/**
	 * @brief After a committed solve: c1 from the balance of h c with the
	 *        solute fluxes formed once per face from c1 itself.
	 */
	void Reconstruct();
	/**
	 * @brief The solute through the edge faces of the open sides, outward,
	 *        of the step as solved (_edgeOut): the face fluxes of Reconstruct
	 *        (_G2), or, without the advection, the cell's own c times the
	 *        volume flux and the diffusive flux, as the rows took them.
	 */
	void EdgeOutflow();
	/**
	 * @brief The theta-weighted solute flux through the face between row lo
	 *        and its neighbor on the high side along axis a, seen from the
	 *        real cell st, whose face d it is.
	 */
	double SoluteFlux(Site &st, int d, const std::vector<Site> &svec) const;
	/** @brief D(c_f) h_f / dn^2 of a face, from the two cells' film and concentration. */
	double KFace(double ha, double hb, double ca, double cb, double dn) const
	{
		if(!_swDiffusion) return 0;
		return _NonDim->Df() * _NonDim->Diffusivity(0.5 * (ca + cb)) * (0.5 * (ha + hb)) / (dn * dn);
	}
	/** @brief The concentration a limited scheme is frozen at, in row i. */
	double LimiterC(int i, const Site &s) const
	{
		return (_limiterC && _limiterC->size() == _nRows) ? (*_limiterC)[i] : _field->CC(const_cast<Site&>(s));
	}
	int Index(const Site &s) const { return T::Index(_G, s); }

	const Vector<double> *_limiterC;
	const WMConcStep *_step;

	Grid *_G;
	Field *_field;
	Field _field_ini;
	C_Fluid *_prop;
	NonDimensionalNumbers *_NonDim;
	Boundary *_boundary;
	//  The manufactured-solution source of the concentration (verification only)
	C_ManufacturedSource _mms;

	bool _swConvection;
	bool _swDiffusion;

	C_Interpolation_base4 _scheme;		// the interpolation scheme

	//  The change of the step, its right-hand side, and the guess the next
	//  solve of the step starts from
	Vector<double> _d, _rhsD, _dGuess;
	//  The theta-weighted solute flux of each face (Reconstruct), as F of
	//  ConcStep.h is indexed
	std::vector<double> _G2[T::nAxis];
	//  The solute through the open edges of the committed solve (EdgeOutflow)
	double _edgeOut;

private:
	WellMixed_ConcEquationBase(const WellMixed_ConcEquationBase&);
	void operator=(const WellMixed_ConcEquationBase&);
};


//==================================================================
//  Implementation
//==================================================================

template<class T>
void WellMixed_ConcEquationBase<T>::Build(Input &inp, C_Fluid *prop, C_TimeSteps *ts, Field *ff,
                                          NonDimensionalNumbers *nonDim, Boundary *bc)
{
	_timeStep = ts;
	_field = ff;
	_prop = prop;
	_boundary = bc;
	_NonDim = nonDim;
	_G = ff->GetGridPtr();

	//  The face interpolation of the advection (Numerics.Convection.Concentration)
	_scheme.Build(C_Interpolation_base4::FromKeyOr(inp.Numerics.Convection.Concentration,
	                                               C_Interpolation_base4::Upwind_3rd));

	_swConvection = inp.Physics.ConcField.EnableConvection;
	_swDiffusion  = inp.Physics.ConcField.EnableDiffusion;

	_field_ini = *_field;
	_mms.Build(inp, *nonDim, T::nAxis == 2);

	//  Build the equation solver. This calls SetMatrixA, so it comes last.
	C_ImplicitSolver::Build(ff->CC(), *ts, C_ImplicitSolver::CrankNicholson, false);
	_VecB = 0;

	//  The preconditioner, the scaling and the convergence criteria
	this->_bcg.Configure(inp.Numerics.ConcEquation.LinearSolver, (int)T::nAxis);
	this->_bcg.SetConcRole(true);
}

template<class T>
void WellMixed_ConcEquationBase<T>::SetFlagBC()
{
	for(Site s = _G->template First<Site>(); s != _G->template Last<Site>(); ++s)
		if(_G->CheckGhost(s) != Grid::Inside) this->SetAsBoundary(Index(s));
}

//------------------------------------------------------------------
//  The boundary rows (written back on every call: Apply sets them and
//  their right-hand side b in _Rhs), and zeros on the stencil of every
//  other row, which fixes the structure Solve fills
//------------------------------------------------------------------
template<class T>
void WellMixed_ConcEquationBase<T>::SetMatrixA()
{
	_boundary->Apply(_crsA, _Rhs, BoundaryCondition::Conc, _field_ini);
	for(int ii = _G->imin(); ii <= _G->imax(); ii++)
	{
		Site st = _G->template Of<Site>(ii);
		const int irow = Index(st);
		if(_crsA.IsBC(irow)) continue;
		for(int k = 0; k < T::nWide; ++k)
		{
			const int idx = T::WideIdx(k);
			if(T::IsDiag(idx)) _crsA.Diag(irow) = 0.0;
			else               _crsA(irow, Index(T::StencilSite(_G, st, idx))) = 0.0;
		}
	}
}

template<class T>
int WellMixed_ConcEquationBase<T>::Solve(int &iter, double &resid)
{
	iter = 0; resid = 0;
	if(!_step || !_step->h0 || !_step->h1 || !_step->F0 || !_step->F1)
	{
		_failCause = "the concentration equation was solved without the film of the step";
		return -1;
	}
	YUZEN_PROF("c.solve");
	const WMConcStep &st = *_step;
	_edgeOut = 0;
	SetMatrixA();

	const double dt = _timeStep->DT().now;
	const double th = _theta;
	const double tmid = _timeStep->NowTime() + 0.5 * dt;
	if(_d.size() != _nRows)      { _d.resize(_nRows); _d = 0; }
	if(_dGuess.size() != _nRows) { _dGuess.resize(_nRows); _dGuess = 0; }
	if(_rhsD.size() != _nRows)   _rhsD.resize(_nRows);
	const Vector<double> &c0 = _phi_old;

	//  The boundary rows as A has them; d satisfies them with b - L c0
	for(int i = 0; i < _nRows; i++)
	{
		if(!_crsA.IsBC(i)) continue;
		double v = 0;
		for(int k = 0; k < _crsA.NonZero(i); k++)
		{
			const double a = _crsA.Element(i, k);
			_crsIA.Element(i, k) = a;
			v += a * c0[_crsA.Jcol(i, k)];
		}
		_rhsD[i] = _Rhs[i] - v;
	}

	std::vector<Site> svec(T::nStencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = _G->imin(); ii <= _G->imax(); ii++)
		{
			Site s = _G->template Of<Site>(ii);
			const int i = Index(s);
			if(_crsA.IsBC(i)) continue;
			_G->SetSite(svec, s);

			double M[T::nStencil];
			for(int k = 0; k < T::nStencil; ++k) M[k] = 0;
			const double h1 = (*st.h1)[i], h0 = (*st.h0)[i];
			//  The supply: at the local c, theta d of it is the new level's
			//  (c_S - c0 = theta d); at the nozzle's c_S, all of it is known
			//  (ConcStep.h). Either way the S c0 the film's step holds is
			//  taken out with c0 (h1 - h0)
			const double S = st.src ? (*st.src)[i] : 0.0;
			double diag = h1 + th*dt*((st.r1 ? (*st.r1)[i] : 0.0) - (st.srcSupplied ? 0.0 : S));
			double R = st.evp ? dt * (*st.evp)[i] * c0[i] : 0.0;
			if(st.srcSupplied) R += dt * S * (st.srcConc - c0[i]);
			if(_mms.Enabled())
			{
				double x, y;
				T::SiteXY(_G, s, x, y);
				R += dt * (th*h1 + (1.0 - th)*h0) * _mms.SCnd(x, y, tmid);
			}

			for(int d = 0; d < T::nDir; ++d)
			{
				const int a = d / 2;
				const bool high = (d % 2 == 1);
				const double sg = high ? 1.0 : -1.0;
				const int adj = T::AdjIdx(d);
				const int ia = Index(svec[adj]);
				const int lo = high ? i : ia;
				const double F1 = st.F1[a][lo], F0 = st.F0[a][lo];

				if(_swConvection)
				{
					double f1[4], f0[4], w1[4], w0[4];
					int col[4];
					for(int q = 0; q < 4; ++q)
					{
						col[q] = T::FaceStencil(d, q);
						const int j = Index(svec[col[q]]);
						f0[q] = c0[j];
						f1[q] = LimiterC(j, svec[col[q]]);
					}
					//  The row's cell is L of its high face and R of its low one
					const int rowPos = high ? C_Interpolation_base4::L : C_Interpolation_base4::R;
					_scheme.FaceWeights(f1, F1, rowPos, w1);
					_scheme.FaceWeights(f0, F0, rowPos, w0);
					double e1 = 0, e0 = 0;
					for(int q = 0; q < 4; ++q)
					{
						M[col[q]] += th*dt*sg*F1*w1[q];
						e1 += w1[q]*(f0[q] - c0[i]);
						e0 += w0[q]*(f0[q] - c0[i]);
					}
					R -= dt*sg*(th*F1*e1 + (1.0 - th)*F0*e0);
				}
				else
				{
					//  Without the advection the flow still thins or thickens
					//  the film under the solute
					diag += th*dt*sg*F1;
				}

				if(_swDiffusion)
				{
					const double dn = T::Spacing(_G, a);
					const double K1 = KFace(_field->HH(s), _field->HH(svec[adj]),
					                        _field->CC(s), _field->CC(svec[adj]), dn);
					const double K0 = KFace(h0, (*st.h0)[ia], c0[i], c0[ia], dn);
					diag   += th*dt*K1;
					M[adj] -= th*dt*K1;
					const double dc0 = high ? (c0[ia] - c0[i]) : (c0[i] - c0[ia]);
					R += dt*sg*(th*K1*dc0 + (1.0 - th)*K0*dc0);
				}
			}
			M[T::CenterIdx()] += diag;

			for(int k = 0; k < T::nWide; ++k)
			{
				const int idx = T::WideIdx(k);
				if(T::IsDiag(idx)) _crsIA.Diag(i) = M[idx];
				else               _crsIA(i, Index(svec[idx])) = M[idx];
			}
			_rhsD[i] = R;
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif

	_d = _dGuess;
	int ch;
	{ YUZEN_PROF("c.linsolve"); ch = _bcg.Solve(_crsIA, _rhsD, _d, iter, resid); }
	if(ch == 0)
	{
		//  The boundary rows exactly, not to the tolerance of the solve
		_crsIA.SatisfyBCRows(_rhsD, _d);
		_dGuess = _d;
	}
	else
	{
		Global::Warning<<"Concentration: BiCGStab did not converge: ret = "<< ch << ", iter = "<<iter
		               << ", resid = "<< std::scientific <<resid<< nl;
		bool finite = true;
		for(int i = 0; i < _d.size() && finite; i++) finite = std::isfinite(_d[i]);
		if(!finite) _d = 0;
		_dGuess = 0;
	}
	for(int i = 0; i < _nRows; i++) _phi_new[i] = c0[i] + _d[i];
	if(ch == 0 && st.commit && _swConvection) Reconstruct();
	if(ch == 0 && st.commit) EdgeOutflow();
	//  The ghosts from the real cells, by the boundary rows (b in _Rhs)
	_crsIA.SatisfyBCRows(_Rhs, _phi_new);
	return ch;
}

//------------------------------------------------------------------
//  The theta-weighted solute flux through face d of the real cell st, in
//  +axis: the face value of c1 (the solution) and of c0 by the scheme
//  itself, a limited one evaluated at that very data (FaceValue), and the
//  diffusive flux. Formed from the four points about the face in the same
//  order whichever cell asks.
//------------------------------------------------------------------
template<class T>
double WellMixed_ConcEquationBase<T>::SoluteFlux(Site &s, int d, const std::vector<Site> &svec) const
{
	const WMConcStep &st = *_step;
	const int a = d / 2;
	const bool high = (d % 2 == 1);
	const int i = Index(s);
	const int ia = Index(svec[T::AdjIdx(d)]);
	const int lo = high ? i : ia, hi = high ? ia : i;
	const double F1 = st.F1[a][lo], F0 = st.F0[a][lo];
	const Vector<double> &c0 = _phi_old, &c1 = _phi_new;
	double f1[4], f0[4];
	for(int q = 0; q < 4; ++q)
	{
		const int j = Index(svec[T::FaceStencil(d, q)]);
		f1[q] = c1[j];
		f0[q] = c0[j];
	}
	double G1 = F1 * _scheme.FaceValue(f1, F1);
	double G0 = F0 * _scheme.FaceValue(f0, F0);
	if(_swDiffusion)
	{
		const double dn = T::Spacing(_G, a);
		Site sl = high ? s : svec[T::AdjIdx(d)], sh = high ? svec[T::AdjIdx(d)] : s;
		const double K1 = KFace(_field->HH(sl), _field->HH(sh), _field->CC(sl), _field->CC(sh), dn);
		const double K0 = KFace((*st.h0)[lo], (*st.h0)[hi], c0[lo], c0[hi], dn);
		G1 -= K1 * (c1[hi] - c1[lo]);
		G0 -= K0 * (c0[hi] - c0[lo]);
	}
	return _theta * G1 + (1.0 - _theta) * G0;
}

template<class T>
void WellMixed_ConcEquationBase<T>::Reconstruct()
{
	const WMConcStep &st = *_step;
	const double dt = _timeStep->DT().now;
	const double th = _theta;
	const double tmid = _timeStep->NowTime() + 0.5 * dt;
	for(int a = 0; a < T::nAxis; ++a) _G2[a].assign((size_t)_nRows, 0.0);

	//  Each face once: the high face of every real cell, and the low face
	//  of a real cell whose neighbor there is a ghost
	std::vector<Site> svec(T::nStencil);
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel firstprivate(svec)
	{
#pragma omp for
#endif
		for(int ii = _G->imin(); ii <= _G->imax(); ii++)
		{
			Site s = _G->template Of<Site>(ii);
			const int i = Index(s);
			if(_crsA.IsBC(i)) continue;
			_G->SetSite(svec, s);
			for(int a = 0; a < T::nAxis; ++a)
			{
				_G2[a][(size_t)i] = SoluteFlux(s, 2*a + 1, svec);
				const int il = Index(svec[T::AdjIdx(2*a)]);
				if(_crsA.IsBC(il)) _G2[a][(size_t)il] = SoluteFlux(s, 2*a, svec);
			}
		}
#if defined(_OPENMP) || defined(OMP)
	}	//END of OMP Parallel
#endif

	Vector<double> c1 = _phi_new;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for if(yuzen::opt::Par((_G->imax() - _G->imin() + 1)))
#endif
	for(int ii = _G->imin(); ii <= _G->imax(); ii++)
	{
		Site s = _G->template Of<Site>(ii);
		const int i = Index(s);
		if(_crsA.IsBC(i)) continue;
		const double h1 = (*st.h1)[i], h0 = (*st.h0)[i];
		if(!(h1 > 0.0)) continue;
		double div = 0;
		for(int a = 0; a < T::nAxis; ++a)
			div += _G2[a][(size_t)i] - _G2[a][(size_t)Index(T::Adj(_G, s, T::Dir(2*a)))];
		const double c0 = _phi_old[i];
		double q = h0 * c0 - dt * div;
		if(st.src) q += dt * (*st.src)[i] * (st.srcSupplied ? st.srcConc : th * c1[i] + (1.0 - th) * c0);
		if(st.r1)  q -= dt * (th * (*st.r1)[i] * c1[i] + (1.0 - th) * (*st.r0)[i] * c0);
		if(_mms.Enabled())
		{
			double x, y;
			T::SiteXY(_G, s, x, y);
			q += dt * (th*h1 + (1.0 - th)*h0) * _mms.SCnd(x, y, tmid);
		}
		_phi_new[i] = q / h1;
	}
}

//------------------------------------------------------------------
//  The solute through the edge faces of the open sides, outward, as the
//  committed solve took it: what the cells' balance of h c lost there.
//  With the advection the face fluxes are those Reconstruct formed the
//  solute from (_G2), so the sum of h c over the cells changes by this and
//  the sources alone, to rounding. In order over the cells (not in
//  parallel), so the sum does not depend on the number of threads.
//------------------------------------------------------------------
template<class T>
void WellMixed_ConcEquationBase<T>::EdgeOutflow()
{
	_edgeOut = 0;
	bool open = false;
	for(int d = 0; d < T::nDir; ++d) open = open || BoundaryCondition::IsOpen(_boundary->TypeOfSide(d));
	if(!open) return;
	const WMConcStep &st = *_step;
	const double th = _theta;
	const Vector<double> &c0 = _phi_old, &c1 = _phi_new;
	std::vector<Site> svec(T::nStencil);
	double q = 0;
	for(int ii = _G->imin(); ii <= _G->imax(); ii++)
	{
		Site s = _G->template Of<Site>(ii);
		const int i = Index(s);
		if(_crsA.IsBC(i)) continue;
		_G->SetSite(svec, s);
		for(int d = 0; d < T::nDir; ++d)
		{
			const int ia = Index(svec[T::AdjIdx(d)]);
			if(!_crsA.IsBC(ia) || !BoundaryCondition::IsOpen(_boundary->TypeOfSide(d))) continue;
			const int a = d / 2;
			const bool high = (d % 2 == 1);
			const int lo = high ? i : ia, hi = high ? ia : i;
			double g;
			if(_swConvection) g = _G2[a][(size_t)lo];
			else
			{
				//  The rows without the advection: the cell's own c times
				//  the volume flux, and the diffusive flux
				const double F1 = st.F1[a][lo], F0 = st.F0[a][lo];
				double G1 = F1 * c1[i], G0 = F0 * c0[i];
				if(_swDiffusion)
				{
					const double dn = T::Spacing(_G, a);
					Site sl = high ? s : svec[T::AdjIdx(d)], sh = high ? svec[T::AdjIdx(d)] : s;
					const double K1 = KFace(_field->HH(sl), _field->HH(sh), _field->CC(sl), _field->CC(sh), dn);
					const double K0 = KFace((*st.h0)[lo], (*st.h0)[hi], c0[lo], c0[hi], dn);
					G1 -= K1 * (c1[hi] - c1[lo]);
					G0 -= K0 * (c0[hi] - c0[lo]);
				}
				g = th * G1 + (1.0 - th) * G0;
			}
			q += high ? g : -g;
		}
	}
	_edgeOut = q;
}

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_WM_CONC_EQUATION_BASE_INCLUDED_
