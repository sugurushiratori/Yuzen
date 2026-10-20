/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CZ_ConcEquation_2D.h : concentration equation of the ConcZ model (resolved across the film thickness)
*
*  CZ_ConcEquation.h
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _CZ_CONC_EQUATION_2D_INCLUDED_
#define _CZ_CONC_EQUATION_2D_INCLUDED_

#include "SemiImplicit.h"
#include "ManufacturedSource.h"
#include "FiniteDifferences.h"

#include "LubFields.h"
#include "NonDimensionalNumbers.h"
#include "InitialCondition.h"
#include "BoundaryCondition.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif


/////////////////////////////////////////
/**
 * @brief Concentration equation of the ConcZ model, resolved through the film thickness.
 *
 * Unlike the well-mixed models this one carries c(x, z, t). The film is
 * mapped onto a fixed computational coordinate zeta in [0,1], so the
 * concentration grid stretches with the film as it thins. The equation is
 * written in conservative form for h c on that grid: each cell exchanges
 * solute with its neighbors through face fluxes, and the volume flux
 * through its zeta faces follows from the discrete geometric conservation
 * law with the film height's rate of change (SetDtDh()) and the film's face
 * fluxes (LinkFilmFluxes()), so a uniform c stays uniform to rounding and
 * the solute of a closed film is conserved to rounding (a limited scheme:
 * to the tolerance of SettleLimiter()).
 *
 * No solute crosses the substrate or the free surface: the total (conormal)
 * flux there is zero, the solvent that evaporates leaves the top cell
 * without its solute, which is what concentrates the solute near the
 * surface.
 */
class CZ_ConcEquation_2D : public C_SemiImplicitEquation
{
public:
	CZ_ConcEquation_2D() {}

	/**
	 * @brief Assemble the equation from the grid, properties, field and boundary condition.
	 * @param inp     the whole input file; the blocks this equation needs are read from it
	 * @param prop    fluid properties (viscosity and the concentration dependence)
	 * @param ts      time-step object supplying dt and the current time
	 * @param ff   ConcZ field: h and s live on the 2-D grid, c on the 3-D (x, zeta) grid
	 * @param nonDim  non-dimensional numbers of the problem
	 * @param bc   boundary condition of the 1-D film, applied at the x ends
	 */
	void Build(Input &inp, C_Fluid *prop, C_TimeSteps *ts, Field_CZ *ff, NonDimensionalNumbers *nonDim, BoundaryCondition2D *bc);

	/**
	 * @brief Take dh/dt from the height equation.
	 * @param v  dh/dt of every column; with the layer fluxes it sets the
	 *           flux through the zeta faces, so it has to be set before
	 *           each solve
	 */
	void SetDtDh(const Vector<double>&v)
	{
		_dtdh = v;
	}

	/**
	 * @brief Solve the step again, in one stage, with the concentration of
	 *        the field at cmid.
	 *
	 * For the passes of the height equation that take the coupling of the
	 * two again: the step goes from the same start (the one TimeAdvance()
	 * set), with every concentration-dependent coefficient formed at cmid
	 * instead of at a predictor of its own. The velocity and dh/dt must
	 * have been set for it.
	 * @return 0 on convergence of the linear solver
	 */
	int SolveAtMidpoint(const Vector<double> &cmid, int &iter, double &resid);
	/**
	 * @brief Advance by one step, in one stage, with the concentration of
	 *        the field at cmid (in place of TimeAdvance(), whose predictor
	 *        stage forms its own midpoint).
	 */
	int AdvanceAtMidpoint(const Vector<double> &cmid, int &iter, double &resid)
	{
		_phi_old2 = _phi_old;
		_phi_old  = _phi_new;
		return SolveAtMidpoint(cmid, iter, resid);
	}

	/**
	 * @brief Take hold of the face fluxes of the film-height equation.
	 *
	 * The volume fluxes of the layers of a column are made to add up to
	 * the flux the film moved with, face by face, so that every cell of
	 * the moving zeta grid keeps the volume the film height gives it (the
	 * geometric conservation law; see CZ_ConcEquation_2D.cpp).
	 * @param qx, qy  the flux of the film through x face ix (between ix-1
	 *                and ix) and y face jy, indexed IL(ix, jy)
	 * @param side    the sideways outflow of 1-D of each column
	 * @param pp      the velocity profile PP at the x faces it goes with
	 * @param src     what each column gains through its free surface
	 *                (manufactured source less evaporation and sideways
	 *                outflow), indexed IL(ix, jy)
	 */
	void LinkFilmFluxes(const std::vector<double> *qx, const std::vector<double> *qy,
	                    const std::vector<double> *side, const std::vector<double> *pp,
	                    const std::vector<double> *src)
	{
		_qFilmX = qx; _qFilmY = qy; _sideSink = side; _ppSide = pp; _colSrc = src;
	}

	/** @brief Set the ghosts of the solution from the end state the field holds (see the .cpp). */
	void RefreshGhosts();

	/**
	 * @brief Keep the state a step starts from (KeepStepStart), and go back
	 *        to it (BackToStepStart) when the step has failed and is to be
	 *        solved again at a shorter one: TimeAdvance() and
	 *        AdvanceAtMidpoint() move the unknowns on by a step before they
	 *        solve, and a retry would otherwise start from the failed answer.
	 */
	void KeepStepStart() { _keepOld2 = _phi_old2; _keepOld = _phi_old; _keepNew = _phi_new; }
	void BackToStepStart() { _phi_old2 = _keepOld2; _phi_old = _keepOld; _phi_new = _keepNew; }

	/** @brief The face interpolation of the advection. */
	const C_Interpolation_base4& Scheme() const {return _scheme;}

	/**
	 * @brief The largest Courant rate of the volume fluxes of the last
	 *        assembly (CourantCheck.h): over the cells, the larger flux of
	 *        each pair of opposite faces over the cell's volume, |Q| / (h dx)
	 *        in x and y and |W| / (h dzeta) through the zeta faces, summed.
	 *        Times dt it is the Courant number. 0 before the first assembly.
	 */
	double CourantRate() const;

	// Grid index to matrix row: (i, j) on the 2-D grid, (i, j, k) on the 3-D one
	inline int IL(const int i, const int j) const {return j*_Nx + i;}
	inline int IL(const int i, const int j, const int k) const {return k*(_Nx*_Ny) +j*_Nx + i;}

private:
	C_Interpolation_base4 _scheme;

	//  Mark the ghost points of the 3-D grid as boundary-condition rows.
	virtual void SetFlagBC();

private:
	// Pointer to grid
	C_Grid2D *_G2;
	C_Grid3D *_G3;
	int _Nx, _Ny, _Nz;

	// Pointers to field data
	Field_CZ *_field;

	// Pointer to fluid property
	C_Fluid *_prop;

	// Nondimensional numbers
	NonDimensionalNumbers *_NonDim;

	// Boundary Condition
	BoundaryCondition2D *_boundary;

	C_FiniteDifferences _fdm;

	//  The source of the method of manufactured solutions (verification only)
	C_ManufacturedSource _mms;

	Vector<double> _dtdh;

	//  The face fluxes of the film-height equation (LinkFilmFluxes)
	const std::vector<double> *_qFilmX = nullptr, *_qFilmY = nullptr;
	const std::vector<double> *_sideSink = nullptr, *_ppSide = nullptr, *_colSrc = nullptr;
	//  The volume fluxes of the cells: through x face ix and y face jy of
	//  layer kz, and through the zeta face below cell kz (kz = _Nz-2: the
	//  free surface), indexed IL(ix, jy, kz)
	std::vector<double> _Qx, _Qy, _W;
	//  The free-surface term of each row (its part of the diagonal that a
	//  uniform concentration does not cancel)
	std::vector<double> _react;
	//  The unknowns at the start of the step (KeepStepStart)
	Vector<double> _keepOld2, _keepOld, _keepNew;
	void BuildVolumeFluxes();
	void AdvWeights(const double f[4], double F, double w[4], double &k, int rowPos) const;
	//  A limited scheme: the step again until the midpoint the limiter is
	//  frozen at settles
	int SettleLimiter(int &iter, double &resid);
	static const int kLimiterPasses = 8;
	void ExplicitPart(const Vector<double> &x, Vector<double> &Ax) const;

	//  Assemble the coefficient matrix, including the metric terms that come out
	//  of the mapping onto the zeta coordinate.
	virtual void SetMatrixA();
	//  A predictor and a midpoint corrector; see the .cpp.
	virtual int Solve(int &iter, double &resid);
	int SolveStage(int &iter, double &resid);
	void ApplyBoundaryRows(Vector<double> &x) const;
	//  The ghost rows: the zeta gradient of the conormal flux at the
	//  substrate and the free surface, and periodic or symmetric at both ends
	//  in x and y. The rows are overwritten directly.
	void SetBoundaryCondition(bool atEnd = false);

private:
	CZ_ConcEquation_2D(const CZ_ConcEquation_2D&);
	void operator=(const CZ_ConcEquation_2D&);

};


#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_CZ_CONC_EQUATION_2D_INCLUDED_
