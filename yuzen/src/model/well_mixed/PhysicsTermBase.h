/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  PhysicsTermBase.h : the common base shared by the physics terms of the
*                      film height equation
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _PHYSICS_TERM_BASE_INCLUDED_
#define _PHYSICS_TERM_BASE_INCLUDED_

#include "grid/YuzenGrid.h"
#include "MyMath.h"
#include "Vector.h"
#include "Grids.h"
#include "GridTraits.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

using namespace yuzen::grid;

class NonDimensionalNumbers;

//--------------------------------------------------------------------------
//  The common base shared by the physics terms (pressure terms and body
//  force terms).
//
//  The film height equation is solved in the form
//
//      ∂h/∂t = -A h - B
//
//  (see numerics/ImplicitSolver.h). A physics term is a part that adds its
//  own contribution to this A and this B; it need know nothing at all about
//  how the time integration is arranged. To add a new term, derive from this
//  class and implement the four functions below.
//
//  What each term implements:
//    Build            read the input, and fix the coefficient _Magnitude
//                     and the finite-difference coefficients
//    SetMatrixA       add into the matrix A (the side that carries the
//                     Crank-Nicolson weight)
//    SetJacobiMatrix  add into the Jacobian used by the Newton iteration
//    CalVelocity      write the velocity this term produces into the field
//                     (for visualization and for the mixed velocity)
//
//  The Jacobian may be approximate: that only slows the convergence of the
//  Newton iteration down, it does not change what the iteration converges
//  to. For a strongly nonlinear term such as the disjoining pressure,
//  however, dropping part of it stops the iteration converging at all. See
//  the comments in DisjoiningPressureKernel.h.
//
//  The difference between the dimensions (1-D / 2-D) is absorbed by the
//  traits T. The number of directions T::nDir is 2 in 1-D (Xm, Xp) and 4 in
//  2-D (Xm, Xp, Ym, Yp); the number of axes T::nAxis is 1 and 2. The
//  finite-difference coefficients are held as arrays indexed by that
//  direction and axis, so a derived class simply looks them up with
//  C(d) / CV(a) and need not be aware of the dimension.
//--------------------------------------------------------------------------
template<class T>
class PhysicsTermBase
{
public:
	typedef T                      Tr;
	typedef typename T::Grid       Grid;
	typedef typename T::Site       Site;
	typedef typename T::Field      Field;
	typedef typename T::Boundary   Boundary;
	typedef typename T::tAdjacent  tAdjacent;
	typedef PhysicsFace::tEdge     tEdge;

public:
	/**
	 * @brief Size the finite-difference coefficient arrays for this dimension.
	 *
	 * Everything else (grid, non-dimensional numbers, magnitude) is left
	 * unset until Build() is called.
	 */
	PhysicsTermBase()
	 : _n_stencil(T::nStencil), _enabled(false), _G(0), _NonDim(0), _boundary(0),
	   _Magnitude(1.0), _Nx(0), _Ny(0)
	{
		for(int d = 0; d < T::nDir;  ++d) _c[d].resize(_n_stencil);
		for(int a = 0; a < T::nAxis; ++a) _cV[a].resize(_n_stencil);
	}
	virtual ~PhysicsTermBase(){}

	/**
	 * @brief Read this term's parameters from the input and assemble its
	 *        finite-difference coefficients.
	 * @param inp_phys  the Physics block of the input file
	 * @param g         the grid
	 * @param nonDim    the set of non-dimensional numbers (Bo, Sp, Dp, As*Ma ...)
	 * @param bc        the boundary condition
	 *
	 * Whether the term is active (_enabled) and the overall coefficient
	 * _Magnitude (sign included) are also decided here.
	 */
	virtual void Build(Input::I_Physics &inp_phys, Grid *g, NonDimensionalNumbers *nonDim, Boundary *bc) = 0;

	/**
	 * @brief Add this term's contribution to the matrix A.
	 * @param crsA   the matrix to add into
	 * @param field  the current field
	 *
	 * Only what is put into A receives the Crank-Nicolson weighting
	 * (putting a term into B makes that term backward Euler).
	 * A disabled term (!Enabled()) must return without doing anything.
	 */
	virtual void SetMatrixA(CRS_Matrix &crsA, const Field &field) = 0;

	/**
	 * @brief Add this term to the Jacobian d(A h)/dh used by the Newton iteration.
	 * @param crsJacobi  the matrix to add into
	 * @param field      the current field
	 */
	virtual void SetJacobiMatrix(CRS_Matrix &crsJacobi, const Field &field) = 0;

	/**
	 * @brief Write the velocity produced by this term into the field.
	 * @param field  the field to write into
	 *
	 * Used for visualization, and as an ingredient of the mixed velocity
	 * (LubFields::Cal_MixedVelocity) that advects the concentration.
	 */
	virtual void CalVelocity(Field &field) = 0;

	/** @brief Whether this term is switched on in the input. */
	bool Enabled() const {return _enabled;}

	/**
	 * @brief Whether this term evaluates A h itself, in conservative form
	 *        (AddConservative), instead of as the product of its entries of
	 *        A with h.
	 *
	 * A term that says so still adds its entries to A (the Jacobian starts
	 * from them), but the height equation takes A h of it from
	 * AddConservative. Each face's flux is formed once, from differences of
	 * the film, and is taken out of one cell and put into the other, so
	 * the sum over the cells is the flux through the edges of the domain
	 * to rounding of the fluxes. The product of the entries with h sums to
	 * that only to rounding of the entries times h -- of the order of
	 * eps |A| |h|, which for a stiff term on a thick film is far larger
	 * than the fluxes (tools/analytic.py nozzle_balance).
	 */
	virtual bool ConservativeForm() const { return false; }
	/**
	 * @brief Add this term's A h, in conservative form, to y (rows that are
	 *        not boundary rows): the coefficients those SetMatrixA last
	 *        assembled A with, applied to x, which must be the film height
	 *        of that field.
	 * @param y      the vector to add into
	 * @param x      the film height (the field SetMatrixA was given)
	 * @param crs    the matrix whose boundary rows are skipped
	 * @param field  the field SetMatrixA was given (for the substrate)
	 */
	virtual void AddConservative(Vector<double> &y, const Vector<double> &x,
	                             const CRS_Matrix &crs, const Field &field) const
	{ (void)y; (void)x; (void)crs; (void)field; }

	/**
	 * @brief The volume flux of this term through face d of the cell s in
	 *        the state field (film and concentration alike), as A carries
	 *        it: along the axis of the face, per unit cell measure (divided
	 *        by the spacing across the face), so that the row of s takes
	 *        sum_d (d odd ? +f_d : -f_d), which is this term's A h to
	 *        rounding.
	 *
	 * The film-height equation forms each face's flux once from these
	 * (HeightEquationBase::FaceVolumeFluxes), and the concentration
	 * equation carries the solute with exactly those fluxes: a uniform
	 * concentration then stays uniform whatever the film does. A disabled
	 * term returns 0. svec is the stencil of s (Grid::SetSite), which the
	 * caller forms once for all the faces and terms of the cell.
	 */
	virtual double FaceFlux(Site &s, int d, const Field &field, const std::vector<Site> &svec)
	{ (void)s; (void)d; (void)field; (void)svec; return 0; }
	/**
	 * @brief A volume sink of the cell s that is not a face flux, as a rate
	 *        per unit cell measure (A h of the diagonal part): the sideways
	 *        outflow of the 1-D centrifugal term. 0 for every other term.
	 */
	virtual double SideRate(Site &s, const Field &field) const
	{ (void)s; (void)field; return 0; }

protected:
	/** @brief Finite-difference coefficients of direction d (0:Xm, 1:Xp, 2:Ym, 3:Yp). */
	Vector<double>&       C(int d)       { return _c[d]; }
	const Vector<double>& C(int d) const { return _c[d]; }

	/** @brief Coefficients of axis a (0:X, 1:Y) used to compute the velocity. */
	Vector<double>&       CV(int a)       { return _cV[a]; }
	const Vector<double>& CV(int a) const { return _cV[a]; }

	const int _n_stencil;	// number of stencil points (5 in 1-D, 13 in 2-D)
	bool _enabled;			// whether the input switches this term on

	//  Finite-difference coefficients, fixed by the mesh and the scheme.
	//  They do not depend on the field.
	Vector<double> _c[T::nDir];		// one per direction
	Vector<double> _cV[T::nAxis];	// one per axis (for the velocity)

	Grid *_G;
	NonDimensionalNumbers *_NonDim;
	Boundary *_boundary;

	//  The coefficient that multiplies this term as a whole. The
	//  non-dimensional numbers (Bo, Sp, Dp, As*Ma ...) go in here, and the
	//  sign is decided here as well.
	double _Magnitude;

	int _Nx, _Ny;			// grid sizes (to recognize the edges; _Ny is unused in 1-D)

	/**
	 * @brief Build the finite-difference coefficients, which are fixed by
	 *        the mesh spacing and the scheme alone.
	 *
	 * They do not depend on the field, so this is called once from Build().
	 */
	virtual void SetLinearCoef() = 0;

	/**
	 * @brief The nonlinear coefficient on a face; for most terms this is the
	 *        mobility h³/μ(c).
	 * @param s      the cell the face belongs to
	 * @param tedge  which face (OO means the cell itself)
	 * @param field  the current field
	 */
	virtual double NonlinearTerm(Site &s, tEdge tedge, const Field &field) = 0;

	/**
	 * @brief The above, differentiated with respect to one point (for the
	 *        Jacobian).
	 * @param s      the cell the face belongs to
	 * @param tedge  which face
	 * @param tadj   the neighboring point the derivative is taken with respect to
	 * @param field  the current field
	 */
	virtual double NonlinearJacobi(Site &s, tEdge tedge, tAdjacent tadj, const Field &field) = 0;

private:
	PhysicsTermBase(const PhysicsTermBase&);
	void operator=(const PhysicsTermBase&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_PHYSICS_TERM_BASE_INCLUDED_
