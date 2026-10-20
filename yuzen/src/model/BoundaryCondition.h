/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  BoundaryCondition.h : application of the boundary conditions (symmetry / periodic / pinning / inlet / outlet)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _BOUNDARY_CONDITION_INCLUDED_
#define _BOUNDARY_CONDITION_INCLUDED_


#include "YuzenUtil.h"
using namespace yuzen::util;
#include "CRS_Matrix.h"
#include "LubFields.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif



//---------------------------------
//   Boundary condition class (abstract)
//---------------------------------
//--------------------------------------------------------------------------
//  Boundary conditions.
//
//  The film-height equation is fourth order, so each boundary needs two
//  conditions. Two ghost layers are therefore placed on each side, and the
//  conditions themselves are written into the "ghost rows" of the matrix
//  (the rows marked by CRS_Matrix::SetAsBC). The physical terms skip those
//  rows, so they are not time-integrated.
//
//  What each kind imposes:
//    Symmetry   symmetry; the first and third derivatives are set to 0 (no flux)
//    Periodic   periodic; the value is made equal to the one on the opposite side
//    Pinning    the film height is held fixed (the contact line is pinned)
//    FreshInlet fluid flows in at a fixed film height and concentration
//    Outlet     free outflow; the higher derivatives are set to 0
//--------------------------------------------------------------------------
class BoundaryCondition
{
public:
	// Kinds of field
	enum tField
	{
		Height			= 0,
		HeightJacobi	= 1,
		Conc			= 2
	};

	// Kinds of boundary condition
	enum tBoundaryCondition
	{
		Symmetry,		//symmetry boundary
		Periodic,		//periodic boundary
		Pinning,		//pinning
		FreshInlet,		//inlet boundary
		Outlet			//outlet boundary
	};

	BoundaryCondition(){}

	virtual ~BoundaryCondition(){}

	//--------------------------------------------------------------
	//  Write one condition equation into row i_diag of the matrix.
	//
	//  The arguments are the same throughout:
	//    i_diag  the row the condition is written into (= the ghost point whose
	//            value is being determined)
	//    col_st  the column the stencil starts from
	//    step    the stride by which the stencil advances (1 in 1D, nx for the
	//            y direction in 2D)
	//    b       the value placed on the right-hand side
	//  The direction of the difference follows from the sign of step, so all
	//  four sides can be written with a single piece of code.
	//--------------------------------------------------------------
	/**
	 * @brief Prescribe the value itself: phi = b.
	 * @param crs     coefficient matrix the condition row is written into
	 * @param rhs     right-hand-side vector the condition row is written into
	 * @param i_diag  matrix row the condition is written into (the ghost point being determined)
	 * @param col_st  first column of the stencil
	 * @param step    column increment of the stencil (1 in 1D, nx for the y direction in 2D)
	 * @param b       right-hand side
	 */
	void Dirichlet(CRS_Matrix &crs, Vector<double>& rhs, int i_diag, int col_st, int step, double b);
	/**
	 * @brief Prescribe the first derivative: dphi/dx = b.
	 * @param crs     coefficient matrix the condition row is written into
	 * @param rhs     right-hand-side vector the condition row is written into
	 * @param i_diag  matrix row the condition is written into
	 * @param col_st  first column of the stencil
	 * @param step    column increment of the stencil; its sign sets the direction the difference is taken in
	 * @param b       right-hand side
	 * @param dx      grid spacing along the stencil
	 */
	void Neumann(CRS_Matrix &crs, Vector<double>& rhs,   int i_diag, int col_st, int step, double b, double dx);
	/**
	 * @brief Prescribe the second derivative: d2phi/dx2 = b.
	 * @param crs     coefficient matrix the condition row is written into
	 * @param rhs     right-hand-side vector the condition row is written into
	 * @param i_diag  matrix row the condition is written into
	 * @param col_st  first column of the stencil
	 * @param step    column increment of the stencil
	 * @param b       right-hand side
	 * @param dx      grid spacing along the stencil
	 */
	void Derive2(CRS_Matrix &crs, Vector<double>& rhs,   int i_diag, int col_st, int step, double b, double dx);
	/**
	 * @brief Prescribe the third derivative: d3phi/dx3 = b, which is the no-flux condition.
	 * @param crs     coefficient matrix the condition row is written into
	 * @param rhs     right-hand-side vector the condition row is written into
	 * @param i_diag  matrix row the condition is written into
	 * @param col_st  first column of the stencil
	 * @param step    column increment of the stencil
	 * @param b       right-hand side
	 * @param dx      grid spacing along the stencil
	 */
	void Derive3(CRS_Matrix &crs, Vector<double>& rhs,   int i_diag, int col_st, int step, double b, double dx);
	/**
	 * @brief Tie a point to another one: phi_i_diag = phi_jcol. Used for the periodic condition.
	 * @param crs     coefficient matrix the condition row is written into
	 * @param rhs     right-hand-side vector the condition row is written into
	 * @param i_diag  matrix row the condition is written into
	 * @param jcol    column of the point whose value is copied
	 */
	void EqualVal(CRS_Matrix &crs, Vector<double>& rhs,  int i_diag, int jcol);

	/**
	 * @brief Look up a boundary-condition kind from its name in the input file.
	 * @param typeStr  "Symmetry", "Periodic", "Pinning", "FreshInlet" or "Outlet"
	 * @return The matching kind. An unrecognized string falls through to Symmetry;
	 *         Input.cpp range-checks the string first, so it should not get here.
	 */
	static tBoundaryCondition FromStr(const std::string typeStr)
	{
		if     (typeStr.compare("Symmetry"  ) == 0) return Symmetry;
		else if(typeStr.compare("Periodic"  ) == 0) return Periodic;
		else if(typeStr.compare("Pinning"   ) == 0) return Pinning;
		else if(typeStr.compare("FreshInlet") == 0) return FreshInlet;
		else if(typeStr.compare("Outlet"    ) == 0) return Outlet;
		else return Symmetry;
	}
	/**
	 * @brief Whether liquid crosses an edge of this kind: a FreshInlet or an
	 *        Outlet (an open edge). What crosses it is counted in the
	 *        balance of History.dat (OutflowedVol, OutflowedResin).
	 */
	static bool IsOpen(tBoundaryCondition t) { return t == FreshInlet || t == Outlet; }

	/**
	 * @brief The condition that fills the substrate's ghost cells for an edge
	 *        of the given kind: Periodic across a periodic edge, Symmetry
	 *        (the substrate mirrored at the wall) on every other.
	 */
	static tBoundaryCondition SubstrateBC(tBoundaryCondition edge)
	{
		return edge == Periodic ? Periodic : Symmetry;
	}

	/**
	 * @brief Whether an edge of this kind is a wall: no liquid and no solute
	 *        crosses it (Symmetry, Pinning). Every flux of the film-height
	 *        equation through its face is zero, whatever the ghost cells
	 *        beyond it hold.
	 */
	static bool IsWall(tBoundaryCondition edge)
	{
		return edge == Symmetry || edge == Pinning;
	}

protected:
	//--------------------------------------------------------------
	//  The indices for one boundary (Xmin/Xmax in 1D; Xmin/Xmax/Ymin/Ymax in
	//  2D). Everything that differs from boundary to boundary is packed in here.
	//
	//    layer 0 = outermost ghost, 1 = inner ghost,
	//    2 and 3 = the 1st and 2nd points of the computational domain
	//
	//  The only difference between the four sides is how the row numbers of the
	//  ghost layers are indexed - IL(0,jy), IL(_Nx-1,jy) or IL(ix,0) - so that
	//  is carried by this table and the boundary handling itself is written
	//  only once.
	//--------------------------------------------------------------
	struct EdgeSlot
	{
		int    L[4];		// row numbers of layers 0..3 on this side
		int    opp[4];		// row numbers of layers 0..3 on the opposite boundary (the partner of a periodic condition)
		int    colStart;	// first column of the 4-point stencil (L[0] on the min side, L[3] on the max side)
		int    step;		// column increment of the 4-point stencil (1 for X, _Nx for Y)
		double hb, cb;		// the height and concentration prescribed by FreshInlet

		EdgeSlot() : colStart(0), step(1), hb(0), cb(0)
		{ L[0]=L[1]=L[2]=L[3] = 0; opp[0]=opp[1]=opp[2]=opp[3] = 0; }
	};

	//--------------------------------------------------------------
	//  Only the symmetry boundary and Pinning-Conc are imposed differently in
	//  1D and 2D. In 1D they are written in derivative form (Derive3 / Neumann
	//  / Derive2), in 2D in mirror form (EqualVal). Mathematically they say the
	//  same thing, but the coefficients differ, so the two are kept apart.
	//--------------------------------------------------------------
	enum tEdgeForm { Form1D, Form2D };

	//  Set up the matrix for one boundary. Everything that differs between the
	//  four sides is packed into EdgeSlot, so this single function is all the
	//  body that is needed. For the constraint on the order of insertion into
	//  the matrix, see the comments where this is defined
	//  (BoundaryCondition.cpp).
	void ApplyEdge(CRS_Matrix &crs, Vector<double>& rhs, tField tf,
	               tBoundaryCondition tbc, const EdgeSlot &e, tEdgeForm form);

private:
	BoundaryCondition(const BoundaryCondition& b);
	void operator=(const BoundaryCondition& b);

};


//---------------------------------
//   Boundary condition class (1D)
//---------------------------------
class BoundaryCondition1D : public BoundaryCondition
{
public:
	/** @brief Default-construct with symmetry at both ends and no grid size; call Build() before use. */
	BoundaryCondition1D() : _Nx(0), _tXmin(Symmetry), _tXmax(Symmetry) {}
	/** @brief Same as the default constructor followed by Build(), with nx counting the ghost layers. */
	BoundaryCondition1D(int nx, tBoundaryCondition xmin, tBoundaryCondition xmax) : _Nx(nx), _tXmin(xmin), _tXmax(xmax) {}

	/**
	 * @brief Record the grid size and the kind of condition to impose at each end.
	 * @param nx    total number of points including the ghost layers, i.e. the grid's nx()
	 * @param xmin  kind of condition at the low-x end
	 * @param xmax  kind of condition at the high-x end
	 */
	void Build(int nx, tBoundaryCondition xmin, tBoundaryCondition xmax)
	{
		_Nx = nx;
		_tXmin = xmin;	_tXmax = xmax;
	}

	/**
	 * @brief Write the boundary conditions into the ghost rows of the matrix.
	 * @param crs        matrix being assembled; its four ghost rows are marked as boundary rows and overwritten
	 * @param rhs        right-hand side, written for the same rows
	 * @param t          which field the rows belong to (Height, HeightJacobi or Conc)
	 * @param ini_Field  field the boundary values are read from, e.g. the height and
	 *                   concentration a FreshInlet is held at
	 *
	 * Two rows per end, because the film-height equation is fourth order. The
	 * order in which entries are inserted into the matrix is fixed; see the note
	 * in BoundaryCondition.cpp before changing it.
	 */
	void Apply(CRS_Matrix &crs, Vector<double>& rhs, tField t, Field_WM1D& ini_Field);

	// The kind of condition imposed at each end, as recorded by Build().
	const tBoundaryCondition& TypeXmin() const {return _tXmin;}
	const tBoundaryCondition& TypeXmax() const {return _tXmax;}
	/** @brief The kind of the edge on side d (0: Xmin, 1: Xmax). */
	tBoundaryCondition TypeOfSide(int d) const {return d == 0 ? _tXmin : _tXmax;}

	/**
	 * @brief Impose the boundary conditions directly on a field, in place.
	 * @param field  field whose height, concentration and substrate shape are overwritten
	 *
	 * Used to make an initial or restarted field consistent with the boundary
	 * conditions before time stepping starts. It solves a small linear system in
	 * which the interior rows are the identity and the ghost rows carry the
	 * conditions. The substrate shape is always treated as symmetric at every
	 * boundary, whatever the height and concentration conditions are.
	 */
	void ApplyToField(Field_WM1D& field);

private:
	int _Nx;
	tBoundaryCondition _tXmin;
	tBoundaryCondition _tXmax;

	BoundaryCondition1D(const BoundaryCondition1D&);
	void operator=(const BoundaryCondition1D&);
};

struct stCornerIndex
{
	/*---Generalized corner indices-------
	OIoi  OIii  ||  IIoi    IIii
	OIoo  OIio  ||  IIoo    IIio
	            ================
	OOoi  OOii      IOoi    IOii
	OOoo  OOio      IOoo    IOio
	------------------------------------------*/
	int OOii, OOio, OOoi, OOoo;
	int IIii, IIio, IIoi, IIoo;
	int IOii, IOio, IOoi, IOoo;
	int OIii, OIio, OIoi, OIoo;
};

//  Coefficients of one row of the 5-point corner stencil. For row i,
//     d*x[i] + cII*x[IIii] + cIo*x[IIio] + coI*x[IIoi] + coo*x[IIoo] = r*cornerVal
struct stCornerRow { double d, cII, cIo, coI, coo, r; };

//---------------------------------
//   Boundary condition class (2D)
//---------------------------------
class BoundaryCondition2D : public BoundaryCondition
{
public:
	/** @brief Default-construct with symmetry on all four sides and no grid size; call Build() before use. */
	BoundaryCondition2D() : _Nx(0), _Ny(0), _tXmin(Symmetry), _tXmax(Symmetry), _tYmin(Symmetry), _tYmax(Symmetry) {}
	/** @brief Construct and Build() in one step; the parameters are those of Build(). */
	BoundaryCondition2D(int nx, int ny, tBoundaryCondition xmin, tBoundaryCondition xmax, tBoundaryCondition ymin, tBoundaryCondition ymax)
		: BoundaryCondition2D() { Build(nx, ny, xmin, xmax, ymin, ymax); }

	/**
	 * @brief Record the grid size and the kind of condition to impose on each of the four sides, then build the corner index tables.
	 * @param nx    total number of points in x including the ghost layers, i.e. the grid's nx()
	 * @param ny    total number of points in y including the ghost layers
	 * @param xmin  kind of condition at the low-x side
	 * @param xmax  kind of condition at the high-x side
	 * @param ymin  kind of condition at the low-y side
	 * @param ymax  kind of condition at the high-y side
	 *
	 * Periodic has to be set on both x sides or neither; setting it on one side
	 * only is a fatal error and aborts the run.
	 */
	void Build(int nx, int ny, tBoundaryCondition xmin, tBoundaryCondition xmax, tBoundaryCondition ymin, tBoundaryCondition ymax)
	{
		_Nx = nx; _Ny = ny;
		if( ( xmin == Periodic || xmax == Periodic) && xmax != xmin)
			Global::FatalError<<"Boundary condition mismatch: Periodic must be set on both the Xmin and the Xmax edge."<<nl, Global::FatalError.Abort();

		_tXmin = xmin;	_tXmax = xmax;	_tYmin = ymin;	_tYmax = ymax;

		SetCornerIndex();
	}

	/**
	 * @brief Write the boundary conditions into the ghost rows of the matrix, edges first and then corners.
	 * @param crs        matrix being assembled; the ghost rows are marked as boundary rows and overwritten
	 * @param rhs        right-hand side, written for the same rows
	 * @param t          which field the rows belong to (Height, HeightJacobi or Conc)
	 * @param ini_Field  field the boundary values are read from
	 *
	 * Two rows per side, because the film-height equation is fourth order. The
	 * four corners are handled separately, since two conditions meet there. The
	 * order in which entries are inserted into the matrix is fixed; see the notes
	 * in BoundaryCondition.cpp before changing it.
	 */
	void Apply(CRS_Matrix &crs, Vector<double>& rhs, tField t, Field_WM2D& ini_Field);

	// The kind of condition imposed on each side, as recorded by Build().
	const tBoundaryCondition& TypeXmin() const {return _tXmin;}
	const tBoundaryCondition& TypeXmax() const {return _tXmax;}
	const tBoundaryCondition& TypeYmin() const {return _tYmin;}
	const tBoundaryCondition& TypeYmax() const {return _tYmax;}
	/** @brief The kind of the edge on side d (0: Xmin, 1: Xmax, 2: Ymin, 3: Ymax). */
	tBoundaryCondition TypeOfSide(int d) const
	{ return d == 0 ? _tXmin : d == 1 ? _tXmax : d == 2 ? _tYmin : _tYmax; }

	/**
	 * @brief Impose the boundary conditions directly on a field, in place.
	 * @param field  field whose height, concentration and substrate shape are overwritten
	 *
	 * Used to make an initial or restarted field consistent with the boundary
	 * conditions before time stepping starts. It solves a small linear system in
	 * which the interior rows are the identity and the ghost rows carry the
	 * conditions. The substrate shape is always treated as symmetric on every
	 * side, whatever the height and concentration conditions are.
	 */
	void ApplyToField(Field_WM2D& field);

private:
	inline int IL(int ix, int jy) {return (_Nx*jy) + ix;}
	int _Nx, _Ny;


	tBoundaryCondition _tXmin;
	tBoundaryCondition _tXmax;
	tBoundaryCondition _tYmin;
	tBoundaryCondition _tYmax;

	void ApplyCornersAll(CRS_Matrix &crs, Vector<double>& rhs, tField t, Field_WM2D& ini_Field);
	void ApplyCorners(CRS_Matrix &crs, Vector<double>& rhs, tField t, double cornerVal,
		tBoundaryCondition tbc1, tBoundaryCondition tbc2,
		const stCornerIndex &ci, const stCornerIndex &ciRO, const stCornerIndex &ciOR, const stCornerIndex &ciRR);

	//--------------------------------------------------------------
	//  Shared parts of the corner handling.
	//
	//  There are 5x5 combinations of the conditions in the X and Y directions,
	//  but there are only two things that ever have to be done at a corner.
	//
	//   (a) Make the four ghost rows copies of four points somewhere -> CopyCorner
	//   (b) Place a 5-point stencil on the four ghost rows           -> CornerStencil
	//--------------------------------------------------------------

	// (a) Make OOii/OOio/OOoi/OOoo copies of s_ii/s_io/s_oi/s_oo respectively.
	void CopyCorner(CRS_Matrix &crs, Vector<double>& rhs, const stCornerIndex &ci,
	                int s_ii, int s_io, int s_oi, int s_oo);

	// Three ways of taking the mirror image, determined by whether the X and
	// the Y direction are folded back.
	//   MirrorXY : fold back in both X and Y (Symmetry/Periodic against each other)
	//   MirrorX  : fold back in X only
	//   MirrorY  : fold back in Y only
	void MirrorXY(CRS_Matrix &crs, Vector<double>& rhs, const stCornerIndex &ci, const stCornerIndex &cs)
	{ CopyCorner(crs, rhs, ci, cs.IIoo, cs.IIoi, cs.IIio, cs.IIii); }
	void MirrorX (CRS_Matrix &crs, Vector<double>& rhs, const stCornerIndex &ci, const stCornerIndex &cs)
	{ CopyCorner(crs, rhs, ci, cs.IOoi, cs.IOoo, cs.IOii, cs.IOio); }
	void MirrorY (CRS_Matrix &crs, Vector<double>& rhs, const stCornerIndex &ci, const stCornerIndex &cs)
	{ CopyCorner(crs, rhs, ci, cs.OIio, cs.OIii, cs.OIoo, cs.OIoi); }

	//  The rows come in the order OOoo, OOio, OOoi, OOii.
	//  When setRhs = false, rhs is left untouched (for HeightJacobi).
	void CornerStencil(CRS_Matrix &crs, Vector<double>& rhs, const stCornerIndex &ci,
	                   const stCornerRow row[4], double cornerVal, bool setRhs);

	stCornerIndex ci_mm{}, ci_mp{}, ci_pm{}, ci_pp{};
	void SetCornerIndex();

	BoundaryCondition2D(const BoundaryCondition2D&);
	void operator=(const BoundaryCondition2D&);
};



#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_BOUNDARY_CONDITION_INCLUDED_
