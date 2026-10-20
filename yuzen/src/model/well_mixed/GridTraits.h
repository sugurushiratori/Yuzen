/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  GridTraits.h : the part of "what differs between dimensions" that does
*                 not depend on the physics-term classes
*
*  The types of the grid, the site, the field and the boundary condition,
*  the map from a site to its row number, and the tables of "directions"
*  and "stencils" that the physics-term kernels use.
*
*  The typedefs of the physics-term classes themselves (LaplacePressure and
*  so on) live in WellMixedTraits.h instead. The physics-term headers include
*  this file, so referring to a physics-term class from here would be
*  circular. Only the face identifiers are taken from PhysicsFace.h.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _GRID_TRAITS_INCLUDED_
#define _GRID_TRAITS_INCLUDED_

#include "LubFields.h"
#include "PhysicsFace.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class BoundaryCondition1D;
class BoundaryCondition2D;

//--------------------------------------------------------------------------
//  Traits that collect "what differs between 1-D and 2-D" into types and
//  static functions.
//
//  A physics-term kernel (LaplacePressureKernel and the others) takes these
//  traits as a template argument, so one body of code serves both 1-D and
//  2-D. Three things mainly change with the dimension, and each is provided
//  as a static function here.
//
//    - types      grid, point, field, boundary condition, physics-term base
//    - directions 1-D has 2 faces (Xm,Xp), 2-D has 4 (Xm,Xp,Ym,Yp)
//    - stencil    1-D has 5 points (M m o p P), 2-D has 9
//
//  ** On the order of summation **
//  The order in which terms are summed when they are added into the matrix
//  affects bit-level reproducibility of the results, so Dir(d) (axis order)
//  and StencilDir(k) (stencil order) are used for different purposes. Which
//  one each kernel uses is stated in that kernel's comments. Swapping them
//  breaks the bit-for-bit agreement of the regression tests.
//--------------------------------------------------------------------------
struct WM1D_Grid
{
	// Types the 1-D physics-term kernels work with.
	typedef C_Grid1D                    Grid;
	typedef Xedge                       Site;
	typedef Field_WM1D                  Field;
	typedef BoundaryCondition1D         Boundary;
	typedef VectorField1D<>            VelField;

	/**
	 * @brief Map a site to its row number in the linear system.
	 * @param g  the grid
	 * @param s  the cell
	 * @return   the row index, i.e. the X index of the cell in 1-D
	 */
	static int Index(Grid *g, const Site &s) { return g->ix(s); }

	//--------------------------------------------------------------
	//  The "directions" and "stencil" of the physics-term kernels.
	//
	//  1-D has only the X direction, 2-D has X and Y. Going through these
	//  two, the body of a physics term (nonlinear term, Jacobian, matrix
	//  assembly) can be written once and used in either dimension.
	//--------------------------------------------------------------
	enum { nDir = 2 };						// number of directions (Xm, Xp)
	typedef PhysicsFace::tEdge      tEdge;
	typedef C_Grid1D::tAdjacent     tAdjacent;
	/**
	 * @brief The d-th direction, in axis order.
	 * @param d  direction number, 0 or 1 in 1-D
	 * @return   Xm for d = 0, Xp for d = 1
	 */
	static tEdge Dir(int d)					// the d-th direction
	{ return d == 0 ? PhysicsFace::Xm : PhysicsFace::Xp; }
	/** @brief The identifier that denotes the cell itself rather than a face. */
	static tEdge DirCenter() { return PhysicsFace::OO; }
	/** @brief Whether e is a face direction (OO and any other value give false). */
	static bool IsFace(tEdge e)
	{ return e == PhysicsFace::Xm || e == PhysicsFace::Xp; }
	/**
	 * @brief The direction whose neighboring point is the stencil entry idx.
	 * @return the direction number, or -1 if no direction has idx as its neighbor
	 */
	static int DirOfAdj(int idx)
	{ for(int d = 0; d < nDir; ++d) if(AdjIdx(d) == idx) return d; return -1; }
	/**
	 * @brief Convert a face identifier back to a direction number.
	 * @return the direction number, or -1 if e is not a face
	 */
	static int DirIndex(tEdge e)
	{ for(int d = 0; d < nDir; ++d) if(Dir(d) == e) return d; return -1; }
	/** @brief Whether the stencil point s is skipped when computing velocities
	 *         (the outermost point). */
	static bool SkipInVel(int s) { return s == C_Grid1D::P; }

	/**
	 * @brief Stencil index of the j-th of the four points that straddle the
	 *        face in direction d, ordered LL, L, R, RR.
	 */
	static int FaceStencil(int d, int j)
	{
		static const int m[2][4] = { {C_Grid1D::M, C_Grid1D::m, C_Grid1D::o, C_Grid1D::p},
		                             {C_Grid1D::m, C_Grid1D::o, C_Grid1D::p, C_Grid1D::P} };
		return m[d][j];
	}
	/**
	 * @brief Stencil index of the j-th of the five points along axis a,
	 *        ordered far negative, near negative, center, near positive,
	 *        far positive.
	 */
	static int AxisStencil(int a, int j)
	{
		(void)a;
		static const int m[5] = {C_Grid1D::M, C_Grid1D::m, C_Grid1D::o, C_Grid1D::p, C_Grid1D::P};
		return m[j];
	}

	/**
	 * @brief The neighboring site on the far side of face e.
	 * @return that neighbor, or s itself when e is OO
	 */
	static Site Adj(Grid *g, Site s, tEdge e)
	{
		switch(e)
		{
		case PhysicsFace::Xm: return g->ShiftXm(s);
		case PhysicsFace::Xp: return g->ShiftXp(s);
		default:                   return s;
		}
	}

	// The stencil points written into the matrix. This ordering is the order
	// in which the entries are inserted into the matrix.
	//   Wide   : the 5 points used by the 4-point interpolation (M, m, o, p, P)
	//   Narrow : the 3 points used by the second derivative (m, o, p)
	enum { nStencil = 5 };                  // number of stencil points (M m o p P)
	enum { nWide = 5, nNarrow = 3 };
	/**
	 * @brief The k-th stencil index of the wide set, in matrix insertion order.
	 * @return one of M, m, o, p, P
	 */
	static int  WideIdx(int k)   { static const int a[5] = {C_Grid1D::M, C_Grid1D::m, C_Grid1D::o, C_Grid1D::p, C_Grid1D::P}; return a[k]; }
	/**
	 * @brief The k-th stencil index of the narrow set, in matrix insertion order.
	 * @return one of m, o, p
	 */
	static int  NarrowIdx(int k) { static const int a[3] = {C_Grid1D::m, C_Grid1D::o, C_Grid1D::p}; return a[k]; }
	// Full : every point used by the fourth derivative (Laplace pressure).
	enum { nFull = 5 };
	/**
	 * @brief The k-th stencil index of the full set, in matrix insertion order.
	 *
	 * Used by the fourth-order operator of the Laplace pressure term, which
	 * touches every point of the stencil.
	 */
	static int  FullIdx(int k)   { static const int a[5] = {C_Grid1D::M, C_Grid1D::m, C_Grid1D::o, C_Grid1D::p, C_Grid1D::P}; return a[k]; }
	/** @brief Whether the stencil index idx is the center point, i.e. the diagonal. */
	static bool IsDiag(int idx)  { return idx == C_Grid1D::o; }
	// Axes (X only in 1-D). Direction d = 2*axis + side, with side 0 the
	// negative side.
	enum { nAxis = 1 };
	/**
	 * @brief Mesh spacing along the given axis.
	 * @param g  the grid
	 * @param axis  axis number; ignored in 1-D, where the spacing is always dx
	 */
	static double Spacing(Grid *g, int axis) { (void)axis; return g->dx(); }
	/** @brief Coordinates of a site; y is 0 in 1-D. */
	static void SiteXY(Grid *g, Site &s, double &x, double &y) { x = g->X(s); y = 0.0; }
	/** @brief Stencil index of the point neighboring the cell across direction d. */
	static int AdjIdx(int d) { return d == 0 ? C_Grid1D::m : C_Grid1D::p; }
	/** @brief Stencil index of the center point. */
	static int CenterIdx()   { return C_Grid1D::o; }

	//  There are two orderings of the directions, and which one is used
	//  affects the numbers.
	//   Dir(d)       : axis order over X,Y (Xm, Xp, Ym, Yp). The summations
	//                  in SetMatrixA follow this order.
	//   StencilDir(k): the order of the points along the stencil (Ym, Xm,
	//                  Xp, Yp in 2-D). The diagonal of SetJacobiMatrix and
	//                  the Too summation follow this order.
	//  A floating-point sum depends on the order it is taken in, so the two
	//  must not be swapped.
	/**
	 * @brief The direction number of the k-th direction in stencil order.
	 * @return k itself in 1-D, where the two orderings coincide
	 */
	static int StencilDir(int k) { return k; }		// the two coincide in 1-D

	//--------------------------------------------------------------
	//  Coordinates of a face. Used to measure the distance from the
	//  center of rotation.
	//--------------------------------------------------------------
	/**
	 * @brief Coordinate of the face te, measured along the axis that face
	 *        faces.
	 */
	static double FaceCoord(Grid *g, Site &s, tEdge te)
	{ return te == PhysicsFace::Xm ? g->X(g->HalfXm(s)) : g->X(g->HalfXp(s)); }
	//  The mean coordinate of the two cell centers that straddle the face
	//  te. On a uniform grid this equals FaceCoord, so in 1-D we simply
	//  return FaceCoord.
	/**
	 * @brief Mean coordinate of the two cell centers that straddle the face te.
	 * @return FaceCoord itself in 1-D, the two being equal on a uniform grid
	 */
	static double CellAvgCoord(Grid *g, Site &s, tEdge te)
	{ return FaceCoord(g, s, te); }

	//--------------------------------------------------------------
	//  How a physics term recognizes the edges of the grid.
	//  Every term makes the same test, so it is kept here once rather
	//  than repeated in each term.
	//--------------------------------------------------------------
	/**
	 * @brief Record the grid sizes a term needs in order to recognize the
	 *        edges of the domain.
	 * @param g  the grid
	 * @param Nx  set to the number of points in X
	 * @param Ny  set to 0 in 1-D
	 */
	static void SetExtent(Grid *g, int &Nx, int &Ny) { Nx = g->nx(); Ny = 0; }
	/**
	 * @brief Whether a velocity may be computed at this site.
	 * @return false in the halo region and in the first row inside it
	 */
	static bool InsideForVelocity(Grid *g, int Nx, int Ny, const Site &s)
	{
		(void)Ny;
		const int ix = g->ix(const_cast<Site&>(s));
		return !(ix < 2 || Nx-2 < ix);
	}
	/**
	 * @brief Whether the face te of site s is the face of a wall edge
	 *        (B::IsWall: Symmetry, Pinning), seen from either side: the
	 *        first cell inside (ix 2, Nx-3) or the inner ghost (ix 1, Nx-2).
	 *
	 * A term whose flux is not a gradient (the upwinded body forces and the
	 * Marangoni stress) carries a non-zero flux through that face when its
	 * face value is formed from the ghosts: a mirrored film is not zero at
	 * the wall. Such a term sets its flux there to zero.
	 */
	template<class B>
	static bool WallFace(Grid *g, const B *bc, int Nx, int Ny, Site &s, tEdge te)
	{
		(void)Ny;
		if(!bc) return false;
		const int ix = g->ix(s);
		if((te == PhysicsFace::Xm && ix == 2) || (te == PhysicsFace::Xp && ix == 1))
			return B::IsWall(bc->TypeXmin());
		if((te == PhysicsFace::Xp && ix == Nx-3) || (te == PhysicsFace::Xm && ix == Nx-2))
			return B::IsWall(bc->TypeXmax());
		return false;
	}

	//--------------------------------------------------------------
	//  Writing the velocity on a face.
	//  A physics term computes only the flux through the face; choosing
	//  the position of that face and the component to write is left to
	//  the traits.
	//--------------------------------------------------------------
	/**
	 * @brief Film height on the negative-side face of the given axis,
	 *        taken as the mean of the cells on either side.
	 */
	static double FaceHeight(Grid *g, const Field &f, Site &s, int axis)
	{
		(void)axis;
		return 0.5*(f.HH(g->ShiftXm(s)) + f.HH(s));
	}
	/**
	 * @brief Dimension-specific suppression of a flux.
	 * @return q unchanged; 1-D has nothing to suppress
	 */
	static double AxisMask(int Ny, int axis, double q)
	{ (void)Ny; (void)axis; return q; }
	/**
	 * @brief Write the velocity u onto the negative-side face of the given axis.
	 */
	static void WriteFaceVel(Grid *g, VelField &v, Site &s, int axis, double u)
	{ (void)axis; v.x(g->HalfXm(s)) = u; }

	/**
	 * @brief The site that the stencil index idx refers to, relative to s.
	 * @return the shifted site, or s itself for the center point
	 */
	static Site StencilSite(Grid *g, Site s, int idx)
	{
		switch(idx)
		{
		case C_Grid1D::M: return g->ShiftXmm(s);
		case C_Grid1D::m: return g->ShiftXm(s);
		case C_Grid1D::p: return g->ShiftXp(s);
		case C_Grid1D::P: return g->ShiftXpp(s);
		default:          return s;
		}
	}
};

struct WM2D_Grid
{
	typedef C_Grid2D                    Grid;
	typedef Zface                       Site;
	typedef Field_WM2D                  Field;
	typedef BoundaryCondition2D         Boundary;
	typedef VectorField2D<>            VelField;

	/**
	 * @brief Map a site to its row number.
	 * @return the linear index of the cell, running over both grid directions
	 */
	static int Index(Grid *g, const Site &s) { return g->IL(s); }

	//--- Directions and stencil of the physics-term kernels (see WM1D_Grid) ---
	enum { nDir = 4 };						// Xm, Xp, Ym, Yp
	typedef PhysicsFace::tEdge      tEdge;
	typedef C_Grid2D::tAdjacent     tAdjacent;
	/**
	 * @brief The d-th direction, in axis order.
	 * @return Xm, Xp, Ym, Yp for d = 0, 1, 2, 3
	 */
	static tEdge Dir(int d)
	{
		switch(d)
		{
		case 0:  return PhysicsFace::Xm;
		case 1:  return PhysicsFace::Xp;
		case 2:  return PhysicsFace::Ym;
		default: return PhysicsFace::Yp;
		}
	}
	/** @brief The identifier that denotes the cell itself rather than a face. */
	static tEdge DirCenter() { return PhysicsFace::OO; }
	/** @brief Whether e is a face direction (OO and any other value give false). */
	static bool IsFace(tEdge e)
	{ return e == PhysicsFace::Xm || e == PhysicsFace::Xp
	      || e == PhysicsFace::Ym || e == PhysicsFace::Yp; }
	/**
	 * @brief The direction whose neighboring point is the stencil entry idx.
	 * @return the direction number, or -1 if no direction has idx as its neighbor
	 */
	static int DirOfAdj(int idx)
	{ for(int d = 0; d < nDir; ++d) if(AdjIdx(d) == idx) return d; return -1; }
	/**
	 * @brief Convert a face identifier back to a direction number.
	 * @return the direction number, or -1 if e is not a face
	 */
	static int DirIndex(tEdge e)
	{ for(int d = 0; d < nDir; ++d) if(Dir(d) == e) return d; return -1; }
	/** @brief Whether the stencil point s is skipped when computing velocities
	 *         (the outermost points in X and in Y). */
	static bool SkipInVel(int s) { return s == C_Grid2D::Po || s == C_Grid2D::oP; }

	/**
	 * @brief Stencil index of the j-th of the five points along axis a,
	 *        ordered far negative, near negative, center, near positive,
	 *        far positive.
	 */
	static int AxisStencil(int a, int j)
	{
		static const int m[2][5] = { {C_Grid2D::Mo, C_Grid2D::mo, C_Grid2D::oo, C_Grid2D::po, C_Grid2D::Po},
		                             {C_Grid2D::oM, C_Grid2D::om, C_Grid2D::oo, C_Grid2D::op, C_Grid2D::oP} };
		return m[a][j];
	}
	/**
	 * @brief Stencil index of the j-th of the four points that straddle the
	 *        face in direction d, ordered LL, L, R, RR.
	 */
	static int FaceStencil(int d, int j)
	{
		static const int m[4][4] = { {C_Grid2D::Mo, C_Grid2D::mo, C_Grid2D::oo, C_Grid2D::po},
		                             {C_Grid2D::mo, C_Grid2D::oo, C_Grid2D::po, C_Grid2D::Po},
		                             {C_Grid2D::oM, C_Grid2D::om, C_Grid2D::oo, C_Grid2D::op},
		                             {C_Grid2D::om, C_Grid2D::oo, C_Grid2D::op, C_Grid2D::oP} };
		return m[d][j];
	}

	/**
	 * @brief The neighboring site on the far side of face e.
	 * @return that neighbor, or s itself when e is OO
	 */
	static Site Adj(Grid *g, Site s, tEdge e)
	{
		switch(e)
		{
		case PhysicsFace::Xm: return g->ShiftXm(s);
		case PhysicsFace::Xp: return g->ShiftXp(s);
		case PhysicsFace::Ym: return g->ShiftYm(s);
		case PhysicsFace::Yp: return g->ShiftYp(s);
		default:                   return s;
		}
	}

	enum { nStencil = 13 };                 // number of stencil points
	enum { nWide = 9, nNarrow = 5 };
	/**
	 * @brief The k-th stencil index of the wide set, in matrix insertion order.
	 */
	static int  WideIdx(int k)
	{ static const int a[9] = {C_Grid2D::oM, C_Grid2D::om, C_Grid2D::Mo, C_Grid2D::mo, C_Grid2D::oo,
	                           C_Grid2D::po, C_Grid2D::Po, C_Grid2D::op, C_Grid2D::oP}; return a[k]; }
	/**
	 * @brief The k-th stencil index of the narrow set, in matrix insertion order.
	 */
	static int  NarrowIdx(int k)
	{ static const int a[5] = {C_Grid2D::om, C_Grid2D::mo, C_Grid2D::oo, C_Grid2D::po, C_Grid2D::op}; return a[k]; }
	enum { nFull = 13 };
	/**
	 * @brief The k-th stencil index of the full set, in matrix insertion order.
	 *
	 * Wider than the 1-D counterpart: the 2-D biharmonic operator also has
	 * cross terms, so the four diagonal points mm, pm, mp, pp are included.
	 */
	static int  FullIdx(int k)
	{ static const int a[13] = {C_Grid2D::oM, C_Grid2D::mm, C_Grid2D::om, C_Grid2D::pm,
	                            C_Grid2D::Mo, C_Grid2D::mo, C_Grid2D::oo, C_Grid2D::po, C_Grid2D::Po,
	                            C_Grid2D::mp, C_Grid2D::op, C_Grid2D::pp, C_Grid2D::oP}; return a[k]; }
	/** @brief Whether the stencil index idx is the center point, i.e. the diagonal. */
	static bool IsDiag(int idx) { return idx == C_Grid2D::oo; }
	enum { nAxis = 2 };
	/**
	 * @brief Mesh spacing along the given axis.
	 * @param g  the grid
	 * @param axis  0 for X, 1 for Y
	 */
	static double Spacing(Grid *g, int axis) { return axis == 0 ? g->dx() : g->dy(); }
	/** @brief Coordinates of a site. */
	static void SiteXY(Grid *g, Site &s, double &x, double &y) { x = g->X(s); y = g->Y(s); }
	/** @brief Stencil index of the point neighboring the cell across direction d. */
	static int AdjIdx(int d)
	{
		switch(d)
		{
		case 0:  return C_Grid2D::mo;
		case 1:  return C_Grid2D::po;
		case 2:  return C_Grid2D::om;
		default: return C_Grid2D::op;
		}
	}
	/** @brief Stencil index of the center point. */
	static int CenterIdx() { return C_Grid2D::oo; }

	//  See the explanation in WM1D_Grid. In 2-D the two orderings differ.
	/**
	 * @brief The direction number of the k-th direction in stencil order.
	 * @return the directions in the order Ym, Xm, Xp, Yp, which differs from
	 *         the axis order used by Dir
	 */
	static int StencilDir(int k)
	{ static const int a[4] = {2, 0, 1, 3}; return a[k]; }	// Ym, Xm, Xp, Yp

	//--- Coordinates of a face (see WM1D_Grid) ---
	/**
	 * @brief Coordinate of the face te, measured along the axis that face
	 *        faces.
	 */
	static double FaceCoord(Grid *g, Site &s, tEdge te)
	{
		switch(te)
		{
		case PhysicsFace::Xm: return g->X(g->HalfXm(s));
		case PhysicsFace::Xp: return g->X(g->HalfXp(s));
		case PhysicsFace::Ym: return g->Y(g->HalfYm(s));
		default:              return g->Y(g->HalfYp(s));
		}
	}
	//  In 2-D the mean of the cell centers is used. On a uniform grid that
	//  is the same value as FaceCoord, but the two do not agree down to
	//  the rounding, so they are kept separate.
	/**
	 * @brief Mean coordinate of the two cell centers that straddle the face te.
	 *
	 * Unlike the 1-D version this really averages the two cell centers; on a
	 * uniform grid that equals FaceCoord, but not down to the last bit.
	 */
	static double CellAvgCoord(Grid *g, Site &s, tEdge te)
	{
		const int ix = g->ix(s), jy = g->iy(s);
		switch(te)
		{
		case PhysicsFace::Xm: return 0.5*(g->XOf<Site>(ix-1,jy) + g->XOf<Site>(ix,jy));
		case PhysicsFace::Xp: return 0.5*(g->XOf<Site>(ix+1,jy) + g->XOf<Site>(ix,jy));
		case PhysicsFace::Ym: return 0.5*(g->YOf<Site>(ix,jy-1) + g->YOf<Site>(ix,jy));
		default:              return 0.5*(g->YOf<Site>(ix,jy+1) + g->YOf<Site>(ix,jy));
		}
	}

	//--- Recognizing the edges of the grid (see WM1D_Grid) ---
	/**
	 * @brief Record the grid sizes a term needs in order to recognize the
	 *        edges of the domain.
	 * @param g  the grid
	 * @param Nx  set to the number of points in X
	 * @param Ny  set to the number of points in Y
	 */
	static void SetExtent(Grid *g, int &Nx, int &Ny) { Nx = g->nx(); Ny = g->ny(); }
	/**
	 * @brief Whether a velocity may be computed at this site.
	 *
	 * The Y check is skipped when Ny is 1, so that a grid with a single cell
	 * in Y is not excluded everywhere.
	 */
	static bool InsideForVelocity(Grid *g, int Nx, int Ny, const Site &s)
	{
		Site &z = const_cast<Site&>(s);
		const int ix = g->ix(z), jy = g->iy(z);
		if(ix < 2 || Nx-2 < ix) return false;
		if(Ny != 1 && (jy < 2 || Ny-2 < jy)) return false;
		return true;
	}
	/** @brief Whether face te of s is the face of a wall edge (see WM1D_Grid). */
	template<class B>
	static bool WallFace(Grid *g, const B *bc, int Nx, int Ny, Site &s, tEdge te)
	{
		if(!bc) return false;
		const int ix = g->ix(s), jy = g->iy(s);
		if((te == PhysicsFace::Xm && ix == 2) || (te == PhysicsFace::Xp && ix == 1))
			return B::IsWall(bc->TypeXmin());
		if((te == PhysicsFace::Xp && ix == Nx-3) || (te == PhysicsFace::Xm && ix == Nx-2))
			return B::IsWall(bc->TypeXmax());
		if((te == PhysicsFace::Ym && jy == 2) || (te == PhysicsFace::Yp && jy == 1))
			return B::IsWall(bc->TypeYmin());
		if((te == PhysicsFace::Yp && jy == Ny-3) || (te == PhysicsFace::Ym && jy == Ny-2))
			return B::IsWall(bc->TypeYmax());
		return false;
	}

	//--- Writing the velocity on a face (see WM1D_Grid) ---
	/**
	 * @brief Film height on the negative-side face of the given axis,
	 *        taken as the mean of the cells on either side.
	 */
	static double FaceHeight(Grid *g, const Field &f, Site &s, int axis)
	{
		return axis == 0 ? 0.5*(f.HH(g->ShiftXm(s)) + f.HH(s))
		                 : 0.5*(f.HH(g->ShiftYm(s)) + f.HH(s));
	}
	//  Ny == 6 marks a grid that is treated as holding only one cell in Y;
	//  the velocity in Y is dropped on such a grid.
	/**
	 * @brief Dimension-specific suppression of a flux.
	 * @return 0 for the Y flux on a grid that is only one cell deep in Y,
	 *         q unchanged otherwise
	 */
	static double AxisMask(int Ny, int axis, double q)
	{ return (axis == 1 && Ny == 6) ? 0.0 : q; }
	/**
	 * @brief Write the velocity u onto the negative-side face of the given axis.
	 */
	static void WriteFaceVel(Grid *g, VelField &v, Site &s, int axis, double u)
	{
		if(axis == 0) v.x(g->HalfXm(s)) = u;
		else          v.y(g->HalfYm(s)) = u;
	}

	/**
	 * @brief The site that the stencil index idx refers to, relative to s.
	 * @return the shifted site, or s itself for the center point
	 */
	static Site StencilSite(Grid *g, Site s, int idx)
	{
		switch(idx)
		{
		case C_Grid2D::oM: return g->ShiftYmm(s);
		case C_Grid2D::mm: return g->ShiftXmYm(s);
		case C_Grid2D::pm: return g->ShiftXpYm(s);
		case C_Grid2D::mp: return g->ShiftXmYp(s);
		case C_Grid2D::pp: return g->ShiftXpYp(s);
		case C_Grid2D::om: return g->ShiftYm(s);
		case C_Grid2D::Mo: return g->ShiftXmm(s);
		case C_Grid2D::mo: return g->ShiftXm(s);
		case C_Grid2D::po: return g->ShiftXp(s);
		case C_Grid2D::Po: return g->ShiftXpp(s);
		case C_Grid2D::op: return g->ShiftYp(s);
		case C_Grid2D::oP: return g->ShiftYpp(s);
		default:           return s;
		}
	}
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_GRID_TRAITS_INCLUDED_
