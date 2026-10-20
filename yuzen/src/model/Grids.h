/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Grids.h : the computational grids (RectGrid with ghost layers and neighbor specifiers added)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _GRID_INCLUDED_
#define _GRID_INCLUDED_


#include "grid/YuzenGrid.h"

#include "Input.h"
#include "CRS_Matrix.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

using namespace yuzen::grid;

//////////////////////////////////////////////////////
// yuzen::grid::RectGrid with just what Yuzen needs added to it.
//
//   - The number of ghost layers, and the extent of the real domain with the
//     ghosts excluded
//   - A test for where in the ghost layers a given site lies (CheckGhost)
//   - Shifts to the two-point and diagonal neighbors (the base class only
//     goes as far as the one-point neighbor)
//   - Conversion to a row number in the matrix (IL)
//
// The base class is yuzen::grid::RectGrid<N>. The inheritance is public
// because the Build of a field (Field) takes a const RectGrid<N>&.
//////////////////////////////////////////////////////

class C_Grid1D : public RectGrid<1>
{
public:
	//  Positions of the five stencil points, reaching two points out on each
	//  side of the site itself (o). A physical term indexes the array laid out
	//  by SetSite() with these enumerators.
	enum tAdjacent	// neighbor cell specifiers
	{
		M = 0,		// two points back
		m = 1,		// one point back
		o = 2,		// the site itself
		p = 3,		// one point forward
		P = 4		// two points forward
	};
	//  Whether the site lies inside the computational domain, on a boundary
	//  itself, or outside it. A physical term skips everything but Inside, and
	//  the boundary conditions fill bXmin / bXmax.
	enum tGhost
	{
		Inside = 0,	// inside the computational domain
		Outside= 1,	// outermost ghost layer (set by the condition formula)
		bXmin   = 2,	// boundary on the Xmin side
		bXmax   = 3	// boundary on the Xmax side
	};

	/** @brief Construct an empty grid with the default two ghost layers per side; call Build() to give it an extent. */
	C_Grid1D() : _num_ghostX(2) {}
	/**
	 * @brief Wrap an already-sized base grid.
	 * @param g  base grid whose point count is taken to include the ghost layers already
	 */
	C_Grid1D(const RectGrid<1> &g);

	virtual ~C_Grid1D(){}

	/** @brief Copy the extent, the point count and the ghost-layer count from another grid. */
	void operator=(const C_Grid1D &g);

	/**
	 * @brief Build the grid from the [Grid] input block.
	 * @param g  input grid block; its Nx is the number of interior points only
	 *
	 * Two ghost layers are added on each side, so nx() comes out as Nx + 4.
	 */
	void Build(const Input::I_Grid &g);

	// Coordinates of the ends of the computational domain, ghost layers excluded.
	// Held separately because the base RectGrid xmin()/xmax() include the ghosts.
	double XminReal() const {return _Xmin_real;}
	double XmaxReal() const {return _Xmax_real;}

	// Index range of the computational domain (half-open: IminReal <= ix < ImaxReal).
	int IminReal() const {return _num_ghostX;}
	int ImaxReal() const {return nx() - _num_ghostX;}

	/**
	 * @brief Number of ghost layers on one side.
	 * @return 2. The film-height equation is fourth order, so two layers are needed.
	 */
	int NumGhostX() const {return _num_ghostX;}

	/**
	 * @brief Classify a site as interior, on a boundary, or outside.
	 * @return Inside, bXmin / bXmax for the boundary points themselves, or Outside for the outermost ghost layer
	 */
	tGhost CheckGhost(AnySite s) const
	{
		if (_num_ghostX <= ix(s) && ix(s) < nx() -_num_ghostX) return Inside;
		else if (ix(s) == _num_ghostX-1    ) return bXmin;
		else if (ix(s) == nx() -_num_ghostX) return bXmax;
		else return Outside;
	}

	/**
	 * @brief Fill svec with the five stencil sites centered on S0.
	 * @param svec  output array, indexed afterwards with the tAdjacent enumerators; must hold at least 5 entries
	 * @param S0    center site
	 *
	 * The inner loop of a physical term calls this once per site and then
	 * reads neighbors through tAdjacent.
	 */
	template<class S> inline void SetSite(std::vector<Site<S> > &svec, Site<S> S0) const
	{
		assert(svec.size() >= (size_t)_n_stencil );

		svec[M] = ShiftXmm(S0);
		svec[m] = ShiftXm (S0);
		svec[o] = S0;
		svec[p] = ShiftXp (S0);
		svec[P] = ShiftXpp(S0);
	}

	// Two-point neighbor shifts; the base grid only provides one-point shifts.
	template<class S> inline Site<S> ShiftXmm(Site<S> s) const { return ShiftXm(ShiftXm(s)); }
	template<class S> inline Site<S> ShiftXpp(Site<S> s) const { return ShiftXp(ShiftXp(s)); }

private:
	int _num_ghostX;
	double _Xmin_real;
	double _Xmax_real;
	static const int _n_stencil;
};

///////////////////////////////////////////////////
//--------------------------------------------------------------------------
//  The 2-D computational grid. Its role is the same as that of C_Grid1D.
//
//  The stencil has 13 points (a cross plus four diagonals). The names are two
//  characters long: the first gives the position in x, the second in y.
//    M/m = two points back / one point back, o = the site itself,
//    p/P = one point forward / two points forward
//  For example: oo = the site itself, po = one point forward in x,
//  mp = one point back in x and one point forward in y.
//
//  There are two ghost layers per side in each of x and y. The corner ghosts
//  (Corner**) are where the boundary conditions of two directions meet, and
//  they need to be handled separately.
//--------------------------------------------------------------------------
class C_Grid2D : public RectGrid<2>
{
public:
	//  Positions of the 13 stencil points (first character is x, second is y).
	enum tAdjacent	// neighbor cell specifiers
	{
		oM = 0,
		mm = 1,
		om = 2,
		pm = 3,
		Mo = 4,
		mo = 5,
		oo = 6,
		po = 7,
		Po = 8,
		mp = 9,
		op = 10,
		pp = 11,
		oP = 12
	};
	//  Whether the site is inside the computational domain, on which boundary,
	//  in a corner, or outside.
	enum tGhost
	{
		Inside = 0,
		Outside= 1,
		bXmin   = 2,
		bXmax   = 3,
		bYmin   = 4,
		bYmax   = 5,
		//  Corner ghosts. Both the x and the y boundary condition apply here.
		CornerMM = 6,
		CornerPM = 7,
		CornerMP = 8,
		CornerPP = 9
	};

	/** @brief Construct an empty grid with the default two ghost layers per side; call Build() to give it an extent. */
	C_Grid2D() : _num_ghostX(2), _num_ghostY(2) {}
	/**
	 * @brief Wrap an already-sized base grid.
	 * @param g  base grid whose point counts are taken to include the ghost layers already
	 */
	C_Grid2D(const RectGrid<2> &g);

public:
	virtual ~C_Grid2D(){}
	/**
	 * @brief Build the grid from the [Grid] input block.
	 * @param g  input grid block; Nx and Ny count interior points only
	 *
	 * Two ghost layers are added on each side in both directions, so nx()
	 * comes out as Nx + 4 and ny() as Ny + 4.
	 */
	void Build(const Input::I_Grid &g);

	/** @brief Copy the extent, the point counts and the ghost-layer counts from another grid. */
	void operator=(const C_Grid2D &g);

	/**
	 * @brief Classify a site as interior, on one of the four boundaries, in a corner, or outside.
	 * @return Inside, bXmin / bXmax / bYmin / bYmax, one of the four Corner** values
	 *         where two boundaries meet, or Outside for the outermost ghost layer
	 */
	tGhost CheckGhost(AnySite s) const
	{
		if ( IsRangeX(s) && IsRangeY(s) ) return Inside;
		else if(IsBoundX0(s) && IsBoundY0(s) ) return CornerMM;
		else if(IsBoundX1(s) && IsBoundY0(s) ) return CornerPM;
		else if(IsBoundX0(s) && IsBoundY1(s) ) return CornerMP;
		else if(IsBoundX1(s) && IsBoundY1(s) ) return CornerPP;
		else if(IsBoundX0(s) && IsRangeY(s) ) return bXmin;
		else if(IsRangeX(s)  && IsBoundY0(s)) return bYmin;
		else if(IsBoundX1(s) && IsRangeY(s) ) return bXmax;
		else if(IsRangeX(s)  && IsBoundY1(s)) return bYmax;
		else return Outside;
	}

	// Ends of the computational domain, ghost layers excluded. The base RectGrid
	// xmin()/xmax() include the ghosts, so these are held separately.
	double XminReal() const {return _Xmin_real;}
	double XmaxReal() const {return _Xmax_real;}
	double YminReal() const {return _Ymin_real;}
	double YmaxReal() const {return _Ymax_real;}

	// Number of ghost layers on one side (2, because the film-height equation is fourth order).
	int NumGhostX() const {return _num_ghostX;}
	int NumGhostY() const {return _num_ghostY;}

	// Index range of the computational domain, half-open in each direction.
	int IminReal() const {return _num_ghostX;}
	int ImaxReal() const {return nx() - _num_ghostX;}
	int JminReal() const {return _num_ghostY;}
	int JmaxReal() const {return ny() -_num_ghostY;}

	/**
	 * @brief Convert a site to its linear index, i.e. its row number in the matrix.
	 * @return nx()*iy + ix, counted over the whole grid including ghost layers
	 */
	inline int IL(AnySite s) const { return nx() * iy(s) + ix(s);}

	/**
	 * @brief Fill svec with the 13 stencil sites centered on S0 (cross plus four diagonals).
	 * @param svec  output array, indexed afterwards with the tAdjacent enumerators; must hold at least 13 entries
	 * @param S0    center site
	 *
	 * On a grid with ny() == 1 only the five x-direction entries are written;
	 * the y and diagonal entries are left untouched.
	 */
	template<class S> inline void SetSite(std::vector<Site<S> > &svec, Site<S> S0) const
	{
		assert(svec.size() >= (size_t)_n_stencil );

		svec[Mo] = ShiftXmm(S0);
		svec[mo] = ShiftXm (S0);
		svec[oo] = S0;
		svec[po] = ShiftXp (S0);
		svec[Po] = ShiftXpp(S0);
		if( ny() != 1)
		{
			svec[oM] = ShiftYmm (S0);
			svec[mm] = ShiftXmYm(S0);
			svec[om] = ShiftYm  (S0);
			svec[pm] = ShiftXpYm(S0);
			svec[mp] = ShiftXmYp(S0);
			svec[op] = ShiftYp  (S0);
			svec[pp] = ShiftXpYp(S0);
			svec[oP] = ShiftYpp (S0);
		}
	}

	// Two-point and diagonal neighbor shifts; the base grid only provides one-point shifts.
	template<class S> inline Site<S> ShiftXmm (Site<S> s) const { return ShiftXm(ShiftXm(s)); }
	template<class S> inline Site<S> ShiftXpp (Site<S> s) const { return ShiftXp(ShiftXp(s)); }
	template<class S> inline Site<S> ShiftYmm (Site<S> s) const { return ShiftYm(ShiftYm(s)); }
	template<class S> inline Site<S> ShiftYpp (Site<S> s) const { return ShiftYp(ShiftYp(s)); }
	template<class S> inline Site<S> ShiftXmYm(Site<S> s) const { return ShiftXm(ShiftYm(s)); }
	template<class S> inline Site<S> ShiftXmYp(Site<S> s) const { return ShiftXm(ShiftYp(s)); }
	template<class S> inline Site<S> ShiftXpYm(Site<S> s) const { return ShiftXp(ShiftYm(s)); }
	template<class S> inline Site<S> ShiftXpYp(Site<S> s) const { return ShiftXp(ShiftYp(s)); }

private:
	int _num_ghostX;
	int _num_ghostY;
	double _Xmin_real;
	double _Xmax_real;
	double _Ymin_real;
	double _Ymax_real;

	static const int _n_stencil;

	bool IsRangeX(AnySite s) const { return (_num_ghostX <= ix(s) && ix(s) < nx() -_num_ghostX);}
	bool IsRangeY(AnySite s) const { return (_num_ghostY <= iy(s) && iy(s) < ny() -_num_ghostY);}
	bool IsBoundX0(AnySite s)const { return (ix(s) == _num_ghostX-1);}
	bool IsBoundX1(AnySite s)const { return (ix(s) == nx() -_num_ghostX);}
	bool IsBoundY0(AnySite s)const { return (iy(s) == _num_ghostY-1);}
	bool IsBoundY1(AnySite s)const { return (iy(s) == ny() -_num_ghostY);}

};

///////////////////////////////////////////////////
//--------------------------------------------------------------------------
//  The 3-D computational grid. The concentration field of the ConcZ model
//  lives on it.
//
//  In plane (x, y) it is the same as C_Grid2D, with the film-thickness
//  direction z added. Since z is normalized by the film thickness, the number
//  of divisions does not change as the film thins. There are six boundary
//  faces, so there are more cases to distinguish: the edges (Edge***) and the
//  corners (Corner***).
//--------------------------------------------------------------------------
class C_Grid3D : public RectGrid<3>
{
public:
	//  Whether the site is inside the computational domain, on which face,
	//  edge or corner of the boundary, or outside.
	enum tGhost
	{
		Inside = 0,
		Outside= 1,
		bXmin   = 2,
		bXmax   = 3,
		bYmin   = 4,
		bYmax   = 5,
		bZmin   = 6,
		bZmax   = 7,
		EdgeMMO = 8,
		EdgePMO = 9,
		EdgeMPO = 10,
		EdgePPO = 11,
		EdgeOMM = 12,
		EdgeOPM = 13,
		EdgeOMP = 14,
		EdgeOPP = 15,
		EdgeMOM = 16,
		EdgePOM = 17,
		EdgeMOP = 18,
		EdgePOP = 19,
		CornerMMM = 20,
		CornerPMM = 21,
		CornerMPM = 22,
		CornerPPM = 23,
		CornerMMP = 24,
		CornerPMP = 25,
		CornerMPP = 26,
		CornerPPP = 27
	};

	/** @brief Construct an empty grid with the default two ghost layers per side; call Build() to give it an extent. */
	C_Grid3D()  : _num_ghostX(2), _num_ghostY(2), _num_ghostZ(2) {}

	/**
	 * @brief Wrap an already-sized base grid.
	 * @param g  base grid whose point counts are taken to include the ghost layers already
	 */
	C_Grid3D(const RectGrid<3> &g);
public:
	virtual ~C_Grid3D(){}

	/** @brief Copy the extent, the point counts and the ghost-layer counts from another grid. */
	void operator=(const C_Grid3D &g);

	/**
	 * @brief Build the grid from the [Grid] input block.
	 * @param g  input grid block; Nx, Ny and Nz count interior points only
	 *
	 * Two ghost layers are added on each side in every direction. The z
	 * direction spans the film thickness in normalized coordinates, so the
	 * number of z divisions does not change as the film thins.
	 */
	void Build(const Input::I_Grid &g);

	/**
	 * @brief Classify a site as interior, on one of the six faces, on an edge, in a corner, or outside.
	 * @return Inside, bXmin ... bZmax for the faces, Edge*** where two boundaries meet,
	 *         Corner*** where three meet, or Outside for the outermost ghost layer
	 */
	tGhost CheckGhost(AnySite s) const
	{
		if (IsRangeX(s) && IsRangeY(s) && IsRangeZ(s)) return Inside;
		else if(IsBoundX0(s) && IsBoundY0(s) && IsBoundZ0(s) ) return CornerMMM;
		else if(IsBoundX1(s) && IsBoundY0(s) && IsBoundZ0(s) ) return CornerPMM;
		else if(IsBoundX0(s) && IsBoundY1(s) && IsBoundZ0(s) ) return CornerMPM;
		else if(IsBoundX1(s) && IsBoundY1(s) && IsBoundZ0(s) ) return CornerPPM;
		else if(IsBoundX0(s) && IsBoundY0(s) && IsBoundZ1(s) ) return CornerMMP;
		else if(IsBoundX1(s) && IsBoundY0(s) && IsBoundZ1(s) ) return CornerPMP;
		else if(IsBoundX0(s) && IsBoundY1(s) && IsBoundZ1(s) ) return CornerMPP;
		else if(IsBoundX1(s) && IsBoundY1(s) && IsBoundZ1(s) ) return CornerPPP;
		else if(IsBoundX0(s) && IsBoundY0(s) && IsRangeZ(s)  ) return EdgeMMO;
		else if(IsBoundX1(s) && IsBoundY0(s) && IsRangeZ(s)  ) return EdgePMO;
		else if(IsBoundX0(s) && IsBoundY1(s) && IsRangeZ(s)  ) return EdgeMPO;
		else if(IsBoundX1(s) && IsBoundY1(s) && IsRangeZ(s)  ) return EdgePPO;
		else if(IsRangeX(s) && IsBoundY0(s) && IsBoundZ0(s)  ) return EdgeOMM;
		else if(IsRangeX(s) && IsBoundY1(s) && IsBoundZ0(s)  ) return EdgeOPM;
		else if(IsRangeX(s) && IsBoundY0(s) && IsBoundZ1(s)  ) return EdgeOMP;
		else if(IsRangeX(s) && IsBoundY1(s) && IsBoundZ1(s)  ) return EdgeOPP;
		else if(IsBoundX0(s) && IsRangeY(s) && IsBoundZ0(s)  ) return EdgeMOM;
		else if(IsBoundX1(s) && IsRangeY(s) && IsBoundZ0(s)  ) return EdgePOM;
		else if(IsBoundX0(s) && IsRangeY(s) && IsBoundZ1(s)  ) return EdgeMOP;
		else if(IsBoundX1(s) && IsRangeY(s) && IsBoundZ1(s)  ) return EdgePOP;
		else if(IsBoundX0(s) && IsRangeY(s)  && IsRangeZ(s) ) return bXmin;
		else if(IsBoundX1(s) && IsRangeY(s)  && IsRangeZ(s) ) return bXmax;
		else if(IsRangeX(s)  && IsBoundY0(s) && IsRangeZ(s) ) return bYmin;
		else if(IsRangeX(s)  && IsBoundY1(s) && IsRangeZ(s) ) return bYmax;
		else if(IsRangeX(s)  && IsRangeY(s)  && IsBoundZ1(s)) return bZmin;
		else if(IsRangeX(s)  && IsRangeY(s)  && IsBoundZ1(s)) return bZmax;
		else return Outside;
	}

	// Ends of the computational domain, ghost layers excluded. The z pair is in
	// coordinates normalized by the film thickness, not in the in-plane length scale.
	double XminReal() const {return _Xmin_real;}
	double XmaxReal() const {return _Xmax_real;}
	double YminReal() const {return _Ymin_real;}
	double YmaxReal() const {return _Ymax_real;}
	double ZminReal() const {return _Zmin_real;}
	double ZmaxReal() const {return _Zmax_real;}

	// Number of ghost layers on one side (2 in every direction).
	int NumGhostX() const {return _num_ghostX;}
	int NumGhostY() const {return _num_ghostY;}
	int NumGhostZ() const {return _num_ghostZ;}

	// Index range of the computational domain, half-open in each direction.
	int IminReal() const {return _num_ghostX;}
	int ImaxReal() const {return nx() - _num_ghostX;}
	int JminReal() const {return _num_ghostY;}
	int JmaxReal() const {return ny() -_num_ghostY;}
	int KminReal() const {return _num_ghostZ;}
	int KmaxReal() const {return nz() -_num_ghostZ;}

	/**
	 * @brief Convert a site to its linear index, i.e. its row number in the matrix.
	 * @return nx()*ny()*iz + nx()*iy + ix, counted over the whole grid including ghost layers
	 */
	inline int IL(AnySite s) const { return (nx() * ny() * iz(s)) + (nx() * iy(s)) + ix(s);}

	// Two-point and diagonal neighbor shifts; the base grid only provides one-point shifts.
	template<class S> inline Site<S> ShiftXmm (Site<S> s) const { return ShiftXm(ShiftXm(s)); }
	template<class S> inline Site<S> ShiftXpp (Site<S> s) const { return ShiftXp(ShiftXp(s)); }
	template<class S> inline Site<S> ShiftYmm (Site<S> s) const { return ShiftYm(ShiftYm(s)); }
	template<class S> inline Site<S> ShiftYpp (Site<S> s) const { return ShiftYp(ShiftYp(s)); }
	template<class S> inline Site<S> ShiftXmYm(Site<S> s) const { return ShiftXm(ShiftYm(s)); }
	template<class S> inline Site<S> ShiftXmYp(Site<S> s) const { return ShiftXm(ShiftYp(s)); }
	template<class S> inline Site<S> ShiftXpYm(Site<S> s) const { return ShiftXp(ShiftYm(s)); }
	template<class S> inline Site<S> ShiftXpYp(Site<S> s) const { return ShiftXp(ShiftYp(s)); }

private:
	int _num_ghostX;
	int _num_ghostY;
	int _num_ghostZ;
	double _Xmin_real;
	double _Xmax_real;
	double _Ymin_real;
	double _Ymax_real;
	double _Zmin_real;
	double _Zmax_real;

	bool IsRangeX(AnySite s) const { return (_num_ghostX <= ix(s) && ix(s) < nx() -_num_ghostX);}
	bool IsRangeY(AnySite s) const { return (_num_ghostY <= iy(s) && iy(s) < ny() -_num_ghostY);}
	bool IsRangeZ(AnySite s) const { return (_num_ghostZ <= iz(s) && iz(s) < nz() -_num_ghostZ);}
	bool IsBoundX0(AnySite s)const { return (ix(s) == _num_ghostX-1);}
	bool IsBoundX1(AnySite s)const { return (ix(s) == nx() -_num_ghostX);}
	bool IsBoundY0(AnySite s)const { return (iy(s) == _num_ghostY-1);}
	bool IsBoundY1(AnySite s)const { return (iy(s) == ny() -_num_ghostY);}
	bool IsBoundZ0(AnySite s)const { return (iz(s) == _num_ghostZ-1);}
	bool IsBoundZ1(AnySite s)const { return (iz(s) == nz() -_num_ghostZ);}
};


#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_GRID_INCLUDED_
