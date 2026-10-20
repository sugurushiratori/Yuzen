/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  RectGrid.h : uniform rectangular grid (1D / 2D / 3D)
*
*  = Indices
*  A point of the grid is named by a single linear index.
*
*      i = 1 + ix*1 + iy*nx + iz*nx*ny        (1 <= i <= nx*ny*nz)
*
*  i = 0 is the invalid value and i = imax+1 the end-of-scan sentinel.
*  ix / iy / iz are 0-origin and are that formula solved backwards.
*
*  = Periodicity
*  Neighbors in the index space are always closed periodically. Non-periodic
*  boundary conditions are expressed by Yuzen placing ghost layers; the grid
*  itself knows nothing about them.
*
*  = Stepping
*  A whole-cell step (ShiftXm and friends) does not change the kind of position.
*  A half-cell step (HalfXm and friends) flips the bit of that axis, so the kind
*  does change. Whether the linear index moves is decided by the bit of the
*  position you are standing on.
*
*      +half cell : bit 0 keeps the index, bit 1 gives +1
*      -half cell : bit 1 keeps the index, bit 0 gives -1
*
*  The easiest way to see this: stepping +half a cell from a point on a grid
*  line (bit 0) still leaves you inside the same cell.
*
*  = Coordinates
*      x = xmin + ix*dx + (dx/2 if the x bit is 1)
*  and likewise for y and z.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _YUZEN_GRID_RECTGRID_INCLUDED_
#define _YUZEN_GRID_RECTGRID_INCLUDED_

#include <vector>
#include <cassert>

#include "Site.h"

namespace yuzen { namespace grid {

template<int DIM>
class RectGrid
{
public:
	enum { Dim = DIM };

	RectGrid()
	 : _nx(0), _ny(1), _nz(1), _imax(0),
	   _xmin(0), _xmax(0), _ymin(0), _ymax(0), _zmin(0), _zmax(0),
	   _dx(0), _dy(0), _dz(0), _dx2(0), _dy2(0), _dz2(0) {}

	//----------------------------------------------------------------
	//  Construction
	//
	//  You give only the extent and the number of divisions. Calling a Build()
	//  that does not match DIM stops the build at the static_assert.
	//
	//  The nx passed here is the total number of points the grid holds, not the
	//  number of cells of the computational domain. How many ghost layers to
	//  place is a concern of the layer above (model/Grids.h), and this layer
	//  knows nothing about it.
	//----------------------------------------------------------------
	/**
	 * @brief Builds a 1D grid over [xmin, xmax] with nx equal divisions.
	 * @param xmin  low end of the x extent
	 * @param xmax  high end of the x extent
	 * @param nx  number of points the grid holds along x, NOT the number of
	 *            computational cells. Ghost layers are added by the caller
	 *            (model/Grids.h); the grid itself knows nothing about them.
	 */
	void Build(double xmin, double xmax, int nx)
	{
		static_assert(DIM == 1, "this Build() is for 1-D grids only");
		SetExtentX(xmin, xmax, nx);
		Allocate();
	}

	/**
	 * @brief Builds a 2D grid over [xmin, xmax] x [ymin, ymax].
	 * @param xmin  low end of the x extent
	 * @param xmax  high end of the x extent
	 * @param nx  number of points along x (not cells; see the 1D overload)
	 * @param ymin  low end of the y extent
	 * @param ymax  high end of the y extent
	 * @param ny  number of points along y
	 */
	void Build(double xmin, double xmax, int nx,
	           double ymin, double ymax, int ny)
	{
		static_assert(DIM == 2, "this Build() is for 2-D grids only");
		SetExtentX(xmin, xmax, nx);
		SetExtentY(ymin, ymax, ny);
		Allocate();
	}

	/**
	 * @brief Builds a 3D grid over [xmin, xmax] x [ymin, ymax] x [zmin, zmax].
	 * @param xmin  low end of the x extent
	 * @param xmax  high end of the x extent
	 * @param nx  number of points along x (not cells; see the 1D overload)
	 * @param ymin  low end of the y extent
	 * @param ymax  high end of the y extent
	 * @param ny  number of points along y
	 * @param zmin  low end of the z extent
	 * @param zmax  high end of the z extent
	 * @param nz  number of points along z
	 */
	void Build(double xmin, double xmax, int nx,
	           double ymin, double ymax, int ny,
	           double zmin, double zmax, int nz)
	{
		static_assert(DIM == 3, "this Build() is for 3-D grids only");
		SetExtentX(xmin, xmax, nx);
		SetExtentY(ymin, ymax, ny);
		SetExtentZ(zmin, zmax, nz);
		Allocate();
	}

	//----------------------------------------------------------------
	//  Shape
	//----------------------------------------------------------------
	//  Shape queries: dimension, point count per axis, extent per axis and mesh
	//  spacing per axis. An axis the grid does not use reports a point count of 1
	//  and an extent / spacing of 0, so the same 1D / 2D / 3D code compiles and
	//  the loop over that axis simply runs once.
	int Dimension() const { return DIM; }

	//  Number of points per axis. An unused axis reports 1, so the same code
	//  serves 1D, 2D and 3D (the loop over that axis just runs once).
	int nx() const { return _nx; }

	int ny() const { return (DIM >= 2) ? _ny : 1; }
	int nz() const { return (DIM >= 3) ? _nz : 1; }

	double xmin() const { return _xmin; }
	double xmax() const { return _xmax; }
	double ymin() const { return (DIM >= 2) ? _ymin : 0.0; }
	double ymax() const { return (DIM >= 2) ? _ymax : 0.0; }
	double zmin() const { return (DIM >= 3) ? _zmin : 0.0; }
	double zmax() const { return (DIM >= 3) ? _zmax : 0.0; }

	double dx() const { return _dx; }
	double dy() const { return (DIM >= 2) ? _dy : 0.0; }
	double dz() const { return (DIM >= 3) ? _dz : 0.0; }

	/**
	 * @brief First and last valid LINEAR index of the grid.
	 *
	 * The linear index range starts at 1, not at 0. Index 0 is reserved as the
	 * "invalid site" sentinel and imax()+1 as the end-of-loop sentinel, and the
	 * tables built by Allocate() are sized for both. All points are covered by
	 *
	 *     for(int i = G.imin(); i <= G.imax(); ++i) ...
	 *
	 * A loop written from 0, or up to imax()-1, is off by one.
	 *
	 * @return imin() is always 1; imax() is nx*ny*nz
	 */
	int imin() const { return 1;     }
	int imax() const { return _imax; }

	//----------------------------------------------------------------
	//  Indices
	//----------------------------------------------------------------
	//  Index conversions. Two different integers can describe the same point and
	//  they must never be mixed up:
	//    - the LINEAR index: one number, 1 .. imax(), the thing a Site carries
	//      and the thing a field array is addressed with;
	//    - the GRID COORDINATES ix / iy / iz: one 0-origin counter per axis.
	//  Of<>() converts from the first, At<>() from the second, and Index() goes
	//  back to the first. Both take plain ints, so the compiler cannot tell them
	//  apart - see the blocks on Of<>() and At<>() below.

	//  Linear index of a site, for indexing a field array.
	//  The linear index used to address a field array
	int Index(AnySite s) const { return s.Index(); }

	/**
	 * @brief Makes a site of type ST from a LINEAR index.
	 *
	 * Counterpart of At<>(). Of<>() takes the linear index i (1 .. imax()) that a
	 * Site already carries; At<>() takes the grid coordinates ix, iy, iz. Both
	 * take plain ints, so handing a grid coordinate to Of<>(), or a linear index
	 * to At<>(), compiles cleanly and silently addresses the wrong point. That
	 * mix-up is the known source of off-by-one bugs here.
	 *
	 *     Cell c = G.Of<Cell>(s.Index());   // from a linear index
	 *     Cell c = G.At<Cell>(ix, iy);      // from grid coordinates
	 *
	 * @param i  linear index in 1 .. imax(); it is neither checked nor wrapped
	 * @return   site of type ST sitting at that linear index
	 */
	template<class ST> ST Of(int i) const { return ST(i); }

	/**
	 * @brief Makes a site of type ST from GRID COORDINATES (ix, iy, iz), 0-origin.
	 *
	 * Counterpart of Of<>(): the arguments here are per-axis counters, not the
	 * linear index. Coordinates outside 0 .. n-1 wrap around periodically, so
	 * At<>() always returns a usable site and can never be used to detect a
	 * coordinate that has run off the end of the axis. The overload has to match
	 * DIM, otherwise the static_assert stops the build.
	 *
	 * The 2-D and 3-D overloads take the further coordinates iy and iz, on the
	 * same terms.
	 *
	 * @param ix  grid coordinate along x, 0-origin, wrapped if out of range
	 * @return    site of type ST at that grid coordinate
	 */
	//  Makes a site from grid coordinates (ix, iy, iz). Out-of-range
	//  coordinates wrap around periodically.
	template<class ST> ST At(int ix) const
	{ static_assert(DIM == 1, "1-D grids only"); return ST(Linear(ix, 0, 0)); }

	template<class ST> ST At(int ix, int iy) const
	{ static_assert(DIM == 2, "2-D grids only"); return ST(Linear(ix, iy, 0)); }

	template<class ST> ST At(int ix, int iy, int iz) const
	{ static_assert(DIM == 3, "3-D grids only"); return ST(Linear(ix, iy, iz)); }

	/**
	 * @brief First site and one-past-the-last site of a scan over type ST.
	 *
	 * Last<ST>() is a sentinel at linear index imax()+1. It is past the end of
	 * the grid and must never be dereferenced or used to read or write a field;
	 * it only terminates the loop.
	 *
	 *     for(Zface s = G.First<Zface>(); s != G.Last<Zface>(); ++s) ...
	 *
	 * Use != against Last<ST>(), not <=.
	 */
	//  First and end of a scan (the end is a sentinel; it cannot be read or written)
	template<class ST> ST First() const { return ST(imin());   }
	template<class ST> ST Last()  const { return ST(_imax + 1); }

	//  Grid coordinates of a site (the inverse of At<>()). The values are looked
	//  up in a table built by Allocate(), so no division happens in inner loops.
	//  A site outside the valid range reports 0 on every axis, which is a
	//  sentinel and not a real coordinate.
	//  Looks up the grid coordinate of a site. The values are baked into a
	//  table, so no division happens. A site out of range gets 0 (a sentinel).
	int ix(AnySite s) const { return Valid(s.Index()) ? _ix[s.Index()] : 0; }

	int iy(AnySite s) const
	{
		static_assert(DIM >= 2, "2-D grids and above only");
		return Valid(s.Index()) ? _iy[s.Index()] : 0;
	}

	int iz(AnySite s) const
	{
		static_assert(DIM >= 3, "3-D grids only");
		return Valid(s.Index()) ? _iz[s.Index()] : 0;
	}

	//----------------------------------------------------------------
	//  Whole-cell step (the kind of position does not change)
	//
	//  Moves to a point of the same kind, as from an Xedge to the next Xedge.
	//  Neighbors are baked into a table too, so no index arithmetic happens.
	//  A step that leaves the domain wraps around periodically (how boundaries
	//  behave is decided by BoundaryCondition in the layer above).
	//----------------------------------------------------------------
	//  Whole-cell neighbor shifts, m = minus (-1 cell), p = plus (+1 cell).
	//  The site type is unchanged: an Xedge shifts to the next Xedge. Neighbors
	//  come from a table, so no index arithmetic is done here. A shift that
	//  leaves the domain wraps around periodically; non-periodic boundaries are
	//  the caller's business (ghost layers plus BoundaryCondition).
	template<class S> Site<S> ShiftXm(Site<S> s) const { return Site<S>(_xm[s.Index()]); }
	template<class S> Site<S> ShiftXp(Site<S> s) const { return Site<S>(_xp[s.Index()]); }

	template<class S> Site<S> ShiftYm(Site<S> s) const
	{ static_assert(DIM >= 2, "2-D grids and above only"); return Site<S>(_ym[s.Index()]); }
	template<class S> Site<S> ShiftYp(Site<S> s) const
	{ static_assert(DIM >= 2, "2-D grids and above only"); return Site<S>(_yp[s.Index()]); }

	template<class S> Site<S> ShiftZm(Site<S> s) const
	{ static_assert(DIM >= 3, "3-D grids only"); return Site<S>(_zm[s.Index()]); }
	template<class S> Site<S> ShiftZp(Site<S> s) const
	{ static_assert(DIM >= 3, "3-D grids only"); return Site<S>(_zp[s.Index()]); }

	//----------------------------------------------------------------
	//  Half-cell step (the bit of that axis flips, so the kind of position changes)
	//
	//  This is the heart of the staggered grid. The kind of position is one bit
	//  per axis (Stagger.h): 0 means on the grid line, 1 means half a cell off
	//  it. Stepping half a cell is flipping that bit, so the return type is
	//  decided by HalfStep<S, Axis>::type. For example
	//      Cell  --HalfXm--> Yface   (only the x bit drops)
	//      Zface --HalfYm--> Node
	//  Whether the index moves is decided by the same bit:
	//      -half cell : bit 1 stays put, bit 0 goes one index back
	//      +half cell : bit 0 stays put, bit 1 goes one index forward
	//  A mix-up does not compile, because the types do not match.
	//----------------------------------------------------------------
	/**
	 * @brief Moves half a cell along an axis; this CHANGES the site type.
	 *
	 * Documents the whole HalfXm / HalfXp / HalfYm / HalfYp / HalfZm / HalfZp
	 * family (m = minus, p = plus). A site type is one stagger bit per axis
	 * (Stagger.h): 0 on the grid line, 1 half a cell off it. Stepping half a cell
	 * flips the bit of that axis, so the RETURN TYPE is not the argument type but
	 * HalfStep<S, Axis>::type. For example
	 *
	 *     Cell  --HalfXm--> Yface     (the x bit drops)
	 *     Zface --HalfYm--> Node
	 *
	 * Whether the linear index moves is decided by the same bit:
	 *     -half step: bit 1 stays put, bit 0 goes one index back
	 *     +half step: bit 0 stays put, bit 1 goes one index forward
	 * A point on a grid line (bit 0) stepped by +half a cell therefore stays
	 * inside the same cell and keeps its index.
	 *
	 * Assigning the result to the wrong site type does not compile, which is the
	 * point: a staggered mix-up is hard to spot at run time.
	 *
	 * @param s  site to step from
	 * @return   site of type HalfStep<S, Axis>::type, half a cell away
	 */
	template<class S> Site<typename HalfStep<S, AxisX>::type> HalfXm(Site<S> s) const
	{
		typedef typename HalfStep<S, AxisX>::type R;
		return Site<R>(S::sx ? s.Index() : _xm[s.Index()]);
	}
	template<class S> Site<typename HalfStep<S, AxisX>::type> HalfXp(Site<S> s) const
	{
		typedef typename HalfStep<S, AxisX>::type R;
		return Site<R>(S::sx ? _xp[s.Index()] : s.Index());
	}

	template<class S> Site<typename HalfStep<S, AxisY>::type> HalfYm(Site<S> s) const
	{
		static_assert(DIM >= 2, "2-D grids and above only");
		typedef typename HalfStep<S, AxisY>::type R;
		return Site<R>(S::sy ? s.Index() : _ym[s.Index()]);
	}
	template<class S> Site<typename HalfStep<S, AxisY>::type> HalfYp(Site<S> s) const
	{
		static_assert(DIM >= 2, "2-D grids and above only");
		typedef typename HalfStep<S, AxisY>::type R;
		return Site<R>(S::sy ? _yp[s.Index()] : s.Index());
	}

	template<class S> Site<typename HalfStep<S, AxisZ>::type> HalfZm(Site<S> s) const
	{
		static_assert(DIM >= 3, "3-D grids only");
		typedef typename HalfStep<S, AxisZ>::type R;
		return Site<R>(S::sz ? s.Index() : _zm[s.Index()]);
	}
	template<class S> Site<typename HalfStep<S, AxisZ>::type> HalfZp(Site<S> s) const
	{
		static_assert(DIM >= 3, "3-D grids only");
		typedef typename HalfStep<S, AxisZ>::type R;
		return Site<R>(S::sz ? _zp[s.Index()] : s.Index());
	}

	//----------------------------------------------------------------
	//  Coordinates
	//
	//  Returns the position already offset by half a cell according to the kind
	//  of point. At the same grid coordinate, Node and Cell differ by dx/2.
	//----------------------------------------------------------------
	//  Physical coordinate of a site. The half-cell offset of the site type is
	//  already included, so the same grid coordinate gives different values for
	//  different types: X(Node) and X(Cell) differ by dx/2. Do not add the
	//  stagger offset again in the caller.
	template<class S> double X(Site<S> s) const
	{ return _xmin + ix(s) * _dx + (S::sx ? _dx2 : 0.0); }

	template<class S> double Y(Site<S> s) const
	{
		static_assert(DIM >= 2, "2-D grids and above only");
		return _ymin + iy(s) * _dy + (S::sy ? _dy2 : 0.0);
	}

	template<class S> double Z(Site<S> s) const
	{
		static_assert(DIM >= 3, "3-D grids only");
		return _zmin + iz(s) * _dz + (S::sz ? _dz2 : 0.0);
	}

	/**
	 * @brief Physical coordinate straight from GRID COORDINATES, without a site.
	 *
	 * Covers XOf / YOf / ZOf. The site type is given explicitly as ST (Zface,
	 * Cell, ...) and supplies the half-cell offset, exactly as X(Site) does; the
	 * int arguments are grid coordinates ix / iy / iz, never a linear index.
	 * Each function also accepts the arguments of the axes it ignores so that a
	 * caller holding (ix, iy, iz) can pass all three to any of them.
	 *
	 * The coordinates are given as ix, iy and iz, all 0-origin; each function
	 * names only the axis it reads and leaves the others unnamed.
	 *
	 * @param ix  grid coordinate along x, 0-origin
	 * @return    coordinate of that point, including the stagger offset of ST
	 */
	//  Straight from grid coordinates. ST is the type of the point (Zface etc.).
	//  The arguments of the unused axes are accepted as well so that the caller
	//  can pass (ix, iy, iz) through unchanged.
	template<class ST> double XOf(int ix) const                 { return CoordX<ST>(ix); }
	template<class ST> double XOf(int ix, int) const            { return CoordX<ST>(ix); }
	template<class ST> double XOf(int ix, int, int) const       { return CoordX<ST>(ix); }

	template<class ST> double YOf(int, int iy) const            { return CoordY<ST>(iy); }
	template<class ST> double YOf(int, int iy, int) const       { return CoordY<ST>(iy); }

	template<class ST> double ZOf(int, int, int iz) const       { return CoordZ<ST>(iz); }

protected:
	template<class ST> double CoordX(int ix) const
	{ return _xmin + ix * _dx + (ST::stagger_type::sx ? _dx2 : 0.0); }

	template<class ST> double CoordY(int iy) const
	{
		static_assert(DIM >= 2, "2-D grids and above only");
		return _ymin + iy * _dy + (ST::stagger_type::sy ? _dy2 : 0.0);
	}

	template<class ST> double CoordZ(int iz) const
	{
		static_assert(DIM >= 3, "3-D grids only");
		return _zmin + iz * _dz + (ST::stagger_type::sz ? _dz2 : 0.0);
	}

	//  Periodic wrap-around
	static int Wrap(int i, int n)
	{
		if(0 <= i && i <= n - 1) return i;
		if(i == n ) return 0;
		if(i == -1) return n - 1;
		if(i >  n ) return i % n;
		return i % n + n;			// i < -1
	}

	//  Converts grid coordinates (ix, iy, iz) into a linear index. Every axis is
	//  wrapped by Wrap(), so a valid index always comes back even for an
	//  out-of-range coordinate. It starts at 1 because 0 is reserved as the
	//  sentinel for the invalid value.
	int Linear(int ix, int iy, int iz) const
	{
		return 1 + Wrap(ix, _nx)
		         + Wrap(iy, _ny) * _nx
		         + Wrap(iz, _nz) * _nx * _ny;
	}

	//  Whether the linear index points at a real point. 0 (the invalid value)
	//  and imax+1 (the end of a scan) are sentinels, so they give false.
	bool Valid(int i) const { return (imin() <= i && i <= _imax); }

private:
	void SetExtentX(double lo, double hi, int n)
	{
		assert(n > 0);
		_xmin = lo; _xmax = hi; _nx = n;
		_dx = (hi - lo) / n;
		_dx2 = _dx * 0.5;
	}
	void SetExtentY(double lo, double hi, int n)
	{
		assert(n > 0);
		_ymin = lo; _ymax = hi; _ny = n;
		_dy = (hi - lo) / n;
		_dy2 = _dy * 0.5;
	}
	void SetExtentZ(double lo, double hi, int n)
	{
		assert(n > 0);
		_zmin = lo; _zmax = hi; _nz = n;
		_dz = (hi - lo) / n;
		_dz2 = _dz * 0.5;
	}

	//  Bakes neighbors and grid coordinates into tables (so that no division
	//  happens in the inner loops)
	void Allocate()
	{
		_imax = _nx * _ny * _nz;
		const int len = _imax + 2;		// room for 0 (invalid) and imax+1 (sentinel)

		_ix.assign(len, 0);
		_xm.assign(len, 0);
		_xp.assign(len, 0);
		if(DIM >= 2) { _iy.assign(len, 0); _ym.assign(len, 0); _yp.assign(len, 0); }
		if(DIM >= 3) { _iz.assign(len, 0); _zm.assign(len, 0); _zp.assign(len, 0); }

		for(int iz = 0; iz < _nz; ++iz)
		for(int iy = 0; iy < _ny; ++iy)
		for(int ix = 0; ix < _nx; ++ix)
		{
			const int i = Linear(ix, iy, iz);
			_ix[i] = ix;
			_xm[i] = Linear(ix - 1, iy, iz);
			_xp[i] = Linear(ix + 1, iy, iz);
			if(DIM >= 2)
			{
				_iy[i] = iy;
				_ym[i] = Linear(ix, iy - 1, iz);
				_yp[i] = Linear(ix, iy + 1, iz);
			}
			if(DIM >= 3)
			{
				_iz[i] = iz;
				_zm[i] = Linear(ix, iy, iz - 1);
				_zp[i] = Linear(ix, iy, iz + 1);
			}
		}
	}

	int _nx, _ny, _nz;
	int _imax;

	double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;
	double _dx, _dy, _dz;
	double _dx2, _dy2, _dz2;

	std::vector<int> _ix, _iy, _iz;
	std::vector<int> _xm, _xp;
	std::vector<int> _ym, _yp;
	std::vector<int> _zm, _zp;
};

typedef RectGrid<1> RectGrid1D;
typedef RectGrid<2> RectGrid2D;
typedef RectGrid<3> RectGrid3D;

} }	// END: namespace yuzen::grid

#endif //_YUZEN_GRID_RECTGRID_INCLUDED_
