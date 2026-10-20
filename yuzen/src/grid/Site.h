/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Site.h : an index that names one point of the grid. The kind of position
*
*  (the Stagger) is carried in the type, so a mix-up such as using a cell-center
*  point where a grid point is expected is a compile error.
*
*  A site holds a single integer, the "linear index". It runs from 1 to imax;
*  0 means "no site", and imax+1 is the sentinel that ends a sweep.
*
*      for(Zface s = G.First<Zface>(); s != G.Last<Zface>(); ++s) ...
*
*  AnySite is an index with the kind dropped. It exists only to be passed to
*  tests that do not depend on the kind (whether the point is in the ghost
*  layer, for instance); it cannot be used to read or write a field.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _YUZEN_GRID_SITE_INCLUDED_
#define _YUZEN_GRID_SITE_INCLUDED_

#include "Stagger.h"

namespace yuzen { namespace grid {

//--------------------------------------------------------------------------
//  A grid point with the kind dropped.
//
//  A grid point that does not carry "which staggered position" in its type.
//  It is used to pass a point to code that does not depend on the kind of
//  position, such as a test for being on the boundary. A Site<S> converts to
//  it implicitly, but not the other way round (the kind cannot be recovered).
//--------------------------------------------------------------------------
class AnySite
{
public:
	AnySite() : _i(0) {}

	/**
	 * @brief Builds an untyped site from a LINEAR index (not a grid coordinate).
	 *
	 * The argument is the linear index i = 1 + ix + iy*nx + iz*nx*ny, not one of
	 * the grid coordinates ix / iy / iz. The constructor is explicit, so an int
	 * never converts to a site on its own.
	 *
	 * @param i  linear index; 0 means "no site", imax+1 is the loop sentinel
	 */
	explicit AnySite(int i) : _i(i) {}

	/**
	 * @brief Linear index used to address a field array.
	 * @return the linear index, in 1 .. imax for a valid site
	 */
	int Index() const { return _i; }

private:
	int _i;
};


//--------------------------------------------------------------------------
//  A grid point that carries its staggered position in the type.
//
//  It holds a single linear index, but the template argument S records which
//  position it is, so passing a Zface where an Xedge is expected does not
//  compile. A mix-up between staggered positions is hard to notice at run time,
//  so it is stopped by the type system instead.
//
//  S is one of AtNode / AtXedge / ... from Stagger.h. Aliases are given below.
//--------------------------------------------------------------------------
template<class S>
class Site
{
public:
	typedef S stagger_type;

	Site() : _i(0) {}

	/**
	 * @brief Builds a typed site from a LINEAR index (not a grid coordinate).
	 *
	 * The argument is the linear index i = 1 + ix + iy*nx + iz*nx*ny, NOT a grid
	 * coordinate. Passing ix where a linear index is expected compiles and gives
	 * a silently wrong site, so prefer building sites through the grid:
	 * RectGrid::Of<S>(i) from a linear index, RectGrid::At<S>(ix, iy) from grid
	 * coordinates. The constructor is explicit for the same reason.
	 *
	 * @param i  linear index; 0 means "no site", imax+1 is the loop sentinel
	 */
	explicit Site(int i) : _i(i) {}

	/**
	 * @brief Linear index used to address a field array.
	 * @return the linear index, in 1 .. imax for a valid site
	 */
	int  Index() const { return _i; }
	int& Index()       { return _i; }

	/**
	 * @brief Drops the stagger type, giving an AnySite.
	 *
	 * Used to pass the site to code that does not depend on the stagger type
	 * (a ghost-layer test, for instance). The conversion is one way: the type
	 * cannot be recovered from an AnySite.
	 */
	operator AnySite() const { return AnySite(_i); }

	//  Iteration and comparison: both walk the linear index of one stagger type.
	//  Use with RectGrid::First<S>() .. Last<S>(). Sites of different types have
	//  different C++ types, so they cannot be compared with each other.
	//  Walks points of the same kind in order of the linear index.
	//  Meant for looping from RectGrid::First<S>() to Last<S>().
	Site& operator++()    { ++_i; return *this; }
	Site  operator++(int) { Site t(*this); ++_i; return t; }

	//  Every comparison is between linear indices. Points of different kinds
	//  have different types, so they cannot be compared at all.
	bool operator==(const Site &s) const { return _i == s._i; }
	bool operator!=(const Site &s) const { return _i != s._i; }
	bool operator< (const Site &s) const { return _i <  s._i; }
	bool operator<=(const Site &s) const { return _i <= s._i; }
	bool operator> (const Site &s) const { return _i >  s._i; }
	bool operator>=(const Site &s) const { return _i >= s._i; }

private:
	int _i;
};


//  Aliases for the positions used most often. The name says which face or edge
//  the point sits at the center of.
//    Node  : intersection of grid lines   Cell : cell center
//    Xedge : center of an edge along x    Xface: center of an x face (a yz face)
typedef Site<AtNode>  Node;
typedef Site<AtXedge> Xedge;
typedef Site<AtYedge> Yedge;
typedef Site<AtZedge> Zedge;
typedef Site<AtXface> Xface;
typedef Site<AtYface> Yface;
typedef Site<AtZface> Zface;
typedef Site<AtCell>  Cell;

} }	// END: namespace yuzen::grid

#endif //_YUZEN_GRID_SITE_INCLUDED_
