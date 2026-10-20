/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Field.h : a quantity laid out on the grid. It carries the kind of position
*
*  (the Stagger) in its type, so indexing a cell-center quantity with a grid
*  point is a compile error.
*
*  Inside it is a single array addressed by the linear index. Entry 0 catches
*  the invalid value and one entry at the end is spare room for the end-of-scan
*  sentinel; neither is ever read or written.
*
*  For scalars and vectors, where the quantity lives is fixed per dimension.
*
*      scalar     1D: Xedge   2D: Zface   3D: Cell      (cell center)
*      vector   x 1D: Node    2D: Yedge   3D: Xface     (face with x normal)
*             y          -        Xedge       Yface
*             z          -            -       Zface
*
*  Each of these is just the staggered-layout convention itself: scalars inside
*  the cell, vectors on the cell boundary along its normal.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _YUZEN_GRID_FIELD_INCLUDED_
#define _YUZEN_GRID_FIELD_INCLUDED_

#include <cstddef>
#include <vector>

#include "RectGrid.h"
#include "FirstTouch.h"

namespace yuzen { namespace grid {

//--------------------------------------------------------------------------
//  A field on the grid.
//
//  S records which staggered position the field lives on, so indexing an Xedge
//  field with a Zface point does not compile. Such a mix-up is hard to notice
//  at run time, so it is stopped by the type system.
//
//  The grid is only held by reference, never owned, so the grid must outlive
//  the field. In practice the model owns one grid and the fields point at it.
//--------------------------------------------------------------------------
template<class GRID, class S, class T = double>
class Field
{
public:
	typedef GRID grid_type;
	typedef S    stagger_type;
	typedef T    value_type;

	Field() : _G(0) {}

	/**
	 * @brief Builds the field on grid g straight away; same as Field() + Build().
	 *
	 * The grid is only referenced, never copied or owned, so g must outlive the
	 * field.
	 */
	explicit Field(const GRID &g, T val = T()) : _G(0) { Build(g, val); }

	/**
	 * @brief Binds the field to a grid, allocates storage and fills it with val.
	 *
	 * Allocates imax()+2 entries, not imax(): the linear index starts at 1, index
	 * 0 is the invalid-site sentinel and imax()+1 the end-of-scan sentinel.
	 * Neither sentinel is meant to be read or written.
	 *
	 * @param g    grid to bind to; held by reference only, so it must outlive
	 *             this field
	 * @param val  initial value for every point
	 */
	//  Binds the grid, allocates storage and fills it with val.
	//  It allocates imax()+2 entries, because the linear index starts at 1 and
	//  0 is the sentinel for an invalid point, imax()+1 the end of a scan.
	void Build(const GRID &g, T val = T())
	{
		_G = &g;
		_body.assign(static_cast<std::size_t>(g.imax()) + 2, val);
		_body[0] = T();
	}

	/**
	 * @brief Whether Build() has run, i.e. whether the field has a grid and storage.
	 */
	//  Whether Build() has been called
	bool Constructed() const { return _G != 0; }

	/**
	 * @brief The grid this field is bound to.
	 * @return reference to the grid; calling this before Build() dereferences a
	 *         null pointer
	 */
	//  The grid this field is bound to
	const GRID& Grid() const { return *_G; }

	/**
	 * @brief Reads or writes the value at a site. This is the normal access path.
	 *
	 * Takes a typed Site\<S\>, so a field living on one stagger type cannot be
	 * indexed with a site of another: the mix-up is a compile error. Contrast
	 * operator[], which takes a raw LINEAR index and checks nothing.
	 *
	 * @param s  site of the same stagger type as the field
	 */
	//  Index by site. This is the normal way to use a field.
	const T& operator()(Site<S> s) const { return _body[s.Index()]; }
	      T& operator()(Site<S> s)       { return _body[s.Index()]; }

	/**
	 * @brief Reads or writes by raw LINEAR index, bypassing the type check.
	 *
	 * For code that walks the storage without caring about the stagger type,
	 * such as file I/O. The argument is the linear index (1 .. imax()), never a
	 * grid coordinate, and nothing verifies it: use operator() with a Site
	 * wherever the site type is known.
	 *
	 * @param i  linear index in 1 .. imax(); 0 and imax()+1 are sentinels
	 */
	//  Direct access by linear index (for I/O and other kind-agnostic code)
	const T& operator[](int i) const { return _body[i]; }
	      T& operator[](int i)       { return _body[i]; }

	/**
	 * @brief Sets every point of the field to val.
	 *
	 * Leaves entry 0, the invalid-site sentinel, untouched.
	 */
	//  Sets every point to val (entry 0, the sentinel, is left alone)
	void Fill(T val)
	{
		for(std::size_t k = 1; k < _body.size(); ++k) _body[k] = val;
	}

private:
	const GRID    *_G;
	std::vector<T, yuzen::util::FirstTouchAllocator<T> > _body;
};


//--------------------------------------------------------------------
//  Scalar fields
//
//  Where a scalar lives is fixed per dimension: in 1D the center of an edge
//  along x, in 2D the center of a z face (= the cell center of the plane), in
//  3D the cell center. The film thickness h and the concentration c live there.
//--------------------------------------------------------------------
template<class T = double> using ScalarField1D = Field<RectGrid<1>, AtXedge, T>;
template<class T = double> using ScalarField2D = Field<RectGrid<2>, AtZface, T>;
template<class T = double> using ScalarField3D = Field<RectGrid<3>, AtCell,  T>;


//--------------------------------------------------------------------
//  Vector fields (each component lives in a different place)
//
//  The layout is staggered, so the x component and the y component live in
//  different places. Gathering them into one type spares the caller from having
//  to keep in mind that the components are fields of different types. x and y
//  are public, so they are indexed directly, as in v.x(nd).
//--------------------------------------------------------------------
template<class T = double>
class VectorField1D
{
public:
	//  The x component. In 1D the velocity lives on the grid lines (Node).
	Field<RectGrid<1>, AtNode, T> x;

	/**
	 * @brief Allocates every component on grid g and fills it with val.
	 *
	 * Constructed() and Grid() just below report on the x component. All
	 * components are always built in this one call, so x stands for all of them.
	 */
	//  Allocates every component in one call
	void Build(const RectGrid<1> &g, T val = T()) { x.Build(g, val); }
	/** @brief Whether Build() has run, as reported by the x component. */
	bool Constructed() const { return x.Constructed(); }
	/** @brief The grid the components are bound to. */
	const RectGrid<1>& Grid() const { return x.Grid(); }
};

template<class T = double>
class VectorField2D
{
public:
	//  The x component lives at the center of an edge along y, the y component
	//  at the center of an edge along x. Each sits at the center of the face
	//  that component crosses.
	Field<RectGrid<2>, AtYedge, T> x;
	Field<RectGrid<2>, AtXedge, T> y;

	/**
	 * @brief Allocates both components on grid g and fills them with val.
	 *
	 * Constructed() and Grid() just below report on the x component, which stands
	 * for all components because they are always built together.
	 */
	void Build(const RectGrid<2> &g, T val = T()) { x.Build(g, val); y.Build(g, val); }
	/** @brief Whether Build() has run, as reported by the x component. */
	bool Constructed() const { return x.Constructed(); }
	/** @brief The grid the components are bound to. */
	const RectGrid<2>& Grid() const { return x.Grid(); }
};

template<class T = double>
class VectorField3D
{
public:
	//  Each component lives at the center of the face it crosses.
	Field<RectGrid<3>, AtXface, T> x;
	Field<RectGrid<3>, AtYface, T> y;
	Field<RectGrid<3>, AtZface, T> z;

	/**
	 * @brief Allocates all three components on grid g and fills them with val.
	 *
	 * Constructed() and Grid() just below report on the x component, which stands
	 * for all components because they are always built together.
	 */
	void Build(const RectGrid<3> &g, T val = T()) { x.Build(g, val); y.Build(g, val); z.Build(g, val); }
	/** @brief Whether Build() has run, as reported by the x component. */
	bool Constructed() const { return x.Constructed(); }
	/** @brief The grid the components are bound to. */
	const RectGrid<3>& Grid() const { return x.Grid(); }
};

} }	// END: namespace yuzen::grid

#endif //_YUZEN_GRID_FIELD_INCLUDED_
