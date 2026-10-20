/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Stagger.h : encodes the "kind of position" on a staggered grid in the type
*
*  When a rectangular box is cut by the grid lines (ix, iy, iz), a quantity can
*  sit in only two places along each axis:
*
*      0 : on the grid line              (x = xmin + ix*dx)
*      1 : halfway between two grid lines (x = xmin + ix*dx + dx/2)
*
*  So three bits cover every kind of position. Those three bits are the Stagger.
*
*      Stagger<0,0,0>  grid point            Node
*      Stagger<1,0,0>  edge along x          Xedge
*      Stagger<0,1,0>  edge along y          Yedge
*      Stagger<0,0,1>  edge along z          Zedge
*      Stagger<0,1,1>  face with x normal    Xface
*      Stagger<1,0,1>  face with y normal    Yface
*      Stagger<1,1,0>  face with z normal    Zface
*      Stagger<1,1,1>  cell center           Cell
*
*  Moving half a cell along an axis flips the bit of that axis. That is, the way
*  a half-cell step changes the kind of position is written as a single bit flip.
*  There is no need to write out a table of offsets by hand.
*
*  In 1D only the x bit carries meaning (Node and Xedge); in 2D the x and y bits
*  carry meaning (Node, Xedge, Yedge, Zface).
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _YUZEN_GRID_STAGGER_INCLUDED_
#define _YUZEN_GRID_STAGGER_INCLUDED_

namespace yuzen { namespace grid {

//  Axis numbers
enum { AxisX = 0, AxisY = 1, AxisZ = 2 };

template<int SX, int SY, int SZ>
struct Stagger
{
	enum { sx = SX, sy = SY, sz = SZ };

	/**
	 * @brief Stagger bit of one axis: 0 = on the grid line, 1 = half a cell off it.
	 * @param a  axis number, one of AxisX / AxisY / AxisZ
	 * @return   the bit of that axis (SX, SY or SZ)
	 */
	static int Bit(int a) { return (a == AxisX) ? SX : (a == AxisY) ? SY : SZ; }
};

typedef Stagger<0,0,0> AtNode;
typedef Stagger<1,0,0> AtXedge;
typedef Stagger<0,1,0> AtYedge;
typedef Stagger<0,0,1> AtZedge;
typedef Stagger<0,1,1> AtXface;
typedef Stagger<1,0,1> AtYface;
typedef Stagger<1,1,0> AtZface;
typedef Stagger<1,1,1> AtCell;

//--------------------------------------------------------------------
//  Kind of position reached by a half-cell step along axis A (the bit of A flips)
//--------------------------------------------------------------------
template<class S, int A> struct HalfStep;

template<class S> struct HalfStep<S, AxisX>
{ typedef Stagger<1 - S::sx, S::sy, S::sz> type; };

template<class S> struct HalfStep<S, AxisY>
{ typedef Stagger<S::sx, 1 - S::sy, S::sz> type; };

template<class S> struct HalfStep<S, AxisZ>
{ typedef Stagger<S::sx, S::sy, 1 - S::sz> type; };

} }	// END: namespace yuzen::grid

#endif //_YUZEN_GRID_STAGGER_INCLUDED_
