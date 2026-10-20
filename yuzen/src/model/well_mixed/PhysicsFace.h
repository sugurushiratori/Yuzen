/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  PhysicsFace.h : the identifiers of the "faces" a physics term refers to
*
*  The boundaries between cells (Xm, Xp, Ym, Yp), plus OO, which denotes the
*  cell itself rather than a face. The 1-D case does not use Ym / Yp.
*
*  These identifiers are kept in a header of their own because both the grid
*  traits (GridTraits.h) and the physics-term base (PhysicsTermBase.h) refer
*  to them; having the traits include the base would be circular.
*
*  The values of this enum are never used as array subscripts. They only
*  select a face in a switch, so neither their order nor their values carry
*  any meaning.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _PHYSICS_FACE_INCLUDED_
#define _PHYSICS_FACE_INCLUDED_

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//--------------------------------------------------------------------------
//  A holder whose only purpose is to carry the face specifiers.
//  Xm / Xp / Ym / Yp are the faces shared with the neighboring cells;
//  OO is the cell itself (not a face).
//--------------------------------------------------------------------------
struct PhysicsFace
{
	/**
	 * @brief Identifier of the face a physics term is currently working on.
	 *
	 * Xm / Xp / Ym / Yp are the faces shared with the neighboring cells;
	 * OO denotes the cell itself rather than a face. The 1-D case never
	 * uses Ym / Yp. The values are never used as array subscripts, so
	 * neither their order nor their numeric values carry meaning.
	 */
	enum tEdge{ Xm, Xp, Ym, Yp, OO };
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_PHYSICS_FACE_INCLUDED_
