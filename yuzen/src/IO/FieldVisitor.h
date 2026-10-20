/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  FieldVisitor.h : base of the visitors over the field class hierarchy
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _FIELD_VISITOR_INCLUDED_
#define _FIELD_VISITOR_INCLUDED_


#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class Field_WM1D;
class Field_WM2D;
class Field_CZ;

//--------------------------------------------------------------------------
//  Base of the visitors over the fields.
//
//  This is the mechanism that keeps the branch on the kind of field (WM1D /
//  WM2D / CZ) on the traversing side rather than in the field classes. The
//  aggregation (FieldInfo) and the probe (PointProbe) derive from it. Adding a
//  new traversal therefore does not touch the field classes either.
//--------------------------------------------------------------------------
class FieldVisitor
{
public:
	FieldVisitor() {}
	virtual ~FieldVisitor() {}

	/**
	 * @brief Handles one field. The overload that runs is chosen by the field,
	 *        through its Accept, so a new traversal is added by implementing
	 *        all three here and nothing in the field classes changes.
	 */
	virtual void Visit(Field_WM1D &) = 0;
	virtual void Visit(Field_WM2D &) = 0;
	virtual void Visit(Field_CZ &) = 0;
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_FIELD_VISITOR_INCLUDED_
