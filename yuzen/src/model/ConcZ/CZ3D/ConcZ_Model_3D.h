/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ConcZ_Model_3D.h : driver for the 3-D version of ConcZ (unimplemented skeleton)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _CONCZ_MODEL_3D_INCLUDED_
#define _CONCZ_MODEL_3D_INCLUDED_

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif


/////////////////////////////////////////////////////////////////////////////
/**
 * @brief Placeholder for the 3-D version of the ConcZ model.
 *
 * Nothing is implemented yet. The class does not derive from ISolverModel and
 * no input combination selects it; it only reserves the name. The constructor
 * declared here has no definition, so the class cannot be instantiated.
 */
class ConcZ_Model_3D
{
public:
	ConcZ_Model_3D();
	~ConcZ_Model_3D(){}

};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_CONCZ_MODEL_INCLUDED_
