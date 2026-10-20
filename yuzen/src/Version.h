/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Version.h : definition of the version number
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _DEF_VERSION_INCLUDED_
#define _DEF_VERSION_INCLUDED_


////////////////////////////////////////
//  The version is MAJOR.MINOR.PATCH.
//
//  The SolverName of an input file names only MAJOR.MINOR, so raising PATCH
//  does not require any input file to be rewritten. Any change that alters
//  how the input is interpreted raises MINOR or above.
//
//  1.0 is the first released version.
#define YUZEN_VERSION_Major	1
#define YUZEN_VERSION_Minor	0
#define YUZEN_VERSION_Rev	0


//////////////////////////////////////////
#define _GRID_VERSION_	"YuzenGrid-1.0"





#endif //_DEF_VERSION_INCLUDED_
