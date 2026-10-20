/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ConcEquation_2D.cpp : the 2-D instance of the well-mixed concentration
*                       equation (ConcEquationBase.h)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "ConcEquation_2D.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

template class WellMixed_ConcEquationBase<WM2D_Grid>;

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif
