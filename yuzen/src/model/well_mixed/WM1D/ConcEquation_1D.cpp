/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ConcEquation_1D.cpp : the 1-D instance of the well-mixed concentration
*                       equation (ConcEquationBase.h)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "ConcEquation_1D.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

template class WellMixed_ConcEquationBase<WM1D_Grid>;

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif
