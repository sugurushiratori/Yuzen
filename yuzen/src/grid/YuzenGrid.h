/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  YuzenGrid.h : the structured-grid and field library of Yuzen. It
*  depends on no external library.
*
*  - Stagger.h  encodes the kind of position in three bits
*  - Site.h     a typed linear index
*  - RectGrid.h uniform rectangular grid (1D / 2D / 3D)
*  - Field.h    a quantity defined on the grid
*
*  The essentials of how to use it:
*
*      RectGrid<2> G;
*      G.Build(0.0, 1.0, nx, 0.0, 1.0, ny);
*
*      ScalarField2D<> h;  h.Build(G);
*
*      for(Zface s = G.First<Zface>(); s != G.Last<Zface>(); ++s)
*          h(s) = ...;
*
*      Zface e = G.ShiftXp(s);      // one cell in +x (same kind of position)
*      Yedge f = G.HalfXp(s);       // half a cell in +x (the kind changes)
*      double xc = G.X(s);
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _YUZEN_GRID_INCLUDED_
#define _YUZEN_GRID_INCLUDED_

#include "Stagger.h"
#include "Site.h"
#include "RectGrid.h"
#include "Field.h"

#endif //_YUZEN_GRID_INCLUDED_
