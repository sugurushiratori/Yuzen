/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CZAsciiLayout.h : the points of a ConcZ field*.dat, shared by the writer
*                    (OutputNeutral::Write(Field_CZ&)) and the ASCII restart
*                    reader (C_Restart::ReadAscii_CZ)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#ifndef _YUZEN_CZ_ASCII_LAYOUT_H_
#define _YUZEN_CZ_ASCII_LAYOUT_H_

#include <vector>

#ifdef _DEF_YUZEN_
namespace yuzen {
#endif

namespace cz_ascii
{
	//  One output point along one axis. At a cell center a == b is that cell;
	//  on a wall face a and b are the ghost and the real cell either side of
	//  it (in increasing index), and the value written there is their mean.
	struct Pt { int a, b; double x; };

	//  The points along one axis whose real cells are the half-open range
	//  [lo, hi): the face at lo, every cell center, the face at hi. center(i)
	//  gives the coordinate of the center of cell i; xlo and xhi those of the
	//  two faces.
	template<class F>
	std::vector<Pt> Points(int lo, int hi, F center, double xlo, double xhi)
	{
		std::vector<Pt> p;
		p.push_back(Pt{lo - 1, lo, xlo});
		for(int i = lo; i < hi; i++) p.push_back(Pt{i, i, center(i)});
		p.push_back(Pt{hi - 1, hi, xhi});
		return p;
	}

	inline bool IsCenter(const Pt &p) { return p.a == p.b; }
}

#ifdef _DEF_YUZEN_
}
#endif

#endif
