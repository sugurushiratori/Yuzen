/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Substrate.cpp : generation of the substrate shape
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <fstream>
#include <cmath>

#include "MyMath.h"
#include "Substrate.h"
#include "GenerateShape.h"


#include "Output.h"


#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif



void C_Substrate::MakeSubstrateProfile(const C_Grid1D &G, ScalarField1D<> &profile, const Input::I_Substrate &inp, double scale)
{
	double xx;
	double mixed_shape;


	C_GenerateShape geneShape;
	geneShape.Build(G, scale, inp);

	for(Xedge xe = G.First<Xedge>(); xe != G.Last<Xedge>(); ++xe)
	{
		xx = G.X(xe) * scale;

		mixed_shape = geneShape.Evaluate(xx);
		profile(xe) = inp.TotalHeight * mixed_shape;
	}
}

void C_Substrate::MakeSubstrateProfile(const C_Grid2D &G, ScalarField2D<> &profile, const Input::I_Substrate &inp, double scale)
{
	double xx, yy;
	double mixed_shape;


	C_GenerateShape geneShape;
	geneShape.Build(G, scale, inp);

	for(Zface zf = G.First<Zface>(); zf != G.Last<Zface>(); ++zf)
	{
		xx = G.X(zf) * scale;
		yy = G.Y(zf) * scale;

		mixed_shape = geneShape.Evaluate(xx,yy);
		profile(zf) = inp.TotalHeight * mixed_shape;
	}
}
