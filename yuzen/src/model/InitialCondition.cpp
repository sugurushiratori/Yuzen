/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  InitialCondition.cpp : implementation of the initial film thickness distribution
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "InitialCondition.h"
#include "MyMath.h"
#include "GenerateShape.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

//------------------------------------------------------------------
//  The initial film must be a film: a finite thickness and substrate at
//  every point, and liquid somewhere. A shape with a zero width or
//  wavelength (or a rim softened over zero, met exactly by a grid point)
//  divides by zero, and a film of zero thickness everywhere -- no precursor
//  and no shape above the substrate -- has no liquid in it at all (it would
//  run to the end with NaN in History.dat). Both stop the run here, before
//  any step.
//------------------------------------------------------------------
namespace {
struct FilmCheck
{
	bool finite = true, liquid = false;
	double x = 0, y = 0;
	void Point(double h, double surfaces, double px, double py)
	{
		if(!(std::isfinite(h) && std::isfinite(surfaces)))
		{
			if(finite) { x = px; y = py; }
			finite = false;
		}
		else if(h > 0) liquid = true;
	}
	void Report() const
	{
		if(!finite)
		{
			Global::FatalError<<"Initial.SurfaceShape / Substrate: the initial film thickness or the"
			                  <<" substrate height is not a finite number at (x, y) = ("<<x<<", "<<y
			                  <<") m. Check the Width, Wavelength and SlopeWidth of the shapes, and"
			                  <<" InnerLength, OuterLength and Alpha of an EdgeBump."<<nl;
			Global::FatalError.Abort();
		}
		if(!liquid)
		{
			Global::FatalError<<"Initial.SurfaceShape: the initial film has no liquid in it (its"
			                  <<" thickness is zero everywhere). Give a PrecursorRatio above 0, or a"
			                  <<" shape that rises above the substrate."<<nl;
			Global::FatalError.Abort();
		}
	}
};
}  // namespace

void C_InitialSurface::Set(const C_Grid1D &G, ScalarField1D<> &HH, ScalarField1D<> &SS, Input::I_Initial::I_SurfaceShape &inpH, Input::I_Substrate &inpS, double scale)
{
	double xx;
	double h_ast = inpH.PrecursorRatio;

	double mixed_shape;
	double hs, ss;


	C_GenerateShape geneShape;
	geneShape.Build(G, scale, inpH);
	FilmCheck check;

	for(Xedge xe = G.First<Xedge>(); xe != G.Last<Xedge>(); ++xe)
	{
		xx = G.X(xe) * scale;

		ss = SS(xe);

		mixed_shape = geneShape.Evaluate(xx);
		hs = inpH.TotalHeight * (mixed_shape * (1.0-h_ast) + h_ast);
		HH(xe) = Math::Max(hs - ss, h_ast*inpH.TotalHeight);
		if(G.CheckGhost(xe) == C_Grid1D::Inside) check.Point(HH(xe), hs + ss, xx, 0);
	}
	check.Report();
}
void C_InitialSurface::Set(const C_Grid2D &G, ScalarField2D<> &HH, ScalarField2D<> &SS, Input::I_Initial::I_SurfaceShape &inpH, Input::I_Substrate &inpS, double scale)
{
	double xx, yy;
	double h_ast = inpH.PrecursorRatio;

	double mixed_shape;
	double hs, ss;


	C_GenerateShape geneShape;
	geneShape.Build(G, scale, inpH);
	FilmCheck check;

	for(Zface zf = G.First<Zface>(); zf != G.Last<Zface>(); ++zf)
	{
		xx = G.X(zf) * scale;
		yy = G.Y(zf) * scale;

		ss = SS(zf);

		mixed_shape = geneShape.Evaluate(xx,yy);
		hs = inpH.TotalHeight * (mixed_shape * (1.0-h_ast) + h_ast);
		HH(zf) = Math::Max(hs - ss, h_ast*inpH.TotalHeight);
		if(G.CheckGhost(zf) == C_Grid2D::Inside) check.Point(HH(zf), hs + ss, xx, yy);
	}
	check.Report();
}


void C_InitialSurface::Set(Field_WM1D& field, Input::I_Initial::I_SurfaceShape &inpH, Input::I_Substrate &inpS, double scale)
{
	Set(field.GetGrid(), field.HH(), field.SS(), inpH,inpS, scale);
}

void C_InitialSurface::Set(Field_WM2D& field, Input::I_Initial::I_SurfaceShape &inpH, Input::I_Substrate &inpS, double scale)
{
	Set(field.GetGrid(), field.HH(), field.SS(), inpH,inpS, scale);
}

void C_InitialSurface::Set(Field_CZ& field, Input::I_Initial::I_SurfaceShape &inpH, Input::I_Substrate &inpS, double scale)
{
	Set(field.GetGrid2(), field.HH(), field.SS(), inpH,inpS, scale);
}
