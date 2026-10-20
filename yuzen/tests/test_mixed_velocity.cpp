/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_mixed_velocity.cpp : pins down how the mixed velocity is assembled
*
*  The governing equation is h_t = -div(Q). When it is expanded inside the
*  -div, the gravity term comes out with a + sign and all the others with a
*  - sign. The _Magnitude of each physical term follows that convention, so
*  the velocity field each term writes already carries its physical
*  direction. The mixed velocity must therefore be the **plain sum** of the
*  individual terms.
*
*  Written as a chain of Math::Diff(-ux, b) calls, since Diff(a,b) = a-b
*  and hence Diff(-u,b) = -u-b, it would compute
*
*      _mixed = LP + SC - VG + HG - CF
*
*  so vertical gravity and the centrifugal force would enter with the wrong
*  sign. The mixed velocity is used for the advection in the concentration
*  equation, which means the solute would be transported in the opposite
*  direction. On top of that, adding one more term would flip the sign of
*  every term after it, so no new physical term could be added.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "doctest/doctest.h"

#include "LubFields.h"

using namespace yuzen;
using namespace yuzen::grid;

namespace {

Input::I_Grid MakeGrid1D(int nx)
{
	Input::I_Grid g;
	g.Dimension = "xh";
	g.Xmin = 0.0; g.Xmax = 1.0;
	g.Ymin = 0.0; g.Ymax = 1.0;
	g.Nx = nx; g.Ny = 1; g.Nz = 1;
	return g;
}

Input::I_Grid MakeGrid2D(int nx, int ny)
{
	Input::I_Grid g = MakeGrid1D(nx);
	g.Dimension = "xyh";
	g.Ny = ny;
	return g;
}

//  Put a different constant into each component
const double V_LP = 1.0;
const double V_SC = 2.0;
const double V_VG = 4.0;
const double V_HG = 8.0;
const double V_CF = 16.0;
const double V_DP = 32.0;
const double V_SUM = V_LP + V_SC + V_VG + V_HG + V_CF + V_DP;

} // namespace


TEST_CASE("Mixed velocity: in 1-D it is the plain sum of the terms")
{
	C_Grid1D G;
	G.Build(MakeGrid1D(8));

	Velocity_WM1D V;
	V.Build(G);

	for(Node nd = G.First<Node>(); nd != G.Last<Node>(); ++nd)
	{
		V.U_LaplacePressure().x(nd)   = V_LP;
		V.U_SolutoCapillarity().x(nd) = V_SC;
		V.U_VerticalGravity().x(nd)   = V_VG;
		V.U_HorizontalGravity().x(nd) = V_HG;
		V.U_Centrifugal().x(nd)       = V_CF;
		V.U_DisjoiningPressure().x(nd)= V_DP;
	}

	V.Cal_MixedVelocity();

	for(Node nd = G.First<Node>(); nd != G.Last<Node>(); ++nd)
		REQUIRE(V.UU().x(nd) == doctest::Approx(V_SUM));
}


TEST_CASE("Mixed velocity: in 1-D signs pass through unchanged")
{
	//  Only gravity is active. The direction written by the physical term
	//  must come out unchanged (a chain of differences flips it here).
	C_Grid1D G;
	G.Build(MakeGrid1D(8));

	Velocity_WM1D V;
	V.Build(G);

	for(Node nd = G.First<Node>(); nd != G.Last<Node>(); ++nd)
		V.U_VerticalGravity().x(nd) = 3.0;

	V.Cal_MixedVelocity();

	for(Node nd = G.First<Node>(); nd != G.Last<Node>(); ++nd)
		REQUIRE(V.UU().x(nd) == doctest::Approx(3.0));

	//  The same holds for the centrifugal force
	Velocity_WM1D W;
	W.Build(G);
	for(Node nd = G.First<Node>(); nd != G.Last<Node>(); ++nd)
		W.U_Centrifugal().x(nd) = -5.0;
	W.Cal_MixedVelocity();
	for(Node nd = G.First<Node>(); nd != G.Last<Node>(); ++nd)
		REQUIRE(W.UU().x(nd) == doctest::Approx(-5.0));
}


TEST_CASE("Mixed velocity: in 2-D it is the plain sum too (x and y components alike)")
{
	C_Grid2D G;
	G.Build(MakeGrid2D(6, 5));

	Velocity_WM2D V;
	V.Build(G);

	for(Yedge ye = G.First<Yedge>(); ye != G.Last<Yedge>(); ++ye)
	{
		V.U_LaplacePressure().x(ye)   = V_LP;
		V.U_SolutoCapillarity().x(ye) = V_SC;
		V.U_VerticalGravity().x(ye)   = V_VG;
		V.U_HorizontalGravity().x(ye) = V_HG;
		V.U_Centrifugal().x(ye)       = V_CF;
		V.U_DisjoiningPressure().x(ye)= V_DP;
	}
	for(Xedge xe = G.First<Xedge>(); xe != G.Last<Xedge>(); ++xe)
	{
		V.U_LaplacePressure().y(xe)   = -V_LP;
		V.U_SolutoCapillarity().y(xe) = -V_SC;
		V.U_VerticalGravity().y(xe)   = -V_VG;
		V.U_HorizontalGravity().y(xe) = -V_HG;
		V.U_Centrifugal().y(xe)       = -V_CF;
		V.U_DisjoiningPressure().y(xe)= -V_DP;
	}

	V.Cal_MixedVelocity();

	for(Yedge ye = G.First<Yedge>(); ye != G.Last<Yedge>(); ++ye)
		REQUIRE(V.UU().x(ye) == doctest::Approx(V_SUM));
	for(Xedge xe = G.First<Xedge>(); xe != G.Last<Xedge>(); ++xe)
		REQUIRE(V.UU().y(xe) == doctest::Approx(-V_SUM));
}
