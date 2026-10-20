/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_grid.cpp : checks the properties of the grid library (src/grid)
*                  directly
*
*  - the linear index and the grid coordinates convert back and forth
*  - neighbors close periodically
*  - a half-cell move followed by its reverse returns to the start
*    (both the site type and the index)
*  - the coordinates follow the staggered layout
*  - a field is zero-initialized, and a sweep visits every point exactly
*    once
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "doctest/doctest.h"

#include "grid/YuzenGrid.h"

using namespace yuzen::grid;

namespace {
const int NX = 9, NY = 7, NZ = 5;
const double X0 = -0.75, X1 = 2.25;
const double Y0 =  1.50, Y1 = 4.50;
const double Z0 = -1.00, Z1 = 1.50;
}


TEST_CASE("Grid: the linear index and the grid coordinates are inverses")
{
	RectGrid<3> G;  G.Build(X0, X1, NX, Y0, Y1, NY, Z0, Z1, NZ);

	CHECK(G.imax() == NX * NY * NZ);

	for(int iz = 0; iz < NZ; ++iz)
	for(int iy = 0; iy < NY; ++iy)
	for(int ix = 0; ix < NX; ++ix)
	{
		Cell c = G.At<Cell>(ix, iy, iz);
		REQUIRE(G.ix(c) == ix);
		REQUIRE(G.iy(c) == iy);
		REQUIRE(G.iz(c) == iz);
	}

	//  It visits every point exactly once
	int n = 0;
	for(Cell c = G.First<Cell>(); c != G.Last<Cell>(); ++c) ++n;
	CHECK(n == G.imax());
}


TEST_CASE("Grid: neighbors close periodically")
{
	RectGrid<2> G;  G.Build(X0, X1, NX, Y0, Y1, NY);

	for(int iy = 0; iy < NY; ++iy)
	{
		Zface lo = G.At<Zface>(0,      iy);
		Zface hi = G.At<Zface>(NX - 1, iy);
		CHECK(G.ShiftXm(lo) == hi);
		CHECK(G.ShiftXp(hi) == lo);
	}
	for(int ix = 0; ix < NX; ++ix)
	{
		Zface lo = G.At<Zface>(ix, 0);
		Zface hi = G.At<Zface>(ix, NY - 1);
		CHECK(G.ShiftYm(lo) == hi);
		CHECK(G.ShiftYp(hi) == lo);
	}

	//  Out-of-range grid coordinates wrap around too
	CHECK(G.At<Zface>(NX, 3)     == G.At<Zface>(0, 3));
	CHECK(G.At<Zface>(-1, 3)     == G.At<Zface>(NX - 1, 3));
	CHECK(G.At<Zface>(2, NY)     == G.At<Zface>(2, 0));
	CHECK(G.At<Zface>(2, -1)     == G.At<Zface>(2, NY - 1));
	CHECK(G.At<Zface>(NX + 2, 3) == G.At<Zface>(2, 3));
}


TEST_CASE("Grid: a one-cell move returns to the start when reversed")
{
	RectGrid<3> G;  G.Build(X0, X1, NX, Y0, Y1, NY, Z0, Z1, NZ);

	for(Cell c = G.First<Cell>(); c != G.Last<Cell>(); ++c)
	{
		REQUIRE(G.ShiftXp(G.ShiftXm(c)) == c);
		REQUIRE(G.ShiftYp(G.ShiftYm(c)) == c);
		REQUIRE(G.ShiftZp(G.ShiftZm(c)) == c);
	}
}


TEST_CASE("Grid: a half-cell move returns to the start, site type included")
{
	RectGrid<2> G;  G.Build(X0, X1, NX, Y0, Y1, NY);

	for(Zface f = G.First<Zface>(); f != G.Last<Zface>(); ++f)
	{
		//  The destination type is determined at compile time
		Yedge wx = G.HalfXm(f);		// from the cell center to the x-normal face
		Xedge wy = G.HalfYm(f);		// from the cell center to the y-normal face
		REQUIRE(G.HalfXp(wx) == f);
		REQUIRE(G.HalfYp(wy) == f);
	}

	for(Node nd = G.First<Node>(); nd != G.Last<Node>(); ++nd)
	{
		Xedge ex = G.HalfXp(nd);
		Yedge ey = G.HalfYp(nd);
		REQUIRE(G.HalfXm(ex) == nd);
		REQUIRE(G.HalfYm(ey) == nd);
	}
}


TEST_CASE("Grid: which way the index moves under a half-cell move")
{
	RectGrid<2> G;  G.Build(X0, X1, NX, Y0, Y1, NY);

	Zface f = G.At<Zface>(4, 3);

	//  A cell center is offset by half a cell in x (bit 1).
	//  Moving half a cell toward -x lands on a cell boundary, so the
	//  index does not move.
	CHECK(G.HalfXm(f).Index() == f.Index());
	//  Moving half a cell toward +x lands on the boundary of the next
	//  cell, so the index advances by one.
	CHECK(G.HalfXp(f).Index() == G.ShiftXp(f).Index());

	Node nd = G.At<Node>(4, 3);
	//  A node is not offset (bit 0). The directions are reversed.
	CHECK(G.HalfXp(nd).Index() == nd.Index());
	CHECK(G.HalfXm(nd).Index() == G.ShiftXm(nd).Index());
}


TEST_CASE("Grid: the coordinates follow the staggered layout")
{
	RectGrid<3> G;  G.Build(X0, X1, NX, Y0, Y1, NY, Z0, Z1, NZ);

	const double dx = (X1 - X0) / NX;
	const double dy = (Y1 - Y0) / NY;
	const double dz = (Z1 - Z0) / NZ;

	CHECK(G.dx() == dx);
	CHECK(G.dy() == dy);
	CHECK(G.dz() == dz);

	for(int iz = 0; iz < NZ; ++iz)
	for(int iy = 0; iy < NY; ++iy)
	for(int ix = 0; ix < NX; ++ix)
	{
		Node  nd = G.At<Node >(ix, iy, iz);
		Cell  cl = G.At<Cell >(ix, iy, iz);
		Xedge xe = G.At<Xedge>(ix, iy, iz);

		//  A node has no offset
		REQUIRE(G.X(nd) == X0 + ix * dx);
		REQUIRE(G.Y(nd) == Y0 + iy * dy);
		REQUIRE(G.Z(nd) == Z0 + iz * dz);

		//  A cell center is offset by half a cell on all three axes
		REQUIRE(G.X(cl) == X0 + ix * dx + dx * 0.5);
		REQUIRE(G.Y(cl) == Y0 + iy * dy + dy * 0.5);
		REQUIRE(G.Z(cl) == Z0 + iz * dz + dz * 0.5);

		//  An x-direction edge is offset by half a cell in x only
		REQUIRE(G.X(xe) == X0 + ix * dx + dx * 0.5);
		REQUIRE(G.Y(xe) == Y0 + iy * dy);
		REQUIRE(G.Z(xe) == Z0 + iz * dz);

		//  Looking it up directly from the grid coordinates gives the same
		REQUIRE(G.XOf<Cell>(ix, iy, iz) == G.X(cl));
		REQUIRE(G.YOf<Cell>(ix, iy, iz) == G.Y(cl));
		REQUIRE(G.ZOf<Cell>(ix, iy, iz) == G.Z(cl));
	}
}


TEST_CASE("Grid: nothing degenerates in a direction holding only one line")
{
	RectGrid<2> G;  G.Build(X0, X1, NX, 0.0, 1.0, 1);

	CHECK(G.ny() == 1);
	for(Zface f = G.First<Zface>(); f != G.Last<Zface>(); ++f)
	{
		REQUIRE(G.iy(f) == 0);
		REQUIRE(G.ShiftYm(f) == f);		// only one line, so it comes back to itself
		REQUIRE(G.ShiftYp(f) == f);
	}
}


TEST_CASE("Field: zero-initialized and readable back by index")
{
	RectGrid<2> G;  G.Build(X0, X1, NX, Y0, Y1, NY);

	ScalarField2D<> h;  h.Build(G);
	CHECK(h.Constructed());
	CHECK(&h.Grid() == &G);

	for(Zface f = G.First<Zface>(); f != G.Last<Zface>(); ++f) REQUIRE(h(f) == 0.0);

	for(Zface f = G.First<Zface>(); f != G.Last<Zface>(); ++f) h(f) = G.X(f);
	for(Zface f = G.First<Zface>(); f != G.Last<Zface>(); ++f) REQUIRE(h(f) == G.X(f));

	//  The same values are visible through the linear index
	for(int i = 1; i <= G.imax(); ++i) REQUIRE(h[i] == h(G.Of<Zface>(i)));

	//  It can be built with an initial value given
	ScalarField2D<> s;  s.Build(G, 3.5);
	for(Zface f = G.First<Zface>(); f != G.Last<Zface>(); ++f) REQUIRE(s(f) == 3.5);
}


TEST_CASE("Field: the components of a vector sit in different places")
{
	RectGrid<2> G;  G.Build(X0, X1, NX, Y0, Y1, NY);

	VectorField2D<> u;  u.Build(G);

	//  The x component sits on the x-normal face (= a y-direction edge),
	//  the y component on the y-normal face (= an x-direction edge)
	Zface f = G.At<Zface>(3, 2);
	u.x(G.HalfXm(f)) = 1.0;
	u.y(G.HalfYm(f)) = 2.0;

	CHECK(u.x(G.At<Yedge>(3, 2)) == 1.0);
	CHECK(u.y(G.At<Xedge>(3, 2)) == 2.0);
}
