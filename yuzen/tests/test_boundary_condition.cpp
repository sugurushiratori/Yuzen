/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_boundary_condition.cpp : test_boundary_condition.cpp  --  regression tests for model/BoundaryCondition
*
*  The four edges (two in 1-D) share one body of code through EdgeSlot,
*  which differs per direction only in its row numbers, the start of its
*  stencil and its column stride. A starting column off by one in one
*  direction would refer outside the matrix (column _Nx), and MatVec would
*  read x[_Nx], which is undefined behavior -- a defect hard to find by eye.
*
*  For every combination of boundary conditions, we check here that
*    - the column indices stay within the range of the matrix
*    - every ghost row carries exactly one equation (a non-zero diagonal)
*    - if the boundary conditions are symmetric, the matrix is symmetric
*      left-right and top-bottom as well
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "doctest/doctest.h"
#include "BoundaryCondition.h"
#include "LubFields.h"
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace yuzen;

namespace {

typedef BoundaryCondition::tBoundaryCondition TBC;

const TBC kTypes[5] = { BoundaryCondition::Symmetry, BoundaryCondition::Periodic,
                        BoundaryCondition::Pinning,  BoundaryCondition::FreshInlet,
                        BoundaryCondition::Outlet };
const char *kName[5] = { "Symmetry", "Periodic", "Pinning", "FreshInlet", "Outlet" };

const BoundaryCondition::tField kFields[3] = {
	BoundaryCondition::Height, BoundaryCondition::HeightJacobi, BoundaryCondition::Conc };

//  Combinations whose corner handling is unimplemented and which abort
bool NotImplemented(int x, int y)
{
	if(x == 2 && (y == 3 || y == 4)) return true;	// Pinning - (Inlet|Outlet)
	if((x == 3 || x == 4) && y == 2) return true;	// (Inlet|Outlet) - Pinning
	return false;
}

Input::I_Grid MakeGrid(const char *dim, int nx, int ny)
{
	Input::I_Grid g;
	g.Dimension = dim;
	g.Xmin = -0.5; g.Xmax = 0.5;
	g.Ymin = -0.5; g.Ymax = 0.5;
	g.Nx = nx; g.Ny = ny; g.Nz = 1;
	return g;
}

//  Set up a state with empty ghost rows (the BC flag only) and interior
//  rows forming the identity matrix
void PrepareMatrix(CRS_Matrix &crs, Vector<double> &rhs, int nx, int ny)
{
	const int n = nx * ny;
	crs.Build(n);
	rhs.resize(n);
	for(int i = 0; i < n; ++i)
	{
		const int ix = i % nx;
		const int jy = (ny > 1) ? i / nx : 0;
		bool ghost = (ix < 2 || ix >= nx - 2);
		if(ny > 1 && (jy < 2 || jy >= ny - 2)) ghost = true;

		if(ghost) crs.SetAsBC(i);
		else      crs.Diag(i) = 1.0;
		rhs[i] = 0.0;
	}
}

//  Are the column indices of every row within [0, N)?
bool ColumnsInRange(const CRS_Matrix &m)
{
	for(int i = 0; i < m.RowSize(); ++i)
		for(int k = 0; k < m.NonZero(i); ++k)
		{
			const int j = m.Jcol(i, k);
			if(j < 0 || j >= m.RowSize()) return false;
		}
	return true;
}

//  Does row i carry an equation? (a non-zero diagonal, or at least two
//  non-zero elements)
bool RowIsConstrained(const CRS_Matrix &m, int i)
{
	if(m.Diag(i) != 0.0) return true;
	int nz = 0;
	for(int k = 0; k < m.NonZero(i); ++k) if(m.Element(i, k) != 0.0) ++nz;
	return nz > 0;
}

} // namespace


TEST_CASE("BoundaryCondition1D: every combination keeps the column indices in range")
{
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;

	Field_WM1D f;
	f.Build(MakeGrid("1D", 12, 1), subs, 1.0);

	C_Grid1D *G = f.GetGridPtr();
	int n = 0;
	for(Xedge xe = G->First<Xedge>(); xe != G->Last<Xedge>(); ++xe, ++n)
	{
		f.HH(xe) = 1.0 + 0.01 * (n % 7);
		f.CC(xe) = 0.5 + 0.01 * (n % 5);
	}
	const int NX = G->nx();

	for(int xa = 0; xa < 5; ++xa)
	for(int xb = 0; xb < 5; ++xb)
	{
		// Build aborts unless the periodic condition is set on both sides
		if((xa == 1) != (xb == 1)) continue;

		for(int fi = 0; fi < 3; ++fi)
		{
			CAPTURE(kName[xa]);
			CAPTURE(kName[xb]);
			CAPTURE(fi);

			CRS_Matrix crs; Vector<double> rhs;
			PrepareMatrix(crs, rhs, NX, 1);

			BoundaryCondition1D bc;
			bc.Build(NX, kTypes[xa], kTypes[xb]);
			bc.Apply(crs, rhs, kFields[fi], f);

			// No column outside the matrix (the starting column of a stencil)
			CHECK(ColumnsInRange(crs));

			// All four ghost rows must carry an equation
			CHECK(RowIsConstrained(crs, 0));
			CHECK(RowIsConstrained(crs, 1));
			CHECK(RowIsConstrained(crs, NX-2));
			CHECK(RowIsConstrained(crs, NX-1));
		}
	}
}


TEST_CASE("BoundaryCondition2D: every combination keeps the column indices in range")
{
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;

	Field_WM2D f;
	f.Build(MakeGrid("2D", 10, 8), subs, 1.0);

	C_Grid2D *G = f.GetGridPtr();
	int n = 0;
	for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf, ++n)
	{
		f.HH(zf) = 1.0 + 0.01 * (n % 7);
		f.CC(zf) = 0.5 + 0.01 * (n % 5);
	}
	const int NX = G->nx();
	const int NY = G->ny();

	int checked = 0;
	for(int xa = 0; xa < 5; ++xa)
	for(int xb = 0; xb < 5; ++xb)
	for(int ya = 0; ya < 5; ++ya)
	for(int yb = 0; yb < 5; ++yb)
	{
		if((xa == 1) != (xb == 1)) continue;
		if((ya == 1) != (yb == 1)) continue;
		if(NotImplemented(xa,ya) || NotImplemented(xa,yb) ||
		   NotImplemented(xb,ya) || NotImplemented(xb,yb)) continue;

		for(int fi = 0; fi < 3; ++fi)
		{
			CRS_Matrix crs; Vector<double> rhs;
			PrepareMatrix(crs, rhs, NX, NY);

			BoundaryCondition2D bc;
			bc.Build(NX, NY, kTypes[xa], kTypes[xb], kTypes[ya], kTypes[yb]);
			bc.Apply(crs, rhs, kFields[fi], f);

			if(!ColumnsInRange(crs))
			{
				CAPTURE(kName[xa]); CAPTURE(kName[xb]);
				CAPTURE(kName[ya]); CAPTURE(kName[yb]); CAPTURE(fi);
				CHECK(ColumnsInRange(crs));
			}
			++checked;
		}
	}
	CHECK(checked > 400);	// no combination has been missed
}


TEST_CASE("BoundaryCondition2D: the sized constructor imposes what Build() does")
{
	//  The constructor taking the sizes and the four kinds is Build() in one
	//  step, the corner index tables included: the two give the same rows.
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;

	Field_WM2D f;
	f.Build(MakeGrid("2D", 10, 8), subs, 1.0);
	C_Grid2D *G = f.GetGridPtr();
	int n = 0;
	for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf, ++n)
	{
		f.HH(zf) = 1.0 + 0.01 * (n % 7);
		f.CC(zf) = 0.5 + 0.01 * (n % 5);
	}
	const int NX = G->nx();
	const int NY = G->ny();

	const int combos[3][4] = { {0, 0, 0, 0}, {4, 3, 1, 1}, {1, 1, 2, 0} };
	for(const auto &c : combos)
	for(int fi = 0; fi < 3; ++fi)
	{
		CAPTURE(kName[c[0]]); CAPTURE(kName[c[1]]);
		CAPTURE(kName[c[2]]); CAPTURE(kName[c[3]]); CAPTURE(fi);

		CRS_Matrix a, b; Vector<double> ra, rb;
		PrepareMatrix(a, ra, NX, NY);
		PrepareMatrix(b, rb, NX, NY);

		BoundaryCondition2D built;
		built.Build(NX, NY, kTypes[c[0]], kTypes[c[1]], kTypes[c[2]], kTypes[c[3]]);
		built.Apply(a, ra, kFields[fi], f);

		BoundaryCondition2D sized(NX, NY, kTypes[c[0]], kTypes[c[1]], kTypes[c[2]], kTypes[c[3]]);
		sized.Apply(b, rb, kFields[fi], f);

		bool same = true;
		for(int i = 0; i < NX * NY && same; ++i)
		{
			if(a.NonZero(i) != b.NonZero(i) || a.Diag(i) != b.Diag(i) || ra[i] != rb[i])
				same = false;
			for(int k = 0; k < a.NonZero(i) && same; ++k)
				if(a.Jcol(i, k) != b.Jcol(i, k) || a.Element(i, k) != b.Element(i, k))
					same = false;
		}
		CHECK(same);
	}
}


TEST_CASE("BoundaryCondition2D: four symmetric edges give a mirror-image matrix")
{
	// When the same condition is given on Xmin/Xmax/Ymin/Ymax, the rows on
	// the Xmin side and those on the Xmax side must coincide under the
	// mirroring ix -> NX-1-ix. This confirms that EdgeSlot, which serves
	// all four edges, leaves no room for a per-direction transcription
	// mistake.
	//
	// Note, however, that the coefficients of Derive3 / Neumann are
	// anti-palindromic (-1,+3,-3,+1), so on the max side, where inward and
	// outward are reversed, the sign of the whole row is flipped.
	// Multiplying a whole row by -1 gives the same equation, so we allow
	// this freedom of sign.
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;

	Field_WM2D f;
	f.Build(MakeGrid("2D", 10, 8), subs, 1.0);

	C_Grid2D *G = f.GetGridPtr();
	for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf)
	{
		f.HH(zf) = 1.0;		// with a uniform field the boundary values match on both sides
		f.CC(zf) = 0.5;
	}
	const int NX = G->nx();
	const int NY = G->ny();

	for(int t = 0; t < 5; ++t)
	{
		if(t == 2) continue;	// Pinning handles the corners asymmetrically
		for(int fi = 0; fi < 3; ++fi)
		{
			CAPTURE(std::string(kName[t]));
			CAPTURE(fi);

			CRS_Matrix crs; Vector<double> rhs;
			PrepareMatrix(crs, rhs, NX, NY);

			BoundaryCondition2D bc;
			bc.Build(NX, NY, kTypes[t], kTypes[t], kTypes[t], kTypes[t]);
			bc.Apply(crs, rhs, kFields[fi], f);

			bool ok = true;
			for(int jy = 2; jy < NY-2 && ok; ++jy)
			for(int k  = 0; k  < 2     && ok; ++k)
			{
				const int rm = NX*jy + k;
				const int rp = NX*jy + (NX-1-k);
				if(crs.NonZero(rm) != crs.NonZero(rp)) { ok = false; break; }

				// Does the mirrored equation of the Xmin-side row agree
				// with the Xmax-side row (up to sign)?
				// The storage order is reversed on the max side, so we
				// look the entries up again by column index and compare.
				double sign = 0.0;
				for(int s = 0; s < crs.NonZero(rm); ++s)
				{
					const int jm = crs.Jcol(rm, s);
					if(jm / NX != jy) { ok = false; break; }		// it should close within the same row
					const int jmir = NX*jy + (NX-1 - (jm % NX));

					const int t2 = crs.ColInd(rp, jmir);
					if(t2 < 0) { ok = false; break; }

					const double a = crs.Element(rm, s);
					const double b = crs.Element(rp, t2);
					if(a == 0.0 && b == 0.0) continue;
					if(sign == 0.0)      sign = (a == b) ? 1.0 : ((a == -b) ? -1.0 : 0.0);
					if(sign == 0.0)      { ok = false; break; }
					if(b != sign * a)    { ok = false; break; }
				}
				if(!ok) break;
				if(sign == 0.0) sign = 1.0;
				if(rhs[rp] != sign * rhs[rm]) { ok = false; break; }
			}
			CHECK(ok);
		}
	}
}


//==========================================================================
//  The values the conditions give the ghost points.
//
//  The checks above look only at the shape of the matrix. What an inlet
//  supplies -- the height and concentration held at its face, and the
//  values at the corners it makes with the other edges -- is in the
//  right-hand side and in the coefficients of the corner stencils, and a
//  run that starts uniform cannot see most of it: in the well-mixed model
//  nothing reads the concentration at a corner ghost, nor the height at
//  the outer corner ghosts. So the values are checked here, on the rows
//  themselves:
//    - a uniform film satisfies every condition, so imposing the conditions
//      on it must give it back at every ghost point, corners included;
//    - an inlet holds the mean of its two face points at the value of the
//      field it was given (the initial field, for the concentration);
//    - mirroring or transposing the field and the edge types together
//      mirrors or transposes every ghost value.
//==========================================================================
namespace {

//  Set the height and the concentration at every point, ghosts included.
void Fill(Field_WM2D &f, double (*h)(int, int), double (*c)(int, int))
{
	C_Grid2D *G = f.GetGridPtr();
	for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf)
	{
		f.HH(zf) = h(G->ix(zf), G->iy(zf));
		f.CC(zf) = c(G->ix(zf), G->iy(zf));
	}
}

//  Solve the rows Apply() writes, with the interior held at `now` and the
//  boundary values taken from `ini`. Every ghost row reads
//  a_rr x_r + sum a_rj x_j = rhs_r and involves only the interior and a few
//  other ghosts, so Gauss-Seidel sweeps settle it to rounding (the two
//  rows of a height inlet or outlet are coupled, hence the many sweeps;
//  on a grid this small they cost nothing).
std::vector<double> Impose(BoundaryCondition2D &bc, BoundaryCondition::tField tf,
                           Field_WM2D &ini, Field_WM2D &now)
{
	C_Grid2D *G = ini.GetGridPtr();
	const int NX = G->nx(), NY = G->ny();
	CRS_Matrix crs; Vector<double> rhs;
	PrepareMatrix(crs, rhs, NX, NY);
	std::vector<double> x(NX*NY, 0.0);
	for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf)
		if(G->CheckGhost(zf) == C_Grid2D::Inside)
			x[G->IL(zf)] = (tf == BoundaryCondition::Conc) ? now.CC(zf) : now.HH(zf);
	bc.Apply(crs, rhs, tf, ini);
	for(int sweep = 0; sweep < 200; ++sweep)
		for(int r = 0; r < NX*NY; ++r)
		{
			if(!crs.IsBC(r)) continue;
			double s = rhs[r], d = 0.0;
			for(int k = 0; k < crs.NonZero(r); ++k)
			{
				const int j = crs.Jcol(r, k);
				if(j == r) d += crs.Element(r, k);
				else       s -= crs.Element(r, k) * x[j];
			}
			x[r] = s / d;
		}
	return x;
}

double Uh(int, int) { return 0.8; }
double Uc(int, int) { return 0.3; }
//  Smooth, with no symmetry of its own, and set in the ghosts too
double Ph(int i, int j) { return 1.0 + 0.05*i - 0.03*j + 0.004*i*j + 0.02*std::sin(0.9*i + 0.4*j); }
double Pc(int i, int j) { return 0.4 - 0.02*i + 0.035*j - 0.003*i*j + 0.015*std::cos(0.7*i - 1.1*j); }
//  Another state of the interior: what the film has become since
double Qh(int i, int j) { return 1.3 - 0.02*i + 0.06*j + 0.01*std::sin(1.7*i*j); }
double Qc(int i, int j) { return 0.55 + 0.01*i - 0.02*j + 0.01*std::cos(1.3*i + 0.2*j); }

//  The kinds that can meet an inlet at a corner (Pinning cannot, and aborts)
const int kIO[4] = { 0, 1, 3, 4 };	// Symmetry, Periodic, FreshInlet, Outlet

bool Valid(int xa, int xb, int ya, int yb)
{
	if((xa == 1) != (xb == 1)) return false;
	if((ya == 1) != (yb == 1)) return false;
	return !(NotImplemented(xa,ya) || NotImplemented(xa,yb) ||
	         NotImplemented(xb,ya) || NotImplemented(xb,yb));
}

} // namespace


TEST_CASE("BoundaryCondition1D: an inlet at either end holds its face at the value it was given")
{
	//  The face value (the cubic interpolation across it) is the mean of the
	//  inner ghost and the first interior point of the field given, at
	//  Xmax as at Xmin (the two ghosts would make the same inlet hold
	//  another value at the other end).
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;
	Field_WM1D f;
	f.Build(MakeGrid("1D", 12, 1), subs, 1.0);
	C_Grid1D *G = f.GetGridPtr();
	const int NX = G->nx();
	for(Xedge xe = G->First<Xedge>(); xe != G->Last<Xedge>(); ++xe)
	{
		f.HH(xe) = Ph(G->ix(xe), 3);
		f.CC(xe) = Pc(G->ix(xe), 3);
	}
	for(int end = 0; end < 2; ++end)
	for(int fi = 0; fi < 3; fi += 2)
	{
		CAPTURE(end); CAPTURE(fi);
		const bool conc = (fi == 2);
		CRS_Matrix crs; Vector<double> rhs;
		PrepareMatrix(crs, rhs, NX, 1);
		BoundaryCondition1D bc;
		bc.Build(NX, end == 0 ? BoundaryCondition::FreshInlet : BoundaryCondition::Outlet,
		             end == 0 ? BoundaryCondition::Outlet : BoundaryCondition::FreshInlet);
		bc.Apply(crs, rhs, kFields[fi], f);
		std::vector<double> x(NX, 0.0);
		for(int i = 2; i < NX-2; ++i) x[i] = conc ? Qc(i, 3) : Qh(i, 3);
		for(int sweep = 0; sweep < 200; ++sweep)
			for(int r = 0; r < NX; ++r)
			{
				if(!crs.IsBC(r)) continue;
				double s = rhs[r], d = 0.0;
				for(int k = 0; k < crs.NonZero(r); ++k)
				{
					const int j = crs.Jcol(r, k);
					if(j == r) d += crs.Element(r, k); else s -= crs.Element(r, k) * x[j];
				}
				x[r] = s / d;
			}
		int L[4];
		for(int k = 0; k < 4; ++k) L[k] = (end == 0) ? k : NX-1-k;
		const double want = conc ? 0.5*(Pc(L[1],3) + Pc(L[2],3)) : 0.5*(Ph(L[1],3) + Ph(L[2],3));
		const double face = (-x[L[0]] + 9.0*x[L[1]] + 9.0*x[L[2]] - x[L[3]]) / 16.0;
		CHECK(face == doctest::Approx(want).epsilon(1e-13));
	}
}


TEST_CASE("BoundaryCondition2D: a uniform film is given back at every ghost point, corners included")
{
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;
	Field_WM2D f;
	f.Build(MakeGrid("2D", 10, 8), subs, 1.0);
	Fill(f, Uh, Uc);
	const int NX = f.GetGridPtr()->nx(), NY = f.GetGridPtr()->ny();

	int checked = 0;
	for(int xa = 0; xa < 5; ++xa) for(int xb = 0; xb < 5; ++xb)
	for(int ya = 0; ya < 5; ++ya) for(int yb = 0; yb < 5; ++yb)
	{
		if(!Valid(xa, xb, ya, yb)) continue;
		BoundaryCondition2D bc;
		bc.Build(NX, NY, kTypes[xa], kTypes[xb], kTypes[ya], kTypes[yb]);
		const bool pinned = (xa == 2 || xb == 2 || ya == 2 || yb == 2);
		for(int fi = 0; fi < 3; fi += 2)	// Height and Conc
		{
			//  Pinning holds the height at zero, not at its value
			if(fi == 0 && pinned) continue;
			const double want = (fi == 0) ? Uh(0,0) : Uc(0,0);
			std::vector<double> x = Impose(bc, kFields[fi], f, f);
			double worst = 0.0; int where = 0;
			for(int r = 0; r < NX*NY; ++r)
				if(std::fabs(x[r] - want) > worst) { worst = std::fabs(x[r] - want); where = r; }
			if(worst > 1e-12)
			{
				const std::string combo = std::string(kName[xa]) + " " + kName[xb] + " " + kName[ya] + " " + kName[yb];
				CAPTURE(combo);
				CAPTURE(fi); CAPTURE(where % NX); CAPTURE(where / NX); CAPTURE(x[where]);
				CHECK(worst <= 1e-12);
			}
			++checked;
		}
	}
	CHECK(checked > 200);	// 237: no combination has been missed
}


TEST_CASE("BoundaryCondition2D: an inlet holds its face at the height and concentration it was given")
{
	//  On each edge in turn, next to each kind of edge. The value at the face
	//  (the cubic interpolation of the four points across it, as Dirichlet()
	//  writes it) is held at the mean of the inner ghost and the first
	//  interior point of the field the condition was built from, whatever the
	//  interior has become since. The outer ghost continues the line through
	//  the two for the concentration, and has a zero third difference for the
	//  height.
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;
	Field_WM2D ini, now;
	ini.Build(MakeGrid("2D", 10, 8), subs, 1.0);
	now.Build(MakeGrid("2D", 10, 8), subs, 1.0);
	Fill(ini, Ph, Pc);
	Fill(now, Qh, Qc);
	C_Grid2D *G = ini.GetGridPtr();
	const int NX = G->nx(), NY = G->ny();

	int checked = 0;
	for(int edge = 0; edge < 4; ++edge)		// Xmin, Xmax, Ymin, Ymax
	for(int o = 0; o < 4; ++o)				// the kind of the other edges
	{
		int t[4] = { kIO[o], kIO[o], kIO[o], kIO[o] };
		t[edge] = 3;
		if(kIO[o] == 1) t[edge ^ 1] = 4;	// the pair with the inlet cannot be periodic: an outlet opposite
		if(!Valid(t[0], t[1], t[2], t[3])) continue;
		BoundaryCondition2D bc;
		bc.Build(NX, NY, kTypes[t[0]], kTypes[t[1]], kTypes[t[2]], kTypes[t[3]]);

		for(int fi = 0; fi < 3; fi += 2)
		{
			const bool conc = (fi == 2);
			std::vector<double> x = Impose(bc, kFields[fi], ini, now);
			const int n = (edge < 2) ? NY : NX;
			for(int m = 2; m < n - 2; ++m)
			{
				int L[4], I[4], J[4];
				for(int k = 0; k < 4; ++k)
				{
					switch(edge)
					{
					case 0:  I[k] = k;        J[k] = m;        break;
					case 1:  I[k] = NX-1-k;   J[k] = m;        break;
					case 2:  I[k] = m;        J[k] = k;        break;
					default: I[k] = m;        J[k] = NY-1-k;   break;
					}
					L[k] = NX*J[k] + I[k];
				}
				const double want = conc ? 0.5*(Pc(I[1],J[1]) + Pc(I[2],J[2]))
				                         : 0.5*(Ph(I[1],J[1]) + Ph(I[2],J[2]));
				const std::string other = kName[kIO[o]];
				CAPTURE(edge); CAPTURE(other); CAPTURE(fi); CAPTURE(m);
				const double face = (-x[L[0]] + 9.0*x[L[1]] + 9.0*x[L[2]] - x[L[3]]) / 16.0;
				CHECK(face == doctest::Approx(want).epsilon(1e-13));
				if(conc) CHECK(std::fabs(x[L[0]] - 2.0*x[L[1]] + x[L[2]]) <= 1e-13);
				else     CHECK(std::fabs(-x[L[0]] + 3.0*x[L[1]] - 3.0*x[L[2]] + x[L[3]]) <= 1e-13);
				++checked;
			}
		}
	}
	CHECK(checked > 150);
}


TEST_CASE("BoundaryCondition2D: mirroring or transposing the film and the edges moves every ghost value with them")
{
	//  With a field that has no symmetry of its own, an edge or a corner
	//  handled on one side only -- a boundary value read from the wrong
	//  points, one row of a corner stencil changed -- makes the image of
	//  the result differ from the result for the image.
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;
	const int nx = 10, ny = 8;

	Field_WM2D f, fX, fY, fT;
	f.Build (MakeGrid("2D", nx, ny), subs, 1.0);
	fX.Build(MakeGrid("2D", nx, ny), subs, 1.0);	// x -> -x
	fY.Build(MakeGrid("2D", nx, ny), subs, 1.0);	// y -> -y
	fT.Build(MakeGrid("2D", ny, nx), subs, 1.0);	// x <-> y
	Fill(f, Ph, Pc);
	C_Grid2D *G = f.GetGridPtr();
	const int NX = G->nx(), NY = G->ny();
	for(int j = 0; j < NY; ++j) for(int i = 0; i < NX; ++i)
	{
		Zface a  = G->At<Zface>(i, j);
		Zface bx = fX.GetGridPtr()->At<Zface>(NX-1-i, j);
		Zface by = fY.GetGridPtr()->At<Zface>(i, NY-1-j);
		Zface bt = fT.GetGridPtr()->At<Zface>(j, i);
		fX.HH(bx) = f.HH(a); fX.CC(bx) = f.CC(a);
		fY.HH(by) = f.HH(a); fY.CC(by) = f.CC(a);
		fT.HH(bt) = f.HH(a); fT.CC(bt) = f.CC(a);
	}

	int checked = 0;
	for(int a = 0; a < 4; ++a) for(int b = 0; b < 4; ++b)
	for(int c = 0; c < 4; ++c) for(int d = 0; d < 4; ++d)
	{
		const int xa = kIO[a], xb = kIO[b], ya = kIO[c], yb = kIO[d];
		if(!Valid(xa, xb, ya, yb)) continue;
		if(xa != 3 && xb != 3 && ya != 3 && yb != 3) continue;	// an inlet somewhere

		BoundaryCondition2D bc, bcX, bcY, bcT;
		bc.Build (NX, NY, kTypes[xa], kTypes[xb], kTypes[ya], kTypes[yb]);
		bcX.Build(NX, NY, kTypes[xb], kTypes[xa], kTypes[ya], kTypes[yb]);
		bcY.Build(NX, NY, kTypes[xa], kTypes[xb], kTypes[yb], kTypes[ya]);
		bcT.Build(NY, NX, kTypes[ya], kTypes[yb], kTypes[xa], kTypes[xb]);

		for(int fi = 0; fi < 3; fi += 2)
		{
			std::vector<double> x  = Impose(bc,  kFields[fi], f,  f);
			std::vector<double> xX = Impose(bcX, kFields[fi], fX, fX);
			std::vector<double> xY = Impose(bcY, kFields[fi], fY, fY);
			std::vector<double> xT = Impose(bcT, kFields[fi], fT, fT);
			double worst = 0.0; int where = 0;
			for(int j = 0; j < NY; ++j) for(int i = 0; i < NX; ++i)
			{
				const double v = x[NX*j + i];
				const double e = std::max(std::max(std::fabs(v - xX[NX*j + (NX-1-i)]),
				                                   std::fabs(v - xY[NX*(NY-1-j) + i])),
				                          std::fabs(v - xT[NY*i + j]));
				if(e > worst) { worst = e; where = NX*j + i; }
			}
			if(worst > 1e-12)
			{
				const std::string combo = std::string(kName[xa]) + " " + kName[xb] + " " + kName[ya] + " " + kName[yb];
				CAPTURE(combo);
				CAPTURE(fi); CAPTURE(where % NX); CAPTURE(where / NX); CAPTURE(worst);
				CHECK(worst <= 1e-12);
			}
			++checked;
		}
	}
	CHECK(checked > 100);
}


TEST_CASE("BoundaryCondition2D: the Newton correction leaves every edge and corner where it is")
{
	//  The film-height equation applies the Height conditions to the right-
	//  hand side once, then the HeightJacobi ones before every Newton solve,
	//  into the same vector. The Jacobi rows ask for a correction that keeps
	//  each condition as it is, so their right-hand side must be zero. An
	//  Inlet-Inlet corner that left it untouched would hold the Height value
	//  r*h written there first: every iteration would move the corner by that
	//  much again, the iteration stagnate whatever the step and the run abort
	//  on MinDt (tools/proptest.py seed 1171).
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;

	Field_WM2D f;
	f.Build(MakeGrid("2D", 10, 8), subs, 1.0);
	C_Grid2D *G = f.GetGridPtr();
	int n = 0;
	for(Zface zf = G->First<Zface>(); zf != G->Last<Zface>(); ++zf, ++n)
	{
		f.HH(zf) = 1.0 + 0.01 * (n % 7);
		f.CC(zf) = 0.5;
	}
	const int NX = G->nx(), NY = G->ny();

	int checked = 0;
	for(int xa = 0; xa < 5; ++xa)
	for(int xb = 0; xb < 5; ++xb)
	for(int ya = 0; ya < 5; ++ya)
	for(int yb = 0; yb < 5; ++yb)
	{
		if((xa == 1) != (xb == 1)) continue;
		if((ya == 1) != (yb == 1)) continue;
		if(NotImplemented(xa,ya) || NotImplemented(xa,yb) ||
		   NotImplemented(xb,ya) || NotImplemented(xb,yb)) continue;

		CRS_Matrix crs; Vector<double> rhs;
		PrepareMatrix(crs, rhs, NX, NY);
		BoundaryCondition2D bc;
		bc.Build(NX, NY, kTypes[xa], kTypes[xb], kTypes[ya], kTypes[yb]);
		bc.Apply(crs, rhs, BoundaryCondition::Height, f);
		CRS_Matrix jac; Vector<double> unused;
		PrepareMatrix(jac, unused, NX, NY);
		bc.Apply(jac, rhs, BoundaryCondition::HeightJacobi, f);

		double worst = 0;
		for(int i = 0; i < NX * NY; ++i)
			if(jac.IsBC(i) && std::fabs(rhs[i]) > worst) worst = std::fabs(rhs[i]);
		if(worst != 0.0)
		{
			CAPTURE(kName[xa]); CAPTURE(kName[xb]);
			CAPTURE(kName[ya]); CAPTURE(kName[yb]);
			CHECK(worst == 0.0);
		}
		++checked;
	}
	CHECK(checked > 100);
}


//==========================================================================
//  CRS_Matrix::SatisfyBCRows on the rows of every condition.
//
//  After each Newton solve the boundary rows are re-satisfied from the
//  interior, so that a mirror ghost (Symmetry) is the exact mirror and a
//  periodic ghost the exact copy: the capillary flux of an edge face is a
//  third difference over them, and a ghost off its mirror by the tolerance
//  of the solve leaked volume (tools/analytic.py nozzle_balance). Starting
//  from ghosts that hold nothing in particular, every boundary row must
//  come out satisfied to rounding of its terms, the Symmetry and Periodic
//  ghosts bit for bit, at every edge and corner, for every field, and the
//  interior must not move.
//==========================================================================
namespace {

//  Largest |rhs_i - sum_j A_ij x_j| over the boundary rows, relative to
//  |rhs_i| + sum_j |A_ij x_j| of that row
double BCRowResidual(const CRS_Matrix &crs, const Vector<double> &rhs, const Vector<double> &x)
{
	double worst = 0;
	for(int i = 0; i < crs.RowSize(); ++i)
	{
		if(!crs.IsBC(i)) continue;
		double s = 0, m = std::fabs(rhs[i]);
		for(int k = 0; k < crs.NonZero(i); ++k)
		{
			const double t = crs.Element(i, k) * x[crs.Jcol(i, k)];
			s += t; m += std::fabs(t);
		}
		if(m > 0) worst = std::max(worst, std::fabs(rhs[i] - s) / m);
	}
	return worst;
}

//  What the ghosts hold before: nothing in particular
double Junk(int r) { return 0.7 + 0.5 * std::sin(12.9898 * r + 78.233); }

} // namespace

TEST_CASE("CRS_Matrix::SatisfyBCRows: every condition in 1-D, the ghosts exact")
{
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;
	Field_WM1D f;
	f.Build(MakeGrid("1D", 12, 1), subs, 1.0);
	C_Grid1D *G = f.GetGridPtr();
	for(Xedge xe = G->First<Xedge>(); xe != G->Last<Xedge>(); ++xe)
	{
		f.HH(xe) = Ph(G->ix(xe), 0);
		f.CC(xe) = Pc(G->ix(xe), 0);
	}
	const int NX = G->nx();

	int checked = 0;
	for(int xa = 0; xa < 5; ++xa)
	for(int xb = 0; xb < 5; ++xb)
	{
		if((xa == 1) != (xb == 1)) continue;
		for(int fi = 0; fi < 3; ++fi)
		{
			CRS_Matrix crs; Vector<double> rhs;
			PrepareMatrix(crs, rhs, NX, 1);
			BoundaryCondition1D bc;
			bc.Build(NX, kTypes[xa], kTypes[xb]);
			bc.Apply(crs, rhs, kFields[fi], f);

			Vector<double> x(NX);
			for(int r = 0; r < NX; ++r) x[r] = crs.IsBC(r) ? Junk(r) : Ph(r, 0);
			const Vector<double> x0 = x;
			crs.SatisfyBCRows(rhs, x);

			const std::string combo = std::string(kName[xa]) + " " + kName[xb];
			CAPTURE(combo); CAPTURE(fi);
			CHECK(BCRowResidual(crs, rhs, x) < 1e-15);
			bool interior = true;
			for(int r = 0; r < NX; ++r) if(!crs.IsBC(r) && x[r] != x0[r]) interior = false;
			CHECK(interior);
			if(xa == 0) { CHECK(x[0] == x[3]);       CHECK(x[1] == x[2]); }
			if(xb == 0) { CHECK(x[NX-1] == x[NX-4]); CHECK(x[NX-2] == x[NX-3]); }
			if(xa == 1)
			{
				CHECK(x[0] == x[NX-4]); CHECK(x[1] == x[NX-3]);
				CHECK(x[NX-1] == x[3]); CHECK(x[NX-2] == x[2]);
			}
			++checked;
		}
	}
	CHECK(checked == 51);
}

TEST_CASE("CRS_Matrix::SatisfyBCRows: every condition in 2-D, edges and corners")
{
	Input::I_Substrate subs;
	subs.TotalHeight = 1.0;
	Field_WM2D f;
	f.Build(MakeGrid("2D", 10, 8), subs, 1.0);
	Fill(f, Ph, Pc);
	C_Grid2D *G = f.GetGridPtr();
	const int NX = G->nx(), NY = G->ny();

	int checked = 0;
	for(int xa = 0; xa < 5; ++xa) for(int xb = 0; xb < 5; ++xb)
	for(int ya = 0; ya < 5; ++ya) for(int yb = 0; yb < 5; ++yb)
	{
		if(!Valid(xa, xb, ya, yb)) continue;
		for(int fi = 0; fi < 3; ++fi)
		{
			CRS_Matrix crs; Vector<double> rhs;
			PrepareMatrix(crs, rhs, NX, NY);
			BoundaryCondition2D bc;
			bc.Build(NX, NY, kTypes[xa], kTypes[xb], kTypes[ya], kTypes[yb]);
			bc.Apply(crs, rhs, kFields[fi], f);

			Vector<double> x(NX*NY);
			for(int r = 0; r < NX*NY; ++r) x[r] = crs.IsBC(r) ? Junk(r) : Ph(r % NX, r / NX);
			const Vector<double> x0 = x;
			crs.SatisfyBCRows(rhs, x);

			const double res = BCRowResidual(crs, rhs, x);
			bool interior = true;
			for(int r = 0; r < NX*NY; ++r) if(!crs.IsBC(r) && x[r] != x0[r]) interior = false;
			//  The mirror and the copy along each edge, outside the corners
			bool exact = true;
			for(int j = 2; j < NY-2; ++j)
			{
				const int o = NX*j;
				if(xa == 0 && (x[o+0] != x[o+3] || x[o+1] != x[o+2])) exact = false;
				if(xb == 0 && (x[o+NX-1] != x[o+NX-4] || x[o+NX-2] != x[o+NX-3])) exact = false;
				if(xa == 1 && (x[o+0] != x[o+NX-4] || x[o+1] != x[o+NX-3])) exact = false;
			}
			for(int i = 2; i < NX-2; ++i)
			{
				if(ya == 0 && (x[i] != x[i+3*NX] || x[i+NX] != x[i+2*NX])) exact = false;
				if(yb == 0 && (x[i+NX*(NY-1)] != x[i+NX*(NY-4)] || x[i+NX*(NY-2)] != x[i+NX*(NY-3)])) exact = false;
				if(ya == 1 && (x[i] != x[i+NX*(NY-4)] || x[i+NX] != x[i+NX*(NY-3)])) exact = false;
			}
			if(res >= 1e-15 || !interior || !exact)
			{
				const std::string combo = std::string(kName[xa]) + " " + kName[xb] + " " + kName[ya] + " " + kName[yb];
				CAPTURE(combo); CAPTURE(fi); CAPTURE(res);
				CHECK(res < 1e-15);
				CHECK(interior);
				CHECK(exact);
			}
			++checked;
		}
	}
	CHECK(checked == 411);	// 137 combinations, three fields each
}
