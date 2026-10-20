/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  EnergyBudget.h : the free energy of the film and its dissipation, for
*                   Output.Energy (Energy.dat), well-mixed model
*
*  The free energy, dimensionless (theory.md section 2.3),
*
*      F = int [ sigma/2 |grad(H + S)|^2 + Bo/2 (H + S)^2 + Dp w(H) ] dA,
*      w(H) = -int_{H*}^{H} Pi(s) ds,
*
*  is the functional whose variation is the driving potential with its sign
*  turned, dF/dH = -phi, phi = sigma lap(H + S) + Dp Pi(H) - Bo (H + S).
*  Discretely F is the sum over the cells of the last two terms times the
*  cell measure and, for the first, the sum over the faces of
*  sigma_f/2 ((f_j - f_i)/dx)^2 times the cell measure, f = H + S: its
*  derivative with respect to H_i is then -phi_i times the cell measure,
*  phi_i with the second difference of LubFields' ConvertPotential, as long
*  as sigma is uniform. The faces are the high face of every real cell and
*  the low face of a real cell next to a ghost on an edge that is not
*  periodic (on a periodic edge the seam is the high face of the last cell).
*
*  The dissipation is what the scheme itself dissipates: the power
*  sum_i phi_i R_i of the potential against the rate R_i the face fluxes of
*  every term give each cell (HeightEquationBase::FluxPower), times the cell
*  measure. For the capillary pressure, the disjoining pressure and the
*  vertical gravity, which share one face mobility (Mobility::Face), it is
*  the sum over the faces of D_f (Delta phi / dx)^2 >= 0, i.e.
*  int M/mu |grad phi|^2 to O(dx^2) (the slip mobility with Physics.Slip).
*  Then dF/dt = -D for the semi-discrete system on a closed (Symmetry) or
*  periodic domain, and with Crank-Nicolson
*
*      (F^{n+1} - F^n)/dt + (D^n + D^{n+1})/2 = O(dt^2).
*
*  F decreases only where the flux is the gradient of phi alone: capillary
*  pressure, vertical gravity and disjoining pressure, a uniform
*  concentration, closed or periodic edges. In-plane gravity, the
*  centrifugal force, the Marangoni flux, evaporation, the nozzle and open
*  edges put energy in or take it out; Energy.dat is written for them too,
*  and D is then the power of the potential against the flux, not a
*  dissipation.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#ifndef YUZEN_ENERGY_BUDGET_H
#define YUZEN_ENERGY_BUDGET_H

#include <cmath>
#include <vector>
#include "LubFields.h"
#include "NonDimensionalNumbers.h"

namespace energy_budget {

using namespace yuzen;

//  w(H) = -int_{H*}^{H} [(H*/s)^n - (H*/s)^m] ds (n, m > 1)
inline double W(double h, const PotentialParams &pp)
{
	if(pp.dp == 0 || !(h > 0)) return 0;
	const double hs = pp.hstar, n = pp.n, m = pp.m;
	const double in = std::pow(hs, n) * (std::pow(h, 1.0 - n) - std::pow(hs, 1.0 - n)) / (1.0 - n);
	const double im = std::pow(hs, m) * (std::pow(h, 1.0 - m) - std::pow(hs, 1.0 - m)) / (1.0 - m);
	return -pp.dp * (in - im);
}

inline double Sigma(double c, const PotentialParams &pp)
{
	return pp.nd ? pp.nd->SurfaceTension(c) : 1.0;
}

//  The energy of one cell without its faces: Bo/2 f^2 + Dp w(H)
inline double CellEnergy(double h, double f, const PotentialParams &pp)
{
	return 0.5 * pp.bo * f * f + W(h, pp);
}

//  The energy of a face between two values f of the free surface, dn apart
inline double FaceEnergy(double fi, double fj, double ci, double cj, double dn, const PotentialParams &pp)
{
	if(!pp.capillary) return 0;
	const double g = (fj - fi) / dn;
	return 0.25 * (Sigma(ci, pp) + Sigma(cj, pp)) * g * g;
}

//--------------------------------------------------------------------------
//  1-D. periodic: the X edges are Periodic (the Y ones do not enter). phi
//  is filled at the real cells, by the linear index of the grid.
//--------------------------------------------------------------------------
inline double FreeEnergy(Field_WM1D &ff, const PotentialParams &pp, bool periodic, bool /*periodicY*/)
{
	const C_Grid1D &G = ff.GetGrid();
	const double dx = G.dx();
	auto v = [&](Xedge e) { return ff.HH(e) + ff.SS(e); };
	double F = 0;
	for(int ii = G.imin(); ii <= G.imax(); ++ii)
	{
		Xedge xe = G.Of<Xedge>(ii);
		if(G.CheckGhost(xe) != C_Grid1D::Inside) continue;
		F += CellEnergy(ff.HH(xe), v(xe), pp);
		const Xedge xp = G.ShiftXp(xe), xm = G.ShiftXm(xe);
		F += FaceEnergy(v(xe), v(xp), ff.CC(xe), ff.CC(xp), dx, pp);
		if(!periodic && G.CheckGhost(xm) != C_Grid1D::Inside)
			F += FaceEnergy(v(xm), v(xe), ff.CC(xm), ff.CC(xe), dx, pp);
	}
	return F * dx;
}

inline void Potential(Field_WM1D &ff, const PotentialParams &pp, std::vector<double> &phi)
{
	const C_Grid1D &G = ff.GetGrid();
	const double dx2 = G.dx() * G.dx();
	phi.assign((size_t)G.imax() + 1, 0.0);
	auto v = [&](Xedge e) { return ff.HH(e) + ff.SS(e); };
	for(int ii = G.imin(); ii <= G.imax(); ++ii)
	{
		Xedge xe = G.Of<Xedge>(ii);
		if(G.CheckGhost(xe) != C_Grid1D::Inside) continue;
		const double lap = (v(G.ShiftXp(xe)) - 2.0 * v(xe) + v(G.ShiftXm(xe))) / dx2;
		phi[(size_t)ii] = yuzen::PotentialAt(ff.HH(xe), v(xe), ff.CC(xe), lap, pp);
	}
}

//--------------------------------------------------------------------------
//  2-D. periodicX, periodicY: the X, the Y edges are Periodic.
//--------------------------------------------------------------------------
inline double FreeEnergy(Field_WM2D &ff, const PotentialParams &pp, bool periodicX, bool periodicY)
{
	const C_Grid2D &G = ff.GetGrid();
	const double dx = G.dx(), dy = G.dy();
	auto v = [&](Zface e) { return ff.HH(e) + ff.SS(e); };
	double F = 0;
	for(int ii = G.imin(); ii <= G.imax(); ++ii)
	{
		Zface zf = G.Of<Zface>(ii);
		if(G.CheckGhost(zf) != C_Grid2D::Inside) continue;
		F += CellEnergy(ff.HH(zf), v(zf), pp);
		const Zface xp = G.ShiftXp(zf), xm = G.ShiftXm(zf), yp = G.ShiftYp(zf), ym = G.ShiftYm(zf);
		F += FaceEnergy(v(zf), v(xp), ff.CC(zf), ff.CC(xp), dx, pp);
		F += FaceEnergy(v(zf), v(yp), ff.CC(zf), ff.CC(yp), dy, pp);
		if(!periodicX && G.CheckGhost(xm) != C_Grid2D::Inside)
			F += FaceEnergy(v(xm), v(zf), ff.CC(xm), ff.CC(zf), dx, pp);
		if(!periodicY && G.CheckGhost(ym) != C_Grid2D::Inside)
			F += FaceEnergy(v(ym), v(zf), ff.CC(ym), ff.CC(zf), dy, pp);
	}
	return F * dx * dy;
}

inline void Potential(Field_WM2D &ff, const PotentialParams &pp, std::vector<double> &phi)
{
	const C_Grid2D &G = ff.GetGrid();
	const double dx2 = G.dx() * G.dx(), dy2 = G.dy() * G.dy();
	phi.assign((size_t)G.imax() + 1, 0.0);
	auto v = [&](Zface e) { return ff.HH(e) + ff.SS(e); };
	for(int ii = G.imin(); ii <= G.imax(); ++ii)
	{
		Zface zf = G.Of<Zface>(ii);
		if(G.CheckGhost(zf) != C_Grid2D::Inside) continue;
		const double vo = v(zf);
		const double lap = (v(G.ShiftXp(zf)) - 2.0 * vo + v(G.ShiftXm(zf))) / dx2
		                 + (v(G.ShiftYp(zf)) - 2.0 * vo + v(G.ShiftYm(zf))) / dy2;
		phi[(size_t)ii] = yuzen::PotentialAt(ff.HH(zf), vo, ff.CC(zf), lap, pp);
	}
}

}	// namespace energy_budget

#endif // YUZEN_ENERGY_BUDGET_H
