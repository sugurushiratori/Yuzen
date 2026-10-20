/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  WellMixed_Model_2D.h : the common part has moved to the template in WellMixedModelBase.h
*
*  All that remains here are the small differences specific to this
*  dimension.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _WELL_MIXED_MODEL_2D_INCLUDED_
#define _WELL_MIXED_MODEL_2D_INCLUDED_

#include "WellMixedModelBase.h"
#include "HeightEquation_2D.h"
#include "ConcEquation_2D.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//--------------------------------------------------------------------------
//  The driver of the 2-D well-mixed model.
//
//  Almost all of its body lives in WellMixed_ModelBase; all that remains
//  here is what differs with the dimension (how the boundary condition is
//  built, the representative length, the volume of a cell, how the probes
//  are located, and the column headings of the output).
//--------------------------------------------------------------------------
class WellMixed_Model_2D : public WellMixed_ModelBase<WM2D_Traits>
{
public:
	/** @brief Name of this model as it appears in the log. */
	virtual const char* Name() const { return "Well-Mixed-2D"; }

	WellMixed_Model_2D() : WellMixed_ModelBase<WM2D_Traits>() {}
	~WellMixed_Model_2D(){}

protected:
	/**
	 * @brief Build the boundary condition. There are four edges, in x and
	 *        in y, so four types are read.
	 */
	virtual void   BuildBoundaryCondition(Input &inp, Grid *G);
	/** @brief Both components, Ux and Uy, are zeroed in the initial velocity. */
	virtual void   ZeroInitialVelocity(Grid *G, const Site &s)
	{
		Site &z = const_cast<Site&>(s);
		_field.Ux(G->HalfXm(z)) = 0;
		_field.Uy(G->HalfYm(z)) = 0;
	}
	/**
	 * @brief The representative length used to cap the time step; in 2-D it
	 *        is the smaller of dx and dy.
	 */
	virtual double MinMeshSize(Grid *G) const { return Math::Min(G->dx(), G->dy()); }
	/**
	 * @brief The size of one cell used in the volume integral; in 2-D it is
	 *        the base area dx*dy.
	 */
	virtual double CellMeasure(Grid *G) const { return G->dx() * G->dy(); }
	/** @brief Log the grid information (Nx, Ny and the ranges of x and y). */
	virtual void   InfoGridExtent() const;
	/**
	 * @brief Find the cell nearest to each probe position given in the
	 *        input, searching in both x and y.
	 */
	virtual void   FindProbes(Grid *G);
	/**
	 * @brief The velocity components listed in the headings of the probe
	 *        output; 2-D has Ux and Uy.
	 */
	virtual void   AddVelocityLabels(std::vector<std::string> &labels, bool dimensioned) const
	{
		labels.push_back(dimensioned ? "Ux [m/s]" : "Ux");
		labels.push_back(dimensioned ? "Uy [m/s]" : "Uy");
	}
	/**
	 * @brief Write the (x, y) coordinates of the probe at the top of its
	 *        output file.
	 */
	virtual void   WriteProbeHeader(std::ostringstream &oss, ProbeList::iterator it) const;
	/**
	 * @brief Read an ASCII restart file; only the name of the function
	 *        called differs between the dimensions.
	 */
	virtual bool   ReadRestartAscii(Input::I_Restart &inp, double coef)
	{ return C_Restart::ReadAscii_WM2D(inp, _field, coef); }

private:
	WellMixed_Model_2D(const WellMixed_Model_2D&);
	void operator=(const WellMixed_Model_2D&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_WELL_MIXED_MODEL_2D_INCLUDED_
