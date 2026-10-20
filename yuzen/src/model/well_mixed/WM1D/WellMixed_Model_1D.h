/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  WellMixed_Model_1D.h : the common part has moved to the template in WellMixedModelBase.h
*
*  All that remains here are the small differences specific to this
*  dimension.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _WELL_MIXED_MODEL_1D_INCLUDED_
#define _WELL_MIXED_MODEL_1D_INCLUDED_

#include "WellMixedModelBase.h"
#include "HeightEquation_1D.h"
#include "ConcEquation_1D.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//--------------------------------------------------------------------------
//  The driver of the 1-D well-mixed model.
//
//  Almost all of its body lives in WellMixed_ModelBase; all that remains
//  here is what differs with the dimension (how the boundary condition is
//  built, the representative length, the volume of a cell, how the probes
//  are located, and the column headings of the output).
//--------------------------------------------------------------------------
class WellMixed_Model_1D : public WellMixed_ModelBase<WM1D_Traits>
{
public:
	/** @brief Name of this model as it appears in the log. */
	virtual const char* Name() const { return "Well-Mixed-1D"; }

	WellMixed_Model_1D() : WellMixed_ModelBase<WM1D_Traits>() {}
	~WellMixed_Model_1D(){}

protected:
	/**
	 * @brief Build the boundary condition. There are boundaries only at the
	 *        two ends in x, so only two types are read.
	 */
	virtual void   BuildBoundaryCondition(Input &inp, Grid *G);
	/** @brief Only Ux is zeroed in the initial velocity. */
	virtual void   ZeroInitialVelocity(Grid *G, const Site &s) { _field.Ux(G->HalfXm(const_cast<Site&>(s))) = 0; }
	/**
	 * @brief The representative length used to cap the time step; in 1-D it
	 *        is dx itself.
	 */
	virtual double MinMeshSize(Grid *G) const { return G->dx(); }
	/**
	 * @brief The size of one cell used in the volume integral; in 1-D it is
	 *        the width dx.
	 */
	virtual double CellMeasure(Grid *G) const { return G->dx(); }
	/** @brief Log the grid information (Nx and the range of x). */
	virtual void   InfoGridExtent() const;
	/**
	 * @brief Find the cell nearest to each probe position given in the
	 *        input, searching in x only.
	 */
	virtual void   FindProbes(Grid *G);
	/**
	 * @brief The velocity components listed in the headings of the probe
	 *        output; 1-D has Ux only.
	 */
	virtual void   AddVelocityLabels(std::vector<std::string> &labels, bool dimensioned) const
	{ labels.push_back(dimensioned ? "Ux [m/s]" : "Ux"); }
	/**
	 * @brief Write the x coordinate of the probe at the top of its output
	 *        file.
	 */
	virtual void   WriteProbeHeader(std::ostringstream &oss, ProbeList::iterator it) const;
	/**
	 * @brief Read an ASCII restart file; only the name of the function
	 *        called differs between the dimensions.
	 */
	virtual bool   ReadRestartAscii(Input::I_Restart &inp, double coef)
	{ return C_Restart::ReadAscii_WM1D(inp, _field, coef); }

private:
	WellMixed_Model_1D(const WellMixed_Model_1D&);
	void operator=(const WellMixed_Model_1D&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_WELL_MIXED_MODEL_1D_INCLUDED_
