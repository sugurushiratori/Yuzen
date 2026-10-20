/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Restart.h : reading and writing restarts (binary VTK / ASCII Neutral)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _DEF_RESTART_INCLUDED_
#define _DEF_RESTART_INCLUDED_

#include "YuzenUtil.h"

#include "Input.h"


#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class Field_WM1D;
class Field_WM2D;
class Field_CZ;

//--------------------------------------------------------------------------
//  Reading a restart.
//
//  Two formats can be read. Both take an earlier output file as the input, as
//  it stands.
//    binary : .vts (written by OutputVTK). Every digit is preserved.
//    ASCII  : .dat (written by OutputNeutral). A person can edit it, so it can
//             also be used to feed an analytical solution in as an initial
//             condition.
//
//  ** Correspondence between rows and grid points **
//  The ASCII order is exactly the order Field_WM1D_Out / Field_WM2D_Out
//  produce: the number of rows is nx()-2 (the two ends are the boundary
//  values), and row i goes into the grid x index i+1. If reading and writing
//  drift apart, the failure is a hard one to spot: the shape moves by one cell
//  every time the output is read back in. tests/test_restart.cpp pins this
//  down.
//
//  The boundary conditions are applied after reading, so the values in the
//  ghost layers are overwritten.
//--------------------------------------------------------------------------
class C_Restart
{
public:
	C_Restart(){}

	/**
	 * @brief Reads a binary restart (.vts) into the field.
	 *
	 * The file is one written earlier by OutputVTK; the header element it looks
	 * for is <YuzenData .../>. There is one overload per kind of field.
	 *
	 * @param inp     the [Restart] settings, which carry the file name
	 * @param field   destination field, overwritten in place
	 * @param H_coef  factor applied to the thicknesses as they are read, to get
	 *                from the units in the file to the ones in the field. Pass
	 *                1e-6/h0 for a file written in dimensional form, 1 for one
	 *                already non-dimensional.
	 * @return false if the file is missing or cannot be interpreted; the field
	 *         is then left as it was.
	 */
	//  Reads from binary (.vts). The header element is <YuzenData>.
	//  H_coef is the conversion factor applied to the film thickness (1e-6/h0
	//  if the file was written in dimensional form).
	static bool Read(Input::I_Restart &inp, Field_WM1D &field, double H_coef);
	static bool Read(Input::I_Restart &inp, Field_WM2D &field, double H_coef);
	static bool Read(Input::I_Restart &inp, Field_CZ &field, double H_coef);

	/**
	 * @brief Reads an ASCII restart (.dat) into the field.
	 *
	 * The file is one written earlier by OutputNeutral, in the row order noted
	 * above; since it is plain text it can also be hand-written, which is how
	 * an analytical solution is fed in as an initial condition.
	 *
	 * @param inp     the Restart block of the input file, giving the file to read
	 * @param field   field the values are read into
	 * @param H_coef  factor applied to the thicknesses as they are read, as for
	 *                the binary Read
	 * @return false if the file is missing, or if its row count does not match
	 *         the grid; the field is then left as it was.
	 */
	//  Reads from ASCII (.dat). Returns false if the number of grid points
	//  does not match.
	static bool ReadAscii_WM1D(Input::I_Restart &inp, Field_WM1D &field, double H_coef);
	static bool ReadAscii_WM2D(Input::I_Restart &inp, Field_WM2D &field, double H_coef);
	static bool ReadAscii_CZ(Input::I_Restart &inp, Field_CZ &field, double H_coef);

private:
	C_Restart(const C_Restart&);
	void operator=(const C_Restart&);
};


#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif


#endif //_DEF_RESTART_INCLUDED_
