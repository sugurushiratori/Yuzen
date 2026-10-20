/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LubFields.h : the lubrication-approximation fields (height, concentration, substrate, velocity) and their output-side copies
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _LUB_FIELDS_INCLUDED_
#define _LUB_FIELDS_INCLUDED_


#include <vector>
#include "Grids.h"
#include "Vector.h"
#include "Substrate.h"


#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

class FieldVisitor;
class C_Output;

//--------------------------------------------------------------------------
//  Common base of the fields.
//
//  Which fields are held depends on the dimension and on the model
//  (Field_WM1D / Field_WM2D / Field_CZ), so the only thing they have in common
//  is that they accept a visitor. Aggregation (FieldInfo) and probing (Probe)
//  are written on the visitor side, which collects the branching over field
//  types in one place.
//--------------------------------------------------------------------------
class LubField
{
public:
	LubField() : _constructed(false) {}
public:
	virtual ~LubField(){}
	/**
	 * @brief Accept a visitor.
	 *
	 * The concrete field types hold different members, so the only thing they
	 * share is this entry point. Aggregation (FieldInfo) and probing (Probe)
	 * are written on the visitor side, which keeps the per-field-type branching
	 * in one place.
	 */
	virtual void Accept(FieldVisitor &) = 0;

protected:
	bool _constructed;
private:
	LubField(const LubField& c);
	void operator=(const LubField& c);
};


/////////////////////////////

//--------------------------------------------------------------------------
//  The 1-D velocity field.
//
//  The velocity of each physical term is held separately, and at the end they
//  are combined into the "mixed" velocity. Keeping the terms apart is not only
//  for visualization: it is needed in order to check which term acts where.
//
//  *** The mixed velocity is the PLAIN SUM of the per-term velocities. ***
//  The governing equation is h_t = -div(Q); when it is expanded inside that
//  -div, only the gravity term comes out with a + sign and the rest with a -.
//  The _Magnitude of each term follows that convention, so the velocity field
//  each term writes already has its physical direction. The combination must
//  therefore be a plain sum. This velocity is also used to advect the
//  concentration, so getting the sign wrong reverses the direction in which
//  the solute is transported.
//  tests/test_mixed_velocity.cpp pins this down.
//--------------------------------------------------------------------------
class Velocity_WM1D
{
public:
	Velocity_WM1D() : _constructed(false) {}
	/** @brief Copy-construct a velocity field of the same shape as vv, allocating its own storage. */
	Velocity_WM1D(const Velocity_WM1D& vv) {Build(vv);}
	~Velocity_WM1D(){}

	/**
	 * @brief Allocate every component on the given grid.
	 * @param g  grid the field is laid out on; a pointer to it is retained, so it must outlive this object
	 */
	void Build(C_Grid1D &g);
	/** @brief Allocate a second velocity field on the same grid as vv and copy its values across. */
	void Build(const Velocity_WM1D& vv);

	void operator=(const Velocity_WM1D& vv);

	// The combined ("mixed") velocity: the effective flow speed, and the one the
	// concentration is advected by. Only valid after Cal_MixedVelocity() has run.
	const VectorField1D<>& UU() const {return _mixed;}
	      VectorField1D<>& UU()       {return _mixed;}
	const double &Ux(Node nd) const {return _mixed.x(nd);}
	      double &Ux(Node nd)       {return _mixed.x(nd);}

	//--------------------------------------------------------------
	//  One velocity field per physical term. Each kernel's CalVelocity() writes
	//  its own. They are written out for visualization (VTK) as they stand.
	//--------------------------------------------------------------
	// One velocity field per physical term. Each kernel's CalVelocity() writes
	// its own; each already carries its physical sign (see the note above), and
	// each is written out for visualization as it stands.
	const VectorField1D<> &U_LaplacePressure() const {return _LaplacePressure;}
	      VectorField1D<> &U_LaplacePressure()       {return _LaplacePressure;}
	const VectorField1D<> &U_SolutoCapillarity() const {return _SolutoCapillarity;}
	      VectorField1D<> &U_SolutoCapillarity()       {return _SolutoCapillarity;}
	// Velocity due to the disjoining pressure (wettability).
	const VectorField1D<> &U_DisjoiningPressure() const {return _DisjoiningPressure;}
	      VectorField1D<> &U_DisjoiningPressure()       {return _DisjoiningPressure;}
	const VectorField1D<> &U_VerticalGravity() const {return _VerticalGravity;}
	      VectorField1D<> &U_VerticalGravity()       {return _VerticalGravity;}
	const VectorField1D<> &U_HorizontalGravity() const {return _HorizontalGravity;}
	      VectorField1D<> &U_HorizontalGravity()       {return _HorizontalGravity;}
	const VectorField1D<> &U_Centrifugal() const {return _Centrifugal;}
	      VectorField1D<> &U_Centrifugal()       {return _Centrifugal;}

	//  Add the per-term velocities together to form the mixed velocity. It must
	//  be a **plain sum** (see the note above). Call it once, after every
	//  kernel's CalVelocity() has been called.
	/**
	 * @brief Combine the per-term velocities into the mixed velocity.
	 *
	 * The combination is a PLAIN SUM of the per-term fields. Each term's
	 * _Magnitude already follows the sign convention of h_t = -div(Q), so each
	 * per-term velocity already points the right way physically; do not write
	 * this as a chain of differences. The result also advects the solute, so a
	 * sign error here reverses solute transport. See the note on this class and
	 * tests/test_mixed_velocity.cpp, which pins the behavior down.
	 *
	 * Call once, after every kernel's CalVelocity() has run.
	 */
	void Cal_MixedVelocity();

private:
	bool _constructed;
	C_Grid1D *_G;
	VectorField1D<> _mixed;
	VectorField1D<> _LaplacePressure;
	VectorField1D<> _SolutoCapillarity;
	VectorField1D<> _DisjoiningPressure;
	VectorField1D<> _VerticalGravity;
	VectorField1D<> _HorizontalGravity;
	VectorField1D<> _Centrifugal;

};

//--------------------------------------------------------------------------
//  The 1-D (well-mixed) field. It holds the film height, the substrate shape,
//  the concentration and the velocity together.
//
//  The concentration is taken to be uniform across the film thickness, so it
//  is a single value living at the same place as the film height.
//  The grid is **owned by this field**, so there is no need to build one
//  separately.
//--------------------------------------------------------------------------
class Field_WM1D : public LubField
{
public:
	Field_WM1D();
	/** @brief Copy-construct a field of the same shape as f, with its own grid and storage. */
	Field_WM1D(const Field_WM1D& f) : LubField() {Build(f);}
	~Field_WM1D(){}

	/**
	 * @brief Build the grid, allocate the fields, and write the substrate shape into them.
	 * @param g       [Grid] input block
	 * @param i_subs  [Substrate] input block, describing the substrate relief
	 * @param scale   in-plane length scale L0 [m], used to non-dimensionalize the substrate shape
	 *
	 * The grid is owned by this object, so there is no need to build one separately.
	 */
	void Build(const Input::I_Grid &g, const Input::I_Substrate &i_subs, double scale);
	/** @brief Allocate a second field with its own copy of ff's grid and copy ff's values across (used for the intermediate Newton state _field_guess). */
	void Build(const Field_WM1D& ff);

	void operator=(const Field_WM1D& ff);

	/** @brief Accept a visitor (aggregation, probing). */
	virtual void Accept(FieldVisitor &v);

	//--------------------------------------------------------------
	//  Field access. In every case:
	//    no argument  ... the field itself (for running loops over it)
	//    pass a site  ... the value at that site (both readable and writable)
	//  Every value is **non-dimensional**. Restoring the dimensions is the
	//  output side's job.
	//--------------------------------------------------------------
	// Field access. Calling with no argument gives the whole field (for loops);
	// calling with a site gives that site's value, readable and writable.
	// Every value here is NON-DIMENSIONAL; restoring dimensions is the output
	// side's job.
	//
	// Film height h, measured from the substrate, so the free surface sits at h + s.
	const ScalarField1D<>& HH() const {return _height;}
	      ScalarField1D<>& HH()       {return _height;}
	const double &HH(Xedge xe) const {return _height(xe);}
	      double &HH(Xedge xe)       {return _height(xe);}

	// Substrate shape s.
	const ScalarField1D<>& SS() const {return _substrate;}
	      ScalarField1D<>& SS()       {return _substrate;}
	const double &SS(Xedge xe) const {return _substrate(xe);}
	      double &SS(Xedge xe)       {return _substrate(xe);}

	// Solute concentration c.
	const ScalarField1D<>& CC() const {return _conc;}
	      ScalarField1D<>& CC()       {return _conc;}
	const double &CC(Xedge xe) const {return _conc(xe);}
	      double &CC(Xedge xe)       {return _conc(xe);}

	// The whole velocity set. Ux() reaches the mixed velocity, which is only
	// meaningful once Velocity_WM1D::Cal_MixedVelocity() has been called.
	const Velocity_WM1D& UU() const {return _velocity;}
	      Velocity_WM1D& UU()       {return _velocity;}
	const double &Ux(Node nd) const {return _velocity.Ux(nd);}
	      double &Ux(Node nd)       {return _velocity.Ux(nd);}

	//--------------------------------------------------------------
	//  Write the solver's solution vector back into the field.
	//  The equation side holds the unknowns as a single vector, and these
	//  bridge the two.
	//--------------------------------------------------------------
	/**
	 * @brief Write the solver's solution vector back into the film-height field.
	 * @param hh  solution vector, indexed by matrix row number (C_Grid1D::IL ordering)
	 */
	void Update_H(Vector<double> &hh);
	/**
	 * @brief Write the solver's solution vector back into the concentration field.
	 * @param cc  solution vector, indexed by matrix row number
	 */
	void Update_C(Vector<double> &cc);
	/**
	 * @brief Recover the concentration from the conserved quantity Hr = h * c, as c = Hr / h.
	 * @param Hr  solution vector holding the conserved quantity, not the concentration itself
	 *
	 * With evaporation present it is more natural to solve for the conserved
	 * quantity, so this is the path used in that case.
	 */
	void Update_C_fromHr(Vector<double> &Hr);

	// The grid this field owns.
	const C_Grid1D& GetGrid() const {return _G;}
	      C_Grid1D* GetGridPtr()    {return &_G;}

private:
	C_Grid1D _G;
	ScalarField1D<> _height;
	ScalarField1D<> _substrate;
	ScalarField1D<> _conc;
	Velocity_WM1D _velocity;
};

//--------------------------------------------------------------------------
//  The output-side 1-D field.
//
//  There are two reasons why the computational field is not written out as it
//  stands.
//    - The grid is staggered, so the components live at different places. In
//      the output they should all be at the same point, so the face values are
//      averaged onto the nodes.
//    - The output should be dimensional (Output.Dimensioned), whereas the
//      computation is non-dimensional.
//  ConvertFromStg() performs that conversion.
//
//  The layout keeps exactly one ghost layer, so the number of rows is nx()-2.
//  The two ends are the boundary values themselves (the average of the two
//  points inside). The ASCII restart reads this layout verbatim, so changing
//  it makes reading and writing disagree
//  (pinned down by tests/test_restart.cpp).
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//  What the driving potential of the field output is made of
//  (Output.FieldData.WritePotential):
//
//      phi = sigma(c) lap(h + s) + Dp Pi(h) - Bo (h + s),
//      Pi(h) = (h*/h)^n - (h*/h)^m,
//
//  dimensionless, the pressure-like terms whose gradient drives the flux
//  through the mobility h^3 / 3 mu (the README's phi; a uniform phi carries no
//  flux). The body forces (in-plane gravity, centrifugal) and the Marangoni
//  stress are not gradients of it and are left out. lap is the second
//  difference the capillary term differences once more at each face.
//--------------------------------------------------------------------------
class NonDimensionalNumbers;

struct PotentialParams
{
	bool capillary = false;			// Physics.LaplacePressure
	double dp = 0, hstar = 0, n = 3, m = 2;	// Dp (0: off), H*, exponents
	double bo = 0;					// Bo (0: off)
	const NonDimensionalNumbers *nd = nullptr;	// for sigma(c)
	double scale = 1;				// sigma0 h0 / L0^2 when dimensioned, else 1
};
/** @brief The potential at a point: sigma(c) lap + Dp Pi(h) - Bo hs, of the
 *         film h, the free surface hs = h + s, the concentration c and the
 *         Laplacian lap of hs given (dimensionless; pp.scale not applied). */
double PotentialAt(double h, double hs, double c, double lap, const PotentialParams &pp);

class Field_WM1D_Out
{
public:
	Field_WM1D_Out();
	/** @brief Declared only; no definition is provided, so this overload cannot be linked. Use the default constructor followed by Build(). */
	Field_WM1D_Out(const Field_WM1D& ff);
	~Field_WM1D_Out(){}

	/**
	 * @brief Lay out the output arrays from a computational field.
	 * @param ff           source field
	 * @param L0           in-plane length scale [m]
	 * @param H0           film height scale [m]
	 * @param U0           velocity scale [m/s]
	 * @param dimensioned  if false the values are left non-dimensional and L0/H0/U0 are not applied
	 *
	 * This only sizes and records the conversion; ConvertFromStg() fills in the values.
	 */
	void Build(const Field_WM1D& ff, double L0, double H0, double U0, bool dimensioned);

	/**
	 * @brief Map an (in-plane, height-direction) index pair to a linear index into the output arrays.
	 * @param ix  output row number, 0 .. nx()-1 (not a grid index)
	 * @param kz  height-direction index, 0 .. nz()-1
	 */
	int IL(int ix, int kz) const {assert(0<=ix && ix <_nx && 0<=kz && kz<_nz); return (_nx)*kz + ix;}

	//--------------------------------------------------------------
	//  Reading from the output layout. The index is not a grid index but the
	//  **output row number** (0 to nx()-1).
	//  The values are dimensional or non-dimensional according to the
	//  dimensioned argument of Build().
	//--------------------------------------------------------------
	// Reading the output arrays. The index is the OUTPUT ROW NUMBER (0 .. nx()-1),
	// not a grid index, and the values are only filled in once ConvertFromStg()
	// has run.
	//   HH = film height h, SS = substrate shape s, HS = free-surface height h + s,
	//   CC = concentration, Ux = mixed velocity, xx = coordinate.
	// UNITS: with dimensioned = false these stay non-dimensional. With
	// dimensioned = true the heights (HH, SS, HS) come out in MICROMETERS, the
	// coordinate xx in MILLIMETERS and the velocities in m/s. CC carries no
	// dimension and is never rescaled.
	double HH(int ix) const { return _hh[ix];}
	double HS(int ix) const { return _hs[ix];}
	double SS(int ix) const { return _ss[ix];}
	double CC(int ix) const { return _cc[ix];}
	double Ux(int ix) const { return _ux[ix];}
	double xx(int i) const {assert(0<=i && i<_nx); return _xx[i];}

	// Per-term velocity contributions, same indexing as above. The suffix names
	// the term: Ca = Laplace (capillary) pressure, Ma = soluto-capillary
	// (Marangoni), Bo = vertical gravity, Pg = in-plane gravity, Sp = centrifugal.
	// A term disabled in the input stays zero.
	double Ux_Ca(int ix) const { return _ux_Ca[ix];}
	double Ux_Ma(int ix) const { return _ux_Ma[ix];}
	double Ux_Bo(int ix) const { return _ux_Bo[ix];}
	double Ux_Pg(int ix) const { return _ux_Pg[ix];}
	double Ux_Sp(int ix) const { return _ux_Sp[ix];}

	// Sizes of the output arrays. nx() is nx()-2 of the grid: one ghost layer is kept.
	int nx() const {return _nx;}
	int nz() const {return _nz;}

	/**
	 * @brief Fill the output arrays from the staggered computational field.
	 *
	 * Velocity components live on faces, so they are averaged onto the nodes to
	 * put every quantity at the same point, and the values are multiplied by the
	 * scales recorded by Build() when dimensioned output was requested. The two
	 * end rows are the boundary values themselves, taken as the average of the
	 * two points inside. The ASCII restart reads this layout verbatim, so
	 * changing it breaks read/write round-tripping (pinned by tests/test_restart.cpp).
	 */
	void ConvertFromStg(Field_WM1D &ff);

	/**
	 * @brief Fill the driving potential (PotentialParams) from the field, at
	 *        the same rows as the others; written as PHI from then on.
	 *        Call after ConvertFromStg().
	 */
	void ConvertPotential(Field_WM1D &ff, const PotentialParams &pp);
	/** @brief Whether the potential has been filled (ConvertPotential). */
	bool HasPotential() const {return _hasPhi;}
	double PHI(int ix) const {return _phi[ix];}
	double* PHI_ptr() {return _phi.data();}

private:
	std::vector<double> _phi;
	bool _hasPhi = false;
	C_Grid1D _G;
	int _nx, _nz;
	std::vector<double> _xx;
	std::vector<double> _hh, _ss;
	std::vector<double> _cc, _hs;
	std::vector<double> _ux;
	std::vector<double> _ux_Ca;
	std::vector<double> _ux_Ma;
	std::vector<double> _ux_Bo;
	std::vector<double> _ux_Pg;
	std::vector<double> _ux_Sp;

	double _L0, _H0, _U0;
	bool _dimensioned;


	double* H_ptr() {return _hh.data();}
	double* S_ptr() {return _ss.data();}
	double* HS_ptr() {return _hs.data();}
	double* C_ptr() {return _cc.data();}
	double* Ux_ptr() {return _ux.data();}

	double* Ux_Ca_ptr() {return _ux_Ca.data();}
	double* Ux_Ma_ptr() {return _ux_Ma.data();}
	double* Ux_Bo_ptr() {return _ux_Bo.data();}
	double* Ux_Pg_ptr() {return _ux_Pg.data();}
	double* Ux_Sp_ptr() {return _ux_Sp.data();}

	void ConvertVelocity(std::vector<double> &ux, VectorField1D<> &vv);
	void ConvertScalar(std::vector<double> &ff, ScalarField1D<> &ss);

	friend class C_Output;
	friend class OutputVTK;
};


////////////////////////////
//--------------------------------------------------------------------------
//  The 2-D velocity field. It is built the same way as the 1-D version, with
//  two components, x and y. That the mixed velocity is the plain sum of the
//  per-term velocities holds here too (see the note on Velocity_WM1D).
//--------------------------------------------------------------------------
class Velocity_WM2D
{
public:
	Velocity_WM2D(): _constructed(false) {};
	/** @brief Declared only; no definition is provided, so this overload cannot be linked. Use the default constructor followed by Build(). */
	Velocity_WM2D(const Velocity_WM2D& vv);
	~Velocity_WM2D(){}

	/**
	 * @brief Allocate every component on the given grid.
	 * @param g  grid the field is laid out on; a pointer to it is retained, so it must outlive this object
	 */
	void Build(C_Grid2D &g);
	/** @brief Allocate a second velocity field on the same grid as vv and copy its values across. */
	void Build(const Velocity_WM2D& vv);

	void operator=(const Velocity_WM2D& vv);

	// The combined ("mixed") velocity: the effective flow speed, and the one the
	// concentration is advected by. The x and y components sit at different
	// places on the staggered grid, hence the different site types. Only valid
	// after Cal_MixedVelocity() has run.
	const VectorField2D<>& UU() const {return _mixed;}
	      VectorField2D<>& UU()       {return _mixed;}
	const double &Ux(Yedge ye) const {return _mixed.x(ye);}
	      double &Ux(Yedge ye)       {return _mixed.x(ye);}
	const double &Uy(Xedge xe) const {return _mixed.y(xe);}
	      double &Uy(Xedge xe)       {return _mixed.y(xe);}

	// One velocity field per physical term. Each kernel's CalVelocity() writes
	// its own; each already carries its physical sign, and each is written out
	// for visualization as it stands.
	const VectorField2D<> &U_LaplacePressure() const {return _LaplacePressure;}
	      VectorField2D<> &U_LaplacePressure()       {return _LaplacePressure;}
	const VectorField2D<> &U_SolutoCapillarity() const {return _SolutoCapillarity;}
	      VectorField2D<> &U_SolutoCapillarity()       {return _SolutoCapillarity;}
	// Velocity due to the disjoining pressure (wettability).
	const VectorField2D<> &U_DisjoiningPressure() const {return _DisjoiningPressure;}
	      VectorField2D<> &U_DisjoiningPressure()       {return _DisjoiningPressure;}
	const VectorField2D<> &U_VerticalGravity() const {return _VerticalGravity;}
	      VectorField2D<> &U_VerticalGravity()       {return _VerticalGravity;}
	const VectorField2D<> &U_HorizontalGravity() const {return _HorizontalGravity;}
	      VectorField2D<> &U_HorizontalGravity()       {return _HorizontalGravity;}
	const VectorField2D<> &U_Centrifugal() const {return _Centrifugal;}
	      VectorField2D<> &U_Centrifugal()       {return _Centrifugal;}

	/**
	 * @brief Combine the per-term velocities into the mixed velocity, component by component.
	 *
	 * As in the 1D case the combination is a PLAIN SUM: each per-term velocity
	 * already points the right way physically, so this must not be written as a
	 * chain of differences. The result also advects the solute, so a sign error
	 * here reverses solute transport. See the note on Velocity_WM1D.
	 *
	 * Call once, after every kernel's CalVelocity() has run.
	 */
	void Cal_MixedVelocity();

private:
	bool _constructed;
	C_Grid2D *_G;
	VectorField2D<> _mixed;
	VectorField2D<> _LaplacePressure;
	VectorField2D<> _SolutoCapillarity;
	VectorField2D<> _DisjoiningPressure;
	VectorField2D<> _VerticalGravity;
	VectorField2D<> _HorizontalGravity;
	VectorField2D<> _Centrifugal;
};

//--------------------------------------------------------------------------
//  The 2-D (well-mixed) field. It is the x-y version of the 1-D one, and what
//  it holds and how it is used are the same. The x and y velocity components
//  live at different places (the grid is staggered).
//--------------------------------------------------------------------------
class Field_WM2D : public LubField
{
public:
	Field_WM2D();
	/** @brief Declared only; no definition is provided, so this overload cannot be linked. Use the default constructor followed by Build(). */
	Field_WM2D(const Field_WM2D& f);
	~Field_WM2D(){}

	/**
	 * @brief Build the grid, allocate the fields, and write the substrate shape into them.
	 * @param g       [Grid] input block
	 * @param i_subs  [Substrate] input block, describing the substrate relief
	 * @param scale   in-plane length scale L0 [m], used to non-dimensionalize the substrate shape
	 *
	 * The grid is owned by this object, so there is no need to build one separately.
	 */
	void Build(const Input::I_Grid &g, const Input::I_Substrate &i_subs, double scale);
	/** @brief Allocate a second field with its own copy of ff's grid and copy ff's values across (used for the intermediate Newton state). */
	void Build(const Field_WM2D& ff);

	void operator=(const Field_WM2D& ff);

	/** @brief Accept a visitor (aggregation, probing). */
	virtual void Accept(FieldVisitor &v);


	//--------------------------------------------------------------
	//  Field access, under the same convention as the 1-D version: no argument
	//  gives the field itself, passing a site gives the value at that site.
	//  Everything is non-dimensional.
	//--------------------------------------------------------------
	// Field access, same convention as the 1D version: no argument gives the
	// whole field, a site gives that site's value. Every value is non-dimensional.
	//
	// Film height h, measured from the substrate, so the free surface sits at h + s.
	const ScalarField2D<>& HH() const {return _height;}
	      ScalarField2D<>& HH()       {return _height;}
	const double &HH(Zface zf) const {return _height(zf);}
	      double &HH(Zface zf)       {return _height(zf);}

	// Substrate shape s.
	const ScalarField2D<>& SS() const {return _substrate;}
	      ScalarField2D<>& SS()       {return _substrate;}
	const double &SS(Zface zf) const {return _substrate(zf);}
	      double &SS(Zface zf)       {return _substrate(zf);}

	// Solute concentration c.
	const ScalarField2D<>& CC() const {return _conc;}
	      ScalarField2D<>& CC()       {return _conc;}
	const double &CC(Zface zf) const {return _conc(zf);}
	      double &CC(Zface zf)       {return _conc(zf);}

	// The whole velocity set. x and y components sit at different places on the
	// staggered grid, which is why Ux() takes a Yedge and Uy() takes an Xedge.
	const Velocity_WM2D& UU() const {return _velocity;}
	      Velocity_WM2D& UU()       {return _velocity;}
	const double &Ux(Yedge ye) const {return _velocity.Ux(ye);}
	      double &Ux(Yedge ye)       {return _velocity.Ux(ye);}
	const double &Uy(Xedge xe) const {return _velocity.Uy(xe);}
	      double &Uy(Xedge xe)       {return _velocity.Uy(xe);}

	// Write the solver's solution vector back into a field; same as the 1D
	// version. The vectors are indexed by matrix row number (C_Grid2D::IL
	// ordering). Update_C_fromHr() takes the conserved quantity Hr = h * c and
	// recovers c as Hr / h, which is the path used when evaporation is on.
	void Update_H(Vector<double> &hh);
	void Update_C(Vector<double> &cc);
	void Update_C_fromHr(Vector<double> &Hr);

	// The grid this field owns.
	const C_Grid2D& GetGrid() const {return _G;}
	      C_Grid2D* GetGridPtr()    {return &_G;}

private:
	C_Grid2D _G;
	ScalarField2D<> _height;
	ScalarField2D<> _substrate;
	ScalarField2D<> _conc;
	Velocity_WM2D _velocity;
};

/**
 * @brief Output-side copy of a 2D field, collocated at the nodes and optionally re-dimensionalized.
 *
 * Same role and the same layout convention as Field_WM1D_Out: one ghost layer
 * is kept on each side, so nx() is the grid's nx() minus 2 and the end rows hold
 * the boundary values themselves.
 */
class Field_WM2D_Out
{
public:
	Field_WM2D_Out();
	/** @brief Declared only; no definition is provided, so this overload cannot be linked. Use the default constructor followed by Build(). */
	Field_WM2D_Out(const Field_WM2D& ff);
	~Field_WM2D_Out(){}

	/**
	 * @brief Lay out the output arrays from a computational field.
	 * @param ff           source field
	 * @param L0           in-plane length scale [m]
	 * @param H0           film height scale [m]
	 * @param U0           velocity scale [m/s]
	 * @param dimensioned  if false the values are left non-dimensional and L0/H0/U0 are not applied
	 *
	 * This only sizes and records the conversion; ConvertFromStg() fills in the values.
	 */
	void Build(const Field_WM2D& ff, double L0, double H0, double U0, bool dimensioned);

	/**
	 * @brief Map output row indices to a linear index into the output arrays.
	 * The three-argument form takes a further height-direction index kz in
	 * 0 .. nz()-1; this two-argument form is that with kz = 0.
	 *
	 * @param ix  in-plane output index, 0 .. nx()-1 (an output row number, not a grid index)
	 * @param jy  in-plane output index, 0 .. ny()-1
	 */
	int IL(int ix, int jy) const {assert(0<=ix && ix <_nx && 0<=jy && jy<_ny); return (_nx* jy) + ix;}
	int IL(int ix, int jy, int kz) const {assert(0<=ix && ix <_nx && 0<=jy && jy<_ny && 0<=kz && kz<_nz); return (_nx*_ny)*kz + (_nx* jy) + ix;}

	// Reading the output arrays, indexed by output row number, not by grid index,
	// and only filled in once ConvertFromStg() has run.
	//   HH = film height h, SS = substrate shape s, HS = free-surface height h + s,
	//   CC = concentration, Ux/Uy = mixed velocity, xx/yy = coordinates.
	// UNITS: with dimensioned = false these stay non-dimensional. With
	// dimensioned = true the heights (HH, SS, HS) come out in MICROMETERS, the
	// coordinates xx and yy in MILLIMETERS and the velocities in m/s. CC carries
	// no dimension and is never rescaled.
	double HH(int ix, int jy) const { return _hh[IL(ix,jy)];}
	double HS(int ix, int jy) const { return _hs[IL(ix,jy)];}
	double SS(int ix, int jy) const { return _ss[IL(ix,jy)];}
	double CC(int ix, int jy) const { return _cc[IL(ix,jy)];}
	double Ux(int ix, int jy) const { return _ux[IL(ix,jy)];}
	double Uy(int ix, int jy) const { return _uy[IL(ix,jy)];}
	double xx(int i) const {assert(0<=i && i<_nx); return _xx[i];}
	double yy(int j) const {assert(0<=j && j<_ny); return _yy[j];}

	// Per-term velocity contributions, same indexing as above. The suffix names
	// the term: Ca = Laplace (capillary) pressure, Ma = soluto-capillary
	// (Marangoni), Bo = vertical gravity, Pg = in-plane gravity, Sp = centrifugal.
	// A term disabled in the input stays zero.
	double Ux_Ca(int ix, int jy) const { return _ux_Ca[IL(ix,jy)];}
	double Uy_Ca(int ix, int jy) const { return _uy_Ca[IL(ix,jy)];}
	double Ux_Ma(int ix, int jy) const { return _ux_Ma[IL(ix,jy)];}
	double Uy_Ma(int ix, int jy) const { return _uy_Ma[IL(ix,jy)];}
	double Ux_Bo(int ix, int jy) const { return _ux_Bo[IL(ix,jy)];}
	double Uy_Bo(int ix, int jy) const { return _uy_Bo[IL(ix,jy)];}
	double Ux_Pg(int ix, int jy) const { return _ux_Pg[IL(ix,jy)];}
	double Uy_Pg(int ix, int jy) const { return _uy_Pg[IL(ix,jy)];}
	double Ux_Sp(int ix, int jy) const { return _ux_Sp[IL(ix,jy)];}
	double Uy_Sp(int ix, int jy) const { return _uy_Sp[IL(ix,jy)];}

	// Sizes of the output arrays; nx() and ny() are the grid's counts minus 2,
	// one ghost layer being kept on each side.
	int nx() const {return _nx;}
	int ny() const {return _ny;}
	int nz() const {return _nz;}

	/**
	 * @brief Fill the output arrays from the staggered computational field.
	 *
	 * The x and y velocity components live on different faces, so both are
	 * averaged onto the nodes to put every quantity at the same point, and the
	 * values are multiplied by the scales recorded by Build() when dimensioned
	 * output was requested. The two end rows in each direction are the boundary
	 * values themselves, taken as the average of the two points inside.
	 */
	void ConvertFromStg(Field_WM2D &ff);

	/** @brief As Field_WM1D_Out::ConvertPotential. */
	void ConvertPotential(Field_WM2D &ff, const PotentialParams &pp);
	bool HasPotential() const {return _hasPhi;}
	double PHI(int ix, int jy) const {return _phi[IL(ix,jy)];}
	double* PHI_ptr() {return _phi.data();}

private:
	std::vector<double> _phi;
	bool _hasPhi = false;
	C_Grid2D _G;
	int _nx, _ny, _nz;
	std::vector<double> _xx, _yy;
	std::vector<double> _hh, _ss;
	std::vector<double> _cc, _hs;
	std::vector<double> _ux, _uy;
	std::vector<double> _ux_Ca, _uy_Ca;
	std::vector<double> _ux_Ma, _uy_Ma;
	std::vector<double> _ux_Bo, _uy_Bo;
	std::vector<double> _ux_Pg, _uy_Pg;
	std::vector<double> _ux_Sp, _uy_Sp;

	double _L0, _H0, _U0;
	bool _dimensioned;


	double* H_ptr() {return _hh.data();}
	double* S_ptr() {return _ss.data();}
	double* HS_ptr() {return _hs.data();}
	double* C_ptr() {return _cc.data();}
	double* Ux_ptr() {return _ux.data();}
	double* Uy_ptr() {return _uy.data();}

	double* Ux_Ca_ptr() {return _ux_Ca.data();}
	double* Uy_Ca_ptr() {return _uy_Ca.data();}
	double* Ux_Ma_ptr() {return _ux_Ma.data();}
	double* Uy_Ma_ptr() {return _uy_Ma.data();}
	double* Ux_Bo_ptr() {return _ux_Bo.data();}
	double* Uy_Bo_ptr() {return _uy_Bo.data();}
	double* Ux_Pg_ptr() {return _ux_Pg.data();}
	double* Uy_Pg_ptr() {return _uy_Pg.data();}
	double* Ux_Sp_ptr() {return _ux_Sp.data();}
	double* Uy_Sp_ptr() {return _uy_Sp.data();}

	void ConvertVelocity(std::vector<double> &ux, std::vector<double> &uy, VectorField2D<> &vv);
	void ConvertScalar(std::vector<double> &ff, ScalarField2D<> &ss);

	friend class C_Output;
	friend class OutputVTK;
};

//--------------------------------------------------------------------------
//  The field of the ConcZ model.
//
//  The film height lives on the in-plane 2-D grid, but the concentration lives
//  on a 3-D grid that is divided in the film-thickness direction as well.
//  The field therefore holds two grids (GetGridPtr2 / GetGridPtr3).
//  The coordinate in the film-thickness direction is normalized by the film
//  thickness, so the number of divisions does not change as the film thins
//  (the grid shrinks with it).
//--------------------------------------------------------------------------
class Field_CZ : public LubField
{
public:
	Field_CZ();
	/** @brief Copy-construct a field of the same shape as ff, with its own pair of grids and its own storage. */
	Field_CZ(const Field_CZ& ff) : LubField() {Build(ff);}
	~Field_CZ(){}

	/**
	 * @brief Build both grids, allocate the fields, and write the substrate shape into them.
	 * @param g       [Grid] input block
	 * @param i_subs  [Substrate] input block, describing the substrate relief
	 * @param scale   in-plane length scale L0 [m], used to non-dimensionalize the substrate shape
	 *
	 * Both the in-plane 2D grid and the 3D grid carrying the concentration are
	 * owned by this object.
	 */
	void Build(const Input::I_Grid &g, const Input::I_Substrate &i_subs, double scale);
	/** @brief Allocate a second field with its own copies of ff's two grids and copy ff's values across (used for the intermediate Newton state). */
	void Build(const Field_CZ& ff);

	void operator=(const Field_CZ& ff);

	/** @brief Accept a visitor (aggregation, probing). */
	virtual void Accept(FieldVisitor &v);

	//--------------------------------------------------------------
	//  Field access.
	//  The film height and the substrate live on the in-plane 2-D grid (Zface);
	//  the concentration and the velocity live on the 3-D grid that is divided
	//  in the film-thickness direction as well (Cell / Xface / Yface / Zface).
	//--------------------------------------------------------------
	// Field access, same convention as the well-mixed fields: no argument gives
	// the whole field, a site gives that site's value. Every value is
	// non-dimensional, and the z coordinate is normalized by the film thickness.
	//
	// Film height h, on the in-plane grid. The free surface sits at h + s.
	const ScalarField2D<>& HH() const {return _height;}
	      ScalarField2D<>& HH()       {return _height;}
	const double &HH(Zface zf) const {return _height(zf);}
	      double &HH(Zface zf)       {return _height(zf);}

	// Substrate shape s, on the in-plane grid.
	const ScalarField2D<>& SS() const {return _substrate;}
	      ScalarField2D<>& SS()       {return _substrate;}
	const double &SS(Zface zf) const {return _substrate(zf);}
	      double &SS(Zface zf)       {return _substrate(zf);}

	// Solute concentration c, resolved through the film thickness too, so it
	// lives on the 3D grid and is read per cell rather than per in-plane site.
	const ScalarField3D<>& CC() const {return _conc;}
	      ScalarField3D<>& CC()       {return _conc;}
	const double &CC(Cell cl) const {return _conc(cl);}
	      double &CC(Cell cl)       {return _conc(cl);}

	// Velocity on the 3D grid. The three components sit on the three face types
	// of the staggered grid, hence the different site types; Uz is the
	// through-thickness component, in the normalized z coordinate.
	const VectorField3D<>& UU() const {return _velocity;}
	      VectorField3D<>& UU()       {return _velocity;}
	const double &Ux(Xface xf) const {return _velocity.x(xf);}
	      double &Ux(Xface xf)       {return _velocity.x(xf);}
	const double &Uy(Yface yf) const {return _velocity.y(yf);}
	      double &Uy(Yface yf)       {return _velocity.y(yf);}
	const double &Uz(Zface zf) const {return _velocity.z(zf);}
	      double &Uz(Zface zf)       {return _velocity.z(zf);}

	// Write the solver's solution vector back into a field. hh is indexed by
	// C_Grid2D::IL row numbers, cc by C_Grid3D::IL row numbers. There is no
	// Update_C_fromHr() here: with the concentration resolved in z there is no
	// depth-integrated conserved quantity to recover it from.
	void Update_H(Vector<double> &hh);
	void Update_C(Vector<double> &cc);

	// The two grids this field owns: the in-plane grid the film height lives on,
	// and the 3D grid the concentration and velocity live on.
	const C_Grid2D& GetGrid2() const {return _G2;}
	      C_Grid2D* GetGridPtr2()    {return &_G2;}
	const C_Grid3D& GetGrid3() const {return _G3;}
	      C_Grid3D* GetGridPtr3()    {return &_G3;}

private:
	C_Grid2D _G2;
	C_Grid3D _G3;
	ScalarField2D<> _height;
	ScalarField2D<> _substrate;
	ScalarField3D<> _conc;
	VectorField3D<> _velocity;
};


#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_LUB_FIELDS_INCLUDED_
