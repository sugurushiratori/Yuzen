/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  FiniteDifferences.h : coefficients of the difference and interpolation schemes (central, upwind, 4-point)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _FINITE_DIFFERENCES_INCLUDED_
#define _FINITE_DIFFERENCES_INCLUDED_


#include "MyMath.h"
#include "Vector.h"
#include <string>

#include "grid/YuzenGrid.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

using namespace yuzen::grid;


///////////////////////////////////////////////
//   Interpolation onto the midpoint of a 4-point stencil
//     @     @   +   @     @
//     LL    L   o   R     RR
//////////////////////////////////////////////
//  This is where the upwinding of the convective term is absorbed. Changing
//  the scheme leaves the calling code untouched; only the coefficient table
//  (_c1 / _c2) is swapped out.
class C_Interpolation_base4
{
public:
	//  Positions of the four stencil points. Seen from the midpoint being
	//  interpolated, LL < L < midpoint < R < RR.
	enum t_iShift{LL=0, L = 1, R = 2, RR = 3};
	enum t_Scheme
	{
		Upwind_1st,		// 1st-upwind scheme.
		Upwind_2nd,		// 2nd-upwind scheme.
		Central_2nd,	// 2nd-central scheme.
		Upwind_3rd,		// 3rd-upwind scheme.
		QUICK,			// QUICK scheme by Leonard.(Quadratic Upstream Interpolation for Convective Kinematics)
		KK,				// Kawamura-Kuwaraha
		Central_4th,	// 4th-central scheme.
		Koren,			// TVD: kappa = 1/3 (= Upwind_3rd) where smooth, limited by Koren (1993)
		VanLeer			// TVD: MUSCL with the smooth limiter of van Leer (1974)
	};
	//  How a limited face value is written as a sum of weights times the
	//  point values, with the limiter frozen (see FaceWeights).
	enum t_FrozenForm
	{
		DownwindForm,	// f_U + psi/2 (f_D - f_U),   psi = Phi / (f_D - f_U)
		UpwindForm		// f_U + chi/2 (f_U - f_UU),  chi = Phi / (f_U - f_UU)
	};
public:
	/** @brief Construct with the 2nd-order central scheme and build its coefficients. */
	C_Interpolation_base4() : _n_stencil(4), _t_scheme(Central_2nd) {SetCoef();}

public:
	virtual ~C_Interpolation_base4() {}

	/**
	 * @brief Select a scheme and rebuild the coefficient table.
	 * @param t  the scheme to use from now on
	 */
	virtual void Build(t_Scheme t) {_t_scheme = t; SetCoef();}

	/** @brief Number of stencil points, which is 4. */
	int NumStencils() const {return _n_stencil;}

	/**
	 * @brief Interpolation coefficients for stencil point i.
	 * @param i  stencil position, one of LL / L / R / RR
	 * @return c1 is the sign-independent (central) part and c2 the
	 *         sign-dependent (upwinding) part.
	 *
	 * With sign as the flow direction (+1 or -1), the interpolated value is
	 *     sum_i ( c1(i) + sign * c2(i) ) * phi_i
	 * Upwinding appears only in c2, so for a central scheme every c2 is
	 * zero. Because the direction never needs an if, the inner loop stays
	 * branch free.
	 */
	inline double c1(int i){return _c1[i];}
	inline double c2(int i){return _c2[i];}

	/** @brief The scheme currently selected. */
	t_Scheme Scheme() const {return _t_scheme;}
	/** @brief Name of the current scheme, for writing to the log. */
	std::string SchemeName();

	//--------------------------------------------------------------
	//  The limited (TVD) schemes
	//
	//  The face value is not a fixed combination of the four points but
	//
	//      f_face = f_U + Phi(a, b) / 2,   a = f_U - f_UU,  b = f_D - f_U,
	//
	//  U the upwind point next to the face, UU the one behind it and D the
	//  downwind one: (UU, U, D) = (LL, L, R) for sign >= 0 and (RR, R, L)
	//  for sign < 0. Phi(a, b) = psi(r) b with r = a / b is the limiter
	//  written in the two differences, which needs no division by b and no
	//  epsilon at a zero gradient: it is 0 unless a and b have the same
	//  sign, and odd in (a, b) together. The mirror image of a stencil
	//  (x -> -x with the flow reversed) gives the same a and b and so the
	//  same face value to the last bit.
	//
	//    Koren:   Phi = sign(b) min(2|a|, (2|b| + |a|)/3, 2|b|)
	//             (psi = max(0, min(2r, (2+r)/3, 2)); in 0.4 < r < 4 the
	//             face value is (-f_UU + 5 f_U + 2 f_D)/6, i.e. Upwind_3rd)
	//    VanLeer: Phi = 2ab / (a + b)
	//             (psi = (r+|r|)/(1+|r|), smooth for r > 0)
	//
	//  Both keep 0 <= psi <= 2 and 0 <= psi/r <= 2, the TVD region of
	//  Sweby (1984).
	//
	//  c1()/c2() hold, for a limited scheme, the linear scheme it reduces
	//  to in smooth regions (Koren: Upwind_3rd; VanLeer: kappa = 0, Fromm),
	//  for GradX and for anything that asks for fixed coefficients.
	//--------------------------------------------------------------
	/** @brief Whether the face value depends on the data (Koren, VanLeer). */
	bool IsLimited() const {return _t_scheme == Koren || _t_scheme == VanLeer;}

	/**
	 * @brief The limiter in two-difference form, Phi(a, b) = psi(a/b) b.
	 * @param a       the upstream difference, f_U - f_UU
	 * @param b       the difference across the face, f_D - f_U
	 * @param dA, dB  [out, optional] dPhi/da and dPhi/db (one-sided at a
	 *                kink of Koren and at a = 0 or b = 0; 0 where Phi = 0)
	 */
	double Limiter(double a, double b, double *dA = 0, double *dB = 0) const;

	/**
	 * @brief The face value from the four points.
	 * @param f     the values at LL, L, R, RR
	 * @param sign  the direction of the flow through the face (its sign only)
	 */
	double FaceValue(const double f[4], double sign) const;

	/**
	 * @brief The face value as weights times the four point values, with
	 *        the limiter frozen at the data f: sum w_j f_j = FaceValue(f).
	 * @param f       the values at LL, L, R, RR
	 * @param sign    the direction of the flow through the face (its sign only)
	 * @param w       [out] the weights of LL, L, R, RR
	 * @param rowPos  the stencil position (L or R) of the cell whose row the
	 *                weights go into, or -1. When it is the upwind point,
	 *                the upwind form is used (w_UU = -chi/2 <= 0, w_U =
	 *                1 + chi/2), otherwise the downwind form (w_U = 1 -
	 *                psi/2, w_D = psi/2 >= 0). With these, the advection
	 *                u (f_R - f_L)/dx of a row is (c_i - c_{i-1}) times a
	 *                coefficient in [0, 2] u/dx, so that frozen matrix is
	 *                an M-matrix and its explicit part positive for
	 *                CFL <= 1. A linear scheme returns c1 + sign c2.
	 */
	void FaceWeights(const double f[4], double sign, int rowPos, double w[4]) const;

	/**
	 * @brief The face value and its exact derivative d f_face / d f_j (for
	 *        the Jacobian). A linear scheme returns c1 + sign c2.
	 */
	double FaceDerivative(const double f[4], double sign, double g[4]) const;

	/**
	 * @brief The scheme named by an input key (Numerics.Convection):
	 *        "KK", "UTOPIA", "Central4", "Koren", "VanLeer", "QUICK",
	 *        "Upwind1", "Upwind2", "Central2".
	 * @return false if the name is not one of them
	 */
	static bool FromKey(const std::string &key, t_Scheme &t);
	/** @brief The scheme named by key, or def when key is empty (Input read
	 *         it and refused anything else). */
	static t_Scheme FromKeyOr(const std::string &key, t_Scheme def)
	{ t_Scheme t = def; if(!key.empty()) FromKey(key, t); return t; }
	/** @brief The input key of a scheme (the inverse of FromKey). */
	static std::string Key(t_Scheme t);
	/** @brief The accepted keys, as a list for a message. */
	static std::string KeyList();

protected:
	const int _n_stencil;
	Vector<double> _c1;	//the sign-independent part (central)
	Vector<double> _c2;	//the sign-dependent part (upwinding), multiplied by Sign(sign)
	t_Scheme _t_scheme;

	/** @brief Fill the coefficient table _c1 / _c2 for the selected scheme. */
	void SetCoef();
private:
	C_Interpolation_base4(const C_Interpolation_base4& c);
	C_Interpolation_base4& operator=(const C_Interpolation_base4& c);
};


////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------
//  Takes the gradient of a field with the selected scheme.
//
//  Even for one and the same GradX, the stencil depends on which kind of
//  point the gradient is evaluated at, so there is an overload per point
//  type. A mismatched type does not compile, so a staggered position cannot
//  be mixed up.
//
//  sign is the flow direction (+1 / -1). Schemes that do not upwind ignore it.
//--------------------------------------------------------------------------
class C_FiniteDifferences : public C_Interpolation_base4
{
public:
	/** @brief Construct with the scheme C_Interpolation_base4 defaults to. */
	C_FiniteDifferences(){}

public:
	/**
	 * @brief Gradient of a scalar field, evaluated with the selected scheme.
	 * @param s     the field to differentiate
	 * @param nd    the node the gradient is evaluated at
	 * @param sign  flow direction, +1 or -1. Ignored by schemes that do not
	 *              upwind.
	 * @return The derivative at the given point.
	 *
	 * The stencil depends on the kind of point the gradient is evaluated at,
	 * so these are overloaded per point type (Node, Xedge, Zface, Cell and
	 * so on). A mismatched staggered position does not compile, so it cannot
	 * be mixed up at run time.
	 */
	//  1D
	double GradX(const ScalarField1D<> &s, const Node  &nd, double sign = 1.0);
	double GradX(const ScalarField1D<> &s, const Xedge &xe, double sign = 1.0);

	//  2D
	double GradX(const ScalarField2D<> &s, const Zface &zf, double sign = 1.0);
	double GradX(const ScalarField2D<> &s, const Yedge &ye, double sign = 1.0);
	double GradY(const ScalarField2D<> &s, const Zface &zf, double sign = 1.0);
	double GradY(const ScalarField2D<> &s, const Xedge &xe, double sign = 1.0);

	//  3D (the ConcZ model)
	double GradX(const ScalarField3D<> &s, const Cell  &cl, double sign = 1.0);
	double GradX(const ScalarField3D<> &s, const Xface &xf, double sign = 1.0);
	double GradY(const ScalarField3D<> &s, const Cell  &cl, double sign = 1.0);
	double GradY(const ScalarField3D<> &s, const Yface &yf, double sign = 1.0);
	double GradZ(const ScalarField3D<> &s, const Cell  &cl, double sign = 1.0);
	double GradZ(const ScalarField3D<> &s, const Zface &zf, double sign = 1.0);

private:
	C_FiniteDifferences(const C_FiniteDifferences& c);
	void operator=(const C_FiniteDifferences& c);
};


#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_FINITE_DIFFERENCES_INCLUDED_
