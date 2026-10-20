/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  MyMath.h : small math utilities - numerical comparison, rounding, constants
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _MY_MATH_INCLUDED_
#define _MY_MATH_INCLUDED_

#include <cmath>
#include <limits>

#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

//--------------------------------------------------------------------------
//  Small math utilities.
//
//  Only what std does not provide, or what is overloaded per type here so
//  that the wrong type cannot be picked up by mistake at the call site.
//  Everything is static and inline (the bodies are in MyMath_I.h).
//--------------------------------------------------------------------------
class Math
{
public:
	/**
	 * @brief Absolute value.
	 *
	 * Overloaded explicitly for int and double so that an implicit
	 * conversion cannot silently cost precision at the call site.
	 */
	static inline int    abs(const int    &val);
	static inline double abs(const double &val);

	//--------------------------------------------------------------
	//  Floating-point comparison, so that == does not have to be written.
	//  The threshold is based on the double epsilon and is applied to the
	//  absolute value, not a relative one, so IsZero() and IsLessEpsilon()
	//  are not suited to quantities of widely different magnitude.
	//--------------------------------------------------------------
	/**
	 * @brief Test whether a value is zero.
	 * @return true when |val| is at or below the smallest normal double.
	 */
	static inline bool IsZero(const double &val);
	/**
	 * @brief Test whether a value is negligible at machine precision.
	 * @return true when |val| is at or below the double epsilon.
	 */
	static inline bool IsLessEpsilon(const double &val);
	/**
	 * @brief Compare two doubles for equality instead of writing ==.
	 * @return true when |v1 - v2| is within one epsilon of either operand
	 *         (Knuth's test), so this one is a relative comparison.
	 */
	static inline bool Equal(const double &v1, const double &v2);

	//  Minimum and maximum (two-argument and three-argument forms)
	/** @brief Smaller / larger of two values, for int and for double. */
	static inline int    Min(const int    &v1, const int    &v2);
	static inline double Min(const double &v1, const double &v2);
	static inline int    Max(const int    &v1, const int    &v2);
	static inline double Max(const double &v1, const double &v2);

	/** @brief Smallest / largest of three values, for int and for double. */
	static inline int    Min(const int    &v1, const int    &v2, const int    &v3);
	static inline double Min(const double &v1, const double &v2, const double &v3);
	static inline int    Max(const int    &v1, const int    &v2, const int    &v3);
	static inline double Max(const double &v1, const double &v2, const double &v3);


	/**
	 * @brief Sign of a value.
	 * @return +1 when val is zero or positive, -1 when it is negative.
	 *
	 * Used to hand the flow direction to the upwind schemes.
	 */
	static inline int    Sign(const int    &val);
	static inline double Sign(const double &val);
	/**
	 * @brief Magnitude of one argument carrying the sign of the other,
	 *        the same as Fortran's SIGN intrinsic.
	 * @param value     the magnitude to return
	 * @param for_sign  the value whose sign is applied
	 * @return |value| when for_sign is zero or positive, -|value| otherwise.
	 */
	static inline int    Sign(const int    &value, const int    &for_sign);
	static inline double Sign(const double &value, const double &for_sign);

	/**
	 * @brief Smooth step going from 0 to 1 across a width eps centered on x0.
	 * @param x    where the step is evaluated
	 * @param x0   center of the transition (the value where it reads 0.5)
	 * @param eps  width of the transition
	 *
	 * Used wherever a discontinuity has to be avoided, such as the in-plane
	 * distribution of evaporation or the time variation of a source.
	 */
	static inline double Sigmoid(const double &x, const double &x0, const double &eps);

	/** @brief Pi, evaluated as 4*atan(1). */
	static inline double pi(){return 4.0*atan(1.0);}

	/**
	 * @brief Returns (a^2 - b^2)/(a + b), which algebraically equals a - b.
	 * @return The difference, or a + b itself when that denominator is zero.
	 *
	 * Note that Diff(-u, b) comes out as -u-b, so do not use this to combine
	 * signed quantities.
	 */
	static inline double Diff(double a, double b);
};

#include "MyMath_I.h"

#endif //_MY_MATH_INCLUDED_
