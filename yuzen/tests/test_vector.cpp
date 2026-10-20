/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_vector.cpp : test_vector.cpp   --   regression tests for MyMath/Vector.h
*
*  There are two points of emphasis:
*   (1) the basic operations agree with hand calculation
*   (2) the reductions (inner product / Sum / Norm) return bit-identical
*       values regardless of the thread count, which the deterministic
*       block reduction guarantees
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "doctest/doctest.h"
#include "Vector.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace {

Vector<double> ramp(int n, double a, double b)
{
	Vector<double> v(n);
	for(int i = 0; i < n; ++i) v[i] = a + b * static_cast<double>(i);
	return v;
}

} // namespace

TEST_CASE("Vector: basic operations")
{
	Vector<double> v(4);
	v[0] = 1.0; v[1] = -2.0; v[2] = 3.0; v[3] = -4.0;

	CHECK(v.size() == 4);
	CHECK(v.Sum()    == doctest::Approx(-2.0));
	CHECK(v.Max()    == doctest::Approx( 3.0));
	CHECK(v.Min()    == doctest::Approx(-4.0));
	CHECK(v.AbsMax() == doctest::Approx( 4.0));
	CHECK(v.Norm()   == doctest::Approx(std::sqrt(30.0)));

	Vector<double> w(4);
	w = 2.0;
	CHECK((v, w) == doctest::Approx(-4.0));		// inner product

	Vector<double> s = v + w;
	CHECK(s[0] == doctest::Approx(3.0));
	CHECK(s[3] == doctest::Approx(-2.0));

	v *= 2.0;
	CHECK(v[2] == doctest::Approx(6.0));
}

TEST_CASE("Vector: resize / zeroize")
{
	Vector<double> v(8);
	for(int i = 0; i < 8; ++i) CHECK(v[i] == 0.0);	// cleared to 0 by default

	v.resize(3);
	CHECK(v.size() == 3);
	for(int i = 0; i < 3; ++i) CHECK(v[i] == 0.0);
}

TEST_CASE("Vector: reductions of an empty vector")
{
	Vector<double> v(0);
	CHECK(v.Sum()  == doctest::Approx(0.0));
	CHECK(v.Norm() == doctest::Approx(0.0));
}

TEST_CASE("Vector: the reductions are bit-identical whatever the thread count")
{
	// Pick a size that is not a multiple of the block count (256), so that
	// the data straddles the block-partition boundaries
	const int N = 100003;
	Vector<double> a = ramp(N,  0.5,  1e-3);
	Vector<double> b = ramp(N, -1.25, 3e-4);

#if defined(_OPENMP)
	const int saved = omp_get_max_threads();
	double dot1 = 0, sum1 = 0, nrm1 = 0;
	omp_set_num_threads(1);
	dot1 = (a, b); sum1 = a.Sum(); nrm1 = a.Norm();

	for(int th = 2; th <= 8; th *= 2)
	{
		omp_set_num_threads(th);
		CAPTURE(th);
		// We require exact equality, not Approx. If this is loosened,
		// "the answer differs between 1 thread and N threads" comes back.
		CHECK((a, b)  == dot1);
		CHECK(a.Sum() == sum1);
		CHECK(a.Norm()== nrm1);
	}
	omp_set_num_threads(saved);
#else
	// In a serial build we only confirm that it agrees with the plain
	// sequential sum over a single partition
	double ref = 0;
	for(int i = 0; i < N; ++i) ref += a[i] * b[i];
	CHECK((a, b) == doctest::Approx(ref));
#endif
}

TEST_CASE("Vector: a small array agrees exactly with the plain sequential sum")
{
	// The contract is that a size at or below the block count (256) is not
	// split into blocks and is summed sequentially
	const int N = 200;
	Vector<double> a = ramp(N, 0.3, 7e-3);
	double ref = 0;
	for(int i = 0; i < N; ++i) ref += a[i];
	CHECK(a.Sum() == ref);
}
