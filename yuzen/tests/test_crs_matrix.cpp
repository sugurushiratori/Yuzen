/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  test_crs_matrix.cpp : tests for the CRS matrix
*
*  test_crs_matrix.cpp  --  regression tests for numerics/CRS_Matrix.h
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "doctest/doctest.h"
#include "CRS_Matrix.h"

TEST_CASE("CRS_Matrix: element access and the non-zero count")
{
	CRS_Matrix m;
	m.Build(4);

	CHECK(m.RowSize() == 4);
	// Right after Build, each row holds only its diagonal element
	CHECK(m.NonZero() == 4);
	for(int i = 0; i < 4; ++i) CHECK(m.Diag(i) == doctest::Approx(0.0));

	m.Diag(0) = 2.0;  m(0,1) = -1.0;
	m(1,0) = -1.0;    m.Diag(1) = 2.0;  m(1,2) = -1.0;
	m(2,1) = -1.0;    m.Diag(2) = 2.0;  m(2,3) = -1.0;
	m(3,2) = -1.0;    m.Diag(3) = 2.0;

	CHECK(m.NonZero() == 10);
	CHECK(m.NonZero(1) == 3);
	CHECK(m(1,2) == doctest::Approx(-1.0));
	// Reading a position that is not stored returns 0 (the const version)
	const CRS_Matrix &cm = m;
	CHECK(cm(0,3) == doctest::Approx(0.0));

	std::pair<int,int> lu = m.NonZ_LU();
	CHECK(lu.first  == 3);	// number of non-zeros in the lower triangle
	CHECK(lu.second == 3);	// number of non-zeros in the upper triangle
}

TEST_CASE("CRS_Matrix: MatVec and operator* agree with the dense result")
{
	const int n = 5;
	CRS_Matrix m;
	m.Build(n);

	// Tridiagonal [-1, 2, -1]
	for(int i = 0; i < n; ++i)
	{
		m.Diag(i) = 2.0;
		if(i > 0)   m(i, i-1) = -1.0;
		if(i < n-1) m(i, i+1) = -1.0;
	}

	Vector<double> x(n);
	for(int i = 0; i < n; ++i) x[i] = static_cast<double>(i + 1);

	// Reference computation with a dense matrix
	Vector<double> ref(n);
	for(int i = 0; i < n; ++i)
	{
		double s = 2.0 * x[i];
		if(i > 0)   s -= x[i-1];
		if(i < n-1) s -= x[i+1];
		ref[i] = s;
	}

	Vector<double> y(n);
	m.MatVec(y, x);
	for(int i = 0; i < n; ++i) { CAPTURE(i); CHECK(y[i] == doctest::Approx(ref[i])); }

	Vector<double> z = m * x;
	for(int i = 0; i < n; ++i) { CAPTURE(i); CHECK(z[i] == doctest::Approx(ref[i])); }
}

TEST_CASE("CRS_Matrix: ClearData zeroes the values and keeps the sparsity structure")
{
	CRS_Matrix m;
	m.Build(3);
	m.Diag(0) = 1.0; m(0,2) = 5.0;
	const int nnz = m.NonZero();

	m.ClearData();
	CHECK(m.NonZero() == nnz);		// the sparsity structure is kept
	CHECK(m.Diag(0) == doctest::Approx(0.0));
	CHECK(m(0,2)    == doctest::Approx(0.0));
}

TEST_CASE("CRS_Matrix: SetAsBC marks a row and IsBC reads it back")
{
	CRS_Matrix m;
	m.Build(3);
	CHECK(m.IsBC(1) == false);
	m.SetAsBC(1);
	CHECK(m.IsBC(1) == true);
	CHECK(m.IsBC(0) == false);
}
