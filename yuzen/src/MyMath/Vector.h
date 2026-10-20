/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Vector.h : one-dimensional array for numerical work
*
*  The storage is left to std::vector, and the element-by-element operations
*  do not open a parallel region of their own. Vector arithmetic is memory
*  bound, so opening a parallel region on the short inner vectors costs more
*  in fork/join than it saves (that would be more than 20 of them per BiCGStab
*  iteration). The right way to parallelize this is to wrap the whole
*  iteration loop in a single parallel region.
*
*  Only the inner product, Sum and Norm are different: they are written as a
*  deterministic fixed-block reduction so that the result is bit-identical
*  whatever the thread count (see the note below).
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _MY_VECTOR_INCLUDED_
#define _MY_VECTOR_INCLUDED_

#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>

#if defined (OMP) || defined(_OPENMP)
#include <omp.h>
#endif

#include "MyMath.h"
#include "FirstTouch.h"

//---------------------------------------------------------------------------
//  Deterministic reduction
//
//  With #pragma omp parallel for reduction(+:x), both the way the partial
//  sums are split up and the order in which they are added depend on the
//  thread count. Floating-point addition is not associative, so changing the
//  thread count changes the low-order bits of an inner product or a norm.
//  That difference is then amplified through the BiCGStab iteration count and
//  through the decisions made by the time-step control, and it surfaces as
//  "one thread and two threads give different results".
//
//  Instead, the array is split here into a fixed number of blocks that does
//  not depend on the thread count. Within a block the sum is sequential, and
//  the blocks themselves are added in a fixed order as well. The order of the
//  additions is therefore always the same, so the value comes out identical
//  bit for bit however many threads are used.
//---------------------------------------------------------------------------
namespace MyMathDetail
{
	const int kReduceBlocks   = 256;	// upper limit on the number of blocks (the parallelism is not raised beyond this)
	const int kMinBlockLength = 512;	// smallest number of elements in one block

	/**
	 * @brief Number of blocks the reduction over n elements is split into.
	 *
	 * The count follows from the array length alone and does not depend on the
	 * thread count. That is what makes the summation order, and with it the
	 * result, independent of how many threads run it. Cutting a short array
	 * into many blocks would cost more in block bookkeeping than it saves, so
	 * a block is never allowed to hold fewer than kMinBlockLength elements.
	 */
	inline int ReduceBlockCount(int n)
	{
		if(n <= kMinBlockLength) return 1;		// the same order as a sequential sum
		const int nb = n / kMinBlockLength;
		return (nb > kReduceBlocks) ? kReduceBlocks : nb;
	}
	/**
	 * @brief First index of block b, when n elements are split into nb blocks.
	 * @return The index n*b/nb, computed in 64 bits so that the product does
	 *         not overflow for a long array.
	 *
	 * The end of block b is the beginning of block b+1, so the blocks cover
	 * the array exactly once with no gap and no overlap.
	 */
	inline int ReduceBlockBegin(int n, int b, int nb) { return (int)((long long)n * b / nb); }
}


template <class T> class Vector
{
public:
	//---------------------------------
	//    Public methods
	//---------------------------------
	/**
	 * @brief Construct an empty vector.
	 * @param zeroize  ignored: elements are always zero initialized.
	 */
	Vector(bool zeroize = true)           : _zeroize(zeroize) {}
	/**
	 * @brief Construct a vector of the given length.
	 * @param size     number of elements; the storage is zero filled
	 * @param zeroize  ignored: elements are always zero initialized.
	 */
	Vector(int size, bool zeroize = true) : _zeroize(zeroize) { resize(size); }
	/**
	 * @brief Copy constructor.
	 * @param v  vector whose contents and settings are copied
	 */
	Vector(const Vector<T> &v)            : _body(v._body), _zeroize(v._zeroize) {}
	~Vector() {}

	/** @brief Number of elements. */
	int  size()   const {return (int)_body.size();}

	/**
	 * @brief Set the length and zero every element.
	 *
	 * The contents are always zero initialized, resizeFast() included. Any
	 * path that could read an uninitialized value would make the results
	 * non-deterministic, so neither of these cuts that corner.
	 */
	void resize    (int size) { assert(size >= 0); _body.assign((size_t)size, T()); }
	void resizeFast(int size) { resize(size); }
	/**
	 * @brief Change the length while keeping the existing contents.
	 * @param size  new length; elements added when growing are zero.
	 */
	void truncate  (int size) { assert(size >= 0); _body.resize((size_t)size, T()); }

	/**
	 * @brief Pointer to the contiguous storage, for handing the data to
	 *        CRS_Matrix and to the linear solvers.
	 * @return Address of the first element, or a null pointer when empty.
	 */
	      T* data()       { return _body.empty() ? 0 : &_body[0]; }
	const T* data() const { return _body.empty() ? 0 : &_body[0]; }

	//---------------------------------
	//    operators
	//---------------------------------
	/** @brief Element-by-element comparison of the contents. */
	bool operator==(const Vector<T>& v) const { return _body == v._body; }
	bool operator!=(const Vector<T>& v) const { return !operator==(v); }

	/** @brief Copy the contents of another vector, taking its length too. */
	void operator= (const Vector<T>& v) { _body = v._body; }
	/**
	 * @brief Fill every existing element with a scalar.
	 *
	 * The length is not changed; use resize() for that.
	 */
	void operator= (const T x)          { std::fill(_body.begin(), _body.end(), x); }

	//---------------------------------
	//    accessors
	//---------------------------------
	//  Element access. The index is checked by assert in debug builds only.
	T& operator[](int i)
	{
		assert((i >=0) && (i < size()));
		return _body[(size_t)i];
	}
	const T& operator[](int i) const
	{
		assert((i >=0) && (i < size()));
		return _body[(size_t)i];
	}
	T&       operator()(int i)       { return operator[](i); }
	const T& operator()(int i) const { return operator[](i); }

	//---------------------------------
	//    arithmetic
	//---------------------------------
	/**
	 * @brief Element-by-element arithmetic, returning a new vector.
	 *
	 * Note that operator* and operator/ on two vectors multiply and divide
	 * element by element; the inner product is operator, instead.
	 */
	Vector<T> operator+ (const Vector<T>& v) const { Vector<T> t(*this); t += v; return t; }
	Vector<T> operator- (const Vector<T>& v) const { Vector<T> t(*this); t -= v; return t; }
	Vector<T> operator* (const Vector<T>& v) const { Vector<T> t(*this); t *= v; return t; }
	Vector<T> operator/ (const Vector<T>& v) const { Vector<T> t(*this); t /= v; return t; }
	Vector<T> operator* (const        T   x) const { Vector<T> t(*this); t *= x; return t; }
	Vector<T> operator/ (const        T   x) const { Vector<T> t(*this); t /= x; return t; }

	/**
	 * @brief Element-by-element arithmetic in place.
	 *
	 * Both operands must have the same length (checked by assert).
	 */
	void operator+=(const Vector<T>& v)
	{
		assert(size() == v.size());
		const int n = size();
		for(int i = 0; i < n; ++i) _body[(size_t)i] += v[i];
	}
	void operator-=(const Vector<T>& v)
	{
		assert(size() == v.size());
		const int n = size();
		for(int i = 0; i < n; ++i) _body[(size_t)i] -= v[i];
	}
	void operator*=(const Vector<T>& v)
	{
		assert(size() == v.size());
		const int n = size();
		for(int i = 0; i < n; ++i) _body[(size_t)i] *= v[i];
	}
	void operator/=(const Vector<T>& v)
	{
		assert(size() == v.size());
		const int n = size();
		for(int i = 0; i < n; ++i) _body[(size_t)i] /= v[i];
	}
	/** @brief Scale every element by a scalar, in place. */
	void operator*=(const T x) { const int n = size(); for(int i = 0; i < n; ++i) _body[(size_t)i] *= x; }
	void operator/=(const T x) { const int n = size(); for(int i = 0; i < n; ++i) _body[(size_t)i] /= x; }

	/**
	 * @brief Inner product of this vector with v.
	 * @return Sum of the element-by-element products.
	 *
	 * The summation uses the deterministic fixed-block reduction described
	 * below, so the result is bit-identical whatever the thread count.
	 */
	T operator, (const Vector<T>& v) const
	{
		assert(size() == v.size());
		const int n = size();
		const int nb = MyMathDetail::ReduceBlockCount(n);
		T partial[MyMathDetail::kReduceBlocks];
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for schedule(static) if(nb > 1 && yuzen::opt::Par(n))
#endif
		for(int b = 0; b < nb; ++b)
		{
			const int st = MyMathDetail::ReduceBlockBegin(n, b,   nb);
			const int ed = MyMathDetail::ReduceBlockBegin(n, b+1, nb);
			T s = 0;
			for(int i = st; i < ed; ++i) s += _body[(size_t)i] * v[i];
			partial[b] = s;
		}
		T ip = 0;
		for(int b = 0; b < nb; ++b) ip += partial[b];
		return ip;
	}

	//--------------------------------------------------------------
	//  Reductions (Sum / Norm / inner product)
	//
	//  *** These are bit-identical whatever the thread count ***
	//  The array is split into a **fixed number** of blocks, a partial sum is
	//  taken over each, and the partial sums are finally added in block order.
	//  The split is decided by the number of elements only, never by the
	//  thread count, so the additions happen in the same order however many
	//  threads run, and the answer agrees bit for bit.
	//  **Do not replace this with #pragma omp parallel for reduction(+:...)**:
	//  that makes the order depend on the thread count again.
	//  tests/test_vector.cpp pins this behavior down.
	//--------------------------------------------------------------
	/**
	 * @brief Sum of all elements.
	 * @return The total, bit-identical whatever the OpenMP thread count.
	 */
	T Sum() const
	{
		const int n = size();
		const int nb = MyMathDetail::ReduceBlockCount(n);
		T partial[MyMathDetail::kReduceBlocks];
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for schedule(static) if(nb > 1 && yuzen::opt::Par(n))
#endif
		for(int b = 0; b < nb; ++b)
		{
			const int st = MyMathDetail::ReduceBlockBegin(n, b,   nb);
			const int ed = MyMathDetail::ReduceBlockBegin(n, b+1, nb);
			T s = 0;
			for(int i = st; i < ed; ++i) s += _body[(size_t)i];
			partial[b] = s;
		}
		T sum = 0;
		for(int b = 0; b < nb; ++b) sum += partial[b];
		return sum;
	}

	/**
	 * @brief Euclidean (2-) norm.
	 * @return sqrt of the sum of squares, summed with the same deterministic
	 *         reduction as Sum(), so it too is bit-identical across threads.
	 */
	T Norm() const
	{
		const int n = size();
		const int nb = MyMathDetail::ReduceBlockCount(n);
		T partial[MyMathDetail::kReduceBlocks];
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for schedule(static) if(nb > 1 && yuzen::opt::Par(n))
#endif
		for(int b = 0; b < nb; ++b)
		{
			const int st = MyMathDetail::ReduceBlockBegin(n, b,   nb);
			const int ed = MyMathDetail::ReduceBlockBegin(n, b+1, nb);
			T s = 0;
			for(int i = st; i < ed; ++i) s += _body[(size_t)i] * _body[(size_t)i];
			partial[b] = s;
		}
		T norm = 0;
		for(int b = 0; b < nb; ++b) norm += partial[b];
		return std::sqrt(norm);
	}

	//  max and min are exactly associative and commutative, so sequential is enough
	/**
	 * @brief Largest absolute value among the elements.
	 *
	 * Used to measure the size of the Newton correction.
	 * Unlike the sums above this runs sequentially, which is enough because
	 * max and min are exactly associative and commutative.
	 */
	T AbsMax() const
	{
		T m = 0;
		const int n = size();
		for(int i = 0; i < n; ++i) m = Math::Max(Math::abs(_body[(size_t)i]), Math::abs(m));
		return m;
	}
	/**
	 * @brief Largest element.
	 * @return The maximum, or -1e120 if the vector is empty (the search
	 *         starts from that sentinel rather than from -infinity).
	 */
	T Max() const
	{
		T m = T(-1e120);
		const int n = size();
		for(int i = 0; i < n; ++i) m = Math::Max(_body[(size_t)i], m);
		return m;
	}
	/**
	 * @brief Smallest element.
	 * @return The minimum, or 1e120 if the vector is empty.
	 */
	T Min() const
	{
		T m = T(1e120);
		const int n = size();
		for(int i = 0; i < n; ++i) m = Math::Min(_body[(size_t)i], m);
		return m;
	}

protected:
	std::vector<T, yuzen::util::FirstTouchAllocator<T> > _body;
	bool           _zeroize;	// recorded only; the storage is always zero initialized
};

#endif //_MY_VECTOR_INCLUDED_
