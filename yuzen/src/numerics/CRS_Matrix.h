/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CRS_Matrix.h : sparse matrix in compressed row storage (CRS)
*
*  The storage is three flat arrays. Nothing is allocated per row.
*
*  This is genuine CRS (the three arrays row_ptr / col_ind / values), which
*  buys the following.
*   - Copying is the equivalent of a memcpy over contiguous memory.
*   - The SpMV accesses memory contiguously, so vectorization and prefetching
*     work.
*   - The matrix can be handed to Eigen through Eigen::Map with no copy at all
*     (LinearSolver.cpp).
*
*  There are two states.
*   Building (_building = true) : one variable-length array per row. The rows
*                                 are independent of each other, so inserting
*                                 from several OMP threads at once is safe
*                                 (this is what the first assembly of the
*                                 physical terms does).
*   Finalized (_building = false): a single flat CRS. No insertion happens.
*  Any read-side call (MatVec, NonZero(), Values() and so on) finalizes the
*  matrix automatically.
*
*  ** The storage order is fixed: insertion order, with the diagonal first in
*  each row. **
*  ScalingDiag and Residual sum in storage order, so changing the order changes
*  the floating-point rounding and the results no longer match bit for bit.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _CRS_MATRIX_INCLUDED_
#define _CRS_MATRIX_INCLUDED_


#include <vector>
#include <cassert>
#include <string>
#include <utility>

#include "Vector.h"
#include "FirstTouch.h"


//--------------------------------------------------------------------------
//  Sparse matrix in compressed row storage (CRS).
//
//  It has two phases.
//    Building phase  : touching an entry that does not exist yet through
//                      operator() inserts it. Within a row the entries end up
//                      in the order they were inserted.
//    Finalized phase : the entries are gathered into flat CRS arrays. The
//                      matrix moves here when Finalize() is called, and the
//                      structure cannot be changed afterwards.
//
//  ** The storage order affects the computed results. **
//  Within a row the entries are in **insertion order**, and the first entry of
//  every row is always the diagonal. This is not a cosmetic detail: it fixes
//  the **order in which the floating-point additions of the matrix-vector
//  product are performed**, so reordering the entries changes the result bit
//  for bit. The bit-for-bit agreement the regression tests rely on depends on
//  this order, so **do not sort the entries into ascending column order.**
//  For libraries that require ascending columns (Eigen and the like) a sorted
//  mirror is built alongside the matrix itself and handed over instead
//  (see Sorted* below).
//--------------------------------------------------------------------------
class CRS_Matrix
{

public:
	/** @brief Construct an empty matrix, already in the building phase. */
	CRS_Matrix() : _n_row(0), _zero(0.0), _building(true) {}

	/**
	 * @brief Fix the number of rows and enter the building phase.
	 * @param n_row  number of rows, which is the number of unknowns
	 */
	void Build(int n_row);

	/** @brief Number of rows, which is the number of unknowns. */
	int RowSize() const {return _n_row;}

	//------------------------------------------------------------------
	//  Number of non-zero entries
	//------------------------------------------------------------------
	/**
	 * @brief Total number of non-zero entries.
	 * @return The count over the whole matrix.
	 *
	 * Calling this moves the matrix into the finalized phase.
	 */
	int NonZero() const { Finalize(); return (int)_colInd.size(); }
	/**
	 * @brief Number of non-zero entries in one row.
	 * @param i_row  row number
	 */
	int NonZero(const int i_row) const { return RowLen(i_row); }

	/**
	 * @brief Number of non-zero entries in the strict lower and upper
	 *        triangles.
	 * @return (lower, upper) counts, used to size the ILU working arrays.
	 */
	std::pair<int, int> NonZ_LU() const;

	/**
	 * @brief Copy only the values from another matrix, leaving the sparsity
	 *        pattern of this one untouched.
	 */
	void CopyDataOnly(const CRS_Matrix& m);
	/**
	 * @brief Set every value to zero, keeping the structure.
	 *
	 * The structure survives and the first entry of each row is guaranteed
	 * to be the diagonal afterwards. Call this before reassembling the
	 * matrix at each step.
	 */
	void ZeroClear();
	/** @brief Discard the structure as well and return to the building phase. */
	void ClearData();

	/**
	 * @brief Mark a row as a boundary-condition row.
	 *
	 * Marked rows are excluded from the time integration and the physical
	 * terms are not accumulated into them; each kernel skips them via IsBC().
	 */
	void SetAsBC(const int i){_flagBC[i] = true;}
	/** @brief Whether a row has been marked as a boundary-condition row. */
	bool    IsBC(const int i) const {return _flagBC[i] != 0;}

	//------------------------------------------------------------------
	//  Diagonal entry. ZeroClear() is what guarantees that the first entry of
	//  each row (storage position 0) is the diagonal.
	//------------------------------------------------------------------
	/** @brief Diagonal entry of a row. */
	      double& Diag(int i)       {assert(RowLen(i) > 0); return ValRef(i, 0);}
	const double& Diag(int i) const {assert(RowLen(i) > 0); return ValRef(i, 0);}

	/**
	 * @brief Checked access by (row, column) number.
	 * @param i_row  row number
	 * @param j_col  column number
	 * @return Reference to the entry; the non-const version inserts it when
	 *         it does not exist yet, which is only possible in the building
	 *         phase. The const version returns a reference to a shared zero.
	 *
	 * This is what the physical terms use when they accumulate into the
	 * matrix.
	 */
	double& operator()(const int i_row, const int j_col)
	{
		const int k = FindCol(i_row, j_col);
		if(k >= 0) return ValRef(i_row, k);
		return Insert(i_row, j_col);
	}
	const double& operator()(const int i_row, const int j_col) const
	{
		const int k = FindCol(i_row, j_col);
		return (k < 0) ? _zero : ValRef(i_row, k);
	}

	/**
	 * @brief Unchecked access by storage position within a row.
	 * @param i_row    row number
	 * @param col_ind  position within the row, not the column number.
	 *                 Position 0 is the diagonal.
	 * @return Reference to the entry.
	 *
	 * No range check is done, so use this only when walking a structure that
	 * is already fixed.
	 */
	      double& Element(const int i_row, const int col_ind)       {return ValRef(i_row, col_ind);}
	const double& Element(const int i_row, const int col_ind) const {return ValRef(i_row, col_ind);}

	/**
	 * @brief Column number of an entry given its position within the row.
	 * @param i_row    the row
	 * @param col_ind  position within the row (0 is the diagonal)
	 */
	int Jcol(const int i_row, const int col_ind) const {return ColRef(i_row, col_ind);}

	/**
	 * @brief Position within the row of a given column number.
	 * @return The storage position, or -1 when that column is not stored.
	 */
	int ColInd(const int i_row, const int j_col) const {return FindCol(i_row, j_col);}

	/**
	 * @brief Matrix-vector product, returning a newly allocated result.
	 *
	 * Note that in the main matrix, which carries a fourth-order difference,
	 * terms of order 1e15 cancel against each other here to produce the
	 * answer. That cancellation is why the Newton correction stops
	 * decreasing below a certain level; it has nothing to do with the
	 * accuracy of the linear solver (docs/PHASE7.md 8.5).
	 */
	Vector<double> operator*(const Vector<double>& v) const;
	/**
	 * @brief Matrix-vector product y = A x into caller-owned storage.
	 * @param y  [out] result, supplied by the caller so that no allocation
	 *           happens here
	 * @param x  vector to multiply
	 *
	 * The same cancellation noted for operator* applies.
	 */
	void MatVec(Vector<double>& y, const Vector<double>& x) const;
	/**
	 * @brief Make x satisfy the boundary rows (IsBC) of A x = rhs exactly:
	 *        each set of boundary rows that refer to each other's unknowns
	 *        solved for them, with the other unknowns as they are: a row on
	 *        its own as x_i = (rhs_i - sum_{j != i} A_ij x_j) / A_ii, a
	 *        coupled set (the two rows of an Outlet or a FreshInlet, some
	 *        corners) by Gaussian elimination. The other rows are left
	 *        alone.
	 *
	 * An iterative solve leaves the boundary rows satisfied only to its
	 * tolerance, and a mirror ghost (Symmetry) a little off its mirror;
	 * the capillary flux of an edge face, a third difference over the
	 * mirror, is then not 0, and on a stiff film it is a leak of volume.
	 * Satisfied this way, the mirror and the periodic copy are exact.
	 */
	void SatisfyBCRows(const Vector<double>& rhs, Vector<double>& x) const;

	//------------------------------------------------------------------
	//  Direct access to the flat CRS arrays.
	//  They are in the form that can be passed straight to
	//  Eigen::Map<Eigen::SparseMatrix<double, Eigen::RowMajor> >.
	//  Calling any of them moves the matrix into the finalized phase.
	//------------------------------------------------------------------
	/**
	 * @brief Pointers to the flat CRS arrays (row_ptr, col_ind, values).
	 * @return Address of the first element, or a null pointer when empty.
	 *
	 * Calling any of these moves the matrix into the finalized phase. The
	 * column numbers are in storage order, so a library that requires them
	 * to be ascending needs SortedColIndData() instead.
	 */
	const int*    RowPtrData() const { Finalize(); return _rowPtr.empty() ? 0 : &_rowPtr[0]; }
	const int*    ColIndData() const { Finalize(); return _colInd.empty() ? 0 : &_colInd[0]; }
	const double* ValuesData() const { Finalize(); return _values.empty() ? 0 : &_values[0]; }
	      double* ValuesData()       { Finalize(); return _values.empty() ? 0 : &_values[0]; }

	//------------------------------------------------------------------
	//  A sorted mirror, for libraries (Eigen and the like) that require the
	//  column numbers within a row to be ascending. The permutation is built
	//  once, when the structure is finalized; after that each call only
	//  gathers the values, which is O(nnz).
	//  The storage order of the matrix itself (insertion order, diagonal
	//  first) is left alone. Changing that order would change the rounding of
	//  the sums in ScalingDiag and Residual.
	//------------------------------------------------------------------
	/**
	 * @brief Column numbers of the sorted mirror, ascending within each row.
	 * @return Array of nnz column numbers, valid until the structure changes.
	 */
	const int*    SortedColIndData() const;
	/**
	 * @brief Values of the sorted mirror, in the same order as
	 *        SortedColIndData().
	 * @return Array of nnz values, gathered fresh on every call so that the
	 *         latest values are picked up.
	 */
	const double* GatherSortedValues() const;	// gathers the latest values on every call
	/** @brief Where in ValuesData() the k-th entry of the column-sorted order sits. */
	const int*    SortedSrcData() const;

	//For Debug
	/**
	 * @brief Write the matrix to a text file, for debugging.
	 * @param wfile  path of the file to write
	 */
	void Write(const std::string &wfile) const;

protected:
	int _n_row;
	double _zero;

	// Building phase (rows are independent, so inserting in parallel with OMP is safe)
	mutable bool _building;
	std::vector<std::vector<int> >    _bInd;
	std::vector<std::vector<double> > _bVal;

	// Finalized phase (flat CRS)
	mutable std::vector<int, yuzen::util::FirstTouchAllocator<int> >       _rowPtr;	// size = _n_row + 1
	mutable std::vector<int, yuzen::util::FirstTouchAllocator<int> >       _colInd;	// size = nnz
	mutable std::vector<double, yuzen::util::FirstTouchAllocator<double> > _values;	// size = nnz

	// char, not std::vector<bool>, whose bit-packed references are not thread safe
	std::vector<char> _flagBC;

	// The sorted mirror (for talking to Eigen)
	mutable std::vector<int, yuzen::util::FirstTouchAllocator<int> >       _sortedCol;	// column numbers in ascending order
	mutable std::vector<int, yuzen::util::FirstTouchAllocator<int> >       _sortedSrc;	// where in _values the value of _sortedCol[k] sits
	mutable std::vector<double, yuzen::util::FirstTouchAllocator<double> > _sortedVal;
	void BuildSortPermutation() const;

	//------------------------------------------------------------------
	//  Internal accessors. They mean the same thing whether the matrix is
	//  still being built or already finalized.
	//------------------------------------------------------------------
	int RowLen(int i) const
	{
		assert(0 <= i && i < _n_row);
		return _building ? (int)_bInd[i].size() : (_rowPtr[i+1] - _rowPtr[i]);
	}
	double& ValRef(int i, int k)
	{
		assert(0 <= k && k < RowLen(i));
		return _building ? _bVal[i][k] : _values[_rowPtr[i] + k];
	}
	const double& ValRef(int i, int k) const
	{
		assert(0 <= k && k < RowLen(i));
		return _building ? _bVal[i][k] : _values[_rowPtr[i] + k];
	}
	int ColRef(int i, int k) const
	{
		assert(0 <= k && k < RowLen(i));
		return _building ? _bInd[i][k] : _colInd[_rowPtr[i] + k];
	}
	// Linear search in storage order. Returns -1 when the column is not found.
	int FindCol(int i_row, int j_col) const
	{
		const int n = RowLen(i_row);
		if(_building)
		{
			const std::vector<int> &ind = _bInd[i_row];
			for(int k = 0; k < n; ++k) if(ind[k] == j_col) return k;
		}
		else
		{
			const int st = _rowPtr[i_row];
			for(int k = 0; k < n; ++k) if(_colInd[st + k] == j_col) return k;
		}
		return -1;
	}
	double& Insert(int i_row, int j_col);

	void Finalize() const;		// per-row arrays -> flat CRS
	//  Flat CRS -> per-row arrays, so that insertion becomes possible again.
	//  **Do not call this from inside a parallel region.** It reallocates the
	//  whole per-row structure, so it is not the insertion into a single row
	//  that the building phase makes thread safe.
	void Open();
};


#endif //_CRS_MATRIX_INCLUDED_
