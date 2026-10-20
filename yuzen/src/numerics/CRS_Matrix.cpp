/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CRS_Matrix.cpp : implementation of the flat CRS matrix
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "Opt.h"
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <fstream>
using namespace std;

#include "YuzenUtil.h"
using namespace yuzen::util;

#include "CRS_Matrix.h"


void CRS_Matrix::Build(int n_row)
{
	assert(n_row > 0);
	_n_row = n_row;
	_zero  = 0.0;

	_building = true;
	_bInd.assign(_n_row, std::vector<int>());
	_bVal.assign(_n_row, std::vector<double>());
	_rowPtr.clear(); _colInd.clear(); _values.clear();

	ZeroClear();

	_flagBC.assign(_n_row, 0);
}


//------------------------------------------------------------------
//  Reduce the structure back to the diagonal alone.
//  This is where the assumption Diag() relies on is established: that the
//  first entry of a row (storage position 0) is the diagonal.
//------------------------------------------------------------------
void CRS_Matrix::ZeroClear()
{
	Open();
	for(int i = 0; i < _n_row; ++i)
	{
		_bInd[i].clear();	_bInd[i].push_back(i);
		_bVal[i].clear();	_bVal[i].push_back(0.0);
	}
}


//------------------------------------------------------------------
//  Set the values to zero while keeping the structure
//------------------------------------------------------------------
void CRS_Matrix::ClearData()
{
	if(_building)
	{
		for(int i = 0; i < _n_row; ++i)
			std::fill(_bVal[i].begin(), _bVal[i].end(), 0.0);
	}
	else
	{
		std::fill(_values.begin(), _values.end(), 0.0);
	}
}


void CRS_Matrix::CopyDataOnly(const CRS_Matrix& m)
{
	assert(m.RowSize() == _n_row);
	Finalize();
	m.Finalize();
	assert(m._colInd.size() == _colInd.size());
	_values = m._values;
}


//------------------------------------------------------------------
//  Append a new column at the end of a row.
//  In the building phase each row has its own independent array, so threads
//  that are working on different rows may call this at the same time safely
//  (this is what the first assembly of the physical terms does).
//------------------------------------------------------------------
double& CRS_Matrix::Insert(int i_row, int j_col)
{
	// Inserting after the matrix has been finalized is rare; it happens when
	// the boundary conditions are applied for the first time, for instance.
	if(!_building) Open();
	_bInd[i_row].push_back(j_col);
	_bVal[i_row].push_back(0.0);
	return _bVal[i_row].back();
}


//------------------------------------------------------------------
//  Per-row arrays -> flat CRS
//  The storage order is carried over unchanged, so that the rounding of the
//  sums does not change.
//------------------------------------------------------------------
void CRS_Matrix::Finalize() const
{
	if(!_building) return;

	CRS_Matrix *self = const_cast<CRS_Matrix*>(this);

	self->_rowPtr.resize(_n_row + 1);
	int nnz = 0;
	for(int i = 0; i < _n_row; ++i)
	{
		self->_rowPtr[i] = nnz;
		nnz += (int)_bInd[i].size();
	}
	self->_rowPtr[_n_row] = nnz;

	self->_colInd.resize(nnz);
	self->_values.resize(nnz);
	for(int i = 0; i < _n_row; ++i)
	{
		const int st = _rowPtr[i];
		const int n  = (int)_bInd[i].size();
		for(int k = 0; k < n; ++k)
		{
			self->_colInd[st + k] = _bInd[i][k];
			self->_values[st + k] = _bVal[i][k];
		}
	}

	// Release the storage that was only needed while building
	std::vector<std::vector<int> >().swap(self->_bInd);
	std::vector<std::vector<double> >().swap(self->_bVal);

	self->_building = false;
}


//------------------------------------------------------------------
//  Flat CRS -> per-row arrays (back to a state that accepts insertion).
//  **Do not call this from inside a parallel region.** It reallocates the
//  whole per-row structure, which is not what the building phase makes
//  thread safe; only insertion into separate rows is.
//------------------------------------------------------------------
void CRS_Matrix::Open()
{
	if(_building) return;

	_bInd.assign(_n_row, std::vector<int>());
	_bVal.assign(_n_row, std::vector<double>());
	for(int i = 0; i < _n_row; ++i)
	{
		const int st = _rowPtr[i], ed = _rowPtr[i+1];
		_bInd[i].assign(_colInd.begin() + st, _colInd.begin() + ed);
		_bVal[i].assign(_values.begin() + st, _values.begin() + ed);
	}
	decltype(_rowPtr)().swap(_rowPtr);
	decltype(_colInd)().swap(_colInd);
	decltype(_values)().swap(_values);

	_building = true;
}


//------------------------------------------------------------------
//  Number of off-diagonal non-zero entries in the lower and upper triangles
//------------------------------------------------------------------
std::pair<int,int> CRS_Matrix::NonZ_LU() const
{
	Finalize();
	int nzL = 0, nzU = 0;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for reduction(+:nzL,nzU)
#endif
	for(int i_row = 0; i_row < _n_row; i_row++)
	{
		const int st = _rowPtr[i_row], ed = _rowPtr[i_row+1];
		for(int k = st; k < ed; ++k)
		{
			const int j_col = _colInd[k];
			if(0 <= j_col && j_col < _n_row)
			{
				if     (j_col < i_row) nzL++;
				else if(j_col > i_row) nzU++;
			}
		}
	}
	return std::pair<int,int>(nzL, nzU);
}


//------------------------------------------------------------------
//  Satisfy the boundary rows exactly (see the declaration)
//------------------------------------------------------------------
void CRS_Matrix::SatisfyBCRows(const Vector<double>& rhs, Vector<double>& x) const
{
	Finalize();

	//  The boundary rows, grouped into the sets that refer to each other's
	//  unknowns (union-find over the entries of a boundary row that sit in
	//  another boundary row's column). A Symmetry or Periodic ghost is a set
	//  of its own; the two rows of an Outlet, a FreshInlet or a 1-D pinned
	//  concentration are coupled, and so are some corners. Taken one row
	//  after the other (Gauss-Seidel, as this did) such a pair is not
	//  solved: the error of the pair shrinks by a factor 3 per sweep
	//  (Outlet, FreshInlet), or grows by 27 (Neumann with the second
	//  derivative, 1-D Pinning of the concentration).
	std::vector<int> loc(_n_row, -1), rows, parent;
	for(int i = 0; i < _n_row; i++)
		if(IsBC(i)) { loc[i] = (int)rows.size(); rows.push_back(i); parent.push_back((int)parent.size()); }
	if(rows.empty()) return;
	auto find = [&parent](int a) { while(parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; } return a; };
	for(size_t r = 0; r < rows.size(); ++r)
		for(int k = _rowPtr[rows[r]] + 1; k < _rowPtr[rows[r]+1]; ++k)
		{
			const int j = loc[_colInd[k]];
			if(j < 0 || _values[k] == 0.0) continue;
			const int a = find((int)r), b = find(j);
			if(a != b) parent[std::max(a, b)] = std::min(a, b);
		}
	std::vector<std::vector<int> > sets(rows.size());
	for(size_t r = 0; r < rows.size(); ++r) sets[find((int)r)].push_back((int)r);

	std::vector<double> M, b;
	for(size_t c = 0; c < sets.size(); ++c)
	{
		const std::vector<int> &set = sets[c];
		const int m = (int)set.size();
		if(m == 0) continue;
		//  The right-hand side of each row of the set, less its terms in the
		//  unknowns outside the set (the interior, other sets). Those terms
		//  are summed in order of decreasing |A_ij| (equal ones in storage
		//  order), so that a pair that cancels in a relation such as
		//  x0 - 3 x1 + 3 x2 - x3 = 0 cancels first, whichever end of the row
		//  the diagonal sits at. A row whose diagonal is zero is left alone.
		M.assign((size_t)m * m, 0.0);
		b.assign(m, 0.0);
		bool solvable = true;
		for(int p = 0; p < m; ++p)
		{
			const int i = rows[set[p]];
			const int st = _rowPtr[i], ed = _rowPtr[i+1];
			if(ed <= st || _values[st] == 0.0) { solvable = false; break; }
			double s = 0;
			std::vector<char> used(ed - st, 0);
			for(int k = st; k < ed; ++k)
			{
				const int j = loc[_colInd[k]];
				if(j >= 0 && find(j) == (int)c)
				{
					const int q = (int)(std::find(set.begin(), set.end(), j) - set.begin());
					M[(size_t)p * m + q] += _values[k];
					used[k - st] = 1;
				}
			}
			for(int n = st; n < ed; ++n)
			{
				int kb = -1;
				for(int k = st; k < ed; ++k)
					if(!used[k - st] && (kb < 0 || std::fabs(_values[k]) > std::fabs(_values[kb]))) kb = k;
				if(kb < 0) break;
				used[kb - st] = 1;
				s += _values[kb] * x[_colInd[kb]];
			}
			b[p] = rhs[i] - s;
		}
		if(!solvable) continue;
		if(m == 1) { x[rows[set[0]]] = b[0] / M[0]; continue; }

		//  A set of more than one row: Gaussian elimination with partial
		//  pivoting (the sets have a few rows)
		for(int col = 0; col < m && solvable; ++col)
		{
			int pv = col;
			for(int p = col + 1; p < m; ++p)
				if(std::fabs(M[(size_t)p * m + col]) > std::fabs(M[(size_t)pv * m + col])) pv = p;
			if(M[(size_t)pv * m + col] == 0.0) { solvable = false; break; }
			if(pv != col)
			{
				for(int q = 0; q < m; ++q) std::swap(M[(size_t)pv * m + q], M[(size_t)col * m + q]);
				std::swap(b[pv], b[col]);
			}
			for(int p = col + 1; p < m; ++p)
			{
				const double f = M[(size_t)p * m + col] / M[(size_t)col * m + col];
				if(f == 0.0) continue;
				for(int q = col; q < m; ++q) M[(size_t)p * m + q] -= f * M[(size_t)col * m + q];
				b[p] -= f * b[col];
			}
		}
		if(!solvable) continue;
		for(int p = m - 1; p >= 0; --p)
		{
			double s = b[p];
			for(int q = p + 1; q < m; ++q) s -= M[(size_t)p * m + q] * x[rows[set[q]]];
			x[rows[set[p]]] = s / M[(size_t)p * m + p];
		}
	}
}

//------------------------------------------------------------------
//  Sparse matrix-vector product
//------------------------------------------------------------------
void CRS_Matrix::MatVec(Vector<double>& y, const Vector<double>& x) const
{
	assert(y.size() == x.size());
	assert(_n_row == x.size());
	Finalize();

	const int    *rp = &_rowPtr[0];
	const int    *ci = _colInd.empty() ? 0 : &_colInd[0];
	const double *va = _values.empty() ? 0 : &_values[0];
	const double *xx = x.data();
	double       *yy = y.data();

#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for schedule(static) if(yuzen::opt::Par(_n_row))
#endif
	for(int i = 0; i < _n_row; i++)
	{
		double sum = 0.0;
		const int ed = rp[i+1];
		for(int k = rp[i]; k < ed; ++k) sum += va[k] * xx[ci[k]];
		yy[i] = sum;
	}
}


Vector<double> CRS_Matrix::operator*(const Vector<double>& v) const
{
	assert(_n_row == v.size());
	Vector<double> tmp(_n_row);
	MatVec(tmp, v);
	return tmp;
}


void CRS_Matrix::Write(const std::string &wfile) const
{
	std::ofstream fout(wfile.c_str());
	fout<<"irow \t jcol \t m_ind[irow][jcol] \t m_A[irow][jcol] "<<std::endl;
	fout<<scientific;
	for(int irow = 0; irow < RowSize(); irow++)
	{
		for(int jcol = 0; jcol < this->NonZero(irow); jcol++)
			fout<<setw(5)<<irow<<"\t"<<setw(5)<<jcol<<"\t"
			    <<setw(5)<<this->Jcol(irow, jcol)<<"\t"<<this->Element(irow, jcol)<<std::endl;
	}
	fout.close();
}


//------------------------------------------------------------------
//  Build the permutation for the sorted mirror. This is only done when the
//  structure has changed.
//------------------------------------------------------------------
void CRS_Matrix::BuildSortPermutation() const
{
	Finalize();
	const int nnz = (int)_colInd.size();
	if((int)_sortedCol.size() == nnz && !_sortedCol.empty()) return;	// already built

	_sortedCol.resize(nnz);
	_sortedSrc.resize(nnz);
	_sortedVal.resize(nnz);

	std::vector<std::pair<int,int> > tmp;
	for(int i = 0; i < _n_row; ++i)
	{
		const int st = _rowPtr[i], ed = _rowPtr[i+1];
		tmp.clear();
		tmp.reserve(ed - st);
		for(int k = st; k < ed; ++k) tmp.push_back(std::pair<int,int>(_colInd[k], k));
		std::sort(tmp.begin(), tmp.end());
		for(int k = st; k < ed; ++k)
		{
			_sortedCol[k] = tmp[k - st].first;
			_sortedSrc[k] = tmp[k - st].second;
		}
	}
}

const int* CRS_Matrix::SortedColIndData() const
{
	BuildSortPermutation();
	return _sortedCol.empty() ? 0 : &_sortedCol[0];
}

const int* CRS_Matrix::SortedSrcData() const
{
	BuildSortPermutation();
	return _sortedSrc.empty() ? 0 : &_sortedSrc[0];
}

const double* CRS_Matrix::GatherSortedValues() const
{
	BuildSortPermutation();
	const int nnz = (int)_sortedSrc.size();
	//  A gather, every entry on its own: the result does not depend on how
	//  the loop is split between threads.
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for schedule(static) if(yuzen::opt::Par(nnz))
#endif
	for(int k = 0; k < nnz; ++k) _sortedVal[k] = _values[_sortedSrc[k]];
	return _sortedVal.empty() ? 0 : &_sortedVal[0];
}
