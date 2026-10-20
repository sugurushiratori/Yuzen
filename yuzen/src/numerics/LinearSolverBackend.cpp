/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LinearSolverBackend.cpp : the implementations of the linear solver
*
*  Every class here is confined to this file and can only be created from
*  outside through MakeLinearSolverBackend().
*
*    EigenBackend  : Eigen's BiCGSTAB, preconditioned with IncompleteLUT,
*                    Diagonal or Identity. The matrix is referenced through
*                    Eigen::Map, with no copy.
*    AmgclBackend  : AMGCL's BiCGStab, with an algebraic multigrid
*                    preconditioner. The preconditioner itself runs in
*                    parallel under OpenMP as well.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <algorithm>
#include "Profile.h"
#include "Opt.h"
#include <cmath>
#include <cstdlib>
#include <vector>
#include <memory>
#include <tuple>
#include <stdexcept>

#include "YuzenUtil.h"
using namespace yuzen::util;

#include "LinearSolverBackend.h"
#include "CRS_Matrix.h"

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include "Grain.h"

#ifdef YUZEN_WITH_AMGCL
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#endif

namespace {

//==================================================================
//  The Eigen backend
//
//  The matrix is referenced through Eigen::Map onto the sorted mirror of the
//  CRS_Matrix. Eigen's sparse algorithms assume that the column numbers
//  within a row are ascending, so the matrix cannot be handed over in its own
//  storage order. Filling the mirror costs one O(nnz) gather.
//
//  Diagonal scaling is applied to the mirror's values and to a copy of the
//  right-hand side, leaving the matrix itself intact. Scaling by rows does
//  not change the solution.
//==================================================================
typedef Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, int> > ConstSpMap;
typedef Eigen::Map<Eigen::VectorXd>       VecMap;
typedef Eigen::Map<const Eigen::VectorXd> ConstVecMap;

// Extra configuration, per preconditioner.
// The default fillfactor of IncompleteLUT is 10, which on a row with a
// 13-point stencil allows more than 100 entries of fill. ILU is ILU(0):
// it factorizes only on the existing non-zero pattern. Allowing almost no fill therefore puts the amount of work on the
// same footing.
template<class Precond> inline void ConfigPrecond(Precond &, int) {}
template<> inline void ConfigPrecond<Eigen::IncompleteLUT<double,int> >(
		Eigen::IncompleteLUT<double,int> &p, int fill)
{
	p.setFillfactor(fill);
	p.setDroptol(1e-10);
}

template<class Precond>
int RunEigenBiCGSTAB(const ConstSpMap &A, const double *b, double *x, int n,
                     int maxIter, double tol, int fill, int &iter, double &resid)
{
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor, int>, Precond> solver;
	solver.setMaxIterations(maxIter);
	solver.setTolerance(tol);
	ConfigPrecond(solver.preconditioner(), fill);
	{ YUZEN_PROF("lin.precond"); solver.compute(A); }
	if(solver.info() != Eigen::Success)
	{
		iter = 0; resid = 0;
		return 2;		// building the preconditioner failed
	}

	ConstVecMap bb(b, n);
	VecMap      xx(x, n);
	{ YUZEN_PROF("lin.iterate"); xx = solver.solveWithGuess(bb, xx); }

	iter  = (int)solver.iterations();
	resid = solver.error();

	if(solver.info() == Eigen::Success) return 0;
	if(!(resid == resid)) return 2;		// NaN
	return 1;							// hit the iteration limit
}

//------------------------------------------------------------------
//  opt bilu (util/Opt.h): block-Jacobi ILUT. The rows are split into a number of
//  contiguous blocks that follows from n alone (one per kRowsPerBlock rows,
//  at most kMaxBlocks), each diagonal block is factorized by Eigen's
//  IncompleteLUT with the settings of the ILUT preconditioner, and the
//  preconditioner applies the blocks' factors, each to its own rows. The
//  blocks are independent, so they are factorized and applied in parallel,
//  and since the blocking does not depend on the thread count neither does
//  the result. The coupling between blocks is left to the Krylov iteration.
//------------------------------------------------------------------
class BlockILUT
{
public:
	static const int kRowsPerBlock = 1024;
	static const int kMaxBlocks    = 128;

	BlockILUT() : _n(0), _nnz(0), _fill(10), _ok(false) {}

	//  Rebuild the blocking when the pattern changes; refactorize always.
	//  Returns false when a block's factorization failed.
	bool Setup(int n, int nnz, const int *rowPtr, const int *colInd, const double *val, int fill)
	{
		const bool same = (_n == n) && (_nnz == nnz) && (_fill == fill) &&
			std::equal(rowPtr, rowPtr + n + 1, _rowPtr.begin()) &&
			std::equal(colInd, colInd + nnz, _colInd.begin());
		if(!same)
		{
			_n = n; _nnz = nnz; _fill = fill;
			_rowPtr.assign(rowPtr, rowPtr + n + 1);
			_colInd.assign(colInd, colInd + nnz);
			int nb = n / kRowsPerBlock;
			if(nb < 1) nb = 1;
			if(nb > kMaxBlocks) nb = kMaxBlocks;
			_begin.resize((size_t)nb + 1);
			for(int b = 0; b <= nb; ++b) _begin[(size_t)b] = (int)((long long)n * b / nb);
			_blk.assign((size_t)nb, Block());
			for(int b = 0; b < nb; ++b)
			{
				Block &B = _blk[(size_t)b];
				const int r0 = _begin[(size_t)b], r1 = _begin[(size_t)b+1];
				B.ptr.assign(1, 0);
				for(int i = r0; i < r1; ++i)
				{
					for(int k = rowPtr[i]; k < rowPtr[i+1]; ++k)
						if(colInd[k] >= r0 && colInd[k] < r1) { B.col.push_back(colInd[k] - r0); B.src.push_back(k); }
					B.ptr.push_back((int)B.col.size());
				}
				B.val.resize(B.col.size());
				B.ilut.reset();
			}
		}
		const int nb = (int)_blk.size();
		int bad = 0;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for schedule(dynamic) reduction(+:bad)
#endif
		for(int b = 0; b < nb; ++b)
		{
			Block &B = _blk[(size_t)b];
			const int m = _begin[(size_t)b+1] - _begin[(size_t)b];
			for(size_t k = 0; k < B.src.size(); ++k) B.val[k] = val[B.src[k]];
			ConstSpMap A(m, m, (int)B.col.size(), B.ptr.data(), B.col.data(), B.val.data());
			if(!B.ilut)
			{
				B.ilut.reset(new Eigen::IncompleteLUT<double, int>());
				ConfigPrecond(*B.ilut, _fill);
				B.ilut->compute(A);		// analyzePattern and factorize (see SolveILUT)
			}
			else B.ilut->factorize(A);
			if(B.ilut->info() != Eigen::Success) { B.ilut.reset(); bad++; }
		}
		_ok = (bad == 0);
		return _ok;
	}

	//  The preconditioner Eigen's bicgstab() applies
	template<class V> Eigen::VectorXd solve(const V &r) const
	{
		Eigen::VectorXd y(_n);
		const int nb = (int)_blk.size();
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for schedule(dynamic)
#endif
		for(int b = 0; b < nb; ++b)
		{
			const int r0 = _begin[(size_t)b], m = _begin[(size_t)b+1] - r0;
			y.segment(r0, m) = _blk[(size_t)b].ilut->solve(r.segment(r0, m));
		}
		return y;
	}

private:
	struct Block
	{
		std::vector<int> ptr, col, src;
		std::vector<double> val;
		std::shared_ptr<Eigen::IncompleteLUT<double, int> > ilut;
	};
	int _n, _nnz, _fill;
	bool _ok;
	std::vector<int> _rowPtr, _colInd, _begin;
	std::vector<Block> _blk;
};

//------------------------------------------------------------------
//  The values and the right-hand side handed to a backend, each row divided
//  by its diagonal when `scale` is set (a row whose diagonal is below 1e-30
//  in magnitude is copied as it is). Row by row in parallel: every entry is
//  one multiplication of its own, so the result does not depend on the
//  number of threads. With `src`, values are those of the matrix in its
//  storage order and the k-th one handed over is values[src[k]] (the gather
//  into column order done here, in the same loop); with x0, x0 is copied to
//  xOut in the same loop as well (opt regions, util/Opt.h).
//------------------------------------------------------------------
void ScaleRows(const CRS_Matrix &crs, const int *rowPtr, const double *values,
               const double *rhs, int n, bool scale, double *valOut, double *rhsOut,
               const int *src = 0, const double *x0 = 0, double *xOut = 0)
{
	const double epsw = 1.0e-30;
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for schedule(static) if(yuzen::opt::Par(n))
#endif
	for(int i = 0; i < n; ++i)
	{
		if(x0) xOut[i] = x0[i];
		const int st = rowPtr[i], ed = rowPtr[i+1];
		const double d = scale ? crs.Diag(i) : 1.0;
		if(!scale || std::fabs(d) < epsw)
		{
			if(src) for(int k = st; k < ed; ++k) valOut[k] = values[src[k]];
			else    for(int k = st; k < ed; ++k) valOut[k] = values[k];
			rhsOut[i] = rhs[i];
			continue;
		}
		const double inv = 1.0 / d;
		if(src) for(int k = st; k < ed; ++k) valOut[k] = values[src[k]] * inv;
		else    for(int k = st; k < ed; ++k) valOut[k] = values[k] * inv;
		rhsOut[i] = rhs[i] * inv;
	}
}

//------------------------------------------------------------------
//  When a cheaper preconditioner tried first fails, the next tries are
//  skipped for a while: 16 solves after the first failure, doubling after
//  each further one up to 1024, back to 16 after a try that converged. A
//  system on which it does not work (a ConcZ film with dry columns) then
//  costs a few failed tries instead of one per solve. The decision follows
//  from the solves' convergence only (no thread count in it).
//------------------------------------------------------------------
struct TryBackoff
{
	int skip = 0, wait = 16;
	bool Try() { if(skip > 0) { --skip; return false; } return true; }
	void Failed() { skip = wait; wait = (wait < 1024) ? 2 * wait : 1024; }
	void Converged() { wait = 16; }
};

class EigenBackend : public ILinearSolverBackend
{
public:
	EigenBackend() : _ilutFill(10), _reuse(true), _patN(0), _patNnz(0) {}

	virtual const char* Name() const { return "eigen"; }

	virtual void Build(int nRow, int nonZero, int nonZL, int nonZU)
	{
		(void)nonZero; (void)nonZL; (void)nonZU;
		_scaledRhs.resize((size_t)nRow);
		// An escape hatch for tuning. Normally the default of 10 is left alone.
		// A smaller value makes the preconditioner cheaper, but dropping it
		// too far breaks the factorization and the solution comes out NaN
		// (a value of 1 breaks down even on a 1-D problem).
		//  An escape hatch for comparative measurements. Passing 0 makes
		//  every solve call compute(), the factorization, again.
		if(const char *e = std::getenv("YUZEN_ILUT_REUSE"))
			_reuse = (std::atoi(e) != 0);
		if(const char *e = std::getenv("YUZEN_ILUT_FILL"))
		{
			const int v = std::atoi(e);
			if(v >= 1) _ilutFill = v;
		}
	}

	//------------------------------------------------------------------
	//  The solve used when ILUT is the preconditioner.
	//
	//  IncompleteLUT::compute() calls analyzePattern() and then factorize().
	//  What the first of those does is
	//
	//      mat1 = A;  mat2 = A^T;  AtA = mat1 + mat2;  AMD ordering
	//
	//  and it is **determined by the sparsity structure alone; it does not
	//  depend on the values**. The structure, however, never changes after the
	//  first assembly (ClearData zeroes the values but keeps the structure),
	//  so rebuilding this on every solve is wasted work from beginning to end.
	//  Measured on ConcZ, this analyzePattern accounted for roughly 30% of the
	//  total run time.
	//
	//  So when the structure is the same as last time, the solver is reused
	//  and only factorize() is called. The same structure always produces the
	//  same ordering, so **the solution does not change, bit for bit**.
	//------------------------------------------------------------------
	int SolveILUT(const LinearSolverSettings &s, const ConstSpMap &A,
	              const int *rowPtr, const int *colInd, int n, int nnz,
	              const double *b, double *x, int &iter, double &resid)
	{
		const bool samePattern =
			(_ilut.get() != 0) && (_patN == n) && (_patNnz == nnz) &&
			std::equal(rowPtr, rowPtr + n + 1, _patRowPtr.begin()) &&
			std::equal(colInd, colInd + nnz,   _patColInd.begin());

		if(!samePattern)
		{
			_ilut.reset(new IlutSolver());
			ConfigPrecond(_ilut->preconditioner(), _ilutFill);
			//  compute(), not analyzePattern() followed by factorize().
			//  IncompleteLUT (Eigen 3.4.0) neither initializes m_info in its
			//  constructor nor sets it in analyzePattern(), and
			//  IterativeSolverBase::analyzePattern() copies it straight after
			//  calling that -- a load of an uninitialized enum, which UBSan
			//  reports and which is undefined behavior however harmless it is
			//  in practice. compute() runs the same analyzePattern() and
			//  factorize() on the preconditioner, in the same order, and only
			//  then reads its status, so the factorization and therefore the
			//  solution are bit-identical to the two-call form.
			{ YUZEN_PROF("lin.precond"); _ilut->compute(A); }
			_patN = n; _patNnz = nnz;
			_patRowPtr.assign(rowPtr, rowPtr + n + 1);
			_patColInd.assign(colInd, colInd + nnz);
		}
		else
		{
			YUZEN_PROF("lin.precond");
			_ilut->factorize(A);
		}

		_ilut->setMaxIterations(s.maxIteration);
		_ilut->setTolerance(s.epsilon);
		if(_ilut->info() != Eigen::Success)
		{
			_ilut.reset(); iter = 0; resid = 0; return 2;
		}

		ConstVecMap bb(b, n);
		VecMap      xx(x, n);
		{ YUZEN_PROF("lin.iterate"); xx = _ilut->solveWithGuess(bb, xx); }

		iter  = (int)_ilut->iterations();
		resid = _ilut->error();

		if(_ilut->info() == Eigen::Success) return 0;
		if(!(resid == resid)) { _ilut.reset(); return 2; }
		return 1;
	}

	virtual int Solve(const LinearSolverSettings &s, CRS_Matrix &crs,
	                  Vector<double> &rhs, Vector<double> &x, int &iter, double &resid)
	{
		yuzen::grain::SolveScope threads((long long)crs.RowSize());	// util/Grain.h
		const int n   = crs.RowSize();
		const int nnz = crs.NonZero();

		const int    *rowPtr = crs.RowPtrData();
		const int    *colInd = crs.SortedColIndData();
		const double *values = 0;
		//  With opt regions (util/Opt.h) and scaling, the gather into column order
		//  is done by the scaling loop itself
		const bool fuse = s.scaling && yuzen::opt::On(yuzen::opt::Regions);
		if(!fuse) { YUZEN_PROF("lin.prepare"); values = crs.GatherSortedValues(); }

		const double *bptr = rhs.data();

		if(s.scaling)
		{
			YUZEN_PROF("lin.prepare");
			// Divide each row by its diagonal. A and b are divided by the same
			// factor, so the solution is unchanged.
			_scaledVal.resize((size_t)nnz);
			_scaledRhs.resize((size_t)n);
			ScaleRows(crs, rowPtr, fuse ? crs.ValuesData() : values, rhs.data(), n, true,
			          _scaledVal.empty() ? 0 : &_scaledVal[0],
			          _scaledRhs.empty() ? 0 : &_scaledRhs[0],
			          fuse ? crs.SortedSrcData() : 0);
			values = _scaledVal.empty() ? 0 : &_scaledVal[0];
			bptr   = _scaledRhs.empty() ? 0 : &_scaledRhs[0];
		}

		ConstSpMap A(n, n, nnz, rowPtr, colInd, values);

		//--------------------------------------------------------------
		//  The cheaper preconditioners first, the one they replace as the
		//  fall-back: opt concpc (util/Opt.h) tries Jacobi on a
		//  concentration equation, opt bilu block-Jacobi ILUT; a try that
		//  does not converge within kTryIter (300) iterations (or breaks down) is
		//  dropped, x is set back to the initial guess, and the solve goes
		//  on as without the option. Which way a solve went is decided by
		//  its own convergence only, so it does not depend on the thread
		//  count. The iterations of a dropped try are counted in iter.
		//--------------------------------------------------------------
		const int kTryIter = (s.maxIteration < 300) ? s.maxIteration : 300;
		int iterTried = 0;
		Vector<double> x0;
		const bool tryJacobi = s.conc && s.precondition == LinearSolverSettings::tp_ILU
		                    && yuzen::opt::On(yuzen::opt::ConcPC) && _offJacobi.Try();
		const bool tryBlock  = s.precondition == LinearSolverSettings::tp_ILU && yuzen::opt::On(yuzen::opt::BlockILU)
		                    && _offBlock.Try();
		if(tryJacobi || tryBlock) x0 = x;
		if(tryJacobi)
		{
			const int r = RunEigenBiCGSTAB<Eigen::DiagonalPreconditioner<double> >(
			           A, bptr, x.data(), n, kTryIter, s.epsilon, 0, iter, resid);
			if(r == 0) { _offJacobi.Converged(); return 0; }
			_offJacobi.Failed();
			iterTried += iter;
			x = x0;
		}
		if(tryBlock)
		{
			bool ok;
			{ YUZEN_PROF("lin.precond"); ok = _bilu.Setup(n, nnz, rowPtr, colInd, values, _ilutFill); }
			if(ok)
			{
				ConstVecMap bb(bptr, n);
				Eigen::VectorXd xx = VecMap(x.data(), n);
				Eigen::Index it = kTryIter;
				double err = s.epsilon;
				bool ret;
				{ YUZEN_PROF("lin.iterate"); ret = Eigen::internal::bicgstab(A, bb, xx, _bilu, it, err); }
				if(ret && err == err && err <= s.epsilon)
				{
					VecMap(x.data(), n) = xx;
					iter = (int)it + iterTried; resid = err;
					_offBlock.Converged();
					return 0;
				}
				iterTried += (int)it;
			}
			_offBlock.Failed();
			x = x0;
		}
		if(iterTried > 0)
		{
			YUZEN_PROF(tryJacobi ? "lin.fallback.conc" : "lin.fallback.bilu");
			const int r = SolveDefault(s, A, rowPtr, colInd, n, nnz, bptr, x, iter, resid);
			iter += iterTried;
			return r;
		}
		return SolveDefault(s, A, rowPtr, colInd, n, nnz, bptr, x, iter, resid);
	}

	//  The solve with the preconditioner of the input
	int SolveDefault(const LinearSolverSettings &s, const ConstSpMap &A, const int *rowPtr,
	                 const int *colInd, int n, int nnz, const double *bptr, Vector<double> &x,
	                 int &iter, double &resid)
	{
		switch(s.precondition)
		{
		case LinearSolverSettings::tp_None:
			return RunEigenBiCGSTAB<Eigen::IdentityPreconditioner>(
			           A, bptr, x.data(), n, s.maxIteration, s.epsilon, 0, iter, resid);
		case LinearSolverSettings::tp_Jacobi:
			return RunEigenBiCGSTAB<Eigen::DiagonalPreconditioner<double> >(
			           A, bptr, x.data(), n, s.maxIteration, s.epsilon, 0, iter, resid);
		case LinearSolverSettings::tp_ILU:
		default:
			if(_reuse)
				return SolveILUT(s, A, rowPtr, colInd, n, nnz, bptr, x.data(), iter, resid);
			return RunEigenBiCGSTAB<Eigen::IncompleteLUT<double, int> >(
			           A, bptr, x.data(), n, s.maxIteration, s.epsilon, _ilutFill, iter, resid);
		}
	}

private:
	typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor, int>,
	                        Eigen::IncompleteLUT<double, int> > IlutSolver;

	int _ilutFill;		// fill factor of IncompleteLUT
	bool _reuse;		// reuse of the pattern analysis (YUZEN_ILUT_REUSE=0 turns it off)
	// Working storage holding the values and the right-hand side after
	// diagonal scaling, so that the matrix itself is left intact
	std::vector<double, yuzen::util::FirstTouchAllocator<double> > _scaledVal;
	std::vector<double, yuzen::util::FirstTouchAllocator<double> > _scaledRhs;

	//  The solver being reused, and the sparsity structure it analyzed
	std::unique_ptr<IlutSolver> _ilut;
	BlockILUT _bilu;		// opt bilu (util/Opt.h)
	TryBackoff _offJacobi, _offBlock;
	int _patN, _patNnz;
	std::vector<int> _patRowPtr, _patColInd;
};


#ifdef YUZEN_WITH_AMGCL
//==================================================================
//  The AMGCL backend
//
//  In AMGCL's builtin backend both the SpMV and the preconditioner are
//  parallelized with OpenMP. ILU is inherently sequential, so this is the
//  backend to use when the run has to scale in parallel.
//  The preconditioners map across as
//    ILU    -> smoothed-aggregation AMG with SPAI0 (fully parallel)
//    Jacobi -> damped Jacobi
//    None   -> no preconditioner
//==================================================================
class AmgclBackend : public ILinearSolverBackend
{
public:
	virtual const char* Name() const { return "amgcl"; }

	virtual void Build(int nRow, int nonZero, int nonZL, int nonZU)
	{ (void)nRow; (void)nonZero; (void)nonZL; (void)nonZU; }

	virtual int Solve(const LinearSolverSettings &s, CRS_Matrix &crs,
	                  Vector<double> &rhs, Vector<double> &x, int &iter, double &resid)
	{
		yuzen::grain::SolveScope threads((long long)crs.RowSize());	// util/Grain.h
		typedef amgcl::backend::builtin<double> Backend;

		const int n   = crs.RowSize();
		const int nnz = crs.NonZero();

		YUZEN_PROF("lin.amgcl");
		//  The row pointers and column indices are handed over in place
		//  (AMGCL copies the matrix into its own storage while it builds the
		//  preconditioner); the values and the right-hand side are written,
		//  scaled or not, into storage kept from one solve to the next.
		const int *rowPtr = crs.RowPtrData();
		const int *colInd = crs.SortedColIndData();
		{
		YUZEN_PROF("lin.prepare");
		_val.resize((size_t)nnz);
		_b.resize((size_t)n);
		_x.resize((size_t)n);
		if(yuzen::opt::On(yuzen::opt::Regions))
		{
			//  The gather, the scaling and the copy of x in one loop
			ScaleRows(crs, rowPtr, crs.ValuesData(), rhs.data(), n, s.scaling,
			          _val.empty() ? 0 : &_val[0], _b.empty() ? 0 : &_b[0],
			          crs.SortedSrcData(), x.data(), _x.empty() ? 0 : &_x[0]);
		}
		else
		{
		const double *values = crs.GatherSortedValues();
		ScaleRows(crs, rowPtr, values, rhs.data(), n, s.scaling,
		          _val.empty() ? 0 : &_val[0], _b.empty() ? 0 : &_b[0]);
		const double *x0 = x.data();
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for schedule(static)
#endif
		for(int i = 0; i < n; ++i) _x[(size_t)i] = x0[i];
		}
		}
		const auto ptr = amgcl::make_iterator_range(rowPtr, rowPtr + n + 1);
		const auto col = amgcl::make_iterator_range(colInd, colInd + nnz);
		const auto val = amgcl::make_iterator_range(_val.data(), _val.data() + nnz);
		std::vector<double> &b = _b, &xx = _x;

		size_t it = 0;
		double err = 0.0;

		//  opt concpc (util/Opt.h): SPAI0 as the preconditioner of a
		//  concentration equation first (set-up: one pass over the rows);
		//  if it does not converge within 300 iterations, x back to the
		//  initial guess and the preconditioner of the input (as with Eigen)
		size_t itTried = 0;
		bool done = false;
		if(s.conc && s.precondition == LinearSolverSettings::tp_ILU
		&& yuzen::opt::On(yuzen::opt::ConcPC) && _offSpai.Try())
		{
			typedef amgcl::make_solver<
				amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::spai0>,
				amgcl::solver::bicgstab<Backend> > Solver;
			Solver::params prm;
			prm.solver.maxiter = (s.maxIteration < 300) ? s.maxIteration : 300;
			prm.solver.tol     = s.epsilon;
			std::unique_ptr<Solver> solve;
			//  A breakdown (AMGCL throws "Zero rho in BiCGStab") is a try
			//  that did not converge
			try
			{
				{ YUZEN_PROF("lin.precond"); solve.reset(new Solver(std::tie(n, ptr, col, val), prm)); }
				{ YUZEN_PROF("lin.iterate"); std::tie(it, err) = (*solve)(b, xx); }
			}
			catch(const std::exception &) { err = 1e300; }
			if(err <= s.epsilon) { done = true; _offSpai.Converged(); }
			else
			{
				YUZEN_PROF("lin.fallback.conc");
				_offSpai.Failed();
				itTried = it;
				const double *x0 = x.data();
				for(int i = 0; i < n; ++i) xx[(size_t)i] = x0[i];
			}
		}
		if(done) {}
		else if(s.precondition == LinearSolverSettings::tp_None
		|| s.precondition == LinearSolverSettings::tp_Jacobi)
		{
			typedef amgcl::make_solver<
				amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::damped_jacobi>,
				amgcl::solver::bicgstab<Backend> > Solver;
			Solver::params prm;
			prm.solver.maxiter = s.maxIteration;
			prm.solver.tol     = s.epsilon;
			std::unique_ptr<Solver> solve;
			{ YUZEN_PROF("lin.precond"); solve.reset(new Solver(std::tie(n, ptr, col, val), prm)); }
			{ YUZEN_PROF("lin.iterate"); std::tie(it, err) = (*solve)(b, xx); }
		}
		else
		{
			typedef amgcl::make_solver<
				amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0>,
				amgcl::solver::bicgstab<Backend> > Solver;
			Solver::params prm;
			prm.solver.maxiter = s.maxIteration;
			prm.solver.tol     = s.epsilon;
			//  opt amgcoarse (util/Opt.h): coarsen to 500 unknowns before the
			//  direct (sequential) solve of the coarsest level, not 3000 --
			//  only for a system larger than 3000. One that small AMGCL solves
			//  directly, in one iteration; coarsened to 500 the fourth-order
			//  1-D film took 80 BiCGStab iterations a solve and ran 6.5
			//  times slower on one thread (tools/scaling bench1d_n02000)
			if(yuzen::opt::On(yuzen::opt::AmgCoarse) && (size_t)n > prm.precond.coarse_enough)
				prm.precond.coarse_enough = 500;
			std::unique_ptr<Solver> solve;
			{ YUZEN_PROF("lin.precond"); solve.reset(new Solver(std::tie(n, ptr, col, val), prm)); }
			{ YUZEN_PROF("lin.iterate"); std::tie(it, err) = (*solve)(b, xx); }
		}

		{
			YUZEN_PROF("lin.prepare");
			double *xo = x.data();
#if defined(_OPENMP) || defined(OMP)
#pragma omp parallel for schedule(static) if(yuzen::opt::Par(n))
#endif
			for(int i = 0; i < n; ++i) xo[i] = xx[(size_t)i];
		}
		iter  = (int)(it + itTried);
		resid = err;
		if(err <= s.epsilon) return 0;
		if(!(err == err))    return 2;
		return 1;
	}

private:
	//  Storage kept between solves, so that a solve allocates nothing here
	//  (b and x are std::vector: AMGCL's builtin backend takes no other
	//  allocator; the matrix values are the bulk of the data)
	std::vector<double, yuzen::util::FirstTouchAllocator<double> > _val;
	std::vector<double> _b, _x;
	TryBackoff _offSpai;		// opt concpc: SPAI0 tried first
};
#endif //YUZEN_WITH_AMGCL

} // namespace


//==================================================================
//  Factory
//==================================================================
ILinearSolverBackend* MakeLinearSolverBackend(const std::string &name, int inPlane)
{
	if(name == "auto")
	{
		//------------------------------------------------------------
		//  Choose by the dimension of the film (measured on a 112-core Xeon
		//  Max with libomp, up to one socket of threads).
		//
		//  1-D: eigen. On a system that is close to tridiagonal there is
		//       almost nothing for algebraic multigrid to coarsen: at Nx=8000
		//       AMGCL needs about 100 linear iterations a solve where eigen
		//       needs 12, and is 5.7 times slower. Only at ~1e5 cells and on
		//       many threads does amgcl pull ahead (1.5x at Nx=131072 on 56).
		//  2-D: amgcl. The two are level at 128^2; from 256^2 up eigen is
		//       1.2 to 1.6 times faster on one thread and amgcl 1.1 to 1.5
		//       times faster on a socket of threads, the count a run takes by
		//       default.
		//
		//  Where the crossovers sit depends on the machine and the thread
		//  count, and in 2-D the gain is small either way. Putting a threshold
		//  here would only make that threshold itself a machine-dependent
		//  guess, so the choice is made on the dimension alone.
		//  (The measurements are in the manual, sec:parallel:backend.)
		//------------------------------------------------------------
#ifdef YUZEN_WITH_AMGCL
		if(inPlane >= 2) return new AmgclBackend();
		return new EigenBackend();
#else
		(void)inPlane;
		return new EigenBackend();
#endif
	}

	if(name == "amgcl")
	{
#ifdef YUZEN_WITH_AMGCL
		return new AmgclBackend();
#else
		Global::Warning<<"LinearSolver: Backend=\"amgcl\" was requested, but this "
		                 "build has no AMGCL. Falling back to \"eigen\"."<<nl;
		return new EigenBackend();
#endif
	}

	return new EigenBackend();
}
