/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  LinearSolver.cpp : the front door to the linear solver
*
*  This file goes as far as holding the settings and choosing the
*  implementation from the strings in the input file. The solving itself is
*  in LinearSolverBackend.cpp.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <cassert>

#include "YuzenUtil.h"
using namespace yuzen::util;

#include "LinearSolver.h"
#include "CRS_Matrix.h"


LinearSolver::LinearSolver()
 : _set(), _backend(MakeLinearSolverBackend("eigen")), _constructed(false),
   _nRow(0), _nonZero(0), _nonZL(0), _nonZU(0)
{}

LinearSolver::~LinearSolver()
{
	delete _backend;
}


void LinearSolver::SetBackend(const std::string &s, int inPlane)
{
	ILinearSolverBackend *next = MakeLinearSolverBackend(s, inPlane);
	delete _backend;
	_backend = next;
	// When Build() has already run, allocate the working storage again for
	// the implementation that was just selected.
	if(_constructed) _backend->Build(_nRow, _nonZero, _nonZL, _nonZU);
}

void LinearSolver::SetPrecondition(const std::string &s)
{
	if     (s == "ILU")    _set.precondition = LinearSolverSettings::tp_ILU;
	else if(s == "Jacobi") _set.precondition = LinearSolverSettings::tp_Jacobi;
	else                   _set.precondition = LinearSolverSettings::tp_None;
}

std::string LinearSolver::BackendName() const
{
	return _backend->Name();
}


void LinearSolver::Build(int nRow, int nonZero, int nonZL, int nonZU)
{
	_nRow = nRow; _nonZero = nonZero; _nonZL = nonZL; _nonZU = nonZU;
	_backend->Build(nRow, nonZero, nonZL, nonZU);
	_constructed = true;
}


int LinearSolver::Solve(CRS_Matrix &crs, Vector<double> &rhs, Vector<double> &x,
                        int &iter, double &resid)
{
	assert(_constructed);
	return _backend->Solve(_set, crs, rhs, x, iter, resid);
}
