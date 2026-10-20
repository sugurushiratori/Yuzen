/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Yuzen.h : top-level driver of the solver
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _DEF_YUZEN_INCLUDED_
#define _DEF_YUZEN_INCLUDED_

#include "YuzenUtil.h"
using namespace yuzen::util;
#include "Input.h"
#include "CheckParameters.h"
#include "SolverModel.h"

#include <iostream>
#include <fstream>


namespace yuzen
{
	//----------------------------------------------------------------------
	//  Top-level driver of the solver.
	//
	//  All that main has to do is this:
	//      C_Yuzen solver;
	//      solver.ReadInputFile(file);
	//      solver.Run();
	//
	//  It reads the input and checks it for consistency, creates the result
	//  directory and the logs, and then picks a model from the input's
	//  Grid.Dimension and Physics.ConcField.Model and runs it. Choosing the
	//  model is left to MakeSolverModel() (SolverModel.h), so this class does
	//  not know which models exist.
	//----------------------------------------------------------------------
	class C_Yuzen
	{
	public:
		/**
		 * @brief Builds a driver with no input read yet and no model selected.
		 *
		 * Note that the form taking a file is only declared, never defined, so
		 * using it fails at link time. Use the default constructor and then
		 * ReadInputFile().
		 */
		C_Yuzen();
		C_Yuzen(FilePath& file);
		/**
		 * @brief Destroys the driver and deletes the solver model it selected.
		 */
		virtual ~C_Yuzen();

		/**
		 * @brief Reads the input file and checks the parameters for consistency.
		 *
		 * Also derives the result directory from the job name and records it.
		 * Call this before Run().
		 *
		 * @param file  path of the input file
		 * @return      false if the file could not be read or the parameters are
		 *              inconsistent, in which case a warning has been issued
		 */
		//  Reads the input file and checks the values for consistency.
		//  On a mismatch it issues a warning and returns false.
		bool ReadInputFile(const FilePath& file);

		/**
		 * @brief Runs the whole computation: set-up, time integration, shutdown.
		 *
		 * Creates the result directory and the log files, picks the solver model
		 * from the input (Grid.Dimension and Physics.ConcField.Model) through
		 * MakeSolverModel(), runs it to completion, then closes the job log.
		 * Requires a successful ReadInputFile() first.
		 *
		 * @return false if set-up failed, for example because the result
		 *         directory already exists and overwriting is not allowed, if
		 *         the run diverged, or if it was stopped by SIGINT or SIGTERM
		 *         (Interrupt in util/YuzenUtil.h; the caller turns that into
		 *         exit code 128 + the signal)
		 */
		//  Creates the result directory and the logs, picks a model, runs it through
		bool Run();

		/**
		 * @brief Whether the restart file named in the input (Restart.FileName,
		 *        when Restart.Enable) is there, relative to the working directory.
		 * @return the reason it cannot be read, empty when it can (or when
		 *         Restart.Enable is false)
		 */
		std::string RestartFileProblem() const;
		/** @brief The grid cells of the run, as util/ThreadCap.h counts them. */
		long long GridCells() const;

		// static members

		/**
		 * @brief Solver version as "Major.Minor.Rev" (see Version.h).
		 */
		static std::string Version();		// version of the solver

		/**
		 * @brief Version series as "Major.Minor", without the patch number.
		 *
		 * This is what an input file names in SolverName, so that a patch
		 * release does not invalidate existing input files.
		 */
		static std::string VersionSeries();	// the unit of input-file compatibility

		/**
		 * @brief Multi-line ASCII banner printed at start-up and into the logs.
		 */
		static std::string Logo();			// logo printed at start-up

		/**
		 * @brief Version string of the grid library the solver was built against.
		 */
		static std::string VersionGrid();	// version of the grid library

	private:
		FilePath _inputFile;
		Input _inputParam;

		std::ofstream _jobInfoStream;
		std::ofstream _warningStream;
		std::ofstream _errorStream;
		// Fatal errors go to both FatalErrors.txt and standard error. This is
		// built after _errorStream and held by pointer so that it is not
		// destroyed before the stream it writes into.
		TeeStream    *_errorTee;

		ISolverModel *_model;	// the model chosen from the input

		//  Creates the result directory, copies the input file into it, and
		//  redirects the message destinations (JobInfo / Info / Warning /
		//  FatalError) to files inside it. Returns false, with the path and the
		//  reason on standard error, when the directory cannot be used: it
		//  exists and overwriting is not allowed, it is not a directory, it
		//  cannot be made or it cannot be written to.
		bool CreateDirAndFiles();

		//  Start-up sequence: CreateDirAndFiles(), then record the version
		//  information, then pick a model with MakeSolverModel() and Build()
		//  it. Called at the top of Run().
		bool Initialize();

		//  Closes out the job log
		void Finalize();

	private:
		C_Yuzen(const C_Yuzen&);
		void operator=(const C_Yuzen&);

	};
}

#endif //_DEF_YUZEN_INCLUDED_
