/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Yuzen.cpp : top-level driver. Reads the input, creates the result directory, picks a model and runs it
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "ThreadCap.h"
#include "Opt.h"
#include "Grain.h"
using namespace yuzen::util;
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "Yuzen.h"

using namespace yuzen;
using namespace Global;

C_Yuzen::C_Yuzen() : _errorTee(0), _model(0)
{

}

C_Yuzen::~C_Yuzen()
{
	delete _model;
	// Point FatalError back at standard error before the tee is destroyed,
	// so that it is never left holding a stream that no longer exists.
	if(_errorTee)
	{
		FatalError.Build(std::cerr, "--> FATAL ERROR", messageStream::FATAL);
		delete _errorTee;
		_errorTee = 0;
	}
}


bool C_Yuzen::ReadInputFile(const FilePath &rfile)
{
	_inputFile = rfile;
	Sout<<"# Reading input file : "<<_inputFile<<" ... ";
	if (_inputParam.Read(_inputFile) != 0)
	{
		WarningIn("ReadInputFile() is failed.");
		return false;
	}
	Sout<<" finished."<<nl;

	Sout<<"### Check consistency #####"<<nl;
	if(!CheckParameters::Check(_inputParam))
	{
		Global::FatalError<<"Input mismatch found. Program will abort."<<nl;
		Global::FatalError.Abort();
	}

	// Tell JobInfo which directory the results go into
	FilePath dir = _inputParam.JobName;
	JobInfo.SetResultDir(dir);

	return true;
}

bool C_Yuzen::CreateDirAndFiles()
{
	FilePath resDir = JobInfo.ResultDir();

	//--------------------------------------------------------------
	//  Every way the result directory can be unusable is found here,
	//  before anything is computed, and stops the run with the path, the
	//  key that names it and the reason. (A directory that exists but
	//  cannot be written to used to run to the end with exit code 0 and
	//  no output at all.) None of these messages can go to
	//  FatalErrors.txt, which is in the directory that failed.
	//--------------------------------------------------------------
	Sout<<nl<<"# Create the result directory .... ";
	auto refuse = [&](const std::string &what) {
		Sout<<nl;
		std::cerr<<"--> FATAL ERROR: the result directory \""<<resDir<<"\" (JobName) "<<what<<nl;
		return false;
	};
	std::error_code ec;
	const std::filesystem::path dirPath(static_cast<const std::string&>(resDir));
	if(std::filesystem::exists(dirPath, ec))
	{
		if(!std::filesystem::is_directory(dirPath, ec))
			return refuse("exists and is not a directory. Remove it or change JobName.");
		if(_inputParam.Output.Overwrite == false)
			return refuse("already exists, and Output.Overwrite is false. Rename or remove it,"
			              " or set Output.Overwrite = true.");
	}
	else
	{
		std::filesystem::create_directories(dirPath, ec);
		if(ec || !std::filesystem::is_directory(dirPath))
			return refuse("cannot be created: " + (ec ? ec.message() : std::string("unknown reason")) + ".");
	}
	{
		std::string why;
		if(!System::CanWriteIn(resDir, why))
			return refuse("cannot be written to: " + why + ".");
	}
	Sout<<"OK."<<nl;

	//Keep a copy of the input file in the result directory
	Sout<<nl<<"# Copy the input file to result directory ... ";
	{
		//  The input file may already be the copy (run from inside the result
		//  directory of an earlier run): nothing to do then.
		const std::filesystem::path in(static_cast<const std::string&>(_inputFile));
		const std::filesystem::path copy = dirPath / in.filename();
		if(std::filesystem::exists(copy, ec) && std::filesystem::equivalent(in, copy, ec))
			Sout<<"OK (it is the input file itself)."<<nl;
		else if(System::Copy(_inputFile, resDir)) Sout<<"OK."<<nl;
		else return refuse("could not take a copy of the input file " + _inputFile + ".");
	}

	// Set up JobInfo
	std::string file = "JobInfo.txt";
	_jobInfoStream.open((resDir / file).c_str());
	if(!_jobInfoStream.good()) return refuse("could not take JobInfo.txt: " + System::LastIOError() + ".");
	JobInfo.Build(_jobInfoStream, "", messageStream::INFO);

	//--------------------------------------------------------------
	//  Where the various messages are sent
	//
	//  Info stays on standard output: its content is the progress log
	//  itself. If that log is wanted as a file, redirecting at run time
	//  (yuzen case.toml > run.log) is the straightforward way, and it still
	//  allows watching the run in a terminal.
	//
	//  Warning goes to Warnings.txt. Fatal errors go to both FatalErrors.txt
	//  and standard error: the file alone leaves nothing on the terminal, and
	//  the terminal alone leaves no record in the result directory.
	//--------------------------------------------------------------
	Info.Build(std::cout, "", messageStream::INFO);

	// Set up Warning (change the output stream)
	file = "Warnings.txt";
	_warningStream.open((resDir / file).c_str());
	Warning.Build(_warningStream, "", messageStream::WARNING);

	// Set up FatalError (to both the file and standard error)
	file = "FatalErrors.txt";
	_errorStream.open((resDir / file).c_str());
	if(_errorStream.good())
	{
		delete _errorTee;
		_errorTee = new TeeStream(_errorStream, std::cerr);
		FatalError.Build(*_errorTee, "--> FATAL ERROR", messageStream::FATAL);
	}

	//--------------------------------------------------------------
	//  Do not write the logo into Warnings.txt or FatalErrors.txt.
	//
	//  Either file is then empty after a run where nothing went wrong, so
	//  content in it means that something really did happen.
	//--------------------------------------------------------------

	//  The algorithms YUZEN_OPT turned off (util/Opt.h), when it is set
	if(const char *e = std::getenv("YUZEN_OPT"))
	{
		std::string bad;
		const unsigned f = yuzen::opt::Parse(e, &bad);
		Info<<" YUZEN_OPT : "<<yuzen::opt::Names(f)<<nl;
		if(!bad.empty())
			Warning<<"YUZEN_OPT: not recognized, ignored: "<<bad<<nl;
	}

#if defined(_OPENMP) || defined(OMP)
	//  The thread limit of util/ThreadCap.h, set here, before the model
	//  allocates and first touches its arrays. A request above it that
	//  OMP_NUM_THREADS made has stopped the run already (Initialize).
	{
		const int requested = omp_get_max_threads();
		const char *capEnv  = std::getenv("YUZEN_THREAD_CAP");
		const int perThread = yuzen::util::CellsPerThreadFromEnv(capEnv);
		const int socket    = yuzen::util::CoresPerSocket();
		const bool sizeTerm = yuzen::util::SizeTermApplies(capEnv, yuzen::grain::KmpRuntime());
		const long long cells = GridCells();
		const int used      = yuzen::util::CappedThreads(cells, requested, perThread, socket, sizeTerm);
		if(used != requested) omp_set_num_threads(used);
		JobInfo.SetThreads(requested, perThread);

		Info<<" Running with OpenMP: "<<nl;
		Info<<" Number of CPUs:   "<<omp_get_num_procs() <<nl;
		Info<<" OMP_MAX_THREADS : "<< omp_get_max_threads();
		if(used != requested && sizeTerm)
			Info<<" (lowered from the default "<<requested<<": "<<cells<<" cells; at most "
			    <<yuzen::util::ThreadLimit(cells, perThread, socket)<<" threads, sec:parallel:cap)";
		else if(used != requested)
			Info<<" (lowered from the default "<<requested<<": the cores of one socket, sec:parallel:cap)";
		Info<<nl;
	}
#endif
	//  The thread count of linear solves by their size (util/Grain.h)
	yuzen::grain::Start();
#if defined(_OPENMP) || defined(OMP)
	Info<<" OpenMP runtime  : "<<yuzen::grain::RuntimeName()<<"; linear solves: "
	    <<yuzen::grain::Names(yuzen::grain::Get())<<nl;
#endif
	if(const char *e = std::getenv("YUZEN_GRAIN"))
	{
		std::string bad;
		yuzen::grain::Parse(e, 0, &bad);
		if(!bad.empty())
			Warning<<"YUZEN_GRAIN: not recognized, ignored: "<<bad<<nl;
	}
#if !(defined(_OPENMP) || defined(OMP))
	Info<<" Running in serial. "<<nl;
#endif


	return true;
}

long long C_Yuzen::GridCells() const
{
	const Input::I_Grid &g = _inputParam.Grid;
	long long cells = g.Nx;
	if(g.Dimension == "xyh") cells *= g.Ny;
	if(_inputParam.Physics.ConcField.Model == "CZ") cells *= g.Nz;
	return cells;
}

std::string C_Yuzen::RestartFileProblem() const
{
	if(!_inputParam.Restart.Enable) return "";
	const std::string &f = _inputParam.Restart.FileName;
	std::error_code ec;
	const std::filesystem::path p(f);
	if(!std::filesystem::exists(p, ec))       return "does not exist";
	if(std::filesystem::is_directory(p, ec))  return "is a directory, not a field file";
	if(!std::ifstream(f.c_str()).good())      return "cannot be read";
	return "";
}

bool C_Yuzen::Initialize()
{
	//  A restart file that is not there stops the run before the result
	//  directory is made: otherwise the directory is made, set up and left
	//  behind by a run that never starts.
	{
		const std::string why = RestartFileProblem();
		if(!why.empty())
		{
			std::cerr<<"--> FATAL ERROR: the restart file \""<<_inputParam.Restart.FileName
			         <<"\" (Restart.FileName, relative to the working directory) "<<why<<"."<<nl;
			return false;
		}
	}

	//  A thread count asked for (OMP_NUM_THREADS) above the limit of
	//  util/ThreadCap.h stops the run before the result directory is made:
	//  the cores asked for would sit idle, or the run be slower on them.
	//  Left to the OpenMP default, the run uses the limit (CreateDirAndFiles).
#if defined(_OPENMP) || defined(OMP)
	if(std::getenv("OMP_NUM_THREADS") != 0)
	{
		const int requested = omp_get_max_threads();
		const char *capEnv  = std::getenv("YUZEN_THREAD_CAP");
		const int perThread = yuzen::util::CellsPerThreadFromEnv(capEnv);
		const int socket    = yuzen::util::CoresPerSocket();
		const bool sizeTerm = yuzen::util::SizeTermApplies(capEnv, yuzen::grain::KmpRuntime());
		const long long cells = GridCells();
		const int limit = yuzen::util::CappedThreads(cells, requested, perThread, socket, sizeTerm);
		if(limit < requested && !sizeTerm)
		{
			std::cerr<<"--> FATAL ERROR: OMP_NUM_THREADS="<<requested<<" is more threads than one socket has ("
			         <<socket<<" cores); across sockets the run is slower. Run with OMP_NUM_THREADS="<<socket
			         <<" and give the other socket another case, or set YUZEN_THREAD_CAP=0 to run on "
			         <<requested<<" all the same."<<nl;
			return false;
		}
		if(limit < requested)
		{
			std::cerr<<"--> FATAL ERROR: OMP_NUM_THREADS="<<requested<<" is more threads than this run can use: "
			         <<cells<<" grid cells";
			if(socket >= yuzen::util::kMinThreads) std::cerr<<", "<<socket<<" cores a socket";
			std::cerr<<", so at most "<<limit<<" threads (one per "<<perThread<<" cells, at least "
			         <<yuzen::util::kMinThreads
			         <<(socket >= yuzen::util::kMinThreads ? ", within one socket" : "")
			         <<"). Run with OMP_NUM_THREADS="<<limit
			         <<" and give the other cores another case, or set YUZEN_THREAD_CAP=0 to run on "
			         <<requested<<" all the same."<<nl;
			return false;
		}
	}
#endif

	if(CreateDirAndFiles() == false) return false;

	//Opening information for JobInfo
	JobInfo<<Logo()<<nl;

	JobInfo<<"###### Program Info ######\n";
	JobInfo.Add("ProgramVersion",    Version(),              "Program Version");
	JobInfo.Add("ProgramBuilt",      System::BuildDate(),    "Program Built");
	JobInfo.Add("Compiler",          System::CompilerName(), "Compiler");
	JobInfo<<"\n";
	std::string RB_Ver = yuzen::config::Block::libraryName() + "-" + yuzen::config::Block::libraryVersion();
	JobInfo.Add("InputReader",  RB_Ver,                 "Ver. <input reader>");
	JobInfo.Add("GridVersion",       VersionGrid(),          "Ver. <grid>");
	JobInfo.Add("UtilVersion",       yuzen::util::Version(),      "Ver. <util>");


	//  Which model is solved follows from two settings in the input. Choosing
	//  the actual object is left to MakeSolverModel() (model/SolverModel.cpp).
	ISolverModel *m = MakeSolverModel(_inputParam.Physics.ConcField.Model,
	                                  _inputParam.Grid.Dimension);
	if(m == 0)
	{
		FatalError<<"Incorrect model: "<<_inputParam.Physics.ConcField.Model.c_str()
		          <<", "<<_inputParam.Grid.Dimension<<nl;
		Global::FatalError.Abort();
		return false;
	}
	delete _model;
	_model = m;

	//  Build first, then report. ConcZ decides its name from the grid it has
	//  just built (ConcZ-2D or ConcZ-3D), so asking before Build() would
	//  give the 1-D film name.
	_model->Build(_inputParam);
	Info<<"Model: "<<_model->Name()<<nl;

	JobInfo.Start();


	return true;
}

void C_Yuzen::Finalize()
{
	if(Interrupt::Requested()) JobInfo.End(std::string("interrupted (") + Interrupt::Name() + ")");
	else JobInfo.End();
}

bool C_Yuzen::Run()
{
	//  From here on Ctrl-C (SIGINT) and SIGTERM stop the time loop after
	//  the step it is in, with the output of the state reached written as
	//  at the end of a run, so that it can be restarted from.
	Interrupt::Install();

	if(Initialize() == false) return false;

	_model->Run();

	Finalize();

	//  Report failure if the run ended by diverging. The caller turns this
	//  into the exit code, so a batch job or a test can tell it from a success.
	return !_model->Diverged() && !Interrupt::Requested();
}
