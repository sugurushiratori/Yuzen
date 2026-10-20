/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  mainProg.cpp : entry point
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#include "YuzenUtil.h"
using namespace yuzen::util;
#include "Yuzen.h"
#include "Grain.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(_WIN32)
#include <windows.h>
#endif


namespace {

void Usage(std::ostream &os, const char *argv0)
{
	os<<"usage: "<<argv0<<" [--check-input] <InputFile>"<<std::endl
	  <<"       "<<argv0<<" --help | --version"<<std::endl;
}

//  A run that cannot go on: the reason through FatalError (the terminal and,
//  once the result directory is set up, FatalErrors.txt), then exit code 1
int Failed(const std::string &why)
{
	Global::FatalError<<std::endl<<"--> FATAL ERROR: "<<why<<std::endl;
	Global::JobInfo.Abort();
	return 1;
}

//  The whole run, from the command line to the exit code
int Main(int argc, char*argv[])
{
	//  --help and --version answer and stop before the banner, so that
	//  their output is only what was asked for.
	const std::vector<std::string> args(argv + 1, argv + argc);
	for(const std::string &a : args)
	{
		if(a == "--help" || a == "-h")
		{
			Usage(std::cout, argv[0]);
			std::cout<<"  Runs the case described by <InputFile> (TOML). The results go into the"<<std::endl
			         <<"  directory named by its JobName, relative to the working directory."<<std::endl
			         <<"  --check-input  read and check <InputFile>, run nothing"<<std::endl
			         <<"  --help, -h     this text"<<std::endl
			         <<"  --version      the version of this build"<<std::endl;
			return 0;
		}
		if(a == "--version")
		{
			std::cout<<"Yuzen-"<<yuzen::C_Yuzen::Version()<<std::endl;
			return 0;
		}
	}

	//  --check-input reads the file and checks it -- unknown keys, ranges,
	//  accepted values, the combinations CheckParameters refuses -- and stops
	//  there: nothing is built, written or run. Every failure exits with 1
	//  from inside ReadInputFile(), as it does for a run.
	bool checkOnly = false;
	std::vector<std::string> files;
	for(const std::string &a : args)
	{
		if(a == "--check-input") checkOnly = true;
		else if(a.size() > 1 && a[0] == '-')
		{
			std::cerr<<"--> FATAL ERROR: unknown option "<<a<<std::endl;
			Usage(std::cerr, argv[0]);
			return 1;
		}
		else files.push_back(a);
	}
	if(files.size() != 1)
	{
		if(files.size() > 1)
			std::cerr<<"--> FATAL ERROR: one input file is read per run, "<<files.size()<<" were given"<<std::endl;
		Usage(std::cerr, argv[0]);
		return 1;
	}
	const FilePath rfile = files[0];

	yuzen::C_Yuzen solver;

	Sout<<yuzen::C_Yuzen::Logo();
	Sout<<" Yuzen-"<<yuzen::C_Yuzen::Version();
	//  The same stamp as JobInfo.txt: __DATE__/__TIME__ here would be those
	//  of this file, compiled at another moment than the library.
	Sout<<" ("<<System::BuildDate()<<")"<<std::endl;
	// Guarded so that this still compiles in a build without OpenMP
#if defined(_OPENMP)
	Sout<<" OpenMP: enabled,  omp_get_max_threads = " << omp_get_max_threads() <<std::endl;
#else
	Sout<<" OpenMP: disabled (serial build)" <<std::endl;
#endif
	Sout<<std::endl;

	//  The input file is looked at before the parser is: a parser handed a
	//  directory, or a file it cannot open, says little about why.
	{
		std::error_code ec;
		const std::filesystem::path p(static_cast<const std::string&>(rfile));
		const char *why = 0;
		if(!std::filesystem::exists(p, ec))            why = "does not exist";
		else if(std::filesystem::is_directory(p, ec))  why = "is a directory, not a file";
		else
		{
			std::ifstream f(rfile.c_str());
			if(!f.good()) why = "cannot be read (no permission?)";
		}
		if(why)
		{
			std::cerr<<"--> FATAL ERROR: the input file \""<<rfile<<"\" "<<why<<"."<<std::endl;
			return 1;
		}
	}

	if(!solver.ReadInputFile(rfile)) return 1;
	if(checkOnly)
	{
		//  The restart file is named relative to the working directory of the
		//  run, which need not be this one: noted, not refused.
		const std::string why = solver.RestartFileProblem();
		if(!why.empty())
			Sout<<"Note: the restart file (Restart.FileName) "<<why
			    <<" relative to this working directory; a run started here would stop."<<std::endl;
		Sout<<"Input check passed: "<<rfile<<std::endl;
		return 0;
	}

	//  The return value of Run() is the exit code: a set-up that failed or a
	//  run that diverged must not end with 0, so that a batch job can tell
	//  those cases apart from a success. A run stopped by SIGINT or SIGTERM
	//  ends, after writing its output, with 128 + the signal number, which is
	//  what the shell reports for a process killed by that signal.
	//  Running out of memory -- a grid too large for the machine -- throws
	//  std::bad_alloc wherever the allocation happens. Caught here, while
	//  the solver and its FatalErrors.txt still exist, it ends the run with
	//  a message there and on the terminal, "abort" in JobInfo.txt and exit
	//  code 1, instead of std::terminate and a core dump.
	bool ok = false;
	try
	{
		ok = solver.Run();
	}
	catch(const std::bad_alloc &)
	{
		return Failed("out of memory. The grid (Grid.Nx, Ny, Nz) is too large for the memory"
		              " of this machine.");
	}
	catch(const std::length_error &e)
	{
		return Failed(std::string("an array too large to allocate (") + e.what()
		              + "). The grid (Grid.Nx, Ny, Nz) is too large.");
	}
	if(Interrupt::Requested())
	{
		Sout<<"Stopped by "<<Interrupt::Name()<<"."<<std::endl;
		return 128 + Interrupt::Signal();
	}
	return ok ? 0 : 1;
}

} // namespace


int main(int argc, char*argv[])
{
#if defined(_WIN32)
	//  The paths in the messages are UTF-8 (yuzen.manifest); so is the console
	SetConsoleOutputCP(CP_UTF8);
#endif
	//  Hot teams for the LLVM and Intel OpenMP runtimes, before their start
	yuzen::grain::PrepareRuntime();
	//  Anything else thrown out of the run, as a message and exit code 1
	try
	{
		return Main(argc, argv);
	}
	catch(const std::exception &e)
	{
		std::cerr<<std::endl<<"--> FATAL ERROR: "<<e.what()<<std::endl;
	}
	return 1;
}
