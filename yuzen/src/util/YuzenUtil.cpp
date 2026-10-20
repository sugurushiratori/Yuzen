/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  YuzenUtil.cpp : implementation of the logging, paths and job information
*
*  YuzenUtil.cpp
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include "YuzenUtil.h"
#include "Opt.h"
#include "Grain.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <csignal>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <iomanip>

#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace yuzen { namespace util {

std::string Version() { return "YuzenUtil-1.0"; }

//==================================================================
//  FilePath
//==================================================================
void FilePath::Normalize()
{
	std::replace(begin(), end(), '\\', '/');
}

FilePath::tType FilePath::Type() const
{
	std::error_code ec;
	const fs::path p(static_cast<const std::string&>(*this));
	if(fs::is_symlink(p, ec))   return LINK;
	if(fs::is_directory(p, ec)) return DIRECTORY;
	if(fs::is_regular_file(p, ec)) return FILE;
	return UNDEFINED;
}

std::string FilePath::Name() const
{
	return fs::path(static_cast<const std::string&>(*this)).filename().string();
}

FilePath FilePath::Path() const
{
	return FilePath(fs::path(static_cast<const std::string&>(*this)).parent_path().string());
}

FilePath FilePath::WithoutExt() const
{
	const fs::path p(static_cast<const std::string&>(*this));
	return FilePath((p.parent_path() / p.stem()).string());
}

std::string FilePath::Ext() const
{
	return fs::path(static_cast<const std::string&>(*this)).extension().string();
}

namespace {
FilePath join(const std::string &a, const std::string &b)
{
	if(a.empty()) return FilePath(b);
	if(b.empty()) return FilePath(a);
	std::string s = a;
	if(s[s.size()-1] != '/' && s[s.size()-1] != '\\') s += '/';
	s += b;
	return FilePath(s);
}
} // namespace

FilePath operator/(const std::string &a, const std::string &b) { return join(a, b); }
FilePath operator/(const char *a,        const std::string &b) { return join(a ? a : "", b); }
FilePath operator/(const std::string &a, const char *b)        { return join(a, b ? b : ""); }


//==================================================================
//  System
//==================================================================
bool System::Exists(const FilePath &p)
{
	std::error_code ec;
	return fs::exists(fs::path(static_cast<const std::string&>(p)), ec);
}
bool System::IsDir(const FilePath &p)
{
	std::error_code ec;
	return fs::is_directory(fs::path(static_cast<const std::string&>(p)), ec);
}
bool System::IsFile(const FilePath &p)
{
	std::error_code ec;
	return fs::is_regular_file(fs::path(static_cast<const std::string&>(p)), ec);
}
bool System::MkDir(const FilePath &p)
{
	std::error_code ec;
	const fs::path path(static_cast<const std::string&>(p));
	if(fs::exists(path, ec)) return fs::is_directory(path, ec);
	return fs::create_directories(path, ec);
}
bool System::Copy(const FilePath &from, const FilePath &to)
{
	std::error_code ec;
	const fs::path src(static_cast<const std::string&>(from));
	fs::path dst(static_cast<const std::string&>(to));
	if(fs::is_directory(dst, ec)) dst /= src.filename();
	fs::copy_file(src, dst, fs::copy_options::overwrite_existing, ec);
	return !ec;
}

std::string System::HostName()
{
#if defined(_WIN32)
	const char *h = std::getenv("COMPUTERNAME");
	return h ? h : "unknown";
#else
	char buf[256] = {0};
	if(::gethostname(buf, sizeof(buf) - 1) == 0) return buf;
	return "unknown";
#endif
}

std::string System::OS_Name()
{
#if defined(_WIN32)
	return "Windows";
#elif defined(__APPLE__)
	return "macOS";
#elif defined(__linux__)
	return "Linux";
#else
	return "unknown";
#endif
}

std::string System::CompilerName()
{
	std::ostringstream oss;
#if defined(__clang__)
	oss << "clang-" << __clang_major__ << "." << __clang_minor__;
#elif defined(__INTEL_COMPILER)
	oss << "icc-" << __INTEL_COMPILER;
#elif defined(__GNUC__)
	oss << "gcc-" << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#elif defined(_MSC_VER)
	oss << "msvc-" << _MSC_VER;
#else
	oss << "unknown";
#endif
	return oss.str();
}

std::string System::BuildDate()
{
	std::ostringstream oss;
	oss << "build at: " << __TIME__ << " on " << __DATE__;
	return oss.str();
}

std::string System::CPU_Info()
{
#if defined(__linux__)
	std::ifstream fin("/proc/cpuinfo");
	std::string line;
	while(std::getline(fin, line))
	{
		if(line.compare(0, 10, "model name") == 0)
		{
			const size_t c = line.find(':');
			if(c != std::string::npos)
			{
				std::string s = line.substr(c + 1);
				while(!s.empty() && s[0] == ' ') s.erase(0, 1);
				return s;
			}
		}
	}
#endif
	return "unknown";
}

std::string System::UserName()
{
	const char *u = std::getenv("USER");
	if(!u) u = std::getenv("USERNAME");
	return u ? u : "";
}

long System::GetPID()
{
#if defined(_WIN32)
	return (long)::_getpid();
#else
	return (long)::getpid();
#endif
}

double System::MemoryTotalMB()
{
#if defined(__linux__)
	std::ifstream fin("/proc/meminfo");
	std::string key; long kb = 0;
	while(fin >> key)
	{
		if(key == "MemTotal:") { fin >> kb; return kb / 1024.0; }
		std::string rest; std::getline(fin, rest);
	}
#endif
	return 0.0;
}

double System::MemoryUsedMB()
{
#if defined(__linux__)
	std::ifstream fin("/proc/self/statm");
	long pages = 0, resident = 0;
	if(fin >> pages >> resident) return resident * (double)::sysconf(_SC_PAGESIZE) / (1024.0 * 1024.0);
#endif
	return 0.0;
}

FilePath System::Cwd()
{
	std::error_code ec;
	return FilePath(fs::current_path(ec).string());
}

bool System::ReplaceFile(const std::string &tmp, const std::string &path, std::string &why)
{
	std::error_code ec;
	//  std::filesystem::rename replaces an existing file on POSIX and on
	//  Windows alike (where std::rename would refuse to)
	fs::rename(fs::path(tmp), fs::path(path), ec);
	if(!ec) return true;
	//  ... except with some older standard libraries on Windows: there the
	//  file is removed first (it is then absent for a moment, never cut short)
	std::error_code ec2;
	if(fs::exists(fs::path(path), ec2) && fs::remove(fs::path(path), ec2))
	{
		fs::rename(fs::path(tmp), fs::path(path), ec);
		if(!ec) return true;
	}
	why = ec.message();
	return false;
}

std::string System::LastIOError()
{
	const int e = errno;
	if(e != 0) return std::strerror(e);
	return "an input/output error (is the disk full, or the directory not writable?)";
}

bool System::CanWriteIn(const FilePath &dir, std::string &why)
{
	const std::string probe = static_cast<const std::string&>(dir / ".yuzen_write_test");
	errno = 0;
	std::FILE *f = std::fopen(probe.c_str(), "wb");
	if(!f) { why = LastIOError(); return false; }
	const bool ok = std::fputs("yuzen\n", f) >= 0 && std::fflush(f) == 0;
	if(!ok) why = LastIOError();
	std::fclose(f);
	std::remove(probe.c_str());
	return ok;
}


//==================================================================
//  Interrupt
//==================================================================
namespace {
volatile std::sig_atomic_t g_signal = 0;

extern "C" void yuzenOnSignal(int s)
{
	g_signal = s;
	//  The second Ctrl-C ends the process at once
	std::signal(s, SIG_DFL);
}
} // namespace

void Interrupt::Install()
{
	g_signal = 0;
	std::signal(SIGINT,  yuzenOnSignal);
	std::signal(SIGTERM, yuzenOnSignal);
#if defined(SIGXFSZ)
	std::signal(SIGXFSZ, SIG_IGN);
#endif
}

int Interrupt::Signal() { return (int)g_signal; }

const char* Interrupt::Name()
{
	switch((int)g_signal)
	{
	case SIGINT:  return "SIGINT";
	case SIGTERM: return "SIGTERM";
	default:      return "a signal";
	}
}


//==================================================================
//  messageStream
//==================================================================
messageStream::messageStream() : _streamPtr(&std::cout), _title(""), _mtype(INFO) {}

messageStream::messageStream(std::ostream &os, const std::string &title, mesType t)
	: _streamPtr(&os), _title(title), _mtype(t) {}

void messageStream::Build(std::ostream &os, const std::string &title, mesType t)
{
	_streamPtr = &os;
	_title     = title;
	_mtype     = t;
}

messageStream::operator std::ostream&()
{
	return *_streamPtr;
}

std::ostream& messageStream::operator()(const char *funcName, const char *srcName, int srcLine)
{
	return operator()(std::string(funcName ? funcName : ""), srcName, srcLine);
}

std::ostream& messageStream::operator()(const std::string &funcName, const char *srcName, int srcLine)
{
	std::ostream &os = *_streamPtr;
	if(!_title.empty()) os << _title << " ";
	os << "in [" << funcName << "] : [" << (srcName ? srcName : "") << ":" << srcLine << "]:: ";
	return os;
}


//==================================================================
//  Error
//==================================================================
Error::Error(const std::string &title, bool swAbort, bool swThrow)
	: messageStream(std::cerr, title, FATAL), _abort(swAbort), _throwExceptions(swThrow) {}

std::ostream& Error::operator()(const char *funcName, const char *srcName, int line)
{
	return messageStream::operator()(funcName, srcName, line);
}
std::ostream& Error::operator()(const std::string &funcName, const char *srcName, int line)
{
	return messageStream::operator()(funcName, srcName, line);
}

void Error::Abort()
{
	// Write to this logger's own stream. If the tee onto both FatalErrors.txt
	// and stderr has been set up, the record of the abort reaches both as well.
	static_cast<std::ostream&>(*this)
		<< std::endl << "####### Aborting ############" << std::endl;
	Global::JobInfo.Abort();
	if(_throwExceptions) throw *this;
	if(_abort) std::exit(EXIT_FAILURE);
}


//==================================================================
//  C_JobInfo
//==================================================================
bool C_JobInfo::constructed = false;

namespace {

double nowSeconds()
{
	using clock = std::chrono::steady_clock;
	return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

std::string dateString()
{
	const std::time_t t = std::time(0);
	std::tm tmv{};
#if defined(_WIN32)
	localtime_s(&tmv, &t);
#else
	localtime_r(&t, &tmv);
#endif
	char buf[64];
	std::strftime(buf, sizeof(buf), "%b %d %Y", &tmv);
	return buf;
}

std::string timeString()
{
	const std::time_t t = std::time(0);
	std::tm tmv{};
#if defined(_WIN32)
	localtime_s(&tmv, &t);
#else
	localtime_r(&t, &tmv);
#endif
	char buf[64];
	std::strftime(buf, sizeof(buf), "%H:%M:%S", &tmv);
	return buf;
}

} // namespace

C_JobInfo::C_JobInfo() : _resultDir("not given"), _startSeconds(nowSeconds()) {}

void C_JobInfo::Add(std::string key, std::string str, const std::string &mesTitle)
{
	_infoList[key] = str;
	_streamPtr->setf(std::ios::left);
	*_streamPtr << std::setw(25) << (mesTitle + " :") << "\t" << str << std::endl;
}
void C_JobInfo::Add(std::string key, int val, const std::string &mesTitle)
{
	std::ostringstream oss; oss << val;
	_infoList[key] = oss.str();
	_streamPtr->setf(std::ios::left);
	*_streamPtr << std::setw(25) << (mesTitle + " :") << "\t" << val << std::endl;
}
void C_JobInfo::Add(std::string key, double val, const std::string &mesTitle)
{
	std::ostringstream oss; oss << val;
	_infoList[key] = oss.str();
	_streamPtr->setf(std::ios::left);
	*_streamPtr << std::setw(25) << (mesTitle + " :") << "\t" << val << std::endl;
}

void C_JobInfo::MachineInfo()
{
	*_streamPtr << std::endl << std::endl << "###### Machine Info ######" << std::endl;
	Add("OS_Name",  System::OS_Name(),  "Operation System");
	Add("HostName", System::HostName(), "Host Name");
	Add("CPUName",  System::CPU_Info(), "CPU Name");
#if defined(_OPENMP)
	Add("NumProc",          omp_get_num_procs(),   "Number of processors");
	Add("OMP_MAX_THREADS ", omp_get_max_threads(), "OMP_MAX_THREADS");
	if(_threadsRequested > 0)
	{
		Add("OMP_THREADS_REQUESTED", _threadsRequested, "OMP_THREADS_REQUESTED");
		Add("OMP_THREAD_CAP", _cellsPerThread, "OMP_THREAD_CAP");
	}
#endif
	Add("YUZEN_OPT", yuzen::opt::Names(yuzen::opt::Flags()), "YUZEN_OPT");
	Add("OMP_RUNTIME", yuzen::grain::RuntimeName(), "OMP_RUNTIME");
	Add("YUZEN_GRAIN", yuzen::grain::Names(yuzen::grain::Get()), "YUZEN_GRAIN");
	Add("TotalMemory", System::MemoryTotalMB(), "Total Memory [MB]");
	Add("UsedMemory",  System::MemoryUsedMB(),  "Used  Memory [MB]");
	*_streamPtr << std::endl;
	Add("UserName", System::UserName(),      "User Name");
	Add("PID",      (int)System::GetPID(),   "Process ID");
}

void C_JobInfo::Start()
{
	MachineInfo();
	*_streamPtr << std::endl << std::endl << "###### Starting Info ######" << std::endl;
	Add("StartDate", dateString(),     "Start Date");
	Add("StartTime", timeString(),     "Start Time");
	Add("StartDir",  System::Cwd(),    "Starting Directory");
	Add("ResultDir", _resultDir,       "Result Directory");
	_startSeconds = nowSeconds();
	constructed = true;
}

void C_JobInfo::End(const std::string &terminationType)
{
	if(constructed)
	{
		*_streamPtr << std::endl << std::endl << "###### Ending Info ######" << std::endl;
		Add("ElapsedTime", nowSeconds() - _startSeconds, "Elapsed Time [s]");
		Add("EndDate",     dateString(),                 "End Date");
		Add("EndTime",     timeString(),                 "End Time");
		Add("Termination", terminationType,              "Termination Type");
	}
	constructed = false;
}

void C_JobInfo::End()   { End("normal"); }
void C_JobInfo::Exit()  { End("exit");   }
void C_JobInfo::Abort() { End("abort");  }

const std::string C_JobInfo::operator[](const std::string &key)
{
	tInfoList::iterator it = _infoList.find(key);
	if(it == _infoList.end())
	{
		WarningIn("C_JobInfo::operator[]") << "key: " << key << " was not found." << std::endl;
		return "Key Not Found";
	}
	return it->second;
}


//==================================================================
//  Global instances
//==================================================================
namespace Global {
	messageStream Info(std::cout, "", messageStream::INFO);
	messageStream Warning(std::cerr, "--> WARNING", messageStream::WARNING);
	Error         FatalError("--> FATAL ERROR", true, false);
	C_JobInfo     JobInfo;
}

} } // namespace yuzen::util
