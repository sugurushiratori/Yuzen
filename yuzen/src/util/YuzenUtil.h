/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  YuzenUtil.h : logging, paths, filesystem and job information
*
*  What this layer provides is only
*    nl / Sout / Info / Warning / FatalError / FilePath / JobInfo /
*    System / sprintf_s
*  The path and filesystem work itself is left to C++17 <filesystem>.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _YUZEN_UTIL_INCLUDED_
#define _YUZEN_UTIL_INCLUDED_

#include <chrono>
#include <cstdio>
#include <exception>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#if defined(_OPENMP)
#include <omp.h>
#endif

//------------------------------------------------------------------
//  Map the MSVC "safe" string functions onto snprintf with a variadic macro.
//  The buffer length is passed along too, so the length check still applies.
//------------------------------------------------------------------
#ifndef _MSC_VER
#define sprintf_s(buf, ...) std::snprintf((buf), sizeof(buf), __VA_ARGS__)
#define sscanf_s            std::sscanf
#endif

#define Sout std::cout
#define Serr std::cerr
#define nl   std::endl

namespace yuzen { namespace util {

/**
 * @brief Returns the version tag of this utility layer.
 */
std::string Version();

//==================================================================
//  Perf : a clock used for performance measurement only. It takes no part in
//  the numerics.
//
//  The Elapsed Time in JobInfo covers building the grid, setting the initial
//  field and writing the initial output as well, so on a large grid the
//  set-up stops being negligible and that figure cannot be used as a measure
//  of parallel scaling. The model drivers use this clock instead to time the
//  time-step loop on its own, and write that under ###### Performance ######
//  in JobInfo.txt.
//==================================================================
namespace Perf {
	/** @brief Wall-clock seconds from a monotonic clock. Differences only. */
	inline double NowSeconds()
	{
		return std::chrono::duration<double>(
		           std::chrono::steady_clock::now().time_since_epoch()).count();
	}
}

//==================================================================
//  FilePath
//
//  This derives from std::string so that callers can pass it around either as
//  a std::string or as a FilePath. The work itself is left to std::filesystem.
//==================================================================
class FilePath : public std::string
{
public:
	/** What the path points at on disk. UNDEFINED means it does not exist. */
	enum tType { UNDEFINED, FILE, DIRECTORY, LINK };

	/**
	 * @brief Builds a path. The std::string and char* forms normalize the text,
	 *        turning backslashes into slashes; the copy form does not (the
	 *        source is already normalized). A null char* gives an empty path.
	 */
	FilePath() {}
	FilePath(const FilePath &p)    : std::string(p) {}
	FilePath(const std::string &s) : std::string(s) { Normalize(); }
	FilePath(const char *s)        : std::string(s ? s : "") { Normalize(); }

	/**
	 * @brief Assigns a new path, normalizing separators for the std::string and
	 *        char* forms as the corresponding constructors do.
	 */
	const FilePath& operator=(const FilePath &p)    { std::string::operator=(p); return *this; }
	const FilePath& operator=(const std::string &s) { std::string::operator=(s); Normalize(); return *this; }
	const FilePath& operator=(const char *s)        { std::string::operator=(s ? s : ""); Normalize(); return *this; }

	/**
	 * @brief Asks the filesystem what this path currently is.
	 * @return FILE, DIRECTORY, LINK, or UNDEFINED if nothing is there.
	 */
	tType       Type() const;
	/**
	 * @brief Last element of the path, extension included.
	 */
	std::string Name() const;		// last element, extension included
	/**
	 * @brief Parent directory of this path.
	 */
	FilePath    Path() const;		// parent directory
	/**
	 * @brief The path with the extension stripped, directory part kept.
	 */
	FilePath    WithoutExt() const;
	/**
	 * @brief Extension of the last element, including the leading dot
	 *        (empty when there is none).
	 */
	std::string Ext() const;		// extension, including the "."

private:
	void Normalize();				// turns backslashes into slashes
};

/**
 * @brief Joins two path fragments, inserting a single separator between them.
 *        An empty fragment gives back the other one unchanged.
 */
FilePath operator/(const std::string &a, const std::string &b);
FilePath operator/(const char *a,        const std::string &b);
FilePath operator/(const std::string &a, const char *b);


//==================================================================
//  System : only the methods that are actually used
//==================================================================
class System
{
public:
	// Existence and type queries. They never throw: a path that cannot be
	// examined (no permission, broken link, ...) simply reports false.
	static bool Exists(const FilePath &p);
	static bool IsDir (const FilePath &p);
	static bool IsFile(const FilePath &p);
	/**
	 * @brief Creates a directory, together with any missing parents.
	 * @return true if the directory now exists, including the case where it
	 *         already existed before the call.
	 */
	static bool MkDir (const FilePath &p);
	/**
	 * @brief Copies a single file, overwriting the destination if present.
	 * @param from  source file
	 * @param to  destination file; if it names an existing directory, the file
	 *            is copied into it under its original name.
	 * @return true on success.
	 */
	static bool Copy  (const FilePath &from, const FilePath &to);

	// Machine and job identification, written into JobInfo.txt.
	// Values that cannot be obtained on the current platform come back as
	// "unknown" (strings) or 0 (numbers) rather than failing.
	/** @brief Name of the machine this job is running on. */
	static std::string HostName();
	/** @brief Operating system family: Windows, macOS, Linux or unknown. */
	static std::string OS_Name();
	/** @brief Compiler name and version this binary was built with. */
	static std::string CompilerName();
	/** @brief Build timestamp baked in at compile time, as "build at: <time> on <date>". */
	static std::string BuildDate();
	/** @brief CPU model name (read from /proc/cpuinfo; Linux only). */
	static std::string CPU_Info();
	/** @brief Name of the user who started the job, taken from the environment. */
	static std::string UserName();
	/** @brief Process ID of the running job. */
	static long        GetPID();
	/** @brief Total physical memory in MB (Linux only; 0 elsewhere). */
	static double      MemoryTotalMB();
	/** @brief Resident memory of this process in MB (Linux only; 0 elsewhere). */
	static double      MemoryUsedMB();
	/** @brief Current working directory. */
	static FilePath    Cwd();

	// Writing files safely. See OutputVTK::Write for how they are used.
	/**
	 * @brief Puts a file written under a temporary name in the place of path,
	 *        replacing what was there.
	 *
	 * A file is written whole under tmp first and only then renamed, so that
	 * a run killed or out of disk space while writing leaves at worst a
	 * stray tmp file, never a file under its real name that looks whole and
	 * is not.
	 * @param tmp  the file written whole under its temporary name
	 * @param path the name it is to have
	 * @param why  set to the reason on failure
	 * @return true on success.
	 */
	static bool ReplaceFile(const std::string &tmp, const std::string &path, std::string &why);
	/**
	 * @brief The reason the last file operation failed, from errno ("No space
	 *        left on device", "Permission denied", ...), or a generic text
	 *        when errno was not set.
	 */
	static std::string LastIOError();
	/**
	 * @brief Whether a file can be created in the directory dir: one is
	 *        created, written and removed again.
	 * @param dir  the directory to try
	 * @param why  set to the reason on failure
	 */
	static bool CanWriteIn(const FilePath &dir, std::string &why);
};


//==================================================================
//  Interrupt : SIGINT (Ctrl-C) and SIGTERM end the run cleanly
//==================================================================
/**
 * @brief Catches SIGINT and SIGTERM so that the time loop can stop after the
 *        step it is in and write its output, as at the end of a run.
 *
 * The first signal only sets a flag, which the time loops test after every
 * step. The handler then puts the default action back, so a second signal
 * ends the process at once. SIGXFSZ (a file grown past `ulimit -f`) is
 * ignored, so that the write fails and is reported instead of killing the
 * process silently.
 */
class Interrupt
{
public:
	/** @brief Installs the handlers. Call once, before the time loop. */
	static void Install();
	/** @brief The signal received, 0 if none. */
	static int  Signal();
	/** @brief True once SIGINT or SIGTERM has been received. */
	static bool Requested() { return Signal() != 0; }
	/** @brief The name of the signal received ("SIGINT", "SIGTERM"). */
	static const char* Name();
};


//==================================================================
//  TeeStream : duplicates one output into two streams
//
//  This exists for the fatal errors. Writing to the file alone leaves nothing
//  on the terminal, and writing to the terminal alone leaves no record in the
//  result directory. Neither is acceptable, so both are written.
//==================================================================
/**
 * @brief An output stream that duplicates everything written to it into two
 *        underlying streams.
 *
 * Used for the fatal-error log, which has to reach both the terminal and the
 * result directory. Neither stream is owned, so both have to outlive the tee.
 * Flushing and sync are forwarded to both.
 */
class TeeStream : public std::ostream
{
public:
	/**
	 * @brief Builds a tee over two streams.
	 * @param a  first destination; neither stream is owned, and both have to
	 *           outlive this object
	 * @param b  second destination, on the same terms
	 */
	TeeStream(std::ostream &a, std::ostream &b) : std::ostream(&_buf), _buf(a, b) {}

private:
	class Buf : public std::streambuf
	{
	public:
		Buf(std::ostream &a, std::ostream &b) : _a(a.rdbuf()), _b(b.rdbuf()) {}
	protected:
		virtual int overflow(int c)
		{
			if(c == EOF) return !EOF;
			const int ra = _a->sputc(static_cast<char>(c));
			const int rb = _b->sputc(static_cast<char>(c));
			return (ra == EOF || rb == EOF) ? EOF : c;
		}
		virtual int sync()
		{
			const int ra = _a->pubsync();
			const int rb = _b->pubsync();
			return (ra == 0 && rb == 0) ? 0 : -1;
		}
	private:
		std::streambuf *_a;
		std::streambuf *_b;
	};
	Buf _buf;
};


//==================================================================
//  messageStream : a logger whose destination can be swapped out
//==================================================================
class messageStream
{
public:
	/** Severity of the messages this stream carries. */
	enum mesType { INFO, WARNING, FATAL };

	/**
	 * @brief Builds a logger targeting std::cout with an empty title and INFO
	 *        severity.
	 */
	messageStream();
	/**
	 * @brief Builds a logger, setting destination, title and severity at once.
	 * @param os     destination stream. It is held by pointer and is not owned,
	 *               so it has to outlive the logger.
	 * @param title  label printed ahead of each message
	 * @param t      severity the messages carry
	 */
	messageStream(std::ostream &os, const std::string &title, mesType t);
	virtual ~messageStream() {}

	/**
	 * @brief Retargets an already constructed logger at another stream.
	 * @param os     destination stream; not owned, and has to outlive the logger.
	 * @param title  label printed ahead of each message
	 * @param t      severity the messages carry
	 */
	void Build(std::ostream &os, const std::string &title, mesType t);

	// Severity and title of this logger.
	mesType            Type()  const { return _mtype; }
	const std::string& Title() const { return _title; }

	/**
	 * @brief Writes the message header and returns the stream to go on writing to.
	 *
	 * Normally reached through the InfoIn(fn) / WarningIn(fn) / FatalErrorIn(fn)
	 * macros, which fill in the source file and line for the caller.
	 *
	 * @param funcName  name of the reporting function, printed in the header
	 * @param srcName   source file name, usually __FILE__
	 * @param srcLine   source line, usually __LINE__
	 * @return the underlying stream, positioned after the header.
	 */
	// Used through InfoIn(fn) / WarningIn(fn). Writes the header, returns the stream.
	std::ostream& operator()(const char *funcName, const char *srcName, int srcLine = 0);
	std::ostream& operator()(const std::string &funcName, const char *srcName, int srcLine = 0);

	/**
	 * @brief Gives the underlying stream with no header written.
	 *        Use this to append to a message, or to log without source location.
	 */
	operator std::ostream&();
	std::ostream& operator()() { return operator std::ostream&(); }

protected:
	std::ostream *_streamPtr;
	std::string   _title;
	mesType       _mtype;
};

/**
 * @brief Streams a value into a logger without writing a header first.
 * @return the underlying stream, so further << calls chain as usual.
 */
template<class T>
std::ostream& operator<<(messageStream &ms, const T &v)
{
	std::ostream &os = static_cast<std::ostream&>(ms);
	os << v;
	return os;
}
/**
 * @brief Applies a stream manipulator (std::endl, nl, ...) to a logger.
 */
// For manipulators such as std::endl
inline std::ostream& operator<<(messageStream &ms, std::ostream& (*f)(std::ostream&))
{
	std::ostream &os = static_cast<std::ostream&>(ms);
	return f(os);
}


//==================================================================
//  Error : abnormal termination
//==================================================================
class Error : public std::exception, public messageStream
{
public:
	/**
	 * @brief Builds a fatal-error logger writing to std::cerr.
	 * @param title    label printed ahead of each message
	 * @param swAbort  if true, Abort() ends the process with EXIT_FAILURE
	 * @param swThrow  if true, Abort() throws *this instead of (or before)
	 *                 exiting, which is what the unit tests rely on
	 */
	Error(const std::string &title = "FATAL ERROR", bool swAbort = true, bool swThrow = false);
	virtual ~Error() throw() {}

	/**
	 * @brief Writes the error header and returns the stream to go on writing to.
	 *        Reaching this does not yet stop the run; call Abort() for that.
	 *        Normally used through the FatalErrorIn(fn) macro.
	 */
	std::ostream& operator()(const char *funcName, const char *srcName, int line = 0);
	std::ostream& operator()(const std::string &funcName, const char *srcName, int line = 0);

	/**
	 * @brief Ends the run: closes out the job information, then throws or exits
	 *        according to the flags given to the constructor. Normally does not
	 *        return to the caller.
	 */
	void Abort();
	/** @brief Exception message: the title this error was built with. */
	virtual const char* what() const throw() { return _title.c_str(); }

private:
	bool _abort;
	bool _throwExceptions;
};


//==================================================================
//  C_JobInfo : the job information written into JobInfo.txt
//==================================================================
class C_JobInfo : public messageStream
{
public:
	typedef std::map<std::string, std::string> tInfoList;

	/**
	 * @brief Builds an empty job record targeting std::cout.
	 *
	 * The result directory reads "not given" until SetResultDir() is called,
	 * and the elapsed-time clock is restarted by Start().
	 */
	C_JobInfo();
	virtual ~C_JobInfo() {}

	/**
	 * @brief Records one entry and echoes it to the job information stream.
	 * @param key       name the entry is stored under and looked up by
	 * @param str       the value, in the textual form it is stored and printed in
	 * @param mesTitle  label used in the printed line; the stored value keeps
	 *                  its own textual form (numbers go through operator<<)
	 */
	void Add(std::string key, std::string str, const std::string &mesTitle);
	void Add(std::string key, int         val, const std::string &mesTitle);
	void Add(std::string key, double      val, const std::string &mesTitle);

	/**
	 * @brief Opens the job record: writes the machine information and the start
	 *        date, time, working directory and result directory, and starts the
	 *        elapsed-time clock.
	 */
	void Start();
	/** @brief Writes the machine section (OS, host, CPU, thread count, memory, user, PID). */
	void MachineInfo();
	/**
	 * @brief Records the thread count OMP_NUM_THREADS asked for and the cells
	 *        per thread of the cap (0: no cap), for the machine section that
	 *        Start() writes (OMP_THREADS_REQUESTED, OMP_THREAD_CAP); the
	 *        count used is OMP_MAX_THREADS there.
	 */
	void SetThreads(int requested, int cellsPerThread)
	{ _threadsRequested = requested; _cellsPerThread = cellsPerThread; }
	/**
	 * @brief Closes the job record with the elapsed time, end date/time and the
	 *        reason it ended: End() normal, Exit() exit, Abort() abort.
	 *        Only the first of these has an effect; later calls do nothing.
	 */
	void End();
	void Exit();
	void Abort();
	/** @brief Drops all recorded entries. The already written output is untouched. */
	void Clear() { _infoList.clear(); }

	/**
	 * @brief Reads back a recorded entry by key.
	 * @return the stored text, or "Key Not Found" (with a warning) if the key
	 *         was never added.
	 */
	const std::string operator[](const std::string &key);

	// Directory the results of this job go into, as reported in the job record.
	void            SetResultDir(const FilePath &dir) { _resultDir = dir; }
	const FilePath& ResultDir() const { return _resultDir; }

	// True only between Start() and the first End()/Exit()/Abort().
	static bool constructed;

	/** @brief As End(), with the termination type given ("interrupted", ...). */
	void End(const std::string &terminationType);

private:

	tInfoList _infoList;
	int _threadsRequested = 0;	// 0: SetThreads() was not called
	int _cellsPerThread = 0;
	FilePath  _resultDir;
	double    _startSeconds;
};


//==================================================================
//  Global instances
//==================================================================
// The single logger instances the whole solver writes through.
// Info goes to std::cout, Warning and FatalError to std::cerr, and JobInfo
// collects what ends up in JobInfo.txt.
namespace Global {
	extern messageStream Info;
	extern messageStream Warning;
	extern Error         FatalError;
	extern C_JobInfo     JobInfo;
}

} } // namespace yuzen::util

// Message entry points. Each writes a header carrying fn plus this file and
// line, then returns a stream to append the message to. notImplemented(fn)
// additionally aborts the run.
#define InfoIn(fn)         yuzen::util::Global::Info(fn, __FILE__, __LINE__)
#define WarningIn(fn)      yuzen::util::Global::Warning(fn, __FILE__, __LINE__)
#define FatalErrorIn(fn)   yuzen::util::Global::FatalError(fn, __FILE__, __LINE__)
#define notImplemented(fn) FatalErrorIn(fn) << "Not implemented", yuzen::util::Global::FatalError.Abort()

#endif //_YUZEN_UTIL_INCLUDED_
