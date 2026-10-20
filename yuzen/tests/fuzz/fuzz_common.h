/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  fuzz_common.h : shared plumbing of the libFuzzer harnesses
*
*  The input layer refuses a bad file by printing the reason and calling
*  exit(1) (Input::Error, Input::Read, Error::Abort). That is the intended
*  behavior, not a defect, but libFuzzer reports any exit() of the target as
*  a crash. The harnesses are linked with -Wl,--wrap=exit, so every call of
*  exit() from yuzencore lands in __wrap_exit below, which -- while a target
*  is running -- jumps back to the harness instead of ending the process.
*  Objects on the frames skipped are leaked, so run with -detect_leaks=0
*  and ASAN_OPTIONS=detect_leaks=0 (CMakeLists.txt here).
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#ifndef YUZEN_FUZZ_COMMON_H
#define YUZEN_FUZZ_COMMON_H

#include <csetjmp>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>

namespace yzfuzz {

inline std::jmp_buf &JumpBuffer() { static std::jmp_buf b; return b; }
inline bool &Trapping() { static bool t = false; return t; }
inline int &LastExitCode() { static int c = 0; return c; }

//  A file of this process's own, in $YUZEN_FUZZ_TMP (or $TMPDIR, or /tmp).
inline std::string TempPath(const std::string &stem, const std::string &ext)
{
	const char *d = std::getenv("YUZEN_FUZZ_TMP");
	if(!d) d = std::getenv("TMPDIR");
	if(!d) d = "/tmp";
	return std::string(d) + "/" + stem + "_" + std::to_string((long)getpid()) + ext;
}

inline bool WriteFile(const std::string &path, const uint8_t *data, size_t size)
{
	FILE *f = std::fopen(path.c_str(), "wb");
	if(!f) return false;
	if(size) std::fwrite(data, 1, size, f);
	std::fclose(f);
	return true;
}

}  // namespace yzfuzz

//  Run STMT; an exit() inside it comes back here, with the code in
//  yzfuzz::LastExitCode(). Evaluates to true when STMT returned normally.
#define YZFUZZ_TRAP_EXIT(STMT)                                    \
	([&]() -> bool {                                              \
		if(setjmp(yzfuzz::JumpBuffer()) == 0)                     \
		{                                                         \
			yzfuzz::Trapping() = true;                            \
			STMT;                                                 \
			yzfuzz::Trapping() = false;                           \
			return true;                                          \
		}                                                         \
		yzfuzz::Trapping() = false;                               \
		return false;                                             \
	}())

#endif
