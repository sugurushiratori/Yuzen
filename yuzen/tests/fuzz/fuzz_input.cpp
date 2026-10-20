/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  fuzz_input.cpp : libFuzzer harness of the input file -- what
*  `yuzen --check-input FILE` does, in process: ConfigBlock (the toml++
*  loader, or the .inp lexer and parser), every block and key of Input::Read
*  with its range checks and the unknown-key report, then CheckParameters.
*
*  Built twice: fuzz_input_toml (the bytes are a .toml) and fuzz_input_inp
*  (the bytes are a block-structured .inp), selected by YZFUZZ_EXT.
*
*  A refusal is an exit() with a non-zero code and is fine; a crash, a
*  sanitizer report, an uncaught exception, a hang, or an exit(0) is not.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "fuzz_common.h"

#include "Input.h"
#include "CheckParameters.h"

#ifndef YZFUZZ_EXT
#define YZFUZZ_EXT ".toml"
#endif

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
	static const std::string path = yzfuzz::TempPath("yz_fuzz_input", YZFUZZ_EXT);
	if(!yzfuzz::WriteFile(path, data, size)) return 0;

	yuzen::Input *inp = new yuzen::Input();
	const bool returned = YZFUZZ_TRAP_EXIT(
		if(inp->Read(path) == 0) (void)yuzen::CheckParameters::Check(*inp));
	if(!returned && yzfuzz::LastExitCode() == 0)
	{
		std::fprintf(stderr, "input refused with exit code 0\n");
		std::abort();
	}
	//  After an exit() the object is whole (the exit is never taken in the
	//  middle of changing a container), so it can be released.
	delete inp;
	return 0;
}
