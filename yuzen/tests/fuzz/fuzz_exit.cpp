/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  fuzz_exit.cpp : the exit() wrapper of the fuzz harnesses (see fuzz_common.h)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "fuzz_common.h"

extern "C" void __real_exit(int code);

extern "C" void __wrap_exit(int code)
{
	if(yzfuzz::Trapping())
	{
		yzfuzz::Trapping() = false;
		yzfuzz::LastExitCode() = code;
		std::longjmp(yzfuzz::JumpBuffer(), 1);
	}
	__real_exit(code);
}
