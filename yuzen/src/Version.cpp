/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Version.cpp : version numbers and the start-up logo
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <cstdio>

#include "Version.h"
#include "Yuzen.h"

#include <string>
#include <sstream>

std::string yuzen::C_Yuzen::Logo()
{
	std::ostringstream ss;
	ss<<"\n";
	ss<<"############################################\n";
	ss<<"#  __   __                                 #\n";
	ss<<"#  \\ \\ / /   _ _______ _ __                #\n";
	ss<<"#   \\ V / | | |_  / _ \\ '_ \\               #\n";
	ss<<"#    | || |_| |/ /  __/ | | |              #\n";
	ss<<"#    |_| \\__,_/___\\___|_| |_|              #\n";
	ss<<"#                                          #\n";
	ss<<"#  thin liquid film flow solver            #\n";
	ss<<"#  (lubrication approximation)             #\n";
	ss<<"############################################\n";
	return ss.str();
}

std::string yuzen::C_Yuzen::Version()
{
	std::ostringstream ss;
	ss<<YUZEN_VERSION_Major<<"."<<YUZEN_VERSION_Minor<<"."<<YUZEN_VERSION_Rev;
	return ss.str();
}

std::string yuzen::C_Yuzen::VersionSeries()
{
	// This is as far as the SolverName of an input file goes.
	// Raising PATCH leaves existing input files valid as they are.
	std::ostringstream ss;
	ss<<YUZEN_VERSION_Major<<"."<<YUZEN_VERSION_Minor;
	return ss.str();
}

std::string yuzen::C_Yuzen::VersionGrid()
{
	std::ostringstream ss;
	ss<<_GRID_VERSION_;
	return ss.str();
}
