/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Input.cpp : reading every parameter of the input file and checking its range
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#include "Input.h"
#include "FiniteDifferences.h"
#include "Yuzen.h"
#include <cctype>
#include <cmath>
#include <string>

using namespace yuzen;

// --------------------------------
// Constructor
// --------------------------------
Input::Input()
{
}
// --------------------------------
// private functions
// --------------------------------
int Input::CheckBlockSize(Block *parent, string name)
{
	int size = parent->getBlockSize(name.c_str());
	return size;
}
template<typename Type> bool Input::CheckVectorSize(vector<Type> &vec, const size_t size_must_be)
{
	return (vec.size() == size_must_be);
}
bool Input::CheckRange(const int    val, const int    min, const    int max)
{
	return (min <= val && val <= max);
}
bool Input::CheckRange(const double val, const double min, const double max)
{
	return (min <= val && val <= max);
}
bool Input::CheckRange(vector<int>    &vec, const int    min, const int    max)
{
	for(size_t i = 0; i<vec.size(); i++)
	{
		if(!(min <= vec[i] && vec[i] <= max)) return false;
	}
	return true;
}
bool Input::CheckRange(vector<double> &vec, const double min, const double max)
{
	for(size_t i = 0; i<vec.size(); i++)
	{
		if(!(min <= vec[i] && vec[i] <= max)) return false;
	}
	return true;
}
//  Values are matched exactly, case included. A value that differs from a
//  candidate only in case is still refused, but the message names the
//  spelling that would have been accepted.
static std::string CaseHint(const string &str, const vector<string> &candidates)
{
	auto lower = [](string s) { for(char &c : s) c = (char)std::tolower((unsigned char)c); return s; };
	for(size_t i = 0; i < candidates.size(); i++)
		if(lower(str) == lower(candidates[i]))
			return " Did you mean `" + candidates[i] + "'? Values are case-sensitive.";
	return "";
}
bool Input::CheckCandidates(const string& str, vector<string> &candidates)
{
	_caseHint.clear();
	for(size_t i = 0; i<candidates.size(); i++)
	{
		if(str == candidates[i]) return true;
	}
	_caseHint = CaseHint(str, candidates);
	return false;
}
bool Input::CheckCandidates(const vector<string>& str, vector<string> &candidates)
{
	bool check;
	for(size_t i = 0; i<str.size(); i++)
	{
		check = false;
		for(size_t j = 0; j<candidates.size(); j++)
		{
			if(str[i] == candidates[j]) check = true;
		}
		if(check == false) return false;
	}
	return true;
}
void Input::Error(const string &mes)
{
	std::cerr<<"Error in the input: "<<mes<<_caseHint<<endl;
	exit(1);
}
// --------------------------------
// Public functions

//------------------------------------------------------------------
//  Reading the input. There is one function per top-level block.
//------------------------------------------------------------------
int Input::Read(const std::string &filename)
{
	try
	{
		root_block.reset(Block::newread(filename.c_str()));

		ReadRoot();
		ReadGrid();
		ReadBoundary();
		ReadPhysics();
		ReadSubstrate();
		ReadInitial();
		ReadFluid();
		ReadOutput();
		ReadNumerics();
		ReadRestart();

		//  Everything above asks the tree for the keys and blocks it knows,
		//  whether they are there or not, so what is left unasked is what this
		//  version does not know: a misspelled key, a key in the wrong block, a
		//  block with a misspelled name. Such an entry would be read by nobody,
		//  and the setting it was meant to change would keep its default in
		//  silence -- alpha for Alpha would leave an EdgeBump at Alpha = 1 --
		//  so it stops the run before anything is computed.
		std::vector<std::string> unknown;
		root_block->collectUnknown("", unknown);
		if(!unknown.empty())
		{
			std::cerr<<endl<<"Error in the input: "<<filename<<" contains "<<unknown.size()
			         <<(unknown.size() == 1 ? " entry" : " entries")
			         <<" that this version of Yuzen does not read, so whatever "
			         <<(unknown.size() == 1 ? "it was" : "they were")
			         <<" meant to set would keep its default:"<<endl;
			for(size_t i = 0; i < unknown.size(); ++i) std::cerr<<"  "<<unknown[i]<<endl;
			std::cerr<<"Correct the spelling or remove the entry. The manual's input"
			         <<" reference lists every key each block accepts."<<endl;
			exit(1);
		}
	}
	catch(std::runtime_error &e)
	{
		std::cerr<<"Error in the input file "<<filename<<":"<<endl;
		std::cerr<<e.what()<<endl; exit(1);
	}
	return 0;
}
int Input::Write(const std::string &filename)
{
	root_block->write(filename);
	return 0;
}


//------------------------------------------------------------------
//  Block [Root]
//------------------------------------------------------------------
void Input::ReadRoot()
{
	Block *sb = 0;
	int bnum = 0;
	vector<string> cand; cand.clear();
	(void)sb; (void)bnum;

// ####### Get block [root] #####
sb = root_block.get();

// Get Variable [SolverName]
sb->getin("SolverName",  SolverName);
//------------------------------------------------------------------
//  The SolverName values that are accepted.
//
//  The unit is MAJOR.MINOR. Raising PATCH still lets the input file through.
//
//  When MINOR is raised, **add the new series to this list without removing
//  the old one**, unless the interpretation of the input has changed. Removing
//  it rejects every input file written up to that point. Drop an old series
//  only when the interpretation really did change (the meaning of a key, or a
//  default), and record that in the revision history.
//------------------------------------------------------------------
cand.clear();
cand.push_back("Yuzen-" + C_Yuzen::VersionSeries());	// the current series
// cand.push_back("Yuzen-1.0");	// <- when MINOR is raised, keep the past series in this form
 if(CheckCandidates(SolverName,  cand) == false)
 {
	 std::string msg = "Invalid value. SolverName must be one of [";
	 for(size_t i = 0; i < cand.size(); ++i)
		 msg += (i ? ", `" : "`") + cand[i] + "'";
	 Error(msg + "]. ");
 }

// Get Variable [InputFormatVersion]
//  The version of the input format the file was written for. Optional, since
//  there has only ever been one; checked when present, so that a file written
//  for a later format is refused rather than half-understood.
{
	std::string ifv;
	sb->getin("InputFormatVersion", ifv, std::string("1.0"));
	if(ifv != "1.0")
		Error("Invalid value. InputFormatVersion must be `1.0', the only input format this version reads. ");
}

// Get Variable [LengthScaleType]
sb->getin("LengthScaleType",  LengthScaleType, "capillary");
//Set candidates.
cand.clear();
cand.push_back("custom");
cand.push_back("capillary");
 if(CheckCandidates(LengthScaleType,  cand) == false)
	 Error("Invalid value. LengthScaleType must be one of [`custom', `capillary']. ");

// Get Variable [LengthScale]
sb->getin("LengthScale",  LengthScale, 1.0);

// Get Variable [JobName]
sb->getin("JobName",  JobName);
}

//------------------------------------------------------------------
//  Block [Grid]
//------------------------------------------------------------------
void Input::ReadGrid()
{
	Block *sb = 0;
	int bnum = 0;
	vector<string> cand; cand.clear();
	(void)sb; (void)bnum;

	// ####### Get block [Grid] #####
	// Check BlockSize.
	bnum = root_block->getBlockSize("Grid");
	if(bnum != 1) Error("Block [Grid] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
	sb = root_block->getBlock("Grid");

	//------------------------------------------------------------------
	//  Upper bound on the number of grid points.
	//
	//  A per-axis bound of 5000 would make a 1D run of more than 5000 points
	//  impossible. That is not enough to measure parallel scaling (the
	//  reduction would only get 5000/512 = 9 blocks, so there would be work for
	//  just 9 threads). On the other hand the linear index of the grid is an
	//  int, so the total number of points does need a bound. The per-axis bound
	//  is therefore large and the total number of points is what guards us.
	//
	//  The internal grid adds two virtual cells per axis, so the number of
	//  points allocated is (Nx+4)(Ny+4)(Nz+4). A bound of 5e7 fits comfortably
	//  inside an int (2.1e9). In practice memory runs out first (5e7 points is
	//  tens of GB).
	//------------------------------------------------------------------
	const int kMaxN = 1000000;
	const long long kMaxPoints = 50000000LL;

	// Get Variable [Nx]
	sb->getin("Nx",  Grid.Nx);
	 if(CheckRange(Grid.Nx, 1, kMaxN) == false)
		 Error("Invalid value. Grid.Nx must be between [1:1000000]. ");

	// Get Variable [Ny]
	sb->getin("Ny",  Grid.Ny, 1);
	 if(CheckRange(Grid.Ny, 1, kMaxN) == false)
		 Error("Invalid value. Grid.Ny must be between [1:1000000]. ");

	// Get Variable [Nz]
	sb->getin("Nz",  Grid.Nz, 1);
	 if(CheckRange(Grid.Nz, 1, kMaxN) == false)
		 Error("Invalid value. Grid.Nz must be between [1:1000000]. ");

	 if((long long)(Grid.Nx + 4) * (long long)(Grid.Ny + 4)
	    * (long long)(Grid.Nz + 4) > kMaxPoints)
		 Error("Invalid value. (Grid.Nx+4)*(Grid.Ny+4)*(Grid.Nz+4) must not exceed 50000000. ");

	// Get Variable [Xmin]
	sb->getin("Xmin",  Grid.Xmin);

	// Get Variable [Xmax]
	sb->getin("Xmax",  Grid.Xmax);
	if(!(Grid.Xmax > Grid.Xmin))
		Error("Invalid value. Grid.Xmax must exceed Grid.Xmin. ");

	// Get Variable [Ymin]
	sb->getin("Ymin",  Grid.Ymin, 0);

	// Get Variable [Ymax]
	sb->getin("Ymax",  Grid.Ymax, 1);

	// Get Variable [Dimension]
	sb->getin("Dimension",  Grid.Dimension);
	//Set candidates.
	cand.clear();
	cand.push_back("xh");
	cand.push_back("xyh");
	 if(CheckCandidates(Grid.Dimension,  cand) == false)
		 Error("Invalid value. Grid.Dimension must be one of [`xh',`xyh']. ");

	// ####### END: block [Grid] #####
}

//------------------------------------------------------------------
//  Block [Boundary]
//------------------------------------------------------------------
void Input::ReadBoundary()
{
	Block *sb = 0;
	Block absent;	// stands in for an optional block the file leaves out
	int bnum = 0;
	vector<string> cand; cand.clear();
	(void)sb; (void)bnum;

	// ####### Get block [Boundary] #####
	// Check BlockSize.
	bnum = root_block->getBlockSize("Boundary");
	if(bnum != 1) Error("Block [Boundary] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
	sb = root_block->getBlock("Boundary");


		// ####### Get block [XmaxBC] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Boundary")->getBlockSize("XmaxBC");
		if(bnum != 1) Error("Block [XmaxBC] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Boundary")->getBlock("XmaxBC");

		// Get Variable [Type]
		sb->getin("Type",  Boundary.XmaxBC.Type, "Symmetry");
		//Set candidates.
		cand.clear();
		cand.push_back("Symmetry");
		cand.push_back("Periodic");
		cand.push_back("Pinning");
		cand.push_back("FreshInlet");
		cand.push_back("Outlet");
		 if(CheckCandidates(Boundary.XmaxBC.Type,  cand) == false)
			 Error("Invalid value. Boundary.XmaxBC.Type must be one of [`Symmetry',`Periodic',`Pinning',`FreshInlet',`Outlet']. ");

		// ####### END: block [XmaxBC] #####

		// ####### Get block [XminBC] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Boundary")->getBlockSize("XminBC");
		if(bnum != 1) Error("Block [XminBC] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Boundary")->getBlock("XminBC");

		// Get Variable [Type]
		sb->getin("Type",  Boundary.XminBC.Type, "Symmetry");
		//Set candidates.
		cand.clear();
		cand.push_back("Symmetry");
		cand.push_back("Periodic");
		cand.push_back("Pinning");
		cand.push_back("FreshInlet");
		cand.push_back("Outlet");
		 if(CheckCandidates(Boundary.XminBC.Type,  cand) == false)
			 Error("Invalid value. Boundary.XminBC.Type must be one of [`Symmetry',`Periodic',`Pinning',`FreshInlet',`Outlet']. ");

		// ####### END: block [XminBC] #####

		// ####### Get block [YmaxBC] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Boundary")->getBlockSize("YmaxBC");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [YmaxBC] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Boundary")->getBlock("YmaxBC") : &absent;

			// Get Variable [Type]
			sb->getin("Type",  Boundary.YmaxBC.Type, "Symmetry");
			//Set candidates.
			cand.clear();
			cand.push_back("Symmetry");
			cand.push_back("Periodic");
			cand.push_back("Pinning");
			cand.push_back("FreshInlet");
			cand.push_back("Outlet");
			 if(CheckCandidates(Boundary.YmaxBC.Type,  cand) == false)
				 Error("Invalid value. Boundary.YmaxBC.Type must be one of [`Symmetry',`Periodic',`Pinning',`FreshInlet',`Outlet']. ");

		}
		// ####### END: block [YmaxBC] #####

		// ####### Get block [YminBC] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Boundary")->getBlockSize("YminBC");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [YminBC] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Boundary")->getBlock("YminBC") : &absent;

			// Get Variable [Type]
			sb->getin("Type",  Boundary.YminBC.Type, "Symmetry");
			//Set candidates.
			cand.clear();
			cand.push_back("Symmetry");
			cand.push_back("Periodic");
			cand.push_back("Pinning");
			cand.push_back("FreshInlet");
			cand.push_back("Outlet");
			 if(CheckCandidates(Boundary.YminBC.Type,  cand) == false)
				 Error("Invalid value. Boundary.YminBC.Type must be one of [`Symmetry',`Periodic',`Pinning',`FreshInlet',`Outlet']. ");

		}
		// ####### END: block [YminBC] #####
	// ####### END: block [Boundary] #####
}

//------------------------------------------------------------------
//  Block [Physics]
//------------------------------------------------------------------
void Input::ReadPhysics()
{
	Block *sb = 0;
	Block absent;	// stands in for an optional block the file leaves out
	int bnum = 0;
	vector<string> cand; cand.clear();
	(void)sb; (void)bnum;

	// ####### Get block [Physics] #####
	// Check BlockSize.
	bnum = root_block->getBlockSize("Physics");
	if(bnum != 1) Error("Block [Physics] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
	sb = root_block->getBlock("Physics");


		// ####### Get block [ConcField] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Physics")->getBlockSize("ConcField");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [ConcField] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Physics")->getBlock("ConcField") : &absent;

			// Get Variable [EnableConvection]
			sb->getin("EnableConvection",  Physics.ConcField.EnableConvection, true);

			// Get Variable [EnableDiffusion]
			sb->getin("EnableDiffusion",  Physics.ConcField.EnableDiffusion, true);

			// Get Variable [Model]
			//  Required when the block is there. Without the block no
			//  concentration equation is built, which is what "None" says.
			if(bnum == 0)
				Physics.ConcField.Model = "None";
			else
				sb->getin("Model",  Physics.ConcField.Model);
			//Set candidates.
			cand.clear();
			cand.push_back("WM");
			cand.push_back("CZ");
			cand.push_back("None");
			 if(CheckCandidates(Physics.ConcField.Model,  cand) == false)
				 Error("Invalid value. Physics.ConcField.Model must be one of [`WM',`CZ',`None']. ");

		}
		// ####### END: block [ConcField] #####

		// ####### Get block [LaplacePressure] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Physics")->getBlockSize("LaplacePressure");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [LaplacePressure] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Physics")->getBlock("LaplacePressure") : &absent;

			// Get Variable [Enable]
			sb->getin("Enable",  Physics.LaplacePressure.Enable, true);

		}
		//  An absent Physics sub-block is a term that is off.
		if(bnum == 0)
			Physics.LaplacePressure.Enable = false;
		// ####### END: block [LaplacePressure] #####

		// ####### Get block [SolutoCapillarity] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Physics")->getBlockSize("SolutoCapillarity");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [SolutoCapillarity] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Physics")->getBlock("SolutoCapillarity") : &absent;

			// Get Variable [Enable]
			sb->getin("Enable",  Physics.SolutoCapillarity.Enable, true);

		}
		//  An absent Physics sub-block is a term that is off.
		if(bnum == 0)
			Physics.SolutoCapillarity.Enable = false;
		// ####### END: block [SolutoCapillarity] #####

		// ####### Get block [DisjoiningPressure] #####
		// The disjoining-pressure term, which handles wettability (the contact line)
		try
		{
			bnum = root_block->getBlock("Physics")->getBlockSize("DisjoiningPressure");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [DisjoiningPressure] may appear at most once (found " + std::to_string(bnum) + ").");
		if( bnum == 1)
		{
			Block *dp = root_block->getBlock("Physics")->getBlock("DisjoiningPressure");
			dp->getin("Enable",       Physics.DisjoiningPressure.Enable, 0);
			dp->getin("ContactAngle", Physics.DisjoiningPressure.ContactAngle, 0);
			//  (n, m) = (3, 2) is the default pair, as when the block is
			//  absent, so a block that enables the term without them passes
			//  n > m > 1.
			dp->getin("ExponentN",    Physics.DisjoiningPressure.ExponentN, 3.0);
			dp->getin("ExponentM",    Physics.DisjoiningPressure.ExponentM, 2.0);
		}
		else
		{
			Physics.DisjoiningPressure.Enable = false;
			Physics.DisjoiningPressure.ContactAngle = 0.0;
			Physics.DisjoiningPressure.ExponentN = 3.0;
			Physics.DisjoiningPressure.ExponentM = 2.0;
		}
		if(Physics.DisjoiningPressure.Enable)
		{
			if(CheckRange(Physics.DisjoiningPressure.ContactAngle, 0.0, 90.0) == false)
				Error("Invalid value. Physics.DisjoiningPressure.ContactAngle must be between [0:90] deg. ");
			if( !(Physics.DisjoiningPressure.ExponentN > Physics.DisjoiningPressure.ExponentM
			      && Physics.DisjoiningPressure.ExponentM > 1.0) )
				Error("Physics.DisjoiningPressure: ExponentN > ExponentM > 1 is required.");
		}

		// ####### Get block [VerticalGravity] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Physics")->getBlockSize("VerticalGravity");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [VerticalGravity] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Physics")->getBlock("VerticalGravity") : &absent;

			// Get Variable [Acceleration]
			sb->getin("Acceleration",  Physics.VerticalGravity.Acceleration, 9.80665);

			// Get Variable [Enable]
			sb->getin("Enable",  Physics.VerticalGravity.Enable, false);

		}
		// ####### END: block [VerticalGravity] #####

		// ####### Get block [HorizontalGravity] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Physics")->getBlockSize("HorizontalGravity");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [HorizontalGravity] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Physics")->getBlock("HorizontalGravity") : &absent;

			// Get Variable [Enable]
			sb->getin("Enable",  Physics.HorizontalGravity.Enable, false);

			// Get VecVariable [Acceleration].
			// Set Default Values.
			Physics.HorizontalGravity.Acceleration.push_back(0);
			Physics.HorizontalGravity.Acceleration.push_back(0);
			sb->getin("Acceleration", Physics.HorizontalGravity.Acceleration, Physics.HorizontalGravity.Acceleration);
			// Check the size of vector
			if(CheckVectorSize( Physics.HorizontalGravity.Acceleration, 2) == false )
			{ Error("Acceleration must be an array of exactly 2 numbers.");}

		}
		// ####### END: block [HorizontalGravity] #####

		// ####### Get block [Centrifugal] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Physics")->getBlockSize("Centrifugal");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [Centrifugal] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Physics")->getBlock("Centrifugal") : &absent;

			// Get Variable [RotationRate]
			sb->getin("RotationRate",  Physics.Centrifugal.RotationRate, 0);

			// Get Variable [Enable]
			sb->getin("Enable",  Physics.Centrifugal.Enable, false);

			// Get VecVariable [RotationAxis].
			// Set Default Values.
			Physics.Centrifugal.RotationAxis.push_back(0);
			Physics.Centrifugal.RotationAxis.push_back(0);
			sb->getin("RotationAxis", Physics.Centrifugal.RotationAxis, Physics.Centrifugal.RotationAxis);
			// Check the size of vector
			if(CheckVectorSize( Physics.Centrifugal.RotationAxis, 2) == false )
			{ Error("RotationAxis must be an array of exactly 2 numbers.");}

		}
		// ####### END: block [Centrifugal] #####

		// ####### Get block [Slip] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Physics")->getBlockSize("Slip");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [Slip] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Physics")->getBlock("Slip") : &absent;

			// Get Variable [Enable]
			sb->getin("Enable",  Physics.Slip.Enable, false);

			// Get Variable [Length]
			sb->getin("Length",  Physics.Slip.Length, 0.0);
		}
		// ####### END: block [Slip] #####

		// ####### Get block [Evaporation] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Physics")->getBlockSize("Evaporation");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [Evaporation] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Physics")->getBlock("Evaporation") : &absent;

			// Get Variable [Coef]
			sb->getin("Coef",  Physics.Evaporation.Coef, 1);

			// Get Variable [NegPowExponent]
			sb->getin("NegPowExponent",  Physics.Evaporation.NegPowExponent, 0.25);

			// Get Variable [NegPowRmax]
			sb->getin("NegPowRmax",  Physics.Evaporation.NegPowRmax, 1);

			// Get Variable [PowerN]
			sb->getin("PowerN",  Physics.Evaporation.PowerN, 0.5);

			// Get Variable [Enable]
			sb->getin("Enable",  Physics.Evaporation.Enable, true);

			// Get Variable [SpatialFunc]
			sb->getin("SpatialFunc",  Physics.Evaporation.SpatialFunc, "uniform");
			//Set candidates.
			cand.clear();
			cand.push_back("uniform");
			cand.push_back("negPow");
			 if(CheckCandidates(Physics.Evaporation.SpatialFunc,  cand) == false)
				 Error("Invalid value. Physics.Evaporation.SpatialFunc must be one of [`uniform',`negPow']. ");

			// Get Variable [ConcFunc]
			sb->getin("ConcFunc",  Physics.Evaporation.ConcFunc, "Raoult");
			//Set candidates.
			cand.clear();
			cand.push_back("Raoult");
			cand.push_back("PowerN");
			 if(CheckCandidates(Physics.Evaporation.ConcFunc,  cand) == false)
				 Error("Invalid value. Physics.Evaporation.ConcFunc must be one of [`Raoult',`PowerN']. ");

			// Get Variable [EndConcType]
			sb->getin("EndConcType",  Physics.Evaporation.EndConcType, "mean");
			//Set candidates.
			cand.clear();
			cand.push_back("mean");
			cand.push_back("max");
			cand.push_back("min");
			if(CheckCandidates(Physics.Evaporation.EndConcType,  cand) == false)
				 Error("Invalid value. Physics.Evaporation.EndConcType must be one of [`mean',`max',`min']. ");

			// Get Variable [EndConcValue]
			sb->getin("EndConcValue",  Physics.Evaporation.EndConcValue, 0.95);

			// Get Variable [ConcR]
			sb->getin("ConcR",  Physics.Evaporation.ConcR, 0.0);

		}
		//  An absent Physics sub-block is a term that is off.
		if(bnum == 0)
			Physics.Evaporation.Enable = false;
		// ####### END: block [Evaporation] #####

		// ####### Get block [Source] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Physics")->getBlockSize("Source");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [Source] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Physics")->getBlock("Source") : &absent;

			// Get Variable [Radius]
			sb->getin("Radius",  Physics.Source.Radius, 1);

			// Get Variable [InletRate]
			sb->getin("InletRate",  Physics.Source.InletRate, 0);

			// Get Variable [Conc]
			sb->getin("Conc",  Physics.Source.Conc, 0);

			// Get Variable [ConcType]: the concentration the supplied liquid
			// arrives at, that of the film where it lands ("local") or Conc
			// ("supplied")
			sb->getin("ConcType",  Physics.Source.ConcType, "local");
			cand.clear();
			cand.push_back("local");
			cand.push_back("supplied");
			if(CheckCandidates(Physics.Source.ConcType,  cand) == false)
				 Error("Invalid value. Physics.Source.ConcType must be one of [`local',`supplied']. ");

			// Get Variable [TimeConstant]
			sb->getin("TimeConstant",  Physics.Source.TimeConstant, 0);

			// Get Variable [ArmLength]
			sb->getin("ArmLength",  Physics.Source.ArmLength, 1);

			// Get Variable [RotationRate]
			sb->getin("RotationRate",  Physics.Source.RotationRate, 0);

			// Get Variable [StartAngle]
			sb->getin("StartAngle",  Physics.Source.StartAngle, 0);

			// Get Variable [RampTime]
			sb->getin("RampTime",  Physics.Source.RampTime, 0);

			// Get Variable [PlateauTime]
			sb->getin("PlateauTime",  Physics.Source.PlateauTime, 1e+010);

			// Get Variable [SubRegions]
			sb->getin("SubRegions",  Physics.Source.SubRegions, 1);
			if(Physics.Source.SubRegions < 1)
				Error("Invalid value. Physics.Source.SubRegions must be 1 or larger. ");

			// Get Variable [Enable]
			sb->getin("Enable",  Physics.Source.Enable, true);

			// Get VecVariable [ArmCenter].
			// Set Default Values.
			Physics.Source.ArmCenter.push_back(0);
			Physics.Source.ArmCenter.push_back(0);
			sb->getin("ArmCenter", Physics.Source.ArmCenter, Physics.Source.ArmCenter);
			// Check the size of vector
			if(CheckVectorSize( Physics.Source.ArmCenter, 2) == false )
			{ Error("ArmCenter must be an array of exactly 2 numbers.");}

		}
		//  An absent Physics sub-block is a term that is off.
		if(bnum == 0)
			Physics.Source.Enable = false;
		// ####### END: block [Source] #####
	// ####### END: block [Physics] #####
}

//------------------------------------------------------------------
//  Block [Substrate]
//------------------------------------------------------------------
void Input::ReadSubstrate()
{
	Block *sb = 0;
	Block absent;	// stands in for an optional block the file leaves out
	int bnum = 0;
	vector<string> cand; cand.clear();
	(void)sb; (void)bnum;

	// ####### Get block [Substrate] #####
	// Check BlockSize.
	try
	{
		bnum = root_block->getBlockSize("Substrate");
	}
	catch(const BlockNotFound&)
	{
		bnum = 0;
	}
	if(bnum > 1) Error("Block [Substrate] may appear at most once (found " + std::to_string(bnum) + ").");
	//  Absent, the block reads as an empty one: every key takes its default.
	{
		sb = (bnum == 1) ? root_block->getBlock("Substrate") : &absent;

		// Get Variable [TotalHeight]
		sb->getin("TotalHeight",  Substrate.TotalHeight, 1);


			// ####### Get block [RectBump] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Substrate")->getBlockSize("RectBump");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : RectBump.
			Substrate.RectBump.resize(bnum);	 //allocate memory.
			for(size_t i_001 = 0; i_001 < Substrate.RectBump.size(); i_001++)
			{
				try
				{
					sb = root_block->getBlock("Substrate")->getBlock("RectBump", i_001);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Substrate.RectBump[i_001].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Substrate.RectBump[i_001].Base, 0);

				// Get Variable [Height]
				sb->getin("Height",  Substrate.RectBump[i_001].Height, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Substrate.RectBump[i_001].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Substrate.RectBump[i_001].CreateType,  cand) == false)
					 Error("Invalid value. Substrate.RectBump block " + std::to_string(i_001 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Center].
				// Set Default Values.
				Substrate.RectBump[i_001].Center.push_back(0);
				Substrate.RectBump[i_001].Center.push_back(0);
				sb->getin("Center", Substrate.RectBump[i_001].Center, Substrate.RectBump[i_001].Center);
				// Check the size of vector
				if(CheckVectorSize( Substrate.RectBump[i_001].Center, 2) == false )
				{ Error("Center must be an array of exactly 2 numbers.");}

				// Get VecVariable [Width].
				// Set Default Values.
				Substrate.RectBump[i_001].Width.push_back(1);
				Substrate.RectBump[i_001].Width.push_back(1);
				sb->getin("Width", Substrate.RectBump[i_001].Width, Substrate.RectBump[i_001].Width);
				// Check the size of vector
				if(CheckVectorSize( Substrate.RectBump[i_001].Width, 2) == false )
				{ Error("Width must be an array of exactly 2 numbers.");}

				// Get VecVariable [SlopeWidth].
				// Set Default Values.
				Substrate.RectBump[i_001].SlopeWidth.push_back(0.05);
				Substrate.RectBump[i_001].SlopeWidth.push_back(0.05);
				sb->getin("SlopeWidth", Substrate.RectBump[i_001].SlopeWidth, Substrate.RectBump[i_001].SlopeWidth);
				// Check the size of vector
				if(CheckVectorSize( Substrate.RectBump[i_001].SlopeWidth, 2) == false )
				{ Error("SlopeWidth must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [RectBump] #####

			// ####### Get block [EllipseBump] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Substrate")->getBlockSize("EllipseBump");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : EllipseBump.
			Substrate.EllipseBump.resize(bnum);	 //allocate memory.
			for(size_t i_002 = 0; i_002 < Substrate.EllipseBump.size(); i_002++)
			{
				try
				{
					sb = root_block->getBlock("Substrate")->getBlock("EllipseBump", i_002);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Substrate.EllipseBump[i_002].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Substrate.EllipseBump[i_002].Base, 0);

				// Get Variable [Height]
				sb->getin("Height",  Substrate.EllipseBump[i_002].Height, 1);

				// Get Variable [SlopeWidth]
				sb->getin("SlopeWidth",  Substrate.EllipseBump[i_002].SlopeWidth, 0.05);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Substrate.EllipseBump[i_002].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Substrate.EllipseBump[i_002].CreateType,  cand) == false)
					 Error("Invalid value. Substrate.EllipseBump block " + std::to_string(i_002 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Center].
				// Set Default Values.
				Substrate.EllipseBump[i_002].Center.push_back(0);
				Substrate.EllipseBump[i_002].Center.push_back(0);
				sb->getin("Center", Substrate.EllipseBump[i_002].Center, Substrate.EllipseBump[i_002].Center);
				// Check the size of vector
				if(CheckVectorSize( Substrate.EllipseBump[i_002].Center, 2) == false )
				{ Error("Center must be an array of exactly 2 numbers.");}

				// Get VecVariable [Width].
				// Set Default Values.
				Substrate.EllipseBump[i_002].Width.push_back(1);
				Substrate.EllipseBump[i_002].Width.push_back(1);
				sb->getin("Width", Substrate.EllipseBump[i_002].Width, Substrate.EllipseBump[i_002].Width);
				// Check the size of vector
				if(CheckVectorSize( Substrate.EllipseBump[i_002].Width, 2) == false )
				{ Error("Width must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [EllipseBump] #####

			// ####### Get block [PeriodicWave] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Substrate")->getBlockSize("PeriodicWave");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : PeriodicWave.
			Substrate.PeriodicWave.resize(bnum);	 //allocate memory.
			for(size_t i_003 = 0; i_003 < Substrate.PeriodicWave.size(); i_003++)
			{
				try
				{
					sb = root_block->getBlock("Substrate")->getBlock("PeriodicWave", i_003);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Substrate.PeriodicWave[i_003].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Substrate.PeriodicWave[i_003].Base, 0);

				// Get Variable [Amplitude]
				sb->getin("Amplitude",  Substrate.PeriodicWave[i_003].Amplitude, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Substrate.PeriodicWave[i_003].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Substrate.PeriodicWave[i_003].CreateType,  cand) == false)
					 Error("Invalid value. Substrate.PeriodicWave block " + std::to_string(i_003 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Wavelength].
				// Set Default Values.
				Substrate.PeriodicWave[i_003].Wavelength.push_back(1);
				Substrate.PeriodicWave[i_003].Wavelength.push_back(1);
				sb->getin("Wavelength", Substrate.PeriodicWave[i_003].Wavelength, Substrate.PeriodicWave[i_003].Wavelength);
				// Check the size of vector
				if(CheckVectorSize( Substrate.PeriodicWave[i_003].Wavelength, 2) == false )
				{ Error("Wavelength must be an array of exactly 2 numbers.");}

				// Get VecVariable [Phase].
				// Set Default Values.
				Substrate.PeriodicWave[i_003].Phase.push_back(0);
				Substrate.PeriodicWave[i_003].Phase.push_back(0);
				sb->getin("Phase", Substrate.PeriodicWave[i_003].Phase, Substrate.PeriodicWave[i_003].Phase);
				// Check the size of vector
				if(CheckVectorSize( Substrate.PeriodicWave[i_003].Phase, 2) == false )
				{ Error("Phase must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [PeriodicWave] #####

			// ####### Get block [FlatPlane] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Substrate")->getBlockSize("FlatPlane");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : FlatPlane.
			Substrate.FlatPlane.resize(bnum);	 //allocate memory.
			for(size_t i_004 = 0; i_004 < Substrate.FlatPlane.size(); i_004++)
			{
				try
				{
					sb = root_block->getBlock("Substrate")->getBlock("FlatPlane", i_004);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Substrate.FlatPlane[i_004].Order, 0);

				// Get Variable [Height]
				sb->getin("Height",  Substrate.FlatPlane[i_004].Height, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Substrate.FlatPlane[i_004].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Substrate.FlatPlane[i_004].CreateType,  cand) == false)
					 Error("Invalid value. Substrate.FlatPlane block " + std::to_string(i_004 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Center].
				// Set Default Values.
				Substrate.FlatPlane[i_004].Center.push_back(0);
				Substrate.FlatPlane[i_004].Center.push_back(0);
				sb->getin("Center", Substrate.FlatPlane[i_004].Center, Substrate.FlatPlane[i_004].Center);
				// Check the size of vector
				if(CheckVectorSize( Substrate.FlatPlane[i_004].Center, 2) == false )
				{ Error("Center must be an array of exactly 2 numbers.");}

				// Get VecVariable [Angle].
				// Set Default Values.
				Substrate.FlatPlane[i_004].Angle.push_back(0);
				Substrate.FlatPlane[i_004].Angle.push_back(0);
				sb->getin("Angle", Substrate.FlatPlane[i_004].Angle, Substrate.FlatPlane[i_004].Angle);
				// Check the size of vector
				if(CheckVectorSize( Substrate.FlatPlane[i_004].Angle, 2) == false )
				{ Error("Angle must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [FlatPlane] #####

			// ####### Get block [MexicanHat] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Substrate")->getBlockSize("MexicanHat");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : MexicanHat.
			Substrate.MexicanHat.resize(bnum);	 //allocate memory.
			for(size_t i_005 = 0; i_005 < Substrate.MexicanHat.size(); i_005++)
			{
				try
				{
					sb = root_block->getBlock("Substrate")->getBlock("MexicanHat", i_005);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Substrate.MexicanHat[i_005].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Substrate.MexicanHat[i_005].Base, 0);

				// Get Variable [Height]
				sb->getin("Height",  Substrate.MexicanHat[i_005].Height, 1);

				// Get Variable [Beta]
				sb->getin("Beta",  Substrate.MexicanHat[i_005].Beta, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Substrate.MexicanHat[i_005].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Substrate.MexicanHat[i_005].CreateType,  cand) == false)
					 Error("Invalid value. Substrate.MexicanHat block " + std::to_string(i_005 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Center].
				// Set Default Values.
				Substrate.MexicanHat[i_005].Center.push_back(0);
				Substrate.MexicanHat[i_005].Center.push_back(0);
				sb->getin("Center", Substrate.MexicanHat[i_005].Center, Substrate.MexicanHat[i_005].Center);
				// Check the size of vector
				if(CheckVectorSize( Substrate.MexicanHat[i_005].Center, 2) == false )
				{ Error("Center must be an array of exactly 2 numbers.");}

				// Get VecVariable [Width].
				// Set Default Values.
				Substrate.MexicanHat[i_005].Width.push_back(1);
				Substrate.MexicanHat[i_005].Width.push_back(1);
				sb->getin("Width", Substrate.MexicanHat[i_005].Width, Substrate.MexicanHat[i_005].Width);
				// Check the size of vector
				if(CheckVectorSize( Substrate.MexicanHat[i_005].Width, 2) == false )
				{ Error("Width must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [MexicanHat] #####

			// ####### Get block [SmoothConvex] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Substrate")->getBlockSize("SmoothConvex");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : SmoothConvex.
			Substrate.SmoothConvex.resize(bnum);	 //allocate memory.
			for(size_t i_006 = 0; i_006 < Substrate.SmoothConvex.size(); i_006++)
			{
				try
				{
					sb = root_block->getBlock("Substrate")->getBlock("SmoothConvex", i_006);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Substrate.SmoothConvex[i_006].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Substrate.SmoothConvex[i_006].Base, 0);

				// Get Variable [Height]
				sb->getin("Height",  Substrate.SmoothConvex[i_006].Height, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Substrate.SmoothConvex[i_006].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Substrate.SmoothConvex[i_006].CreateType,  cand) == false)
					 Error("Invalid value. Substrate.SmoothConvex block " + std::to_string(i_006 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Center].
				// Set Default Values.
				Substrate.SmoothConvex[i_006].Center.push_back(0);
				Substrate.SmoothConvex[i_006].Center.push_back(0);
				sb->getin("Center", Substrate.SmoothConvex[i_006].Center, Substrate.SmoothConvex[i_006].Center);
				// Check the size of vector
				if(CheckVectorSize( Substrate.SmoothConvex[i_006].Center, 2) == false )
				{ Error("Center must be an array of exactly 2 numbers.");}

				// Get VecVariable [Width].
				// Set Default Values.
				Substrate.SmoothConvex[i_006].Width.push_back(1);
				Substrate.SmoothConvex[i_006].Width.push_back(1);
				sb->getin("Width", Substrate.SmoothConvex[i_006].Width, Substrate.SmoothConvex[i_006].Width);
				// Check the size of vector
				if(CheckVectorSize( Substrate.SmoothConvex[i_006].Width, 2) == false )
				{ Error("Width must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [SmoothConvex] #####
	}
	// ####### END: block [Substrate] #####
}

//------------------------------------------------------------------
//  Block [Initial]
//------------------------------------------------------------------
void Input::ReadInitial()
{
	Block *sb = 0;
	int bnum = 0;
	vector<string> cand; cand.clear();
	(void)sb; (void)bnum;

	// ####### Get block [Initial] #####
	// Check BlockSize.
	bnum = root_block->getBlockSize("Initial");
	if(bnum != 1) Error("Block [Initial] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
	sb = root_block->getBlock("Initial");


		// ####### Get block [Concentration] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Initial")->getBlockSize("Concentration");
		if(bnum != 1) Error("Block [Concentration] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Initial")->getBlock("Concentration");

		// Get Variable [Value]
		sb->getin("Value",  Initial.Concentration.Value);
		 if(CheckRange(Initial.Concentration.Value, 0.0, 1.0) == false)
			 Error("Invalid value. Initial.Concentration.Value must be between [0:1]. ");

		// Get Variable [Type]
		sb->getin("Type",  Initial.Concentration.Type, "uniform");
		//Set candidates.
		cand.clear();
		cand.push_back("uniform");
		 if(CheckCandidates(Initial.Concentration.Type,  cand) == false)
			 Error("Invalid value. Initial.Concentration.Type must be one of [`uniform']. ");

		// ####### END: block [Concentration] #####

		// ####### Get block [SurfaceShape] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Initial")->getBlockSize("SurfaceShape");
		if(bnum != 1) Error("Block [SurfaceShape] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Initial")->getBlock("SurfaceShape");

		// Get Variable [TotalHeight]
		sb->getin("TotalHeight",  Initial.SurfaceShape.TotalHeight, 1);

		// Get Variable [PrecursorRatio]
		sb->getin("PrecursorRatio",  Initial.SurfaceShape.PrecursorRatio, 0);


			// ####### Get block [RectBump] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlockSize("RectBump");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : RectBump.
			Initial.SurfaceShape.RectBump.resize(bnum);	 //allocate memory.
			for(size_t i_007 = 0; i_007 < Initial.SurfaceShape.RectBump.size(); i_007++)
			{
				try
				{
					sb = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlock("RectBump", i_007);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Initial.SurfaceShape.RectBump[i_007].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Initial.SurfaceShape.RectBump[i_007].Base, 0);

				// Get Variable [Height]
				sb->getin("Height",  Initial.SurfaceShape.RectBump[i_007].Height, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Initial.SurfaceShape.RectBump[i_007].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Initial.SurfaceShape.RectBump[i_007].CreateType,  cand) == false)
					 Error("Invalid value. Initial.SurfaceShape.RectBump block " + std::to_string(i_007 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Center].
				// Set Default Values.
				Initial.SurfaceShape.RectBump[i_007].Center.push_back(0);
				Initial.SurfaceShape.RectBump[i_007].Center.push_back(0);
				sb->getin("Center", Initial.SurfaceShape.RectBump[i_007].Center, Initial.SurfaceShape.RectBump[i_007].Center);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.RectBump[i_007].Center, 2) == false )
				{ Error("Center must be an array of exactly 2 numbers.");}

				// Get VecVariable [Width].
				// Set Default Values.
				Initial.SurfaceShape.RectBump[i_007].Width.push_back(1);
				Initial.SurfaceShape.RectBump[i_007].Width.push_back(1);
				sb->getin("Width", Initial.SurfaceShape.RectBump[i_007].Width, Initial.SurfaceShape.RectBump[i_007].Width);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.RectBump[i_007].Width, 2) == false )
				{ Error("Width must be an array of exactly 2 numbers.");}

				// Get VecVariable [SlopeWidth].
				// Set Default Values.
				Initial.SurfaceShape.RectBump[i_007].SlopeWidth.push_back(0.05);
				Initial.SurfaceShape.RectBump[i_007].SlopeWidth.push_back(0.05);
				sb->getin("SlopeWidth", Initial.SurfaceShape.RectBump[i_007].SlopeWidth, Initial.SurfaceShape.RectBump[i_007].SlopeWidth);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.RectBump[i_007].SlopeWidth, 2) == false )
				{ Error("SlopeWidth must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [RectBump] #####

			// ####### Get block [EllipseBump] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlockSize("EllipseBump");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : EllipseBump.
			Initial.SurfaceShape.EllipseBump.resize(bnum);	 //allocate memory.
			for(size_t i_008 = 0; i_008 < Initial.SurfaceShape.EllipseBump.size(); i_008++)
			{
				try
				{
					sb = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlock("EllipseBump", i_008);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Initial.SurfaceShape.EllipseBump[i_008].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Initial.SurfaceShape.EllipseBump[i_008].Base, 0);

				// Get Variable [Height]
				sb->getin("Height",  Initial.SurfaceShape.EllipseBump[i_008].Height, 1);

				// Get Variable [SlopeWidth]
				sb->getin("SlopeWidth",  Initial.SurfaceShape.EllipseBump[i_008].SlopeWidth, 0.05);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Initial.SurfaceShape.EllipseBump[i_008].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Initial.SurfaceShape.EllipseBump[i_008].CreateType,  cand) == false)
					 Error("Invalid value. Initial.SurfaceShape.EllipseBump block " + std::to_string(i_008 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Center].
				// Set Default Values.
				Initial.SurfaceShape.EllipseBump[i_008].Center.push_back(0);
				Initial.SurfaceShape.EllipseBump[i_008].Center.push_back(0);
				sb->getin("Center", Initial.SurfaceShape.EllipseBump[i_008].Center, Initial.SurfaceShape.EllipseBump[i_008].Center);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.EllipseBump[i_008].Center, 2) == false )
				{ Error("Center must be an array of exactly 2 numbers.");}

				// Get VecVariable [Width].
				// Set Default Values.
				Initial.SurfaceShape.EllipseBump[i_008].Width.push_back(1);
				Initial.SurfaceShape.EllipseBump[i_008].Width.push_back(1);
				sb->getin("Width", Initial.SurfaceShape.EllipseBump[i_008].Width, Initial.SurfaceShape.EllipseBump[i_008].Width);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.EllipseBump[i_008].Width, 2) == false )
				{ Error("Width must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [EllipseBump] #####

			// ####### Get block [PeriodicWave] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlockSize("PeriodicWave");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : PeriodicWave.
			Initial.SurfaceShape.PeriodicWave.resize(bnum);	 //allocate memory.
			for(size_t i_009 = 0; i_009 < Initial.SurfaceShape.PeriodicWave.size(); i_009++)
			{
				try
				{
					sb = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlock("PeriodicWave", i_009);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Initial.SurfaceShape.PeriodicWave[i_009].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Initial.SurfaceShape.PeriodicWave[i_009].Base, 0);

				// Get Variable [Amplitude]
				sb->getin("Amplitude",  Initial.SurfaceShape.PeriodicWave[i_009].Amplitude, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Initial.SurfaceShape.PeriodicWave[i_009].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Initial.SurfaceShape.PeriodicWave[i_009].CreateType,  cand) == false)
					 Error("Invalid value. Initial.SurfaceShape.PeriodicWave block " + std::to_string(i_009 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Wavelength].
				// Set Default Values.
				Initial.SurfaceShape.PeriodicWave[i_009].Wavelength.push_back(1);
				Initial.SurfaceShape.PeriodicWave[i_009].Wavelength.push_back(1);
				sb->getin("Wavelength", Initial.SurfaceShape.PeriodicWave[i_009].Wavelength, Initial.SurfaceShape.PeriodicWave[i_009].Wavelength);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.PeriodicWave[i_009].Wavelength, 2) == false )
				{ Error("Wavelength must be an array of exactly 2 numbers.");}

				// Get VecVariable [Phase].
				// Set Default Values.
				Initial.SurfaceShape.PeriodicWave[i_009].Phase.push_back(0);
				Initial.SurfaceShape.PeriodicWave[i_009].Phase.push_back(0);
				sb->getin("Phase", Initial.SurfaceShape.PeriodicWave[i_009].Phase, Initial.SurfaceShape.PeriodicWave[i_009].Phase);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.PeriodicWave[i_009].Phase, 2) == false )
				{ Error("Phase must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [PeriodicWave] #####

			// ####### Get block [FlatPlane] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlockSize("FlatPlane");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : FlatPlane.
			Initial.SurfaceShape.FlatPlane.resize(bnum);	 //allocate memory.
			for(size_t i_010 = 0; i_010 < Initial.SurfaceShape.FlatPlane.size(); i_010++)
			{
				try
				{
					sb = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlock("FlatPlane", i_010);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Initial.SurfaceShape.FlatPlane[i_010].Order, 0);

				// Get Variable [Height]
				sb->getin("Height",  Initial.SurfaceShape.FlatPlane[i_010].Height, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Initial.SurfaceShape.FlatPlane[i_010].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Initial.SurfaceShape.FlatPlane[i_010].CreateType,  cand) == false)
					 Error("Invalid value. Initial.SurfaceShape.FlatPlane block " + std::to_string(i_010 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Center].
				// Set Default Values.
				Initial.SurfaceShape.FlatPlane[i_010].Center.push_back(0);
				Initial.SurfaceShape.FlatPlane[i_010].Center.push_back(0);
				sb->getin("Center", Initial.SurfaceShape.FlatPlane[i_010].Center, Initial.SurfaceShape.FlatPlane[i_010].Center);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.FlatPlane[i_010].Center, 2) == false )
				{ Error("Center must be an array of exactly 2 numbers.");}

				// Get VecVariable [Angle].
				// Set Default Values.
				Initial.SurfaceShape.FlatPlane[i_010].Angle.push_back(0);
				Initial.SurfaceShape.FlatPlane[i_010].Angle.push_back(0);
				sb->getin("Angle", Initial.SurfaceShape.FlatPlane[i_010].Angle, Initial.SurfaceShape.FlatPlane[i_010].Angle);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.FlatPlane[i_010].Angle, 2) == false )
				{ Error("Angle must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [FlatPlane] #####

			// ####### Get block [MexicanHat] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlockSize("MexicanHat");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : MexicanHat.
			Initial.SurfaceShape.MexicanHat.resize(bnum);	 //allocate memory.
			for(size_t i_011 = 0; i_011 < Initial.SurfaceShape.MexicanHat.size(); i_011++)
			{
				try
				{
					sb = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlock("MexicanHat", i_011);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Initial.SurfaceShape.MexicanHat[i_011].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Initial.SurfaceShape.MexicanHat[i_011].Base, 0);

				// Get Variable [Height]
				sb->getin("Height",  Initial.SurfaceShape.MexicanHat[i_011].Height, 1);

				// Get Variable [Beta]
				sb->getin("Beta",  Initial.SurfaceShape.MexicanHat[i_011].Beta, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Initial.SurfaceShape.MexicanHat[i_011].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Initial.SurfaceShape.MexicanHat[i_011].CreateType,  cand) == false)
					 Error("Invalid value. Initial.SurfaceShape.MexicanHat block " + std::to_string(i_011 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Center].
				// Set Default Values.
				Initial.SurfaceShape.MexicanHat[i_011].Center.push_back(0);
				Initial.SurfaceShape.MexicanHat[i_011].Center.push_back(0);
				sb->getin("Center", Initial.SurfaceShape.MexicanHat[i_011].Center, Initial.SurfaceShape.MexicanHat[i_011].Center);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.MexicanHat[i_011].Center, 2) == false )
				{ Error("Center must be an array of exactly 2 numbers.");}

				// Get VecVariable [Width].
				// Set Default Values.
				Initial.SurfaceShape.MexicanHat[i_011].Width.push_back(1);
				Initial.SurfaceShape.MexicanHat[i_011].Width.push_back(1);
				sb->getin("Width", Initial.SurfaceShape.MexicanHat[i_011].Width, Initial.SurfaceShape.MexicanHat[i_011].Width);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.MexicanHat[i_011].Width, 2) == false )
				{ Error("Width must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [MexicanHat] #####

			// ####### Get block [SmoothConvex] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlockSize("SmoothConvex");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : SmoothConvex.
			Initial.SurfaceShape.SmoothConvex.resize(bnum);	 //allocate memory.
			for(size_t i_012 = 0; i_012 < Initial.SurfaceShape.SmoothConvex.size(); i_012++)
			{
				try
				{
					sb = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlock("SmoothConvex", i_012);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Initial.SurfaceShape.SmoothConvex[i_012].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Initial.SurfaceShape.SmoothConvex[i_012].Base, 0);

				// Get Variable [Height]
				sb->getin("Height",  Initial.SurfaceShape.SmoothConvex[i_012].Height, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Initial.SurfaceShape.SmoothConvex[i_012].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Initial.SurfaceShape.SmoothConvex[i_012].CreateType,  cand) == false)
					 Error("Invalid value. Initial.SurfaceShape.SmoothConvex block " + std::to_string(i_012 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get VecVariable [Center].
				// Set Default Values.
				Initial.SurfaceShape.SmoothConvex[i_012].Center.push_back(0);
				Initial.SurfaceShape.SmoothConvex[i_012].Center.push_back(0);
				sb->getin("Center", Initial.SurfaceShape.SmoothConvex[i_012].Center, Initial.SurfaceShape.SmoothConvex[i_012].Center);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.SmoothConvex[i_012].Center, 2) == false )
				{ Error("Center must be an array of exactly 2 numbers.");}

				// Get VecVariable [Width].
				// Set Default Values.
				Initial.SurfaceShape.SmoothConvex[i_012].Width.push_back(1);
				Initial.SurfaceShape.SmoothConvex[i_012].Width.push_back(1);
				sb->getin("Width", Initial.SurfaceShape.SmoothConvex[i_012].Width, Initial.SurfaceShape.SmoothConvex[i_012].Width);
				// Check the size of vector
				if(CheckVectorSize( Initial.SurfaceShape.SmoothConvex[i_012].Width, 2) == false )
				{ Error("Width must be an array of exactly 2 numbers.");}

			}
			// ####### END: block [SmoothConvex] #####

			// ####### Get block [EdgeBump] #####
			// Check BlockSize.
			try
			{
				bnum = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlockSize("EdgeBump");
			}
			catch(const BlockNotFound&)
			{
				bnum = 0;
			}
			//------- Loop for block : EdgeBump.
			Initial.SurfaceShape.EdgeBump.resize(bnum);	 //allocate memory.
			for(size_t i_013 = 0; i_013 < Initial.SurfaceShape.EdgeBump.size(); i_013++)
			{
				try
				{
					sb = root_block->getBlock("Initial")->getBlock("SurfaceShape")->getBlock("EdgeBump", i_013);
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}

				// Get Variable [Order]
				sb->getin("Order",  Initial.SurfaceShape.EdgeBump[i_013].Order, 0);

				// Get Variable [Base]
				sb->getin("Base",  Initial.SurfaceShape.EdgeBump[i_013].Base, 0);

				// Get Variable [Height]
				sb->getin("Height",  Initial.SurfaceShape.EdgeBump[i_013].Height, 1);

				// Get Variable [Alpha]
				sb->getin("Alpha",  Initial.SurfaceShape.EdgeBump[i_013].alpha, 1);

				// Get Variable [InnerLength]
				sb->getin("InnerLength",  Initial.SurfaceShape.EdgeBump[i_013].Lin, 1);

				// Get Variable [OuterLength]
				sb->getin("OuterLength",  Initial.SurfaceShape.EdgeBump[i_013].Lout, 1);

				// Get Variable [CreateType]
				sb->getin("CreateType",  Initial.SurfaceShape.EdgeBump[i_013].CreateType, "add");
				//Set candidates.
				cand.clear();
				cand.push_back("add");
				cand.push_back("remove");
				cand.push_back("stack");
				 if(CheckCandidates(Initial.SurfaceShape.EdgeBump[i_013].CreateType,  cand) == false)
					 Error("Invalid value. Initial.SurfaceShape.EdgeBump block " + std::to_string(i_013 + 1) + ": CreateType must be one of [`add',`remove',`stack']. ");

				// Get Variable [Location]
				sb->getin("Location",  Initial.SurfaceShape.EdgeBump[i_013].Location, "Xmax");
				//Set candidates.
				cand.clear();
				cand.push_back("Xmin");
				cand.push_back("Xmax");
				cand.push_back("Ymin");
				cand.push_back("Ymax");
				 if(CheckCandidates(Initial.SurfaceShape.EdgeBump[i_013].Location,  cand) == false)
					 Error("Invalid value. Initial.SurfaceShape.EdgeBump block " + std::to_string(i_013 + 1) + ": Location must be one of [`Xmin',`Xmax',`Ymin',`Ymax']. ");

			}
			// ####### END: block [EdgeBump] #####
		// ####### END: block [SurfaceShape] #####
	// ####### END: block [Initial] #####
}

//------------------------------------------------------------------
//  Block [Fluid]
//------------------------------------------------------------------
void Input::ReadFluid()
{
	Block *sb = 0;
	int bnum = 0;
	vector<string> cand; cand.clear();
	(void)sb; (void)bnum;

	// ####### Get block [Fluid] #####
	// Check BlockSize.
	bnum = root_block->getBlockSize("Fluid");
	if(bnum != 1) Error("Block [Fluid] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
	sb = root_block->getBlock("Fluid");

	// Get Variable [Name]
	sb->getin("Name",  Fluid.Name, "CR-EX");


		// ####### Get block [Density] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Fluid")->getBlockSize("Density");
		if(bnum != 1) Error("Block [Density] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Fluid")->getBlock("Density");

		// Get Variable [Value]
		sb->getin("Value",  Fluid.Density.Value, 1070);

		// ####### END: block [Density] #####

		// ####### Get block [SurfaceTension] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Fluid")->getBlockSize("SurfaceTension");
		if(bnum != 1) Error("Block [SurfaceTension] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Fluid")->getBlock("SurfaceTension");

		// Get Variable [BaseValue]
		sb->getin("BaseValue",  Fluid.SurfaceTension.BaseValue, 0.028);

		// Get Variable [BaseConc]
		sb->getin("BaseConc",  Fluid.SurfaceTension.BaseConc, 0.5);

		// Get Variable [Coef]
		sb->getin("Coef",  Fluid.SurfaceTension.Coef, 0.014);

		// ####### END: block [SurfaceTension] #####

		// ####### Get block [Viscosity] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Fluid")->getBlockSize("Viscosity");
		if(bnum != 1) Error("Block [Viscosity] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Fluid")->getBlock("Viscosity");

		// Get Variable [BaseValue]
		sb->getin("BaseValue",  Fluid.Viscosity.BaseValue, 0.115);

		// Get Variable [BaseConc]
		sb->getin("BaseConc",  Fluid.Viscosity.BaseConc, 0.5);

		// Get Variable [Exponent]
		sb->getin("Exponent",  Fluid.Viscosity.Exponent, 8);

		// ####### END: block [Viscosity] #####

		// ####### Get block [Diffusivity] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Fluid")->getBlockSize("Diffusivity");
		if(bnum != 1) Error("Block [Diffusivity] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Fluid")->getBlock("Diffusivity");

		// Get Variable [BaseValue]
		sb->getin("BaseValue",  Fluid.Diffusivity.BaseValue, 9.78e-011);

		// Get Variable [BaseConc]
		sb->getin("BaseConc",  Fluid.Diffusivity.BaseConc, 0.5);

		// Get Variable [Exponent]
		sb->getin("Exponent",  Fluid.Diffusivity.Exponent, 6.94);

		// ####### END: block [Diffusivity] #####

		// ####### Get block [MolarMass] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Fluid")->getBlockSize("MolarMass");
		if(bnum != 1) Error("Block [MolarMass] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Fluid")->getBlock("MolarMass");

		// Get Variable [Solvent]
		sb->getin("Solvent",  Fluid.MolarMass.Solvent, 106);

		// Get Variable [Resin]
		sb->getin("Resin",  Fluid.MolarMass.Resin, 2200);

		// ####### END: block [MolarMass] #####
	// ####### END: block [Fluid] #####
}

//------------------------------------------------------------------
//  Block [Output]
//------------------------------------------------------------------
void Input::ReadOutput()
{
	Block *sb = 0;
	Block absent;	// stands in for an optional block the file leaves out
	int bnum = 0;
	vector<string> cand; cand.clear();
	(void)sb; (void)bnum;

	// ####### Get block [Output] #####
	// Check BlockSize.
	bnum = root_block->getBlockSize("Output");
	if(bnum != 1) Error("Block [Output] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
	sb = root_block->getBlock("Output");

	// Get Variable [Overwrite]
	sb->getin("Overwrite",  Output.Overwrite, false);

	// Get Variable [Dimensioned]
	sb->getin("Dimensioned",  Output.Dimensioned, true);


		// ####### Get block [FieldData] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Output")->getBlockSize("FieldData");
		if(bnum != 1) Error("Block [FieldData] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Output")->getBlock("FieldData");

		// Get Variable [TimeInterval]; whether it was given at all, since
		// it cannot be given with Times or a LogSchedule (CheckParameters)
		sb->getin("TimeInterval",  Output.FieldData.TimeInterval, 0.01);
		Output.FieldData.HasTimeInterval = sb->has("TimeInterval");

		// Get Variable [NeutralFile]
		sb->getin("NeutralFile",  Output.FieldData.NeutralFile, true);

		// Get Variable [VTKFile]
		sb->getin("VTKFile",  Output.FieldData.VTKFile, true);

		// Get Variable [WritePotential]
		sb->getin("WritePotential",  Output.FieldData.WritePotential, false);

		// Output at listed times: [Times] (and [ExactTimes]), or the block
		// [LogSchedule]; either replaces TimeInterval. CheckParameters
		// refuses both with it, and checks the times against EndTime.
		Output.FieldData.Times.assign(0, 0.0);
		sb->getin("Times", Output.FieldData.Times, Output.FieldData.Times);
		if(sb->has("Times") && Output.FieldData.Times.empty())
			Error("Output.FieldData.Times is empty: list at least one time, or leave the key out.");
		sb->getin("ExactTimes", Output.FieldData.ExactTimes, false);

			// ####### Get block [LogSchedule] (optional) #####
			try { bnum = sb->getBlockSize("LogSchedule"); } catch(const BlockNotFound&) { bnum = 0; }
			if(bnum > 1) Error("Block [Output.FieldData.LogSchedule] may appear only once (found " + std::to_string(bnum) + ").");
			Output.FieldData.LogSchedule.Present = (bnum == 1);
			if(bnum == 1)
			{
				Block *lb = sb->getBlock("LogSchedule");
				lb->getin("Start", Output.FieldData.LogSchedule.Start);
				lb->getin("End", Output.FieldData.LogSchedule.End);
				lb->getin("PerDecade", Output.FieldData.LogSchedule.PerDecade);
				lb->getin("ExactTimes", Output.FieldData.LogSchedule.ExactTimes, false);
				if(!(Output.FieldData.LogSchedule.Start > 0.0))
					Error("Output.FieldData.LogSchedule.Start must be positive (given " + std::to_string(Output.FieldData.LogSchedule.Start) + ").");
				if(!(Output.FieldData.LogSchedule.End > Output.FieldData.LogSchedule.Start))
					Error("Output.FieldData.LogSchedule.End must be after Start.");
				if(Output.FieldData.LogSchedule.PerDecade < 1)
					Error("Output.FieldData.LogSchedule.PerDecade must be 1 or more (given " + std::to_string(Output.FieldData.LogSchedule.PerDecade) + ").");
				if(!Output.FieldData.Times.empty())
					Error("Output.FieldData: give Times or a LogSchedule, not both.");
				//  Start 10^(i / PerDecade), i = 0, 1, ...: each from i, not by
				//  repeated multiplication, so that the times do not drift;
				//  End itself is kept when it falls on one to 1e-9
				const double decades = std::log10(Output.FieldData.LogSchedule.End / Output.FieldData.LogSchedule.Start);
				const long long last = (long long)std::floor(decades * Output.FieldData.LogSchedule.PerDecade + 1e-9);
				for(long long i = 0; i <= last; i++)
					Output.FieldData.Times.push_back(Output.FieldData.LogSchedule.Start
					                                 * std::pow(10.0, (double)i / Output.FieldData.LogSchedule.PerDecade));
				Output.FieldData.ExactTimes = Output.FieldData.LogSchedule.ExactTimes;
			}
			else
			{
				Output.FieldData.LogSchedule.Start = 0;
				Output.FieldData.LogSchedule.End = 0;
				Output.FieldData.LogSchedule.PerDecade = 0;
				Output.FieldData.LogSchedule.ExactTimes = false;
			}
			// ####### END: block [LogSchedule] #####

		// ####### END: block [FieldData] #####

		// ####### Get block [GeneralHistory] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Output")->getBlockSize("GeneralHistory");
		if(bnum != 1) Error("Block [GeneralHistory] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Output")->getBlock("GeneralHistory");

		// Get Variable [TimeInterval]
		sb->getin("TimeInterval",  Output.GeneralHistory.TimeInterval, 0.1);

		// ####### END: block [GeneralHistory] #####

		// ####### Get block [CalcInfo] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Output")->getBlockSize("CalcInfo");
		if(bnum != 1) Error("Block [CalcInfo] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Output")->getBlock("CalcInfo");

		// Get Variable [TimeInterval]
		sb->getin("TimeInterval",  Output.CalcInfo.TimeInterval, 0.1);

		// ####### END: block [CalcInfo] #####

		// ####### Get block [ProbeList] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Output")->getBlockSize("ProbeList");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [ProbeList] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Output")->getBlock("ProbeList") : &absent;

			// Get Variable [TimeInterval]
			sb->getin("TimeInterval",  Output.ProbeList.TimeInterval, 0.1);


				// ####### Get block [Probe] #####
				// Check BlockSize.
				try
				{
					bnum = root_block->getBlock("Output")->getBlock("ProbeList")->getBlockSize("Probe");
				}
				catch(const BlockNotFound&)
				{
					bnum = 0;
				}
				//------- Loop for block : Probe.
				Output.ProbeList.Probe.resize(bnum);	 //allocate memory.
				for(size_t i_014 = 0; i_014 < Output.ProbeList.Probe.size(); i_014++)
				{
					try
					{
						sb = root_block->getBlock("Output")->getBlock("ProbeList")->getBlock("Probe", i_014);
					}
					catch(const BlockNotFound&)
					{
						bnum = 0;
					}

					// Get Variable [X]
					sb->getin("X",  Output.ProbeList.Probe[i_014].px);

					// Get Variable [Y]
					sb->getin("Y",  Output.ProbeList.Probe[i_014].py, 0);

					// Get Variable [Z]
					sb->getin("Z",  Output.ProbeList.Probe[i_014].pz, 0);

				}
				// ####### END: block [Probe] #####
		}
		// ####### END: block [ProbeList] #####

		// ####### Get block [Energy] #####
		// Check BlockSize.
		try
		{
			bnum = root_block->getBlock("Output")->getBlockSize("Energy");
		}
		catch(const BlockNotFound&)
		{
			bnum = 0;
		}
		if(bnum > 1) Error("Block [Energy] may appear at most once (found " + std::to_string(bnum) + ").");
		//  Absent, the block reads as an empty one: every key takes its default.
		{
			sb = (bnum == 1) ? root_block->getBlock("Output")->getBlock("Energy") : &absent;

			// Get Variable [Enable]
			sb->getin("Enable",  Output.Energy.Enable, false);
		}
		// ####### END: block [Energy] #####
	// ####### END: block [Output] #####
}

//------------------------------------------------------------------
//  Block [Numerics]
//------------------------------------------------------------------
void Input::ReadNumerics()
{
	Block *sb = 0;
	int bnum = 0;
	vector<string> cand; cand.clear();
	(void)sb; (void)bnum;

	// ####### Get block [Numerics] #####
	// Check BlockSize.
	bnum = root_block->getBlockSize("Numerics");
	if(bnum != 1) Error("Block [Numerics] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
	sb = root_block->getBlock("Numerics");

	// Get Variable [MaxSteps]
	//  0, the default, sets no limit: the run is bounded by EndTime alone.
	sb->getin("MaxSteps",  Numerics.MaxSteps, 0);
	if(Numerics.MaxSteps < 0)
		Error("Invalid value. Numerics.MaxSteps must be 0 (no limit) or larger.");

	// Get Variable [InitialDt]
	sb->getin("InitialDt",  Numerics.InitialDt);

	// Get Variable [EndTime]
	sb->getin("EndTime",  Numerics.EndTime, 0);

	//------------------------------------------------------------------
	//  The steady-state test.
	//
	//  It looks at the rate of change per unit time, which does not depend on
	//  the time step. Looking at the change per step instead would silently
	//  declare a run whose time step has collapsed to be "steady".
	//  SteadyStateCriterion <= 0 switches the steady-state test off.
	//------------------------------------------------------------------
	sb->getin("SteadyStateCriterion", Numerics.SteadyStateCriterion, 1.0e-3);
	sb->getin("SteadyStateSteps",     Numerics.SteadyStateSteps, 10);
	if(!CheckRange(Numerics.SteadyStateSteps, 1, 1000000))
	{
		Error("Numerics.SteadyStateSteps must be in [1:1000000].");
	}


		// ####### Get block [TimeStepControl] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Numerics")->getBlockSize("TimeStepControl");
		if(bnum != 1) Error("Block [TimeStepControl] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Numerics")->getBlock("TimeStepControl");

		// Get Variable [KeepSteps]
		sb->getin("KeepSteps",  Numerics.TimeStepControl.KeepSteps, 50);
		if(Numerics.TimeStepControl.KeepSteps < 1)
			Error("Invalid value. Numerics.TimeStepControl.KeepSteps must be 1 or larger. ");

		// Get Variable [MaxDtFactor]
		sb->getin("MaxDtFactor",  Numerics.TimeStepControl.MaxDtFactor, 10000);

		//------------------------------------------------------------------
		//  Lower bound of the weight used when the time-step test measures the
		//  error.
		//
		//  The error is made relative by dividing by the local film thickness.
		//  The thin regions (the precursor film) then weigh tens of times more
		//  heavily than the bulk, so ScaleFloor puts in the lower bound "where
		//  the film is thinner than this fraction of the maximum, do not judge
		//  it any more strictly than that".
		//    0 (the default) : divide by the local film thickness
		//    1               : normalize by the global maximum alone
		//
		//  Note: this is **not** what holds the time step down in runs with a
		//     contact line; the Newton convergence test is (see the comment
		//     on NewtonEpsilonRel). Moving this from 0 to 1 does not shift the
		//     time step by even one order of magnitude there. It may matter
		//     in other runs, such as ones drying out.
		//------------------------------------------------------------------
		sb->getin("ScaleFloor", Numerics.TimeStepControl.ScaleFloor, 0.0);
		if(!CheckRange(Numerics.TimeStepControl.ScaleFloor, 0.0, 1.0))
		{
			Error("Numerics.TimeStepControl.ScaleFloor must be in [0:1].");
		}

		// Get Variable [MinDt]
		sb->getin("MinDt",  Numerics.TimeStepControl.MinDt, 1e-015);

		// Get Variable [Criterion]
		sb->getin("Criterion",  Numerics.TimeStepControl.Criterion, 0.01);

		// Get Variable [ChangeFactor]
		sb->getin("ChangeFactor",  Numerics.TimeStepControl.ChangeFactor, 0.5);
		//  The step is multiplied by it to shrink and divided by it to grow, so
		//  outside (0, 1) the two exchange roles or the step becomes zero.
		if(!(Numerics.TimeStepControl.ChangeFactor > 0.0 && Numerics.TimeStepControl.ChangeFactor < 1.0))
			Error("Invalid value. Numerics.TimeStepControl.ChangeFactor must be between (0:1), both ends excluded. ");

		// Get Variable [Enable]
		sb->getin("Enable",  Numerics.TimeStepControl.Enable, false);

		// ####### END: block [TimeStepControl] #####

		// ####### Get block [HeightEquation] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Numerics")->getBlockSize("HeightEquation");
		if(bnum != 1) Error("Block [HeightEquation] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Numerics")->getBlock("HeightEquation");

		// Get Variable [NewtonMaxIterations]
		sb->getin("NewtonMaxIterations",  Numerics.HeightEquation.NewtonMaxIterations, 30);

		// Get Variable [NewtonMinIterations]
		// The default is 1 (leave on the convergence test alone). A larger
		//          value makes every step take at least that many iterations.
		sb->getin("NewtonMinIterations",  Numerics.HeightEquation.NewtonMinIterations, 1);
		 if(CheckRange(Numerics.HeightEquation.NewtonMinIterations, 1, 100) == false)
			 Error("Invalid value. Numerics.HeightEquation.NewtonMinIterations must be between [1:100]. ");

					// Get Variable [JacobianUpdateInterval]
					// Interval at which the Jacobian is rebuilt. 1 (the default) means every iteration.
					sb->getin("JacobianUpdateInterval",  Numerics.HeightEquation.JacobianUpdateInterval, 1);
					 if(CheckRange(Numerics.HeightEquation.JacobianUpdateInterval, 1, 100) == false)
						 Error("Invalid value. Numerics.HeightEquation.JacobianUpdateInterval must be between [1:100]. ");

		// Get Variable [NewtonEpsilon]
		sb->getin("NewtonEpsilon",  Numerics.HeightEquation.NewtonEpsilon, 1e-010);

		//------------------------------------------------------------------
		//  The relative convergence test.
		//
		//  After falling quadratically, the Newton correction hits a floor set
		//  by the accuracy of the linear solver (around 1e-10 at this size) and
		//  stops there. Setting NewtonEpsilon below that floor (to 1e-12, say) means
		//  it can never be reached: the iteration spins all the way to
		//  NewtonMaxIterations, is then reported as "failed to converge", and the time
		//  step gets cut as well. Measured, that alone made the run 100 times
		//  slower.
		//
		//  NewtonEpsilonRel says "by how many orders of magnitude the correction
		//  must fall below the first correction of that step to count as
		//  converged". It depends neither on the scale of the solution nor on
		//  the time step. At 0 or below, only the absolute test is made.
		//------------------------------------------------------------------
		sb->getin("NewtonEpsilonRel", Numerics.HeightEquation.NewtonEpsilonRel, 1e-6);
		if(!CheckRange(Numerics.HeightEquation.NewtonEpsilonRel, 0.0, 1.0))
		{
			Error("Numerics.HeightEquation.NewtonEpsilonRel must be in [0:1].");
		}


		// ####### Get block [ManufacturedSolution] (optional, verification only) #####
		Numerics.ManufacturedSolution.Enable = false;
		Numerics.ManufacturedSolution.Mean = 0;
		Numerics.ManufacturedSolution.Amplitude = 0;
		Numerics.ManufacturedSolution.Omega = 0;
		Numerics.ManufacturedSolution.Wavelength.assign(2, 1.0);
		Numerics.ManufacturedSolution.ConcMean = 0.5;
		Numerics.ManufacturedSolution.ConcAmplitude = 0;
		Numerics.ManufacturedSolution.ConcOmega = 0;
		Numerics.ManufacturedSolution.ConcPhase = 0;
		Numerics.ManufacturedSolution.ConcWavelength.assign(2, 1.0);
		Numerics.ManufacturedSolution.ConcStratAmplitude = 0;
		Numerics.ManufacturedSolution.ConcStratOmega = 0;
		Numerics.ManufacturedSolution.ConcStratWavelength.assign(2, 1.0);
		Numerics.ManufacturedSolution.SubstrateMean = 0;
		Numerics.ManufacturedSolution.SubstrateAmplitude = 0;
		Numerics.ManufacturedSolution.SubstrateWavelength.assign(2, 1.0);
		{
			int nms = 0;
			try { nms = root_block->getBlock("Numerics")->getBlockSize("ManufacturedSolution"); }
			catch(const BlockNotFound&) { nms = 0; }
			if(nms > 1) Error("Block [ManufacturedSolution] may appear at most once (found " + std::to_string(nms) + ").");
			if(nms == 1)
			{
				Block *mb = root_block->getBlock("Numerics")->getBlock("ManufacturedSolution");
				mb->getin("Enable",    Numerics.ManufacturedSolution.Enable, false);
				mb->getin("Mean",      Numerics.ManufacturedSolution.Mean, 0.0);
				mb->getin("Amplitude", Numerics.ManufacturedSolution.Amplitude, 0.0);
				mb->getin("Omega",     Numerics.ManufacturedSolution.Omega, 0.0);
				mb->getin("Wavelength", Numerics.ManufacturedSolution.Wavelength, Numerics.ManufacturedSolution.Wavelength);
				if(CheckVectorSize(Numerics.ManufacturedSolution.Wavelength, 2) == false)
					Error("Numerics.ManufacturedSolution.Wavelength must be an array of exactly 2 numbers.");
				mb->getin("ConcMean",      Numerics.ManufacturedSolution.ConcMean, 0.5);
				mb->getin("ConcAmplitude", Numerics.ManufacturedSolution.ConcAmplitude, 0.0);
				mb->getin("ConcOmega",     Numerics.ManufacturedSolution.ConcOmega, 0.0);
				mb->getin("ConcPhase",     Numerics.ManufacturedSolution.ConcPhase, 0.0);
				mb->getin("ConcWavelength", Numerics.ManufacturedSolution.ConcWavelength, Numerics.ManufacturedSolution.ConcWavelength);
				if(CheckVectorSize(Numerics.ManufacturedSolution.ConcWavelength, 2) == false)
					Error("Numerics.ManufacturedSolution.ConcWavelength must be an array of exactly 2 numbers.");
				mb->getin("ConcStratAmplitude", Numerics.ManufacturedSolution.ConcStratAmplitude, 0.0);
				mb->getin("ConcStratOmega",     Numerics.ManufacturedSolution.ConcStratOmega, 0.0);
				mb->getin("ConcStratWavelength", Numerics.ManufacturedSolution.ConcStratWavelength, Numerics.ManufacturedSolution.ConcStratWavelength);
				if(CheckVectorSize(Numerics.ManufacturedSolution.ConcStratWavelength, 2) == false)
					Error("Numerics.ManufacturedSolution.ConcStratWavelength must be an array of exactly 2 numbers.");
				mb->getin("SubstrateMean",      Numerics.ManufacturedSolution.SubstrateMean, 0.0);
				mb->getin("SubstrateAmplitude", Numerics.ManufacturedSolution.SubstrateAmplitude, 0.0);
				mb->getin("SubstrateWavelength", Numerics.ManufacturedSolution.SubstrateWavelength, Numerics.ManufacturedSolution.SubstrateWavelength);
				if(CheckVectorSize(Numerics.ManufacturedSolution.SubstrateWavelength, 2) == false)
					Error("Numerics.ManufacturedSolution.SubstrateWavelength must be an array of exactly 2 numbers.");
				if(Numerics.ManufacturedSolution.Enable
				   && !(Numerics.ManufacturedSolution.ConcMean - std::fabs(Numerics.ManufacturedSolution.ConcAmplitude)
				          - std::fabs(Numerics.ManufacturedSolution.ConcStratAmplitude) > 0
				        && Numerics.ManufacturedSolution.ConcMean + std::fabs(Numerics.ManufacturedSolution.ConcAmplitude)
				          + std::fabs(Numerics.ManufacturedSolution.ConcStratAmplitude) < 1))
					Error("Numerics.ManufacturedSolution needs the concentration inside (0, 1).");
				if(Numerics.ManufacturedSolution.Enable
				   && !(Numerics.ManufacturedSolution.Mean > Numerics.ManufacturedSolution.Amplitude
				        && Numerics.ManufacturedSolution.Amplitude >= 0))
					Error("Numerics.ManufacturedSolution needs Mean > Amplitude >= 0.");
			}
		}

		// ####### Get block [Convection] (optional) #####
		//  The face interpolation of the convective terms. Absent, each model
		//  keeps the scheme it always used.
		Numerics.Convection.Concentration = "";
		Numerics.Convection.Height = "";
		{
			int ncv = 0;
			try { ncv = root_block->getBlock("Numerics")->getBlockSize("Convection"); }
			catch(const BlockNotFound&) { ncv = 0; }
			if(ncv > 1) Error("Block [Convection] may appear at most once (found " + std::to_string(ncv) + ").");
			if(ncv == 1)
			{
				Block *cb = root_block->getBlock("Numerics")->getBlock("Convection");
				cb->getin("Concentration", Numerics.Convection.Concentration, std::string(""));
				cb->getin("Height",        Numerics.Convection.Height,        std::string(""));
			}
			Numerics.Convection.ConcentrationGiven = !Numerics.Convection.Concentration.empty();
			Numerics.Convection.HeightGiven        = !Numerics.Convection.Height.empty();
			C_Interpolation_base4::t_Scheme tsc;
			{
				//  The scheme names, for a hint on a value given in the wrong case.
				vector<string> schemes;
				const string kl = C_Interpolation_base4::KeyList();
				for(size_t a = kl.find('`'); a != string::npos; a = kl.find('`', a + 1))
				{
					const size_t b = kl.find('\'', a);
					schemes.push_back(kl.substr(a + 1, b - a - 1));
					a = b;
				}
				const string &cv = Numerics.Convection.Concentration, &hv = Numerics.Convection.Height;
				if(Numerics.Convection.ConcentrationGiven && !C_Interpolation_base4::FromKey(cv, tsc))
					_caseHint = CaseHint(cv, schemes);
				else if(Numerics.Convection.HeightGiven && !C_Interpolation_base4::FromKey(hv, tsc))
					_caseHint = CaseHint(hv, schemes);
			}
			if(Numerics.Convection.ConcentrationGiven
			   && !C_Interpolation_base4::FromKey(Numerics.Convection.Concentration, tsc))
				Error("Invalid value. Numerics.Convection.Concentration = `" + Numerics.Convection.Concentration
				      + "' is not a scheme Yuzen has; it must be one of [" + C_Interpolation_base4::KeyList() + "]. ");
			if(Numerics.Convection.HeightGiven
			   && !C_Interpolation_base4::FromKey(Numerics.Convection.Height, tsc))
				Error("Invalid value. Numerics.Convection.Height = `" + Numerics.Convection.Height
				      + "' is not a scheme Yuzen has; it must be one of [" + C_Interpolation_base4::KeyList() + "]. ");
			//  The defaults, chosen by the comparison of the schemes (manual,
			//  sec:verification:schemes): UTOPIA for the concentration of both
			//  models, the Koren-limited scheme for the upwinded fluxes of the
			//  film height.
			if(!Numerics.Convection.ConcentrationGiven)
				Numerics.Convection.Concentration = "UTOPIA";
			if(!Numerics.Convection.HeightGiven)
				Numerics.Convection.Height = "Koren";
		}

			// ####### Get block [LinearSolver] #####
			// Check BlockSize.
			bnum = root_block->getBlock("Numerics")->getBlock("HeightEquation")->getBlockSize("LinearSolver");
			if(bnum != 1) Error("Block [LinearSolver] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
			sb = root_block->getBlock("Numerics")->getBlock("HeightEquation")->getBlock("LinearSolver");

			// Get Variable [MaxIterations]
			sb->getin("MaxIterations",  Numerics.HeightEquation.LinearSolver.MaxIterations, 5000);
			if(Numerics.HeightEquation.LinearSolver.MaxIterations < 1)
				Error("Invalid value. Numerics.HeightEquation.LinearSolver.MaxIterations must be 1 or larger. ");

			// Get Variable [Epsilon]
			sb->getin("Epsilon",  Numerics.HeightEquation.LinearSolver.Epsilon, 1e-013);

			// Get Variable [RowScaling]
			sb->getin("RowScaling",  Numerics.HeightEquation.LinearSolver.RowScaling, true);

			// Get Variable [PreconditionType]
			sb->getin("PreconditionType",  Numerics.HeightEquation.LinearSolver.PreconditionType, "ILU");
			//Set candidates.
			cand.clear();
			cand.push_back("None");
			cand.push_back("Jacobi");
			cand.push_back("ILU");
			 if(CheckCandidates(Numerics.HeightEquation.LinearSolver.PreconditionType,  cand) == false)
				 Error("Invalid value. Numerics.HeightEquation.LinearSolver.PreconditionType must be one of [`None',`Jacobi',`ILU']. ");

			// Get Variable [Backend]
			sb->getin("Backend",  Numerics.HeightEquation.LinearSolver.Backend, "auto");
			cand.clear();
			cand.push_back("auto");
			cand.push_back("eigen");
			cand.push_back("amgcl");
			//  Any other name is rejected here rather than silently treated as
			//  something else.
			 if(CheckCandidates(Numerics.HeightEquation.LinearSolver.Backend,  cand) == false)
				 Error("Invalid value. Numerics.HeightEquation.LinearSolver.Backend must be one of [`auto',`eigen',`amgcl']. ");

			// ####### END: block [LinearSolver] #####
		// ####### END: block [HeightEquation] #####

		// ####### Get block [ConcEquation] #####
		// Check BlockSize.
		bnum = root_block->getBlock("Numerics")->getBlockSize("ConcEquation");
		if(bnum != 1) Error("Block [ConcEquation] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
		sb = root_block->getBlock("Numerics")->getBlock("ConcEquation");


			// ####### Get block [LinearSolver] #####
			// Check BlockSize.
			bnum = root_block->getBlock("Numerics")->getBlock("ConcEquation")->getBlockSize("LinearSolver");
			if(bnum != 1) Error("Block [LinearSolver] is required and must appear exactly once (found " + std::to_string(bnum) + ").");
			sb = root_block->getBlock("Numerics")->getBlock("ConcEquation")->getBlock("LinearSolver");

			// Get Variable [MaxIterations]
			sb->getin("MaxIterations",  Numerics.ConcEquation.LinearSolver.MaxIterations, 5000);
			if(Numerics.ConcEquation.LinearSolver.MaxIterations < 1)
				Error("Invalid value. Numerics.ConcEquation.LinearSolver.MaxIterations must be 1 or larger. ");

			// Get Variable [Epsilon]
			sb->getin("Epsilon",  Numerics.ConcEquation.LinearSolver.Epsilon, 1e-013);

			// Get Variable [RowScaling]
			sb->getin("RowScaling",  Numerics.ConcEquation.LinearSolver.RowScaling, true);

			// Get Variable [PreconditionType]
			sb->getin("PreconditionType",  Numerics.ConcEquation.LinearSolver.PreconditionType, "ILU");
			//Set candidates.
			cand.clear();
			cand.push_back("None");
			cand.push_back("Jacobi");
			cand.push_back("ILU");
			 if(CheckCandidates(Numerics.ConcEquation.LinearSolver.PreconditionType,  cand) == false)
				 Error("Invalid value. Numerics.ConcEquation.LinearSolver.PreconditionType must be one of [`None',`Jacobi',`ILU']. ");

			// Get Variable [Backend]
			sb->getin("Backend",  Numerics.ConcEquation.LinearSolver.Backend, "auto");
			cand.clear();
			cand.push_back("auto");
			cand.push_back("eigen");
			cand.push_back("amgcl");
			//  Any other name is rejected here rather than silently treated as
			//  something else.
			 if(CheckCandidates(Numerics.ConcEquation.LinearSolver.Backend,  cand) == false)
				 Error("Invalid value. Numerics.ConcEquation.LinearSolver.Backend must be one of [`auto',`eigen',`amgcl']. ");

			// ####### END: block [LinearSolver] #####
		// ####### END: block [ConcEquation] #####
	// ####### END: block [Numerics] #####
}

//------------------------------------------------------------------
//  Block [Restart]
//------------------------------------------------------------------
void Input::ReadRestart()
{
	Block *sb = 0;
	Block absent;	// stands in for an optional block the file leaves out
	int bnum = 0;
	vector<string> cand; cand.clear();
	(void)sb; (void)bnum;

	// ####### Get block [Restart] #####
	// Check BlockSize.
	try
	{
		bnum = root_block->getBlockSize("Restart");
	}
	catch(const BlockNotFound&)
	{
		bnum = 0;
	}
	if(bnum > 1) Error("Block [Restart] may appear at most once (found " + std::to_string(bnum) + ").");
	//  Absent, the block reads as an empty one: every key takes its default.
	{
		sb = (bnum == 1) ? root_block->getBlock("Restart") : &absent;

		// Get Variable [FileIndex]
		sb->getin("FileIndex",  Restart.FileIndex, 0);

		// Get Variable [RestartTime]
		sb->getin("RestartTime",  Restart.RestartTime, 0);

		// Get Variable [Enable]
		sb->getin("Enable",  Restart.Enable, false);

		// Get Variable [ASCII]
		sb->getin("ASCII",  Restart.ASCII, false);

		// Get Variable [FileName]
		sb->getin("FileName",  Restart.FileName, "restartfile");

	}
	// ####### END: block [Restart] #####
}
