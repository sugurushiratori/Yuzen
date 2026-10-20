/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Input.h : definition of every parameter of the input file
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#ifndef _DEF_Input_INCLUDED_
#define _DEF_Input_INCLUDED_

#include <memory>
#include "ConfigBlock.h"
using namespace yuzen::config;
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <stdexcept>

using std::string;
using std::list;
using std::vector;
using std::cout;
using std::endl;

namespace yuzen{
class Input
{
private:
	// The tree structure of the input file. Input is copied by value
	// (NonDimensionalNumbers.cpp / Evaporation.cpp), so ownership is shared.
	std::shared_ptr<Block> root_block;
	int CheckBlockSize(Block *b, string s);
	template<typename Type> bool CheckVectorSize(vector<Type> &vec, const size_t size_must_be);
	bool CheckRange(const int    val, const int    min, const    int max);
	bool CheckRange(const double val, const double min, const double max);
	bool CheckRange(vector<int>    &vec, const int    min, const int    max);
	bool CheckRange(vector<double> &vec, const double min, const double max);
	bool CheckCandidates(const string& str, vector<string> &candidates);
	bool CheckCandidates(const vector<string>& str, vector<string> &candidates);
	void Error(const string &mes);
	//  Set by a failed CheckCandidates() when the value differs from one of
	//  the candidates only in case; Error() appends it.
	string _caseHint;
public:
	Input();
	virtual ~Input(){}

	/**
	 * @brief Reads an input file and fills in every parameter of this object.
	 *
	 * The format (.inp or .toml) is picked from the extension by the reader
	 * layer. Ranges and allowed values are checked while reading; anything out
	 * of range is reported and ends the run rather than being returned as an
	 * error. Items left out of the file get their built-in defaults. Checks
	 * that span several parameters are not done here - see CheckParameters.
	 *
	 * @param filename  path of the input file
	 * @return 0. A failed read does not return; it terminates the process.
	 */
	//  Reads the input file (.inp / .toml) and fills in every parameter.
	//  Ranges are checked while reading; anything out of range prints a warning
	//  and aborts the run. Items left out get their defaults.
	int Read(const std::string &filename);

private:
	// Read(), at 1,918 lines, split up one top-level block at a time
	void ReadRoot();
	void ReadGrid();
	void ReadBoundary();
	void ReadPhysics();
	void ReadSubstrate();
	void ReadInitial();
	void ReadFluid();
	void ReadOutput();
	void ReadNumerics();
	void ReadRestart();
public:

	/**
	 * @brief Writes the parsed input file back out in TOML form.
	 *
	 * What is written is the block tree as it was read, not the members of this
	 * object, so defaults that were filled in for missing items do not appear.
	 * Read has to have been called first.
	 *
	 * @param filename  path of the file to write
	 * @return 0.
	 */
	int Write(const std::string &filename);

// --------------------------------
// User-defined block Structures
// --------------------------------
	string SolverName;
	string JobName;
	string LengthScaleType;
	double LengthScale;

	/** @brief Mirrors the [Grid] block: domain extent and number of cells. */
	class I_Grid
	{
	public:
		string Dimension;
		double Xmin;
		double Xmax;
		double Ymin;
		double Ymax;
		int Nx;
		int Ny;
		int Nz;
	};
	I_Grid Grid;
	/**
	 * @brief Mirrors the [Boundary] block: one sub-block per domain edge
	 *        (XmaxBC, XminBC, YmaxBC, YminBC).
	 */
	class I_Boundary
	{
	public:
		class I_XmaxBC
		{
		public:
			string Type;
		};
		I_XmaxBC XmaxBC;
		class I_XminBC
		{
		public:
			string Type;
		};
		I_XminBC XminBC;
		class I_YmaxBC
		{
		public:
			string Type;
		};
		I_YmaxBC YmaxBC;
		class I_YminBC
		{
		public:
			string Type;
		};
		I_YminBC YminBC;
	};
	I_Boundary Boundary;
	/**
	 * @brief Mirrors the [Physics] block: which effects are switched on and the
	 *        coefficients they need.
	 *
	 * Each effect is its own sub-block with an Enable flag, so terms can be
	 * turned on and off one at a time.
	 */
	class I_Physics
	{
	public:
		class I_ConcField
		{
		public:
			string Model;
			bool EnableConvection;
			bool EnableDiffusion;
		};
		I_ConcField ConcField;
		class I_LaplacePressure
		{
		public:
			bool Enable;
		};
		I_LaplacePressure LaplacePressure;
		class I_SolutoCapillarity
		{
		public:
			bool Enable;
		};
		I_SolutoCapillarity SolutoCapillarity;
		//  Disjoining pressure (wettability / contact line)
		class I_DisjoiningPressure
		{
		public:
			bool Enable;
			double ContactAngle;	// equilibrium contact angle [deg]
			double ExponentN;		// exponent on the repulsive side
			double ExponentM;		// exponent on the attractive side
		};
		I_DisjoiningPressure DisjoiningPressure;
		class I_VerticalGravity
		{
		public:
			bool Enable;
			double Acceleration;
		};
		I_VerticalGravity VerticalGravity;
		class I_HorizontalGravity
		{
		public:
			bool Enable;
			vector<double> Acceleration;
		};
		I_HorizontalGravity HorizontalGravity;
		class I_Centrifugal
		{
		public:
			bool Enable;
			double RotationRate;
			vector<double> RotationAxis;
		};
		I_Centrifugal Centrifugal;
		//  Navier slip at the substrate (Mobility.h)
		class I_Slip
		{
		public:
			//  Initialized here: a test that fills an Input by hand and does
			//  not know the block gets no slip, not an indeterminate flag
			bool Enable = false;
			double Length = 0.0;	// slip length lambda [m], >= 0
		};
		I_Slip Slip;
		class I_Evaporation
		{
		public:
			bool Enable;
			double Coef;
			string SpatialFunc;
			string ConcFunc;
			double NegPowExponent;
			double NegPowRmax;
			double PowerN;
			double ConcR;
			string EndConcType;
			double EndConcValue;
		};
		I_Evaporation Evaporation;
		class I_Source
		{
		public:
			bool Enable;
			double Radius;
			double InletRate;
			double Conc;
			string ConcType;
			double TimeConstant;
			double ArmLength;
			vector<double> ArmCenter;
			double RotationRate;
			double StartAngle;
			double RampTime;
			double PlateauTime;
			int SubRegions;
		};
		I_Source Source;
	};
	I_Physics Physics;
	/**
	 * @brief Mirrors the [Substrate] block: the topography the film sits on.
	 *
	 * The shape sub-blocks (RectBump, EllipseBump, ...) may each appear any
	 * number of times, which is why they are held as vectors; they are combined
	 * in the order given by their Order field. See C_GenerateShape.
	 */
	class I_Substrate
	{
	public:
		double TotalHeight;
		class I_RectBump
		{
		public:
			int Order;
			string CreateType;
			double Base;
			double Height;
			vector<double> Center;
			vector<double> Width;
			vector<double> SlopeWidth;
		};
		std::vector<I_RectBump> RectBump;
		class I_EllipseBump
		{
		public:
			int Order;
			string CreateType;
			double Base;
			double Height;
			vector<double> Center;
			vector<double> Width;
			double SlopeWidth;
		};
		std::vector<I_EllipseBump> EllipseBump;
		class I_PeriodicWave
		{
		public:
			int Order;
			string CreateType;
			double Base;
			double Amplitude;
			vector<double> Wavelength;
			vector<double> Phase;
		};
		std::vector<I_PeriodicWave> PeriodicWave;
		class I_FlatPlane
		{
		public:
			int Order;
			string CreateType;
			double Height;
			vector<double> Center;
			vector<double> Angle;
		};
		std::vector<I_FlatPlane> FlatPlane;
		class I_MexicanHat
		{
		public:
			int Order;
			string CreateType;
			double Base;
			double Height;
			vector<double> Center;
			vector<double> Width;
			double Beta;
		};
		std::vector<I_MexicanHat> MexicanHat;
		class I_SmoothConvex
		{
		public:
			int Order;
			string CreateType;
			double Base;
			double Height;
			vector<double> Center;
			vector<double> Width;
		};
		std::vector<I_SmoothConvex> SmoothConvex;
	};
	I_Substrate Substrate;
	/**
	 * @brief Mirrors the [Initial] block: the state the run starts from.
	 *
	 * Concentration sets the initial solute field; SurfaceShape builds the
	 * initial film thickness out of repeatable shape sub-blocks, the same way
	 * Substrate builds the topography.
	 */
	class I_Initial
	{
	public:
		class I_Concentration
		{
		public:
			string Type;
			double Value;
		};
		I_Concentration Concentration;
		class I_SurfaceShape
		{
		public:
			double TotalHeight;
			double PrecursorRatio;
			class I_RectBump
			{
			public:
				int Order;
				string CreateType;
				double Base;
				double Height;
				vector<double> Center;
				vector<double> Width;
				vector<double> SlopeWidth;
			};
			std::vector<I_RectBump> RectBump;
			class I_EllipseBump
			{
			public:
				int Order;
				string CreateType;
				double Base;
				double Height;
				vector<double> Center;
				vector<double> Width;
				double SlopeWidth;
			};
			std::vector<I_EllipseBump> EllipseBump;
			class I_PeriodicWave
			{
			public:
				int Order;
				string CreateType;
				double Base;
				double Amplitude;
				vector<double> Wavelength;
				vector<double> Phase;
			};
			std::vector<I_PeriodicWave> PeriodicWave;
			class I_FlatPlane
			{
			public:
				int Order;
				string CreateType;
				double Height;
				vector<double> Center;
				vector<double> Angle;
			};
			std::vector<I_FlatPlane> FlatPlane;
			class I_MexicanHat
			{
			public:
				int Order;
				string CreateType;
				double Base;
				double Height;
				vector<double> Center;
				vector<double> Width;
				double Beta;
			};
			std::vector<I_MexicanHat> MexicanHat;
			class I_SmoothConvex
			{
			public:
				int Order;
				string CreateType;
				double Base;
				double Height;
				vector<double> Center;
				vector<double> Width;
			};
			std::vector<I_SmoothConvex> SmoothConvex;
			class I_EdgeBump
			{
			public:
				int Order;
				string CreateType;
				double Base;
				double Height;
				string Location;
				double alpha;
				double Lin;
				double Lout;
			};
			std::vector<I_EdgeBump> EdgeBump;
		};
		I_SurfaceShape SurfaceShape;
	};
	I_Initial Initial;
	/**
	 * @brief Mirrors the [Fluid] block: material properties of the liquid.
	 *
	 * Surface tension, viscosity and diffusivity are given as a value at a
	 * reference concentration plus a coefficient or exponent describing how
	 * they change with concentration.
	 */
	class I_Fluid
	{
	public:
		string Name;
		class I_Density
		{
		public:
			double Value;
		};
		I_Density Density;
		class I_SurfaceTension
		{
		public:
			double BaseValue;
			double BaseConc;
			double Coef;
		};
		I_SurfaceTension SurfaceTension;
		class I_Viscosity
		{
		public:
			double BaseValue;
			double BaseConc;
			double Exponent;
		};
		I_Viscosity Viscosity;
		class I_Diffusivity
		{
		public:
			double BaseValue;
			double BaseConc;
			double Exponent;
		};
		I_Diffusivity Diffusivity;
		class I_MolarMass
		{
		public:
			double Solvent;
			double Resin;
		};
		I_MolarMass MolarMass;
	};
	I_Fluid Fluid;
	/**
	 * @brief Mirrors the [Output] block: what is written, how often, and in
	 *        which units.
	 *
	 * Dimensioned selects dimensional or non-dimensional values in the written
	 * files. The sub-blocks set the output interval of the field files, the
	 * history file, the calculation-information file and the probes.
	 */
	class I_Output
	{
	public:
		bool Overwrite;
		bool Dimensioned;
		class I_FieldData
		{
		public:
			double TimeInterval;
			bool NeutralFile;
			bool VTKFile;
			//  Output at listed times instead of every TimeInterval:
			//  Times, or a LogSchedule (resolved into Times by the reader)
			std::vector<double> Times;		// [s], ascending, > 0
			bool ExactTimes = false;		// land a step on each
			bool HasTimeInterval = false;		// TimeInterval was given
			bool WritePotential = false;		// the driving potential as PHI
			class I_LogSchedule
			{
			public:
				bool Present = false;
				double Start = 0;		// [s]
				double End = 0;			// [s]
				int PerDecade = 0;
				bool ExactTimes = false;
			};
			I_LogSchedule LogSchedule;
			bool Listed() const {return !Times.empty();}
		};
		I_FieldData FieldData;
		class I_GeneralHistory
		{
		public:
			double TimeInterval;
		};
		I_GeneralHistory GeneralHistory;
		class I_CalcInfo
		{
		public:
			double TimeInterval;
		};
		I_CalcInfo CalcInfo;
		class I_ProbeList
		{
		public:
			double TimeInterval;
			class I_Probe
			{
			public:
				double px;
				double py;
				double pz;
			};
			std::vector<I_Probe> Probe;
		};
		I_ProbeList ProbeList;
		//  Energy.dat: the free energy and its dissipation (well-mixed),
		//  written with History.dat
		class I_Energy
		{
		public:
			bool Enable = false;
		};
		I_Energy Energy;
	};
	I_Output Output;
	/**
	 * @brief Mirrors the [Numerics] block: time stepping, stopping conditions
	 *        and solver settings.
	 *
	 * TimeStepControl holds the adaptive time step, HeightEquation the Newton
	 * iteration for the film-thickness equation, and both HeightEquation and
	 * ConcEquation carry their own LinearSolver sub-block.
	 */
	class I_Numerics
	{
	public:
		double InitialDt;
		long long MaxSteps;		// 0: no limit
		double EndTime;
		//  Threshold of the steady-state test. Once the relative rate of change
		//  of the film thickness per unit (non-dimensional) time has stayed
		//  below this value, the run is taken to be steady and stops.
		//  At 0 or below, the steady-state test is not made at all.
		double SteadyStateCriterion;
		int SteadyStateSteps;
		class I_TimeStepControl
		{
		public:
			bool Enable;
			double MaxDtFactor;
			double MinDt;
			double Criterion;
			//  Lower bound of the weight used when measuring the error, as a
			//  fraction of the global maximum. 0 = divide by the local film
			//  thickness; 1 = divide by the global maximum alone.
			double ScaleFloor;
			int KeepSteps;
			double ChangeFactor;
		};
		I_TimeStepControl TimeStepControl;
		class I_HeightEquation
		{
		public:
			double NewtonEpsilon;
			//  Relative convergence test, as a ratio against the first
			//  correction of that step. At 0 or below it is switched off and
			//  only the absolute test is made.
			double NewtonEpsilonRel;
			int NewtonMaxIterations;
			int NewtonMinIterations;				// minimum number of iterations; 1 means the convergence test alone
			int JacobianUpdateInterval;	// modified Newton method; 1 means every iteration
			class I_LinearSolver
			{
			public:
				int MaxIterations;
				double Epsilon;
				bool RowScaling;
				string PreconditionType;
				// Picks the implementation of the linear solver ("auto" | "eigen" | "amgcl")
				string Backend;
			};
			I_LinearSolver LinearSolver;
		};
		I_HeightEquation HeightEquation;
		class I_ConcEquation
		{
		public:
			class I_LinearSolver
			{
			public:
				int MaxIterations;
				double Epsilon;
				bool RowScaling;
				string PreconditionType;
				// Picks the implementation of the linear solver ("auto" | "eigen" | "amgcl")
				string Backend;
			};
			I_LinearSolver LinearSolver;
		};
		I_ConcEquation ConcEquation;
		//  Verification only: add to the film-height equation the source that
		//  makes h = Mean + Amplitude cos(kx x - Omega t) cos(ky y) an exact
		//  solution, with k = 2 pi / Wavelength (ky = 0 in 1-D). The run is then
		//  compared with that h by tools/mms.py. Lengths in meters, Omega in
		//  rad/s. See ManufacturedSource.h for which terms it covers.
		class I_ManufacturedSolution
		{
		public:
			bool Enable;
			double Mean;
			double Amplitude;
			double Omega;
			vector<double> Wavelength;
			//  The manufactured concentration and substrate
			double ConcMean;
			double ConcAmplitude;
			double ConcOmega;
			double ConcPhase;
			vector<double> ConcWavelength;
			//  ConcZ only: a part that varies through the film,
			//  ConcStratAmplitude cos(2 pi x/Lx - ConcStratOmega t)
			//  cos(2 pi y/Ly) cos(pi zeta), zeta = (z - s)/h
			double ConcStratAmplitude;
			double ConcStratOmega;
			vector<double> ConcStratWavelength;
			double SubstrateMean;
			double SubstrateAmplitude;
			vector<double> SubstrateWavelength;
		};
		I_ManufacturedSolution ManufacturedSolution;
		//  [Numerics.Convection] (optional): the face interpolation of the
		//  convective terms. Concentration is that of the advection of c,
		//  Height that of the upwinded fluxes of the film-height equation
		//  (in-plane gravity, centrifugal force, Marangoni) of the
		//  well-mixed model. After reading they hold the scheme in use, the
		//  model's default where the key is absent (Concentration: "UTOPIA";
		//  Height: "Koren"). HeightGiven records whether
		//  Height was written (CZ has no upwinded term in its height
		//  equation and refuses it).
		class I_Convection
		{
		public:
			//  Empty (as in an Input built in code rather than read) means
			//  the model's default
			I_Convection() : ConcentrationGiven(false), HeightGiven(false) {}
			string Concentration;
			string Height;
			bool ConcentrationGiven;
			bool HeightGiven;
		};
		I_Convection Convection;
	};
	I_Numerics Numerics;
	/**
	 * @brief Mirrors the [Restart] block: which previous result to continue from.
	 *
	 * ASCII selects the .dat form over the binary .vts one, RestartTime is the
	 * time the continued run starts at, and FileIndex is the number the output
	 * file numbering picks up from. See C_Restart.
	 */
	class I_Restart
	{
	public:
		bool Enable;
		bool ASCII;
		string FileName;
		double RestartTime;
		int FileIndex;
	};
	I_Restart Restart;
};


}	// end namespace yuzen

#endif	//_DEF_Input_INCLUDED_
