/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  Output.h : writing the results out (VTK / text / time series / ParaView PVD)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#ifndef _DEF_OUTPUT_INCLUDED_
#define _DEF_OUTPUT_INCLUDED_

#include <sstream>

#include "YuzenUtil.h"
using namespace yuzen::util;
#include "LubFields.h"

#include "TimeSteps.h"

#include "grid/YuzenGrid.h"


#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

using namespace yuzen::grid;

class NonDimensionalNumbers;


//-----------------------
//   C_Output_Interface
//-----------------------
//--------------------------------------------------------------------------
//  The part common to the output classes.
//
//  It holds the result directory, the file name head, the extension and the
//  sequence number. The sequence number can be set from outside with SetFNum so
//  that a restart can carry it over. The dimensional / non-dimensional switch
//  (Output.Dimensioned) lives here too, along with the reference to the
//  NonDimensionalNumbers used for the conversion.
//--------------------------------------------------------------------------
/**
 * @brief Stops the run after an output file could not be written.
 *
 * The writers report which file and why (no space left on the device, a file
 * size limit, no permission) through FatalError and return false; the model
 * calls this, which says what is left on disk and aborts with exit code 1. A
 * run that cannot write its results is not continued: it would end with exit
 * code 0 and empty or truncated files.
 */
void StopOnOutputError();

class C_Output_Interface
{
public:
	/**
	 * @brief Path of the file the next Write would produce.
	 *
	 * For the sequence writers this is the result directory, the file name head,
	 * the current sequence number and the extension; it changes as the sequence
	 * advances. For the single-file writers it is fixed.
	 */
	//  Path of the file that is written next
	virtual std::string FileName() const = 0;

protected:
	C_Output_Interface() : _timeStep(0), _NonDim(0), _dimensioned(false), _constructed(false), _dir(""), _head(""), _ext(""), _fnum(0)  {}
	virtual ~C_Output_Interface(){;}

	void Build(C_TimeSteps &tt, NonDimensionalNumbers &nd, bool dimensioned, const FilePath &resultDir, const std::string &fhead, const std::string &ext, int fst = 0)
	{
		_timeStep = &tt;
		_NonDim = &nd;
		_dimensioned = dimensioned;
		_dir = resultDir;
		_head = fhead;
		_ext = ext;
		_fnum = fst;
	}

	C_TimeSteps *_timeStep;
	NonDimensionalNumbers *_NonDim;
	bool _dimensioned;
	bool _constructed;
	FilePath _dir;
	std::string _head;
	std::string _ext;
	int _fnum;

private:
	C_Output_Interface(const C_Output_Interface&);
	void operator=(const C_Output_Interface&);

};

//-----------------------
//   C_PVD_File
//-----------------------
//--------------------------------------------------------------------------
//  The ParaView time-series file (.pvd).
//  It lists the individual .vts files with their time stamps; opening it gives
//  an animation.
//--------------------------------------------------------------------------
class C_PVD_File : public C_Output_Interface
{
public:
	C_PVD_File() :  _header("field"), _delta_t(1), _time(0), _start_time(0), _fnum_st(0) {}

	virtual ~C_PVD_File() {}
	/**
	 * @brief Creates the .pvd file and writes its opening tags.
	 *
	 * @param tt        time-step object the output times are taken from
	 * @param nd        non-dimensional numbers, used to convert back to dimensional values
	 * @param dimensioned  write dimensional values instead of non-dimensional ones, converting through nd
	 * @param resultDir directory the file is written into
	 * @param fname     name of the .pvd file itself (extension added here)
	 * @param head      file name head of the .vts files it will point at
	 * @param st_time   time stamp given to the first entry
	 * @param dt        time added per entry; the .pvd assumes the field files
	 *                  are written at a fixed interval
	 * @param fst       sequence number the first entry refers to
	 */
	void Build(C_TimeSteps &tt, NonDimensionalNumbers &nd, bool dimensioned, const FilePath &resultDir, const std::string &fname, const std::string &head, double st_time, double dt, int fst);
	/**
	 * @brief Appends one entry, then advances the sequence number by one and
	 *        the time stamp by dt. Call it once for every .vts written.
	 */
	//  Call this once a .vts has been written. It appends one line carrying
	//  the time stamp and the sequence number.
	bool operator++();
	bool operator++(int){ return operator++();}
	/**
	 * @brief Appends one entry stamped with the time actually reached.
	 *
	 * \c operator++ stamps each entry with the nominal schedule,
	 * \c start+n*dt, which is right while the schedule is being kept and
	 * wrong for an entry written off it --- the final output of a run, whose
	 * state sits at whatever time the last step landed on. Prefer this.
	 *
	 * @param timeNow  the time of the state just written, on the same scale
	 *                 as \c st_time was given to Build()
	 */
	bool Add(double timeNow);
	/**
	 * @brief Continue the collection of the run a restart continues.
	 *
	 * Build() starts the file afresh. When a run is continued exactly (from
	 * a field*.state) into the result directory it was written in, the
	 * entries that file already held up to the time restarted from are put
	 * back, so that the .pvd goes on listing the whole run. Entries after
	 * it belong to the part of the earlier run that is being done again,
	 * and are dropped. Which is which is told by the field file an entry
	 * points at, kept when numbered below the one set by SetFNum() (call
	 * that first); the six-digit time stamp decides only for an entry
	 * pointing at some other file.
	 * @param tmax  the latest time stamp kept, on the scale of the file
	 */
	bool Resume(double tmax);
	/**
	 * @brief Sets the sequence number the next entry refers to.
	 *
	 * A restarted run numbers its field files on from Restart.FileIndex, and
	 * the .pvd has to point at those files, not at field00000.vts onwards.
	 */
	void SetFNum(int n){_fnum = n;}
	/**
	 * @brief Writes the file a last time. The file is valid XML, its
	 *        closing tags included, after every entry: it is rewritten whole
	 *        (under a temporary name, then renamed) each time an entry is
	 *        added, so ParaView opens the .pvd of a run still going or killed.
	 * @return false if the file could not be written (reported through
	 *         FatalError), as for Add() and Resume().
	 */
	bool Close();

	virtual std::string FileName() const {return _fpath;}

private:
	std::string _header;
	double _delta_t;
	double _time;
	double _start_time;
	std::vector<std::string> _previous;	// the entries the file held before Build()
	std::vector<std::string> _entries;	// the entries the file holds now
	int _fnum_st;
	FilePath _fpath;

	bool Rewrite();

	C_PVD_File(const C_PVD_File&);
	void operator=(const C_PVD_File&);

};

//-----------------------
//   OutputVTK
//-----------------------
//--------------------------------------------------------------------------
//  Writes a VTK structured grid (.vts). The values are embedded as raw binary.
//  ParaView can open it, and the binary restart reads this format too.
//--------------------------------------------------------------------------
class OutputVTK  : public C_Output_Interface
{
public:
	OutputVTK(){}
	virtual ~OutputVTK(){}

	/**
	 * @brief Sets up the writer. No file is created until the first Write.
	 *
	 * @param tt         time-step object the output times are taken from
	 * @param nd         non-dimensional numbers, used to convert back to dimensional values
	 * @param resultDir  directory the file is written into
	 * @param dimensioned  write dimensional values instead of non-dimensional
	 *                     ones, converting through nd
	 * @param fhead        file name head; the sequence number and ".vts" are
	 *                     appended to it
	 * @param fst          sequence number the first file gets
	 */
	void Build(C_TimeSteps &tt, NonDimensionalNumbers &nd, bool dimensioned, const FilePath &resultDir, const std::string &fhead, int fst);
	/**
	 * @brief Writes one .vts file for this field and advances the sequence
	 *        number. There is one overload per kind of field.
	 * The file is written under its name with ".tmp" added and renamed when
	 * it is whole, so a run stopped while writing never leaves a truncated
	 * file under the real name.
	 *
	 * @return false if the writer was never built or the file could not be
	 *         written (the reason is reported through FatalError; the caller
	 *         stops the run with StopOnOutputError()); in that case the
	 *         sequence number does not move.
	 */
	//  Writes one file and advances the sequence number by one. One overload
	//  per kind of field.
	bool Write(Field_WM1D_Out &ff);
	bool Write(Field_WM2D_Out &ff);
	bool Write(Field_CZ &ff);
	/**
	 * @brief Sets the sequence number of the next file, so a restarted run
	 *        carries on numbering from where the previous one stopped.
	 */
	//  Sets the sequence number from outside (so a restart numbers on from
	//  where it left off)
	void SetFNum(int n){_fnum = n;}

	virtual std::string FileName() const;

private:
	OutputVTK(const OutputVTK&);
	void operator=(const OutputVTK&);
};


//-----------------------
//   OutputNeutral
//-----------------------
//--------------------------------------------------------------------------
//  Writes the text format (.dat). Tab separated, one point per line.
//
//  ** This format is also the input of the ASCII restart. **
//  Changing the column order or the number of digits therefore puts writing
//  and reading out of step. The digits are written in full, with
//  setprecision(17). At the default 6 digits the pressure is destroyed once
//  the film thickness is put through the second difference
//  (docs/PHASE7.md section 7.7).
//--------------------------------------------------------------------------
class OutputNeutral : public C_Output_Interface
{
public:
	OutputNeutral(){}
	virtual ~OutputNeutral(){;}

	/**
	 * @brief Sets up the writer. No file is created until the first Write.
	 *
	 * @param tt         time-step object the output times are taken from
	 * @param nd         non-dimensional numbers, used to convert back to dimensional values
	 * @param resultDir  directory the file is written into
	 * @param dimensioned  write dimensional values instead of non-dimensional
	 *                     ones, converting through nd
	 * @param fhead        file name head; the sequence number and ".dat" are
	 *                     appended to it
	 * @param fst          sequence number the first file gets
	 */
	void Build(C_TimeSteps &tt, NonDimensionalNumbers &nd, bool dimensioned, const FilePath &resultDir, const std::string &fhead, int fst);
	/**
	 * @brief Writes one .dat file for this field and advances the sequence
	 *        number. There is one overload per kind of field.
	 *
	 * The file produced here is also what C_Restart reads back as an ASCII
	 * restart, so the column order and the number of digits are part of the
	 * format rather than a presentation choice.
	 *
	 * The file is written under its name with ".tmp" added and renamed when
	 * it is whole, so a run stopped while writing never leaves a truncated
	 * file under the real name.
	 *
	 * @return false if the writer was never built or the file could not be
	 *         written (the reason is reported through FatalError; the caller
	 *         stops the run with StopOnOutputError()); in that case the
	 *         sequence number does not move.
	 */
	//  Writes one file and advances the sequence number by one. One overload
	//  per kind of field.
	bool Write(Field_WM1D_Out &ff);
	bool Write(Field_WM2D_Out &ff);
	bool Write(Field_CZ &ff);
	/**
	 * @brief Sets the sequence number of the next file, so a restarted run
	 *        carries on numbering from where the previous one stopped.
	 */
	//  Sets the sequence number from outside (so a restart numbers on from
	//  where it left off)
	void SetFNum(int n){_fnum = n;}

	virtual std::string FileName() const;

private:
	OutputNeutral(const OutputNeutral&);
	void operator=(const OutputNeutral&);
};

//-----------------------
//   OutputTimeline
//-----------------------
//--------------------------------------------------------------------------
//  Writes a time series (History.dat / CalcInfo.dat / the probes).
//  It writes the header line first and then appends one line per step.
//--------------------------------------------------------------------------
class OutputTimeline : public C_Output_Interface
{
public:
	OutputTimeline() : _num_vals(0), _labels(0), _fpath("") {}
	virtual ~OutputTimeline() {}

	/**
	 * @brief Creates the file and writes the column header line.
	 *
	 * An existing file of the same name is truncated, so Build is called once
	 * per run and not per output.
	 *
	 * @param tt       time-step object the output times are taken from
	 * @param nd       non-dimensional numbers, used to convert back to dimensional values
	 * @param dimensioned  write dimensional values instead of non-dimensional ones, converting through nd
	 * @param resultDir  directory the file is written into
	 * @param fhead    file name; ".dat" is appended, with no sequence number
	 * @param num      number of value columns after the time column
	 * @param labels   column labels, num of them
	 * @param pre_str  optional line written above the header, used for comments
	 */
	void Build(C_TimeSteps &tt, NonDimensionalNumbers &nd, bool dimensioned, const FilePath &resultDir, const std::string &fhead,
			 int num, std::vector<std::string> &labels, std::string pre_str = "");
	/**
	 * @brief Appends one line to the file.
	 * @param oss  the whole line, time first and then the values, already
	 *             assembled tab-separated by the caller. The newline is added
	 *             here.
	 * @return false if the line could not be written (reported through
	 *         FatalError; the caller stops the run with StopOnOutputError()).
	 */
	//  Appends one line. The caller assembles the time and the values tab
	//  separated and passes them in.
	bool Write(const ostringstream& oss);
	virtual std::string FileName() const {return _fpath;}
	/**
	 * @brief Continue the file of the run a restart continues.
	 *
	 * Build() truncates the file. When a run is continued exactly (from a
	 * field*.state) into the result directory it was written in, the rows
	 * the file held up to the time restarted from are put back under the
	 * new header; rows after it belong to the part of the earlier run that
	 * is being done again, and are dropped.
	 * A last row without its newline -- the one being written when the
	 * earlier run was killed or ran out of space -- is dropped.
	 * @param tmax  the latest time kept, on the scale of the time column
	 */
	bool Resume(double tmax);

private:
	int _num_vals;
	std::vector<std::string> _labels;
	FilePath _fpath;
	std::vector<std::string> _previous;	// the lines the file held before Build()

	OutputTimeline(const OutputTimeline&);
	void operator=(const OutputTimeline&);
};



#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_DEF_OUTPUT_INCLUDED_
