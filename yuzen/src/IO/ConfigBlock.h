/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ConfigBlock.h : the input-file reading layer
*
*  It offers only these six: newread / getin / getBlock / getBlockSize / write /
*  BlockNotFound.
*
*  Two formats are supported. The right one is picked from the extension.
*
*   .toml : read with toml++ (a public header-only library, MIT). Preferred for
*           anything new.
*   anything else : the .inp format, read by the self-contained parser in
*           ConfigBlock.cpp.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _CONFIG_BLOCK_INCLUDED_
#define _CONFIG_BLOCK_INCLUDED_

#include "YuzenUtil.h"
using namespace yuzen::util;
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace yuzen { namespace config {

//------------------------------------------------------------------
//  Exceptions
//------------------------------------------------------------------
/** @brief Thrown when a block of the requested name (or index) is not in the tree. */
//  There is no block of the given name. When a block is optional (it may or may
//  not be present), the caller catches this and uses its default instead.
class BlockNotFound : public std::runtime_error
{
public:
	explicit BlockNotFound(const std::string &what) : std::runtime_error(what) {}
};

/**
 * @brief Thrown when a variable of the requested name is not in the block, or
 *        when its text cannot be converted to the requested type.
 */
//  There is no variable of the given name. The getin overloads that take a
//  default do not throw; they leave the default in place.
class VariableNotFound : public std::runtime_error
{
public:
	explicit VariableNotFound(const std::string &what) : std::runtime_error(what) {}
};


//------------------------------------------------------------------
//  A single value. The conversion to a type happens when it is read out.
//------------------------------------------------------------------
/**
 * @brief One entry of a block, kept as text until getin() converts it.
 *
 * Holding the text rather than a parsed number is what lets the same entry be
 * read as an int, a double or a string, decided by the type getin() is handed.
 */
struct Value
{
	bool                     isVector;
	std::vector<std::string> items;		// strings have had their quotes stripped
	bool                     isString;	// whether it was a string literal (for writing back out)

	Value() : isVector(false), isString(false) {}
};


//------------------------------------------------------------------
//  The block tree
//------------------------------------------------------------------
//--------------------------------------------------------------------------
//  The block tree. One input file becomes one such tree.
//
//  Blocks can nest, and several blocks may share the same name (listing a
//  number of EllipseBump blocks to stack shapes on top of each other, for
//  instance). Values are held as text; the conversion to a type happens when
//  getin() reads them out.
//--------------------------------------------------------------------------
class Block
{
public:
	Block() {}
	/** @brief Deletes this block and, recursively, every block below it. */
	~Block();

	/**
	 * @brief Reads an input file and builds the block tree for it.
	 *
	 * The format is picked from the extension: ".toml" is parsed with toml++,
	 * anything else with the self-contained .inp parser.
	 *
	 * @param filename  path of the input file
	 * @return the root block. Ownership passes to the caller; deleting the root
	 *         deletes the whole tree.
	 */
	// Reads the file, telling .toml from .inp by the extension. Ownership of the
	// result passes to the caller.
	static Block* newread(const std::string &filename);

	//--- Blocks ---
	/**
	 * @brief How many blocks of this name this block directly contains.
	 * @param name  name of the child block
	 * @return the count, which is at least 1.
	 * @throws BlockNotFound if there is no block of that name at all.
	 */
	//  How many blocks of that name there are (none if 0)
	int    getBlockSize(const std::string &name) const;
	/**
	 * @brief Looks up one child block by name.
	 * @param name   name of the child block
	 * @param index  which one, when the same name appears more than once,
	 *               counted in input-file order from 0
	 * @return a pointer into the tree. Ownership stays with the tree, so the
	 *         caller must not delete it, and it dies with the root.
	 * @throws BlockNotFound if the name is absent or the index is out of range.
	 */
	//  Looks up the block at position index. Throws BlockNotFound if absent.
	//  What comes back is a pointer into the tree; ownership does not pass.
	Block* getBlock(const std::string &name, int index = 0) const;

	//--- Variables (no default: throws VariableNotFound if absent) ---
	/**
	 * @brief Reads a required variable of this block into v.
	 *
	 * Values are kept as text and converted here, so the target type is decided
	 * by the type of v rather than by how the value was written. A scalar read
	 * out of a list-valued entry takes the first item; a list read into a
	 * vector replaces whatever the vector held.
	 *
	 * @param name  variable name as written in the input file
	 * @param v     destination; left untouched when the call throws
	 * @throws VariableNotFound if the variable is not in this block, or if
	 *         its text is not a value of the type asked for: a number that
	 *         does not parse or is not finite, a non-integral value for an
	 *         int, a word that is not a boolean, an empty list, or a list of
	 *         several values where one is wanted. The message names the key.
	 */
	//  The type is decided by the argument. A failed text -> type conversion
	//  also throws, rather than give 0, or false, without a word.
	void getin(const std::string &name, int &v) const;
	void getin(const std::string &name, long long &v) const;
	void getin(const std::string &name, double &v) const;
	void getin(const std::string &name, bool &v) const;
	void getin(const std::string &name, std::string &v) const;
	void getin(const std::string &name, std::vector<double> &v) const;
	void getin(const std::string &name, std::vector<int> &v) const;
	void getin(const std::string &name, std::vector<std::string> &v) const;

	//--- Variables (with a default: v is left alone if absent) ---
	/**
	 * @brief Reads an optional variable, falling back to def when it is absent.
	 *
	 * Use this for anything the input file may leave out. Nothing is thrown
	 * for an absent variable; one that is present but does not convert throws
	 * VariableNotFound as above. The caller does not have to pre-load v,
	 * since def is written into it.
	 *
	 * @param name  name of the variable in this block
	 * @param v     receives the value read, or def when the variable is absent
	 * @param def   value stored in v when the variable is not in this block
	 */
	//  Use these for entries that may be left out. The caller does not have to
	//  put a default into v beforehand (def goes in as it is).
	void getin(const std::string &name, int &v, int def) const;
	void getin(const std::string &name, long long &v, long long def) const;
	void getin(const std::string &name, double &v, double def) const;
	void getin(const std::string &name, bool &v, bool def) const;
	void getin(const std::string &name, std::string &v, const std::string &def) const;
	void getin(const std::string &name, std::vector<double> &v, const std::vector<double> &def) const;
	void getin(const std::string &name, std::vector<int> &v, const std::vector<int> &def) const;
	void getin(const std::string &name, std::vector<std::string> &v, const std::vector<std::string> &def) const;

	/**
	 * @brief Whether this block holds a variable of this name (for a key
	 *        whose absence means something its default cannot say, such
	 *        as Output.FieldData.TimeInterval beside a list of times).
	 *        Asks nothing: the key must still be read with getin().
	 */
	bool has(const std::string &name) const { return _vars.count(name) != 0; }

	//--- Writing out (TOML format) ---
	/**
	 * @brief Writes this tree back out as TOML.
	 *
	 * Numbers are written with
	 * setprecision(17) so a round trip does not lose digits. A file that cannot
	 * be opened is silently skipped.
	 *
	 * @param filename  path of the file to write
	 */
	//  Writes the tree that was read back out as TOML, unchanged.
	void write(const std::string &filename) const;

	//--- Construction (used by the parsers) ---
	/** @brief Stores a variable under this name, replacing any previous value. */
	void        setValue(const std::string &name, const Value &v) { _vars[name] = v; }
	/**
	 * @brief Appends a child block of this name and returns it for filling in.
	 *
	 * Repeating a name is allowed; each call appends another instance rather
	 * than replacing the previous one.
	 *
	 * @return a pointer into the tree, owned by this block.
	 */
	Block*      addBlock(const std::string &name);
	// Name this block was declared under ("" for the root).
	const std::string& name() const { return _name; }
	/** @brief Sets the name this block was declared under; called by the parser. */
	void        setName(const std::string &n) { _name = n; }
	/** @brief Dotted path of this block from the root ("" for the root), for messages. */
	const std::string& path() const { return _path; }

	// Name and version of the reader layer, for the job information listing.
	/**
	 * @brief Lists every entry of the tree that the reader never asked for.
	 *
	 * Every getin() and every block lookup is recorded, found or not, so once
	 * the input has been read, an entry nobody asked for is one the program
	 * does not know: a misspelled key (\c alpha for \c Alpha), a key of the
	 * wrong block, or a block with a misspelled name. Unreported, such an entry
	 * would simply never be read, and the key it was meant to set would keep
	 * its default without a word.
	 *
	 * @param prefix  dotted path of this block, "" for the root
	 * @param out     one line per unknown entry is appended, giving its full
	 *                path and, where one is close, the name that was probably
	 *                meant
	 */
	void collectUnknown(const std::string &prefix, std::vector<std::string> &out) const;

private:
	//  Every key name the reader asked any block for, with the blocks it was
	//  asked of, so that a key put in the wrong block can be pointed to the
	//  right one.
	void collectAsked(const std::string &prefix,
	                  std::map<std::string, std::set<std::string> > &asked) const;
	void collectUnknown(const std::string &prefix, std::vector<std::string> &out,
	                    const std::map<std::string, std::set<std::string> > &asked) const;
public:

	/** @brief Name of the library that parsed the input, for the job information listing. */
	static std::string libraryName();
	/** @brief Version of the library that parsed the input, for the job information listing. */
	static std::string libraryVersion();

private:
	const Value* find(const std::string &name) const;
	void writeTo(std::ostream &os, const std::string &prefix) const;

	std::string _name;
	std::string _path;
	//  The value of a variable that is present, for conversion: one item,
	//  else VariableNotFound naming the key (an empty list, or a list where
	//  one value is wanted, does not silently give the default or the first
	//  item).
	const std::string& scalarItem(const std::string &name, const Value &v) const;
	std::string keyPath(const std::string &name) const
	{ return _path.empty() ? name : _path + "." + name; }
	std::map<std::string, Value>                _vars;
	std::map<std::string, std::vector<Block*> > _blocks;
	std::vector<std::string>                    _blockOrder;	// keeps the order for writing out

	//  What the reader has asked this block for, present or not. Mutable because
	//  asking is a const operation on the tree; see collectUnknown().
	mutable std::set<std::string>               _askedVars;
	mutable std::set<std::string>               _askedBlocks;

	Block(const Block&);
	Block& operator=(const Block&);
};

} }  // namespace yuzen::config

#endif //_CONFIG_BLOCK_INCLUDED_
