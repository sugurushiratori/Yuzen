/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ConfigBlock.cpp : lexer and parser for the block-structured input (.inp), and
*                    the implementation of the .toml loader
*
*  ConfigBlock.cpp
*
*  The .inp parser and the .toml loader.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>

#include "ConfigBlock.h"

#include <toml++/toml.hpp>

namespace yuzen { namespace config {

//==================================================================
//  Block
//==================================================================
Block::~Block()
{
	for(std::map<std::string, std::vector<Block*> >::iterator it = _blocks.begin();
	    it != _blocks.end(); ++it)
	{
		for(size_t i = 0; i < it->second.size(); ++i) delete it->second[i];
	}
}

Block* Block::addBlock(const std::string &name)
{
	if(_blocks.find(name) == _blocks.end()) _blockOrder.push_back(name);
	Block *b = new Block();
	b->setName(name);
	b->_path = keyPath(name);
	_blocks[name].push_back(b);
	return b;
}

int Block::getBlockSize(const std::string &name) const
{
	_askedBlocks.insert(name);
	std::map<std::string, std::vector<Block*> >::const_iterator it = _blocks.find(name);
	if(it == _blocks.end()) throw BlockNotFound("Block not found: [" + name + "]");
	return (int)it->second.size();
}

Block* Block::getBlock(const std::string &name, int index) const
{
	_askedBlocks.insert(name);
	std::map<std::string, std::vector<Block*> >::const_iterator it = _blocks.find(name);
	if(it == _blocks.end()) throw BlockNotFound("Block not found: [" + name + "]");
	if(index < 0 || index >= (int)it->second.size())
		throw BlockNotFound("Block index out of range: [" + name + "]");
	return it->second[(size_t)index];
}

const Value* Block::find(const std::string &name) const
{
	_askedVars.insert(name);
	std::map<std::string, Value>::const_iterator it = _vars.find(name);
	return (it == _vars.end()) ? 0 : &it->second;
}

namespace {

//  Edit distance, case-insensitive, for suggesting what an unknown name was
//  probably meant to be.
size_t Distance(const std::string &a, const std::string &b)
{
	std::vector<size_t> prev(b.size() + 1), cur(b.size() + 1);
	for(size_t j = 0; j <= b.size(); ++j) prev[j] = j;
	for(size_t i = 1; i <= a.size(); ++i)
	{
		cur[0] = i;
		for(size_t j = 1; j <= b.size(); ++j)
		{
			const bool same = std::tolower((unsigned char)a[i-1]) == std::tolower((unsigned char)b[j-1]);
			cur[j] = std::min(std::min(prev[j] + 1, cur[j-1] + 1), prev[j-1] + (same ? 0 : 1));
		}
		prev.swap(cur);
	}
	return prev[b.size()];
}

std::string Suggest(const std::string &name, const std::set<std::string> &known)
{
	std::string best;
	size_t bestd = (size_t)-1;
	for(std::set<std::string>::const_iterator it = known.begin(); it != known.end(); ++it)
	{
		const size_t d = Distance(name, *it);
		if(d < bestd) { bestd = d; best = *it; }
	}
	//  Only a near miss is worth naming: a third of the length, at least one.
	const size_t limit = std::max<size_t>(1, name.size() / 3);
	return (bestd <= limit) ? best : std::string();
}

} // namespace

void Block::collectAsked(const std::string &prefix,
                         std::map<std::string, std::set<std::string> > &asked) const
{
	const std::string here = prefix.empty() ? std::string() : prefix + ".";
	for(std::set<std::string>::const_iterator it = _askedVars.begin(); it != _askedVars.end(); ++it)
		asked[*it].insert(prefix.empty() ? std::string("(top level)") : "[" + prefix + "]");
	for(std::map<std::string, std::vector<Block*> >::const_iterator bt = _blocks.begin(); bt != _blocks.end(); ++bt)
		for(size_t i = 0; i < bt->second.size(); ++i)
			bt->second[i]->collectAsked(here + bt->first, asked);
}

void Block::collectUnknown(const std::string &prefix, std::vector<std::string> &out) const
{
	std::map<std::string, std::set<std::string> > asked;
	collectAsked(prefix, asked);
	collectUnknown(prefix, out, asked);
}

void Block::collectUnknown(const std::string &prefix, std::vector<std::string> &out,
                           const std::map<std::string, std::set<std::string> > &asked) const
{
	const std::string here = prefix.empty() ? std::string() : prefix + ".";
	for(std::map<std::string, Value>::const_iterator it = _vars.begin(); it != _vars.end(); ++it)
	{
		if(_askedVars.count(it->first)) continue;
		std::string line = "key   " + here + it->first;
		const std::string s = Suggest(it->first, _askedVars);
		std::map<std::string, std::set<std::string> >::const_iterator where = asked.find(it->first);
		if(!s.empty()) line += "   (did you mean " + here + s + "?)";
		else if(where != asked.end())
		{
			//  The right name in the wrong block, PrecursorRatio under
			//  [Initial] rather than [Initial.SurfaceShape]: say where it goes.
			line += "   (this key is read in ";
			for(std::set<std::string>::const_iterator w = where->second.begin(); w != where->second.end(); ++w)
				line += (w == where->second.begin() ? "" : ", ") + *w;
			line += ")";
		}
		out.push_back(line);
	}
	for(size_t k = 0; k < _blockOrder.size(); ++k)
	{
		const std::string &bn = _blockOrder[k];
		std::map<std::string, std::vector<Block*> >::const_iterator bt = _blocks.find(bn);
		if(bt == _blocks.end()) continue;
		if(!_askedBlocks.count(bn))
		{
			std::string line = "block [" + here + bn + "]";
			const std::string s = Suggest(bn, _askedBlocks);
			if(!s.empty()) line += "   (did you mean [" + here + s + "]?)";
			out.push_back(line);
			continue;
		}
		for(size_t i = 0; i < bt->second.size(); ++i)
			bt->second[i]->collectUnknown(here + bn, out, asked);
	}
}

//------------------------------------------------------------------
//  Type conversion
//
//  Each conversion takes the whole text or refuses it. What strtod makes of
//  a text is not enough: "fast" for a viscosity would be 0, "1e-3x" 1e-3,
//  "nan" would pass every range check written as a comparison, a NaN for an
//  int would be cast to int (undefined behavior), "2.7" for Nx would be 3,
//  and any word that is not one of the spellings of true would be false --
//  all without a word (yuzen/tests/fuzz feeds the reader such texts).
//------------------------------------------------------------------
namespace {

std::string trimmed(const std::string &s)
{
	size_t a = 0, b = s.size();
	while(a < b && std::isspace((unsigned char)s[a])) ++a;
	while(b > a && std::isspace((unsigned char)s[b-1])) --b;
	return s.substr(a, b - a);
}

std::string quoted(const std::string &s)
{
	return "`" + (s.size() > 60 ? s.substr(0, 60) + "..." : s) + "'";
}

bool toBool(const std::string &s, bool &out)
{
	std::string t;
	const std::string u = trimmed(s);
	for(size_t i = 0; i < u.size(); ++i) t += (char)std::tolower((unsigned char)u[i]);
	if(t == "true" || t == "yes" || t == "on" || t == "1" || t == ".true.")   { out = true;  return true; }
	if(t == "false" || t == "no" || t == "off" || t == "0" || t == ".false.") { out = false; return true; }
	return false;
}

bool toDouble(const std::string &s, double &out)
{
	const std::string t = trimmed(s);
	if(t.empty()) return false;
	char *end = 0;
	const double d = std::strtod(t.c_str(), &end);
	if(end != t.c_str() + t.size() || !std::isfinite(d)) return false;
	out = d;
	return true;
}

bool toInt(const std::string &s, int &out)
{
	// Read via a double so that exponent notation such as "1e3" is still
	// picked up as an integer
	double d = 0;
	if(!toDouble(s, d) || d != std::floor(d)) return false;
	//  Cutting a value that does not fit in an int down with (int) is undefined
	//  behavior, and in practice MaxSteps = 4000000000 (> INT_MAX) came out
	//  negative. Saturate instead.
	if(d >= (double)std::numeric_limits<int>::max()) out = std::numeric_limits<int>::max();
	else if(d <= (double)std::numeric_limits<int>::min()) out = std::numeric_limits<int>::min();
	else out = (int)d;
	return true;
}

//  A count that may exceed an int (Numerics.MaxSteps). Read via a double
//  like toInt, so whole numbers are exact up to 2^53; saturated beyond the
//  range of long long.
bool toInt64(const std::string &s, long long &out)
{
	double d = 0;
	if(!toDouble(s, d) || d != std::floor(d)) return false;
	if(d >= 9.2e18) out = std::numeric_limits<long long>::max();
	else if(d <= -9.2e18) out = std::numeric_limits<long long>::min();
	else out = (long long)d;
	return true;
}

const char *kindOf(const bool *)   { return "true or false"; }
const char *kindOf(const long long *) { return "a whole number"; }
const char *kindOf(const int *)    { return "a whole number"; }
const char *kindOf(const double *) { return "a finite number"; }

} // namespace

const std::string& Block::scalarItem(const std::string &name, const Value &v) const
{
	if(v.items.size() != 1)
		throw VariableNotFound("Variable [" + keyPath(name) + "]: "
		                       + (v.items.empty() ? std::string("an empty value or list")
		                                          : std::string("a list of several values"))
		                       + " where one value is wanted.\n");
	return v.items[0];
}

#define CONFIG_SCALAR_GETTER(TYPE, CONV)                                            \
	void Block::getin(const std::string &name, TYPE &v) const                        \
	{                                                                                \
		if(!find(name)) throw VariableNotFound("Variable not found: [" + name + "]"); \
		TYPE r = TYPE();                                                             \
		getin(name, r, r);                                                           \
		v = r;                                                                       \
	}                                                                                \
	void Block::getin(const std::string &name, TYPE &v, TYPE def) const               \
	{                                                                                \
		const Value *p = find(name);                                                 \
		if(!p) { v = def; return; }                                                  \
		const std::string &item = scalarItem(name, *p);                              \
		TYPE r = TYPE();                                                             \
		if(!CONV(item, r))                                                           \
			throw VariableNotFound("Variable [" + keyPath(name) + "] = " + quoted(item) \
			                       + " is not " + kindOf((const TYPE*)0) + ".\n");   \
		v = r;                                                                       \
	}

CONFIG_SCALAR_GETTER(int,    toInt)
typedef long long long_long;
CONFIG_SCALAR_GETTER(long_long, toInt64)
CONFIG_SCALAR_GETTER(double, toDouble)
CONFIG_SCALAR_GETTER(bool,   toBool)
#undef CONFIG_SCALAR_GETTER

void Block::getin(const std::string &name, std::string &v) const
{
	const Value *p = find(name);
	if(!p) throw VariableNotFound("Variable not found: [" + name + "]");
	v = scalarItem(name, *p);
}
void Block::getin(const std::string &name, std::string &v, const std::string &def) const
{
	const Value *p = find(name);
	v = p ? scalarItem(name, *p) : def;
}

#define CONFIG_VECTOR_GETTER(TYPE, CONV)                                             \
	void Block::getin(const std::string &name, std::vector<TYPE> &v) const            \
	{                                                                                 \
		if(!find(name)) throw VariableNotFound("Variable not found: [" + name + "]"); \
		getin(name, v, v);                                                            \
	}                                                                                 \
	void Block::getin(const std::string &name, std::vector<TYPE> &v,                  \
	                  const std::vector<TYPE> &def) const                             \
	{                                                                                 \
		const Value *p = find(name);                                                  \
		if(!p) { v = def; return; }                                                   \
		std::vector<TYPE> r;                                                          \
		for(size_t i = 0; i < p->items.size(); ++i)                                   \
		{                                                                             \
			TYPE x = TYPE();                                                          \
			if(!CONV(p->items[i], x))                                                 \
				throw VariableNotFound("Variable [" + keyPath(name) + "]: item "      \
				                       + std::to_string(i + 1) + " = " + quoted(p->items[i]) \
				                       + " is not " + kindOf((const TYPE*)0) + ".\n"); \
			r.push_back(x);                                                           \
		}                                                                             \
		v = r;                                                                        \
	}

CONFIG_VECTOR_GETTER(double, toDouble)
CONFIG_VECTOR_GETTER(int,    toInt)
#undef CONFIG_VECTOR_GETTER

void Block::getin(const std::string &name, std::vector<std::string> &v) const
{
	const Value *p = find(name);
	if(!p) throw VariableNotFound("Variable not found: [" + name + "]");
	v = p->items;
}
void Block::getin(const std::string &name, std::vector<std::string> &v,
                  const std::vector<std::string> &def) const
{
	const Value *p = find(name);
	v = p ? p->items : def;
}


//==================================================================
//  The .inp parser
//
//  Grammar:
//    comment        // ... to end of line   /  # ... to end of line
//    block          block Name  ...  end     (case insensitive)
//    variable       Name = value
//    value          number / "string" / true,false / ( a, b, ... )
//==================================================================
namespace {

struct Lexer
{
	const std::string &src;
	size_t pos;
	int    line;
	bool   openString;		// the last string ran to the end of the file
	explicit Lexer(const std::string &s) : src(s), pos(0), line(1), openString(false) {}

	void skip()
	{
		for(;;)
		{
			while(pos < src.size() && std::isspace((unsigned char)src[pos]))
			{
				if(src[pos] == '\n') ++line;
				++pos;
			}
			if(pos + 1 < src.size() && src[pos] == '/' && src[pos+1] == '/')
			{
				while(pos < src.size() && src[pos] != '\n') ++pos;
				continue;
			}
			if(pos < src.size() && src[pos] == '#')
			{
				while(pos < src.size() && src[pos] != '\n') ++pos;
				continue;
			}
			break;
		}
	}
	bool eof() { skip(); return pos >= src.size(); }
	char peek() { skip(); return pos < src.size() ? src[pos] : '\0'; }

	// Returns one identifier, number or punctuation mark. A string comes back
	// with its quotes stripped.
	std::string next(bool *wasString = 0)
	{
		if(wasString) *wasString = false;
		skip();
		if(pos >= src.size()) return std::string();

		const char c = src[pos];
		if(c == '"' || c == '\'')
		{
			const char q = c;
			++pos;
			std::string out;
			while(pos < src.size() && src[pos] != q)
			{
				if(src[pos] == '\\' && pos + 1 < src.size()) ++pos;
				if(src[pos] == '\n') ++line;
				out += src[pos++];
			}
			if(pos < src.size()) ++pos;		// the closing quote
			else openString = true;
			if(wasString) *wasString = true;
			return out;
		}
		if(c == '(' || c == ')' || c == ',' || c == '=' || c == '[' || c == ']')
		{
			++pos;
			return std::string(1, c);
		}
		std::string out;
		while(pos < src.size())
		{
			const char d = src[pos];
			if(std::isspace((unsigned char)d) || d == '(' || d == ')' || d == ','
			   || d == '=' || d == '"' || d == '\'' || d == '[' || d == ']') break;
			if(d == '/' && pos + 1 < src.size() && src[pos+1] == '/') break;
			if(d == '#') break;
			out += d;
			++pos;
		}
		return out;
	}
};

std::string lower(const std::string &s)
{
	std::string t;
	for(size_t i = 0; i < s.size(); ++i) t += (char)std::tolower((unsigned char)s[i]);
	return t;
}

//  The deepest nesting of blocks either reader builds. The input itself
//  goes four deep; the parsers recurse once per level, so without a limit a
//  file of a hundred thousand nested blocks would run them out of stack.
const int kMaxDepth = 64;

void parseInpBlock(Lexer &lx, Block *blk, const std::string &file, int depth = 0)
{
	if(depth > kMaxDepth)
	{
		std::ostringstream oss;
		oss << file << ":" << lx.line << ": blocks nested more than " << kMaxDepth << " deep";
		throw BlockNotFound(oss.str());
	}
	for(;;)
	{
		if(lx.eof()) return;					// end of the root
		const size_t save = lx.pos;
		bool wasStr = false;
		const std::string tok = lx.next(&wasStr);
		if(tok.empty()) return;

		const std::string lt = lower(tok);
		if(lt == "end")  return;
		if(lt == "block" || lt == "root")
		{
			const std::string bname = lx.next();
			if(bname.empty())
			{
				std::ostringstream oss;
				oss << file << ":" << lx.line << ": block name expected";
				throw BlockNotFound(oss.str());
			}
			parseInpBlock(lx, blk->addBlock(bname), file, depth + 1);
			continue;
		}

		// A variable: Name = value
		const std::string eq = lx.next();
		if(eq != "=")
		{
			std::ostringstream oss;
			oss << file << ":" << lx.line << ": '=' expected after '" << tok << "'";
			throw BlockNotFound(oss.str());
		}
		(void)save;

		Value val;
		const char c = lx.peek();
		const int vline = lx.line;				// where the value starts
		if(c == '(' || c == '[')
		{
			const char close = (c == '(') ? ')' : ']';
			lx.next();							// the opening bracket
			val.isVector = true;
			for(;;)
			{
				const char p = lx.peek();
				if(p == close) { lx.next(); break; }
				if(p == '\0')
				{
					std::ostringstream oss;
					oss << file << ":" << lx.line << ": unterminated vector for '" << tok << "'";
					throw BlockNotFound(oss.str());
				}
				bool q = false;
				const std::string item = lx.next(&q);
				if(item == ",") continue;
				if(item == std::string(1, close)) break;
				val.items.push_back(item);
				if(q) val.isString = true;
			}
		}
		else
		{
			bool q = false;
			val.items.push_back(lx.next(&q));
			val.isString = q;
		}
		//  A string without its closing quote would take the rest of the
		//  file as its value, and every key after it would keep its default
		if(lx.openString)
		{
			std::ostringstream oss;
			oss << file << ":" << vline << ": unterminated string for '" << tok << "'";
			throw BlockNotFound(oss.str());
		}
		blk->setValue(tok, val);
	}
}

Block* readInp(const std::string &filename)
{
	std::ifstream fin(filename.c_str(), std::ios::binary);
	if(!fin) throw BlockNotFound("Cannot open input file: " + filename);
	std::ostringstream ss;
	ss << fin.rdbuf();
	const std::string src = ss.str();

	Block *root = new Block();
	root->setName("root");
	try
	{
		Lexer lx(src);
		parseInpBlock(lx, root, filename);
	}
	catch(...)
	{
		delete root;
		throw;
	}
	return root;
}


//------------------------------------------------------------------
//  The .toml loader (toml++)
//------------------------------------------------------------------
void fillFromToml(const toml::table &tbl, Block *blk, int depth = 0);

//  A value of a kind the input has no use for (a date or a time, an array
//  inside an array) is refused: dropped without a word, the key would be
//  neither read nor reported as unknown, and its setting would keep the
//  default.
[[noreturn]] void unsupported(const Block *blk, const std::string &key, const char *what)
{
	const std::string path = blk->path().empty() ? key : blk->path() + "." + key;
	throw BlockNotFound("Variable [" + path + "]: " + what + " is not a value Yuzen reads"
	                    " (numbers, strings, booleans and arrays of them are).\n");
}

void setScalar(Block *blk, const std::string &key, const toml::node &n)
{
	Value v;
	if(const toml::value<std::string> *s = n.as_string())
	{
		v.items.push_back(s->get());
		v.isString = true;
	}
	else if(const toml::value<bool> *b = n.as_boolean())
	{
		v.items.push_back(b->get() ? "true" : "false");
	}
	else if(const toml::value<int64_t> *i = n.as_integer())
	{
		std::ostringstream oss; oss << i->get();
		v.items.push_back(oss.str());
	}
	else if(const toml::value<double> *d = n.as_floating_point())
	{
		std::ostringstream oss; oss << std::setprecision(17) << d->get();
		v.items.push_back(oss.str());
	}
	else unsupported(blk, key, "a date or time");
	blk->setValue(key, v);
}

void setArray(Block *blk, const std::string &key, const toml::array &arr)
{
	Value v;
	v.isVector = true;
	for(size_t i = 0; i < arr.size(); ++i)
	{
		const toml::node *n = arr.get(i);
		if(!n) continue;
		if(const toml::value<std::string> *s = n->as_string())
		{
			v.items.push_back(s->get()); v.isString = true;
		}
		else if(const toml::value<bool> *b = n->as_boolean())
		{
			v.items.push_back(b->get() ? "true" : "false");
		}
		else if(const toml::value<int64_t> *iv = n->as_integer())
		{
			std::ostringstream oss; oss << iv->get(); v.items.push_back(oss.str());
		}
		else if(const toml::value<double> *d = n->as_floating_point())
		{
			std::ostringstream oss; oss << std::setprecision(17) << d->get();
			v.items.push_back(oss.str());
		}
		else unsupported(blk, key, n->is_table() ? "a table inside an array of values"
		                          : n->is_array() ? "an array inside an array" : "a date or time");
	}
	blk->setValue(key, v);
}

void fillFromToml(const toml::table &tbl, Block *blk, int depth)
{
	if(depth > kMaxDepth)
		throw BlockNotFound("[" + blk->path() + "]: tables nested more than "
		                    + std::to_string(kMaxDepth) + " deep\n");
	for(const auto &kv : tbl)
	{
		const std::string key(kv.first.str());
		const toml::node &n = kv.second;

		if(const toml::table *sub = n.as_table())
		{
			fillFromToml(*sub, blk->addBlock(key), depth + 1);
		}
		else if(const toml::array *arr = n.as_array())
		{
			if(arr->is_array_of_tables())
			{
				// A repeated [[Name]] becomes several blocks of the same name
				for(size_t i = 0; i < arr->size(); ++i)
				{
					if(const toml::table *t = arr->get(i)->as_table())
						fillFromToml(*t, blk->addBlock(key), depth + 1);
				}
			}
			else setArray(blk, key, *arr);
		}
		else setScalar(blk, key, n);
	}
}

Block* readToml(const std::string &filename)
{
	Block *root = new Block();
	root->setName("root");
	try
	{
		const toml::table tbl = toml::parse_file(filename);
		fillFromToml(tbl, root);
	}
	catch(const toml::parse_error &e)
	{
		delete root;
		std::ostringstream oss;
		oss << "TOML parse error in " << filename << ": " << e.description();
		throw BlockNotFound(oss.str());
	}
	catch(...)
	{
		delete root;
		throw;
	}
	return root;
}

} // namespace


Block* Block::newread(const std::string &filename)
{
	const size_t dot = filename.find_last_of('.');
	std::string ext = (dot == std::string::npos) ? "" : lower(filename.substr(dot));
	if(ext == ".toml") return readToml(filename);
	return readInp(filename);
}


//==================================================================
//  Writing out in TOML format
//==================================================================
namespace {

bool looksNumeric(const std::string &s)
{
	if(s.empty()) return false;
	char *end = 0;
	std::strtod(s.c_str(), &end);
	return end && *end == '\0';
}

void emitValue(std::ostream &os, const Value &v)
{
	const bool quote = v.isString || (!v.items.empty() && !looksNumeric(v.items[0])
	                                  && lower(v.items[0]) != "true" && lower(v.items[0]) != "false");
	if(v.isVector) os << "[";
	for(size_t i = 0; i < v.items.size(); ++i)
	{
		if(i) os << ", ";
		if(quote) os << '"' << v.items[i] << '"';
		else      os << v.items[i];
	}
	if(v.isVector) os << "]";
}

} // namespace

void Block::writeTo(std::ostream &os, const std::string &prefix) const
{
	for(std::map<std::string, Value>::const_iterator it = _vars.begin(); it != _vars.end(); ++it)
	{
		os << it->first << " = ";
		emitValue(os, it->second);
		os << "\n";
	}
	for(size_t i = 0; i < _blockOrder.size(); ++i)
	{
		const std::string &bn = _blockOrder[i];
		std::map<std::string, std::vector<Block*> >::const_iterator bt = _blocks.find(bn);
		if(bt == _blocks.end()) continue;
		const std::string path = prefix.empty() ? bn : (prefix + "." + bn);
		const bool repeated = bt->second.size() > 1;
		for(size_t k = 0; k < bt->second.size(); ++k)
		{
			os << "\n" << (repeated ? "[[" : "[") << path << (repeated ? "]]" : "]") << "\n";
			bt->second[k]->writeTo(os, path);
		}
	}
}

void Block::write(const std::string &filename) const
{
	std::ofstream fout(filename.c_str());
	if(!fout) return;
	fout << "# Yuzen input (TOML) - auto generated\n";
	fout << std::setprecision(17);
	writeTo(fout, "");
}

std::string Block::libraryName()    { return "ConfigBlock"; }
std::string Block::libraryVersion() { return "1.0 (.inp + toml++)"; }

} }  // namespace yuzen::config
