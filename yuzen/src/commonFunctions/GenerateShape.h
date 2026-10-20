/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  GenerateShape.h : assembles the shape settings of the input into a set of shape functions
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _GENERATE_SHAPE_INCLUDED_
#define _GENERATE_SHAPE_INCLUDED_

#include <utility>
#include <vector>
#include <list>

#include "Input.h"
#include "Grids.h"
#include "ShapeFunctions.h"


#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

///////////
//--------------------------------------------------------------------------
//  Assembles the shape settings listed in the input into a sequence of shape
//  functions.
//
//  The input describes a shape by listing blocks. There may be several blocks
//  of the same kind, so a single shape is named by a pair (tShapeSet) of "the
//  kind" and "which one of that kind it is". The order of application follows
//  from the Order of each block.
//
//  Evaluating the assembled sequence in order and combining the results gives
//  the initial film thickness (or the substrate topography). Writing that into
//  the field itself is done by C_InitialSurface.
//--------------------------------------------------------------------------
class C_GenerateShape
{
public:
	//  The kinds of shape available. Each corresponds to a block name in the input.
	enum tShapeFunc
	{
		RectBump = 0,		// rectangular bump (the side walls can be rounded)
		EllipseBump = 1,	// elliptical bump
		FlatPlane = 2,		// tilted plane
		MexicanHat = 3,		// Mexican hat
		SmoothConvex = 4,	// smooth convex cap
		PeriodicWave = 5,	// periodic wave
		EdgeBump = 6		// raised rim along an edge
	};

	//  A pair naming one shape (its kind, plus its number within that kind)
	class tShapeSet
	{
	public:
		tShapeSet() : _tShapeFunc(RectBump), _index(0) {}


		/**
		 * @brief Names one shape.
		 * @param tt  which kind of shape it is
		 * @param ii  position within the shapes of that kind, counted in
		 *           input-file order from 0
		 */
		tShapeSet(tShapeFunc tt, int ii) : _tShapeFunc(tt), _index(ii) {}
		// Which kind of shape, and which one of that kind. Both are readable
		// and writable; the non-const forms hand out a reference to the member.
		//  which kind of shape
		const tShapeFunc& Type() const {return _tShapeFunc;}
		      tShapeFunc& Type()       {return _tShapeFunc;}
		//  number within that kind (in the order listed in the input)
		const int& Index() const {return _index;}
		      int& Index()       {return _index;}

	private:
		tShapeFunc _tShapeFunc;
		int _index;

		// The default copy constructor and assignment operator are allowed
		//tShapeSet(const tShapeSet&);
		//void operator=(const tShapeSet&);

	};

public:
	C_GenerateShape() : _Xmax(1), _Xmin(0), _Ymax(1), _Ymin(0)
	{
		_shapeOrder.resize(0);
		_vecRectBump.resize(0);
		_vecEllipseBump.resize(0);
		_vecFlatPlane.resize(0);
		_vecMexicanHat.resize(0);
		_vecSmoothConvex.resize(0);
		_vecPeriodicWave.resize(0);
		_vecEdgeBump.resize(0);
	};

	/**
	 * @brief Reads the shape blocks of the input and builds the list of shape
	 *        functions, sorted into the order they are to be applied.
	 *
	 * Overloaded twice over: once per grid dimension, and once per input type,
	 * because the initial film thickness and the substrate topography are
	 * described by different blocks even though the shapes are the same.
	 * Calling Build again replaces the previous list.
	 *
	 * @param g      grid, used only for the extent of the domain
	 * @param scale  in-plane length scale. The lengths in the input file are
	 *               dimensional and are divided by this on the way in, so the
	 *               shape functions work in non-dimensional coordinates.
	 * @param inp    the input block the shapes are read from: the initial film
	 *               thickness block, or the substrate block
	 */
	//--------------------------------------------------------------
	//  Reads the shape settings of the input and builds the sequence of
	//  shape functions, arranged in the order they are applied.
	//  scale is the in-plane length scale (the input is written dimensionally).
	//  Overloaded because the film thickness and the substrate use
	//  different input types.
	//--------------------------------------------------------------
	void Build(const C_Grid1D &g, double scale, const Input::I_Initial::I_SurfaceShape &inp);
	void Build(const C_Grid1D &g, double scale, const Input::I_Substrate &inp);

	void Build(const C_Grid2D &g, double scale, const Input::I_Initial::I_SurfaceShape &inp);
	void Build(const C_Grid2D &g, double scale, const Input::I_Substrate &inp);

	/**
	 * @brief Height at an in-plane point, with every shape applied in order.
	 *
	 * Each shape is evaluated in turn and merged into the running value the way
	 * its CreateType says: Add keeps whichever of the two is larger, Remove
	 * whichever is smaller, and Stack adds the shape on top of what is there.
	 * Order therefore matters. With no shapes built the result is 0.
	 *
	 * @param px  in-plane coordinate, non-dimensional. This one-argument form
	 *            is the 1-D case and evaluates at py = 0.
	 * @return the combined height, non-dimensional.
	 */
	double Evaluate(double px);

	/**
	 * @brief Height at an in-plane point, with every shape applied in order.
	 *
	 * The two-dimensional form of Evaluate(double); the shapes are combined
	 * in exactly the same way.
	 *
	 * @param px  first in-plane coordinate, non-dimensional
	 * @param py  second in-plane coordinate, non-dimensional
	 * @return the combined height, non-dimensional.
	 */
	double Evaluate(double px, double py);


private:
	double _Xmax, _Xmin, _Ymax, _Ymin;

	std::vector<tShapeSet> _shapeOrder;

	std::vector<C_RectBump>		_vecRectBump;
	std::vector<C_EllipseBump>	_vecEllipseBump;
	std::vector<C_FlatPlane>	_vecFlatPlane;
	std::vector<C_MexicanHat>	_vecMexicanHat;
	std::vector<C_SmoothConvex>	_vecSmoothConvex;
	std::vector<C_PeriodicWave>	_vecPeriodicWave;
	std::vector<C_EdgeBump>		_vecEdgeBump;

	void SortOrder();
	C_ShapeFunc_base::tCreateType GetType(std::string type);

	void Build(double scale, const Input::I_Initial::I_SurfaceShape &inp);
	void Build(double scale, const Input::I_Substrate &inp);

	int GetOrder(tShapeSet);

private:
	C_GenerateShape(const C_GenerateShape&);
	void operator=(const C_GenerateShape&);

};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_GENERATE_SHAPE_INCLUDED_
