/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ShapeFunctions.h : shape functions - rectangular bump, elliptical bump, plane and so on
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _SHAPE_FUNCTIONS_INCLUDED_
#define _SHAPE_FUNCTIONS_INCLUDED_

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif



//-------------------
//   Common part
//-------------------
//--------------------------------------------------------------------------
//  Base of the shape functions.
//
//  The initial film thickness and the substrate topography are built by
//  listing these and superposing them. Each one on its own is no more than
//  "a function that returns a height given an in-plane coordinate (px, py)".
//
//  How they are combined is decided by tCreateType:
//    Add    add to the value already there
//    Remove subtract from the value already there (cut a groove)
//    Stack  stack on top of the value already there (take the underlying
//           height as the reference)
//  They are applied in the order they are listed, so changing the order
//  changes the result.
//--------------------------------------------------------------------------
class C_ShapeFunc_base
{
public:
	//  how this shape is combined with what is already there
	enum tCreateType
	{
		Add,	// add
		Remove,	// subtract
		Stack	// stack on top of what is underneath
	};
	C_ShapeFunc_base() : _order(0), _CreateType(Add), _base(0), _amp(1), _cntX(0), _cntY(0) {}
	//  The shapes are held by value in per-type vectors, so nothing is deleted
	//  through this base today. The class has virtual functions all the same,
	//  so give it a virtual destructor: it costs nothing here and keeps a
	//  future base-pointer delete from being undefined behavior.
	virtual ~C_ShapeFunc_base() {}

	// Application order and combination rule, both fixed by Build. C_GenerateShape
	// sorts the shapes by Order and then merges each one according to Type.
	//  the order in which this shape is applied (the Order of the input)
	int Order(){return _order;}
	//  how this shape is combined with what is already there
	tCreateType Type(){return _CreateType;}

protected:
	int  _order;			// order of application
	tCreateType _CreateType;	// how this shape is combined
	double _base;			// pedestal (the height this shape holds at minimum)
	double _amp;			// amplitude
	double _cntX, _cntY;	// in-plane coordinates of the center

	//  height of this shape at the in-plane coordinate (px, py)
	virtual double operator()(double px, double py) const = 0;

private:
	// The default copy constructor and assignment operator are allowed
	//C_ShapeFunc_base(const C_ShapeFunc_base&);
	//void operator=(const C_ShapeFunc_base&);
};

//-------------------
//   RectBump
//-------------------
class C_RectBump : public C_ShapeFunc_base
{
public:
	C_RectBump() : C_ShapeFunc_base(), _widthX(1), _widthY(1), _slopeX(0.1), _slopeY(0.1) {}
	/**
	 * @brief Sets up a rectangular bump with rounded side walls.
	 *
	 * The first six arguments are the ones every shape takes: order is the
	 * position in the application order, type the combination rule, base the
	 * height the shape holds away from the feature, amp the height of the
	 * feature above base, and (cntX, cntY) its center. All lengths are
	 * non-dimensional.
	 *
	 * @param order  position of this shape in the application order
	 * @param type   rule by which this shape is combined with the ones already applied
	 * @param base   height the shape holds away from the feature
	 * @param amp    height of the feature above base
	 * @param cntX   x coordinate of the center of the feature
	 * @param cntY   y coordinate of the center of the feature
	 * @param widthX  width of the bump in x, measured between the mid-points of
	 *                the two side walls
	 * @param widthY  the same width in y
	 * @param slopeX  width of the rounded side wall in x; larger values give a
	 *                more gradual wall, and a value approaching zero gives a step
	 * @param slopeY  the same wall width in y
	 */
	void Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double widthX, double widthY, double slopeX, double slopeY);
	/** @brief Height of this shape at the in-plane point (px, py). */
	virtual double operator()(double px, double py) const;

private:
	double _widthX, _widthY;
	double _slopeX, _slopeY;

	// The default copy constructor and assignment operator are allowed
	//	C_RectBump(const C_RectBump&);
	//void operator=(const C_RectBump&);
};

//-------------------
//   EllipseBump
//-------------------
class C_EllipseBump : public C_ShapeFunc_base
{
public:
	C_EllipseBump() : C_ShapeFunc_base(), _widthX(1), _widthY(1), _slope(0.1) {}
	/**
	 * @brief Sets up an elliptical bump with a rounded rim.
	 *
	 * The first six arguments are the common ones, as for C_RectBump.
	 *
	 * @param order  position of this shape in the application order
	 * @param type   rule by which this shape is combined with the ones already applied
	 * @param base   height the shape holds away from the feature
	 * @param amp    height of the feature above base
	 * @param cntX   x coordinate of the center of the feature
	 * @param cntY   y coordinate of the center of the feature
	 * @param widthX  semi-axis of the ellipse in x; the rim sits where
	 *                ((px-cntX)/widthX)^2 + ((py-cntY)/widthY)^2 = 1/4
	 * @param widthY  the same semi-axis in y
	 * @param slope   width of the rounded rim, in units of that same radius
	 */
	void Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double widthX, double widthY, double slope);
	/** @brief Height of this shape at the in-plane point (px, py). */
	virtual double operator()(double px, double py) const;

private:
	double _widthX, _widthY;
	double _slope;

	// The default copy constructor and assignment operator are allowed
	//C_EllipseBump(const C_EllipseBump&);
	//void operator=(const C_EllipseBump&);
};

//----------------------
//     FlatPlane
//----------------------
class C_FlatPlane : public C_ShapeFunc_base
{
public:
	C_FlatPlane() : C_ShapeFunc_base(), _angleX(0), _angleY(0) {}
	/**
	 * @brief Sets up a tilted plane.
	 *
	 * The first six arguments are the common ones, as for C_RectBump, except
	 * that here amp is the height of the plane at (cntX, cntY) and base takes
	 * no part in the result.
	 *
	 * @param order  position of this shape in the application order
	 * @param type   rule by which this shape is combined with the ones already applied
	 * @param base   height the shape holds away from the feature
	 * @param amp    height of the feature above base
	 * @param cntX   x coordinate of the center of the feature
	 * @param cntY   y coordinate of the center of the feature
	 * @param angleX  tilt in the x direction, in radians; the height changes by
	 *                tan(angleX) per unit of px
	 * @param angleY  the same tilt in y
	 */
	void Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double angleX, double angleY);
	/** @brief Height of this shape at the in-plane point (px, py). */
	virtual double operator()(double px, double py) const;

private:
	double _angleX, _angleY;

	// The default copy constructor and assignment operator are allowed
	//C_FlatPlane(const C_FlatPlane&);
	//void operator=(const C_FlatPlane&);
};

//---------------------------------------------
//   Mexican-hat (1-(x/L)^2) * exp(-b*(x/L)^2)
//---------------------------------------------
class C_MexicanHat : public C_ShapeFunc_base
{
public:
	C_MexicanHat(): C_ShapeFunc_base(), _Lx(1), _Ly(1), _beta(0.1) {}
	/**
	 * @brief Sets up a Mexican-hat profile: a central peak ringed by a trough.
	 *
	 * The first six arguments are the common ones, as for C_RectBump.
	 *
	 * @param order  position of this shape in the application order
	 * @param type   rule by which this shape is combined with the ones already applied
	 * @param base   height the shape holds away from the feature
	 * @param amp    height of the feature above base
	 * @param cntX   x coordinate of the center of the feature
	 * @param cntY   y coordinate of the center of the feature
	 * @param Lx    length scale in x. The profile crosses zero at a scaled radius
	 *              of 1, so this sets where the trough starts. Zero in either
	 *              direction makes the shape evaluate to 0.
	 * @param Ly    the same length scale in y
	 * @param beta  decay of the Gaussian envelope; larger values pull the
	 *              trough in and make it shallower
	 */
	void Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double Lx, double Ly, double beta);
	/** @brief Height of this shape at the in-plane point (px, py). */
	virtual double operator()(double px, double py) const;

private:
	double _Lx, _Ly;
	double _beta;

	// The default copy constructor and assignment operator are allowed
	//C_MexicanHat(const C_MexicanHat&);
	//void operator=(const C_MexicanHat&);
};

//--------------------------------
// Smooth convex (1-(x/L)^2)^2
//--------------------------------
class C_SmoothConvex : public C_ShapeFunc_base
{
public:
	C_SmoothConvex(): C_ShapeFunc_base(), _Lx(1), _Ly(1) {}
	/**
	 * @brief Sets up a smooth convex cap of finite support.
	 *
	 * The first six arguments are the common ones, as for C_RectBump. Outside
	 * the ellipse the shape is exactly base, and it meets that level with zero
	 * slope, so there is no kink at the edge.
	 *
	 * @param order  position of this shape in the application order
	 * @param type   rule by which this shape is combined with the ones already applied
	 * @param base   height the shape holds away from the feature
	 * @param amp    height of the feature above base
	 * @param cntX   x coordinate of the center of the feature
	 * @param cntY   y coordinate of the center of the feature
	 * @param Lx  semi-axis in x of the ellipse the cap sits on
	 * @param Ly  the same semi-axis in y
	 */
	void Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double Lx, double Ly);
	/** @brief Height of this shape at the in-plane point (px, py). */
	virtual double operator()(double px, double py) const;

private:
	double _Lx, _Ly;

	// The default copy constructor and assignment operator are allowed
	//C_SmoothConvex(const C_SmoothConvex&);
	//void operator=(const C_SmoothConvex&);
};

//--------------------------------
//   PeriodicWave
//--------------------------------
class C_PeriodicWave : public C_ShapeFunc_base
{
public:
	C_PeriodicWave(): C_ShapeFunc_base(), _Lx(1), _Ly(1) {}
	/**
	 * @brief Sets up a cosine corrugation, periodic over the whole domain.
	 *
	 * The first six arguments are the common ones, as for C_RectBump; here
	 * (cntX, cntY) is the point where the wave is at its crest, which is how
	 * the phase of the input is expressed.
	 *
	 * @param order  position of this shape in the application order
	 * @param type   rule by which this shape is combined with the ones already applied
	 * @param base   height the shape holds away from the feature
	 * @param amp    height of the feature above base
	 * @param cntX   x coordinate of the center of the feature
	 * @param cntY   y coordinate of the center of the feature
	 * @param Lx  wavelength in x
	 * @param Ly  wavelength in y
	 */
	void Build(int order, tCreateType type, double base, double amp, double cntX, double cntY,
		double Lx, double Ly);
	/** @brief Height of this shape at the in-plane point (px, py). */
	virtual double operator()(double px, double py) const;

private:
	double _Lx, _Ly;

	// The default copy constructor and assignment operator are allowed
	//C_PeriodicWave(const C_PeriodicWave&);
	//void operator=(const C_PeriodicWave&);
};

//--------------------------------
//   EdgeBump
//--------------------------------
class C_EdgeBump : public C_ShapeFunc_base
{
public:
	/** Which edge of the domain the ridge is attached to. */
	enum tLocation
	{
		Xmin, Xmax, Ymin, Ymax
	};

	C_EdgeBump(): C_ShapeFunc_base(), _alpha(0.5), _Lin(1), _Lout(1), _rs(0), _location(Xmax) {}
	/**
	 * @brief Sets up the raised rim that forms along a domain edge.
	 *
	 * The shape is one-dimensional: it varies only along the direction normal
	 * to the chosen edge, and is uniform along it. Two profiles are used, one
	 * on the inner side of rs and one on the outer, blended smoothly across it.
	 * Unlike the other shapes this one takes no center; rs plays that part.
	 *
	 * @param order  position in the application order (input Order)
	 * @param type   how this shape merges with what is already there
	 * @param base   height the shape holds away from the ridge
	 * @param amp    height at the crest of the ridge
	 * @param alpha  shape of the inner profile; larger values spread the decay
	 *               over a wider region
	 * @param Lin    length scale of the inner profile, toward the interior
	 * @param Lout   length scale of the outer profile, toward the edge. It is
	 *               also the distance from the edge at which the crest is put.
	 * @param rs     position of the crest along the normal direction, set by
	 *               the caller to Lout inside the chosen edge
	 * @param loc    which edge this ridge belongs to
	 */
	void Build(int order, tCreateType type, double base, double amp,
		double alpha, double Lin, double Lout, double rs, tLocation loc);
	/** @brief Height of this shape at the in-plane point (px, py). */
	virtual double operator()(double px, double py) const;


private:
	double _alpha, _Lin, _Lout, _rs;
	tLocation _location;

	// The default copy constructor and assignment operator are allowed
	//C_EdgeBump(const C_EdgeBump&);
	//void operator=(const C_EdgeBump&);
};



#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_SHAPE_FUNCTIONS_INCLUDED_
