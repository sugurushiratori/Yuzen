/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  GenerateShape.cpp : builds the initial film thickness and substrate shape from the input settings
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/


#include "GenerateShape.h"
#include "MyMath.h"

#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif


void C_GenerateShape::Build(const C_Grid1D &g, double scale, const Input::I_Initial::I_SurfaceShape &inp)
{
	_Xmin = g.XminReal() * scale;
	_Xmax = g.XmaxReal() * scale;
	_Ymin = 0;
	_Ymax = 1;

	Build(scale, inp);
}

void C_GenerateShape::Build(const C_Grid2D &g, double scale, const Input::I_Initial::I_SurfaceShape &inp)
{
	_Xmin = g.XminReal() * scale;
	_Xmax = g.XmaxReal() * scale;
	_Ymin = g.YminReal() * scale;
	_Ymax = g.YmaxReal() * scale;
	Build(scale, inp);
}

void C_GenerateShape::Build(double scale, const Input::I_Initial::I_SurfaceShape &inp)
{
	// RectBump
	for(unsigned int i = 0; i< inp.RectBump.size(); i++)
	{
		Input::I_Initial::I_SurfaceShape::I_RectBump i_obj = inp.RectBump[i];
		C_RectBump c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Height, i_obj.Center[0], i_obj.Center[1],
			i_obj.Width[0], i_obj.Width[1], i_obj.SlopeWidth[0], i_obj.SlopeWidth[1]);
		_vecRectBump.push_back(c_obj);
		tShapeSet p(RectBump, i);
		_shapeOrder.push_back(p);
	}

	// EllipseBump
	for(unsigned int i = 0; i< inp.EllipseBump.size(); i++)
	{
		Input::I_Initial::I_SurfaceShape::I_EllipseBump i_obj = inp.EllipseBump[i];
		C_EllipseBump c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Height, i_obj.Center[0], i_obj.Center[1],
			i_obj.Width[0], i_obj.Width[1], i_obj.SlopeWidth);
		_vecEllipseBump.push_back(c_obj);
		tShapeSet p(EllipseBump, i);
		_shapeOrder.push_back(p);
	}

	// FlatPlane
	for(unsigned int i = 0; i< inp.FlatPlane.size(); i++)
	{
		Input::I_Initial::I_SurfaceShape::I_FlatPlane i_obj = inp.FlatPlane[i];
		C_FlatPlane c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), 0, i_obj.Height, i_obj.Center[0], i_obj.Center[1],
			i_obj.Angle[0], i_obj.Angle[1]);
		_vecFlatPlane.push_back(c_obj);
		tShapeSet p(FlatPlane, i);
		_shapeOrder.push_back(p);
	}

	// MexicanHat
	for(unsigned int i = 0; i< inp.MexicanHat.size(); i++)
	{
		Input::I_Initial::I_SurfaceShape::I_MexicanHat i_obj = inp.MexicanHat[i];
		C_MexicanHat c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Height, i_obj.Center[0], i_obj.Center[1],
			i_obj.Width[0], i_obj.Width[1], i_obj.Beta);
		_vecMexicanHat.push_back(c_obj);
		tShapeSet p(MexicanHat, i);
		_shapeOrder.push_back(p);
	}

	// SmoothConvex
	for(unsigned int i = 0; i< inp.SmoothConvex.size(); i++)
	{
		Input::I_Initial::I_SurfaceShape::I_SmoothConvex i_obj = inp.SmoothConvex[i];
		C_SmoothConvex c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Height, i_obj.Center[0], i_obj.Center[1],
			i_obj.Width[0], i_obj.Width[1]);
		_vecSmoothConvex.push_back(c_obj);
		tShapeSet p(SmoothConvex, i);
		_shapeOrder.push_back(p);
	}

	// PeriodicWave
	for(unsigned int i = 0; i< inp.PeriodicWave.size(); i++)
	{
		Input::I_Initial::I_SurfaceShape::I_PeriodicWave i_obj = inp.PeriodicWave[i];
		C_PeriodicWave c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Amplitude, i_obj.Phase[0], i_obj.Phase[1],
			i_obj.Wavelength[0], i_obj.Wavelength[1]);
		_vecPeriodicWave.push_back(c_obj);
		tShapeSet p(PeriodicWave, i);
		_shapeOrder.push_back(p);
	}

	// EdgeBump
	double rs = 0;
	C_EdgeBump::tLocation loc = C_EdgeBump::Xmax;
	for(unsigned int i = 0; i< inp.EdgeBump.size(); i++)
	{
		Input::I_Initial::I_SurfaceShape::I_EdgeBump i_obj = inp.EdgeBump[i];
		if(i_obj.Location.compare("Xmin") == 0){loc = C_EdgeBump::Xmin; rs = _Xmin + i_obj.Lout;}
		if(i_obj.Location.compare("Xmax") == 0){loc = C_EdgeBump::Xmax; rs = _Xmax - i_obj.Lout;}
		if(i_obj.Location.compare("Ymin") == 0){loc = C_EdgeBump::Ymin; rs = _Ymin + i_obj.Lout;}
		if(i_obj.Location.compare("Ymax") == 0){loc = C_EdgeBump::Ymax; rs = _Ymax - i_obj.Lout;}
		C_EdgeBump c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Height,
			i_obj.alpha, i_obj.Lin, i_obj.Lout, rs, loc);
		_vecEdgeBump.push_back(c_obj);
		tShapeSet p(EdgeBump, i);
		_shapeOrder.push_back(p);
	}

	SortOrder();
}

void C_GenerateShape::Build(const C_Grid1D &g, double scale, const Input::I_Substrate &inp)
{
	_Xmin = g.XminReal() * scale;
	_Xmax = g.XmaxReal() * scale;
	_Ymin = 0;
	_Ymax = 0;
	Build(scale, inp);
}

void C_GenerateShape::Build(const C_Grid2D &g, double scale, const Input::I_Substrate &inp)
{
	_Xmin = g.XminReal() * scale;
	_Xmax = g.XmaxReal() * scale;
	_Ymin = g.YminReal() * scale;
	_Ymax = g.YmaxReal() * scale;

	Build(scale, inp);
}

void C_GenerateShape::Build(double scale, const Input::I_Substrate &inp)
{
	// RectBump
	for(unsigned int i = 0; i< inp.RectBump.size(); i++)
	{
		Input::I_Substrate::I_RectBump i_obj = inp.RectBump[i];
		C_RectBump c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Height, i_obj.Center[0], i_obj.Center[1],
			i_obj.Width[0], i_obj.Width[1], i_obj.SlopeWidth[0], i_obj.SlopeWidth[1]);
		_vecRectBump.push_back(c_obj);
		tShapeSet p(RectBump, i);
		_shapeOrder.push_back(p);
	}

	// EllipseBump
	for(unsigned int i = 0; i< inp.EllipseBump.size(); i++)
	{
		Input::I_Substrate::I_EllipseBump i_obj = inp.EllipseBump[i];
		C_EllipseBump c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Height, i_obj.Center[0], i_obj.Center[1],
			i_obj.Width[0], i_obj.Width[1], i_obj.SlopeWidth);
		_vecEllipseBump.push_back(c_obj);
		tShapeSet p(EllipseBump, i);
		_shapeOrder.push_back(p);
	}

	// FlatPlane
	for(unsigned int i = 0; i< inp.FlatPlane.size(); i++)
	{
		Input::I_Substrate::I_FlatPlane i_obj = inp.FlatPlane[i];
		C_FlatPlane c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), 0, i_obj.Height, i_obj.Center[0], i_obj.Center[1],
			i_obj.Angle[0], i_obj.Angle[1]);
		_vecFlatPlane.push_back(c_obj);
		tShapeSet p(FlatPlane, i);
		_shapeOrder.push_back(p);
	}

	// MexicanHat
	for(unsigned int i = 0; i< inp.MexicanHat.size(); i++)
	{
		Input::I_Substrate::I_MexicanHat i_obj = inp.MexicanHat[i];
		C_MexicanHat c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Height, i_obj.Center[0], i_obj.Center[1],
			i_obj.Width[0], i_obj.Width[1], i_obj.Beta);
		_vecMexicanHat.push_back(c_obj);
		tShapeSet p(MexicanHat, i);
		_shapeOrder.push_back(p);
	}

	// SmoothConvex
	for(unsigned int i = 0; i< inp.SmoothConvex.size(); i++)
	{
		Input::I_Substrate::I_SmoothConvex i_obj = inp.SmoothConvex[i];
		C_SmoothConvex c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Height, i_obj.Center[0], i_obj.Center[1],
			i_obj.Width[0], i_obj.Width[1]);
		_vecSmoothConvex.push_back(c_obj);
		tShapeSet p(SmoothConvex, i);
		_shapeOrder.push_back(p);
	}

	// PeriodicWave
	for(unsigned int i = 0; i< inp.PeriodicWave.size(); i++)
	{
		Input::I_Substrate::I_PeriodicWave i_obj = inp.PeriodicWave[i];
		C_PeriodicWave c_obj;
		c_obj.Build(i_obj.Order, GetType(i_obj.CreateType), i_obj.Base, i_obj.Amplitude, i_obj.Phase[0], i_obj.Phase[1],
			i_obj.Wavelength[0], i_obj.Wavelength[1]);
		_vecPeriodicWave.push_back(c_obj);
		tShapeSet p(PeriodicWave, i);
		_shapeOrder.push_back(p);
	}

	SortOrder();

}

C_ShapeFunc_base::tCreateType C_GenerateShape::GetType(std::string type)
{
	C_ShapeFunc_base::tCreateType tt = C_ShapeFunc_base::Add;
	if(type.compare("add")==0) tt =  C_ShapeFunc_base::Add;
	else if(type.compare("remove")==0) tt =  C_ShapeFunc_base::Remove;
	else if(type.compare("stack")==0) tt =  C_ShapeFunc_base::Stack;

	return tt;
}

int C_GenerateShape::GetOrder(tShapeSet set)
{
	int order = 0;
	switch(set.Type())
	{
	default: order = 0; break;	// unknown shape kind; 0 keeps this from indexing _vecRectBump[] out of range
	case RectBump:		order = _vecRectBump[set.Index()].Order(); break;
	case EllipseBump:	order = _vecEllipseBump[set.Index()].Order(); break;
	case FlatPlane:		order = _vecFlatPlane[set.Index()].Order(); break;
	case MexicanHat:	order = _vecMexicanHat[set.Index()].Order(); break;
	case SmoothConvex:	order = _vecSmoothConvex[set.Index()].Order(); break;
	case PeriodicWave:	order = _vecPeriodicWave[set.Index()].Order(); break;
	case EdgeBump:		order = _vecEdgeBump[set.Index()].Order(); break;
	}
	return order;
}

void C_GenerateShape::SortOrder()
{
	tShapeSet set_tmp;
	int order_i, order_j;
	for(unsigned int i = 0; i< _shapeOrder.size(); i++)
	{
		for(unsigned int j = i+1; j< _shapeOrder.size(); j++)
		{
			order_i = GetOrder(_shapeOrder[i]);
			order_j = GetOrder(_shapeOrder[j]);

			if(order_i > order_j)
			{
				set_tmp			= _shapeOrder[i];
				_shapeOrder[i]	= _shapeOrder[j];
				_shapeOrder[j]	= set_tmp;
			}
		}
	}
}

double C_GenerateShape::Evaluate(double px)
{
	return Evaluate(px, 0);
}

double C_GenerateShape::Evaluate(double px, double py)
{
	double mixed_val = 0;
	double val = 0;
	C_ShapeFunc_base::tCreateType type = C_ShapeFunc_base::Add;

	tShapeSet set;
	for(unsigned int i = 0; i< _shapeOrder.size(); i++)
	{
		set = _shapeOrder[i];
		switch(set.Type())
		{
		default: val = 0; break;	// unknown shape kind (as above)
		case RectBump:
			type = _vecRectBump[set.Index()].Type();
			val  = _vecRectBump[set.Index()](px,py);
			break;
		case EllipseBump:
			type = _vecEllipseBump[set.Index()].Type();
			val  = _vecEllipseBump[set.Index()](px,py);
			break;
		case FlatPlane:
			type = _vecFlatPlane[set.Index()].Type();
			val  = _vecFlatPlane[set.Index()](px,py);
			break;
		case MexicanHat:
			type = _vecMexicanHat[set.Index()].Type();
			val  = _vecMexicanHat[set.Index()](px,py);
			break;
		case SmoothConvex:
			type = _vecSmoothConvex[set.Index()].Type();
			val  = _vecSmoothConvex[set.Index()](px,py);
			break;
		case PeriodicWave:
			type = _vecPeriodicWave[set.Index()].Type();
			val  = _vecPeriodicWave[set.Index()](px,py);
			break;
		case EdgeBump:
			type = _vecEdgeBump[set.Index()].Type();
			val  = _vecEdgeBump[set.Index()](px,py);
			break;
		}

		switch(type)
		{
		default:
		case C_ShapeFunc_base::Add:
			mixed_val = Math::Max(mixed_val, val);
			break;
		case C_ShapeFunc_base::Remove:
			mixed_val = Math::Min(mixed_val, val);
			break;
		case C_ShapeFunc_base::Stack:
			mixed_val += val;
			break;
		}
	}

	return mixed_val;
}
