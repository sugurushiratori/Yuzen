/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  FiniteDifferences.cpp : implementation of the difference and interpolation scheme coefficients
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#include <fstream>
#include <cmath>

#include "FiniteDifferences.h"



#ifdef _DEF_YUZEN_
using namespace yuzen;
#endif

///////////////////////////////////////////////////////////////////////
void C_Interpolation_base4::SetCoef()
{
	_c1.resize(_n_stencil);
	_c2.resize(_n_stencil);
	double denom = 1;
	switch(_t_scheme)
	{
	case Upwind_1st:	// 1st-upwind scheme.
		_c1[L] =  1; _c1[R] = 1;
		_c2[L] =  1; _c2[R] = -1;
		denom = _c1.Sum();
		break;
	case Central_2nd:	// 2nd-central scheme.
		_c1[L] = 1; _c1[R] = 1;
		denom = _c1.Sum();
		break;
	case Central_4th:	// 4th-central scheme.
		_c1[LL] = -1; _c1[L] = 9; _c1[R] =  9; _c1[RR] = -1;
		denom = _c1.Sum();
		break;
	case Upwind_2nd:	// 2nd-upwind scheme.
		_c1[LL] = -1; _c1[L] = 3; _c1[R] =  3; _c1[RR] = -1;
		_c2[LL] = -1; _c2[L] = 3; _c2[R] = -3; _c2[RR] =  1;
		denom = _c1.Sum();
		break;
	case Upwind_3rd:	// 3rd-upwind scheme. (UTOPIA)
		_c1[LL] = -1; _c1[L] = 7; _c1[R] =  7; _c1[RR] = -1;
		_c2[LL] = -1; _c2[L] = 3; _c2[R] = -3; _c2[RR] =  1;
		denom = _c1.Sum();
		break;
	case QUICK:		// QUICK scheme by Leonard.(Quadratic Upstream Interpolation for Convective Kinematics)
		_c1[LL] = -1; _c1[L] = 9; _c1[R] =  9; _c1[RR] = -1;
		_c2[LL] = -1; _c2[L] = 3; _c2[R] = -3; _c2[RR] =  1;
		denom = _c1.Sum();
		break;
	case KK:		// Kuwahara-Kawamura scheme.
		_c1[LL] = -1; _c1[L] = 7; _c1[R] =  7; _c1[RR] = -1;
		_c2[LL] = -1; _c2[L] = 3; _c2[R] = -3; _c2[RR] =  1;
		denom = _c1.Sum();
		_c2 *= 3.0;
		break;
	case Koren:		// limited; where smooth it is Upwind_3rd
		_c1[LL] = -1; _c1[L] = 7; _c1[R] =  7; _c1[RR] = -1;
		_c2[LL] = -1; _c2[L] = 3; _c2[R] = -3; _c2[RR] =  1;
		denom = _c1.Sum();
		break;
	case VanLeer:	// limited; linearized about r = 1 it is kappa = 0 (Fromm)
		_c1[LL] = -1; _c1[L] = 5; _c1[R] =  5; _c1[RR] = -1;
		_c2[LL] = -1; _c2[L] = 3; _c2[R] = -3; _c2[RR] =  1;
		denom = _c1.Sum();
		break;
	}
	_c1 /= denom;
	_c2 /= denom;
}

///////////////////////////////////////////////////////////////////////
std::string C_Interpolation_base4::SchemeName()
{
	std::string name("");
	switch(_t_scheme)
	{
	case Upwind_1st:	name = "Upwind-1st";  break;
	case Central_2nd:	name = "Central-2nd";  break;
	case Central_4th:	name = "Central-4th";  break;
	case Upwind_2nd:	name = "Upwind-2nd";  break;
	case Upwind_3rd:	name = "Upwind-3rd";  break;
	case QUICK:			name = "QUICK";  break;
	case KK:			name = "Kawamura-Kuwahara";  break;
	case Koren:			name = "Koren (TVD)";  break;
	case VanLeer:		name = "van Leer (TVD)";  break;
	}
	return name;
}

///////////////////////////////////////////////////////////////////////
//  The limited schemes
///////////////////////////////////////////////////////////////////////
double C_Interpolation_base4::Limiter(double a, double b, double *dA, double *dB) const
{
	double phi = 0, pa = 0, pb = 0;
	//  Zero unless a and b have the same (strict) sign. Written with
	//  comparisons rather than a*b > 0, which underflows.
	const bool same = (a > 0 && b > 0) || (a < 0 && b < 0);
	if(same)
	{
		switch(_t_scheme)
		{
		case Koren:
		{
			const double s  = (b > 0) ? 1.0 : -1.0;
			const double aa = std::fabs(a), ab = std::fabs(b);
			const double t1 = 2.0*aa, t2 = (2.0*ab + aa)/3.0, t3 = 2.0*ab;
			//  At a tie the branch listed first is taken, the same for (a, b)
			//  and (-a, -b)
			if(t1 <= t2 && t1 <= t3)	{ phi = s*t1; pa = 2.0;     pb = 0.0; }
			else if(t2 <= t3)			{ phi = s*t2; pa = 1.0/3.0; pb = 2.0/3.0; }
			else						{ phi = s*t3; pa = 0.0;     pb = 2.0; }
			break;
		}
		case VanLeer:
		{
			//  In the ratios a/(a+b) and b/(a+b), both in (0, 1) since a and
			//  b have the same sign: 2ab/(a+b) and b^2/(a+b)^2 as they stood
			//  underflowed to 0 and 0/0 = NaN for |a|, |b| below ~1e-162 and
			//  overflowed above ~1e154.
			const double sum = a + b;			// never 0: a and b have the same sign
			const double ra = a/sum, rb = b/sum;
			phi = 2.0*a*rb;
			pa  = 2.0*rb*rb;
			pb  = 2.0*ra*ra;
			break;
		}
		default:
			break;
		}
	}
	if(dA) *dA = pa;
	if(dB) *dB = pb;
	return phi;
}

namespace {
//  The stencil positions of UU, U and D for the direction of the flow
inline void UpwindPositions(double sign, int &uu, int &u, int &d)
{
	if(sign >= 0) { uu = C_Interpolation_base4::LL; u = C_Interpolation_base4::L; d = C_Interpolation_base4::R; }
	else          { uu = C_Interpolation_base4::RR; u = C_Interpolation_base4::R; d = C_Interpolation_base4::L; }
}
}

double C_Interpolation_base4::FaceValue(const double f[4], double sign) const
{
	if(!IsLimited())
	{
		double v = 0;
		for(int j = 0; j < 4; ++j) v += (_c1[j] + Math::Sign(sign)*_c2[j]) * f[j];
		return v;
	}
	int uu, u, d;
	UpwindPositions(sign, uu, u, d);
	return f[u] + 0.5*Limiter(f[u] - f[uu], f[d] - f[u]);
}

void C_Interpolation_base4::FaceWeights(const double f[4], double sign, int rowPos, double w[4]) const
{
	if(!IsLimited())
	{
		for(int j = 0; j < 4; ++j) w[j] = _c1[j] + Math::Sign(sign)*_c2[j];
		return;
	}
	int uu, u, d;
	UpwindPositions(sign, uu, u, d);
	const double a = f[u] - f[uu], b = f[d] - f[u];
	const double phi = Limiter(a, b);
	for(int j = 0; j < 4; ++j) w[j] = 0;
	if(rowPos == u)
	{	//  Upwind form: f_U + chi/2 (f_U - f_UU)
		const double chi = (a != 0) ? phi / a : 0.0;
		w[u]  = 1.0 + 0.5*chi;
		w[uu] = -0.5*chi;
	}
	else
	{	//  Downwind form: f_U + psi/2 (f_D - f_U)
		const double psi = (b != 0) ? phi / b : 0.0;
		w[u] = 1.0 - 0.5*psi;
		w[d] = 0.5*psi;
	}
}

double C_Interpolation_base4::FaceDerivative(const double f[4], double sign, double g[4]) const
{
	if(!IsLimited())
	{
		//  The value is FaceValue's own. Summed here from g[] it would be the
		//  same arithmetic, but a compiler that contracts across statements
		//  (clang -ffp-contract=fast) can fuse the two loops differently, and
		//  the two would disagree in the last bit.
		for(int j = 0; j < 4; ++j) g[j] = _c1[j] + Math::Sign(sign)*_c2[j];
		return FaceValue(f, sign);
	}
	int uu, u, d;
	UpwindPositions(sign, uu, u, d);
	double pa, pb;
	const double phi = Limiter(f[u] - f[uu], f[d] - f[u], &pa, &pb);
	for(int j = 0; j < 4; ++j) g[j] = 0;
	//  a = f_U - f_UU, b = f_D - f_U
	g[uu] = -0.5*pa;
	g[u]  = 1.0 + 0.5*(pa - pb);
	g[d]  = 0.5*pb;
	return f[u] + 0.5*phi;
}

namespace {
struct KeyEntry { const char *key; C_Interpolation_base4::t_Scheme t; };
const KeyEntry kKeys[] = {
	{"KK",       C_Interpolation_base4::KK},
	{"UTOPIA",   C_Interpolation_base4::Upwind_3rd},
	{"Central4", C_Interpolation_base4::Central_4th},
	{"Koren",    C_Interpolation_base4::Koren},
	{"VanLeer",  C_Interpolation_base4::VanLeer},
	{"QUICK",    C_Interpolation_base4::QUICK},
	{"Upwind1",  C_Interpolation_base4::Upwind_1st},
	{"Upwind2",  C_Interpolation_base4::Upwind_2nd},
	{"Central2", C_Interpolation_base4::Central_2nd}};
const int kNumKeys = sizeof(kKeys) / sizeof(kKeys[0]);
}

bool C_Interpolation_base4::FromKey(const std::string &key, t_Scheme &t)
{
	for(int i = 0; i < kNumKeys; ++i)
		if(key == kKeys[i].key) { t = kKeys[i].t; return true; }
	return false;
}

std::string C_Interpolation_base4::Key(t_Scheme t)
{
	for(int i = 0; i < kNumKeys; ++i)
		if(t == kKeys[i].t) return kKeys[i].key;
	return "";
}

std::string C_Interpolation_base4::KeyList()
{
	std::string s;
	for(int i = 0; i < kNumKeys; ++i)
	{
		if(i) s += ", ";
		s += std::string("`") + kKeys[i].key + "'";
	}
	return s;
}

/////////////////////////////////////////////////////////////////////////
double C_FiniteDifferences::GradX(const ScalarField1D<> &s, const Node  &nd, double sign)
{
	double gR, gL;
	Xedge xeR, xeL;
	const RectGrid<1> &G = s.Grid();
	xeR = G.HalfXp(nd);
	xeL = G.HalfXm(nd);
	gR = GradX(s, xeR, sign);
	gL = GradX(s, xeL, sign);
	return 0.5*(gR+gL);
}
double C_FiniteDifferences::GradX(const ScalarField1D<> &s, const Xedge &xe, double sign)
{
	const RectGrid<1> &G = s.Grid();
	Vector<double> ff(5);
	double dx = G.dx();
	ff[0]  = s(G.ShiftXm(G.ShiftXm(xe)));
	ff[1]  = s(G.ShiftXm(xe));
	ff[2]  = s(xe);
	ff[3]  = s(G.ShiftXp(xe));
	ff[4]  = s(G.ShiftXp(G.ShiftXp(xe)));

	double fR = 0, fL = 0;
	for(int i=0;i<4;i++)
	{
		// fR takes ff[i+1] and fL ff[i], as in the other eleven
		//          GradX/GradY/GradZ overloads; the other way round inverts
		//          the sign.
		fR += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i+1];
		fL += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i];
	}
	return (fR-fL)/dx;
}

double C_FiniteDifferences::GradX(const ScalarField2D<> &s, const Zface &zf, double sign)
{
	const RectGrid<2> &G = s.Grid();
	Vector<double> ff(5);
	double dx = G.dx();
	ff[0]  = s(G.ShiftXm(G.ShiftXm(zf)));
	ff[1]  = s(G.ShiftXm(zf));
	ff[2]  = s(zf);
	ff[3]  = s(G.ShiftXp(zf));
	ff[4]  = s(G.ShiftXp(G.ShiftXp(zf)));

	double fR = 0, fL = 0;
	for(int i=0;i<4;i++)
	{
		fR += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i+1];
		fL += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i];
	}
	return (fR-fL)/dx;
}
double C_FiniteDifferences::GradX(const ScalarField2D<> &s, const Yedge &ye, double sign)
{
	double gR, gL;
	Zface zfR, zfL;
	const RectGrid<2> &G = s.Grid();
	zfR = G.HalfXp(ye);
	zfL = G.HalfXm(ye);
	gR = GradX(s, zfR, sign);
	gL = GradX(s, zfL, sign);
	return 0.5*(gR+gL);
}
double C_FiniteDifferences::GradY(const ScalarField2D<> &s, const Zface &zf, double sign)
{
	const RectGrid<2> &G = s.Grid();
	Vector<double> ff(5);
	double dy = G.dy();
	ff[0]  = s(G.ShiftYm(G.ShiftYm(zf)));
	ff[1]  = s(G.ShiftYm(zf));
	ff[2]  = s(zf);
	ff[3]  = s(G.ShiftYp(zf));
	ff[4]  = s(G.ShiftYp(G.ShiftYp(zf)));

	double fR = 0, fL = 0;
	for(int i=0;i<4;i++)
	{
		fR += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i+1];
		fL += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i];
	}
	return (fR-fL)/dy;
}
double C_FiniteDifferences::GradY(const ScalarField2D<> &s, const Xedge &xe, double sign)
{
	double gR, gL;
	Zface zfR, zfL;
	const RectGrid<2> &G = s.Grid();
	zfR = G.HalfYp(xe);
	zfL = G.HalfYm(xe);
	gR = GradY(s, zfR, sign);
	gL = GradY(s, zfL, sign);
	return 0.5*(gR+gL);
}

double C_FiniteDifferences::GradX(const ScalarField3D<> &s, const Cell  &cl, double sign)
{
	const RectGrid<3> &G = s.Grid();
	Vector<double> ff(5);
	double dx = G.dx();
	ff[0]  = s(G.ShiftXm(G.ShiftXm(cl)));
	ff[1]  = s(G.ShiftXm(cl));
	ff[2]  = s(cl);
	ff[3]  = s(G.ShiftXp(cl));
	ff[4]  = s(G.ShiftXp(G.ShiftXp(cl)));

	double fR = 0, fL = 0;
	for(int i=0;i<4;i++)
	{
		fR += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i+1];
		fL += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i];
	}
	return (fR-fL)/dx;
}
double C_FiniteDifferences::GradX(const ScalarField3D<> &s, const Xface &xf, double sign)
{
	double gR, gL;
	Cell clR, clL;
	const RectGrid<3> &G = s.Grid();
	clR = G.HalfXp(xf);
	clL = G.HalfXm(xf);
	gR = GradX(s, clR, sign);
	gL = GradX(s, clL, sign);
	return 0.5*(gR+gL);
}
double C_FiniteDifferences::GradY(const ScalarField3D<> &s, const Cell  &cl, double sign)
{
	const RectGrid<3> &G = s.Grid();
	Vector<double> ff(5);
	double dy = G.dx();
	ff[0]  = s(G.ShiftYm(G.ShiftYm(cl)));
	ff[1]  = s(G.ShiftYm(cl));
	ff[2]  = s(cl);
	ff[3]  = s(G.ShiftYp(cl));
	ff[4]  = s(G.ShiftYp(G.ShiftYp(cl)));

	double fR = 0, fL = 0;
	for(int i=0;i<4;i++)
	{
		fR += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i+1];
		fL += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i];
	}
	return (fR-fL)/dy;
}
double C_FiniteDifferences::GradY(const ScalarField3D<> &s, const Yface &yf, double sign)
{
	double gR, gL;
	Cell clR, clL;
	const RectGrid<3> &G = s.Grid();
	clR = G.HalfYp(yf);
	clL = G.HalfYm(yf);
	gR = GradY(s, clR, sign);
	gL = GradY(s, clL, sign);
	return 0.5*(gR+gL);
}

double C_FiniteDifferences::GradZ(const ScalarField3D<> &s, const Cell  &cl, double sign)
{
	const RectGrid<3> &G = s.Grid();
	Vector<double> ff(5);
	double dz = G.dz();
	ff[0]  = s(G.ShiftZm(G.ShiftZm(cl)));
	ff[1]  = s(G.ShiftZm(cl));
	ff[2]  = s(cl);
	ff[3]  = s(G.ShiftZp(cl));
	ff[4]  = s(G.ShiftZp(G.ShiftZp(cl)));

	double fR = 0, fL = 0;
	for(int i=0;i<4;i++)
	{
		fR += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i+1];
		fL += (_c1[i] +Math::Sign(sign)*_c2[i])*ff[i];
	}
	return (fR-fL)/dz;
}
double C_FiniteDifferences::GradZ(const ScalarField3D<> &s, const Zface &zf, double sign)
{
	double gR, gL;
	Cell clR, clL;
	const RectGrid<3> &G = s.Grid();
	clR = G.HalfZp(zf);
	clL = G.HalfZm(zf);
	gR = GradZ(s, clR, sign);
	gL = GradZ(s, clL, sign);
	return 0.5*(gR+gL);
}
