/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  MyMath_I.h : bodies of the inline functions (included from MyMath.h)
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

inline int    Math::abs(const int    &val)
{
	return ((val < 0) ? -val : val);
}

inline double Math::abs(const double &val)
{
	return ((val < 0) ? -val : val);
}

inline bool Math::IsZero(const double &val)
{
#ifdef _WIN32
	if(abs(val) <= DBL_MIN )return true;
#else
	if(abs(val) <= std::numeric_limits<double>::min() )return true;
#endif
	else return false;
}

inline bool Math::IsLessEpsilon(const double &val)
{
	return ( abs(val) <= std::numeric_limits<double>::epsilon() ? true : false);
}

inline bool Math::Equal(const double &v1, const double &v2)
{
	// Knuth's method
	double criterion = std::numeric_limits<double>::epsilon();
	if(abs(v1-v2) <= criterion *abs(v1) || abs(v1-v2) <= criterion *abs(v2)) return true;
	else return false;
}


inline int    Math::Min(const int    &v1, const int    &v2)
{
	return (v1 < v2 ? v1 : v2);
}

inline double Math::Min(const double &v1, const double &v2)
{
	return (v1 < v2 ? v1 : v2);
}

inline int    Math::Max(const int    &v1, const int    &v2)
{
	return (v1 > v2 ? v1 : v2);
}

inline double Math::Max(const double &v1, const double &v2)
{
	return (v1 > v2 ? v1 : v2);
}

inline int    Math::Min(const int    &v1, const int    &v2, const int    &v3)
{
	int v23 = Min(v2,v3);
	return (v1 < v23 ? v1 : v23);
}

inline double Math::Min(const double &v1, const double &v2, const double &v3)
{
	double v23 = Min(v2,v3);
	return (v1 < v23 ? v1 : v23);
}

inline int    Math::Max(const int    &v1, const int    &v2, const int    &v3)
{
	int v23 = Max(v2,v3);
	return (v1 > v23 ? v1 : v23);
}

inline double Math::Max(const double &v1, const double &v2, const double &v3)
{
	double v23 = Max(v2,v3);
	return (v1 > v23 ? v1 : v23);
}


inline int    Math::Sign(const int    &v1)
{
	return (v1 >= 0 ? 1 : -1);
}

inline double Math::Sign(const double &v1)
{
	return (v1 >= 0.0 ? 1.0 : -1.0);
}

inline int    Math::Sign(const int    &value,const int    &for_sign)
{
	return (for_sign >= 0.0 ? abs(value) : -abs(value));
}

inline double Math::Sign(const double &value, const double &for_sign)
{
	return (for_sign >= 0.0 ? abs(value) : -abs(value));
}

inline double Math::Sigmoid(const double &xx, const double &x0, const double &eps)
{
	return 1.0/(1.0 + exp(-(xx-x0)/(eps)));
}

inline double Math::Diff(double a, double b)
{
	double denom = a + b;
	if(IsZero(denom)) return denom;
	double cc = (a*a) - (b*b);
	return cc / denom;
}
