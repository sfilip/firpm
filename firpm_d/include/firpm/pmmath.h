#ifndef __PMMATH_H__
#define __PMMATH_H__

/**
 * @file pmmath.h
 * @author Graeme Smecher
 * @date 20 Nov 2019
 * @brief Math templates that dispatch to appropriate type specializations
 */

// firpm_d
// Copyright (C) 2019  G. Smecher
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>

#include <cmath>
#ifdef HAVE_MPFR
	#include "mpreal.h"
#endif

namespace pmmath {
	template<typename T> T sin(T);
	template<typename T> T cos(T);
	template<typename T> T tan(T);

	template<typename T> T asin(T);
	template<typename T> T acos(T);
	template<typename T> T atan(T);

	template<typename T> T log(T);
	template<typename T> T exp(T);
	template<typename T> T sqrt(T);
	template<typename T> T pow(T, T);

	template<typename T> T fabs(T);
	template<typename T> T fmax(T, T);
	template<typename T> T fmin(T, T);

	template<typename T> bool signbit(T);
	template<typename T> bool isfinite(T);
	template<typename T> bool isnan(T);

	template<typename T> long round(T);

	/* Specializations: double precision */
	template<> inline double sin<double>(double x) { return std::sin(x); };
	template<> inline double cos<double>(double x) { return std::cos(x); };
	template<> inline double tan<double>(double x) { return std::tan(x); };

	template<> inline double asin<double>(double x) { return std::asin(x); };
	template<> inline double acos<double>(double x) { return std::acos(x); };
	template<> inline double atan<double>(double x) { return std::atan(x); };

	template<> inline double log<double>(double x) { return std::log(x); };
	template<> inline double exp<double>(double x) { return std::exp(x); };
	template<> inline double sqrt<double>(double x) { return std::sqrt(x); };
	template<> inline double pow<double>(double x, double y) { return std::pow(x, y); };

	template<> inline double fabs<double>(double x) { return std::fabs(x); };
	template<> inline double fmax<double>(double x, double y) { return std::fmax(x, y); };
	template<> inline double fmin<double>(double x, double y) { return std::fmin(x, y); };

	template<> inline bool signbit<double>(double x) { return std::signbit(x); };
	template<> inline bool isfinite<double>(double x) { return std::isfinite(x); };
	template<> inline bool isnan<double>(double x) { return std::isnan(x); };

	template<> inline long round<double>(double x) { return std::round(x); };

	/* Specializations: long double precision */
	template<> inline long double sin<long double>(long double x) { return sinl(x); };
	template<> inline long double cos<long double>(long double x) { return cosl(x); };
	template<> inline long double tan<long double>(long double x) { return tanl(x); };

	template<> inline long double asin<long double>(long double x) { return asinl(x); };
	template<> inline long double acos<long double>(long double x) { return acosl(x); };
	template<> inline long double atan<long double>(long double x) { return atanl(x); };

	template<> inline long double log<long double>(long double x) { return logl(x); };
	template<> inline long double exp<long double>(long double x) { return expl(x); };
	template<> inline long double sqrt<long double>(long double x) { return sqrtl(x); };
	template<> inline long double pow<long double>(long double x, long double y) { return powl(x, y); };

	template<> inline long double fabs<long double>(long double x) { return fabsl(x); };
	template<> inline long double fmax<long double>(long double x, long double y) { return fmaxl(x, y); };
	template<> inline long double fmin<long double>(long double x, long double y) { return fminl(x, y); };

	template<> inline bool signbit<long double>(long double x) { return std::signbit(x); };
	template<> inline bool isfinite<long double>(long double x) { return std::isfinite(x); };
	template<> inline bool isnan<long double>(long double x) { return std::isnan(x); };

	template<> inline long round<long double>(long double x) { return std::round(x); };

	/* Specialization: multiple precision mpreal */
#ifdef HAVE_MPFR
	template<> inline mpfr::mpreal sin<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::sin(x); };
	template<> inline mpfr::mpreal cos<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::cos(x); };
	template<> inline mpfr::mpreal tan<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::tan(x); };

	template<> inline mpfr::mpreal asin<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::asin(x); };
	template<> inline mpfr::mpreal acos<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::acos(x); };
	template<> inline mpfr::mpreal atan<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::atan(x); };

	template<> inline mpfr::mpreal log<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::log(x); };
	template<> inline mpfr::mpreal exp<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::exp(x); };
	template<> inline mpfr::mpreal sqrt<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::sqrt(x); };
	template<> inline mpfr::mpreal pow<mpfr::mpreal>(mpfr::mpreal x, mpfr::mpreal y) { return mpfr::pow(x, y); };

	template<> inline mpfr::mpreal fabs<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::abs(x); };
	template<> inline mpfr::mpreal fmax<mpfr::mpreal>(mpfr::mpreal x, mpfr::mpreal y) { return mpfr::max(x, y); };
	template<> inline mpfr::mpreal fmin<mpfr::mpreal>(mpfr::mpreal x, mpfr::mpreal y) { return mpfr::min(x, y); };

	template<> inline bool signbit<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::signbit(x); };
	template<> inline bool isfinite<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::isfinite(x); };
	template<> inline bool isnan<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::isnan(x); };

	template<> inline long round<mpfr::mpreal>(mpfr::mpreal x) { return mpfr::round(x).toLong(); };
#endif

}

#endif
