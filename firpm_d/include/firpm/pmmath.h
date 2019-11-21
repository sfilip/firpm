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

namespace pmmath {
	template<typename T> constexpr T sin(T);
	template<typename T> constexpr T cos(T);
	template<typename T> constexpr T tan(T);

	template<typename T> constexpr T asin(T);
	template<typename T> constexpr T acos(T);
	template<typename T> constexpr T atan(T);

	template<typename T> constexpr T log(T);
	template<typename T> constexpr T exp(T);
	template<typename T> constexpr T sqrt(T);
	template<typename T> constexpr T pow(T, T);

	template<typename T> constexpr T fabs(T);
	template<typename T> constexpr T fmax(T, T);
	template<typename T> constexpr T fmin(T, T);

	/* Specializations: double precision */
	template<> constexpr double sin<double>(double x) { return std::sin(x); };
	template<> constexpr double cos<double>(double x) { return std::cos(x); };
	template<> constexpr double tan<double>(double x) { return std::tan(x); };

	template<> constexpr double asin<double>(double x) { return std::asin(x); };
	template<> constexpr double acos<double>(double x) { return std::acos(x); };
	template<> constexpr double atan<double>(double x) { return std::atan(x); };

	template<> constexpr double log<double>(double x) { return std::log(x); };
	template<> constexpr double exp<double>(double x) { return std::exp(x); };
	template<> constexpr double sqrt<double>(double x) { return std::sqrt(x); };
	template<> constexpr double pow<double>(double x, double y) { return std::pow(x, y); };

	template<> constexpr double fabs<double>(double x) { return std::fabs(x); };
	template<> constexpr double fmax<double>(double x, double y) { return std::fmax(x, y); };
	template<> constexpr double fmin<double>(double x, double y) { return std::fmin(x, y); };

	/* Specializations: long double precision */
	template<> constexpr long double sin<long double>(long double x) { return sinl(x); };
	template<> constexpr long double cos<long double>(long double x) { return cosl(x); };
	template<> constexpr long double tan<long double>(long double x) { return tanl(x); };

	template<> constexpr long double asin<long double>(long double x) { return asinl(x); };
	template<> constexpr long double acos<long double>(long double x) { return acosl(x); };
	template<> constexpr long double atan<long double>(long double x) { return atanl(x); };

	template<> constexpr long double log<long double>(long double x) { return logl(x); };
	template<> constexpr long double exp<long double>(long double x) { return expl(x); };
	template<> constexpr long double sqrt<long double>(long double x) { return sqrtl(x); };
	template<> constexpr long double pow<long double>(long double x, long double y) { return powl(x, y); };

	template<> constexpr long double fabs<long double>(long double x) { return fabsl(x); };
	template<> constexpr long double fmax<long double>(long double x, long double y) { return fmaxl(x, y); };
	template<> constexpr long double fmin<long double>(long double x, long double y) { return fminl(x, y); };
}

#endif
