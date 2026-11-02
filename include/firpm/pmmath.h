#ifndef __PMMATH_H__
#define __PMMATH_H__

/**
 * @file pmmath.h
 * @author Graeme Smecher
 * @date 20 Nov 2019
 * @brief Math templates that dispatch to appropriate type specializations
 */

// firpm
// Copyright (C) 2019  G. Smecher

#include <cmath>
#include <numbers>
#ifdef HAVE_MPFR
	#include "mpreal.h"
#endif

namespace pm {
	namespace pmmath {
		/* Unqualified calls resolve by argument-dependent lookup: std:: for
		 * double and long double, mpfr:: for mpreal. */
		template<typename T> T sin(T x) { using std::sin; return sin(x); };
		template<typename T> T cos(T x) { using std::cos; return cos(x); };
		template<typename T> T tan(T x) { using std::tan; return tan(x); };
		template<typename T> T asin(T x) { using std::asin; return asin(x); };
		template<typename T> T acos(T x) { using std::acos; return acos(x); };
		template<typename T> T atan(T x) { using std::atan; return atan(x); };
		template<typename T> T log(T x) { using std::log; return log(x); };
		template<typename T> T exp(T x) { using std::exp; return exp(x); };
		template<typename T> T sqrt(T x) { using std::sqrt; return sqrt(x); };
		template<typename T> T fabs(T x) { using std::fabs; return fabs(x); };
		template<typename T> T pow(T x, T y) { using std::pow; return pow(x, y); };
		template<typename T> T fmax(T x, T y) { using std::fmax; return fmax(x, y); };
		template<typename T> T fmin(T x, T y) { using std::fmin; return fmin(x, y); };
		template<typename T> bool signbit(T x) { using std::signbit; return signbit(x); };
		template<typename T> bool isfinite(T x) { using std::isfinite; return isfinite(x); };
		template<typename T> bool isnan(T x) { using std::isnan; return isnan(x); };
		template<typename T> long round(T x) { using std::round; return static_cast<long>(round(x)); };
		template<typename T> long floor(T x) { using std::floor; return static_cast<long>(floor(x)); };
		/* x = m * 2^e with |m| in [0.5, 1): mpreal supplies frexp(x, int*)
		 * and ldexp(x, mp_exp_t), both found by ADL. */
		template<typename T> T frexp(T x, int* e) { using std::frexp; return frexp(x, e); };
		template<typename T> T ldexp(T x, int e) { using std::ldexp; return ldexp(x, e); };

		template<typename T> T const_pi(void) { return std::numbers::pi_v<T>; };
#ifdef HAVE_MPFR
		/* Computed at the current default precision. */
		template<> inline mpfr::mpreal const_pi<mpfr::mpreal>(void) { return mpfr::const_pi(); };
#endif
	} // namespace pmmath
} // namespace pm

#endif
