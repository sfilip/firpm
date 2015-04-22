/**
 * @file barycentric.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief Functions used to perform barycentric Lagrange interpolation
 * during the execution of the Parks-McClellan algorithm
 *
 */

//    firpm_mp
//    Copyright (C) 2015  S. Filip
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>




#ifndef BARYCENTRIC_H_
#define BARYCENTRIC_H_

#include "util.h"
#include "cheby.h"
#include "band.h"

/*! Procedure which computes the weights used in
 * the evaluation of the barycentric interpolation
 * formulas (see [Berrut&Trefethen2004] and [Pachon&Trefethen2009]
 * for the implementation ideas)
 * @param[out] w the computed weights
 * @param[in] x the interpolation points
 * @param[in] prec numerical precision used for the computations
 */
void barycentricWeights(std::vector<mpfr::mpreal>& w,
        std::vector<mpfr::mpreal>& x, mp_prec_t prec = 165ul);

/*! Determines the current reference error according to the
 * barycentric formula (internally it also computes the barycentric weights)
 * @param[out] delta the value of the current reference error
 * @param[in] x the current reference set (i.e. interpolation points)
 * @param[in] bands information relating to the ideal frequency response of
 * the filter
 * @param[in] prec numerical precision used for the computations
 */
void computeDelta(mpfr::mpreal &delta, std::vector<mpfr::mpreal>& x,
        std::vector<Band> &bands, mp_prec_t prec = 165ul);

/*! Determines the current reference error according to the
 * barycentric formula
 * @param[out] delta the value of the current reference error
 * @param[in] w the barycentric weights associated with the current reference
 * set
 * @param[in] x the current reference set (i.e. interpolation points)
 * @param[in] bands information relating to the ideal frequency response of
 * @param[in] prec numerical precision used for the computations
 */

void computeDelta(mpfr::mpreal &delta, std::vector<mpfr::mpreal>& w,
        std::vector<mpfr::mpreal>& x, std::vector<Band> &bands,
        mp_prec_t prec = 165ul);

/*! Computes the filter response at the current reference set
 * @param[out] C the vector of frequency responses at the reference set
 * @param[in] delta the current reference error
 * @param[in] x the current reference vector
 * @param[in] bands frequency band information for the ideal filter
 * @param[in] prec numerical precision used for the computations
 */
void computeC(std::vector<mpfr::mpreal> &C, mpfr::mpreal &delta,
        std::vector<mpfr::mpreal> &x, std::vector<Band> &bands,
        mp_prec_t prec = 165ul);

/*! Computes the frequency response of the current filter
 * @param[out] Pc the frequency response amplitude value at the current node
 * @param[in] xVal the current frequency node where we do our computation
 * (the point is given in the \f$\left[-1,1\right]\f$ interval,
 * and not the initial \f$\left[0,\pi\right]\f$)
 * @param[in] x the current reference set
 * @param[in] C the frequency responses at the current reference set
 * @param[in] w the current barycentric weights
 * @param[in] prec numerical precision used for the computations
 */
void computeApprox(mpfr::mpreal &Pc, const mpfr::mpreal &xVal,
        std::vector<mpfr::mpreal> &x, std::vector<mpfr::mpreal> &C,
        std::vector<mpfr::mpreal> & w,
        mp_prec_t prec = 165ul);

/*! Computes the approximation error at a given node using the current set of
 * reference points
 * @param[out] error the requested error value
 * @param[in] xVal the current frequency node where we do our computation
 * @param[in] delta the current reference error
 * @param[in] x the current reference set
 * @param[in] C the frequency response values at the x nodes
 * @param[in] w the barycentric weights
 * @param[in] bands frequency band information for the ideal filter
 * @param[in] prec numerical precision used for the computations
 */
void computeError(mpfr::mpreal &error, const mpfr::mpreal &xVal,
        mpfr::mpreal &delta, std::vector<mpfr::mpreal> &x,
        std::vector<mpfr::mpreal> &C, std::vector<mpfr::mpreal> &w,
        std::vector<Band> &bands,
        mp_prec_t prec = 165ul);

/*! The ideal frequency response and weight information at the given frequency
 * node (it can be in the \f$\left[-1,1\right]\f$ interval,
 * and not the initial \f$\left[0,\pi\right]\f$, the difference is made with
 * information from the bands parameter)
 * @param[out] D ideal frequency response
 * @param[out] W weight value for the current point
 * @param[in] xVal the current frequency node where we do our computation
 * @param[in] bands frequency band information for the ideal filter
 */
void computeIdealResponseAndWeight(mpfr::mpreal &D, mpfr::mpreal &W,
        const mpfr::mpreal &xVal, std::vector<Band> &bands);


#endif // BARYCENTRIC_H
