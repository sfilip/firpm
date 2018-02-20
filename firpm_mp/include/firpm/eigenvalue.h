/**
 * @file eigenvalue.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief Functions used to generate eigenvalues of colleague matrices
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



#ifndef EIGENVALUE_H_
#define EIGENVALUE_H_

#include "util.h"
#include "cheby.h"
#include "barycentric.h"

/** Eigen matrix container for mpfr::mpreal values */
typedef Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> MatrixXq;
/** Eigen vector container for mpfr::mpreal values */
typedef Eigen::Matrix<std::complex<mpfr::mpreal>, Eigen::Dynamic, 1> VectorXcq;

/*! Funtion that generates the colleague matrix for a Chebyshev interpolant
 * expressed using the basis of Chebyshev polynomials of the first kind
 * (WITH or WITHOUT balancing in the vein of [Parlett&Reinsch1969] "Balancing a Matrix for
 * Calculation of Eigenvalues and Eigenvectors")
 * @param[out] C the corresponding colleague matrix
 * @param[in] a the coefficients of the Chebyshev interpolant
 * @param[in] withBalancing perform a balancing operation on the colleague matrix C
 * @param[in] prec the working precision used for the computations
 */
void FIRPM_EXPORT generateColleagueMatrix1stKind(MatrixXq& C,
        std::vector<mpfr::mpreal>& a,
        bool withBalancing,
        mp_prec_t prec = 165ul);

void FIRPM_EXPORT generateColleagueMatrix2ndKind(MatrixXq& C,
        std::vector<mpfr::mpreal>& a,
        bool withBalancing,
        mp_prec_t prec = 165ul);

/*! Funtion that generates the colleague matrix for a Chebyshev interpolant
 * expressed using the basis of Chebyshev polynomials of the second kind
 * (WITH or WITHOUT balancing in the vein of [Parlett&Reinsch1969] "Balancing a Matrix for
 * Calculation of Eigenvalues and Eigenvectors")
 * @param[out] C the corresponding colleague matrix
 * @param[in] a the coefficients of the Chebyshev interpolant
 * @param[in] withBalancing perform a balancing operation on the colleague matrix C
 * @param[in] prec the working precision used for the computations
 */
void FIRPM_EXPORT generateColleagueMatrix2ndKindWithBalancing(MatrixXq& C,
        std::vector<mpfr::mpreal>& a,
        mp_prec_t prec = 165ul);

/*! Function that computes the eigenvalues of a given matrix
 * @param[out] eigenvalues the computed eigenvalues
 * @param[in] C the corresponding matrix
 */
void FIRPM_EXPORT determineEigenvalues(VectorXcq &eigenvalues,
        MatrixXq &C);

/*! Function that computes the real values located inside
 * an interval \f$[a,b]\f$ from a vector of complex values
 * @param[out] realValues the set of real values inside \f$[a,b]\f$
 * @param[in] complexValues the set of complex values to search from
 * @param[in] a left side of the closed interval
 * @param[in] b right side of the closed interval
 */
void FIRPM_EXPORT getRealValues(std::vector<mpfr::mpreal> &realValues,
        VectorXcq &complexValues,
        mpfr::mpreal &a, mpfr::mpreal &b);


#endif
