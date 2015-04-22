/**
 * @file eigenvalue.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief Functions used to generate eigenvalues of colleague matrices
 *
 */

#ifndef EIGENVALUE_H
#define EIGENVALUE_H

#include "util.h"
#include "cheby.h"
#include "barycentric.h"

/** Eigen matrix container for double values */
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXq;
/** Eigen vector container for double values */
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> VectorXcq;

/*! Funtion that generates the colleague matrix for a Chebyshev interpolant
 * expressed using the basis of Chebyshev polynomials of the first kind
 * (WITHOUT balancing)
 * @param[out] C the corresponding colleague matrix
 * @param[in] a the coefficients of the Chebyshev interpolant
 */
void generateColleagueMatrix1stKind(MatrixXq& C,
        std::vector<double>& a);

/*! Funtion that generates the colleague matrix for a Chebyshev interpolant
 * expressed using the basis of Chebyshev polynomials of the first kind
 * (WITH balancing)
 * @param[out] C the corresponding colleague matrix
 * @param[in] a the coefficients of the Chebyshev interpolant
 */

void generateColleagueMatrix1stKindWithBalancing(MatrixXq& C,
        std::vector<double>& a);

/*! Funtion that generates the colleague matrix for a Chebyshev interpolant
 * expressed using the basis of Chebyshev polynomials of the second kind
 * (WITHOUT balancing)
 * @param[out] C the corresponding colleague matrix
 * @param[in] a the coefficients of the Chebyshev interpolant
 */
void generateColleagueMatrix2ndKind(MatrixXq& C,
        std::vector<double>& a);

/*! Funtion that generates the colleague matrix for a Chebyshev interpolant
 * expressed using the basis of Chebyshev polynomials of the second kind
 * (WITH balancing)
 * @param[out] C the corresponding colleague matrix
 * @param[in] a the coefficients of the Chebyshev interpolant
 */
void generateColleagueMatrix2ndKindWithBalancing(MatrixXq& C,
        std::vector<double>& a);

/*! Function that computes the eigenvalues of a given matrix
 * @param[out] eigenvalues the computed eigenvalues
 * @param[in] C the corresponding matrix
 */
void determineEigenvalues(VectorXcq &eigenvalues,
        MatrixXq &C);

/*! Function that computes the real values located inside
 * an interval \f$[a,b]\f$ from a vector of complex values
 * @param[out] realValues the set of real values inside \f$[a,b]\f$
 * @param[in] complexValues the set of complex values to search from
 * @param[in] a left side of the closed interval
 * @param[in] b right side of the closed interval
 */
void getRealValues(std::vector<double> &realValues,
        VectorXcq &complexValues,
        double &a, double &b);


#endif
