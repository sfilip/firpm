/**
 * @file cheby.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief Useful functions for working with Chebyshev interpolants
 *
 */

//    firpm_d
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

#ifndef CHEBY_H_
#define CHEBY_H_

#include "util.h"

/*! Computes the cosines of the elements of a vector
 * @param[in] in the vector to process
 * @param[out] out the vector containing the cosines of the elements from the
 * vector in
 */
void applyCos(std::vector<double>& out,
        std::vector<double> const& in);

/*! Does a change of variable from the interval \f$\left[-1, 1\right]\f$ to the
 * interval \f$\left[a, b\right]\f$ on the elements of a given input vector
 * @param[in] in the vector of elements in \f$\left[-1, 1\right]\f$
 * @param[out] out the vector of elements after the change of variable
 * @param[in] a left bound of the desired interval
 * @param[in] b right bound of the desired interval
 */
void changeOfVariable(std::vector<double>& out,
        std::vector<double> const& in,
        double& a, double& b);

/*! The Clenshaw algorithm which evaluates the value of a Chebyshev interpolant
 *  (CI) at a certain point
 *  @param[out] result reference to the variable that will contain the
 *  evaluation of the CI at the specified point (parameter x)
 *  @param[in] p a vector containing the coefficients of the CI
 *  @param[in] x the point at which we want to compute the value
 *  of the CI
 *  @param[in] a left bound of the interval where the CI is considered
 *  @param[in] b right bound of the interval where the CI is considered
 */
void evaluateClenshaw(double &result, std::vector<double> &p,
                        double &x, double &a, double &b);

/*! The Clenshaw algorithm which evaluates the value of a CI implicitly
 *  considered on the \f$\left[-1,1\right]\f$ interval
 *  @param[out] result reference to the variable that will contain the
 *  evaluation of the CI at the specified point (parameter x)
 *  @param[in] p a vector containing the coefficients of the CI
 *  @param[in] x the point at which we want to compute the value
 *  of the CI
 */
void evaluateClenshaw(double &result, std::vector<double> &p,
                        double &x);

/*! The Clenshaw algorithm which evaluates the value of a CI expressed
 *  using a basis consisting of Chebyshev polynomials of the second kind.
 *  The working interval is considered to be \f$\left[-1,1\right]\f$.
 *  @param[out] result reference to the variable that will contain the
 *  evaluation of the CI at the specified point (parameter x)
 *  @param[in] p a vector containing the coefficients of the CI
 *  @param[in] x the point at which we want to compute the value
 *  of the CI
 */
void evaluateClenshaw2ndKind(double &result, std::vector<double> &p,
                        double &x);


/*! Function that generates equidistant nodes inside the
 * \f$\left[0,\pi\right]\f$ interval, meaning values of the
 * form \f$\frac{i\pi}{n}\f$, where \f$0\leq i\leq n\f$
 * @param[out] v the vector that will contain the equi-distributed points
 * @param[in] n the number of points - 1 which will be computed
 */
void generateEquidistantNodes(std::vector<double>& v, std::size_t n);

/*! Function that generates the Chebyshev nodes of the second kind
 * \f$\mu_k=\cos\left(\frac{k\pi}{n}\right), k=0,\ldots,n\f$.
 * @param[out] v the vector that will contain the Chebyshev nodes
 * @param[in] n the number of points - 1 which will be computed
 */
void generateChebyshevPoints(std::vector<double>& v, std::size_t n);


/*! This function computes the values of the coefficients of the CI when
 * Chebyshev nodes of the second kind are used
 * @param[out] c vector used to hold the values of the computed coefficients
 * @param[in] fv vector that holds the value of the current function to
 * approximate at the Chebyshev nodes of the second kind scaled to the
 * appropriate interval (in our case it will be \f$\left[0,\pi\right]\f$)
 * @param[in] n degree of the CI
 */
void generateChebyshevCoefficients(std::vector<double>& c,
                std::vector<double>& fv, std::size_t n);

/*! Function that generates the coefficients of the derivative of a given CI
 *  @param[out] derivC the vector of coefficients of the derivative of the CI
 *  @param[in] c the vector of coefficients of the CI whose derivative we
 *  want to compute
 */
void derivativeCoefficients1stKind(std::vector<double>& derivC,
                        std::vector<double>& c);

/*! Function that generates the coefficients of the derivative of a given
 *  CI, but expressed in the orthogonal basis of Chebyshev polynomials of
 *  the second kind.
 *  @param[out] derivC the vector of coefficients of the derivative of the CI
 *  @param[in] c the vector of coefficients of the CI whose derivative we
 *  want to compute
 */
void derivativeCoefficients2ndKind(std::vector<double>& derivC,
        std::vector<double>& c);

#endif /* CHEBY_H_ */
