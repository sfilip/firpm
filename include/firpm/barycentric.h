/**
 * @file barycentric.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief Functions used to perform barycentric Lagrange interpolation
 * during the execution of the Parks-McClellan algorithm
 *
 */

//    firpm
//    Copyright (C) 2015 - 2019  S. Filip
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

#ifndef __PMBARYCENTRIC_H__
#define __PMBARYCENTRIC_H__

#include "util.h"
#include "cheby.h"
#include "band.h"

namespace pm {

        /*! Procedure which computes the weights used in
        * the evaluation of the barycentric interpolation
        * formulas (see [Berrut&Trefethen2004] and [Pachon&Trefethen2009]
        * for the implementation ideas)
        * @param[in] x the interpolation points
        * @return the computed weights
        */
        template<typename T>
        std::vector<T> baryweights(std::vector<T> const& x);

        /*! Determines the current reference error according to the
        * barycentric formula (internally it also computes the barycentric weights)
        * @param[in] x the current reference set (i.e., interpolation points)
        * @param[in] bands information relating to the ideal frequency response of
        * the filter
        * @return delta the value of the current reference error
        */
        template<typename T>
        T compdelta(std::vector<T>& x,
                std::vector<band_t<T>> const& bands);

        /*! Determines the current reference error according to the
        * barycentric formula
        * @param[in] w the barycentric weights associated with the current reference
        * set
        * @param[in] x the current reference set (i.e. interpolation points)
        * @param[in] bands information relating to the ideal frequency response of
        * the filter
        * @return the value of the current reference error
        */

        template<typename T>
        T compdelta(std::vector<T>& w,
                std::vector<T>& x, std::vector<band_t<T>> const& bands);

        /*! Computes the filter response at the current reference set
        * @param[in] delta the current reference error
        * @param[in] x the current reference vector
        * @param[in] bands frequency band information for the ideal filter
        * @return the vector of frequency responses at the reference set
        */
        template<typename T>
        std::vector<T> compc(T& delta,
                std::vector<T>& x, std::vector<band_t<T>> const& bands);

        /*! Computes the frequency response of the current filter
        * @param[in] xVal the current frequency node where we do our computation
        * (the point is given in the \f$\left[-1,1\right]\f$ interval,
        * and not the initial \f$\left[0,\pi\right]\f$)
        * @param[in] x the current reference set
        * @param[in] C the frequency responses at the current reference set
        * @param[in] w the current barycentric weights
        * @return the frequency response amplitude value at the current node
        */
        template<typename T>
        T approx(T const& xVal,
                std::vector<T> const& x,
                std::vector<T> const& C,
                std::vector<T> const& w);

        /*! Computes the approximation error at a given node using the current set of
        * reference points
        * @param[in] xVal the current frequency node where we do our computation
        * @param[in] delta the current reference error
        * @param[in] x the current reference set
        * @param[in] C the frequency response values at the x nodes
        * @param[in] w the barycentric weights
        * @param[in] bands frequency band information for the ideal filter
        * @return error the requested error value
        */
        template<typename T>
        T comperror(T const& xVal,
                T& delta,
                std::vector<T> const& x,
                std::vector<T> const& C,
                std::vector<T> const& w,
                std::vector<band_t<T>> const& bands);

        /*! The ideal frequency response and weight information at the given frequency
        * node (it can be in the \f$\left[-1,1\right]\f$ interval,
        * and not the initial \f$\left[0,\pi\right]\f$, the difference is made with
        * information from the bands parameter)
        * @param[in] xVal the current frequency node where we do our computation
        * @param[in] bands frequency band information for the ideal filter
        * @return a (D, W) tuple containing ideal frequency response and weight
        */
        template<typename T>
        std::pair<T, T> idealvals(T const& xVal, std::vector<band_t<T>>& bands);

} // namespace pm

#endif
