/**
 * @file barycentric.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief Functions used to perform barycentric Lagrange interpolation
 * during the execution of the Parks-McClellan algorithm
 *
 */

//    firpm
//    Copyright (C) 2015 - 2024  S. Filip

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
        * @param[out] w the computed weights
        * @param[in] x the interpolation points
        */
        template<typename T>
        void baryweights(std::vector<T>& w,
                std::vector<T>& x);

        /*! Determines the current reference error according to the
        * barycentric formula (internally it also computes the barycentric weights)
        * @param[out] delta the value of the current reference error
        * @param[in] x the current reference set (i.e., interpolation points)
        * @param[in] bands information relating to the ideal frequency response of
        * the filter
        */
        template<typename T>
        void compdelta(T& delta, std::vector<T>& x,
                std::vector<band_t<T>>& bands);

        /*! Determines the current reference error according to the
        * barycentric formula
        * @param[out] delta the value of the current reference error
        * @param[in] w the barycentric weights associated with the current reference
        * set
        * @param[in] x the current reference set (i.e. interpolation points)
        * @param[in] bands information relating to the ideal frequency response of
        * the filter
        */

        template<typename T>
        void compdelta(T& delta, std::vector<T>& w,
                std::vector<T>& x, std::vector<band_t<T>>& bands);

        /*! Computes the filter response at the current reference set
        * @param[out] C the vector of frequency responses at the reference set
        * @param[in] delta the current reference error
        * @param[in] x the current reference vector
        * @param[in] bands frequency band information for the ideal filter
        */
        template<typename T>
        void compc(std::vector<T>& C, T& delta,
                std::vector<T>& x, std::vector<band_t<T>>& bands);

        /*! Computes the frequency response of the current filter
        * @param[out] Pc the frequency response amplitude value at the current node
        * @param[in] xVal the current frequency node where we do our computation
        * (the point is given in the \f$\left[-1,1\right]\f$ interval,
        * and not the initial \f$\left[0,\pi\right]\f$)
        * @param[in] x the current reference set
        * @param[in] C the frequency responses at the current reference set
        * @param[in] w the current barycentric weights
        */
        template<typename T>
        void approx(T& Pc, T const& xVal,
                std::vector<T>& x, std::vector<T>& C,
                std::vector<T>& w);

        /*! Computes the approximation error at a given node using the current set of
        * reference points
        * @param[out] error the requested error value
        * @param[in] xVal the current frequency node where we do our computation
        * @param[in] delta the current reference error
        * @param[in] x the current reference set
        * @param[in] C the frequency response values at the x nodes
        * @param[in] w the barycentric weights
        * @param[in] bands frequency band information for the ideal filter
        */
        template<typename T>
        void comperror(T& error, T const& xVal,
                T& delta, std::vector<T>& x,
                std::vector<T>& C, std::vector<T>& w,
                std::vector<band_t<T>>& bands);

        /*! The ideal frequency response and weight information at the given frequency
        * node (it can be in the \f$\left[-1,1\right]\f$ interval,
        * and not the initial \f$\left[0,\pi\right]\f$, the difference is made with
        * information from the bands parameter)
        * @param[out] D ideal frequency response
        * @param[out] W weight value for the current point
        * @param[in] xVal the current frequency node where we do our computation
        * @param[in] bands frequency band information for the ideal filter
        */
        template<typename T>
        void idealvals(T& D, T& W,
                T const& xVal, std::vector<band_t<T>>& bands);

} // namespace pm

#endif