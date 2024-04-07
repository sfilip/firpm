/**
 * @file band.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief Utilities for managing frequency bands
 * */

//    firpm
//    Copyright (C) 2015 - 2024  S. Filip

#ifndef __PMBAND_H__
#define __PMBAND_H__

#include "util.h"

namespace pm {
    /**
     * The interval type being considered (in \f$\left[0, \pi\right]\f$
     * or in \f$\left[-1, 1\right]\f$)
     */
    enum space_t {
        FREQ,           /**< not done the change of variable
                        (i.e., we are in \f$\left[0, \pi\right]\f$) */
        CHEBY           /**< done the change of variable
                        (i.e., we are in \f$\left[-1, 1\right]\f$) */
    };

    /**
     * @brief A data type encapsulating information relevant to a
     * frequency band
     *
     * Contains important information concerning a frequency band
     * used during the execution of the exchange algorithm.
     */
    template<typename T>
    struct band_t {
        space_t space;          /**< the space in which we are working */
        std::function<T(space_t, T)> amplitude;
                                /**< the ideal amplitude for this band */
        T start;           /**< the left bound of the band */
        T stop;            /**< the right bound of the band */
        std::function<T(space_t, T)> weight;
                                /**< weight function value on the band */
        std::size_t xs;         /**< number of interpolation points taken in the band */
        std::vector<T> part;    /**< partition points (if any) inside the band */
    };

    /**
     * Gives the direction in which the change of variable is performed
     */
    enum convdir_t {
        FROMFREQ,               /**< apply the change of variable \f$y=\cos(x)\f$*/
        TOFREQ                  /**< apply the change of variable \f$y=\arccos(x)\f$*/
    };

    /*! Performs the change of variable on the set of bands of interest
    * @param[out] out output frequency bands
    * @param[in]  in input frequency bands
    * @param[in]  direction the direction in which the change of
    * variable is performed
    */
    template<typename T>
    void bandconv(std::vector<band_t<T>>& out, std::vector<band_t<T>>& in,
            convdir_t direction);

} // namespace pm

#endif