/**
 * @file band.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief Utilities for managing frequency bands
 * */

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
    * @param[in]  in input frequency bands
    * @param[in]  direction the direction in which the change of
    * variable is performed
    * @return output frequency bands
    */
    template<typename T>
    std::vector<band_t<T>> bandconv(std::vector<band_t<T>> const &in,
            convdir_t direction);

} // namespace pm

#endif
