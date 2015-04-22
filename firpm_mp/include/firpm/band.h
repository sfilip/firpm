/**
 * @file band.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief Utilities for managing frequency bands
 * */


#ifndef BAND_H_
#define BAND_H_

#include "util.h"

/**
 * The interval type being considered (in \f$\left[0, \pi\right]\f$ or in \f$\left[-1, 1\right]\f$)
 */
enum BandSpace {
    FREQ,           /**< not done the change of variable (i.e. we are in \f$\left[0, \pi\right]\f$) */
    CHEBY           /**< done the change of variable (i.e. we are in \f$\left[-1, 1\right]\f$) */
};

/**
 * @brief A data type encapsulating information relevant to a frequency band
 *
 * Contains important information concerning a frequency band used during the
 * execution of the exchange algorithm.
 */
struct Band {
    BandSpace space;        /**< the space in which we are working */
    std::function<mpfr::mpreal(BandSpace, mpfr::mpreal)> amplitude;
                            /**< the ideal amplitude for this band */
    mpfr::mpreal start;           /**< the left bound of the band */
    mpfr::mpreal stop;            /**< the right bound of the band */
    std::function<mpfr::mpreal(BandSpace, mpfr::mpreal)> weight;
                            /**< weight function value on the band */
    std::size_t extremas;   /**< number of interpolation points taken in the band */
};

/**
 * Gives the direction in which the change of variable is performed
 */
enum ConversionDirection {
    FROMFREQ,               /**< apply the change of variable \f$y=\cos(x)\f$*/
    TOFREQ                  /**< apply the change of variable \f$y=\arccos(x)\f$*/
};

/*! Performs the change of variable on the set of bands of interest
 * @param[out] out output frequency bands
 * @param[in]  in input frequency bands
 * @param[in]  direction the direction in which the change of variable is performed
 * @param[in]  prec the working precision used to perform the computations
 */
void bandConversion(std::vector<Band> &out, std::vector<Band> &in,
        ConversionDirection direction, mp_prec_t prec = 165ul);

#endif /* BAND_H_ */
