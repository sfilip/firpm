/**
 * @file firpm.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief The file to include for accessing the library functions
 *
 * This file contains the central routines for constructing minimax FIR filters.
 */


//    firpm
//    Copyright (C) 2015-2020  S. Filip
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

#ifndef __PMPM_H__
#define __PMPM_H__

#include "util.h"
#include "cheby.h"
#include "barycentric.h"

namespace pm {
    /** @enum filter_t marker to distinguish
     * between the two categories of filters (digital
     * differentiators and Hilbert transformers) that
     * can be constructed using type III and IV FIR
     * filters. */
    enum class filter_t {
        FIR_DIFFERENTIATOR,     /**< marker for constructing digital differentiators */
        FIR_HILBERT             /**< marker for constructing Hilbert transformers */
    };

    /** @enum init_t flag representing the
     * initialization strategies that can be used
     * at the lowest level of the scaling approach 
     */
    enum class init_t {
        UNIFORM,                /**< uniform initialization marker */
        SCALING,                /**< reference scaling-based initialization */
        AFP                     /**< AFP algorithm-based initialization */
    };

    /** @enum status_t code to distinguish the
     * various states in which the Parks-McClellan
     * algorithm execution finished in. */
    enum class status_t {
        STATUS_SUCCESS,                     /**< successful execution */
        STATUS_FREQUENCY_INVALID_INTERVAL,  /**< invalid frequency inputs */
        STATUS_AMPLITUDE_VECTOR_MISMATCH,   /**< amplitude/frequency vector sizes mismatch */
        STATUS_AMPLITUDE_DISCONTINUITY,     /**< discontinuous amplitude values detected */
        STATUS_WEIGHT_NEGATIVE,             /**< negative weight value detected */
        STATUS_WEIGHT_VECTOR_MISMATCH,      /**< weight/frequency vector sizes mismatch */
        STATUS_WEIGHT_DISCONTINUITY,        /**< discontinuous weight value detected */
        STATUS_SCALING_INVALID,             /**< failure in performing valid reference scaling */
        STATUS_AFP_INVALID,                 /**< numerical failure in performing AFP initialization */ 
        STATUS_COEFFICIENT_SET_INVALID,     /**< invalid final coefficient set */
        STATUS_EXCHANGE_FAILURE,            /**< runtime error in producing a valid reference set */
        STATUS_CONVERGENCE_WARNING,         /**< successful execution, but with convengence warnings */
        STATUS_UNKNOWN_FAILURE              /**< unknown runtime failure */
    };

    /**
     * @brief The type of the object returned by the Parks-McClellan algorithm.
     *
     * Utility object which contains useful information about the filter computed by the
     * Parks-McClellan algorithm
     */
    template<typename T>
    struct pmoutput_t
    {
        std::vector<T> h;           /**< the final filter coefficients*/
        std::vector<T> x;           /**< the reference set used to generate the final
                                    filter (values are in \f$[-1,1]\f$ and NOT 
                                    \f$[0,\pi]\f$)*/
        std::size_t iter;           /**< number of iterations that were necessary to
                                    achieve convergence*/
        T delta;                    /**< the final reference error */
        T q;                        /**< convergence parameter value */
        status_t status;            /**< status code for the output object */
    };

    /*! An implementation of the uniform initialization approach for
    * starting the Parks-McClellan algorithm
    * @param[out] omega the initial set of references to be computed
    * @param[in] B the frequency bands of interest (i.e., stopbands and
    * passbands for example)
    * @param[in] n the size of the reference set
    */

    template<typename T>
    void uniform(std::vector<T>& omega,
            std::vector<band_t<T>>& B, std::size_t n);

    /*! An implementation of the reference scaling approach mentioned
    * in section 4 of the article.
    * @param[out] status saves diagnostic information in case the routine does not
    * complete successfully
    * @param[out] nx the reference set obtained from scaling the initial set x
    * @param[out] ncbands contains information about the bands of interest 
    * corresponding to the new reference (i.e., how many reference points are 
    * inside each band). The bands are given inside \f$[-1,1]\f$ (i.e., the CHEBY 
    * band space)
    * @param[out] nfbands contains information about the bands of interest corresponding 
    * to the new reference (i.e., how many reference points are inside each band). 
    * The bands are given inside \f$[0,\pi]\f$ (i.e., the FREQ band space)
    * @param[in] nxs the size of the new reference set
    * @param[in] x the input reference on which will be used to perform the scaling
    * @param[in] cbands band information for the filter to which the x reference 
    * corresponds to. The bands are given inside \f$[-1,1]\f$ (i.e., the CHEBY 
    * band space)
    * @param[in] fbands band information for the filter to which the x reference 
    * corresponds to. The bands are given inside \f$[0,\pi]\f$ (i.e., the FREQ 
    * band space)
    */
    template<typename T>
    void refscaling(status_t& status,
            std::vector<T>& nx, std::vector<band_t<T>>& ncbands,
            std::vector<band_t<T>>& nfbands, std::size_t nxs,
            std::vector<T>& x, std::vector<band_t<T>>& cbands,
            std::vector<band_t<T>>& fbands);

    /*! An internal routine which implements the exchange algorithm for designing 
    * FIR filters
    * @param[in] x the initial reference set
    * @param[in] cbands band information for the filter to which the x reference 
    * corresponds to. The bands are given inside \f$[-1,1]\f$ (i.e., the CHEBY band space)
    * @param[in] eps convergence parameter threshold (i.e., quantizes the number of 
    * significant digits of the minimax error that are accurate at the end of the 
    * final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last 
    * iteration. Since there is no explicit filter type (I to IV) given as input, 
    * the h vector of the output will correspond to the coefficients of the frequency 
    * response \f$H_d(\omega)=\sum_{k=0}^{n}h_k\cos(\omega k)\f$
    *
    * The following piece of code shows how to use this function to construct a lowpass 
    * filter. It also shows how to write frequency band specifications and how to use 
    * the uniform initialization routine. It is equivalent to the sample code from 
    * firpm. <b>Important:</b> in order to determine the final coefficients of the 
    * transfer function of the filter, some preprocessing specific to the type of the 
    * filter (I to IV) has to be done on the elements of <tt>output.h</tt> (see for 
    * example the source code of the <tt>firpm</tt> functions on how this is done for 
    * each type of filter).
    * @see firpm
    * @code
    * // frequency band specification
    * std::vector<band_t<double>> fbands(2);
    * double pi = M_PI;
    *
    * fbands[0].start = 0;
    * fbands[0].stop = pi * 0.4;
    * fbands[0].weight = [] (space_t, double) -> double {return 1.0; };
    * fbands[0].space = space_t::FREQ;
    * fbands[0].amplitude = [](space_t, double) -> double { return 1.0; };
    *
    * fbands[1].start = pi * 0.5;
    * fbands[1].stop = pi;
    * fbands[1].weight = [] (space_t, double) -> double {return 10.0; };
    * fbands[1].space = BandSpace::FREQ;
    * fbands[1].amplitude = [](space_t, double) -> double { return 10.0; };
    * std::size_t degree = 100;  // filter degree
    *
    * // reference initialization code
    * std::vector<band_t<double>> cbands;
    * std::vector<double> omega;
    * std::vector<double> x;
    * uniform(omega, fbands, degree + 2u);
    * // apply the change of variable y = cos(x) so that we are working inside [-1, 1]
    * cos(x, omega);
    * bandconv(cbands, fbands, convdir_t::FROMFREQ);
    * // apply the exchange algorithm
    * pmoutput_t<double> output = exchange(x, cbands);
    * @endcode
    */

    template<typename T>
    pmoutput_t<T> exchange(std::vector<T>& x,
            std::vector<band_t<T>>& cbands,
            double eps = 0.01,
            std::size_t nmax = 4u, 
            unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type I and II FIR filters.
    * This routine is the most general and can be set to use any of the three
    * proposed initialization types. If the default parameters are used, then
    * it will use uniform initialization.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final 
    * transfer function. For even n, the filter will be type I, while for odd 
    * n the type is II.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] eps convergence parameter threshold (i.e., quantizes the number 
    * of significant digits of the minimax error that are accurate at the end of 
    * the final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] strategy initialization strategy. Can be UNIFORM, SCALING or AFP
    * @param[in] depth in case the SCALING initialization strategy is used, 
    * specifies the number of scaling levels to use (by default, the value is set
    * to 1, meaning a filter of length approximatively n/2 is used to construct
    * the initial reference for the requested n coefficient filter)
    * @param[in] rstrategy in case SCALING is used, specifies how to initialize
    * the smallest length filter used to perform reference scaling (UNIFORM by
    * default)
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last 
    * iteration. The h vector of the output contains the coefficients corresponding 
    * to the transfer function of the final filter (in this case, for types I and II,
    * the values are symmetrical to the middle coefficient(s))
    *
    * An example of how to use this function is given below. It designs a degree 
    * \f$100\f$ type I lowpass filter, with passband \f$[0, 0.4\pi]\f$ and stopband 
    * \f$[0.5\pi, \pi]\f$. It has unit weight inside the passband and weight 10
    * inside the stopband. More examples, including code on how to use the customized
    * reference scaling and AFP versions <tt>firpmRS, firpmAFP</tt>, are provided inside 
    * the test files.
    * @code
    * pmoutput_t<double> output = firpm<double>(200, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    * @endcode
    */

    template<typename T>
    pmoutput_t<T> firpm(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                double eps = 0.01,
                std::size_t nmax = 4u,
                init_t strategy = init_t::UNIFORM,
                std::size_t depth = 0u,
                init_t rstrategy = init_t::UNIFORM,
                unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type I and II FIR filters. This routine uses 
    * reference scaling by default and is just a wrapper over the <tt>firpm</tt>.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final transfer function. 
    * For even n, the filter will be type I, while for odd n the type is II.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] eps convergence parameter threshold (i.e., quantizes the number of 
    * significant digits of the minimax error that are accurate at the end of the 
    * final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] depth how many times should reference scaling be applied 
    * recursively (default value is 1)
    * @param[in] rstrategy  what initialization strategy to use at the lowest level 
    * (uniform or AFP-based) if strategy is reference scaling 
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last iteration. 
    * The h vector of the output contains the coefficients corresponding to the transfer 
    * function of the final filter (in this case, for types I and II, the values are 
    * symmetrical to the middle coefficient(s))
    */

    template<typename T>
    pmoutput_t<T> firpmRS(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                double eps = 0.01,
                std::size_t nmax = 4u,
                std::size_t depth = 1u,
                init_t rstrategy = init_t::UNIFORM,
                unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type I and II FIR filters. This routine 
    * uses AFP-based initialization and is just a wrapper over <tt>firpm</tt>.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final transfer function. 
    * For even n, the filter will be type I, while for odd n the type is II.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] eps convergence parameter threshold (i.e., quantizes the number of 
    * significant digits of the minimax error that are accurate at the end of the 
    * final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last iteration. 
    * The h vector of the output contains the coefficients corresponding to the transfer 
    * function of the final filter (in this case, for types I and II, the values are 
    * symmetrical to the middle coefficient(s))
    */

    template<typename T>
    pmoutput_t<T> firpmAFP(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                double eps = 0.01,
                std::size_t nmax = 4u,
                unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type III and IV FIR filters. 
    * This routine is the most general and can be set to use any of the three
    * proposed initialization types. If the default parameters are used, then
    * it will use uniform initialization.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final 
    * transfer function. For even n, the filter will be type III, while for odd 
    * n the type is IV.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] type denotes the type of filter we want to design: digital 
    * differentiator of Hilbert transformer
    * @param[in] eps convergence parameter threshold (i.e., quantizes the number 
    * of significant digits of the minimax error that are accurate at the end of 
    * the final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] strategy the initialization strategy
    * @param[in] depth how many times should reference scaling be applied 
    * recursively (default value is 1)
    * @param[in] rstrategy  what initialization strategy to use at the lowest level 
    * (uniform or AFP-based) if strategy is reference scaling 
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last iteration. 
    * The h vector of the output contains the coefficients corresponding to the 
    * transfer function of the final filter (in this case, for types III and IV, 
    * the values are antisymmetrical to the middle coefficient(s))*/

    template<typename T>
    pmoutput_t<T> firpm(std::size_t n,
            std::vector<T>const& f,
            std::vector<T>const& a,
            std::vector<T>const& w,
            filter_t type,
            double eps = 0.01,
            std::size_t nmax = 4,
            init_t strategy = init_t::UNIFORM,
            std::size_t depth = 0u,
            init_t rstrategy = init_t::UNIFORM,
            unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type III and IV FIR filters. 
    * This routine uses reference scaling.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final 
    * transfer function. For even n, the filter will be type III, while for 
    * odd n the type is IV.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] type denotes the type of filter we want to design: digital 
    * differentiator of Hilbert transformer
    * @param[in] eps convergence parameter threshold (i.e., quantizes the 
    * number of significant digits of the minimax error that are accurate 
    * at the end of the final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] depth how many times should reference scaling be applied 
    * recursively (default value is 1)
    * @param[in] rstrategy  what initialization strategy to use at the 
    * lowest level (uniform or AFP-based)
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last 
    * iteration. The h vector of the output contains the coefficients corresponding 
    * to the transfer function of the final filter (in this case, for types III and 
    * IV, the values are antisymmetrical to the middle coefficient(s))*/

    template<typename T>
    pmoutput_t<T> firpmRS(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                filter_t type,
                double eps = 0.01,
                std::size_t nmax = 4u,
                std::size_t depth = 1u,
                init_t rstrategy = init_t::UNIFORM,
                unsigned long prec = 165ul);

    /*! Parks-McClellan routine for implementing type III and IV FIR filters.
    * This routine uses AFP-based initialization.
    * @param[in] n \f$n+1\f$ denotes the number of coefficients of the final 
    * transfer function. For even n, the filter will be type III, while for 
    * odd n the type is IV.
    * @param[in] f vector denoting the frequency ranges of each band of interest
    * @param[in] a the ideal amplitude at each point of f
    * @param[in] w the wight function value on each band
    * @param[in] type denotes the type of filter we want to design: digital 
    * differentiator of Hilbert transformer
    * @param[in] eps convergence parameter threshold (i.e quantizes the number 
    * of significant digits of the minimax error that are accurate at the end 
    * of the final iteration)
    * @param[in] nmax the degree used by the CPR method on each subinterval
    * @param[in] prec the numerical precision of the MPFR type (will be disregarded for
    * the double and long double instantiations of the functions)
    * @return information pertaining to the polynomial computed at the last 
    * iteration. The h vector of the output contains the coefficients corresponding 
    * to the transfer function of the final filter (in this case, for types III and 
    * IV, the values are antisymmetrical to the middle coefficient(s))*/

    template<typename T>
    pmoutput_t<T> firpmAFP(std::size_t n,
                std::vector<T>const& f,
                std::vector<T>const& a,
                std::vector<T>const& w,
                filter_t type,
                double eps = 0.01,
                std::size_t nmax = 4u,
                unsigned long prec = 165ul);

} // namespace pm

#endif
