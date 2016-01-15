/**
 * @file pm.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief The routines the user should call to construct FIR filters
 *
 * This file contains the central routines for constructing minimax FIR filters.
 */


//    firpm_d
//    Copyright (C) 2015-2016  S. Filip
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

#ifndef PM_H
#define PM_H

#include "util.h"
#include "cheby.h"
#include "barycentric.h"
#include "eigenvalue.h"

/** utility type for storing interval endpoints */
typedef std::pair<double, double> Interval;

/** @enum ftype marker to distinguish
 * between the two categories of filters (digital
 * differentiators and Hilbert transformers) that
 * can be constructed using type III and IV FIR
 * filters. */
enum class ftype {
    FIR_DIFFERENTIATOR,     /**< marker for constructing digital differentiators */
    FIR_HILBERT             /**< marker for constructing Hilbert transformers */
};

/** @enum RootSolver flag representing the
 * initialization strategis that can be used
 * at the lowest level of the scaling approach
 * */
enum class RootSolver {
    UNIFORM,                /**< uniform initialization marker */
    AFP                     /**< AFP algorithm-based initialization */
};


/**
 * @brief The type of the object returned by the Parks-McClellan algorithm.
 *
 * Utility object which contains useful information about the filter computed by the
 * Parks-McClellan algorithm
 */
struct PMOutput
{
    std::vector<double> h;    /**< the final filter coefficients*/
    std::vector<double> x;    /**< the reference set used to
                                      generate the final filter (values are in \f$[-1,1]\f$
                                      and NOT \f$[0,\pi]\f$)*/
    std::size_t iter;               /**< number of iterations that were necessary to achieve
                                      convergence*/
    double delta;             /**< the final reference error */
    double Q;                 /**< convergence parameter value */
};

/*! An implementation of the uniform initialization approach for starting the Parks-McClellan
 * algorithm
 * @param[out] omega the initial set of references to be computed
 * @param[in] B the frequency bands of interest (i.e. stopbands and passbands for example)
 */

void initUniformExtremas(std::vector<double>& omega,
        std::vector<Band>& B);

/*! An implementation of the reference scaling approach mentioned in section 4 of the article.
 * @param[out] newX the reference set obtained from scaling the initial set x
 * @param[out] newChebyBands contains information about the bands of interest corresponding to the
 * new reference (i.e. how many reference points are inside each band). The bands are given inside
 * \f$[-1,1]\f$ (i.e. the CHEBY band space)
 * @param[out] newFreqBands contains information about the bands of interest corresponding to the
 * new reference (i.e. how many reference points are inside each band). The bands are given inside
 * \f$[0,\pi]\f$ (i.e. the FREQ band space)
 * @param[in] newXSize the size of the new reference set
 * @param[in] x the input reference on which will be used to perform the scaling
 * @param[in] chebyBands band information for the filter to which the x reference corresponds to.
 * The bands are given inside \f$[-1,1]\f$ (i.e. the CHEBY band space)
 * @param[in] freqBands band information for the filter to which the x reference corresponds to.
 * The bands are given inside \f$[0,\pi]\f$ (i.e. the FREQ band space)
 */
void referenceScaling(std::vector<double>& newX, std::vector<Band>& newChebyBands,
        std::vector<Band>& newFreqBands, std::size_t newXSize,
        std::vector<double>& x, std::vector<Band>& chebyBands,
        std::vector<Band>& freqBands);

/*! An internal routine which implements the exchange algorithm for designing FIR filters
 * @param[in] x the initial reference set
 * @param[in] chebyBands band information for the filter to which the x reference corresponds to.
 * The bands are given inside \f$[-1,1]\f$ (i.e. the CHEBY band space)
 * @param[in] epsT convergence parameter threshold (i.e quantizes the number of significant digits
 * of the minimax error that are accurate at the end of the final iteration)
 * @param[in] Nmax the degree used by the CPR method on each subinterval
 * @return information pertaining to the polynomial computed at the last iteration. Since there
 * is no explicit filter type (I to IV) given as input, the h vector of the output will correspond
 * to the coefficients of the frequency response \f$H_d(\omega)=\sum_{k=0}^{n}h_k\cos(\omega k)\f$
 *
 * The following piece of code shows how to use this function to construct a lowpass filter. It also
 * shows how to write frequency band specifications and how to use the uniform initialization routine.
 * It is equivalent to the sample code from firpm. <b>Important:</b> in order to determine the final coefficients
 * of the transfer function of the filter, some preprocessing specific to the type of the filter (I to IV) has to
 * be done on the elements of <tt>output.h</tt> (see for example the source code of the <tt>firpm</tt> functions on
 * how this is done for each type of filter).
 * @see firpm
 * @code
 * // frequency band specification
 * std::vector<Band> freqBands(2);
 * double pi = M_PI;
 *
 * freqBands[0].start = 0;
 * freqBands[0].stop = pi * 0.4;
 * freqBands[0].weight = [] (BandSpace, double) -> double {return double(1); };
 * freqBands[0].space = BandSpace::FREQ;
 * freqBands[0].amplitude = [](BandSpace, double) -> double { return double(1); };
 *
 * freqBands[1].start = pi * 0.5;
 * freqBands[1].stop = pi;
 * freqBands[1].weight = [] (BandSpace, double) -> double {return double(10); };
 * freqBands[1].space = BandSpace::FREQ;
 * freqBands[1].amplitude = [](BandSpace, double) -> double { return double(0); };
 * std::size_t degree = 100;  // filter degree
 *
 * // reference initialization code
 * std::vector<Band> chebyBands;
 * std::vector<double> omega(degree + 2u);
 * std::vector<double> x(degree + 2u);
 * initUniformExtremas(omega, freqBands, prec);
 * // apply the change of variable y = cos(x) so that we are working inside [-1, 1]
 * applyCos(x, omega);
 * bandConversion(chebyBands, freqBands, ConversionDirection::FROMFREQ);
 * // apply the exchange algorithm
 * PMOutput output = exchange(x, chebyBands);
 * @endcode
 */

PMOutput exchange(std::vector<double>& x,
        std::vector<Band>& chebyBands,
        double epsT = 0.01,
        int Nmax = 4);

/*! Parks-McClellan routine for implementing type I and II FIR filters. This routine uses uniform
 * initialization.
 * @param[in] N \f$N+1\f$ denotes the number of coefficients of the final transfer function. For even n, the
 * filter will be type I, while for odd n the type is II.
 * @param[in] f vector denoting the frequency ranges of each band of interest
 * @param[in] a the ideal amplitude at each point of f
 * @param[in] w the wight function value on each band
 * @param[in] epsT convergence parameter threshold (i.e quantizes the number of significant digits
 * of the minimax error that are accurate at the end of the final iteration)
 * @param[in] Nmax the degree used by the CPR method on each subinterval
 * @return information pertaining to the polynomial computed at the last iteration. The h vector of the
 * output contains the coefficients corresponding to the transfer function of the final filter
 * (in this case, for types I and II, the values are symmetrical to the middle coefficient(s))
 *
 * An example of how to use this function is given below. It designs a degree \f$n=100\f$ type I lowpass filter, with
 * passband \f$[0, 0.4\pi]\f$ and stopband \f$[0.5\pi, \pi]\f$. It has unit weight inside the passband and weight 10
 * inside the stopband. More examples, including code on how to use the reference scaling versions <tt>firpmRS</tt>,
 * are provided inside the test files.
 * @code
 * PMOutput output = firpm(200, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
 * @endcode
 */

PMOutput firpm(std::size_t N,
        std::vector<double>const& f,
        std::vector<double>const& a,
        std::vector<double>const& w,
        double epsT = 0.01,
        int Nmax = 4);


/*! Parks-McClellan routine for implementing type I and II FIR filters. This routine uses reference scaling by default.
 * @param[in] N \f$N+1\f$ denotes the number of coefficients of the final transfer function. For even n, the
 * filter will be type I, while for odd n the type is II.
 * @param[in] f vector denoting the frequency ranges of each band of interest
 * @param[in] a the ideal amplitude at each point of f
 * @param[in] w the wight function value on each band
 * @param[in] epsT convergence parameter threshold (i.e quantizes the number of significant digits
 * of the minimax error that are accurate at the end of the final iteration)
 * @param[in] depth how many times should reference scaling be applied recursively (default value is 1)
 * @param[in] Nmax the degree used by the CPR method on each subinterval
 * @param[in] root  what initialization strategy to use at the lowest level (uniform or AFP-based)
 * @return information pertaining to the polynomial computed at the last iteration. The h vector of the
 * output contains the coefficients corresponding to the transfer function of the final filter
 * (in this case, for types I and II, the values are symmetrical to the middle coefficient(s))*/

PMOutput firpmRS(std::size_t N,
        std::vector<double>const& f,
        std::vector<double>const& a,
        std::vector<double>const& w,
        double epsT = 0.01,
        std::size_t depth = 1u,
        int Nmax = 4,
        RootSolver root = RootSolver::UNIFORM);

/*! Parks-McClellan routine for implementing type I and II FIR filters. This routine uses AFP-based
 * initialization.
 * @param[in] N \f$N+1\f$ denotes the number of coefficients of the final transfer function. For even n, the
 * filter will be type I, while for odd n the type is II.
 * @param[in] f vector denoting the frequency ranges of each band of interest
 * @param[in] a the ideal amplitude at each point of f
 * @param[in] w the wight function value on each band
 * @param[in] epsT convergence parameter threshold (i.e quantizes the number of significant digits
 * of the minimax error that are accurate at the end of the final iteration)
 * @param[in] Nmax the degree used by the CPR method on each subinterval
 * @return information pertaining to the polynomial computed at the last iteration. The h vector of the
 * output contains the coefficients corresponding to the transfer function of the final filter
 * (in this case, for types I and II, the values are symmetrical to the middle coefficient(s))
 */

PMOutput firpmAFP(std::size_t N,
        std::vector<double>const& f,
        std::vector<double>const& a,
        std::vector<double>const& w,
        double epsT = 0.01,
        int Nmax = 4);


/*! Parks-McClellan routine for implementing type III and IV FIR filters. This routine uses uniform
 * initialization.
 * @param[in] N \f$N+1\f$ denotes the number of coefficients of the final transfer function. For even n, the
 * filter will be type III, while for odd n the type is IV.
 * @param[in] f vector denoting the frequency ranges of each band of interest
 * @param[in] a the ideal amplitude at each point of f
 * @param[in] w the wight function value on each band
 * @param[in] type denotes the type of filter we want to design: digital differentiator of Hilbert transformer
 * @param[in] epsT convergence parameter threshold (i.e quantizes the number of significant digits
 * of the minimax error that are accurate at the end of the final iteration)
 * @param[in] Nmax the degree used by the CPR method on each subinterval
 * @return information pertaining to the polynomial computed at the last iteration. The h vector of the
 * output contains the coefficients corresponding to the transfer function of the final filter
 * (in this case, for types III and IV, the values are antisymmetrical to the middle coefficient(s))*/

PMOutput firpm(std::size_t N,
        std::vector<double>const& f,
        std::vector<double>const& a,
        std::vector<double>const& w,
        ftype type,
        double epsT = 0.01,
        int Nmax = 4);

/*! Parks-McClellan routine for implementing type III and IV FIR filters. This routine uses reference scaling.
 * @param[in] N \f$N+1\f$ denotes the number of coefficients of the final transfer function. For even n, the
 * filter will be type III, while for odd n the type is IV.
 * @param[in] f vector denoting the frequency ranges of each band of interest
 * @param[in] a the ideal amplitude at each point of f
 * @param[in] w the wight function value on each band
 * @param[in] type denotes the type of filter we want to design: digital differentiator of Hilbert transformer
 * @param[in] epsT convergence parameter threshold (i.e quantizes the number of significant digits
 * of the minimax error that are accurate at the end of the final iteration)
 * @param[in] depth how many times should reference scaling be applied recursively (default value is 1)
 * @param[in] root  what initialization strategy to use at the lowest level (uniform or AFP-based)
 * @param[in] Nmax the degree used by the CPR method on each subinterval
 * @return information pertaining to the polynomial computed at the last iteration. The h vector of the
 * output contains the coefficients corresponding to the transfer function of the final filter
 * (in this case, for types III and IV, the values are antisymmetrical to the middle coefficient(s))*/
PMOutput firpmRS(std::size_t N,
        std::vector<double>const& f,
        std::vector<double>const& a,
        std::vector<double>const& w,
        ftype type,
        double epsT = 0.01,
        std::size_t depth = 1u,
        int Nmax = 4,
        RootSolver root = RootSolver::UNIFORM
        );

/*! Parks-McClellan routine for implementing type III and IV FIR filters. This routine uses AFP-based
 * initialization.
 * @param[in] N \f$N+1\f$ denotes the number of coefficients of the final transfer function. For even n, the
 * filter will be type III, while for odd n the type is IV.
 * @param[in] f vector denoting the frequency ranges of each band of interest
 * @param[in] a the ideal amplitude at each point of f
 * @param[in] w the wight function value on each band
 * @param[in] type denotes the type of filter we want to design: digital differentiator of Hilbert transformer
 * @param[in] epsT convergence parameter threshold (i.e quantizes the number of significant digits
 * of the minimax error that are accurate at the end of the final iteration)
 * @param[in] Nmax the degree used by the CPR method on each subinterval
 * @return information pertaining to the polynomial computed at the last iteration. The h vector of the
 * output contains the coefficients corresponding to the transfer function of the final filter
 * (in this case, for types III and IV, the values are antisymmetrical to the middle coefficient(s))*/

PMOutput firpmAFP(std::size_t N,
        std::vector<double>const& f,
        std::vector<double>const& a,
        std::vector<double>const& w,
        ftype type,
        double epsT = 0.01,
        int Nmax = 4);



#endif
