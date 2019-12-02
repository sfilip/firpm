//    firpm_d
//    Copyright (C) 2015-2019  S. Filip
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

#include "firpm/pm.h"
#include "firpm/band.h"
#include "firpm/barycentric.h"
#include "firpm/pmmath.h"
#include <set>
#include <fstream>
#include <sstream>

template<typename T>
using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using VectorXd = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T>
void chebvand(MatrixXd<T>& A, std::size_t degree,
        std::vector<T>& meshPoints,
        std::function<T(T)>& weightFunction)
{

    A.resize(degree + 1u, meshPoints.size());
    for(std::size_t i{0u}; i < meshPoints.size(); ++i)
    {
        T pointWeight = weightFunction(meshPoints[i]);
        A(0u, i) = 1;
        A(1u, i) = meshPoints[i];
        for(std::size_t j{2u}; j <= degree; ++j)
            A(j, i) = meshPoints[i] * A(j - 1u, i) * 2 - A(j - 2u, i);
        for(std::size_t j{0u}; j <= degree; ++j)
            A(j, i) *= pointWeight;
    }
}

// approximate Fekete points
template<typename T>
void afp(std::vector<T>& points, MatrixXd<T>& A,
         std::vector<T>& mesh)
{
    VectorXd<T> b = VectorXd<T>::Ones(A.rows());
    b(0) = 2;
    VectorXd<T> y = A.colPivHouseholderQr().solve(b);
    points.clear();

    for(Eigen::Index i{0}; i < y.rows(); ++i)
        if(y(i) != 0.0)
            points.push_back(mesh[i]);
    std::sort(points.begin(), points.end(),
            [](const T& lhs,
               const T& rhs) {
                return lhs < rhs;
            });
}

template<typename T>
void countBand(std::vector<band_t<T>>& cb, std::vector<T>& x)
{
    for(auto& it : cb)
        it.xs = 0u;
    std::size_t bandIt = 0u;
    for(std::size_t i{0u}; i < x.size(); ++i)
    {
        while(bandIt < cb.size() && cb[bandIt].stop < x[i])
            bandIt++;
        ++cb[bandIt].xs;
    }
}

template<typename T>
void wam(std::vector<T>& wam, std::vector<band_t<T>>& cb,
         std::size_t deg)
{
    std::vector<T> cp;
    equipts(cp, deg + 2u);
    cos(cp, cp);
    std::sort(begin(cp), end(cp));
    for(std::size_t i{0u}; i < cb.size(); ++i)
    {
        if(cb[i].start != cb[i].stop)
        {
            std::vector<T> bufferNodes;
            chgvar(bufferNodes, cp, cb[i].start, cb[i].stop);
            bufferNodes[0] = cb[i].start;
            bufferNodes[bufferNodes.size()-1u] = cb[i].stop;
            for(auto& it : bufferNodes)
                wam.push_back(it);
        }
        else
            wam.push_back(cb[i].start);
    }
}

template<typename T>
void uniform(std::vector<T>& omega,
        std::vector<band_t<T>>& B, std::size_t n)
{
    T avgDist = 0;
    omega.resize(n);

    std::vector<T> bandwidths(B.size());
    std::vector<std::size_t> nonPointBands;
    for(std::size_t i{0u}; i < B.size(); ++i) {
        bandwidths[i] = B[i].stop - B[i].start;
        if(bandwidths[i] > 0.0)
        {
            nonPointBands.push_back(i);
            avgDist += bandwidths[i];
        }
        B[i].xs = 1u;
    }
    if(nonPointBands.empty())
    {
        std::cerr << "ERROR: All intervals are points!\n";
        exit(EXIT_FAILURE);
    }

    avgDist /= (omega.size() - B.size());
    std::size_t npSize = nonPointBands.size();

    B[nonPointBands[npSize - 1u]].xs = omega.size() - (B.size() - npSize);
    T buffer;
    buffer = bandwidths[nonPointBands[0]] / avgDist;

    if (npSize > 1) {
        B[nonPointBands[0]].xs = pmmath::round(buffer) + 1;
        B[nonPointBands[npSize - 1u]].xs -= B[nonPointBands[0]].xs;
    }

    for(std::size_t i{1u}; i < npSize - 1u; ++i) {
        buffer = bandwidths[nonPointBands[i]] / avgDist;
        B[nonPointBands[i]].xs = pmmath::round(buffer) + 1;
        B[nonPointBands[npSize - 1u]].xs -= B[nonPointBands[i]].xs;
    }


    std::size_t startIndex = 0ul;
    for(std::size_t i{0ul}; i < B.size(); ++i) {
        if(B[i].xs > 1u)
            buffer = bandwidths[i] / (B[i].xs - 1);
        omega[startIndex] = B[i].start;
        omega[startIndex + B[i].xs - 1] = B[i].stop;
        for(std::size_t j{1ul}; j < B[i].xs - 1; ++j)
            omega[startIndex + j] = omega[startIndex + j - 1] + buffer;
        startIndex += B[i].xs;
    }
}

template<typename T>
void referenceScaling(std::vector<T>& newX, std::vector<band_t<T>>& newChebyBands,
        std::vector<band_t<T>>& newFreqBands, std::size_t newXSize,
        std::vector<T>& x, std::vector<band_t<T>>& chebyBands,
        std::vector<band_t<T>>& freqBands)
{
    std::vector<std::size_t> newDistribution(chebyBands.size());
    for(std::size_t i{0u}; i < chebyBands.size(); ++i)
        newDistribution[i] = 0u;
    std::size_t multipointBands = 0u;
    std::size_t offset = 0u;
    int twoInt = 0;
    for(std::size_t i{0u}; i < chebyBands.size(); ++i)
    {
        newX.push_back(x[offset]);
        if(chebyBands[i].xs > 2u)
        {
            ++multipointBands;
            for(std::size_t j{1u}; j < chebyBands[i].xs - 2u; ++j)
            {
                newX.push_back((x[offset + j] + x[offset + j + 1]) / 2);
                newX.push_back(x[offset + j]);
            }
            newX.push_back(x[offset + chebyBands[i].xs - 2u]);
            newX.push_back(x[offset + chebyBands[i].xs - 1u]);
            twoInt += 2;
        }
        else if(chebyBands[i].xs == 2u)
        {
            ++multipointBands;
            ++twoInt;
            newX.push_back(x[offset + 1u]);
        }
        offset += chebyBands[i].xs;
    }
    int threeInt = newXSize - newX.size() - twoInt;
    offset = 0u;
    for(std::size_t i{0u}; i < chebyBands.size(); ++i)
    {
            if(chebyBands[i].xs > 1u)
            {
                if(threeInt > 0)
                {
                    newX.push_back(x[offset] + (x[offset + 1] - x[offset]) / 3);
                    T secondValue = x[offset] + (x[offset + 1] - x[offset]) / 3
                        + (x[offset + 1] - x[offset]) / 3;
                    newX.push_back(secondValue);
                    threeInt--;
                    twoInt--;
                }
                else if (twoInt > 0)
                {
                    newX.push_back((x[offset] + x[offset + 1]) / 2);
                    twoInt--;
                }
            }
        offset += chebyBands[i].xs;
    }
    offset = 0;
    for(std::size_t i{0u}; i < chebyBands.size(); ++i)
    {
            if(chebyBands[i].xs > 2u)
            {
                if(threeInt > 0)
                {
                    newX.push_back(x[offset + chebyBands[i].xs - 2u] +
                            (x[offset + chebyBands[i].xs - 1u] -
                                x[offset + chebyBands[i].xs - 2u]) / 3);
                    T secondValue = x[offset + chebyBands[i].xs - 2u] +
                        (x[offset + chebyBands[i].xs - 1u] -
                            x[offset + chebyBands[i].xs - 2u]) / 3 +
                        (x[offset + chebyBands[i].xs - 1u] -
                            x[offset + chebyBands[i].xs - 2u]) / 3;
                    newX.push_back(secondValue);
                    threeInt--;
                    twoInt--;
                }
                else if (twoInt > 0)
                {
                    newX.push_back((x[offset + chebyBands[i].xs - 2u] +
                                x[offset + chebyBands[i].xs - 1u]) / 2);
                    twoInt--;
                }
            }
        offset += chebyBands[i].xs;
    }
    if(newXSize > newX.size())
    {
        std::cerr << "ERROR: Failed to do reference scaling\n";
        exit(EXIT_FAILURE);
    }
    newX.resize(newXSize);
    std::sort(newX.begin(), newX.end());
    std::size_t total = 0u;
    for(std::size_t i{0u}; i < newX.size(); ++i)
    {
            for(std::size_t j{0u}; j < chebyBands.size(); ++j)
                if(newX[i] >= chebyBands[j].start && newX[i] <= chebyBands[j].stop)
                {
                    newDistribution[j]++;
                    ++total;
                }
    }
    if(total != newXSize)
    {
        std::cout << "ERROR: Failed to find reference scaling distribution!\n";
        exit(EXIT_FAILURE);
    }


    for (std::size_t i{0u}; i < chebyBands.size(); ++i)
    {
        newFreqBands[freqBands.size() - 1u - i].xs = newDistribution[i];
        newChebyBands[i].xs = newDistribution[i];
    }
}

template<typename T>
void split(std::vector<std::pair<T, T>>& subIntervals,
        std::vector<band_t<T>>& chebyBands,
        std::vector<T> &x)
{
    std::size_t bandOffset{0u};
    for(std::size_t i{0u}; i < chebyBands.size(); ++i)
    {
        T middleValA{0.0}, middleValB{0.0};
        if(chebyBands[i].xs == 0u) {
            subIntervals.push_back(
                std::make_pair(chebyBands[i].start, 
                               chebyBands[i].stop));
        } else {
            if (bandOffset < x.size() 
                && x[bandOffset] > chebyBands[i].start
                && x[bandOffset] < chebyBands[i].stop)
            {
                middleValA = x[bandOffset];
                subIntervals.push_back(
                    std::make_pair(chebyBands[i].start, middleValA));
                if (chebyBands[i].xs == 1u) {
                    subIntervals.push_back(
                        std::make_pair(middleValA, 
                                       chebyBands[i].stop));
                }
            } else {
                middleValA = chebyBands[i].start;
                if(chebyBands[i].xs == 1u) {
                    subIntervals.push_back(std::make_pair(
                        middleValA, chebyBands[i].stop
                    ));
                }
            }
        }
        if(bandOffset < x.size() 
           && chebyBands[i].xs > 1)
        {
            for(std::size_t j{bandOffset};
                j < bandOffset + chebyBands[i].xs - 1u; ++j)
            {
                middleValB = x[j + 1];
                subIntervals.push_back(std::make_pair(middleValA, middleValB));
                middleValA = middleValB;
            }
            if(middleValA != chebyBands[i].stop)
                subIntervals.push_back(
                    std::make_pair(middleValA, chebyBands[i].stop));
        }
        bandOffset += chebyBands[i].xs;
    }
}

template<typename T>
void extremaSearch(T& convergenceOrder,
        T& delta, std::vector<T>& eigenExtrema,
        std::vector<T>& x, std::vector<band_t<T>>& chebyBands,
        std::size_t Nmax, unsigned long prec)
{
#ifdef HAVE_MPFR
    mpfr_prec_t prevPrec = mpfr::mpreal::get_default_prec();
    mpfr::mpreal::set_default_prec(prec);
#endif

    // 1.   Split the initial [-1, 1] interval in subintervals
    //      in order that we can use a reasonable size matrix
    //      eigenvalue solver on the subintervals
    std::vector<std::pair<T, T>> subIntervals;
    std::pair<T, T> dom{std::make_pair(-1.0, 1.0)};

    split(subIntervals, chebyBands, x);
    //std::cout << "Number of subintervals: "
    //    << subIntervals.size() << std::endl;

    // 2.   Compute the barycentric variables (i.e. weights)
    //      needed for the current iteration

    std::vector<T> w(x.size());
    baryweights(w, x);


    compdelta(delta, w, x, chebyBands);
    //std::cout << "delta = " << delta << std::endl;

    std::vector<T> C(x.size());
    compc(C, delta, x, chebyBands);

    // 3.   Use an eigenvalue solver on each subinterval to find the
    //      local extrema that are located inside the frequency bands
    std::vector<T> chebyNodes(Nmax + 1u);
    equipts(chebyNodes, Nmax + 1u);
    cos(chebyNodes, chebyNodes);

    std::vector<std::pair<T, T>> potentialExtrema;
    std::vector<T> pEx;
    T extremaErrorValueLeft;
    T extremaErrorValueRight;
    T extremaErrorValue;
    comperror(extremaErrorValue, chebyBands[0].start,
            delta, x, C, w, chebyBands);
    potentialExtrema.push_back(std::make_pair(
            chebyBands[0].start, extremaErrorValue));


    for (std::size_t i{0u}; i < chebyBands.size() - 1u; ++i)
    {
        comperror(extremaErrorValueLeft, chebyBands[i].stop,
                delta, x, C, w, chebyBands);
        comperror(extremaErrorValueRight, chebyBands[i + 1].start,
                delta, x, C, w, chebyBands);
        bool sgnLeft = pmmath::signbit(extremaErrorValueLeft);
        bool sgnRight = pmmath::signbit(extremaErrorValueRight);
        if (sgnLeft != sgnRight) {
            potentialExtrema.push_back(std::make_pair(
                    chebyBands[i].stop, extremaErrorValueLeft));
            potentialExtrema.push_back(std::make_pair(
                    chebyBands[i + 1u].start, extremaErrorValueRight));
        } else {
            T abs1 = pmmath::fabs(extremaErrorValueLeft);
            T abs2 = pmmath::fabs(extremaErrorValueRight);
            if(abs1 > abs2)
                potentialExtrema.push_back(std::make_pair(
                        chebyBands[i].stop, extremaErrorValueLeft));
            else
                potentialExtrema.push_back(std::make_pair(
                        chebyBands[i + 1u].start, extremaErrorValueRight));
        }
    }
    comperror(extremaErrorValue,
            chebyBands[chebyBands.size() - 1u].stop,
            delta, x, C, w, chebyBands);
    potentialExtrema.push_back(std::make_pair(
            chebyBands[chebyBands.size() - 1u].stop,
            extremaErrorValue));


    std::vector<std::vector<T>> pExs(subIntervals.size());

    #pragma omp parallel for
    for (std::size_t i = 0u; i < subIntervals.size(); ++i)
    {
        #ifdef HAVE_MPFR
            mpfr_prec_t prevPrec = mpfr::mpreal::get_default_prec();
            mpfr::mpreal::set_default_prec(prec);
        #endif
        // find the Chebyshev nodes scaled to the current subinterval  
        std::vector<T> siCN(Nmax + 1u);
        chgvar(siCN, chebyNodes, subIntervals[i].first,
                subIntervals[i].second);

        // compute the Chebyshev interpolation function values on the
        // current subinterval
        std::vector<T> fx(Nmax + 1u);
        for (std::size_t j{0u}; j < fx.size(); ++j)
            comperror(fx[j], siCN[j], delta, x, C, w,
                    chebyBands);

        // compute the values of the CI coefficients and those of its
        // derivative
        std::vector<T> c(Nmax + 1u);
        chebcoeffs(c, fx);
        std::vector<T> dc(Nmax);
        diffcoeffs(dc, c);

        // solve the corresponding eigenvalue problem and determine the
        // local extrema situated in the current subinterval
        std::vector<T> eigenRoots;
        roots(eigenRoots, dc, dom);
        if(!eigenRoots.empty()) {
            chgvar(eigenRoots, eigenRoots,
                    subIntervals[i].first, subIntervals[i].second);
            for (std::size_t j{0u}; j < eigenRoots.size(); ++j)
                pExs[i].push_back(eigenRoots[j]);
        }
        pExs[i].push_back(subIntervals[i].first);
        pExs[i].push_back(subIntervals[i].second);
        #ifdef HAVE_MPFR
            mpfr::mpreal::set_default_prec(prevPrec);
        #endif
    }

    for(std::size_t i{0u}; i < pExs.size(); ++i)
        for(std::size_t j{0u}; j < pExs[i].size(); ++j)
            pEx.push_back(pExs[i][j]);

    std::size_t startingOffset = potentialExtrema.size();
    potentialExtrema.resize(potentialExtrema.size() + pEx.size());
    #pragma omp parallel for
    for(std::size_t i = 0u; i < pEx.size(); ++i)
    {
        #ifdef HAVE_MPFR
            mpfr_prec_t prevPrec = mpfr::mpreal::get_default_prec();
            mpfr::mpreal::set_default_prec(prec);
        #endif
        T valBuffer;
        comperror(valBuffer, pEx[i],
                delta, x, C, w, chebyBands);
        potentialExtrema[startingOffset + i] = std::make_pair(pEx[i], valBuffer);
        #ifdef HAVE_MPFR
            mpfr::mpreal::set_default_prec(prevPrec);
        #endif
    }

    // sort list of potential extrema in increasing order
    std::sort(potentialExtrema.begin(), potentialExtrema.end(),
            [](const std::pair<T, T>& lhs,
               const std::pair<T, T>& rhs) {
                return lhs.first < rhs.first;
            });

    eigenExtrema.clear();
    std::size_t extremaIt{0u};
    std::vector<std::pair<T, T>> alternatingExtrema;
    T minError = INT_MAX;
    T maxError = INT_MIN;
    T absError;

    while (extremaIt < potentialExtrema.size())
    {
        std::pair<T, T> maxErrorPoint;
        maxErrorPoint = potentialExtrema[extremaIt];
        while(extremaIt < potentialExtrema.size() - 1u &&
            (pmmath::signbit(maxErrorPoint.second) ==
             pmmath::signbit(potentialExtrema[extremaIt + 1u].second)))
        {
            ++extremaIt;
            if (pmmath::fabs(maxErrorPoint.second) < pmmath::fabs(potentialExtrema[extremaIt].second))
                maxErrorPoint = potentialExtrema[extremaIt];
        }
        if (pmmath::isfinite(maxErrorPoint.second)) {
            alternatingExtrema.push_back(maxErrorPoint);
        }
        ++extremaIt;
    }
    std::vector<std::pair<T, T>> bufferExtrema;
    //std::cout << "Alternating extrema: " << x.size() << " | "
    //    << alternatingExtrema.size() << std::endl;

    if(alternatingExtrema.size() < x.size())
    {
        std::cerr << "WARNING: The exchange algorithm did not converge.\n";
        std::cerr << "TRIGGER: Not enough alternating extrema!\n"
            << "POSSIBLE CAUSE: Nmax too small\n";
        convergenceOrder = 2.0;
        return;
    }
    else if (alternatingExtrema.size() > x.size())
    {
        std::size_t remSuperfluous = alternatingExtrema.size() - x.size();
        if (remSuperfluous % 2 != 0)
        {
            if(remSuperfluous == 1u)
            {
                std::vector<T> x1, x2;
                x1.push_back(alternatingExtrema[0u].first);
                for(std::size_t i{1u}; i < alternatingExtrema.size() - 1u; ++i)
                {
                    x1.push_back(alternatingExtrema[i].first);
                    x2.push_back(alternatingExtrema[i].first);
                }
                x2.push_back(alternatingExtrema[alternatingExtrema.size() - 1u].first);
                T delta1, delta2;
                compdelta(delta1, x1, chebyBands);
                compdelta(delta2, x2, chebyBands);
                delta1 = pmmath::fabs(delta1);
                delta2 = pmmath::fabs(delta2);
                std::size_t sIndex{1u};
                if(delta1 > delta2)
                    sIndex = 0u;
                for(std::size_t i{sIndex}; i < alternatingExtrema.size() + sIndex - 1u; ++i)
                    bufferExtrema.push_back(alternatingExtrema[i]);
                alternatingExtrema = bufferExtrema;
                bufferExtrema.clear();
            }
            else
            {
                T abs1 = pmmath::fabs(alternatingExtrema[0u].second);
                T abs2 = pmmath::fabs(alternatingExtrema[alternatingExtrema.size() - 1u].second);
                std::size_t sIndex = 0u;
                if (abs1 < abs2)
                    sIndex = 1u;
                for(std::size_t i{sIndex}; i < alternatingExtrema.size() + sIndex - 1u; ++i)
                    bufferExtrema.push_back(alternatingExtrema[i]);
                alternatingExtrema = bufferExtrema;
                bufferExtrema.clear();
            }
        }


        while (alternatingExtrema.size() > x.size())
        {
            std::size_t toRemoveIndex{0u};
            T minValToRemove = pmmath::fmin(pmmath::fabs(alternatingExtrema[0u].second),
                                              pmmath::fabs(alternatingExtrema[1u].second));
            T removeBuffer;
            for (std::size_t i{1u}; i < alternatingExtrema.size() - 1u; ++i)
            {
                removeBuffer = pmmath::fmin(pmmath::fabs(alternatingExtrema[i].second),
                                   pmmath::fabs(alternatingExtrema[i + 1u].second));
                if (removeBuffer < minValToRemove)
                {
                    minValToRemove = removeBuffer;
                    toRemoveIndex  = i;
                }
            }
            for (std::size_t i{0u}; i < toRemoveIndex; ++i)
                bufferExtrema.push_back(alternatingExtrema[i]);
            for (std::size_t i{toRemoveIndex + 2u}; i < alternatingExtrema.size(); ++i)
                bufferExtrema.push_back(alternatingExtrema[i]);
            alternatingExtrema = bufferExtrema;
            bufferExtrema.clear();
        }


    }
    if (alternatingExtrema.size() < x.size())
    {
        std::cerr << "Trouble!\n";
        exit(EXIT_FAILURE);
    }

    //std::cout << "After removal: " << alternatingExtrema.size() << std::endl;
    for (auto& it : alternatingExtrema)
    {
        eigenExtrema.push_back(it.first);
        absError = pmmath::fabs(it.second);
        minError = pmmath::fmin(minError, absError);
        maxError = pmmath::fmax(maxError, absError);
    }

    //std::cout << "Min error = " << minError << std::endl;
    //std::cout << "Max error = " << maxError << std::endl;
    convergenceOrder = (maxError - minError) / maxError;
    //std::cout << "Convergence order = " << convergenceOrder << std::endl;
    // update the extrema count in each frequency band
    std::size_t bIndex{0u};
    for(std::size_t i{0u}; i < chebyBands.size(); ++i)
    {
        chebyBands[i].xs = 0u;
    }
    for(auto &it : eigenExtrema)
    {
        if(chebyBands[bIndex].start <= it && it <= chebyBands[bIndex].stop)
        {
            ++chebyBands[bIndex].xs;
        }
        else
        {
            ++bIndex;
            ++chebyBands[bIndex].xs;
        }
    }
#ifdef HAVE_MPFR
    mpfr::mpreal::set_default_prec(prevPrec);
#endif
}


// REMARK: remember that this routine assumes that the information
// pertaining to the reference x and the frequency bands (i.e. the
// number of reference values inside each band) is given at the
// beginning of the execution
template<typename T>
pmoutput_t<T> exchange(std::vector<T>& x,
        std::vector<band_t<T>>& chebyBands, double eps,
        std::size_t Nmax, unsigned long prec)
{
    pmoutput_t<T> output;

    std::size_t degree = x.size() - 2u;
    std::sort(x.begin(), x.end(),
            [](const T& lhs,
               const T& rhs) {
                return lhs < rhs;
            });
    std::vector<T> startX{x};
    std::cout.precision(20);

    output.q = 1;
    output.iter = 0u;
    do {
        ++output.iter;
        //std::cout << "*********ITERATION " << output.iter << " **********\n";
        extremaSearch(output.q, output.delta,
                output.x, startX, chebyBands, Nmax, prec);
        startX = output.x;
        if(output.q > 1.0)
            break;
    } while (output.q > eps && output.iter <= 100u);

    if(pmmath::isnan(output.delta) || pmmath::isnan(output.q))
        std::cerr << "WARNING: The exchange algorithm did not converge.\n"
            << "TRIGGER: numerical instability\n"
            << "POSSIBLE CAUSES: poor starting reference and/or "
            << "a too small value for Nmax.\n";

    if(output.iter >= 101u)
        std::cerr << "WARNING: The exchange algorithm did not converge.\n"
            << "TRIGGER: exceeded iteration threshold of 100\n"
            << "POSSIBLE CAUSES: poor starting reference and/or "
            << "a too small value for Nmax.\n";


    output.h.resize(degree + 1u);
    std::vector<T> finalC(output.x.size());
    std::vector<T> finalAlpha(output.x.size());
    baryweights(finalAlpha, output.x);
    T finalDelta = output.delta;
    output.delta = pmmath::fabs(output.delta);
    //std::cout << "MINIMAX delta = " << output.delta << std::endl;
    compc(finalC, finalDelta, output.x, chebyBands);
    std::vector<T> finalChebyNodes(degree + 1);
    equipts(finalChebyNodes, degree + 1);
    cos(finalChebyNodes, finalChebyNodes);
    std::vector<T> fv(degree + 1);

    for (std::size_t i{0u}; i < fv.size(); ++i)
        approx(fv[i], finalChebyNodes[i], output.x,
                finalC, finalAlpha);

    chebcoeffs(output.h, fv);

    return output;
}

template<typename T>
void parseSpecification(std::vector<T> const &f,
            std::vector<T> const &a,
            std::vector<T> const &w)
{
    if(f.size() != a.size()) {
        std::cerr << "ERROR: Frequency and amplitude vector sizes"
            << " do not match!\n";
        exit(EXIT_FAILURE);
    }

    if(f.size() % 2 != 0) {
        std::cerr << "ERROR: Frequency band edges must come in pairs!\n";
        exit(EXIT_FAILURE);
    }

    if(f.size() != w.size() * 2u) {
        std::cerr << "ERROR: Weight vector size does not match the"
            << " the number of frequency bands in the specification!\n";
        exit(EXIT_FAILURE);
    }

    for(std::size_t i{0u}; i < f.size() - 1u; ++i) {
        if(f[i] == f[i + 1u] && (a[i] != a[i + 1u])) {
            std::cerr << "ERROR: Adjacent bands with discontinuities"
                << " are not allowed!\n";
            exit(EXIT_FAILURE);
        }
        if(f[i] > f[i + 1u]) {
            std::cerr << "ERROR: Frequency vector entries must be in "
                << "nondecreasing order!\n";
            exit(EXIT_FAILURE);
        }
    }
    for(std::size_t i{0u}; i < w.size(); ++i) {
        if(w[i] <= 0.0) {
            std::cerr << "ERROR: Band weights must be positive!\n";
            exit(EXIT_FAILURE);
        }
    }

    if(f[0u] < 0.0 || f[f.size() - 1u] > 1.0) {
        std::cerr << "ERROR: Normalized frequency band edges must be "
            << "between 0 and 1!\n";
        exit(EXIT_FAILURE);
    }
}

template<typename T>
pmoutput_t<T> firpm(std::size_t n,
            std::vector<T>const &f,
            std::vector<T>const &a,
            std::vector<T>const &w,
            double eps,
            std::size_t nmax,
            init_t strategy,
            std::size_t depth,
            init_t rstrategy, 
            unsigned long prec)
{
    parseSpecification(f, a, w);
    std::vector<T> h;
    std::vector<band_t<T>> fbands(w.size());
    std::vector<band_t<T>> cbands;
    if(n % 2 != 0) {
        if(f[f.size()-1u] == 1 && a[a.size()-1u] != 0) {
            std::cout << "WARNING: gain at Nyquist frequency different from 0.\n"
                << "Increasing the number of taps by one and passing to a "
                << "type I filter" << std::endl;
            ++n;
        }
    }
    std::size_t deg = n / 2;
    if(n % 2 == 0) {            // type I filter
        for(std::size_t i{0u}; i < fbands.size(); ++i) {
            fbands[i].start = pmmath::const_pi<T>() * f[2u*i];
            fbands[i].stop  = pmmath::const_pi<T>() * f[2u*i+1u];
            fbands[i].space = space_t::FREQ;
            fbands[i].amplitude = [i, &a, &fbands](space_t space, T x) -> T {
                if(a[2u*i] != a[2u*i+1u]) {
                    if(space == space_t::CHEBY)
                        x = pmmath::acos(x);
                    return ((x-fbands[i].start) * a[2u*i+1u] - 
                            (x-fbands[i].stop) * a[2u*i]) /
                            (fbands[i].stop - fbands[i].start);
                }
                return a[2u*i];
            };
            fbands[i].weight = [i, &w](space_t, T x) -> T {
                return w[i];
            };
        }
    } else {                    // type II filter
        for(std::size_t i{0u}; i < fbands.size(); ++i) {
            fbands[i].start = pmmath::const_pi<T>() * f[2u*i];
            if(f[2u*i + 1u] == 1.0) {
                if(f[2u*i] < 0.9999)
                    fbands[i].stop = pmmath::const_pi<T>() * T(0.9999);
                else
                    fbands[i].stop = pmmath::const_pi<T>() * ((f[2u*i] + 1) / 2);
            } 
            else 
                fbands[i].stop = pmmath::const_pi<T>() * f[2u*i+1u];
            fbands[i].space = space_t::FREQ;
            fbands[i].amplitude = [i, &a, &fbands](space_t space, T x) -> T {
                if(a[2u*i] != a[2u*i+1u]) {
                    if(space == space_t::CHEBY)
                        x = pmmath::acos(x);
                    return (((x-fbands[i].start) * a[2u*i+1u] -
                            (x-fbands[i].stop) * a[2u*i]) /
                            (fbands[i].stop - fbands[i].start)) / pmmath::cos(x/2);
                }
                if(space == space_t::FREQ)
                    return a[2u*i] / pmmath::cos(x/2);
                else
                    return a[2u*i] / pmmath::sqrt((x+1)/2);
            };
            fbands[i].weight = [i, &w](space_t space, T x) -> T {
                if(space == space_t::FREQ)
                    return pmmath::cos(x/2) * w[i];
                else
                    return pmmath::sqrt((x+1)/2) * w[i];
            };
        }
    }

    pmoutput_t<T> output;
    std::vector<T> x;
    bandconv(cbands, fbands, convdir_t::FROMFREQ);
    std::function<T(T)> wf = [&cbands](T x) -> T {
        for(std::size_t i{0u}; i < cbands.size(); ++i)
            if(cbands[i].start <= x && x <= cbands[i].stop)
                return cbands[i].weight(space_t::CHEBY, x);
        // this should never execute
        return 1.0;
    };

    if(fbands.size() > (deg + 2u) / 4)
        strategy = init_t::AFP;

    switch(strategy) {
        case init_t::UNIFORM: 
        {
            if (fbands.size() <= (deg + 2u) / 4) { 
                std::vector<T> omega;
                uniform(omega, fbands, deg + 2u);
                cos(x, omega);
                bandconv(cbands, fbands, convdir_t::FROMFREQ);
            } else {
                // use AFP strategy for very small degrees (wrt nb of bands)
                std::vector<T> mesh;
                wam(mesh, cbands, deg);
                MatrixXd<T> A;
                chebvand(A, deg+1u, mesh, wf);
                afp(x, A, mesh);
                if(x.size() != deg + 2u) {
                    std::cerr << "ERROR: AFP strategy failed to produce a valid starting reference\n"
                        << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix\n";
                    exit(EXIT_FAILURE);
                }
                countBand(cbands, x);
            }
            output = exchange(x, cbands, eps, nmax, prec);
        } break;
        case init_t::SCALING: 
        {
            std::vector<std::size_t> sdegs(depth+1u);
            sdegs[depth] = deg;
            for(int i{(int)depth-1}; i >= 0; --i)
                sdegs[i] = sdegs[i+1]/2;

            if(rstrategy == init_t::UNIFORM) {
                if (fbands.size() <= (sdegs[0] + 2u) / 4) {
                    std::vector<T> omega;
                    uniform(omega, fbands, sdegs[0]+2u);
                    cos(x, omega);
                    bandconv(cbands, fbands, convdir_t::FROMFREQ);
                } else {
                    // use AFP strategy for very small degrees (wrt nb of bands)
                    std::vector<T> mesh;
                    wam(mesh, cbands, sdegs[0]);
                    MatrixXd<T> A;
                    chebvand(A, sdegs[0]+1u, mesh, wf);
                    afp(x, A, mesh);
                    if(x.size() != sdegs[0] + 2u) {
                        std::cerr << "ERROR: AFP strategy failed to produce a valid starting reference\n"
                            << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix\n";
                        exit(EXIT_FAILURE);
                    }
                    countBand(cbands, x);                    
                }
                output = exchange(x, cbands, eps, nmax, prec);
            } else { // AFP-based strategy
                std::vector<T> mesh;
                wam(mesh, cbands, sdegs[0]);
                MatrixXd<T> A;
                chebvand(A, sdegs[0]+1u, mesh, wf);
                afp(x, A, mesh);
                if(x.size() != sdegs[0] + 2u) {
                    std::cerr << "ERROR: AFP strategy failed to produce a valid starting reference\n"
                        << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix\n";
                    exit(EXIT_FAILURE);
                }
                countBand(cbands, x);
                output = exchange(x, cbands, eps, nmax, prec);
            }
            for(std::size_t i{1u}; i <= depth && output.q <= 0.5; ++i) {
                x.clear();
                referenceScaling(x, cbands, fbands, sdegs[i]+2u, 
                                 output.x, cbands, fbands);
                output = exchange(x, cbands, eps, nmax, prec);
            }
        } break;
        default: { // AFP-based initialization
            std::vector<T> mesh;
            wam(mesh, cbands, deg);
            MatrixXd<T> A;
            chebvand(A, deg+1u, mesh, wf);
            afp(x, A, mesh);
            if(x.size() != deg + 2u) {
                std::cerr << "ERROR: AFP strategy failed to produce a valid starting reference\n"
                    << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix\n";
                exit(EXIT_FAILURE);
            }
            countBand(cbands, x);
            output = exchange(x, cbands, eps, nmax, prec);
        }
    }

    h.resize(n+1u);
    if(n % 2 == 0) {
        h[deg] = output.h[0];
        for(std::size_t i{0u}; i < deg; ++i)
            h[i] = h[n-i] = output.h[deg-i] / 2u;
    } else {
        h[0] = h[n] = output.h[deg] / 4u;
        h[deg] = h[deg+1u] = (output.h[0] * 2 + output.h[1]) / 4u;
        for(std::size_t i{2u}; i < deg + 1u; ++i)
            h[deg+1u-i] = h[deg+i] = (output.h[i-1u] + output.h[i]) / 4u;
    }
    output.h = h;

    return output;
}

template<typename T>
pmoutput_t<T> firpmRS(std::size_t n,
            std::vector<T>const &f,
            std::vector<T>const &a,
            std::vector<T>const &w,
            double eps,
	        std::size_t nmax,
            std::size_t depth, 
            init_t rstrategy,
            unsigned long prec)
{
    return firpm<T>(n, f, a, w, eps, nmax, init_t::SCALING, depth, rstrategy, prec);
}

template<typename T>
pmoutput_t<T> firpmAFP(std::size_t n,
            std::vector<T>const &f,
            std::vector<T>const &a,
            std::vector<T>const &w,
            double eps, std::size_t nmax,
            unsigned long prec)
{
    return firpm<T>(n, f, a, w, eps, nmax, init_t::AFP, 0u, init_t::AFP, prec);
}

template<typename T>
pmoutput_t<T> firpm(std::size_t n,
            std::vector<T>const& f,
            std::vector<T>const& a,
            std::vector<T>const& w,
            filter_t type,
            double eps,
            std::size_t nmax,
            init_t strategy,
            std::size_t depth,
            init_t rstrategy,
            unsigned long prec)
{
    parseSpecification(f, a, w);
    std::vector<T> h;
    std::vector<band_t<T>> fbands(w.size());
    std::vector<band_t<T>> cbands;
    std::vector<T> fn{f};
    std::size_t deg = n / 2u;
    T sFactor = a[1] / (f[1] * pmmath::const_pi<T>());

    if(n % 2 == 0) { // TYPE III
        if(f[0u] == 0.0) {
            if(fn[1u] < 1e-5)
                fn[0u] = fn[1u] / 2;
            else
                fn[0u] = 1e-5;                    
        }
        if(f[f.size() - 1u] == 1.0) {
            if(f[f.size() - 2u] > 0.9999)
                fn[f.size() - 1u] = (f[f.size() - 2u] + 1.0) / 2;
            else
                fn[f.size() - 1u] = 0.9999;
        }
        --deg;
        fbands[0u].start = pmmath::const_pi<T>() * fn[0u];
        fbands[0u].stop  = pmmath::const_pi<T>() * fn[1u];
        fbands[0u].space = space_t::FREQ;
        if(type == filter_t::FIR_DIFFERENTIATOR) {
            fbands[0u].weight = [&w](space_t space, T x) -> T {
                if(space == space_t::FREQ)
                    return (pmmath::sin(x) / x) * w[0u];
                else
                    return (pmmath::sqrt(T(1.0) - x * x) / pmmath::acos(x)) * w[0u];
            };
            fbands[0u].amplitude = [sFactor](space_t space, T x) -> T {
                if(space == space_t::FREQ)
                    return (x / pmmath::sin(x)) * sFactor;
                else 
                    return (pmmath::acos(x) / pmmath::sqrt(T(1.0) - x * x)) * sFactor;
            };
        } else { // FIR_HILBERT
            fbands[0u].weight = [&w](space_t space, T x) -> T {
                if(space == space_t::FREQ)
                    return pmmath::sin(x) * w[0u];
                else
                    return pmmath::sqrt(T(1.0) - x * x) * w[0u];
            };
            fbands[0u].amplitude = [&a, &fbands](space_t space, T x) -> T {
                if(space == space_t::CHEBY)
                    x = pmmath::acos(x);
                if(a[0u] != a[1u]) 
                    return (((x - fbands[0].start) * a[1u] -
                            (x - fbands[0].stop) * a[0u]) /
                            (fbands[0u].stop - fbands[0u].start)) / pmmath::sin(x);
                return a[0u] / pmmath::sin(x);
            };
        }
        for(std::size_t i{1u}; i < fbands.size(); ++i) {
            fbands[i].start = pmmath::const_pi<T>() * fn[2u * i];
            fbands[i].stop  = pmmath::const_pi<T>() * fn[2u * i + 1u];
            fbands[i].space = space_t::FREQ;
            if(type == filter_t::FIR_DIFFERENTIATOR) {
                fbands[i].weight = [&w, i](space_t space, T x) -> T {
                    if(space == space_t::FREQ)
                        return pmmath::sin(x) * w[i];
                    else
                        return pmmath::sqrt(T(1.0) - x * x) * w[i];
                };
                fbands[i].amplitude = [&fbands, &a, i](space_t space, T x) -> T {
                    if(a[2u * i] != a[2u * i + 1u]) {
                        if(space == space_t::CHEBY) {
                            x = pmmath::acos(x);
                            return ((x - fbands[i].start) * a[2u * i + 1u] -
                                    (x - fbands[i].stop) * a[2u * i]) /
                                    (fbands[i].stop - fbands[i].start);
                        }
                    }
                    return a[2u * i];
                };
            } else { // FIR_HILBERT
                fbands[i].weight = [&w, i](space_t space, T x) -> T {
                    if(space == space_t::FREQ)
                        return pmmath::sin(x) * w[i];
                    else
                        return pmmath::sqrt(T(1.0) - x * x) * w[i];
                };
                fbands[i].amplitude = [&fbands, &a, i](space_t space, T x) -> T {
                    if(space == space_t::CHEBY)
                        x = pmmath::acos(x);
                    if(a[2u * i] != a[2u * i + 1u])
                        return (((x - fbands[i].start) * a[2u * i + 1u] -
                                (x - fbands[i].stop) * a[2u * i]) / 
                                (fbands[i].stop - fbands[i].start)) / pmmath::sin(x);
                    return a[2u * i] / pmmath::sin(x);
                };
            }
        }
    } else { // TYPE IV
        if(f[0u] == 0.0) {
            if(fn[1u] < 1e-5)
                fn[0u] = fn[1u] / 2;
            else
                fn[0u] = 1e-5;
        }

        fbands[0u].start = pmmath::const_pi<T>() * fn[0u];
        fbands[0u].stop  = pmmath::const_pi<T>() * fn[1u];
        fbands[0u].space = space_t::FREQ;
        if(type == filter_t::FIR_DIFFERENTIATOR) {
            fbands[0u].weight = [&w](space_t space, T x) -> T {
                if(space == space_t::FREQ)
                    return (pmmath::sin(x / 2) / x) * w[0u];
                else
                    return (pmmath::sin(pmmath::acos(x) / 2) / pmmath::acos(x)) * w[0u];
            };
            fbands[0u].amplitude = [sFactor](space_t space, T x) -> T {
                if(space == space_t::FREQ) 
                    return (x / pmmath::sin(x / 2)) * sFactor;
                else 
                    return (pmmath::acos(x) / pmmath::sin(pmmath::acos(x) / 2)) * sFactor;
            };
        } else { // FIR_HILBERT
            fbands[0u].weight = [&w](space_t space, T x) -> T {
                if(space == space_t::FREQ)
                    return pmmath::sin(x / 2) * w[0u];
                else
                    return pmmath::sin(pmmath::acos(x) / 2) * w[0u];
            };
            fbands[0u].amplitude = [&fbands, &a](space_t space, T x) -> T {
                if(space == space_t::CHEBY)
                    x = pmmath::acos(x);
                if(a[0u] != a[1u]) 
                    return (((x - fbands[0].start) * a[1u] -
                            (x - fbands[0].stop) * a[0u]) /
                            (fbands[0].stop - fbands[0].start)) / pmmath::sin(x / 2);
                return a[0u] / pmmath::sin(x / 2);
            };
        }
        for(std::size_t i{1u}; i < fbands.size(); ++i) {
            fbands[i].start = pmmath::const_pi<T>() * fn[2u * i];
            fbands[i].stop  = pmmath::const_pi<T>() * fn[2u * i + 1u];
            fbands[i].space = space_t::FREQ;
            if(type == filter_t::FIR_DIFFERENTIATOR) {
                fbands[i].weight = [&w, i](space_t space, T x) -> T {
                    if(space == space_t::FREQ)
                        return pmmath::sin(x / 2) * w[i];
                    else
                        return (pmmath::sin(pmmath::acos(x) / 2)) * w[i];
                };
                fbands[i].amplitude = [&fbands, &a, i](space_t space, T x) -> T {
                    if(a[2u * i] != a[2u * i + 1u]) {
                        if(space == space_t::CHEBY)
                            x = pmmath::acos(x);
                        return ((x - fbands[i].start) * a[2u * i + 1u] -
                                (x - fbands[i].stop) * a[2u * i]) /
                                (fbands[i].stop - fbands[i].start);
                    }
                    return a[2u * i];
                };
            } else { // FIR_HILBERT
                fbands[i].weight = [&w, i](space_t space, T x) -> T {
                    if(space == space_t::FREQ)
                        return pmmath::sin(x / 2) * w[i];
                    else
                        return pmmath::sin(pmmath::acos(x) / 2) * w[i];
                };
                fbands[i].amplitude = [&fbands, &a, i](space_t space, T x) -> T {
                    if(space == space_t::CHEBY)
                        x = pmmath::acos(x);
                    if(a[2u * i] != a[2u * i + 1u]) 
                        return (((x - fbands[i].start) * a[2u * i + 1u] -
                                (x - fbands[i].stop) * a[2u * i]) /
                                (fbands[i].stop - fbands[i].start)) / pmmath::sin(x / 2);
                    return a[2u * i] / pmmath::sin(x / 2);
                };
            }
        }
    }

    pmoutput_t<T> output;
    std::vector<T> x;
    bandconv(cbands, fbands, convdir_t::FROMFREQ);
    std::function<T(T)> wf = [&cbands](T x) -> T {
        for(std::size_t i{0u}; i < cbands.size(); ++i)
            if(cbands[i].start <= x && x <= cbands[i].stop)
                return cbands[i].weight(space_t::CHEBY, x);
        // this should never execute
        return 1.0;
    };

    if(fbands.size() > (deg + 2u) / 4)
        strategy = init_t::AFP;

    switch(strategy) {
        case init_t::UNIFORM: 
        {
            if (fbands.size() <= (deg + 2u) / 4) { 
                std::vector<T> omega;
                uniform(omega, fbands, deg + 2u);
                cos(x, omega);
                bandconv(cbands, fbands, convdir_t::FROMFREQ);
            } else {
                // use AFP strategy for very small degrees (wrt nb of bands)
                std::vector<T> mesh;
                wam(mesh, cbands, deg);
                MatrixXd<T> A;
                chebvand(A, deg+1u, mesh, wf);
                afp(x, A, mesh);
                if(x.size() != deg + 2u) {
                    std::cerr << "ERROR: AFP strategy failed to produce a valid starting reference\n"
                        << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix\n";
                    exit(EXIT_FAILURE);
                }
                countBand(cbands, x);
            }
            output = exchange(x, cbands, eps, nmax, prec);
        } break;
        case init_t::SCALING: 
        {
            std::vector<std::size_t> sdegs(depth+1u);
            sdegs[depth] = deg;
            for(int i{(int)depth-1}; i >= 0; --i)
                sdegs[i] = sdegs[i+1]/2;

            if(rstrategy == init_t::UNIFORM) {
                if (fbands.size() <= (sdegs[0] + 2u) / 4) {
                    std::vector<T> omega;
                    uniform(omega, fbands, sdegs[0]+2u);
                    cos(x, omega);
                    bandconv(cbands, fbands, convdir_t::FROMFREQ);
                } else {
                    // use AFP strategy for very small degrees (wrt nb of bands)
                    std::vector<T> mesh;
                    wam(mesh, cbands, sdegs[0]);
                    MatrixXd<T> A;
                    chebvand(A, sdegs[0]+1u, mesh, wf);
                    afp(x, A, mesh);
                    if(x.size() != sdegs[0] + 2u) {
                        std::cerr << "ERROR: AFP strategy failed to produce a valid starting reference\n"
                            << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix\n";
                        exit(EXIT_FAILURE);
                    }
                    countBand(cbands, x);                    
                }
                output = exchange(x, cbands, eps, nmax);
            } else { // AFP-based strategy
                std::vector<T> mesh;
                wam(mesh, cbands, sdegs[0]);
                MatrixXd<T> A;
                chebvand(A, sdegs[0]+1u, mesh, wf);
                afp(x, A, mesh);
                if(x.size() != sdegs[0] + 2u) {
                    std::cerr << "ERROR: AFP strategy failed to produce a valid starting reference\n"
                        << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix\n";
                    exit(EXIT_FAILURE);
                }
                countBand(cbands, x);
                output = exchange(x, cbands, eps, nmax);
            }
            for(std::size_t i{1u}; i <= depth && output.q <= 0.5; ++i) {
                x.clear();
                referenceScaling(x, cbands, fbands, sdegs[i]+2u, 
                                 output.x, cbands, fbands);
                output = exchange(x, cbands, eps, nmax, prec);
            }
        } break;
        default: { // AFP-based initialization
            std::vector<T> mesh;
            wam(mesh, cbands, deg);
            MatrixXd<T> A;
            chebvand(A, deg+1u, mesh, wf);
            afp(x, A, mesh);
            if(x.size() != deg + 2u) {
                std::cerr << "ERROR: AFP strategy failed to produce a valid starting reference\n"
                    << "POSSIBLE CAUSE: badly conditioned Chebyshev Vandermonde matrix\n";
                exit(EXIT_FAILURE);
            }
            countBand(cbands, x);
            output = exchange(x, cbands, eps, nmax, prec);
        }
    }

    h.resize(n + 1u);
    if(n % 2 == 0)
    {
        h[deg + 1u] = 0;
        h[deg] = (output.h[0u] * 2.0 - output.h[2]) / 4u;
        h[deg + 2u] = -h[deg];
        h[1u] = output.h[deg - 1u] / 4;
        h[2u * deg + 1u] = -h[1u];
        h[0u] =  output.h[deg] / 4;
        h[2u * (deg + 1u)] = -h[0u];
        for(std::size_t i{2u}; i < deg; ++i)
        {
            h[deg + 1u - i] = (output.h[i - 1u] - output.h[i + 1u]) / 4;
            h[deg + 1u + i] = -h[deg + 1u - i];
        }
    } else {
        ++deg;
        h[deg - 1u] = (output.h[0u] * 2.0 - output.h[1u]) / 4;
        h[deg] = -h[deg - 1u];
        h[0u] = output.h[deg - 1u] / 4;
        h[2u * deg - 1u] = -h[0u];
        for(std::size_t i{2u}; i < deg; ++i)
        {
            h[deg - i] = (output.h[i - 1u] - output.h[i]) / 4;
            h[deg + i - 1u] = -h[deg - i];
        }
    }

    output.h = h;
    return output;
}

template<typename T>
pmoutput_t<T> firpmRS(std::size_t n,
            std::vector<T>const &f,
            std::vector<T>const &a,
            std::vector<T>const &w,
            filter_t type,
	        double eps,
            std::size_t nmax, std::size_t depth,
            init_t rstrategy,
            unsigned long prec)
{
    return firpm<T>(n, f, a, w, type, eps, nmax,
                 init_t::SCALING, depth, rstrategy, prec);
}

template<typename T>
pmoutput_t<T> firpmAFP(std::size_t n,
            std::vector<T>const &f,
            std::vector<T>const &a,
            std::vector<T>const &w,
            filter_t type,
	        double eps,
            std::size_t nmax,
            unsigned long prec)
{
    return firpm<T>(n, f, a, w, type, eps,
                 nmax, init_t::AFP, 0u, init_t::AFP, prec);
}

/* Explicit instantiations, since template code is not in header */

/* double precision */
template void uniform<double>(std::vector<double>& omega,
            std::vector<band_t<double>>& B, std::size_t n);

template void referenceScaling<double>(std::vector<double>& newX,
	        std::vector<band_t<double>>& newChebyBands,
            std::vector<band_t<double>>& newFreqBands,
	        std::size_t newXSize,
            std::vector<double>& x,
	        std::vector<band_t<double>>& chebyBands,
            std::vector<band_t<double>>& freqBands);

template pmoutput_t<double> exchange<double>(std::vector<double>& x,
            std::vector<band_t<double>>& chebyBands, 
            double eps, 
            std::size_t nmax,
            unsigned long prec);

template pmoutput_t<double> firpm<double>(std::size_t n,
            std::vector<double>const &f,
            std::vector<double>const &a,
            std::vector<double>const &w,
            double eps,
            std::size_t nmax,
            init_t strategy,
            std::size_t depth,
            init_t rstrategy,
            unsigned long prec);

template pmoutput_t<double> firpm<double>(std::size_t n,
            std::vector<double>const& f,
            std::vector<double>const& a,
            std::vector<double>const& w,
            filter_t type,
            double eps,
            std::size_t nmax,
            init_t strategy,
            std::size_t depth,
            init_t rstrategy,
            unsigned long prec);

template pmoutput_t<double> firpmRS<double>(std::size_t n,
            std::vector<double>const &f,
            std::vector<double>const &a,
            std::vector<double>const &w,
            double eps, std::size_t nmax,
            std::size_t depth,
            init_t rstrategy,
            unsigned long prec);

template pmoutput_t<double> firpmRS<double>(std::size_t n,
            std::vector<double>const &f,
            std::vector<double>const &a,
            std::vector<double>const &w,
            filter_t type,
            double eps,
            std::size_t nmax,
            std::size_t depth,
            init_t rstrategy,
            unsigned long prec);

template pmoutput_t<double> firpmAFP<double>(std::size_t n,
            std::vector<double>const &f,
            std::vector<double>const &a,
            std::vector<double>const &w,
            double eps, std::size_t nmax,
            unsigned long prec);

template pmoutput_t<double> firpmAFP<double>(std::size_t n,
            std::vector<double>const &f,
            std::vector<double>const &a,
            std::vector<double>const &w,
            filter_t type,
            double eps,
            std::size_t nmax,
            unsigned long prec);

/* long double precision */
template void uniform<long double>(std::vector<long double>& omega,
            std::vector<band_t<long double>>& B, std::size_t n);

template void referenceScaling<long double>(std::vector<long double>& newX,
	        std::vector<band_t<long double>>& newChebyBands,
            std::vector<band_t<long double>>& newFreqBands,
	        std::size_t newXSize,
            std::vector<long double>& x,
	        std::vector<band_t<long double>>& chebyBands,
            std::vector<band_t<long double>>& freqBands);

template pmoutput_t<long double> exchange<long double>(std::vector<long double>& x,
            std::vector<band_t<long double>>& chebyBands, 
            double eps, 
            std::size_t nmax,
            unsigned long prec);

template pmoutput_t<long double> firpm<long double>(std::size_t n,
            std::vector<long double>const &f,
            std::vector<long double>const &a,
            std::vector<long double>const &w,
            double eps,
            std::size_t nmax,
            init_t strategy,
            std::size_t depth,
            init_t rstrategy,
            unsigned long prec);

template pmoutput_t<long double> firpm<long double>(std::size_t n,
            std::vector<long double>const& f,
            std::vector<long double>const& a,
            std::vector<long double>const& w,
            filter_t type,
            double eps,
            std::size_t nmax,
            init_t strategy,
            std::size_t depth,
            init_t rstrategy,
            unsigned long prec);

template pmoutput_t<long double> firpmRS<long double>(std::size_t n,
            std::vector<long double>const &f,
            std::vector<long double>const &a,
            std::vector<long double>const &w,
            double eps, std::size_t nmax,
            std::size_t depth,
            init_t rstrategy,
            unsigned long prec);

template pmoutput_t<long double> firpmRS<long double>(std::size_t n,
            std::vector<long double>const &f,
            std::vector<long double>const &a,
            std::vector<long double>const &w,
            filter_t type,
            double eps,
            std::size_t nmax,
            std::size_t depth,
            init_t rstrategy,
            unsigned long prec);

template pmoutput_t<long double> firpmAFP<long double>(std::size_t n,
            std::vector<long double>const &f,
            std::vector<long double>const &a,
            std::vector<long double>const &w,
            double eps, std::size_t nmax,
            unsigned long prec);

template pmoutput_t<long double> firpmAFP<long double>(std::size_t n,
            std::vector<long double>const &f,
            std::vector<long double>const &a,
            std::vector<long double>const &w,
            filter_t type,
            double eps,
            std::size_t nmax,
            unsigned long prec);

/* multiple precision mpreal */
#ifdef HAVE_MPFR
    template void uniform<mpfr::mpreal>(std::vector<mpfr::mpreal>& omega,
                std::vector<band_t<mpfr::mpreal>>& B, std::size_t n);

    template void referenceScaling<mpfr::mpreal>(std::vector<mpfr::mpreal>& newX,
                std::vector<band_t<mpfr::mpreal>>& newChebyBands,
                std::vector<band_t<mpfr::mpreal>>& newFreqBands,
                std::size_t newXSize,
                std::vector<mpfr::mpreal>& x,
                std::vector<band_t<mpfr::mpreal>>& chebyBands,
                std::vector<band_t<mpfr::mpreal>>& freqBands);

    template pmoutput_t<mpfr::mpreal> exchange<mpfr::mpreal>(std::vector<mpfr::mpreal>& x,
                std::vector<band_t<mpfr::mpreal>>& chebyBands, 
                double eps, 
                std::size_t nmax, unsigned long prec);

    template pmoutput_t<mpfr::mpreal> firpm<mpfr::mpreal>(std::size_t n,
                std::vector<mpfr::mpreal>const &f,
                std::vector<mpfr::mpreal>const &a,
                std::vector<mpfr::mpreal>const &w,
                double eps,
                std::size_t nmax,
                init_t strategy,
                std::size_t depth,
                init_t rstrategy,
                unsigned long prec);

    template pmoutput_t<mpfr::mpreal> firpm<mpfr::mpreal>(std::size_t n,
                std::vector<mpfr::mpreal>const& f,
                std::vector<mpfr::mpreal>const& a,
                std::vector<mpfr::mpreal>const& w,
                filter_t type,
                double eps,
                std::size_t nmax,
                init_t strategy,
                std::size_t depth,
                init_t rstrategy,
                unsigned long prec);

    template pmoutput_t<mpfr::mpreal> firpmRS<mpfr::mpreal>(std::size_t n,
                std::vector<mpfr::mpreal>const &f,
                std::vector<mpfr::mpreal>const &a,
                std::vector<mpfr::mpreal>const &w,
                double eps, std::size_t nmax,
                std::size_t depth,
                init_t rstrategy,
                unsigned long prec);

    template pmoutput_t<mpfr::mpreal> firpmRS<mpfr::mpreal>(std::size_t n,
                std::vector<mpfr::mpreal>const &f,
                std::vector<mpfr::mpreal>const &a,
                std::vector<mpfr::mpreal>const &w,
                filter_t type,
                double eps,
                std::size_t nmax,
                std::size_t depth,
                init_t rstrategy,
                unsigned long prec);

    template pmoutput_t<mpfr::mpreal> firpmAFP<mpfr::mpreal>(std::size_t n,
                std::vector<mpfr::mpreal>const &f,
                std::vector<mpfr::mpreal>const &a,
                std::vector<mpfr::mpreal>const &w,
                double eps, std::size_t nmax,
                unsigned long prec);

    template pmoutput_t<mpfr::mpreal> firpmAFP<mpfr::mpreal>(std::size_t n,
                std::vector<mpfr::mpreal>const &f,
                std::vector<mpfr::mpreal>const &a,
                std::vector<mpfr::mpreal>const &w,
                filter_t type,
                double eps,
                std::size_t nmax,
                unsigned long prec);
#endif