//    firpm_d
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
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "firpm/barycentric.h"
#include "firpm/pmmath.h"

namespace pm {

    template<typename T>
    std::vector<T> baryweights(std::vector<T> const& x)
    {
        std::vector<T> w(x.size());
        if(x.size() > 500u)
        {
            for(std::size_t i{0u}; i < x.size(); ++i)
            {
                T one = 1;
                T denom = 0.0;
                T xi = x[i];
                for(std::size_t j{0u}; j < x.size(); ++j)
                {
                    if (j != i) {
                        denom += pmmath::log(((xi - x[j] > 0) ? (xi - x[j]) : (x[j] - xi)));
                        one *= ((xi - x[j] > 0) ? 1 : -1);
                    }
                }
                w[i] = one / pmmath::exp(denom + pmmath::log(2.0)* (x.size() - 1));
            }
        }
        else
        {
            std::size_t step = (x.size() - 2) / 15 + 1;
            T one = 1u;
            for(std::size_t i{0u}; i < x.size(); ++i)
            {
                T denom = 1.0;
                T xi = x[i];
                for(std::size_t j{0u}; j < step; ++j)
                {
                    for(std::size_t k{j}; k < x.size(); k += step)
                        if (k != i)
                            denom *= ((xi - x[k]) * 2);
                }
                w[i] = one / denom;
            }
        }
        return w;
    }

    template<typename T>
    std::pair<T, T> idealvals(T const& x, std::vector<band_t<T>> const& bands)
    {
        T D, W;
        for (auto &it : bands) {
            if (x >= it.start && x <= it.stop) {
                D = it.amplitude(it.space, x);
                W = it.weight(it.space, x);
                return std::pair(D, W);
            }
        }
        return std::pair(0, 0);
    }

    template<typename T>
    T compdelta(std::vector<T>& x,
            std::vector<band_t<T>> const& bands)
    {
        std::vector<T> w = baryweights(x);

        T num, denom, buffer;
        num = denom = 0;
        for (std::size_t i{0u}; i < w.size(); ++i)
        {
            auto [D, W] = idealvals(x[i], bands);
            buffer = w[i];
            num += buffer * D;
            buffer = w[i] / W;
            if (i % 2 == 0)
                buffer = -buffer;
            denom += buffer;
        }

        return num / denom;
    }

    template<typename T>
    T compdelta(std::vector<T>& w,
            std::vector<T>& x, std::vector<band_t<T>> const& bands)
    {
        T num, denom, buffer;
        num = denom = 0;
        for (std::size_t i{0u}; i < w.size(); ++i)
        {
            auto [D, W] = idealvals(x[i], bands);
            buffer = w[i];
            num += buffer * D;
            buffer = w[i] / W;
            if (i % 2 == 0)
                buffer = -buffer;
            denom += buffer;
        }
        return num / denom;
    }


    template<typename T>
    std::vector<T> compc(T& delta, std::vector<T>& omega,
            std::vector<band_t<T>> const& bands)
    {
        std::vector<T> C(omega.size());

        for (std::size_t i{0u}; i < omega.size(); ++i)
        {
            auto [D, W] = idealvals(omega[i], bands);
            if (i % 2 != 0)
                W = -W;
            C[i] = D + (delta / W);
        }

        return C;
    }

    template<typename T>
    T approx(T const& omega,
            std::vector<T> const& x,
            std::vector<T> const& C,
            std::vector<T> const& w)
    {
        T num, denom;
        T buff;
        num = denom = 0;

        T Pc = omega;
        std::size_t r = x.size();
        for (std::size_t i{0u}; i < r; ++i)
        {
            if (Pc == x[i])
                return C[i];

            buff = w[i] / (Pc - x[i]);
            num += buff * C[i];
            denom += buff;
        }
        return num / denom;
    }

    template<typename T>
    T comperror(T const& xVal,
            T& delta,
            std::vector<T> const& x,
            std::vector<T> const& C,
            std::vector<T> const& w,
            std::vector<band_t<T>> const& bands)
    {
        for (std::size_t i{0u}; i < x.size(); ++i)
        {
            if (xVal == x[i]) {
                if (i % 2 == 0)
                    return delta;
                else
                    return -delta;
            }
        }

        auto [D, W] = idealvals(xVal, bands);
        T error = approx(xVal, x, C, w);
        error -= D;
        error *= W;
        return error;
    }

    /* Template instantiations */

    /* double precision */

    template std::vector<double> baryweights<double>(std::vector<double> const& x);

    template double compdelta<double>(std::vector<double>& x,
            std::vector<band_t<double>> const& bands);

    template double compdelta<double>(std::vector<double>& w,
            std::vector<double>& x,
            std::vector<band_t<double>> const& bands);

    template std::vector<double> compc<double>(double& delta,
            std::vector<double>& x, std::vector<band_t<double>> const& bands);

    template double approx<double>(
            double const& xVal,
            std::vector<double> const& x,
            std::vector<double> const& C,
            std::vector<double> const& w);

    template double comperror<double>(double const& xVal,
            double& delta,
            std::vector<double> const& x,
            std::vector<double> const& C,
            std::vector<double> const& w,
            std::vector<band_t<double>> const& bands);

    /* long double precision */

    template std::vector<long double> baryweights<long double>(std::vector<long double> const& x);

    template long double compdelta<long double>(
            std::vector<long double>& x, std::vector<band_t<long double>> const& bands);

    template long double compdelta<long double>(
            std::vector<long double>& w, std::vector<long double>& x,
            std::vector<band_t<long double>> const& bands);

    template std::vector<long double> compc<long double>(long double& delta,
            std::vector<long double> &x, std::vector<band_t<long double>> const& bands);

    template long double approx<long double>(long double const& xVal,
            std::vector<long double> const& x,
            std::vector<long double> const& C,
            std::vector<long double> const& w);

    template long double comperror<long double>(long double const& xVal,
            long double& delta,
            std::vector<long double> const& x,
            std::vector<long double> const& C,
            std::vector<long double> const& w,
            std::vector<band_t<long double>> const& bands);

#ifdef HAVE_MPFR
    // separate implementation for the MPFR version; it is much faster and
    // the higher precision should usually compensate for any eventual
    // ill-conditioning
    template<> std::vector<mpfr::mpreal> baryweights<mpfr::mpreal>(std::vector<mpfr::mpreal> const& x)
    {
        std::vector<mpfr::mpreal> w(x.size());
        std::size_t step = (x.size() - 2u) / 15 + 1;
        mpfr::mpreal one = 1u;
        for(std::size_t i{0u}; i < x.size(); ++i)
        {
            mpfr::mpreal denom = 1.0;
            mpfr::mpreal xi = x[i];
            for(std::size_t j{0u}; j < step; ++j)
            {
                for(std::size_t k{j}; k < x.size(); k += step)
                    if (k != i)
                        denom *= ((xi - x[k]) << 1);
            }
            w[i] = one / denom;
        }
        return w;
    }
    template mpfr::mpreal compdelta<mpfr::mpreal>(
            std::vector<mpfr::mpreal>& x, std::vector<band_t<mpfr::mpreal>> const& bands);

    template mpfr::mpreal compdelta<mpfr::mpreal>(
            std::vector<mpfr::mpreal>& w, std::vector<mpfr::mpreal>& x,
            std::vector<band_t<mpfr::mpreal>> const& bands);

    template std::vector<mpfr::mpreal> compc<mpfr::mpreal>(mpfr::mpreal& delta,
            std::vector<mpfr::mpreal>& x, std::vector<band_t<mpfr::mpreal>> const& bands);

    template mpfr::mpreal approx<mpfr::mpreal>(mpfr::mpreal const& xVal,
            std::vector<mpfr::mpreal> const& x,
            std::vector<mpfr::mpreal> const& C,
            std::vector<mpfr::mpreal> const& w);

    template mpfr::mpreal comperror<mpfr::mpreal>(mpfr::mpreal const& xVal,
            mpfr::mpreal& delta,
            std::vector<mpfr::mpreal> const& x,
            std::vector<mpfr::mpreal> const& C,
            std::vector<mpfr::mpreal> const& w,
            std::vector<band_t<mpfr::mpreal>> const& bands);
#endif

} // namespace pm
