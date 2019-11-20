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

template<typename T>
void baryweights(std::vector<T>& w,
        std::vector<T>& x)
{
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
                    denom += log(((xi - x[j] > 0) ? (xi - x[j]) : (x[j] - xi)));
                    one *= ((xi - x[j] > 0) ? 1 : -1);
                }
            }
            w[i] = one / exp(denom + log(2.0)* (x.size() - 1));
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
}


template<typename T>
void idealvals(T &D, T &W,
        const T &x, std::vector<band_t<T>> &bands)
{
    for (auto &it : bands) {
        if (x >= it.start && x <= it.stop) {
            D = it.amplitude(it.space, x);
            W = it.weight(it.space, x);
            return;
        }
    }
}

template<typename T>
void compdelta(T &delta, std::vector<T>& x,
        std::vector<band_t<T>> &bands)
{
    std::vector<T> w(x.size());
    baryweights(w, x);

    T num, denom, D, W, buffer;
    num = denom = D = W = 0;
    for (std::size_t i{0u}; i < w.size(); ++i)
    {
        idealvals(D, W, x[i], bands);
        buffer = w[i];
        num += buffer * D;
        buffer = w[i] / W;
        if (i % 2 == 0)
            buffer = -buffer;
        denom += buffer;
    }

    delta = num / denom;
}

template<typename T>
void compdelta(T &delta, std::vector<T>& w,
        std::vector<T>& x, std::vector<band_t<T>> &bands)
{
    T num, denom, D, W, buffer;
    num = denom = D = W = 0;
    for (std::size_t i{0u}; i < w.size(); ++i)
    {
        idealvals(D, W, x[i], bands);
        buffer = w[i];
        num += buffer * D;
        buffer = w[i] / W;
        if (i % 2 == 0)
            buffer = -buffer;
        denom += buffer;
    }
    delta = num / denom;
}


template<typename T>
void compc(std::vector<T> &C, T &delta,
        std::vector<T> &omega, std::vector<band_t<T>> &bands)
{
    T D, W;
    D = W = 0;
    for (std::size_t i{0u}; i < omega.size(); ++i)
    {
        idealvals(D, W, omega[i], bands);
        if (i % 2 != 0)
            W = -W;
        C[i] = D + (delta / W);
    }
}

template<typename T>
void approx(T &Pc, const T &omega,
        std::vector<T> &x, std::vector<T> &C,
        std::vector<T> &w)
{
    T num, denom;
    T buff;
    num = denom = 0;

    Pc = omega;
    std::size_t r = x.size();
    for (std::size_t i{0u}; i < r; ++i)
    {
        if (Pc == x[i]) {
            Pc = C[i];
            return;
        }
        buff = w[i] / (Pc - x[i]);
        num += buff * C[i];
        denom += buff;
    }
    Pc = num / denom;
}

template<typename T>
void comperror(T &error, const T &xVal,
        T &delta, std::vector<T> &x,
        std::vector<T> &C, std::vector<T> &w,
        std::vector<band_t<T>> &bands)
{
    for (std::size_t i{0u}; i < x.size(); ++i)
    {
        if (xVal == x[i]) {
            if (i % 2 == 0)
                error = delta;
            else
                error = -delta;
            return;
        }
    }

    T D, W;
    D = W = 0;
    idealvals(D, W, xVal, bands);
    approx(error, xVal, x, C, w);
    error -= D;
    error *= W;
}

/* Template instantiations */
template void baryweights<double>(std::vector<double>& w,
		std::vector<double>& x);

template void compdelta<double>(double &delta,
		std::vector<double>& x, std::vector<band_t<double>> &bands);

template void compdelta<double>(double &delta,
		std::vector<double>& w, std::vector<double>& x,
		std::vector<band_t<double>> &bands);

template void compc<double>(std::vector<double> &C, double &delta,
		std::vector<double> &x, std::vector<band_t<double>> &bands);

template void approx<double>(double &Pc, const double &xVal,
		std::vector<double> &x, std::vector<double> &C,
		std::vector<double> & w);

template void comperror<double>(double &error, const double &xVal,
		double &delta, std::vector<double> &x,
		std::vector<double> &C, std::vector<double> &w,
		std::vector<band_t<double>> &bands);

