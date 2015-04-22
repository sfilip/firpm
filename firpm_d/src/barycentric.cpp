//    firpm_d
//    Copyright (C) 2015  S. Filip
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

void barycentricWeights(std::vector<double>& w,
        std::vector<double>& x)
{
    if(x.size() > 500u)
    {
        for(std::size_t i = 0u; i < x.size(); ++i)
        {
            double one = 1;
            double denom = 0.0;
            double xi = x[i];
            for(std::size_t j = 0u; j < x.size(); ++j)
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
        double one = 1u;
        for(std::size_t i = 0u; i < x.size(); ++i)
        {
            double denom = 1.0;
            double xi = x[i];
            for(std::size_t j = 0u; j < step; ++j)
            {
                for(std::size_t k = j; k < x.size(); k += step)
                    if (k != i)
                        denom *= ((xi - x[k]) * 2);
            }
            w[i] = one / denom;
        }
    }
}


void computeIdealResponseAndWeight(double &D, double &W,
        const double &xVal, std::vector<Band> &bands)
{
    for (auto &it : bands) {
        if (xVal >= it.start && xVal <= it.stop) {
            D = it.amplitude(it.space, xVal);
            W = it.weight(it.space, xVal);
            return;
        }
    }
}

void computeDelta(double &delta, std::vector<double>& x,
        std::vector<Band> &bands)
{
    std::vector<double> w(x.size());
    barycentricWeights(w, x);

    double num, denom, D, W, buffer;
    num = denom = D = W = 0;
    for (std::size_t i = 0u; i < w.size(); ++i)
    {
        computeIdealResponseAndWeight(D, W, x[i], bands);
        buffer = w[i];
        num += buffer * D;
        buffer = w[i] / W;
        if (i % 2 == 0)
            buffer = -buffer;
        denom += buffer;
    }

    delta = num / denom;
}

void computeDelta(double &delta, std::vector<double>& w,
        std::vector<double>& x, std::vector<Band> &bands)
{
    double num, denom, D, W, buffer;
    num = denom = D = W = 0;
    for (std::size_t i = 0u; i < w.size(); ++i)
    {
        computeIdealResponseAndWeight(D, W, x[i], bands);
        buffer = w[i];
        num += buffer * D;
        buffer = w[i] / W;
        if (i % 2 == 0)
            buffer = -buffer;
        denom += buffer;
    }
    delta = num / denom;
}


void computeC(std::vector<double> &C, double &delta,
        std::vector<double> &omega, std::vector<Band> &bands)
{
    double D, W;
    D = W = 0;
    for (std::size_t i = 0u; i < omega.size(); ++i)
    {
        computeIdealResponseAndWeight(D, W, omega[i], bands);
        if (i % 2 != 0)
            W = -W;
        C[i] = D + (delta / W);
    }
}

void computeApprox(double &Pc, const double &omega,
        std::vector<double> &x, std::vector<double> &C,
        std::vector<double> &w)
{
    double num, denom;
    double buff;
    num = denom = 0;

    Pc = omega;
    std::size_t r = x.size();
    for (std::size_t i = 0u; i < r; ++i)
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

void computeError(double &error, const double &xVal,
        double &delta, std::vector<double> &x,
        std::vector<double> &C, std::vector<double> &w,
        std::vector<Band> &bands)
{
    for (std::size_t i = 0u; i < x.size(); ++i)
    {
        if (xVal == x[i]) {
            if (i % 2 == 0)
                error = delta;
            else
                error = -delta;
            return;
        }
    }

    double D, W;
    D = W = 0;
    computeIdealResponseAndWeight(D, W, xVal, bands);
    computeApprox(error, xVal, x, C, w);
    error -= D;
    error *= W;
}
