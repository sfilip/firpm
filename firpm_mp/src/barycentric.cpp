//    firpm_mp
//    Copyright (C) 2015 -2019  S. Filip
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



#include "firpm/barycentric.h"


void baryweights(std::vector<mpfr::mpreal>& w,
        std::vector<mpfr::mpreal>& x, mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    std::size_t step = (x.size() - 2u) / 15 + 1;
    mpreal one = 1u;
    for(std::size_t i{0u}; i < x.size(); ++i)
    {
        mpreal denom = 1.0;
        mpreal xi = x[i];
        for(std::size_t j{0u}; j < step; ++j)
        {
            for(std::size_t k{j}; k < x.size(); k += step)
                if (k != i)
                    denom *= ((xi - x[k]) << 1);
        }
        w[i] = one / denom;
    }

    mpreal::set_default_prec(prevPrec);
}

void idealvals(mpfr::mpreal &D, mpfr::mpreal &W,
        const mpfr::mpreal &x, std::vector<band_t> &bands,
        mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    for (auto &it : bands) {
        if (x >= it.start && x <= it.stop) {

            D = it.amplitude(it.space, x);
            W = it.weight(it.space, x);
            return;
        }
    }

    mpreal::set_default_prec(prevPrec);
}

void compdelta(mpfr::mpreal &delta, std::vector<mpfr::mpreal>& x,
        std::vector<band_t> &bands, mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);


    std::vector<mpfr::mpreal> w(x.size());
    baryweights(w, x);

    mpfr::mpreal num, denom, D, W, buffer;
    num = denom = D = W = 0;
    for (std::size_t i{0u}; i < w.size(); ++i)
    {
        idealvals(D, W, x[i], bands);
        buffer = w[i];
        num = mpfr::fma(buffer, D, num);
        buffer = w[i] / W;
        if (i % 2 == 0)
            buffer = -buffer;
        denom += buffer;
    }
    delta = num / denom;

    mpreal::set_default_prec(prevPrec);
}

void compdelta(mpfr::mpreal &delta, std::vector<mpfr::mpreal>& w,
        std::vector<mpfr::mpreal>& x, std::vector<band_t> &bands,
        mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    mpfr::mpreal num, denom, D, W, buffer;
    num = denom = D = W = 0;
    for (std::size_t i{0u}; i < w.size(); ++i)
    {
        idealvals(D, W, x[i], bands);
        buffer = w[i];
        num = mpfr::fma(buffer, D, num);
        buffer = w[i] / W;
        if (i % 2 == 0)
            buffer = -buffer;
        denom += buffer;
    }
    delta = num / denom;

    mpreal::set_default_prec(prevPrec);
}


void compc(std::vector<mpfr::mpreal> &C, mpfr::mpreal &delta,
        std::vector<mpfr::mpreal> &omega, std::vector<band_t> &bands,
        mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    mpfr::mpreal D, W;
    D = W = 0;
    for (std::size_t i{0u}; i < omega.size(); ++i)
    {
        idealvals(D, W, omega[i], bands);
        if (i % 2 != 0)
            W = -W;
        C[i] = D + (delta / W);
    }

    mpreal::set_default_prec(prevPrec);
}

void approx(mpfr::mpreal &Pc, const mpfr::mpreal &omega,
        std::vector<mpfr::mpreal> &x, std::vector<mpfr::mpreal> &C,
        std::vector<mpfr::mpreal> &w, mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    mpfr::mpreal num, denom;
    mpfr::mpreal buff;
    num = denom = 0;

    Pc = omega;
    std::size_t r = x.size();
    for (std::size_t i{0u}; i < r; ++i)
    {
        if (Pc == x[i]) {
            Pc = C[i];
            mpreal::set_default_prec(prevPrec);
            return;
        }
        buff = w[i] / (Pc - x[i]);
        num = fma(buff, C[i], num);
        denom += buff;
    }
    Pc = num / denom;

    mpreal::set_default_prec(prevPrec);
}

void comperror(mpfr::mpreal &error, const mpfr::mpreal &xVal,
        mpfr::mpreal &delta, std::vector<mpfr::mpreal> &x,
        std::vector<mpfr::mpreal> &C, std::vector<mpfr::mpreal> &w,
        std::vector<band_t> &bands, mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    for (std::size_t i{0u}; i < x.size(); ++i)
    {
        if (xVal == x[i]) {
            if (i % 2 == 0)
                error = delta;
            else
                error = -delta;
            mpreal::set_default_prec(prevPrec);
            return;
        }
    }

    mpfr::mpreal D, W;
    D = W = 0;
    idealvals(D, W, xVal, bands);
    approx(error, xVal, x, C, w);
    error -= D;
    error *= W;
    
    mpreal::set_default_prec(prevPrec);
}
