//    firpm
//    Copyright (C) 2015 - 2024  S. Filip

#include "firpm/barycentric.h"
#include "firpm/pmmath.h"

namespace pm {

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
    }


    template<typename T>
    void idealvals(T& D, T& W,
            T const& x, std::vector<band_t<T>>& bands)
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
    void compdelta(T& delta, std::vector<T>& x,
            std::vector<band_t<T>>& bands)
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
    void compdelta(T& delta, std::vector<T>& w,
            std::vector<T>& x, std::vector<band_t<T>>& bands)
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
    void compc(std::vector<T>& C, T& delta,
            std::vector<T>& omega, std::vector<band_t<T>>& bands)
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
    void approx(T& Pc, T const& omega,
            std::vector<T>& x, std::vector<T>& C,
            std::vector<T>& w)
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
    void comperror(T& error, T const& xVal,
            T& delta, std::vector<T>& x,
            std::vector<T>& C, std::vector<T>& w,
            std::vector<band_t<T>>& bands)
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

    /* double precision */

    template void baryweights<double>(std::vector<double>& w,
            std::vector<double>& x);

    template void compdelta<double>(double& delta,
            std::vector<double>& x, std::vector<band_t<double>>& bands);

    template void compdelta<double>(double& delta,
            std::vector<double>& w, std::vector<double>& x,
            std::vector<band_t<double>>& bands);

    template void compc<double>(std::vector<double>& C, double& delta,
            std::vector<double>& x, std::vector<band_t<double>>& bands);

    template void approx<double>(double& Pc, double const& xVal,
            std::vector<double>& x, std::vector<double>& C,
            std::vector<double>& w);

    template void comperror<double>(double& error, double const& xVal,
            double& delta, std::vector<double>& x,
            std::vector<double>& C, std::vector<double>& w,
            std::vector<band_t<double>>& bands);

    /* long double precision */

    template void baryweights<long double>(std::vector<long double>& w,
            std::vector<long double>& x);

    template void compdelta<long double>(long double& delta,
            std::vector<long double>& x, std::vector<band_t<long double>>& bands);

    template void compdelta<long double>(long double& delta,
            std::vector<long double>& w, std::vector<long double>& x,
            std::vector<band_t<long double>>& bands);

    template void compc<long double>(std::vector<long double>& C, long double& delta,
            std::vector<long double> &x, std::vector<band_t<long double>>& bands);

    template void approx<long double>(long double& Pc, long double const& xVal,
            std::vector<long double>& x, std::vector<long double>& C,
            std::vector<long double>& w);

    template void comperror<long double>(long double& error, long double const& xVal,
            long double& delta, std::vector<long double>& x,
            std::vector<long double>& C, std::vector<long double>& w,
            std::vector<band_t<long double>>& bands);

#ifdef HAVE_MPFR
    // separate implementation for the MPFR version; it is much faster and
    // the higher precision should usually compensate for any eventual
    // ill-conditioning
    template<> void baryweights<mpfr::mpreal>(std::vector<mpfr::mpreal>& w,
            std::vector<mpfr::mpreal>& x)
    {
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
    }
    template void compdelta<mpfr::mpreal>(mpfr::mpreal& delta,
            std::vector<mpfr::mpreal>& x, std::vector<band_t<mpfr::mpreal>>& bands);

    template void compdelta<mpfr::mpreal>(mpfr::mpreal& delta,
            std::vector<mpfr::mpreal>& w, std::vector<mpfr::mpreal>& x,
            std::vector<band_t<mpfr::mpreal>>& bands);

    template void compc<mpfr::mpreal>(std::vector<mpfr::mpreal>& C, mpfr::mpreal& delta,
            std::vector<mpfr::mpreal>& x, std::vector<band_t<mpfr::mpreal>>& bands);

    template void approx<mpfr::mpreal>(mpfr::mpreal& Pc, mpfr::mpreal const& xVal,
            std::vector<mpfr::mpreal>& x, std::vector<mpfr::mpreal>& C,
            std::vector<mpfr::mpreal>& w);

    template void comperror<mpfr::mpreal>(mpfr::mpreal& error, mpfr::mpreal const& xVal,
            mpfr::mpreal& delta, std::vector<mpfr::mpreal>& x,
            std::vector<mpfr::mpreal>& C, std::vector<mpfr::mpreal>& w,
            std::vector<band_t<mpfr::mpreal>>& bands);
#endif

} // namespace pm