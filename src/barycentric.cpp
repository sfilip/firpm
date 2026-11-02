//    firpm
//    Copyright (C) 2015 - 2024  S. Filip

#include "firpm/barycentric.h"
#include "firpm/pmmath.h"

#include <algorithm>
#include <numeric>

namespace pm {

    template<typename T>
    void baryweights(std::vector<T>& w,
            std::vector<T>& x)
    {
        // w[i] = 1 / prod_{j != i} 2 (x[i] - x[j]).  The product can leave
        // the exponent range of T, so its exponent is carried separately.
        #pragma omp parallel for
        for(std::ptrdiff_t i = 0; i < std::ssize(x); ++i)
        {
            T m = 1;
            int e = 0;
            for(std::ptrdiff_t j = 0; j < std::ssize(x); ++j)
            {
                if (j != i) {
                    int k;
                    m = pmmath::frexp(m * ((x[i] - x[j]) * 2), &k);
                    e += k;
                }
            }
            w[i] = pmmath::ldexp(1 / m, -e);
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
        compdelta(delta, w, x, bands);
    }

    // delta = sum_i w_i D_i / sum_i (-1)^i w_i / W_i.
    //
    // The barycentric weights of n >= 2 distinct nodes sum to zero, so the
    // numerator also equals sum_i w_i (D_i - c) for any constant c, and its
    // rounding error is bounded by eps * sum_i |w_i| |D_i - c|.  The
    // |w|-weighted median of the D_i minimises that bound.  Summing w_i D_i
    // directly cancels catastrophically when delta is much smaller than the
    // D_i, as at the starting reference of a filter whose reference set lies
    // nearly all at one level.
    template<typename T>
    void compdelta(T& delta, std::vector<T>& w,
            std::vector<T>& x, std::vector<band_t<T>>& bands)
    {
        std::size_t n = w.size();
        std::vector<T> D(n), W(n);
        for (std::size_t i{0u}; i < n; ++i)
            idealvals(D[i], W[i], x[i], bands);

        std::vector<std::size_t> order(n);
        std::iota(order.begin(), order.end(), 0u);
        std::ranges::sort(order, {}, [&D](std::size_t i) { return D[i]; });
        T total = 0;
        for (auto const& wi : w)
            total += pmmath::fabs(wi);
        T mass = 0;
        T c = 0;
        for (auto i : order) {
            c = D[i];
            mass += pmmath::fabs(w[i]);
            if (mass * 2 >= total)
                break;
        }

        T num = 0, denom = 0;
        for (std::size_t i{0u}; i < n; ++i)
        {
            num += w[i] * (D[i] - c);
            T buffer = w[i] / W[i];
            denom += (i % 2 == 0) ? -buffer : buffer;
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
    void compref(std::vector<T>& D, std::vector<T>& S,
            std::vector<T>& omega, std::vector<band_t<T>>& bands)
    {
        T W = 0;
        for (std::size_t i{0u}; i < omega.size(); ++i)
        {
            idealvals(D[i], W, omega[i], bands);
            if (i % 2 != 0)
                W = -W;
            S[i] = 1 / W;
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

    // Error evaluation for the exchange.  P(xVal) - D = sum_i b_i (C_i - D)
    // / sum_i b_i with b_i = w_i / (xVal - x_i), since the barycentric
    // formula reproduces constants; subtracting D inside the sum rather than
    // from P(xVal) keeps the error accurate when it is far smaller than D.
    // The reference values C_i = D_i + delta s_i arrive from compref() as D_i
    // and s_i = (-1)^i / W_i, so that C_i - D = (D_i - D) + delta s_i is
    // formed without first rounding delta s_i into D_i.
    template<typename T>
    void comperror(T& error, T const& xVal,
            T& delta, std::vector<T>& x,
            std::vector<T>& D, std::vector<T>& S, std::vector<T>& w,
            std::vector<band_t<T>>& bands)
    {
        for (std::size_t i{0u}; i < x.size(); ++i)
        {
            if (xVal == x[i]) {
                error = (i % 2 == 0) ? delta : -delta;
                return;
            }
        }

        T Dx, W;
        Dx = W = 0;
        idealvals(Dx, W, xVal, bands);

        T num = 0, denom = 0;
        for (std::size_t i{0u}; i < x.size(); ++i)
        {
            T buff = w[i] / (xVal - x[i]);
            num += buff * ((D[i] - Dx) + delta * S[i]);
            denom += buff;
        }
        error = W * (num / denom);
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

    template void compref<double>(std::vector<double>& D, std::vector<double>& S,
            std::vector<double>& x, std::vector<band_t<double>>& bands);

    template void approx<double>(double& Pc, double const& xVal,
            std::vector<double>& x, std::vector<double>& C,
            std::vector<double>& w);

    template void comperror<double>(double& error, double const& xVal,
            double& delta, std::vector<double>& x,
            std::vector<double>& C, std::vector<double>& w,
            std::vector<band_t<double>>& bands);

    template void comperror<double>(double& error, double const& xVal,
            double& delta, std::vector<double>& x,
            std::vector<double>& D, std::vector<double>& S, std::vector<double>& w,
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

    template void compref<long double>(std::vector<long double>& D, std::vector<long double>& S,
            std::vector<long double>& x, std::vector<band_t<long double>>& bands);

    template void approx<long double>(long double& Pc, long double const& xVal,
            std::vector<long double>& x, std::vector<long double>& C,
            std::vector<long double>& w);

    template void comperror<long double>(long double& error, long double const& xVal,
            long double& delta, std::vector<long double>& x,
            std::vector<long double>& C, std::vector<long double>& w,
            std::vector<band_t<long double>>& bands);

    template void comperror<long double>(long double& error, long double const& xVal,
            long double& delta, std::vector<long double>& x,
            std::vector<long double>& D, std::vector<long double>& S, std::vector<long double>& w,
            std::vector<band_t<long double>>& bands);

#ifdef HAVE_MPFR
    // The generic template serves mpreal: its worker threads start at the
    // default precision, but every operation takes the larger precision of
    // its operands, so the weights come out at the precision of x.
    template void baryweights<mpfr::mpreal>(std::vector<mpfr::mpreal>& w,
            std::vector<mpfr::mpreal>& x);
    template void compdelta<mpfr::mpreal>(mpfr::mpreal& delta,
            std::vector<mpfr::mpreal>& x, std::vector<band_t<mpfr::mpreal>>& bands);

    template void compdelta<mpfr::mpreal>(mpfr::mpreal& delta,
            std::vector<mpfr::mpreal>& w, std::vector<mpfr::mpreal>& x,
            std::vector<band_t<mpfr::mpreal>>& bands);

    template void compc<mpfr::mpreal>(std::vector<mpfr::mpreal>& C, mpfr::mpreal& delta,
            std::vector<mpfr::mpreal>& x, std::vector<band_t<mpfr::mpreal>>& bands);

    template void compref<mpfr::mpreal>(std::vector<mpfr::mpreal>& D, std::vector<mpfr::mpreal>& S,
            std::vector<mpfr::mpreal>& x, std::vector<band_t<mpfr::mpreal>>& bands);

    template void approx<mpfr::mpreal>(mpfr::mpreal& Pc, mpfr::mpreal const& xVal,
            std::vector<mpfr::mpreal>& x, std::vector<mpfr::mpreal>& C,
            std::vector<mpfr::mpreal>& w);

    template void comperror<mpfr::mpreal>(mpfr::mpreal& error, mpfr::mpreal const& xVal,
            mpfr::mpreal& delta, std::vector<mpfr::mpreal>& x,
            std::vector<mpfr::mpreal>& C, std::vector<mpfr::mpreal>& w,
            std::vector<band_t<mpfr::mpreal>>& bands);

    template void comperror<mpfr::mpreal>(mpfr::mpreal& error, mpfr::mpreal const& xVal,
            mpfr::mpreal& delta, std::vector<mpfr::mpreal>& x,
            std::vector<mpfr::mpreal>& D, std::vector<mpfr::mpreal>& S, std::vector<mpfr::mpreal>& w,
            std::vector<band_t<mpfr::mpreal>>& bands);
#endif

} // namespace pm