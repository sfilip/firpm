//    firpm_d
//    Copyright (C) 2015 - 2019 S. Filip
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

#include "firpm/cheby.h"
#include "firpm/pmmath.h"

template<typename T>
using MatrixXd = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using VectorXcd = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>;

template<typename T>
void balance(MatrixXd<T> &A)
{
    std::size_t n = A.rows();

    T rNorm;      // row norm
    T cNorm;      // column norm
    bool converged = false;

    T g, f, s;
    while(!converged)
    {
        converged = true;
        for(std::size_t i{0u}; i < n; ++i)
        {
            rNorm = cNorm = 0.0;
            for(std::size_t j{0u}; j < n; ++j)
            {
                if(i == j)
                    continue;
                cNorm += pmmath::fabs(A(j, i));
                rNorm += pmmath::fabs(A(i, j));
            }
            if((cNorm == 0.0) || (rNorm == 0))
                continue;

            g = rNorm / 2.0;
            f = 1.0;
            s = cNorm + rNorm;


            while(pmmath::isfinite(cNorm) && cNorm < g)
            {
                f *= 2.0;
                cNorm *= 4.0;
            }

            g = rNorm * 2.0;

            while(pmmath::isfinite(cNorm) && cNorm > g)
            {
                f /= 2.0;
                cNorm /= 4.0;
            }

            if((rNorm + cNorm) < s * f * 0.95)
            {
                converged = false;
                g = T(1.0) / f;
                // multiply by D^{-1} on the left
                A.row(i) *= g;
                // multiply by D on the right
                A.col(i) *= f;
            }

        }
    }
}

template<typename T>
MatrixXd<T> colleague(std::vector<T> const &c,
                   chebkind_t kind,
                   bool bal)
{
    std::vector<T> a{c};
    std::size_t n = c.size() - 1u;
    MatrixXd<T> C(n, n);

    for(std::size_t i{0u}; i < n; ++i)
        for(std::size_t j{0u}; j < n; ++j)
            C(i, j) = 0;

    T denom = -1;
    denom /= a[n];
    denom /= 2;
    for(std::size_t i{0u}; i < c.size() - 1u; ++i)
        a[i] *= denom;
    a[n - 2u] += 0.5;

    for (std::size_t i{0u}; i < n - 1; ++i)
        C(i, i + 1) = C(i + 1, i) = 0.5;
    switch(kind) {
        case FIRST: C(n - 2, n - 1) = 1.0;  break;
        default:    C(n - 2, n - 1) = 0.5;  break;
    }
    for(std::size_t i{0u}; i < n; ++i)
        C(i, 0) = a[n - i - 1u];

    if(bal)
        balance(C);    

    return C;
}

template<typename T>
void cos(std::vector<T>& out,
        std::vector<T> const& in)
{
    out.resize(in.size());
    for(std::size_t i{0u}; i < in.size(); ++i)
        out[i] = pmmath::cos(in[i]);
}

template<typename T>
void chgvar(std::vector<T>& out,
        std::vector<T> const& in,
        T& a, T& b)
{
    out.resize(in.size());
    for(std::size_t i{0u}; i < in.size(); ++i)
        out[i] = (b + a) / 2 + in[i] * (b - a) / 2;
}

template<typename T>
void clenshaw(T &result, const std::vector<T> &p,
        const T &x, const T &a, const T &b)
{
    T bn1, bn2, bn;
    T buffer;

    bn1 = 0;
    bn2 = 0;

    // compute the value of (2*x - b - a)/(b - a) 
    // in the temporary variable buffer
    buffer = (x * 2 - b - a) / (b - a);

    int n = (int)p.size() - 1;
    for(int k{n}; k >= 0; --k) {
        bn = buffer * 2;
        bn = bn * bn1 - bn2 + p[k];
        // update values
        bn2 = bn1;
        bn1 = bn;
    }

    // set the value for the result
    // (i.e., the CI value at x)
    result = bn1 - buffer * bn2;
}

template<typename T>
void clenshaw(T &result,
              const std::vector<T> &p,
              const T &x,
              chebkind_t kind)
{
    T bn1, bn2, bn;

    int n = (int)p.size() - 1;
    bn2 = 0;
    bn1 = p[n];
    for(int k{n - 1}; k >= 1; --k) {
        bn = x * 2;
        bn = bn * bn1 - bn2 + p[k];
        // update values
        bn2 = bn1;
        bn1 = bn;
    }

    if(kind == FIRST)
        result = x * bn1 - bn2 + p[0];
    else
        result = (x * 2) * bn1 - bn2 + p[0];
}

template<typename T>
void equipts(std::vector<T>& v, std::size_t n)
{
    v.resize(n);
    // store the points in the vector v as 
    // v[i] = i * pi / (n-1)
    for(std::size_t i{0u}; i < n; ++i) {
        v[i] = M_PI * i;
        v[i] /= (n-1);
    }
}

// this function computes the values of the coefficients of 
// the CI when Chebyshev nodes of the second kind are used
template<typename T>
void chebcoeffs(std::vector<T>& c,
        std::vector<T>& fv)
{
    std::size_t n = fv.size();
    std::vector<T> v(n);
    equipts(v, n);

    T buffer;

    // halve the first and last coefficients
    T oldValue1 = fv[0];
    T oldValue2 = fv[n-1u];
    fv[0] /= 2;
    fv[n-1u] /= 2;

    for(std::size_t i{0u}; i < n; ++i) {
        // compute the actual value at the Chebyshev
        // node cos(i * pi / n)
        buffer = pmmath::cos(v[i]);
        clenshaw(c[i], fv, buffer, FIRST);

        if(i == 0u || i == n-1u) {
            c[i] /= (n-1u);
        } else {
            c[i] *= 2;
            c[i] /= (n-1u);
        }
    }
    fv[0] = oldValue1;
    fv[n-1u] = oldValue2;

}

// function that generates the coefficients of the 
// derivative of a given CI
template<typename T>
void diffcoeffs(std::vector<T>& dc,
                std::vector<T>& c,
                chebkind_t kind)
{
    dc.resize(c.size()-1);
    switch(kind) {
        case FIRST: {
            int n = c.size() - 2;
            dc[n] = c[n + 1] * (2 * (n + 1));
            dc[n - 1] = c[n] * (2 * n);
            for(int i{n - 2}; i >= 0; --i) {
                dc[i] = 2 * (i + 1);
                dc[i] = dc[i] * c[i + 1] + dc[i + 2];
            }
            dc[0] /= 2;
        };
        break;
        default: {
            int n = c.size() - 1;
            for(int i{n}; i > 0; --i)
                dc[i - 1] = c[i] * i;
        }
        break;
    }
}

template<typename T>
void roots(std::vector<T>& r, std::vector<T>& c,
           std::pair<T, T> const &dom,
           chebkind_t kind,
           bool balance)
{
    r.clear();
    for(auto &it : c)
        if(!pmmath::isfinite(it))
            return;

    MatrixXd<T> C = colleague(c, kind, balance);
    Eigen::EigenSolver<MatrixXd<T>> es(C);
    VectorXcd<T> eigs = es.eigenvalues();

    T threshold = 1e-20;
    for(Eigen::Index i{0}; i < eigs.size(); ++i) {
        if(pmmath::fabs(eigs(i).imag()) < threshold)
            if(dom.first < eigs(i).real() && 
               dom.second >= eigs(i).real()) {
                r.push_back(eigs(i).real());
            }
    }

    std::sort(begin(r), end(r));
}

/* Explicit instantiation */

/* double precision */

template void cos<double>(std::vector<double>& out,
        std::vector<double> const& in);

template void chgvar<double>(std::vector<double>& out,
        std::vector<double> const& in,
        double& a, double& b);

template void equipts<double>(std::vector<double>& v, std::size_t n);


template void chebcoeffs<double>(std::vector<double>& c,
                std::vector<double>& fv);

template void diffcoeffs<double>(std::vector<double>& dc,
                std::vector<double>& c,
                chebkind_t kind);

template void roots<double>(std::vector<double>& r, std::vector<double>& c,
           std::pair<double, double> const &dom,
           chebkind_t kind, bool balance);

/* long double precision */

template void cos<long double>(std::vector<long double>& out,
        std::vector<long double> const& in);

template void chgvar<long double>(std::vector<long double>& out,
        std::vector<long double> const& in,
        long double& a, long double& b);

template void equipts<long double>(std::vector<long double>& v, std::size_t n);


template void chebcoeffs<long double>(std::vector<long double>& c,
                std::vector<long double>& fv);

template void diffcoeffs<long double>(std::vector<long double>& dc,
                std::vector<long double>& c,
                chebkind_t kind);

template void roots<long double>(std::vector<long double>& r, std::vector<long double>& c,
           std::pair<long double, long double> const &dom,
           chebkind_t kind, bool balance);

#ifdef HAVE_MPFR
    template void cos<mpfr::mpreal>(std::vector<mpfr::mpreal>& out,
            std::vector<mpfr::mpreal> const& in);

    template void chgvar<mpfr::mpreal>(std::vector<mpfr::mpreal>& out,
            std::vector<mpfr::mpreal> const& in,
            mpfr::mpreal& a, mpfr::mpreal& b);

    template void equipts<mpfr::mpreal>(std::vector<mpfr::mpreal>& v, std::size_t n);


    template void chebcoeffs<mpfr::mpreal>(std::vector<mpfr::mpreal>& c,
                    std::vector<mpfr::mpreal>& fv);

    template void diffcoeffs<mpfr::mpreal>(std::vector<mpfr::mpreal>& dc,
                    std::vector<mpfr::mpreal>& c,
                    chebkind_t kind);

    template void roots<mpfr::mpreal>(std::vector<mpfr::mpreal>& r, std::vector<mpfr::mpreal>& c,
            std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
            chebkind_t kind, bool balance);
#endif