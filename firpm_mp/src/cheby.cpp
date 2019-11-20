//    firpm_mp
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

/** Eigen matrix container for mpfr::mpreal values */
typedef Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> MatrixXmp;
/** Eigen vector container for mpfr::mpreal values */
typedef Eigen::Matrix<std::complex<mpfr::mpreal>, Eigen::Dynamic, 1> VectorXcmp;

void balance(MatrixXmp& A, mp_prec_t prec)
{
    using mpfr::mpreal;
    mpfr_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    std::size_t n = A.rows();

    mpfr::mpreal rNorm;      // row norm
    mpfr::mpreal cNorm;      // column norm
    bool converged = false;
    mpfr::mpreal one = mpfr::mpreal(1.0);

    mpfr::mpreal g, f, s;
    while(!converged)
    {
        converged = true;
        for(std::size_t i = 0u; i < n; ++i)
        {
            rNorm = cNorm = 0.0;
            for(std::size_t j = 0u; j < n; ++j)
            {
                if(i == j)
                    continue;
                cNorm += mpfr::fabs(A(j, i));
                rNorm += mpfr::fabs(A(i, j));
            }
            if((cNorm == 0.0) || (rNorm == 0))
                continue;

            g = rNorm >> 1u;
            f = 1.0;
            s = cNorm + rNorm;


            while(mpfr::isfinite(cNorm) && cNorm < g)
            {
                f <<= 1u;
                cNorm <<= 2u;
            }

            g = rNorm << 1u;

            while(mpfr::isfinite(cNorm) && cNorm > g)
            {
                f >>= 1u;
                cNorm >>= 2u;
            }

            if((rNorm + cNorm) < s * f * 0.95)
            {
                converged = false;
                g = one / f;
                // multiply by D^{-1} on the left
                A.row(i) *= g;
                // multiply by D on the right
                A.col(i) *= f;
            }

        }
    }

    mpreal::set_default_prec(prevPrec);
}

MatrixXmp colleague(std::vector<mpfr::mpreal> const &c,
                   chebkind_t kind,
                   bool bal, mp_prec_t prec)
{
    using mpfr::mpreal;
    mpfr_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    std::vector<mpfr::mpreal> a{c};
    std::size_t n = c.size() - 1u;
    MatrixXmp C(n, n);

    for(std::size_t i{0u}; i < n; ++i)
        for(std::size_t j{0u}; j < n; ++j)
            C(i, j) = 0;

    mpfr::mpreal denom = -1;
    denom /= a[n];
    denom >>= 1;
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
        balance(C, prec);   

    mpreal::set_default_prec(prevPrec); 

    return C;
}

void cos(std::vector<mpfr::mpreal>& out,
        std::vector<mpfr::mpreal> const& in,
        mp_prec_t prec)
{
    using mpfr::mpreal;
    mpfr_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);
    
    out.resize(in.size());
    for (std::size_t i{0u}; i < in.size(); ++i)
        out[i] = mpfr::cos(in[i]);

    mpreal::set_default_prec(prevPrec); 
}

void chgvar(std::vector<mpfr::mpreal>& out,
        std::vector<mpfr::mpreal> const& in,
        mpfr::mpreal& a, mpfr::mpreal& b,
        mp_prec_t prec)
{
    using mpfr::mpreal;
    mpfr_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    out.resize(in.size());
    for (std::size_t i{0u}; i < in.size(); ++i)
        out[i] = mpfr::fma((b - a) / 2, in[i], (b + a) / 2);

    mpreal::set_default_prec(prevPrec); 
}



void clenshaw(mpfr::mpreal &result, const std::vector<mpfr::mpreal> &p,
        const mpfr::mpreal &x, const mpfr::mpreal &a, const mpfr::mpreal &b,
        mp_prec_t prec)
{
    using mpfr::mpreal;
    mpfr_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    mpfr::mpreal bn1, bn2, bn;
    mpfr::mpreal buffer;

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

    mpreal::set_default_prec(prevPrec); 
}

void clenshaw(mpfr::mpreal &result, 
              const std::vector<mpfr::mpreal> &p,
              const mpfr::mpreal &x,
              chebkind_t kind,
              mp_prec_t prec)
{
    using mpfr::mpreal;
    mpfr_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    mpfr::mpreal bn1, bn2, bn;

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

    mpreal::set_default_prec(prevPrec); 
}

void equipts(std::vector<mpfr::mpreal>& v, std::size_t n,
             mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);
    mpreal pi = mpfr::const_pi(prec);

    v.resize(n);
    // store the points in the vector v as 
    // v[i] = i * pi / (n-1)
    for(std::size_t i{0u}; i < n; ++i) {
        v[i] = pi * i;
        v[i] /= (n-1);
    }

    mpreal::set_default_prec(prevPrec);
}


// this function computes the values of the coefficients of 
// the CI when Chebyshev nodes of the second kind are used
void chebcoeffs(std::vector<mpfr::mpreal>& c,
        std::vector<mpfr::mpreal>& fv, 
        mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    std::size_t n = fv.size();
    std::vector<mpfr::mpreal> v(n);
    equipts(v, n);

    mpfr::mpreal buffer;

    // halve the first and last coefficients
    mpfr::mpreal oldValue1 = fv[0];
    mpfr::mpreal oldValue2 = fv[n-1u];
    fv[0u] >>= 1;
    fv[n-1u] >>= 1;

    for(std::size_t i{0u}; i < n; ++i) {
        // compute the actual value at the Chebyshev
        // node cos(i * pi / n)
        buffer = mpfr::cos(v[i]);
        clenshaw(c[i], fv, buffer, FIRST, prec);

        if(i == 0u || i == n-1u) {
            c[i] /= (n-1u);
        } else {
            c[i] <<= 1;
            c[i] /= (n-1u);
        }
    }
    fv[0u] = oldValue1;
    fv[n-1u] = oldValue2;

    mpreal::set_default_prec(prevPrec);
}

// function that generates the coefficients of the 
// derivative of a given CI
void diffcoeffs(std::vector<mpfr::mpreal>& dc,
                std::vector<mpfr::mpreal>& c,
                chebkind_t kind,
                mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

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
            dc[0] >>= 1;
        };
        break;
        default: {
            int n = c.size() - 1;
            for(int i{n}; i > 0; --i)
                dc[i - 1] = c[i] * i;
        }
        break;
    }

    mpreal::set_default_prec(prevPrec);
}

void roots(std::vector<mpfr::mpreal>& r, std::vector<mpfr::mpreal>& c,
           std::pair<mpfr::mpreal, mpfr::mpreal> const &dom,
           chebkind_t kind,
           bool balance,
           mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    r.clear();
    for(auto &it : c)
        if(!mpfr::isfinite(it))
            return;

    MatrixXmp C = colleague(c, kind, balance, prec);
    Eigen::EigenSolver<MatrixXmp> es(C);
    VectorXcmp eigs = es.eigenvalues();

    mpfr::mpreal threshold = 1e-20;
    for(Eigen::Index i{0}; i < eigs.size(); ++i) {
        if(mpfr::fabs(eigs(i).imag()) < threshold)
            if(dom.first < eigs(i).real() && 
               dom.second >= eigs(i).real()) {
                r.push_back(eigs(i).real());
            }
    }

    std::sort(begin(r), end(r));

    mpreal::set_default_prec(prevPrec);
}