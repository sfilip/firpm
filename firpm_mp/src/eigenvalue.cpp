//    firpm_mp
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
//    along with this program.  If not, see <http://www.gnu.org/licenses/>



#include "firpm/eigenvalue.h"
#include <algorithm>

void balance(MatrixXq& A)
{
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
}


void generateColleagueMatrix1stKind(MatrixXq& C,
        std::vector<mpfr::mpreal>& a,
        bool withBalancing,
        mp_prec_t prec)
{
    using mpfr::mpreal;
    mpfr_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);
    std::vector<mpfr::mpreal> c = a;

    std::size_t n = a.size() - 1;
    // construct the initial matrix

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            C(i, j) = 0;


    mpreal denom = -1;
    denom /= c[n];
    denom >>= 1;
    for(std::size_t i = 0u; i < a.size() - 1; ++i)
        c[i] *= denom;
    c[n - 2] += 0.5;

    for (std::size_t i = 0u; i < n - 1; ++i)
        C(i, i + 1) = C(i + 1, i) = 0.5;
    C(n - 2, n - 1) = 1;

    for(std::size_t i = 0u; i < n; ++i)
    {
        C(i, 0) = c[n - i - 1];
    }

    if(withBalancing)
        balance(C);
    mpreal::set_default_prec(prevPrec);

}

void determineEigenvalues(VectorXcq &eigenvalues,
        MatrixXq &C)
{
    Eigen::EigenSolver<MatrixXq> es(C);
    eigenvalues = es.eigenvalues();
}


void getRealValues(std::vector<mpfr::mpreal> &roots,
        VectorXcq &eigenValues,
        mpfr::mpreal &a, mpfr::mpreal &b)
{
    using mpfr::mpreal;
    mpreal threshold = 10;
    mpfr_pow_si(threshold.mpfr_ptr(),
            threshold.mpfr_srcptr(), -20, GMP_RNDN);
    for (int i = 0; i < eigenValues.size(); ++i)
    {
        mpreal imagValue = mpfr::abs(eigenValues(i).imag());
        if(mpfr::abs(eigenValues(i).imag()) < threshold) {
            if(a <= eigenValues(i).real() && b >= eigenValues(i).real()) {
                roots.push_back(eigenValues(i).real());
            }
        }
    }
    std::sort(roots.begin(), roots.end());
}



void generateColleagueMatrix2ndKind(MatrixXq& C,
        std::vector<mpfr::mpreal>& a,
        bool withBalancing,
        mp_prec_t prec)
{
    using mpfr::mpreal;
    mpfr_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);
    std::vector<mpfr::mpreal> c = a;

    std::size_t n = a.size() - 1;
    // construct the initial matrix

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            C(i, j) = 0;


    mpreal denom = -1;
    denom /= c[n];
    denom >>= 1;
    for(std::size_t i = 0u; i < a.size() - 1; ++i)
        c[i] *= denom;
    c[n - 2] += 0.5;

    for (std::size_t i = 0u; i < n - 1; ++i)
        C(i, i + 1) = C(i + 1, i) = 0.5;
    C(n - 2, n - 1) = 0.5;

    for(std::size_t i = 0u; i < n; ++i)
    {
        C(i, 0) = c[n - i - 1];
    }

    if(withBalancing)
        balance(C);
    mpreal::set_default_prec(prevPrec);

}
