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
//    along with this program.  If not, see <http://www.gnu.org/licenses/>

#include "firpm/eigenvalue.h"
#include <algorithm>

void balance(MatrixXq& A)
{
    std::size_t n = A.rows();

    double rNorm;      // row norm
    double cNorm;      // column norm
    bool converged = false;

    double g, f, s;
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
                cNorm += fabs(A(j, i));
                rNorm += fabs(A(i, j));
            }
            if((cNorm == 0.0) || (rNorm == 0))
                continue;

            g = rNorm / 2.0;
            f = 1.0;
            s = cNorm + rNorm;


            while(std::isfinite(cNorm) && cNorm < g)
            {
                f *= 2.0;
                cNorm *= 4.0;
            }

            g = rNorm * 2.0;

            while(std::isfinite(cNorm) && cNorm > g)
            {
                f /= 2.0;
                cNorm /= 4.0;
            }

            if((rNorm + cNorm) < s * f * 0.95)
            {
                converged = false;
                g = 1.0 / f;
                // multiply by D^{-1} on the left
                A.row(i) *= g;
                // multiply by D on the right
                A.col(i) *= f;
            }

        }
    }
}


void generateColleagueMatrix1stKind(MatrixXq& C,
        std::vector<double>& a, bool withBalancing)
{
    std::vector<double> c = a;

    std::size_t n = a.size() - 1;
    // construct the initial matrix

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            C(i, j) = 0;


    double denom = -1;
    denom /= c[n];
    denom /= 2;
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
}

void generateColleagueMatrix2ndKind(MatrixXq& C,
        std::vector<double>& a, bool withBalancing)
{
    std::vector<double> c = a;

    std::size_t n = a.size() - 1;
    // construct the initial matrix

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            C(i, j) = 0;


    double denom = -1;
    denom /= c[n];
    denom /= 2;
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
}

void determineEigenvalues(VectorXcq &eigenvalues,
        MatrixXq &C)
{
    Eigen::EigenSolver<MatrixXq> es(C);
    eigenvalues = es.eigenvalues();
}


void getRealValues(std::vector<double> &realValues,
        VectorXcq &complexValues,
        double &a, double &b)
{
    double threshold = 10;
    threshold = pow(10, -20);
    for (int i = 0; i < complexValues.size(); ++i)
    {
        double imagValue = fabs(complexValues(i).imag());
        if(imagValue < threshold) {
            if(a <= complexValues(i).real() && b >= complexValues(i).real()) {
                realValues.push_back(complexValues(i).real());
            }
        }
    }
    std::sort(realValues.begin(), realValues.end());
}
