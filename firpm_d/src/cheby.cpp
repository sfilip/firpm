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

/** Eigen matrix container for double values */
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
/** Eigen vector container for double values */
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> VectorXcd;

void balance(MatrixXd &A)
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

MatrixXd colleague(std::vector<double> const &c,
                   chebkind_t kind,
                   bool bal)
{
    std::vector<double> a = c;
    std::size_t n = c.size() - 1;
    MatrixXd C(n, n);

    for(std::size_t i{0u}; i < n; ++i)
        for(std::size_t j{0u}; j < n; ++j)
            C(i, j) = 0;

    double denom = -1;
    denom /= a[n];
    denom /= 2;
    for(std::size_t i{0u}; i < c.size() - 1u; ++i)
        a[i] *= denom;
    a[n - 2u] += 0.5;

    for (std::size_t i = 0u; i < n - 1; ++i)
        C(i, i + 1) = C(i + 1, i) = 0.5;
    switch(kind) {
        case FIRST: C(n - 2, n - 1) = 1;    break;
        default:    C(n - 2, n - 1) = 0.5;  break;
    }
    for(std::size_t i{0u}; i < n; ++i)
        C(i, 0) = a[n - i - 1];

    if(bal)
        balance(C);    

    return C;
}

void cos(std::vector<double>& out,
        std::vector<double> const& in)
{
    out.resize(in.size());
    for(std::size_t i{0u}; i < in.size(); ++i)
        out[i] = cosl(in[i]);
}

void chgvar(std::vector<double>& out,
        std::vector<double> const& in,
        double& a, double& b)
{
    out.resize(in.size());
    for(std::size_t i{0u}; i < in.size(); ++i)
        out[i] = (b + a) / 2 + in[i] * (b - a) / 2;
}

void clenshaw(double &result, const std::vector<double> &p,
        const double &x, const double &a, const double &b)
{
    double bn1, bn2, bn;
    double buffer;

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

void clenshaw(double &result, 
              const std::vector<double> &p,
              const double &x,
              chebkind_t kind)
{
    double bn1, bn2, bn;

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

void equipts(std::vector<double>& v, std::size_t n)
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
void chebcoeffs(std::vector<double>& c,
        std::vector<double>& fv)
{
    std::size_t n = fv.size();
    std::vector<double> v(n);
    equipts(v, n);

    double buffer;

    // halve the first and last coefficients
    double oldValue1 = fv[0];
    double oldValue2 = fv[n-1u];
    fv[0] /= 2;
    fv[n-1u] /= 2;

    for(std::size_t i{0u}; i < n; ++i) {
        // compute the actual value at the Chebyshev
        // node cos(i * pi / n)
        buffer = cosl(v[i]);
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
void diffcoeffs(std::vector<double>& dc,
                std::vector<double>& c,
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
            for(int i{n}; i > 0u; --i)
            dc[i - 1] = c[i] * i;
        }
        break;
    }
}

void roots(std::vector<double>& r, std::vector<double>& c,
           std::pair<double, double> const &dom,
           chebkind_t kind,
           bool balance)
{
    r.clear();
    for(auto &it : c)
        if(!std::isfinite(it))
            return;

    MatrixXd C = colleague(c, kind, balance);
    Eigen::EigenSolver<MatrixXd> es(C);
    VectorXcd eigs = es.eigenvalues();

    double threshold = 1e-20;
    for(std::size_t i{0u}; i < eigs.size(); ++i) {
        if(fabs(eigs(i).imag()) < threshold)
            if(dom.first < eigs(i).real() && 
               dom.second >= eigs(i).real()) {
                r.push_back(eigs(i).real());
            }
    }

    std::sort(begin(r), end(r));
}