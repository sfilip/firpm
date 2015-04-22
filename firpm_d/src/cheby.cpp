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

#include "firpm/cheby.h"

void applyCos(std::vector<double>& out,
        std::vector<double> const& in)
{
    for (std::size_t i{0u}; i < in.size(); ++i)
        out[i] = cosl(in[i]);
}

void changeOfVariable(std::vector<double>& out,
        std::vector<double> const& in,
        double& a, double& b)
{
    for (std::size_t i{0u}; i < in.size(); ++i)
        out[i] = (b + a) / 2 + in[i] * (b - a) / 2;
}

void evaluateClenshaw(double &result, std::vector<double> &p,
        double &x, double &a, double &b)
{
    double bn1, bn2, bn;
    double buffer;

    bn1 = 0;
    bn2 = 0;

    // compute the value of (2*x - b - a)/(b - a) in the temporary
    // variable buffer
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
    // (i.e. the CI value at x)
    result = bn1 - buffer * bn2;
}

void evaluateClenshaw(double &result, std::vector<double> &p,
                            double &x)
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

    result = x * bn1 - bn2 + p[0];
}

void evaluateClenshaw2ndKind(double &result, std::vector<double> &p,
                            double &x)
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

    result = (x * 2) * bn1 - bn2 + p[0];
}


void generateEquidistantNodes(std::vector<double>& v, std::size_t n)
{
    // store the points in the vector v as v[i] = i * pi / n
    for(std::size_t i{0u}; i <= n; ++i) {
        v[i] = M_PI * i;
        v[i] /= n;
    }
}

void generateChebyshevPoints(std::vector<double>& x, std::size_t n)
{
    // n is the number of points - 1
    x.reserve(n + 1u);
    if(n > 0u)
    {
        for(int k = n; k >= -n; k -= 2)
            x.push_back(sin(M_PI * k / (n * 2)));
    }
    else
    {
        x.push_back(0);
    }
}

// this function computes the values of the coefficients of the CI when
// Chebyshev nodes of the second kind are used
void generateChebyshevCoefficients(std::vector<double>& c,
        std::vector<double>& fv, std::size_t n)
{
    std::vector<double> v(n + 1);
    generateEquidistantNodes(v, n);

    double buffer;

    // halve the first and last coefficients
    double oldValue1 = fv[0];
    double oldValue2 = fv[n];
    fv[0] /= 2;
    fv[n] /= 2;

    for(std::size_t i{0u}; i <= n; ++i) {
        buffer = cos(v[i]);         // compute the actual value at the Chebyshev
                                    // node cos(i * pi / n)

        evaluateClenshaw(c[i], fv, buffer);

        if(i == 0u || i == n) {
            c[i] /= n;
        } else {
            c[i] *= 2;
            c[i] /= n;
        }
    }
    fv[0] = oldValue1;
    fv[n] = oldValue2;

}

// function that generates the coefficients of the derivative of a given CI
void derivativeCoefficients1stKind(std::vector<double>& derivC,
                                        std::vector<double>& c)
{
    int n = c.size() - 2;
    derivC[n] = c[n + 1] * (2 * (n + 1));
    derivC[n - 1] = c[n] * (2 * n);
    for(int i{n - 2}; i >= 0; --i) {
        derivC[i] = 2 * (i + 1);
        derivC[i] = derivC[i] * c[i + 1] + derivC[i + 2];
    }
    derivC[0] /= 2;
}

// use the formula (T_n(x))' = n * U_{n-1}(x)
void derivativeCoefficients2ndKind(std::vector<double>& derivC,
        std::vector<double>& c)
{
    std::size_t n = c.size() - 1;
    for(std::size_t i{n}; i > 0u; --i)
        derivC[i - 1] = c[i] * i;
}

