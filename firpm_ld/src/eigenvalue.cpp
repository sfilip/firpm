#include "firpm/eigenvalue.h"
#include "firpm/balancing.h"
#include <algorithm>


void generateColleagueMatrix1stKind(MatrixXq& C,
        std::vector<long double>& a)
{
    std::vector<long double> c = a;

    std::size_t n = a.size() - 1;
    // construct the initial matrix

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            C(i, j) = 0;


    long double denom = -1;
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
}

void generateColleagueMatrix1stKindWithBalancing(MatrixXq& C,
        std::vector<long double>& a)
{
    std::size_t n = a.size() - 1;


    long double* A = new long double[n * n];

    std::vector<long double> c = a;

    // construct the initial matrix

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            A[n * i + j] = 0;


    long double denom = -1;
    denom /= c[n];
    denom /= 2;
    for(std::size_t i = 0u; i < a.size() - 1; ++i)
        c[i] *= denom;
    c[n - 2] += 0.5;

    for (std::size_t i = 0u; i < n - 1; ++i)
        A[n * i + i + 1] = A[(i + 1) * n + i] = 0.5;
    A[n * (n - 1) - 1] = 1;

    for(std::size_t i = 0u; i < n; ++i)
    {
        A[n * i] = c[n - i - 1];
    }


    balance(A, n);

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            C(i, j) = A[i * n + j];

    delete[] A;
}


void determineEigenvalues(VectorXcq &eigenvalues,
        MatrixXq &C)
{
    Eigen::EigenSolver<MatrixXq> es(C);
    eigenvalues = es.eigenvalues();
}


void getRealValues(std::vector<long double> &realValues,
        VectorXcq &complexValues,
        long double &a, long double &b)
{
    long double threshold = 10;
    threshold = pow(10, -20);
    for (int i = 0; i < complexValues.size(); ++i)
    {
        long double imagValue = fabs(complexValues(i).imag());
        if(imagValue < threshold) {
            if(a <= complexValues(i).real() && b >= complexValues(i).real()) {
                realValues.push_back(complexValues(i).real());
            }
        }
    }
    std::sort(realValues.begin(), realValues.end());
}


void generateColleagueMatrix2ndKind(MatrixXq& C,
        std::vector<long double>& a)
{
    std::vector<long double> c = a;

    std::size_t n = a.size() - 1;
    // construct the initial matrix

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            C(i, j) = 0;


    long double denom = -1;
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
}



void generateColleagueMatrix2ndKindWithBalancing(MatrixXq& C,
        std::vector<long double>& a)
{
    std::size_t n = a.size() - 1;


    long double* A = new long double[n * n];

    std::vector<long double> c = a;

    // construct the initial matrix

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            A[n * i + j] = 0;


    long double denom = -1;
    denom /= c[n];
    denom /= 2;
    for(std::size_t i = 0u; i < a.size() - 1; ++i)
        c[i] *= denom;
    c[n - 2] += 0.5;

    for (std::size_t i = 0u; i < n - 1; ++i)
        A[n * i + i + 1] = A[(i + 1) * n + i] = 0.5;
    A[n * (n - 1) - 1] = 0.5;

    for(std::size_t i = 0u; i < n; ++i)
    {
        A[n * i] = c[n - i - 1];
    }

    balance(A, n);

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            C(i, j) = A[i * n + j];

    delete[] A;
}

