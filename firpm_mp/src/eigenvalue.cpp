#include "firpm/eigenvalue.h"
#include "firpm/balancing.h"
#include <algorithm>


void generateColleagueMatrix1stKind(MatrixXq& C,
        std::vector<mpfr::mpreal>& a, mp_prec_t prec)
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
    mpreal::set_default_prec(prevPrec);

}

void generateColleagueMatrix1stKindWithBalancing(MatrixXq& C,
        std::vector<mpfr::mpreal>& a, mp_prec_t prec)
{
    using mpfr::mpreal;
    mpfr_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);
    std::size_t n = a.size() - 1;


    mpfr::mpreal* A = new mpfr::mpreal[n * n];

    std::vector<mpfr::mpreal> c = a;

    // construct the initial matrix

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            A[n * i + j] = 0;


    mpreal denom = -1;
    denom /= c[n];
    denom >>= 1;
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
        std::vector<mpfr::mpreal>& a, mp_prec_t prec)
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
    mpreal::set_default_prec(prevPrec);

}



void generateColleagueMatrix2ndKindWithBalancing(MatrixXq& C,
        std::vector<mpfr::mpreal>& a, mp_prec_t prec)
{
    using mpfr::mpreal;
    mpfr_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);
    std::size_t n = a.size() - 1;


    mpfr::mpreal* A = new mpfr::mpreal[n * n];

    std::vector<mpfr::mpreal> c = a;

    // construct the initial matrix

    for (std::size_t i = 0u; i < n; ++i)
        for (std::size_t j = 0u; j < n; ++j)
            A[n * i + j] = 0;


    mpreal denom = -1;
    denom /= c[n];
    denom >>= 1;
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
    mpreal::set_default_prec(prevPrec);
}
