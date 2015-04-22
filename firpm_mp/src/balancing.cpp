#include "firpm/balancing.h"

void swapRows(mpfr::mpreal *A, int row1, int row2, int ncols)
{
    mpfr::mpreal buffer;
    mpfr::mpreal* row1A = A + ncols * row1;
    mpfr::mpreal* row2A = A + ncols * row2;

    for(int j = 0; j < ncols; ++j)
    {
        buffer = *row1A;
        *row1A++ = *row2A;
        *row2A++ = buffer;
    }
}

void swapCols(mpfr::mpreal *A, int col1, int col2,
        int nrows, int ncols)
{
    mpfr::mpreal buffer;
    mpfr::mpreal* col1A = A + col1;
    mpfr::mpreal* col2A = A + col2;
    for(int i = 0; i < nrows; col1A += ncols, col2A += ncols, ++i)
    {
        buffer = *col1A;
        *col1A = *col2A;
        *col2A = buffer;
    }
}

void multiplyRowByScalar(mpfr::mpreal *A, const mpfr::mpreal &x,
        int row, int ncols)
{
    mpfr::mpreal* rowA = A + ncols * row;
    for (int j = 0; j < ncols; ++j)
        *rowA++ *= x;
}

int searchRows(mpfr::mpreal *A, int p[], int n)
{
    mpfr::mpreal *pA;
    int i, j;
    int bottomBlock;

    // search for rows isolating an eigenvalue and push them down
    bottomBlock = n - 1;

    pA = A + bottomBlock * n;
    for(i = bottomBlock; i >= 0; --i, pA -= n) {
        for(j = 0; j < i; ++j)
            if(*(pA + j) != 0.0)
                break;
        if(j < i)
            continue;
        for(j = i + 1; j < n; ++j)
            if(*(pA + j) != 0.0)
                break;
        if(j < n)
            continue;
        swapRows(A, i, bottomBlock, n);
        swapCols(A, i, bottomBlock, n, n);
        p[bottomBlock] = i;
        --bottomBlock;
    }
    return bottomBlock;
}

int searchCols(mpfr::mpreal *A, int p[], int bottomBlock, int n)
{
    mpfr::mpreal *pA;
    int i,j;
    int topBlock;

    // search for columns isolating an eigenvalue and push them left
    topBlock = 0;
    pA = A + 1;
    for(i = 0; i < bottomBlock; ++i, pA = A + i) {
       for(j = 0; j < i; ++j, pA += n)
            if(*pA != 0.0)
                break;
        if(j < i)
            continue;
        pA += n;
        for(j = i + 1; j < n; ++j, pA += n)
            if (*pA != 0.0)
                break;
        if(j < n)
            continue;
        if(i > topBlock) {
            swapRows(A, i, topBlock, n);
            swapCols(A, i, topBlock, n, n);
            p[topBlock] = i;
            --i;
        }
        ++topBlock;
    }
    return topBlock;

}


mpfr::mpreal computeExponent(const mpfr::mpreal &x)
{
    mpfr::mpreal y = 1;
    int k = 0;

    k = mpfr::floor(mpfr::log2(x >> 1) >> 1).toLong();
    y <<= k;
    while(x <= ((y * y) >> 1))
        y >>= 1;
    while(x > ((y * y) << 1))
        y <<= 1;
    return y;

}


void balance(mpfr::mpreal *A, int n)
{
    int topBlock, bottomBlock, converged;
    mpfr::mpreal *pA, *pC, *pTop;
    mpfr::mpreal colNorm, rowNorm;
    mpfr::mpreal x, y, z;
    int *p = new int[n + 2];

    // keep track of the row/column permutations
    for(int i = 0; i < n; ++i)
        p[i] = i;

    // search for rows isolating an eigenvalue and push them down
    bottomBlock = searchRows(A, p, n);
    p[n + 1] = bottomBlock;

    // search for columns isolating an eigenvalue and push them left
    topBlock = searchCols(A, p, bottomBlock, n);
    p[n] = topBlock;

    // upper triangular matrix; no need to perform balancing
    if (topBlock >= bottomBlock)
    {
        delete[] p;
        return;
    }

    // balance the submatrix block
    converged = 0;
    while(!converged) {
        converged = 1;
        pTop = A + topBlock * n;
        pA = pTop;
        for (int i = topBlock; i <= bottomBlock; ++i, pA += n) {
            colNorm = 0.0;
            rowNorm = 0.0;
            pC = pTop + i;
            for(int j = topBlock; j <= bottomBlock; ++j, pC += n) {
                if(j == i)
                    continue;
                colNorm += mpfr::abs(*pC);
                rowNorm += mpfr::abs(*(pA + j));
            }
            x = rowNorm / colNorm;
            y = computeExponent(x);
            z = mpfr::mpreal(1.0) / y;
            if((colNorm * y + rowNorm * z) < (colNorm + rowNorm) * 0.95) {
                pC = A + i;
                for(int j = 0; j <= bottomBlock; ++j, pC += n)
                    *pC *= y;
                for(int j = topBlock; j < n; ++j)
                    *(pA + j) *= z;
                converged = 0;
            }
        }
    }
    delete[] p;
}
