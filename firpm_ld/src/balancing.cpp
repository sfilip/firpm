#include "firpm/balancing.h"

void swapRows(long double *A, int row1, int row2, int ncols)
{
    long double buffer;
    long double* row1A = A + ncols * row1;
    long double* row2A = A + ncols * row2;

    for(int j = 0; j < ncols; ++j)
    {
        buffer = *row1A;
        *row1A++ = *row2A;
        *row2A++ = buffer;
    }
}

void swapCols(long double *A, int col1, int col2,
        int nrows, int ncols)
{
    long double buffer;
    long double* col1A = A + col1;
    long double* col2A = A + col2;
    for(int i = 0; i < nrows; col1A += ncols, col2A += ncols, ++i)
    {
        buffer = *col1A;
        *col1A = *col2A;
        *col2A = buffer;
    }
}

void multiplyRowByScalar(long double *A, const long double &x,
        int row, int ncols)
{
    long double* rowA = A + ncols * row;
    for (int j = 0; j < ncols; ++j)
        *rowA++ *= x;
}

int searchRows(long double *A, int p[], int n)
{
    long double *pA;
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

int searchCols(long double *A, int p[], int bottomBlock, int n)
{
    long double *pA;
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

/*long double computeExponent(const long double &x)
{
    long double y = 1;
    int k = 0;

    k = floor(log2(x * 0.5L) / 2);
    y = powl(2.0L, k);
    while ( x <= ( y * y * 0.5L ) ) y /= 2.0L;
    while ( x > ((y * y) * 2.0L) ) y *= 2.0L;
    return y;
}*/

long double computeExponent(const long double &z)
{
    union { long double dbl; short word[5]; } x, y;
    long double y2;
    const short exp_mask = 0x7FF0;                  // Exponent mask
    const short bias     = 0x3FF0;                  // Exponent of 1.0
    short offset;

    x.dbl = z;
    y.dbl = 1.0;
    offset = x.word[4] - bias + 0x0010;
    if (offset >= 0) y.word[4] = ((offset >> 1) & exp_mask) + bias;
    else y.word[4] = bias - (((-offset) >> 1) & exp_mask);
    y2 = y.dbl * y.dbl;
    if ( x.dbl <= (y2 * 0.5L ) ) y.word[4] -= 0x0010;
    if ( x.dbl > (2.0L * y2 ) ) y.word[4] += 0x0010;

    return y.dbl;
}




void balance(long double *A, int n)
{
    int topBlock, bottomBlock, converged;
    long double *pA, *pC, *pTop;
    long double colNorm, rowNorm;
    long double x, y, z;
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
                colNorm += fabs(*pC);
                rowNorm += fabs(*(pA + j));
            }
            x = rowNorm / colNorm;
            y = computeExponent(x);
            z = 1.0L / y;
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
