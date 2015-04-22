#include "firpm/balancing.h"

void swapRows(double *A, int row1, int row2, int ncols)
{
    double buffer;
    double* row1A = A + ncols * row1;
    double* row2A = A + ncols * row2;

    for(int j = 0; j < ncols; ++j)
    {
        buffer = *row1A;
        *row1A++ = *row2A;
        *row2A++ = buffer;
    }
}

void swapCols(double *A, int col1, int col2,
        int nrows, int ncols)
{
    double buffer;
    double* col1A = A + col1;
    double* col2A = A + col2;
    for(int i = 0; i < nrows; col1A += ncols, col2A += ncols, ++i)
    {
        buffer = *col1A;
        *col1A = *col2A;
        *col2A = buffer;
    }
}

void multiplyRowByScalar(double *A, const double &x,
        int row, int ncols)
{
    double* rowA = A + ncols * row;
    for (int j = 0; j < ncols; ++j)
        *rowA++ *= x;
}

int searchRows(double *A, int p[], int n)
{
    double *pA;
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

int searchCols(double *A, int p[], int bottomBlock, int n)
{
    double *pA;
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


double computeExponent(const double &z)
{
    union { double dbl; short word[4]; } x, y;
    double y2;
    const short exp_mask = 0x7FF0;                  // Exponent mask
    const short bias     = 0x3FF0;                  // Exponent of 1.0
    short offset;

    x.dbl = z;
    y.dbl = 1.0;
    offset = x.word[3] - bias + 0x0010;
    if (offset >= 0) y.word[3] = ((offset >> 1) & exp_mask) + bias;
    else y.word[3] = bias - (((-offset) >> 1) & exp_mask);
    y2 = y.dbl * y.dbl;
    if ( x.dbl <= (y2 * 0.5L ) ) y.word[3] -= 0x0010;
    if ( x.dbl > (2.0L * y2 ) ) y.word[3] += 0x0010;

    return y.dbl;
}

/*double computeExponent(const double &x)
{
   double y;
   int k;

   k = log(0.5 * x) / log(4.0);
   y = pow(2.0, k);
   while ( x <= ( y * y * 0.5) ) y /= 2.0;
   while ( x > (2.0 * y * y) ) y *= 2.0;
   return y;
}*/



void balance(double *A, int n)
{
    int topBlock, bottomBlock, converged;
    double *pA, *pC, *pTop;
    double colNorm, rowNorm;
    double x, y, z;
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
            z = double(1.0) / y;
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
