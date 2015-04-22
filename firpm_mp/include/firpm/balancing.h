/**
 * @file balancing.h
 * @author Silviu Filip
 * @date 12 March 2015
 * @brief Utility function for balancing a square matrix
 *
 */

#ifndef BALANCING_H_
#define BALANCING_H_

#include "util.h"

/*! Routine implementing the matrix balancing algorithm from B. N. Parlett and C. Reinsch (1969)
 * "Balancing a matrix for calculation of eigenvalues and eigenvectors" to be used as a
 * preprocessing step for performing eigenvalue computations.
 * @param[out] A the address of the matrix to be balanced (will be modified by the routine)
 * @param[in] n the dimension of the matrix
 */
void balance(mpfr::mpreal *A, int n);

#endif
