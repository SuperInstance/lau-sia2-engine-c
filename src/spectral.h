#ifndef SIA2_SPECTRAL_H
#define SIA2_SPECTRAL_H

#ifdef SIA2_STATIC
#define SIA2_API static
#else
#define SIA2_API extern
#endif

#include "common.h"
#include "arena.h"
#include <math.h>
#include <string.h>

/*
 * Power iteration eigendecomposition.
 * Finds the top k eigenvalues/eigenvectors of a symmetric n×n matrix.
 */

typedef struct {
    double value;
    double vector[SIA2_MAX_DIM];
} sia2_eigenmode;

typedef struct {
    sia2_eigenmode modes[SIA2_MAX_DIM];
    int n_modes;
} sia2_spectrum;

/*
 * Power iteration for dominant eigenpair.
 * A is n×n row-major symmetric.  Returns SIA2_OK or SIA2_ERR_CONVERGE.
 */
SIA2_API int sia2_power_iteration(const double *A, int n,
                                  double *eigenvalue, double *eigenvector,
                                  int max_iters, double tol);

/*
 * Deflated eigendecomposition — top k modes via repeated power iteration
 * with Hotelling deflation.
 */
SIA2_API int sia2_eigendecompose(const double *A, int n, int k,
                                 sia2_spectrum *spec,
                                 int max_iters, double tol);

/*
 * Rank eigenmodes by magnitude (descending).
 */
SIA2_API void sia2_rank_eigenmodes(sia2_spectrum *spec);

/*
 * Find the weakest eigenmode (smallest |eigenvalue|).
 */
SIA2_API int sia2_weakest_mode(const sia2_spectrum *spec);

/*
 * Compute spectral gap: |λ_1| - |λ_2| for top two modes.
 */
SIA2_API double sia2_spectral_gap(const sia2_spectrum *spec);

#ifdef SIA2_STATIC
#include "spectral.c"
#endif

#endif /* SIA2_SPECTRAL_H */
