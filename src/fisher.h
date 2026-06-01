#ifndef SIA2_FISHER_H
#define SIA2_FISHER_H

#ifdef SIA2_STATIC
#define SIA2_API_FIS static
#else
#define SIA2_API_FIS extern
#endif

#include "common.h"

/*
 * Fisher information matrix and natural gradient for SIA².
 * Enables curvature-aware improvement steps.
 */

typedef struct {
    double matrix[SIA2_MAX_DIM * SIA2_MAX_DIM];
    int    dim;
} sia2_fisher;

typedef struct {
    double direction[SIA2_MAX_DIM];
    double step_size;
    int    dim;
} sia2_natural_grad;

/*
 * Compute Fisher information matrix from score vectors.
 * score_vectors: m × dim matrix (m samples, dim parameters).
 * F = (1/m) Σ ∇log p · ∇log p^T
 */
SIA2_API_FIS int sia2_fisher_compute(const double *score_vectors, int m, int dim,
                                     sia2_fisher *F);

/*
 * Approximate natural gradient:  Δθ ≈ F^(-1) ∇L
 * Uses diagonal approximation for stability (inverts only the diagonal).
 * grad: dimension dim.  Result stored in ng->direction.
 */
SIA2_API_FIS int sia2_natural_gradient(const sia2_fisher *F,
                                       const double *grad,
                                       double damping, sia2_natural_grad *ng);

/*
 * Fisher diagonal only — lighter weight.
 */
SIA2_API_FIS int sia2_fisher_diagonal(const double *score_vectors, int m, int dim,
                                      double *diag);

/*
 * Compute effective rank: number of eigenvalues > threshold * max_eigenvalue.
 */
SIA2_API_FIS int sia2_fisher_effective_rank(const sia2_fisher *F, double threshold);

#ifdef SIA2_STATIC
#include "fisher.c"
#endif
#endif /* SIA2_FISHER_H */
