#include "fisher.h"
#include <math.h>
#include <string.h>

int sia2_fisher_compute(const double *score_vectors, int m, int dim,
                        sia2_fisher *F)
{
    int i, j, k;
    if (!score_vectors || !F || m < 1 || dim < 1 || dim > SIA2_MAX_DIM)
        return SIA2_ERR_NULL;

    F->dim = dim;
    memset(F->matrix, 0, dim * dim * sizeof(double));

    /* F = (1/m) Σ s_i s_i^T */
    for (k = 0; k < m; k++) {
        const double *s = score_vectors + k * dim;
        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                F->matrix[i * dim + j] += s[i] * s[j];
    }
    for (i = 0; i < dim * dim; i++)
        F->matrix[i] /= m;

    return SIA2_OK;
}

int sia2_fisher_diagonal(const double *score_vectors, int m, int dim,
                         double *diag)
{
    int k, i;
    if (!score_vectors || !diag || m < 1) return SIA2_ERR_NULL;
    memset(diag, 0, dim * sizeof(double));

    for (k = 0; k < m; k++) {
        const double *s = score_vectors + k * dim;
        for (i = 0; i < dim; i++)
            diag[i] += s[i] * s[i];
    }
    for (i = 0; i < dim; i++)
        diag[i] /= m;

    return SIA2_OK;
}

int sia2_natural_gradient(const sia2_fisher *F, const double *grad,
                          double damping, sia2_natural_grad *ng)
{
    int i;
    if (!F || !grad || !ng) return SIA2_ERR_NULL;
    int dim = F->dim;
    ng->dim = dim;

    /* Diagonal natural gradient:  Δθ_i = grad_i / (F_ii + λ) */
    for (i = 0; i < dim; i++) {
        double fii = F->matrix[i * dim + i] + damping;
        if (fabs(fii) < SIA2_EPSILON) {
            ng->direction[i] = 0.0;
        } else {
            ng->direction[i] = grad[i] / fii;
        }
    }

    /* Compute step size as L2 norm */
    double norm = 0.0;
    for (i = 0; i < dim; i++) norm += ng->direction[i] * ng->direction[i];
    ng->step_size = sqrt(norm);

    return SIA2_OK;
}

int sia2_fisher_effective_rank(const sia2_fisher *F, double threshold)
{
    int i, rank = 0;
    double max_eig = 0.0;
    int dim;
    if (!F) return 0;
    dim = F->dim;

    /* Find max diagonal element as proxy for max eigenvalue */
    for (i = 0; i < dim; i++) {
        double v = fabs(F->matrix[i * dim + i]);
        if (v > max_eig) max_eig = v;
    }

    if (max_eig < SIA2_EPSILON) return 0;
    for (i = 0; i < dim; i++) {
        if (fabs(F->matrix[i * dim + i]) > threshold * max_eig)
            rank++;
    }
    return rank;
}
