#include "spectral.h"
#include <math.h>
#include <string.h>

/* Normalize vector in-place to unit length.  Returns norm before normalization. */
static double normalize(double *v, int n)
{
    double norm = 0.0;
    int i;
    for (i = 0; i < n; i++) norm += v[i] * v[i];
    norm = sqrt(norm);
    if (norm < SIA2_EPSILON) return 0.0;
    for (i = 0; i < n; i++) v[i] /= norm;
    return norm;
}

/* Matrix-vector multiply:  y = A * x.  A is n×n row-major. */
static void matvec(const double *A, const double *x, double *y, int n)
{
    int i, j;
    for (i = 0; i < n; i++) {
        double s = 0.0;
        for (j = 0; j < n; j++) s += A[i * n + j] * x[j];
        y[i] = s;
    }
}

/* Dot product of two n-vectors. */
static double dot(const double *a, const double *b, int n)
{
    double s = 0.0;
    int i;
    for (i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

int sia2_power_iteration(const double *A, int n,
                         double *eigenvalue, double *eigenvector,
                         int max_iters, double tol)
{
    int i, iter;
    double v[SIA2_MAX_DIM];
    double w[SIA2_MAX_DIM];
    double prev_lambda;

    if (!A || !eigenvalue || !eigenvector || n < 1 || n > SIA2_MAX_DIM)
        return SIA2_ERR_NULL;

    /* Initial guess: unit vector [1,0,...] */
    memset(v, 0, sizeof(v));
    v[0] = 1.0;

    prev_lambda = 0.0;
    for (iter = 0; iter < max_iters; iter++) {
        matvec(A, v, w, n);
        double norm = normalize(w, n);
        if (norm < SIA2_EPSILON) {
            *eigenvalue = 0.0;
            memcpy(eigenvector, v, n * sizeof(double));
            return SIA2_OK;
        }
        *eigenvalue = dot(v, w, n) > 0 ? norm : -norm;

        /* Check convergence */
        if (iter > 0 && fabs(*eigenvalue - prev_lambda) < tol * (fabs(*eigenvalue) + SIA2_EPSILON)) {
            memcpy(eigenvector, w, n * sizeof(double));
            return SIA2_OK;
        }
        prev_lambda = *eigenvalue;
        memcpy(v, w, n * sizeof(double));
    }
    memcpy(eigenvector, v, n * sizeof(double));
    return SIA2_ERR_CONVERGE;
}

int sia2_eigendecompose(const double *A, int n, int k,
                        sia2_spectrum *spec,
                        int max_iters, double tol)
{
    double Adef[SIA2_MAX_DIM * SIA2_MAX_DIM];
    double ev, vec[SIA2_MAX_DIM];
    int m, rc;

    if (!A || !spec || n < 1 || k < 1 || k > n) return SIA2_ERR_NULL;
    if (n > SIA2_MAX_DIM) return SIA2_ERR_SIZE;

    memcpy(Adef, A, n * n * sizeof(double));
    spec->n_modes = 0;

    for (m = 0; m < k; m++) {
        rc = sia2_power_iteration(Adef, n, &ev, vec, max_iters, tol);
        if (rc != SIA2_OK && rc != SIA2_ERR_CONVERGE) return rc;

        spec->modes[m].value = ev;
        memcpy(spec->modes[m].vector, vec, n * sizeof(double));
        spec->n_modes = m + 1;

        /* Hotelling deflation: A' = A - λ v v^T */
        if (fabs(ev) > SIA2_EPSILON) {
            int i, j;
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    Adef[i * n + j] -= ev * vec[i] * vec[j];
        }
    }
    return SIA2_OK;
}

void sia2_rank_eigenmodes(sia2_spectrum *spec)
{
    int i, j;
    if (!spec || spec->n_modes < 2) return;
    /* Simple insertion sort descending by |eigenvalue| */
    for (i = 1; i < spec->n_modes; i++) {
        sia2_eigenmode tmp = spec->modes[i];
        double key = fabs(tmp.value);
        j = i - 1;
        while (j >= 0 && fabs(spec->modes[j].value) < key) {
            spec->modes[j + 1] = spec->modes[j];
            j--;
        }
        spec->modes[j + 1] = tmp;
    }
}

int sia2_weakest_mode(const sia2_spectrum *spec)
{
    int i, weakest;
    double min_mag;
    if (!spec || spec->n_modes == 0) return -1;

    weakest = 0;
    min_mag = fabs(spec->modes[0].value);
    for (i = 1; i < spec->n_modes; i++) {
        double mag = fabs(spec->modes[i].value);
        if (mag < min_mag) { min_mag = mag; weakest = i; }
    }
    return weakest;
}

double sia2_spectral_gap(const sia2_spectrum *spec)
{
    if (!spec || spec->n_modes < 2) return 0.0;
    return fabs(spec->modes[0].value) - fabs(spec->modes[1].value);
}
