#ifndef SIA2_BANACH_H
#define SIA2_BANACH_H

#ifdef SIA2_STATIC
#define SIA2_API_BAN static
#else
#define SIA2_API_BAN extern
#endif

#include "common.h"

/*
 * Banach fixed-point theorem machinery for SIA².
 * Quantifies contraction and predicts convergence of improvement iterations.
 */

typedef struct {
    double ratio;         /* Lipschitz / contraction ratio q */
    int    is_contraction; /* q < 1 ? */
    int    predicted_iters; /* estimate of iterations to converge within eps */
    double residual;       /* ||Φ(x) - x|| at last evaluation */
} sia2_banach_result;

/*
 * Estimate contraction ratio q = sup ||Φ(x) - Φ(y)|| / ||x - y||
 * from samples.  xs, ys are pairs of points; phis_x, phis_y are their images.
 * count = number of pairs.
 */
SIA2_API_BAN double sia2_contraction_ratio(const double *xs, const double *ys,
                                           const double *phis_x, const double *phis_y,
                                           int count, int dim);

/*
 * Full Banach analysis: contraction ratio + convergence prediction.
 */
SIA2_API_BAN sia2_banach_result sia2_banach_analyze(const double *xs,
                                                     const double *ys,
                                                     const double *phis_x,
                                                     const double *phis_y,
                                                     int count, int dim,
                                                     double epsilon);

/*
 * Detect if we're at a fixed point: ||Φ(x) - x|| < tol.
 */
SIA2_API_BAN int sia2_is_fixed_point(const double *x, const double *phi_x,
                                     int dim, double tol);

/*
 * Predict iterations to converge: ceil(log(eps / d0) / log(q)).
 * d0 = initial distance, q = contraction ratio, eps = target precision.
 */
SIA2_API_BAN int sia2_predict_convergence(double q, double d0, double epsilon);

#ifdef SIA2_STATIC
#include "banach.c"
#endif
#endif /* SIA2_BANACH_H */
