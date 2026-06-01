#include "banach.h"
#include <math.h>
#include <float.h>

static double vec_dist(const double *a, const double *b, int dim)
{
    double s = 0.0;
    int i;
    for (i = 0; i < dim; i++) { double d = a[i] - b[i]; s += d * d; }
    return sqrt(s);
}

double sia2_contraction_ratio(const double *xs, const double *ys,
                              const double *phis_x, const double *phis_y,
                              int count, int dim)
{
    double max_ratio = 0.0;
    int i;
    if (!xs || !ys || !phis_x || !phis_y || count < 1) return DBL_MAX;

    for (i = 0; i < count; i++) {
        double d_xy = vec_dist(xs + i * dim, ys + i * dim, dim);
        double d_ph = vec_dist(phis_x + i * dim, phis_y + i * dim, dim);
        if (d_xy > SIA2_EPSILON) {
            double r = d_ph / d_xy;
            if (r > max_ratio) max_ratio = r;
        }
    }
    return max_ratio;
}

int sia2_predict_convergence(double q, double d0, double epsilon)
{
    if (q >= 1.0 || q <= 0.0) return -1;  /* won't converge */
    if (d0 < epsilon) return 0;
    /* n ≥ log(ε/d₀) / log(q)  — note log(q) < 0, log(ε/d₀) < 0 */
    return (int)ceil(log(epsilon / d0) / log(q));
}

int sia2_is_fixed_point(const double *x, const double *phi_x, int dim, double tol)
{
    return vec_dist(x, phi_x, dim) < tol;
}

sia2_banach_result sia2_banach_analyze(const double *xs, const double *ys,
                                       const double *phis_x, const double *phis_y,
                                       int count, int dim, double epsilon)
{
    sia2_banach_result r;
    r.ratio = sia2_contraction_ratio(xs, ys, phis_x, phis_y, count, dim);
    r.is_contraction = (r.ratio < 1.0);
    r.residual = 0.0;

    if (count > 0) {
        /* Residual from last pair */
        r.residual = vec_dist(xs + (count - 1) * dim, phis_x + (count - 1) * dim, dim);
    }

    if (r.is_contraction && r.residual > epsilon) {
        r.predicted_iters = sia2_predict_convergence(r.ratio, r.residual, epsilon);
    } else {
        r.predicted_iters = 0;
    }
    return r;
}
