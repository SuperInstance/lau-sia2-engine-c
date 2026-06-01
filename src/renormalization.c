#include "renormalization.h"
#include <math.h>
#include <string.h>

sia2_rg_class sia2_rg_classify(double beta, double threshold)
{
    double ab = fabs(beta);
    if (ab < threshold * 0.1) return SIA2_RG_TRIVIAL;
    if (ab < threshold)       return SIA2_RG_GAUSSIAN;
    if (ab < 1.0)             return SIA2_RG_NONTRIVIAL;
    return SIA2_RG_UNSTABLE;
}

int sia2_rg_detect_fixed_point(const double *g, int len, double threshold)
{
    int i;
    if (!g || len < 2) return 0;
    for (i = 1; i < len; i++) {
        double dg = g[i] - g[i - 1];
        if (fabs(dg) > threshold) return 0;
    }
    return 1;
}

sia2_rg_result sia2_rg_beta(const double *g, const double *scales,
                            int len, double threshold)
{
    sia2_rg_result r;
    memset(&r, 0, sizeof(r));

    if (!g || !scales || len < 2) {
        r.universality_class = SIA2_RG_TRIVIAL;
        return r;
    }

    /* Central difference for β at midpoint */
    int mid = len / 2;
    if (mid < 1) mid = 1;
    if (mid >= len - 1) mid = len - 2;

    double ds = scales[mid + 1] - scales[mid - 1];
    double dg = g[mid + 1] - g[mid - 1];

    if (fabs(ds) > SIA2_EPSILON)
        r.beta = dg / ds;
    else
        r.beta = 0.0;

    /* Anomalous dimension: γ ≈ dβ/dg at the midpoint */
    if (fabs(dg) > SIA2_EPSILON) {
        /* Approximate from data spread */
        double range = 0.0;
        int i;
        for (i = 0; i < len - 1; i++) range += fabs(g[i + 1] - g[i]);
        r.gamma = (len > 1) ? range / (len - 1) : 0.0;
    }

    r.is_fixed_point = (fabs(r.beta) < threshold);
    r.universality_class = sia2_rg_classify(r.beta, threshold);

    /* Correlation length exponent */
    if (fabs(r.beta) > threshold) {
        r.nu = 1.0 / fabs(r.beta);
    } else {
        r.nu = DBL_MAX; /* diverges at fixed point */
    }

    return r;
}

double sia2_rg_nu_exponent(const double *beta_values, int len)
{
    double sum = 0.0;
    int count = 0, i;
    if (!beta_values || len < 2) return 0.0;

    for (i = 0; i < len; i++) {
        double ab = fabs(beta_values[i]);
        if (ab > SIA2_EPSILON) {
            sum += 1.0 / ab;
            count++;
        }
    }
    return (count > 0) ? sum / count : 0.0;
}
