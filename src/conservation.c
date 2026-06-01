#include "conservation.h"
#include <math.h>
#include <string.h>

int sia2_check_capability(double cap_before, double cap_after, double threshold)
{
    if (fabs(cap_before) < SIA2_EPSILON) return (fabs(cap_after) < threshold) ? 1 : 0;
    double rel = fabs(cap_after - cap_before) / fabs(cap_before);
    return (rel < threshold) ? 1 : 0;
}

int sia2_check_entropy(double entropy_before, double entropy_after, double tolerance)
{
    return (entropy_after <= entropy_before + tolerance) ? 1 : 0;
}

int sia2_check_continuity(const double *f_before, const double *f_after,
                          int n, double threshold)
{
    int i;
    double diff_norm = 0.0, f_norm = 0.0;
    if (!f_before || !f_after || n < 1) return 0;
    for (i = 0; i < n; i++) {
        double d = f_after[i] - f_before[i];
        diff_norm += d * d;
        f_norm += f_before[i] * f_before[i];
    }
    f_norm = sqrt(f_norm);
    if (f_norm < SIA2_EPSILON) return (sqrt(diff_norm) < threshold) ? 1 : 0;
    return (sqrt(diff_norm) / f_norm < threshold) ? 1 : 0;
}

int sia2_check_monotonicity(const double *trajectory, int len)
{
    int i;
    if (!trajectory || len < 2) return 0;
    for (i = 1; i < len; i++) {
        if (trajectory[i] < trajectory[i - 1] - SIA2_EPSILON) return 0;
    }
    return 1;
}

sia2_conservation_result sia2_check_all_conservation(
    double cap_before, double cap_after,
    double ent_before, double ent_after,
    const double *f_before, const double *f_after, int fn_dim,
    const double *traj, int traj_len,
    double tol)
{
    sia2_conservation_result r;
    memset(&r, 0, sizeof(r));

    r.capability_before = cap_before;
    r.capability_after  = cap_after;
    r.entropy_before    = ent_before;
    r.entropy_after     = ent_after;

    int c_cap  = sia2_check_capability(cap_before, cap_after, tol);
    int c_ent  = sia2_check_entropy(ent_before, ent_after, tol);
    int c_cont = sia2_check_continuity(f_before, f_after, fn_dim, tol);
    int c_mono = sia2_check_monotonicity(traj, traj_len);

    if (traj && traj_len >= 2) {
        r.monotonicity_slope = (traj[traj_len - 1] - traj[0]) / (traj_len - 1);
    }
    if (f_before && f_after && fn_dim > 0) {
        double d = 0.0;
        int i;
        for (i = 0; i < fn_dim; i++) {
            double diff = f_after[i] - f_before[i];
            d += diff * diff;
        }
        r.continuity_delta = sqrt(d);
    }

    r.passed = c_cap && c_ent && c_cont && c_mono;
    return r;
}
