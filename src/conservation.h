#ifndef SIA2_CONSERVATION_H
#define SIA2_CONSERVATION_H

#ifdef SIA2_STATIC
#define SIA2_API_CONS static
#else
#define SIA2_API_CONS extern
#endif

#include "common.h"
#include <math.h>

/*
 * Conservation law checking for SIA² improvement trajectories.
 * Ensures improvement doesn't violate structural invariants.
 */

typedef struct {
    int    passed;
    double capability_before;
    double capability_after;
    double entropy_before;
    double entropy_after;
    double continuity_delta;
    double monotonicity_slope;
} sia2_conservation_result;

/*
 * Check capability conservation: |ΔC| / C < threshold.
 * Improvement should not destroy existing capabilities.
 */
SIA2_API_CONS int sia2_check_capability(double cap_before, double cap_after,
                                        double threshold);

/*
 * Check entropy bound: H(after) ≤ H(before) + tolerance.
 * Prevents runaway entropy increase during improvement.
 */
SIA2_API_CONS int sia2_check_entropy(double entropy_before, double entropy_after,
                                     double tolerance);

/*
 * Check continuity: ||f(x+ε) - f(x)|| / ||f(x)|| < threshold.
 */
SIA2_API_CONS int sia2_check_continuity(const double *f_before,
                                        const double *f_after,
                                        int n, double threshold);

/*
 * Check monotonicity: improvement should be non-decreasing in capability.
 */
SIA2_API_CONS int sia2_check_monotonicity(const double *trajectory, int len);

/*
 * Run all conservation checks at once.
 */
SIA2_API_CONS sia2_conservation_result sia2_check_all_conservation(
    double cap_before, double cap_after,
    double ent_before, double ent_after,
    const double *f_before, const double *f_after, int fn_dim,
    const double *traj, int traj_len,
    double tol);

#ifdef SIA2_STATIC
#include "conservation.c"
#endif
#endif /* SIA2_CONSERVATION_H */
