#ifndef SIA2_TRAJECTORY_H
#define SIA2_TRAJECTORY_H

#ifdef SIA2_STATIC
#define SIA2_API_TRAJ static
#else
#define SIA2_API_TRAJ extern
#endif

#include "common.h"

/*
 * Improvement trajectory tracking for SIA².
 * Records and analyzes the evolution of system state across improvement steps.
 */

#define SIA2_TRAJ_MAX_POINTS 256

typedef struct {
    double timestamp;
    double capability;
    double entropy;
    double residual;
    int    step;
} sia2_traj_point;

typedef struct {
    sia2_traj_point points[SIA2_TRAJ_MAX_POINTS];
    int count;
} sia2_trajectory;

/*
 * Initialize a trajectory.
 */
SIA2_API_TRAJ void sia2_trajectory_init(sia2_trajectory *t);

/*
 * Record a point in the trajectory.
 */
SIA2_API_TRAJ int sia2_trajectory_record(sia2_trajectory *t,
                                         double timestamp, double capability,
                                         double entropy, double residual,
                                         int step);

/*
 * Compute convergence rate from last n points.
 * Returns slope of log(residual) vs step.
 */
SIA2_API_TRAJ double sia2_trajectory_convergence_rate(const sia2_trajectory *t,
                                                      int last_n);

/*
 * Detect plateau: residual hasn't decreased significantly in last n points.
 */
SIA2_API_TRAJ int sia2_trajectory_plateau(const sia2_trajectory *t,
                                          int last_n, double threshold);

/*
 * Get the best (lowest residual) point.
 */
SIA2_API_TRAJ int sia2_trajectory_best(const sia2_trajectory *t);

/*
 * Compute total improvement: first to last capability delta.
 */
SIA2_API_TRAJ double sia2_trajectory_total_improvement(const sia2_trajectory *t);

#ifdef SIA2_STATIC
#include "trajectory.c"
#endif
#endif /* SIA2_TRAJECTORY_H */
