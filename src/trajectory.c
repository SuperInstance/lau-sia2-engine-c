#include "trajectory.h"
#include <math.h>
#include <string.h>
#include <float.h>

void sia2_trajectory_init(sia2_trajectory *t)
{
    if (t) memset(t, 0, sizeof(*t));
}

int sia2_trajectory_record(sia2_trajectory *t,
                           double timestamp, double capability,
                           double entropy, double residual, int step)
{
    if (!t) return SIA2_ERR_NULL;
    if (t->count >= SIA2_TRAJ_MAX_POINTS) return SIA2_ERR_SIZE;

    t->points[t->count].timestamp  = timestamp;
    t->points[t->count].capability = capability;
    t->points[t->count].entropy    = entropy;
    t->points[t->count].residual   = residual;
    t->points[t->count].step       = step;
    t->count++;
    return SIA2_OK;
}

double sia2_trajectory_convergence_rate(const sia2_trajectory *t, int last_n)
{
    int start, i, count;
    double sum_x, sum_y, sum_xy, sum_xx;
    if (!t || t->count < 2) return 0.0;
    if (last_n < 2) last_n = 2;
    if (last_n > t->count) last_n = t->count;

    start = t->count - last_n;
    count = last_n;
    sum_x = sum_y = sum_xy = sum_xx = 0.0;

    for (i = 0; i < count; i++) {
        double x = (double)(i);
        double y = log(fabs(t->points[start + i].residual) + SIA2_EPSILON);
        sum_x  += x;
        sum_y  += y;
        sum_xy += x * y;
        sum_xx += x * x;
    }

    double denom = count * sum_xx - sum_x * sum_x;
    if (fabs(denom) < SIA2_EPSILON) return 0.0;
    return (count * sum_xy - sum_x * sum_y) / denom;
}

int sia2_trajectory_plateau(const sia2_trajectory *t, int last_n, double threshold)
{
    int start, i;
    if (!t || t->count < last_n || last_n < 2) return 0;

    start = t->count - last_n;
    for (i = 1; i < last_n; i++) {
        double delta = t->points[start + i - 1].residual - t->points[start + i].residual;
        if (delta > threshold) return 0; /* still improving */
    }
    return 1; /* plateau detected */
}

int sia2_trajectory_best(const sia2_trajectory *t)
{
    int best_idx = 0, i;
    double best_res;
    if (!t || t->count == 0) return -1;

    best_res = t->points[0].residual;
    for (i = 1; i < t->count; i++) {
        if (t->points[i].residual < best_res) {
            best_res = t->points[i].residual;
            best_idx = i;
        }
    }
    return best_idx;
}

double sia2_trajectory_total_improvement(const sia2_trajectory *t)
{
    if (!t || t->count < 2) return 0.0;
    return t->points[t->count - 1].capability - t->points[0].capability;
}
