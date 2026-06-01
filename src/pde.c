#include "pde.h"
#include <math.h>
#include <float.h>

int sia2_heat_solve(double *u, int n, double alpha, double dx, double dt,
                    int steps)
{
    double *buf = (double *)malloc(n * sizeof(double));
    int s, i;
    double r = alpha * dt / (dx * dx);

    if (!u || n < 3 || !buf) {
        free(buf);
        return 0;
    }

    /* Safety: clamp r for stability (CFL) */
    if (r > 0.5) r = 0.45;

    for (s = 0; s < steps; s++) {
        buf[0] = u[0]; /* Dirichlet BC: u[0] fixed */
        for (i = 1; i < n - 1; i++) {
            buf[i] = u[i] + r * (u[i - 1] - 2.0 * u[i] + u[i + 1]);
        }
        buf[n - 1] = u[n - 1]; /* Dirichlet BC: u[n-1] fixed */
        memcpy(u, buf, n * sizeof(double));
    }
    free(buf);
    return steps;
}

sia2_pde_estimate sia2_pde_compute_estimates(const double *u, int n, double dx)
{
    sia2_pde_estimate e;
    int i;
    memset(&e, 0, sizeof(e));

    if (!u || n < 1) return e;

    e.max_value = u[0];
    e.min_value = u[0];
    for (i = 0; i < n; i++) {
        e.energy += u[i] * u[i] * dx;
        if (u[i] > e.max_value) e.max_value = u[i];
        if (u[i] < e.min_value) e.min_value = u[i];
    }
    for (i = 0; i < n - 1; i++) {
        e.total_variation += fabs(u[i + 1] - u[i]);
    }
    return e;
}

int sia2_check_maximum_principle(const double *u0, const double *u, int n)
{
    int i;
    double min0, max0;
    if (!u0 || !u || n < 1) return 0;

    min0 = u0[0]; max0 = u0[0];
    for (i = 1; i < n; i++) {
        if (u0[i] < min0) min0 = u0[i];
        if (u0[i] > max0) max0 = u0[i];
    }
    for (i = 0; i < n; i++) {
        if (u[i] < min0 - SIA2_EPSILON || u[i] > max0 + SIA2_EPSILON)
            return 0;
    }
    return 1;
}

double sia2_cfl_timestep(double alpha, double dx)
{
    if (alpha <= 0.0 || dx <= 0.0) return 0.0;
    return 0.4 * dx * dx / alpha;  /* safety factor 0.4 < 0.5 */
}
