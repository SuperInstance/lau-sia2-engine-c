#ifndef SIA2_PDE_H
#define SIA2_PDE_H

#ifdef SIA2_STATIC
#define SIA2_API_PDE static
#else
#define SIA2_API_PDE extern
#endif

#include "common.h"

/*
 * PDE solvers for SIA² spectral smoothing.
 * Heat equation for diffusion-based improvement regularization.
 */

typedef struct {
    double energy;          /* ∫ u² dx  (discrete L2 energy) */
    double max_value;       /* max(u) */
    double min_value;       /* min(u) */
    double total_variation; /* Σ|u[i+1] - u[i]| */
} sia2_pde_estimate;

/*
 * Solve 1D heat equation:  u_t = α u_xx  using explicit Euler.
 * u: initial state (n points), modified in-place.
 * Returns number of steps actually taken.
 */
SIA2_API_PDE int sia2_heat_solve(double *u, int n,
                                 double alpha, double dx, double dt,
                                 int steps);

/*
 * Compute energy estimates for current state.
 */
SIA2_API_PDE sia2_pde_estimate sia2_pde_compute_estimates(const double *u, int n, double dx);

/*
 * Check maximum principle: min(u₀) ≤ u(x,t) ≤ max(u₀).
 * u0: initial state, u: current state, n: grid points.
 */
SIA2_API_PDE int sia2_check_maximum_principle(const double *u0,
                                              const double *u, int n);

/*
 * Compute CFL-stable time step for heat equation.
 * dt < dx² / (2α).
 */
SIA2_API_PDE double sia2_cfl_timestep(double alpha, double dx);

#ifdef SIA2_STATIC
#include "pde.c"
#endif
#endif /* SIA2_PDE_H */
