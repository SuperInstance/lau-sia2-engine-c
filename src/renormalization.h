#ifndef SIA2_RENORMALIZATION_H
#define SIA2_RENORMALIZATION_H

#ifdef SIA2_STATIC
#define SIA2_API_RG static
#else
#define SIA2_API_RG extern
#endif

#include "common.h"

/*
 * Renormalization group (RG) analysis for SIA².
 * Detects scale-invariant improvement behavior and universality classes.
 */

typedef enum {
    SIA2_RG_TRIVIAL,       /* β ≈ 0: flows to trivial fixed point */
    SIA2_RG_GAUSSIAN,      /* β small: Gaussian-like */
    SIA2_RG_NONTRIVIAL,    /* β significant: interesting fixed point */
    SIA2_RG_UNSTABLE       /* |β| > 1: runaway flow */
} sia2_rg_class;

typedef struct {
    double beta;           /* RG beta function value */
    double gamma;          /* anomalous dimension */
    double nu;             /* correlation length exponent */
    int    is_fixed_point; /* |β| < threshold */
    sia2_rg_class universality_class;
} sia2_rg_result;

/*
 * Compute RG beta function: β(g) = d·g / d·ln(Λ).
 * Uses finite difference from measurements at scales.
 * g: coupling values at different scales (len points).
 * scales: log-scale values (ln Λ).
 */
SIA2_API_RG sia2_rg_result sia2_rg_beta(const double *g, const double *scales,
                                        int len, double threshold);

/*
 * Detect fixed point: |β(g*)| < threshold.
 */
SIA2_API_RG int sia2_rg_detect_fixed_point(const double *g, int len,
                                           double threshold);

/*
 * Classify universality class from beta function.
 */
SIA2_API_RG sia2_rg_class sia2_rg_classify(double beta, double threshold);

/*
 * Compute correlation length exponent: ν = 1 / |β'(g*)|.
 * Approximates β' via finite differences.
 */
SIA2_API_RG double sia2_rg_nu_exponent(const double *beta_values, int len);

#ifdef SIA2_STATIC
#include "renormalization.c"
#endif
#endif /* SIA2_RENORMALIZATION_H */
