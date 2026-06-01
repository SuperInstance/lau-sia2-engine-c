#ifndef SIA2_COMMON_H
#define SIA2_COMMON_H

#define SIA2_VERSION_MAJOR 1
#define SIA2_VERSION_MINOR 0
#define SIA2_VERSION_PATCH 0

/* Return codes */
#define SIA2_OK              0
#define SIA2_ERR_NULL       -1
#define SIA2_ERR_SIZE       -2
#define SIA2_ERR_CONVERGE   -3
#define SIA2_ERR_SINGULAR   -4
#define SIA2_ERR_VIOLATION  -5

/* Limits */
#define SIA2_MAX_DIM        64
#define SIA2_MAX_ITERS      10000
#define SIA2_EPSILON        1e-12
#define SIA2_PI             3.14159265358979323846

#endif /* SIA2_COMMON_H */
