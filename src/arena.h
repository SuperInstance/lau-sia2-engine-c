#ifndef SIA2_ARENA_H
#define SIA2_ARENA_H

/*
 * Arena allocator — no malloc on hot paths.
 * Usage: sia2_arena_init(&a, buf, sizeof(buf));
 *        p = sia2_arena_alloc(&a, n * sizeof(double));
 */

#include <stddef.h>
#include <string.h>
#include <assert.h>

typedef struct {
    unsigned char *base;
    size_t         capacity;
    size_t         offset;
} sia2_arena;

static inline void sia2_arena_init(sia2_arena *a, void *buf, size_t cap)
{
    a->base     = (unsigned char *)buf;
    a->capacity = cap;
    a->offset   = 0;
}

static inline void *sia2_arena_alloc(sia2_arena *a, size_t bytes)
{
    /* align to 8 bytes */
    size_t aligned = (a->offset + 7u) & ~(size_t)7u;
    if (aligned + bytes > a->capacity) return NULL;
    void *p = a->base + aligned;
    a->offset = aligned + bytes;
    memset(p, 0, bytes);
    return p;
}

static inline void sia2_arena_reset(sia2_arena *a) { a->offset = 0; }

static inline size_t sia2_arena_used(const sia2_arena *a) { return a->offset; }

#endif /* SIA2_ARENA_H */
