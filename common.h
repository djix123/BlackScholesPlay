/* common.h - helpers shared by black_scholes.c and heston_calibrate.c. */
#ifndef COMMON_H
#define COMMON_H

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Most time steps a Monte-Carlo path may take; also keeps the count in an int. */
#define MC_MAX_STEPS 1e7

static double norm_cdf(double x)
{
    return 0.5 * erfc(-x / sqrt(2.0));
}

/* Strict number parsing: the whole string must be a finite number. */
static int parse_double(const char *s, double *out)
{
    char *end;
    errno = 0;
    double v = strtod(s, &end);
    if (end == s || *end != '\0' || errno == ERANGE || !isfinite(v))
        return -1;
    *out = v;
    return 0;
}

static int parse_long(const char *s, long *out)
{
    char *end;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (end == s || *end != '\0' || errno == ERANGE)
        return -1;
    *out = v;
    return 0;
}

/* Parse n numeric arguments argv[0..n-1] into out[]; names[] is for errors.
 * Returns 0 on success, -1 (after printing a message) on the first bad one. */
static int parse_args(char **argv, const char *const *names, int n, double *out)
{
    for (int i = 0; i < n; i++) {
        if (parse_double(argv[i], &out[i]) != 0) {
            fprintf(stderr, "Invalid number '%s' for %s\n", argv[i], names[i]);
            return -1;
        }
    }
    return 0;
}

/* Seed drand48 once per process. Set BS_SEED=<unsigned integer> for a
 * reproducible run; otherwise the seed mixes the time, pid and clock. */
static void rng_seed_once(void)
{
    static int done = 0;
    if (done)
        return;
    done = 1;

    unsigned long long s = 0;
    const char *env = getenv("BS_SEED");
    int have = 0;
    if (env && *env) {
        char *end;
        errno = 0;
        s = strtoull(env, &end, 10);
        if (*end == '\0' && errno != ERANGE)
            have = 1;
        else
            fprintf(stderr, "Ignoring invalid BS_SEED '%s'\n", env);
    }
    if (!have)
        s = (unsigned long long)time(NULL) ^
            ((unsigned long long)getpid() << 20) ^ (unsigned long long)clock();
    s *= 0x9E3779B97F4A7C15ULL; /* spread small seeds across all 48 bits */
    s ^= s >> 29;
    unsigned short seed[3] = {(unsigned short)(s & 0xffff),
                              (unsigned short)((s >> 16) & 0xffff),
                              (unsigned short)((s >> 32) & 0xffff)};
    seed48(seed);
}

/* Standard normal via Box-Muller on drand48(). */
static double bm_rand(void)
{
    double u1 = drand48();
    double u2 = drand48();
    while (u1 <= 0.0)
        u1 = drand48();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

#endif
