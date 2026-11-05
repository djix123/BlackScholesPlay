/* heston_calibrate.c
 *
 * Calibrate Heston stochastic-vol parameters (v0, kappa, theta, xi, rho)
 * to an implied-vol surface quoted by relative strike and tenor, and
 * price vanilla options semi-analytically under Heston (characteristic
 * function, "little trap" formulation, Gauss-Legendre quadrature).
 *
 * Modes:
 *   price S K T r q v0 kappa theta xi rho
 *       Print Heston call/put and implied vol for one vanilla.
 *   mc S K T r q v0 kappa theta xi rho PATHS
 *       Price the same vanilla by Heston Monte-Carlo (full-truncation Euler,
 *       antithetic) and compare with the Fourier price in standard errors.
 *   gen OUTFILE S r q v0 kappa theta xi rho [ns nt kmin kmax tmin tmax]
 *       Write an implied-vol surface file for the given Heston params
 *       (default 7x5, relative strikes 0.8..1.2, tenors 0.25..2.0).
 *   fit SURFFILE S r q [v0 kappa theta xi rho] [maxiter]
 *       Calibrate Heston params to the surface (Nelder-Mead on mean
 *       squared implied-vol error). Prints params, per-point errors and
 *       the ready-to-use `black_scholes var ...` argument order.
 *
 * Surface file format ('#' starts a comment):
 *   ns nt
 *   krel_1 ... krel_ns        relative strikes k = K/F
 *   t_1 ... t_nt              tenors in years
 *   vol[strike1][tenor1] ...  ns rows of nt implied vols
 *
 * Build: gcc -O2 -o heston_calibrate heston_calibrate.c -lm
 */

#include <complex.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Fourier integration: composite Gauss-Legendre over panels of width
 * PANEL_WIDTH starting at u = 0, continuing until the integrand has died out
 * (two consecutive panels each contributing < TAIL_TOL in probability units)
 * or UMAX_LIMIT is reached. The cutoff therefore adapts to the tenor: short
 * tenors have slowly decaying characteristic functions and need a far larger
 * cutoff than a fixed u_max would give. */
#define PANEL_NODES 32
#define PANEL_WIDTH 10.0
#define UMAX_LIMIT 40000.0
#define TAIL_TOL 1e-14

#define MAX_DIM 10000          /* max strikes or tenors in a surface file */
#define MAX_POINTS 1000000     /* max ns*nt */
#define VOL_PENALTY 1.0        /* objective error for a non-invertible point */

static double g_gx[PANEL_NODES], g_gw[PANEL_NODES];
static int g_quad_ready = 0;

static double clampd(double x, double lo, double hi)
{
    return x < lo ? lo : (x > hi ? hi : x);
}

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) {
        fprintf(stderr, "Out of memory\n");
        exit(EXIT_FAILURE);
    }
    return p;
}

/* ---------------- strict argument parsing ---------------- */

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

/* Parse n numeric args argv[0..n-1]; names[] is used in the error message. */
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

/* Seed drand48 once per process; BS_SEED=<n> gives a reproducible run. */
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
    s *= 0x9E3779B97F4A7C15ULL;
    s ^= s >> 29;
    unsigned short seed[3] = {(unsigned short)(s & 0xffff),
                              (unsigned short)((s >> 16) & 0xffff),
                              (unsigned short)((s >> 32) & 0xffff)};
    seed48(seed);
}

/* ---------------- Black-Scholes and implied vol ---------------- */

static double norm_cdf(double x)
{
    return 0.5 * erfc(-x / sqrt(2.0));
}

/* Black-Scholes call (is_call=1) or put. */
static double bs_price(int is_call, double S, double K, double T, double r, double q,
                       double sigma)
{
    double dq = exp(-q * T), dr = exp(-r * T);
    if (T <= 0.0 || sigma <= 0.0)
        return is_call ? fmax(S * dq - K * dr, 0.0) : fmax(K * dr - S * dq, 0.0);
    double sq = sigma * sqrt(T);
    double d1 = (log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / sq;
    double d2 = d1 - sq;
    return is_call ? S * dq * norm_cdf(d1) - K * dr * norm_cdf(d2)
                   : K * dr * norm_cdf(-d2) - S * dq * norm_cdf(-d1);
}

/* Implied vol from a (model) CALL price, inverted on the out-of-the-money
 * side: the call for K >= F, otherwise the put obtained by put-call parity.
 * OTM options carry far more vega per unit of price error than deep-ITM calls,
 * so the inversion stays well conditioned in the wings.
 * Returns NaN if the price is outside the no-arbitrage range (<= 0 or at/above
 * the upper bound) or no vol in [1e-6, 10] reproduces it; callers must handle
 * that rather than receive a silently clamped vol. */
static double bs_otm_implied_vol(double call_price, double S, double K, double T,
                                 double r, double q)
{
    double dq = exp(-q * T), dr = exp(-r * T);
    double F = S * exp((r - q) * T);
    int use_call = K >= F;
    double price = use_call ? call_price : call_price - S * dq + K * dr;
    double upper = use_call ? S * dq : K * dr;
    if (!isfinite(price) || !(price > 0.0) || price >= upper)
        return NAN;

    double lo = 1e-6, hi = 10.0;
    if (bs_price(use_call, S, K, T, r, q, lo) > price ||
        bs_price(use_call, S, K, T, r, q, hi) < price)
        return NAN;
    for (int i = 0; i < 100; i++) {
        double mid = 0.5 * (lo + hi);
        if (bs_price(use_call, S, K, T, r, q, mid) < price)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* ---------------- Heston semi-analytic pricer ---------------- */

/* Gauss-Legendre nodes/weights on [-1,1] (Numerical Recipes gauleg). */
static void gauleg(int n, double *x, double *w)
{
    int m = (n + 1) / 2;
    for (int i = 0; i < m; i++) {
        double z = cos(M_PI * (i + 0.75) / (n + 0.5));
        double z1 = 0.0, pp = 0.0;
        do {
            double p1 = 1.0, p2 = 0.0;
            for (int j = 0; j < n; j++) {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > 1e-15);
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
}

static void ensure_quad(void)
{
    if (!g_quad_ready) {
        gauleg(PANEL_NODES, g_gx, g_gw);
        g_quad_ready = 1;
    }
}

/* Heston risk-neutral characteristic function of x_T = ln S_T,
 * "little trap" formulation (Albrecher et al. 2007):
 *   dS/S = (r-q) dt + sqrt(v) dW1
 *   dv   = kappa(theta - v) dt + xi sqrt(v) dW2,   dW1 dW2 = rho dt
 * Accepts complex u; u - i yields the stock-measure CF.
 */
static double complex heston_cf(double complex u, double S, double r, double q,
                                double T, double v0, double kappa, double theta,
                                double xi, double rho)
{
    double complex d = csqrt((rho * xi * I * u - kappa) * (rho * xi * I * u - kappa) +
                             xi * xi * (u * u + I * u));
    if (creal(d) < 0.0)
        d = -d;
    double complex a = kappa - rho * xi * I * u;
    double complex g = (a - d) / (a + d);
    double complex ed = cexp(-d * T);
    double complex C = (kappa * theta / (xi * xi)) *
                       ((a - d) * T - 2.0 * clog((1.0 - g * ed) / (1.0 - g)));
    double complex D = ((a - d) / (xi * xi)) * (1.0 - ed) / (1.0 - g * ed);
    return cexp(C + D * v0 + I * u * (log(S) + (r - q) * T));
}

/* Heston European calls for nK strikes at one maturity:
 *   C = S e^{-qT} P1 - K e^{-rT} P2
 *   P2 = 1/2 + (1/pi) int_0^inf Re[e^{-iu lnK} phi(u) / (iu)] du
 *   P1 = 1/2 + (1/pi) int_0^inf Re[e^{-iu lnK} phi(u-i) / (iu phi(-i))] du
 * with phi(-i) = F = S e^{(r-q)T}. The characteristic function does not depend
 * on K, so it is evaluated once per node and shared across all strikes. */
static void heston_calls(double S, const double *K, int nK, double T, double r, double q,
                         double v0, double kappa, double theta, double xi, double rho,
                         double *calls)
{
    ensure_quad();
    double F = S * exp((r - q) * T);
    double *lnK = xmalloc((size_t)nK * sizeof *lnK);
    double *p1 = xmalloc((size_t)nK * sizeof *p1);
    double *p2 = xmalloc((size_t)nK * sizeof *p2);
    double *d1 = xmalloc((size_t)nK * sizeof *d1);
    double *d2 = xmalloc((size_t)nK * sizeof *d2);
    for (int i = 0; i < nK; i++) {
        lnK[i] = log(K[i]);
        p1[i] = p2[i] = 0.5;
    }

    int quiet = 0;
    for (double a = 0.0; a < UMAX_LIMIT && quiet < 2; a += PANEL_WIDTH) {
        for (int i = 0; i < nK; i++)
            d1[i] = d2[i] = 0.0;
        for (int n = 0; n < PANEL_NODES; n++) {
            double u = a + 0.5 * PANEL_WIDTH * (g_gx[n] + 1.0);
            double w = 0.5 * PANEL_WIDTH * g_gw[n] / M_PI;
            double complex ph_u = heston_cf(u, S, r, q, T, v0, kappa, theta, xi, rho);
            double complex ph_um = heston_cf(u - I, S, r, q, T, v0, kappa, theta, xi, rho);
            double complex f2 = ph_u / (I * u);
            double complex f1 = ph_um / (I * u * F);
            for (int i = 0; i < nK; i++) {
                double complex e = cexp(-I * u * lnK[i]);
                d2[i] += w * creal(e * f2);
                d1[i] += w * creal(e * f1);
            }
        }
        double biggest = 0.0;
        for (int i = 0; i < nK; i++) {
            p1[i] += d1[i];
            p2[i] += d2[i];
            biggest = fmax(biggest, fmax(fabs(d1[i]), fabs(d2[i])));
        }
        quiet = biggest < TAIL_TOL ? quiet + 1 : 0;
    }

    for (int i = 0; i < nK; i++) {
        double c = S * exp(-q * T) * p1[i] - K[i] * exp(-r * T) * p2[i];
        calls[i] = c > 0.0 ? c : 0.0;
    }
    free(lnK);
    free(p1);
    free(p2);
    free(d1);
    free(d2);
}

static double heston_call(double S, double K, double T, double r, double q,
                          double v0, double kappa, double theta, double xi, double rho)
{
    double c;
    heston_calls(S, &K, 1, T, r, q, v0, kappa, theta, xi, rho, &c);
    return c;
}

/* ---------------- implied-vol surface ---------------- */

typedef struct {
    int ns, nt;
    double *krel;  /* relative strikes k = K/F, ns entries */
    double *tenor; /* tenors in years, nt entries */
    double *vol;   /* vol[i*nt + j] = impliedvols[strike i][tenor j] */
} surface;

static void surface_init(surface *s)
{
    s->ns = s->nt = 0;
    s->krel = s->tenor = s->vol = NULL;
}

static void surface_free(surface *s)
{
    free(s->krel);
    free(s->tenor);
    free(s->vol);
    surface_init(s);
}

/* Read the next whitespace-separated token, skipping '#' comments (which may
 * also trail data on a line). Returns 0 on a token, -1 at EOF, -2 if too long. */
static int next_token(FILE *f, char *buf, size_t n)
{
    int c;
    for (;;) {
        c = fgetc(f);
        if (c == EOF)
            return -1;
        if (c == '#') {
            do {
                c = fgetc(f);
            } while (c != EOF && c != '\n');
            continue;
        }
        if (!isspace(c))
            break;
    }
    size_t len = 0;
    while (c != EOF && !isspace(c) && c != '#') {
        if (len + 1 >= n)
            return -2;
        buf[len++] = (char)c;
        c = fgetc(f);
    }
    if (c == '#')
        ungetc(c, f);
    buf[len] = '\0';
    return 0;
}

static int next_dbl(FILE *f, double *out)
{
    char tok[64];
    return next_token(f, tok, sizeof tok) == 0 ? parse_double(tok, out) : -1;
}

static int next_dim(FILE *f, int *out)
{
    char tok[64];
    long v;
    if (next_token(f, tok, sizeof tok) != 0 || parse_long(tok, &v) != 0)
        return -1;
    if (v < 1 || v > MAX_DIM)
        return -1;
    *out = (int)v;
    return 0;
}

/* Load and validate a surface; on failure prints why, frees everything and
 * returns -1. */
static int surface_load(surface *s, const char *path)
{
    surface_init(s);
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Cannot open %s: %s\n", path, strerror(errno));
        return -1;
    }
    const char *why = NULL;
    if (next_dim(f, &s->ns) != 0 || next_dim(f, &s->nt) != 0) {
        why = "bad header: expected 'ns nt' with 1 <= ns, nt <= 10000";
        goto fail;
    }
    if ((long)s->ns * (long)s->nt > MAX_POINTS) {
        why = "surface too large (ns*nt > 1000000)";
        goto fail;
    }
    s->krel = xmalloc((size_t)s->ns * sizeof *s->krel);
    s->tenor = xmalloc((size_t)s->nt * sizeof *s->tenor);
    s->vol = xmalloc((size_t)s->ns * (size_t)s->nt * sizeof *s->vol);
    for (int i = 0; i < s->ns; i++)
        if (next_dbl(f, &s->krel[i]) != 0 || s->krel[i] <= 0.0) {
            why = "bad or missing relative strike (need a number > 0)";
            goto fail;
        }
    for (int j = 0; j < s->nt; j++)
        if (next_dbl(f, &s->tenor[j]) != 0 || s->tenor[j] <= 0.0) {
            why = "bad or missing tenor (need a number > 0)";
            goto fail;
        }
    for (int i = 0; i < s->ns * s->nt; i++)
        if (next_dbl(f, &s->vol[i]) != 0 || s->vol[i] <= 0.0) {
            why = "bad or missing implied vol (need ns*nt numbers > 0)";
            goto fail;
        }
    char extra[64];
    if (next_token(f, extra, sizeof extra) != -1) {
        why = "unexpected data after the last implied vol";
        goto fail;
    }
    fclose(f);
    return 0;

fail:
    fprintf(stderr, "%s: %s\n", path, why);
    fclose(f);
    surface_free(s);
    return -1;
}

/* Model implied vols for every (strike, tenor) of sf under the given Heston
 * parameters, into iv[i*nt + j]. Entries that cannot be inverted are NaN.
 * Returns the number of NaN entries. */
static int model_vols(const surface *sf, double S, double r, double q,
                      double v0, double kappa, double theta, double xi, double rho,
                      double *iv)
{
    double *Ks = xmalloc((size_t)sf->ns * sizeof *Ks);
    double *calls = xmalloc((size_t)sf->ns * sizeof *calls);
    int bad = 0;
    for (int j = 0; j < sf->nt; j++) {
        double T = sf->tenor[j];
        double F = S * exp((r - q) * T);
        for (int i = 0; i < sf->ns; i++)
            Ks[i] = sf->krel[i] * F;
        heston_calls(S, Ks, sf->ns, T, r, q, v0, kappa, theta, xi, rho, calls);
        for (int i = 0; i < sf->ns; i++) {
            double v = bs_otm_implied_vol(calls[i], S, Ks[i], T, r, q);
            iv[i * sf->nt + j] = v;
            if (!isfinite(v))
                bad++;
        }
    }
    free(Ks);
    free(calls);
    return bad;
}

/* ---------------- calibration ---------------- */

typedef struct {
    const surface *surf;
    double S, r, q;
    double *iv; /* scratch, ns*nt */
} cal_ctx;

/* Unconstrained x -> Heston params: positive params via exp, rho via tanh.
 * The exponents are clamped to the ranges below; a fit that ends on a bound is
 * reported by mode_fit(), which compares the final x against these bounds. */
static const double X_LO[4] = {-14.0, -9.0, -14.0, -14.0};
static const double X_HI[4] = {4.0, 4.0, 4.0, 3.0};

static void unpack(const double *x, double *v0, double *kappa, double *theta,
                   double *xi, double *rho)
{
    *v0 = exp(clampd(x[0], X_LO[0], X_HI[0]));
    *kappa = exp(clampd(x[1], X_LO[1], X_HI[1]));
    *theta = exp(clampd(x[2], X_LO[2], X_HI[2]));
    *xi = exp(clampd(x[3], X_LO[3], X_HI[3]));
    *rho = tanh(x[4]);
}

static double objective(const double *x, void *void_ctx)
{
    cal_ctx *c = void_ctx;
    double v0, kappa, theta, xi, rho;
    unpack(x, &v0, &kappa, &theta, &xi, &rho);

    int npts = c->surf->ns * c->surf->nt;
    model_vols(c->surf, c->S, c->r, c->q, v0, kappa, theta, xi, rho, c->iv);
    double sum = 0.0;
    for (int k = 0; k < npts; k++) {
        /* A point that cannot be inverted (model price outside the
         * no-arbitrage range) gets a large fixed error instead of a
         * clamped, misleading vol. */
        double e = isfinite(c->iv[k]) ? c->iv[k] - c->surf->vol[k] : VOL_PENALTY;
        sum += e * e;
    }
    return sum / npts; /* mean squared implied-vol error */
}

/* Nelder-Mead simplex minimization over n variables. x is overwritten with
 * the best point; returns the objective there. *converged is set to 1 if the
 * spread/diameter test below was met, and to 0 if the loop ran out the
 * maxiter budget without meeting it. */
static double nm_minimize(double (*f)(const double *, void *), void *ctx,
                          double *x, int n, int maxiter, int *converged)
{
    int np = n + 1;
    int converged_now = 0; /* set when the spread/diam test below is met */
    double **v = xmalloc((size_t)np * sizeof *v);
    double *fv = xmalloc((size_t)np * sizeof *fv);
    double *cen = xmalloc((size_t)n * sizeof *cen);
    double *tr = xmalloc((size_t)n * sizeof *tr);
    double *ex = xmalloc((size_t)n * sizeof *ex);
    double *xc = xmalloc((size_t)n * sizeof *xc);
    for (int i = 0; i < np; i++)
        v[i] = xmalloc((size_t)n * sizeof **v);

    for (int j = 0; j < n; j++)
        v[0][j] = x[j];
    for (int i = 1; i <= n; i++) {
        for (int j = 0; j < n; j++)
            v[i][j] = x[j];
        v[i][i - 1] += 0.25; /* vertex i perturbs coordinate i-1 */
    }
    for (int i = 0; i < np; i++)
        fv[i] = f(v[i], ctx);

    for (int iter = 0; iter < maxiter; iter++) {
        /* insertion sort, best first */
        for (int i = 1; i < np; i++) {
            double fi = fv[i];
            double *vi = v[i];
            int j = i - 1;
            while (j >= 0 && fv[j] > fi) {
                fv[j + 1] = fv[j];
                v[j + 1] = v[j];
                j--;
            }
            fv[j + 1] = fi;
            v[j + 1] = vi;
        }

        double spread = fv[np - 1] - fv[0];
        double diam = 0.0;
        for (int i = 1; i < np; i++)
            for (int j = 0; j < n; j++) {
                double d = fabs(v[i][j] - v[0][j]);
                if (d > diam)
                    diam = d;
            }
        if (spread < 1e-15 && diam < 1e-10) {
            converged_now = 1;
            break;
        }

        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int i = 0; i < n; i++)
                s += v[i][j];
            cen[j] = s / n;
        }

        int accepted = 0;
        for (int j = 0; j < n; j++)
            tr[j] = cen[j] + (cen[j] - v[np - 1][j]);
        double fr = f(tr, ctx);

        if (fr < fv[0]) {
            for (int j = 0; j < n; j++)
                ex[j] = cen[j] + 2.0 * (cen[j] - v[np - 1][j]);
            double fe = f(ex, ctx);
            if (fe < fr) {
                memcpy(v[np - 1], ex, (size_t)n * sizeof **v);
                fv[np - 1] = fe;
            } else {
                memcpy(v[np - 1], tr, (size_t)n * sizeof **v);
                fv[np - 1] = fr;
            }
            accepted = 1;
        } else if (fr < fv[np - 2]) {
            memcpy(v[np - 1], tr, (size_t)n * sizeof **v);
            fv[np - 1] = fr;
            accepted = 1;
        } else {
            double fc;
            if (fr < fv[np - 1]) {
                for (int j = 0; j < n; j++)
                    xc[j] = cen[j] + 0.5 * (tr[j] - cen[j]);
                fc = f(xc, ctx);
                if (fc <= fr)
                    accepted = 1;
            } else {
                for (int j = 0; j < n; j++)
                    xc[j] = cen[j] + 0.5 * (v[np - 1][j] - cen[j]);
                fc = f(xc, ctx);
                if (fc < fv[np - 1])
                    accepted = 1;
            }
            if (accepted) {
                memcpy(v[np - 1], xc, (size_t)n * sizeof **v);
                fv[np - 1] = fc;
            } else {
                /* shrink towards best */
                for (int i = 1; i < np; i++) {
                    for (int j = 0; j < n; j++)
                        v[i][j] = v[0][j] + 0.5 * (v[i][j] - v[0][j]);
                    fv[i] = f(v[i], ctx);
                }
            }
        }
    }

    *converged = converged_now;

    /* v[0] is best after the last sort (or shrink keeps v[0] best) */
    int best = 0;
    for (int i = 1; i < np; i++)
        if (fv[i] < fv[best])
            best = i;
    memcpy(x, v[best], (size_t)n * sizeof *x);
    double fbest = fv[best];

    for (int i = 0; i < np; i++)
        free(v[i]);
    free(v);
    free(fv);
    free(cen);
    free(tr);
    free(ex);
    free(xc);
    return fbest;
}

/* ---------------- modes ---------------- */

static int usage(const char *prog)
{
    fprintf(stderr,
            "Usage:\n"
            "  %s price S K T r q v0 kappa theta xi rho\n"
            "  %s mc S K T r q v0 kappa theta xi rho PATHS\n"
            "  %s gen OUTFILE S r q v0 kappa theta xi rho [ns nt kmin kmax tmin tmax]\n"
            "  %s fit SURFFILE S r q [v0 kappa theta xi rho] [maxiter]\n"
            "      (the 5 start parameters are all-or-nothing; maxiter may follow either way)\n"
            "\n"
            "Surface file format ('#' comments allowed):\n"
            "  ns nt\n"
            "  krel_1 ... krel_ns     relative strikes k = K/F\n"
            "  t_1 ... t_nt           tenors in years\n"
            "  vol[strike1][t1] ...   ns rows of nt implied vols\n"
            "Set BS_SEED=<n> for a reproducible `mc` run.\n",
            prog, prog, prog, prog);
    return EXIT_FAILURE;
}

static int params_ok(double S, double v0, double kappa, double theta, double xi, double rho)
{
    return S > 0.0 && v0 > 0.0 && kappa > 0.0 && theta > 0.0 && xi > 0.0 &&
           fabs(rho) < 1.0;
}

static void print_feller(double kappa, double theta, double xi)
{
    double ratio = 2.0 * kappa * theta / (xi * xi);
    printf("Feller ratio 2*kappa*theta/xi^2 = %.4f (%s)\n", ratio,
           ratio >= 1.0 ? "satisfied: variance stays positive"
                        : "violated: variance can reach zero");
}

static int mode_price(int argc, char **argv)
{
    if (argc != 12)
        return usage(argv[0]);
    static const char *const names[10] = {"S", "K", "T", "r", "q",
                                          "v0", "kappa", "theta", "xi", "rho"};
    double a[10];
    if (parse_args(argv + 2, names, 10, a) != 0)
        return EXIT_FAILURE;
    double S = a[0], K = a[1], T = a[2], r = a[3], q = a[4];
    double v0 = a[5], kappa = a[6], theta = a[7], xi = a[8], rho = a[9];
    if (K <= 0.0 || T <= 0.0 || !params_ok(S, v0, kappa, theta, xi, rho)) {
        fprintf(stderr, "Invalid input: require S>0, K>0, T>0, v0>0, kappa>0, theta>0, xi>0, |rho|<1\n");
        return EXIT_FAILURE;
    }
    double call = heston_call(S, K, T, r, q, v0, kappa, theta, xi, rho);
    double put = call - S * exp(-q * T) + K * exp(-r * T);
    double iv = bs_otm_implied_vol(call, S, K, T, r, q);
    printf("Heston vanilla  S=%.6g K=%.6g T=%.6g r=%.6g q=%.6g\n", S, K, T, r, q);
    printf("params: v0=%.6g kappa=%.6g theta=%.6g xi=%.6g rho=%.6g\n",
           v0, kappa, theta, xi, rho);
    if (isfinite(iv))
        printf("Call: %.6f  Put: %.6f  Implied vol: %.6f\n", call, put, iv);
    else
        printf("Call: %.6f  Put: %.6f  Implied vol: n/a (price outside the invertible range)\n",
               call, put);
    return EXIT_SUCCESS;
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

/* Monte-Carlo Heston vanilla (full-truncation Euler, antithetic) compared with
 * the Fourier price. This is an independent check of the semi-analytic pricer. */
static int mode_mc(int argc, char **argv)
{
    if (argc != 13)
        return usage(argv[0]);
    static const char *const names[10] = {"S", "K", "T", "r", "q",
                                          "v0", "kappa", "theta", "xi", "rho"};
    double a[10];
    long paths;
    if (parse_args(argv + 2, names, 10, a) != 0)
        return EXIT_FAILURE;
    if (parse_long(argv[12], &paths) != 0 || paths < 2) {
        fprintf(stderr, "Invalid PATHS '%s' (need an integer >= 2)\n", argv[12]);
        return EXIT_FAILURE;
    }
    double S = a[0], K = a[1], T = a[2], r = a[3], q = a[4];
    double v0 = a[5], kappa = a[6], theta = a[7], xi = a[8], rho = a[9];
    if (K <= 0.0 || T <= 0.0 || !params_ok(S, v0, kappa, theta, xi, rho)) {
        fprintf(stderr, "Invalid input: require S>0, K>0, T>0, v0>0, kappa>0, theta>0, xi>0, |rho|<1\n");
        return EXIT_FAILURE;
    }

    rng_seed_once();
    int N = (int)ceil(500.0 * T);
    if (N < 100)
        N = 100;
    double dt = T / N, sdt = sqrt(dt);
    double disc = exp(-r * T);
    double sq1mr2 = sqrt(1.0 - rho * rho);
    double csum = 0.0, csq = 0.0, psum = 0.0, psq = 0.0;

    for (long i = 0; i < paths; i++) {
        double lp = log(S), lm = log(S), vp = v0, vm = v0;
        for (int j = 0; j < N; j++) {
            double z1 = bm_rand(), z2 = bm_rand();
            double dws = sq1mr2 * z1 + rho * z2;
            double up = fmax(vp, 0.0), um = fmax(vm, 0.0);
            lp += (r - q - 0.5 * up) * dt + sqrt(up) * sdt * dws;
            vp += kappa * (theta - up) * dt + xi * sqrt(up) * sdt * z2;
            lm += (r - q - 0.5 * um) * dt - sqrt(um) * sdt * dws;
            vm += kappa * (theta - um) * dt - xi * sqrt(um) * sdt * z2;
        }
        double sp = exp(lp), sm = exp(lm);
        double c = 0.5 * disc * (fmax(sp - K, 0.0) + fmax(sm - K, 0.0));
        double p = 0.5 * disc * (fmax(K - sp, 0.0) + fmax(K - sm, 0.0));
        csum += c;
        csq += c * c;
        psum += p;
        psq += p * p;
    }
    double cm = csum / paths, pm = psum / paths;
    double cse = sqrt(fmax(csq / paths - cm * cm, 0.0) / paths);
    double pse = sqrt(fmax(psq / paths - pm * pm, 0.0) / paths);

    double fc = heston_call(S, K, T, r, q, v0, kappa, theta, xi, rho);
    double fp = fc - S * exp(-q * T) + K * exp(-r * T);
    printf("Heston Monte-Carlo (%ld paths, %d steps, full-truncation Euler, antithetic)\n",
           paths, N);
    printf("Call: %.6f +/- %.6f (Fourier %.6f, diff %+.2f SE)\n", cm, cse, fc,
           cse > 0 ? (cm - fc) / cse : 0.0);
    printf("Put:  %.6f +/- %.6f (Fourier %.6f, diff %+.2f SE)\n", pm, pse, fp,
           pse > 0 ? (pm - fp) / pse : 0.0);
    return EXIT_SUCCESS;
}

static double grid_pt(double lo, double hi, int i, int n)
{
    return n > 1 ? lo + (hi - lo) * i / (n - 1) : lo;
}

static int mode_gen(int argc, char **argv)
{
    if (argc < 11 || argc > 17)
        return usage(argv[0]);
    const char *out = argv[2];
    static const char *const names[8] = {"S", "r", "q", "v0", "kappa", "theta", "xi", "rho"};
    double a[8];
    if (parse_args(argv + 3, names, 8, a) != 0)
        return EXIT_FAILURE;
    double S = a[0], r = a[1], q = a[2];
    double v0 = a[3], kappa = a[4], theta = a[5], xi = a[6], rho = a[7];

    long ns = 7, nt = 5;
    double kmin = 0.8, kmax = 1.2, tmin = 0.25, tmax = 2.0;
    if ((argc > 11 && parse_long(argv[11], &ns) != 0) ||
        (argc > 12 && parse_long(argv[12], &nt) != 0) ||
        (argc > 13 && parse_double(argv[13], &kmin) != 0) ||
        (argc > 14 && parse_double(argv[14], &kmax) != 0) ||
        (argc > 15 && parse_double(argv[15], &tmin) != 0) ||
        (argc > 16 && parse_double(argv[16], &tmax) != 0)) {
        fprintf(stderr, "Invalid number in the optional gen arguments [ns nt kmin kmax tmin tmax]\n");
        return EXIT_FAILURE;
    }
    if (!params_ok(S, v0, kappa, theta, xi, rho) || ns < 1 || ns > MAX_DIM ||
        nt < 1 || nt > MAX_DIM || ns * nt > MAX_POINTS ||
        kmin <= 0.0 || kmax < kmin || tmin <= 0.0 || tmax < tmin) {
        fprintf(stderr, "Invalid gen input\n");
        return EXIT_FAILURE;
    }
    if ((ns == 1 && argc > 14 && kmax != kmin) || (nt == 1 && argc > 16 && tmax != tmin)) {
        fprintf(stderr, "Invalid gen input: ns=1 (nt=1) uses kmin (tmin) only, so kmax (tmax) must equal it\n");
        return EXIT_FAILURE;
    }

    surface sf;
    surface_init(&sf);
    sf.ns = (int)ns;
    sf.nt = (int)nt;
    sf.krel = xmalloc((size_t)ns * sizeof *sf.krel);
    sf.tenor = xmalloc((size_t)nt * sizeof *sf.tenor);
    sf.vol = xmalloc((size_t)(ns * nt) * sizeof *sf.vol);
    for (int i = 0; i < sf.ns; i++)
        sf.krel[i] = grid_pt(kmin, kmax, i, sf.ns);
    for (int j = 0; j < sf.nt; j++)
        sf.tenor[j] = grid_pt(tmin, tmax, j, sf.nt);

    int bad = model_vols(&sf, S, r, q, v0, kappa, theta, xi, rho, sf.vol);
    if (bad) {
        fprintf(stderr, "Cannot write surface: %d point(s) have a model price outside the "
                        "invertible range (e.g. extreme strikes at short tenors):\n", bad);
        for (int i = 0; i < sf.ns; i++)
            for (int j = 0; j < sf.nt; j++)
                if (!isfinite(sf.vol[i * sf.nt + j]))
                    fprintf(stderr, "  krel=%.6g T=%.6g\n", sf.krel[i], sf.tenor[j]);
        surface_free(&sf);
        return EXIT_FAILURE;
    }

    FILE *f = fopen(out, "w");
    if (!f) {
        fprintf(stderr, "Cannot write %s: %s\n", out, strerror(errno));
        surface_free(&sf);
        return EXIT_FAILURE;
    }
    fprintf(f, "# implied vol surface (relative strike x tenor) under Heston\n");
    fprintf(f, "# S=%.10g r=%.10g q=%.10g v0=%.10g kappa=%.10g theta=%.10g xi=%.10g rho=%.10g\n",
            S, r, q, v0, kappa, theta, xi, rho);
    fprintf(f, "%d %d\n", sf.ns, sf.nt);
    for (int i = 0; i < sf.ns; i++)
        fprintf(f, "%.10g%c", sf.krel[i], i + 1 < sf.ns ? ' ' : '\n');
    for (int j = 0; j < sf.nt; j++)
        fprintf(f, "%.10g%c", sf.tenor[j], j + 1 < sf.nt ? ' ' : '\n');
    for (int i = 0; i < sf.ns; i++)
        for (int j = 0; j < sf.nt; j++)
            fprintf(f, "%.10g%c", sf.vol[i * sf.nt + j], j + 1 < sf.nt ? ' ' : '\n');
    int werr = ferror(f);
    if (fclose(f) != 0 || werr) {
        fprintf(stderr, "Error writing %s\n", out);
        surface_free(&sf);
        return EXIT_FAILURE;
    }
    printf("Wrote %d x %d implied-vol surface to %s\n", sf.ns, sf.nt, out);
    surface_free(&sf);
    return EXIT_SUCCESS;
}

static int mode_fit(int argc, char **argv)
{
    if (argc < 6)
        return usage(argv[0]);
    /* Accepted forms: SURFFILE S r q | ... maxiter | ... 5 params | ... 5 params maxiter */
    if (argc != 6 && argc != 7 && argc != 11 && argc != 12)
        return usage(argv[0]);
    const char *path = argv[2];

    static const char *const names3[3] = {"S", "r", "q"};
    double a3[3];
    if (parse_args(argv + 3, names3, 3, a3) != 0)
        return EXIT_FAILURE;
    double S = a3[0], r = a3[1], q = a3[2];

    double v0 = 0.04, kappa = 2.0, theta = 0.04, xi = 0.3, rho = -0.7;
    long maxiter = 2000;
    if (argc >= 11) {
        static const char *const names5[5] = {"v0", "kappa", "theta", "xi", "rho"};
        double a5[5];
        if (parse_args(argv + 6, names5, 5, a5) != 0)
            return EXIT_FAILURE;
        v0 = a5[0];
        kappa = a5[1];
        theta = a5[2];
        xi = a5[3];
        rho = a5[4];
    }
    const char *miter = argc == 7 ? argv[6] : (argc == 12 ? argv[11] : NULL);
    if (miter && parse_long(miter, &maxiter) != 0) {
        fprintf(stderr, "Invalid maxiter '%s'\n", miter);
        return EXIT_FAILURE;
    }
    if (!params_ok(S, v0, kappa, theta, xi, rho) || maxiter < 1 || maxiter > 100000000L) {
        fprintf(stderr, "Invalid fit input: require S>0, v0>0, kappa>0, theta>0, xi>0, |rho|<1, 1<=maxiter<=1e8\n");
        return EXIT_FAILURE;
    }

    surface surf;
    if (surface_load(&surf, path) != 0)
        return EXIT_FAILURE;

    int npts = surf.ns * surf.nt;
    double x[5] = {
        log(v0), log(kappa), log(theta), log(xi),
        atanh(clampd(rho, -0.999, 0.999))
    };
    double *iv = xmalloc((size_t)npts * sizeof *iv);
    cal_ctx ctx = { &surf, S, r, q, iv };
    int converged = 0;
    double fval = nm_minimize(objective, &ctx, x, 5, (int)maxiter, &converged);

    double v0f, kappaf, thetaf, xif, rhof;
    unpack(x, &v0f, &kappaf, &thetaf, &xif, &rhof);
    int bad = model_vols(&surf, S, r, q, v0f, kappaf, thetaf, xif, rhof, iv);

    static const char *const pname[4] = {"v0", "kappa", "theta", "xi"};
    double sse = 0.0, maxe = 0.0;
    int nvalid = 0;
    printf("Calibration of %s (%d strikes x %d tenors), S=%.6g r=%.6g q=%.6g\n",
           path, surf.ns, surf.nt, S, r, q);
    printf("v0=%.6f kappa=%.6f theta=%.6f xi=%.6f rho=%.6f\n",
           v0f, kappaf, thetaf, xif, rhof);
    print_feller(kappaf, thetaf, xif);
    if (!converged)
        printf("warning: stopped at maxiter=%ld without converging; increase maxiter\n", maxiter);
    for (int k = 0; k < 4; k++) {
        if (x[k] <= X_LO[k] || x[k] >= X_HI[k])
            printf("warning: %s ended on its search bound; the fit is not at an interior optimum\n",
                   pname[k]);
    }
    if (fabs(x[4]) > 6.0)
        printf("warning: rho is numerically +/-1\n");
    printf("mean squared vol error %.3e\n\n", fval);
    printf("%-8s %-8s %10s %10s %10s\n", "krel", "T", "mkt_vol", "fit_vol", "err");
    for (int i = 0; i < surf.ns; i++) {
        for (int j = 0; j < surf.nt; j++) {
            double mkt = surf.vol[i * surf.nt + j];
            double fit = iv[i * surf.nt + j];
            if (!isfinite(fit)) {
                printf("%-8.4f %-8.4f %10.6f %10s %10s\n", surf.krel[i], surf.tenor[j],
                       mkt, "n/a", "n/a");
                continue;
            }
            double e = fit - mkt;
            sse += e * e;
            nvalid++;
            if (fabs(e) > maxe)
                maxe = fabs(e);
            printf("%-8.4f %-8.4f %10.6f %10.6f %10.2e\n",
                   surf.krel[i], surf.tenor[j], mkt, fit, e);
        }
    }
    if (nvalid)
        printf("\nRMSE %.3e vol points, max |err| %.3e\n", sqrt(sse / nvalid), maxe);
    if (bad)
        printf("warning: %d of %d points could not be inverted (excluded from RMSE; "
               "they were penalised in the fit)\n", bad, npts);
    printf("\nFeed into the pricer (var-mode arg order v0 kappa theta xi rho):\n");
    printf("  ./black_scholes var <S> <T> <r> <q> <sigma> <Kvar> %.6f %.6f %.6f %.6f %.6f mc <PATHS>\n",
           v0f, kappaf, thetaf, xif, rhof);

    free(iv);
    surface_free(&surf);
    return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "price") == 0)
        return mode_price(argc, argv);
    if (argc >= 2 && strcmp(argv[1], "mc") == 0)
        return mode_mc(argc, argv);
    if (argc >= 2 && strcmp(argv[1], "gen") == 0)
        return mode_gen(argc, argv);
    if (argc >= 2 && strcmp(argv[1], "fit") == 0)
        return mode_fit(argc, argv);
    return usage(argv[0]);
}
