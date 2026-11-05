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
 *   gen OUTFILE S r q v0 kappa theta xi rho [ns nt kmin kmax tmin tmax]
 *       Write an implied-vol surface file for the given Heston params
 *       (default 7x5, relative strikes 0.8..1.2, tenors 0.25..2.0).
 *   fit SURFFILE S r q [v0 kappa theta xi rho] [maxiter]
 *       Calibrate Heston params to the surface (Nelder-Mead on mean
 *       squared implied-vol error). Prints params, per-point errors and
 *       the ready-to-use `black_scholes var ...` argument order.
 *
 * Surface file format ('#' starts a comment line):
 *   ns nt
 *   krel_1 ... krel_ns        relative strikes k = K/F
 *   t_1 ... t_nt              tenors in years
 *   vol[strike1][tenor1] ...  ns rows of nt implied vols
 *
 * Build: gcc -O2 -o heston_calibrate heston_calibrate.c -lm
 */

#include <complex.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NG 128   /* Gauss-Legendre nodes */
#define UMAX 200.0 /* integration cutoff in Fourier space */

static double g_gx[NG], g_gw[NG];
static int g_quad_ready = 0;

static double clampd(double x, double lo, double hi)
{
    return x < lo ? lo : (x > hi ? hi : x);
}

static double norm_cdf(double x)
{
    return 0.5 * erfc(-x / sqrt(2.0));
}

/* Black-Scholes call, used for implied-vol inversion. */
static double bs_call(double S, double K, double T, double r, double q, double sigma)
{
    if (T <= 0.0 || sigma <= 0.0)
        return fmax(S * exp(-q * T) - K * exp(-r * T), 0.0);
    double sq = sigma * sqrt(T);
    double d1 = (log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / sq;
    double d2 = d1 - sq;
    return S * exp(-q * T) * norm_cdf(d1) - K * exp(-r * T) * norm_cdf(d2);
}

static double bs_implied_vol(double price, double S, double K, double T, double r, double q)
{
    double lo = 1e-6, hi = 10.0;
    double intrinsic = fmax(S * exp(-q * T) - K * exp(-r * T), 0.0);
    double cap = S * exp(-q * T);
    if (price <= intrinsic + 1e-12)
        return lo;
    if (price >= cap - 1e-12)
        return hi;
    for (int i = 0; i < 200; i++) {
        double mid = 0.5 * (lo + hi);
        if (bs_call(S, K, T, r, q, mid) < price)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

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
        gauleg(NG, g_gx, g_gw);
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

/* Heston European call:
 *   C = S e^{-qT} P1 - K e^{-rT} P2
 *   P2 = 1/2 + (1/pi) int_0^UMAX Re[e^{-iu lnK} phi(u) / (iu)] du
 *   P1 = 1/2 + (1/pi) int_0^UMAX Re[e^{-iu lnK} phi(u-i) / (iu phi(-i))] du
 * with phi(-i) = F = S e^{(r-q)T}.
 */
static double heston_call(double S, double K, double T, double r, double q,
                          double v0, double kappa, double theta, double xi, double rho)
{
    ensure_quad();
    double lnK = log(K);
    double F = S * exp((r - q) * T);
    double p1 = 0.5, p2 = 0.5;
    for (int i = 0; i < NG; i++) {
        double u = 0.5 * UMAX * (g_gx[i] + 1.0);
        double w = 0.5 * UMAX * g_gw[i];
        double complex e = cexp(-I * u * lnK);
        double complex ph_u = heston_cf(u, S, r, q, T, v0, kappa, theta, xi, rho);
        double complex ph_um = heston_cf(u - I, S, r, q, T, v0, kappa, theta, xi, rho);
        p2 += w * creal(e * ph_u / (I * u)) / M_PI;
        p1 += w * creal(e * ph_um / (I * u * F)) / M_PI;
    }
    double call = S * exp(-q * T) * p1 - K * exp(-r * T) * p2;
    return call > 0.0 ? call : 0.0;
}

/* ---------------- implied-vol surface ---------------- */

typedef struct {
    int ns, nt;
    double *krel;  /* relative strikes k = K/F, ns entries */
    double *tenor; /* tenors in years, nt entries */
    double *vol;   /* vol[i*nt + j] = impliedvols[strike i][tenor j] */
} surface;

static void surface_free(surface *s)
{
    free(s->krel);
    free(s->tenor);
    free(s->vol);
    s->krel = s->tenor = s->vol = NULL;
}

/* Skip whitespace and '#' comment lines, then read one number. */
static int next_int(FILE *f, int *out)
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
        if (isspace(c))
            continue;
        break;
    }
    ungetc(c, f);
    return fscanf(f, "%d", out) == 1 ? 0 : -1;
}

static int next_dbl(FILE *f, double *out)
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
        if (isspace(c))
            continue;
        break;
    }
    ungetc(c, f);
    return fscanf(f, "%lf", out) == 1 ? 0 : -1;
}

static int surface_load(surface *s, const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f)
        return -1;
    if (next_int(f, &s->ns) != 0 || next_int(f, &s->nt) != 0 ||
        s->ns < 1 || s->nt < 1) {
        fclose(f);
        return -1;
    }
    s->krel = malloc((size_t)s->ns * sizeof *s->krel);
    s->tenor = malloc((size_t)s->nt * sizeof *s->tenor);
    s->vol = malloc((size_t)s->ns * (size_t)s->nt * sizeof *s->vol);
    if (!s->krel || !s->tenor || !s->vol) {
        fclose(f);
        return -1;
    }
    for (int i = 0; i < s->ns; i++)
        if (next_dbl(f, &s->krel[i]) != 0 || s->krel[i] <= 0.0) {
            fclose(f);
            return -1;
        }
    for (int j = 0; j < s->nt; j++)
        if (next_dbl(f, &s->tenor[j]) != 0 || s->tenor[j] <= 0.0) {
            fclose(f);
            return -1;
        }
    for (int i = 0; i < s->ns * s->nt; i++)
        if (next_dbl(f, &s->vol[i]) != 0) {
            fclose(f);
            return -1;
        }
    fclose(f);
    return 0;
}

/* ---------------- calibration ---------------- */

typedef struct {
    const surface *surf;
    double S, r, q;
} cal_ctx;

/* Unconstrained x -> Heston params: positive params via exp, rho via tanh. */
static void unpack(const double *x, double *v0, double *kappa, double *theta,
                   double *xi, double *rho)
{
    *v0 = exp(clampd(x[0], -14.0, 4.0));
    *kappa = exp(clampd(x[1], -9.0, 4.0));
    *theta = exp(clampd(x[2], -14.0, 4.0));
    *xi = exp(clampd(x[3], -14.0, 3.0));
    *rho = tanh(x[4]);
}

static double objective(const double *x, void *void_ctx)
{
    cal_ctx *c = void_ctx;
    double v0, kappa, theta, xi, rho;
    unpack(x, &v0, &kappa, &theta, &xi, &rho);

    double sum = 0.0;
    int npts = c->surf->ns * c->surf->nt;
    for (int i = 0; i < c->surf->ns; i++) {
        for (int j = 0; j < c->surf->nt; j++) {
            double T = c->surf->tenor[j];
            double F = c->S * exp((c->r - c->q) * T);
            double K = c->surf->krel[i] * F;
            double price = heston_call(c->S, K, T, c->r, c->q,
                                       v0, kappa, theta, xi, rho);
            double lo = fmax(c->S * exp(-c->q * T) - K * exp(-c->r * T), 0.0);
            double cap = c->S * exp(-c->q * T);
            price = clampd(price, lo, cap);
            double iv = bs_implied_vol(price, c->S, K, T, c->r, c->q);
            double e = iv - c->surf->vol[i * c->surf->nt + j];
            sum += e * e;
        }
    }
    return sum / npts; /* mean squared implied-vol error */
}

/* Nelder-Mead simplex minimization over n variables. x is overwritten with
 * the best point; returns the objective there. */
static double nm_minimize(double (*f)(const double *, void *), void *ctx,
                          double *x, int n, int maxiter)
{
    int np = n + 1;
    double **v = malloc((size_t)np * sizeof *v);
    double *fv = malloc((size_t)np * sizeof *fv);
    double *cen = malloc((size_t)n * sizeof *cen);
    double *tr = malloc((size_t)n * sizeof *tr);
    double *ex = malloc((size_t)n * sizeof *ex);
    double *xc = malloc((size_t)n * sizeof *xc);
    if (!v || !fv || !cen || !tr || !ex || !xc)
        exit(EXIT_FAILURE);
    for (int i = 0; i < np; i++) {
        v[i] = malloc((size_t)n * sizeof **v);
        if (!v[i])
            exit(EXIT_FAILURE);
    }

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
        if (spread < 1e-15 && diam < 1e-10)
            break;

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
            "  %s gen OUTFILE S r q v0 kappa theta xi rho [ns nt kmin kmax tmin tmax]\n"
            "  %s fit SURFFILE S r q [v0 kappa theta xi rho] [maxiter]\n"
            "      (the 5 start parameters are all-or-nothing; maxiter may follow either way)\n"
            "\n"
            "Surface file format ('#' comment lines allowed):\n"
            "  ns nt\n"
            "  krel_1 ... krel_ns     relative strikes k = K/F\n"
            "  t_1 ... t_nt           tenors in years\n"
            "  vol[strike1][t1] ...   ns rows of nt implied vols\n",
            prog, prog, prog);
    return EXIT_FAILURE;
}

static int params_ok(double S, double v0, double kappa, double theta, double xi, double rho)
{
    return S > 0.0 && v0 > 0.0 && kappa > 0.0 && theta > 0.0 && xi > 0.0 &&
           fabs(rho) < 1.0;
}

static int mode_price(int argc, char **argv)
{
    if (argc != 12)
        return usage(argv[0]);
    double S = atof(argv[2]), K = atof(argv[3]), T = atof(argv[4]);
    double r = atof(argv[5]), q = atof(argv[6]);
    double v0 = atof(argv[7]), kappa = atof(argv[8]), theta = atof(argv[9]);
    double xi = atof(argv[10]), rho = atof(argv[11]);
    if (K <= 0.0 || T <= 0.0 || !params_ok(S, v0, kappa, theta, xi, rho)) {
        fprintf(stderr, "Invalid input: require S>0, K>0, T>0, v0>0, kappa>0, theta>0, xi>0, |rho|<1\n");
        return EXIT_FAILURE;
    }
    double call = heston_call(S, K, T, r, q, v0, kappa, theta, xi, rho);
    double put = call - S * exp(-q * T) + K * exp(-r * T);
    double iv = bs_implied_vol(call, S, K, T, r, q);
    printf("Heston vanilla  S=%.6g K=%.6g T=%.6g r=%.6g q=%.6g\n", S, K, T, r, q);
    printf("params: v0=%.6g kappa=%.6g theta=%.6g xi=%.6g rho=%.6g\n",
           v0, kappa, theta, xi, rho);
    printf("Call: %.6f  Put: %.6f  Implied vol: %.6f\n", call, put, iv);
    return EXIT_SUCCESS;
}

static int mode_gen(int argc, char **argv)
{
    if (argc < 11)
        return usage(argv[0]);
    const char *out = argv[2];
    double S = atof(argv[3]), r = atof(argv[4]), q = atof(argv[5]);
    double v0 = atof(argv[6]), kappa = atof(argv[7]), theta = atof(argv[8]);
    double xi = atof(argv[9]), rho = atof(argv[10]);
    int ns = argc > 11 ? atoi(argv[11]) : 7;
    int nt = argc > 12 ? atoi(argv[12]) : 5;
    double kmin = argc > 13 ? atof(argv[13]) : 0.8;
    double kmax = argc > 14 ? atof(argv[14]) : 1.2;
    double tmin = argc > 15 ? atof(argv[15]) : 0.25;
    double tmax = argc > 16 ? atof(argv[16]) : 2.0;
    if (!params_ok(S, v0, kappa, theta, xi, rho) || ns < 1 || nt < 1 ||
        kmin <= 0.0 || kmax < kmin || tmin <= 0.0 || tmax < tmin) {
        fprintf(stderr, "Invalid gen input\n");
        return EXIT_FAILURE;
    }

    FILE *f = fopen(out, "w");
    if (!f) {
        fprintf(stderr, "Cannot write %s\n", out);
        return EXIT_FAILURE;
    }
    fprintf(f, "# implied vol surface (relative strike x tenor) under Heston\n");
    fprintf(f, "# S=%.10g r=%.10g q=%.10g v0=%.10g kappa=%.10g theta=%.10g xi=%.10g rho=%.10g\n",
            S, r, q, v0, kappa, theta, xi, rho);
    fprintf(f, "%d %d\n", ns, nt);
    for (int i = 0; i < ns; i++) {
        double k = ns > 1 ? kmin + (kmax - kmin) * i / (ns - 1) : kmin;
        fprintf(f, "%.10g%c", k, i + 1 < ns ? ' ' : '\n');
    }
    for (int j = 0; j < nt; j++) {
        double t = nt > 1 ? tmin + (tmax - tmin) * j / (nt - 1) : tmin;
        fprintf(f, "%.10g%c", t, j + 1 < nt ? ' ' : '\n');
    }
    for (int i = 0; i < ns; i++) {
        double k = ns > 1 ? kmin + (kmax - kmin) * i / (ns - 1) : kmin;
        for (int j = 0; j < nt; j++) {
            double t = nt > 1 ? tmin + (tmax - tmin) * j / (nt - 1) : tmin;
            double F = S * exp((r - q) * t);
            double K = k * F;
            double price = heston_call(S, K, t, r, q, v0, kappa, theta, xi, rho);
            double iv = bs_implied_vol(price, S, K, t, r, q);
            fprintf(f, "%.10g%c", iv, j + 1 < nt ? ' ' : '\n');
        }
    }
    fclose(f);
    printf("Wrote %d x %d implied-vol surface to %s\n", ns, nt, out);
    return EXIT_SUCCESS;
}

static int mode_fit(int argc, char **argv)
{
    if (argc < 6)
        return usage(argv[0]);
    const char *path = argv[2];
    double S = atof(argv[3]), r = atof(argv[4]), q = atof(argv[5]);
    double v0 = 0.04, kappa = 2.0, theta = 0.04, xi = 0.3, rho = -0.7;
    int maxiter = 2000;
    /* Accepted forms: SURFFILE S r q | ... maxiter | ... 5 params | ... 5 params maxiter */
    if (argc != 6 && argc != 7 && argc != 11 && argc != 12)
        return usage(argv[0]);
    if (argc >= 11) {
        v0 = atof(argv[6]);
        kappa = atof(argv[7]);
        theta = atof(argv[8]);
        xi = atof(argv[9]);
        rho = atof(argv[10]);
    }
    if (argc == 7)
        maxiter = atoi(argv[6]);
    else if (argc == 12)
        maxiter = atoi(argv[11]);
    if (!params_ok(S, v0, kappa, theta, xi, rho) || maxiter < 1) {
        fprintf(stderr, "Invalid fit input: require S>0, v0>0, kappa>0, theta>0, xi>0, |rho|<1, maxiter>=1\n");
        return EXIT_FAILURE;
    }

    surface surf;
    if (surface_load(&surf, path) != 0) {
        fprintf(stderr, "Cannot load surface from %s\n", path);
        return EXIT_FAILURE;
    }

    double x[5] = {
        log(v0), log(kappa), log(theta), log(xi),
        atanh(clampd(rho, -0.999, 0.999))
    };
    cal_ctx ctx = { &surf, S, r, q };
    double fval = nm_minimize(objective, &ctx, x, 5, maxiter);

    double v0f, kappaf, thetaf, xif, rhof;
    unpack(x, &v0f, &kappaf, &thetaf, &xif, &rhof);

    int npts = surf.ns * surf.nt;
    double sse = 0.0, maxe = 0.0;
    printf("Calibration of %s (%d strikes x %d tenors), S=%.6g r=%.6g q=%.6g\n",
           path, surf.ns, surf.nt, S, r, q);
    printf("v0=%.6f kappa=%.6f theta=%.6f xi=%.6f rho=%.6f\n",
           v0f, kappaf, thetaf, xif, rhof);
    printf("mean squared vol error %.3e\n\n", fval);
    printf("%-8s %-8s %10s %10s %10s\n", "krel", "T", "mkt_vol", "fit_vol", "err");
    for (int i = 0; i < surf.ns; i++) {
        for (int j = 0; j < surf.nt; j++) {
            double T = surf.tenor[j];
            double F = S * exp((r - q) * T);
            double K = surf.krel[i] * F;
            double price = heston_call(S, K, T, r, q, v0f, kappaf, thetaf, xif, rhof);
            double lo = fmax(S * exp(-q * T) - K * exp(-r * T), 0.0);
            double cap = S * exp(-q * T);
            price = clampd(price, lo, cap);
            double iv = bs_implied_vol(price, S, K, T, r, q);
            double e = iv - surf.vol[i * surf.nt + j];
            sse += e * e;
            if (fabs(e) > maxe)
                maxe = fabs(e);
            printf("%-8.4f %-8.4f %10.6f %10.6f %10.2e\n",
                   surf.krel[i], T, surf.vol[i * surf.nt + j], iv, e);
        }
    }
    printf("\nRMSE %.6f vol points, max |err| %.6f\n", sqrt(sse / npts), maxe);
    printf("\nFeed into the pricer (var-mode arg order v0 kappa theta xi rho):\n");
    printf("  ./black_scholes var <S> <T> <r> <q> <sigma> <Kvar> %.6f %.6f %.6f %.6f %.6f mc <PATHS>\n",
           v0f, kappaf, thetaf, xif, rhof);

    surface_free(&surf);
    return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "price") == 0)
        return mode_price(argc, argv);
    if (argc >= 2 && strcmp(argv[1], "gen") == 0)
        return mode_gen(argc, argv);
    if (argc >= 2 && strcmp(argv[1], "fit") == 0)
        return mode_fit(argc, argv);
    return usage(argv[0]);
}
