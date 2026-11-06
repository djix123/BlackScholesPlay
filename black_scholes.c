#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

typedef struct {
    double price;
    double delta;
    double gamma;
    double vega;
    double theta;
    double rho;
} bs_result;

typedef struct {
    double call_price;
    double put_price;
    double call_std_err;
    double put_std_err;
} mc_result;

typedef struct {
    double call_price;
    double put_price;
    double call_std_err;
    double put_std_err;
    double mean_variance;
    double mean_variance_std_err;
} var_mc_result;

/* Call (is_call=1) or put (is_call=0) price and greeks.
 * S: spot, K: strike, T: time to expiry in years,
 * r: risk-free rate, q: dividend yield, sigma: volatility.
 * Returns 0 on success, -1 on invalid input.
 */
int black_scholes(bs_result *out, double S, double K, double T,
                  double r, double q, double sigma, int is_call)
{
    if (!out || S <= 0.0 || K <= 0.0 || T < 0.0 || sigma < 0.0)
        return -1;

    if (T == 0.0) {
        /* At expiry the value is intrinsic. */
        out->price = is_call ? fmax(S - K, 0.0) : fmax(K - S, 0.0);
        out->delta = out->gamma = out->vega = out->theta = out->rho = 0.0;
        if (is_call && S > K)
            out->delta = 1.0;
        if (!is_call && S < K)
            out->delta = -1.0;
        return 0;
    }
    if (sigma == 0.0) {
        /* Deterministic forward: the option is worth the discounted payoff at
         * the forward, max(+-(S e^{-qT} - K e^{-rT}), 0). Intrinsic value
         * S-K is NOT a valid floor for a European option when r or q != 0. */
        double disc_q = exp(-q * T), disc_r = exp(-r * T);
        double fwd_val = S * disc_q - K * disc_r; /* PV of (F - K) */
        double sgn = is_call ? 1.0 : -1.0;
        int itm = sgn * fwd_val > 0.0;
        out->price = itm ? sgn * fwd_val : 0.0;
        out->delta = itm ? sgn * disc_q : 0.0;
        out->gamma = out->vega = 0.0;
        out->theta = itm ? sgn * (q * S * disc_q - r * K * disc_r) : 0.0;
        out->rho = itm ? sgn * K * T * disc_r : 0.0;
        return 0;
    }

    double sqrtT = sqrt(T);
    double d1 = (log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    double d2 = d1 - sigma * sqrtT;
    double disc_q = exp(-q * T);
    double disc_r = exp(-r * T);
    double pdf_d1 = exp(-0.5 * d1 * d1) / sqrt(2.0 * M_PI);
    double Nd1 = norm_cdf(d1), Nd2 = norm_cdf(d2);
    double Nmd1 = norm_cdf(-d1), Nmd2 = norm_cdf(-d2);

    if (is_call) {
        out->price = S * disc_q * Nd1 - K * disc_r * Nd2;
        out->delta = disc_q * Nd1;
        out->theta = -S * disc_q * pdf_d1 * sigma / (2.0 * sqrtT)
                    - r * K * disc_r * Nd2 + q * S * disc_q * Nd1;
        out->rho = K * T * disc_r * Nd2;
    } else {
        out->price = K * disc_r * Nmd2 - S * disc_q * Nmd1;
        out->delta = -disc_q * Nmd1;
        out->theta = -S * disc_q * pdf_d1 * sigma / (2.0 * sqrtT)
                    + r * K * disc_r * Nmd2 - q * S * disc_q * Nmd1;
        out->rho = -K * T * disc_r * Nmd2;
    }

    out->gamma = disc_q * pdf_d1 / (S * sigma * sqrtT);
    out->vega = S * disc_q * pdf_d1 * sqrtT;

    return 0;
}

/* Monte-Carlo pricing under risk-neutral GBM:
 *   S_T = S * exp((r - q - sigma^2/2) T + sigma sqrt(T) Z),  Z ~ N(0,1)
 * Exact terminal simulation (single step), antithetic variates.
 * Returns 0 on success, -1 on invalid input.
 */
int black_scholes_mc(mc_result *out, double S, double K, double T,
                     double r, double q, double sigma, long paths)
{
    if (!out || S <= 0.0 || K <= 0.0 || T < 0.0 || sigma < 0.0 || paths <= 0)
        return -1;

    rng_seed_once();

    double drift = (r - q - 0.5 * sigma * sigma) * T;
    double vol = sigma * sqrt(T);
    double disc = exp(-r * T);

    double call_sum = 0.0, put_sum = 0.0;
    double call_sq_sum = 0.0, put_sq_sum = 0.0;
    long n = 0;

    for (long i = 0; i < paths; i++) {
        double z = bm_rand();

        double st_plus = S * exp(drift + vol * z);
        double st_minus = S * exp(drift - vol * z);

        double cp = fmax(st_plus - K, 0.0) * disc;
        double cm = fmax(st_minus - K, 0.0) * disc;
        double pp = fmax(K - st_plus, 0.0) * disc;
        double pm = fmax(K - st_minus, 0.0) * disc;

        call_sum += (cp + cm) * 0.5;
        put_sum += (pp + pm) * 0.5;
        call_sq_sum += (cp + cm) * (cp + cm) * 0.25;
        put_sq_sum += (pp + pm) * (pp + pm) * 0.25;
        n++;
    }

    double call_mean = call_sum / n;
    double put_mean = put_sum / n;
    double call_var = call_sq_sum / n - call_mean * call_mean;
    double put_var = put_sq_sum / n - put_mean * put_mean;
    if (call_var < 0.0) call_var = 0.0;
    if (put_var < 0.0) put_var = 0.0;

    out->call_price = call_mean;
    out->put_price = put_mean;
    out->call_std_err = sqrt(call_var / n);
    out->put_std_err = sqrt(put_var / n);
    return 0;
}

/* Flat implied-vol smile: every strike gets the same volatility.
 * Swap in any sigma(K) here to replicate the swap off a real smile.
 */
static double flat_smile(double K, void *ctx)
{
    (void)K;
    return *(const double *)ctx;
}

/* Model-free fair strike of a continuously-monitored annualized variance
 * swap, statically replicated from OTM vanillas (Carr-Madan):
 *   K_var = (2/T) e^{rT} [ int_0^F P(K)/K^2 dK + int_F^inf C(K)/K^2 dK ],
 *   F = S e^{(r-q)T}.
 * Integrates on a log-strike grid, trapezoidal in ln K. ln K is normal with
 * mean ln F - s^2/2 and standard deviation s = smile(F) sqrt(T), so the grid
 * is asymmetric: [ln F - s^2/2 - 8 s, ln F + 8 s]. That keeps both wings
 * covered (to ~e^-32 of Gaussian tail mass) however large s gets; the strip
 * is still widened to a half-width of 4 whenever the spacing makes that cheap.
 * The puts (K < F) and the calls (K > F) are integrated as two separate
 * segments, each with its own uniform spacing, so ln F is always an endpoint:
 * the integrand has a kink there and trapezoid must not straddle it. Spacing
 * is <= min(0.001, s/100) with no floor; n_points is only a minimum. With
 * s = 0 the log-price is a point mass at F and the strike is exactly 0. The
 * count is computed in double; it returns -1 when the grid would need more
 * than VSS_MAX_POINTS points or when exp() of a grid limit would overflow or
 * underflow.
 * Returns 0 on success, -1 on invalid or unresolvable input.
 */
int variance_swap_strike(double *out, double S, double T, double r, double q,
                         double (*smile)(double, void *), void *smile_ctx,
                         int n_points)
{
    if (!out || !smile || S <= 0.0 || T <= 0.0 || n_points < 10)
        return -1;

    /* Refuse grids this function could not integrate, and grid limits whose
     * exponentials are not representable as doubles (exp() overflows or
     * flushes to zero beyond about +-700 in its argument). */
    const double VSS_MAX_POINTS = 5e7;   /* hard limit, far below INT_MAX */
    const double VSS_WIDE_POINTS = 4e6;  /* budget for widening to +/-4 */
    const double VSS_MAX_LOGK = 700.0;   /* max |ln K| on the grid */

    double lF = log(S) + (r - q) * T;
    if (!(isfinite(lF) && fabs(lF) <= VSS_MAX_LOGK))
        return -1;
    double F = exp(lF);
    double s0 = smile(F, smile_ctx) * sqrt(T);
    if (!(s0 >= 0.0) || !isfinite(s0))
        return -1;

    if (s0 == 0.0) {
        /* Zero-width grid: the log-price is a point mass at F, every OTM
         * vanilla is worthless and the replication is exactly 0. Integrating
         * a strip of arbitrary width instead would leave rounding noise of
         * order 1e-19 from the two K = F endpoints. */
        *out = 0.0;
        return 0;
    }
    /* Asymmetric window around ln F, wide enough for the whole distribution. */
    double lo_off = 0.5 * s0 * s0 + 8.0 * s0; /* ln F down to the put wing */
    double hi_off = 8.0 * s0;                 /* ln F up to the call wing */
    double dx_max = fmin(0.001, s0 / 100.0);
    if (!(dx_max > 0.0) || !isfinite(dx_max))
        return -1;

    /* A half-width of at least 4 is kept only while the s/100 spacing makes it
     * affordable; for tiny s it would cost tens of millions of points. */
    if (hi_off < 4.0 && (fmax(lo_off, 4.0) + 4.0) / dx_max <= VSS_WIDE_POINTS) {
        lo_off = fmax(lo_off, 4.0);
        hi_off = 4.0;
    }

    double lmin = lF - lo_off, lmax = lF + hi_off;
    if (!(lmin > -VSS_MAX_LOGK) || !(lmax < VSS_MAX_LOGK))
        return -1;

    double n_put = ceil(lo_off / dx_max) + 1.0;  /* points on [lmin, ln F] */
    double n_call = ceil(hi_off / dx_max) + 1.0; /* points on [ln F, lmax] */
    double n_total = n_put + n_call - 1.0;       /* ln F belongs to both wings */
    if (!(n_total <= VSS_MAX_POINTS))
        return -1;
    if ((double)n_points < n_put)
        n_points = (int)n_put;
    if ((double)n_points < n_call)
        n_points = (int)n_call;

    /* Trapezoid in ln K over the two wings separately, each with its own
     * uniform spacing and both reaching exactly K = F, where the integrand
     * has a kink. */
    double integral = 0.0;
    for (int seg = 0; seg < 2; seg++) {
        int is_call = seg;
        double a = is_call ? lF : lmin, b = is_call ? lmax : lF;
        double dx = (b - a) / (double)(n_points - 1);
        double prev = 0.0;
        for (int i = 0; i < n_points; i++) {
            double lk = (i == n_points - 1) ? b : a + (double)i * dx;
            double K = exp(lk);
            bs_result res;
            if (black_scholes(&res, S, K, T, r, q, smile(K, smile_ctx), is_call) != 0)
                return -1;
            double g = res.price / K; /* price/K^2 * dK, with dK = K dx */
            if (i > 0)
                integral += 0.5 * (prev + g) * dx;
            prev = g;
        }
    }

    *out = (2.0 / T) * exp(r * T) * integral;
    return 0;
}

/* E[annualized variance] under Heston, continuous monitoring:
 * E[int_0^T v dt]/T = theta + (v0 - theta)(1 - e^{-kappa T})/(kappa T).
 */
static double heston_expected_variance(double v0, double kappa, double theta, double T)
{
    double x = kappa * T;
    double ratio = fabs(x) < 1e-12 ? 1.0 : (1.0 - exp(-x)) / x;
    return theta + (v0 - theta) * ratio;
}

/* Monte-Carlo pricing of annualized variance options under Heston:
 *   dS/S = (r-q) dt + sqrt(v) dW1
 *   dv   = kappa(theta - v) dt + xi sqrt(v) dW2,   dW1 dW2 = rho dt
 * Full-truncation Euler with N = round(252 T) steps (daily fixings).
 * Realized variance per path: V = (1/T) sum (ln S_i/S_{i-1})^2, i.e. annualized
 * by the actual elapsed time so that E[V] matches the analytic Heston E[V]
 * for any T (annualizing by 252/N would bias V when 252 T is not an integer).
 * Antithetic variates. Returns 0 on success, -1 on invalid input or when N
 * would exceed MC_MAX_STEPS.
 */
int heston_variance_mc(var_mc_result *out, double S, double T, double r, double q,
                       double Kvar, double v0, double kappa, double theta,
                       double xi, double rho, long paths)
{
    if (!out || S <= 0.0 || T <= 0.0 || Kvar < 0.0 || v0 < 0.0 ||
        kappa <= 0.0 || theta < 0.0 || xi < 0.0 || fabs(rho) > 1.0 || paths <= 0)
        return -1;

    rng_seed_once();

    double steps = floor(252.0 * T + 0.5);
    if (steps > MC_MAX_STEPS)
        return -1;
    int N = steps < 1.0 ? 1 : (int)steps;
    double dt = T / N;
    double ann = 1.0 / T;
    double sdt = sqrt(dt);
    double disc = exp(-r * T);
    double sq1mr2 = sqrt(fmax(1.0 - rho * rho, 0.0));
    double drift = r - q;

    double call_sum = 0.0, put_sum = 0.0;
    double call_sq_sum = 0.0, put_sq_sum = 0.0;
    double var_sum = 0.0, var_sq_sum = 0.0;
    long n = 0;

    for (long i = 0; i < paths; i++) {
        double vp = v0, vm = v0;
        double accp = 0.0, accm = 0.0;

        for (int j = 0; j < N; j++) {
            double z1 = bm_rand();
            double z2 = bm_rand();
            double dw_s = sq1mr2 * z1 + rho * z2;

            double up = fmax(vp, 0.0);
            double dlnp = (drift - 0.5 * up) * dt + sqrt(up) * sdt * dw_s;
            vp += kappa * (theta - up) * dt + xi * sqrt(up) * sdt * z2;
            accp += dlnp * dlnp;

            double um = fmax(vm, 0.0);
            double dlnm = (drift - 0.5 * um) * dt - sqrt(um) * sdt * dw_s;
            vm += kappa * (theta - um) * dt - xi * sqrt(um) * sdt * z2;
            accm += dlnm * dlnm;
        }

        double Vp = ann * accp, Vm = ann * accm;
        double cp = disc * fmax(Vp - Kvar, 0.0);
        double cm = disc * fmax(Vm - Kvar, 0.0);
        double pp = disc * fmax(Kvar - Vp, 0.0);
        double pm = disc * fmax(Kvar - Vm, 0.0);

        double cpair = 0.5 * (cp + cm);
        double ppair = 0.5 * (pp + pm);
        call_sum += cpair;
        put_sum += ppair;
        call_sq_sum += cpair * cpair;
        put_sq_sum += ppair * ppair;

        double vpair = 0.5 * (Vp + Vm);
        var_sum += vpair;
        var_sq_sum += vpair * vpair;
        n++;
    }

    double call_mean = call_sum / n;
    double put_mean = put_sum / n;
    double call_var = call_sq_sum / n - call_mean * call_mean;
    double put_var = put_sq_sum / n - put_mean * put_mean;
    double mvar = var_sum / n;
    double vvar = var_sq_sum / n - mvar * mvar;
    if (call_var < 0.0) call_var = 0.0;
    if (put_var < 0.0) put_var = 0.0;
    if (vvar < 0.0) vvar = 0.0;

    out->call_price = call_mean;
    out->put_price = put_mean;
    out->call_std_err = sqrt(call_var / n);
    out->put_std_err = sqrt(put_var / n);
    out->mean_variance = mvar;
    out->mean_variance_std_err = sqrt(vvar / n);
    return 0;
}

/* Variance mode:
 *   var S T r q sigma Kvar v0 kappa theta xi rho [mc PATHS]
 * argv[0] is "var"; prog is the program name (for messages).
 */
static int main_variance(const char *prog, int argc, char **argv)
{
    int run_mc = 0;
    long paths = 0;

    if (argc == 14 && strcmp(argv[12], "mc") == 0) {
        run_mc = 1;
        if (parse_long(argv[13], &paths) != 0) {
            fprintf(stderr, "Invalid PATHS '%s'\n", argv[13]);
            return EXIT_FAILURE;
        }
    } else if (argc != 12) {
        fprintf(stderr,
                "Usage: %s var S T r q sigma Kvar v0 kappa theta xi rho [mc PATHS]\n"
                "  S: spot, T: years to expiry, r: risk-free rate, q: dividend yield\n"
                "  sigma: flat vol for the model-free variance-swap replication\n"
                "  Kvar: variance strike (annualized variance, e.g. 0.04 = 20%% vol)\n"
                "  v0: initial variance, kappa: mean-reversion speed, theta: long-run variance\n"
                "  xi: vol of vol, rho: corr(S,v) in [-1,1]\n"
                "  optional 'mc PATHS': also price variance options by Monte-Carlo\n"
                "  (set BS_SEED=<n> for a reproducible Monte-Carlo run)\n",
                prog);
        return EXIT_FAILURE;
    }

    static const char *const names[11] = {"S", "T", "r", "q", "sigma", "Kvar",
                                          "v0", "kappa", "theta", "xi", "rho"};
    double a[11];
    if (parse_args(argv + 1, names, 11, a) != 0)
        return EXIT_FAILURE;
    double S = a[0], T = a[1], r = a[2], q = a[3], sigma = a[4], Kvar = a[5];
    double v0 = a[6], kappa = a[7], theta = a[8], xi = a[9], rho = a[10];

    if (S <= 0.0 || T <= 0.0 || sigma < 0.0 || Kvar < 0.0 || v0 < 0.0 ||
        kappa <= 0.0 || theta < 0.0 || xi < 0.0 || fabs(rho) > 1.0) {
        fprintf(stderr,
                "Invalid input: require S>0, T>0, sigma>=0, Kvar>=0, v0>=0, "
                "kappa>0, theta>=0, xi>=0, |rho|<=1\n");
        return EXIT_FAILURE;
    }

    double smile = sigma;
    double swap_strike;
    if (variance_swap_strike(&swap_strike, S, T, r, q, flat_smile, &smile, 4000) != 0) {
        fprintf(stderr, "Failed to compute variance swap strike\n");
        return EXIT_FAILURE;
    }

    double exp_var = heston_expected_variance(v0, kappa, theta, T);
    double disc = exp(-r * T);
    double det_call = disc * fmax(sigma * sigma - Kvar, 0.0);
    double det_put = disc * fmax(Kvar - sigma * sigma, 0.0);

    printf("Variance mode  S=%.4f T=%.4f r=%.4f q=%.4f Kvar=%.4f\n", S, T, r, q, Kvar);
    printf("Heston: v0=%.4f kappa=%.4f theta=%.4f xi=%.4f rho=%.4f\n",
           v0, kappa, theta, xi, rho);
    if (xi > 0.0) {
        double feller = 2.0 * kappa * theta / (xi * xi);
        printf("Feller ratio 2*kappa*theta/xi^2 = %.4f (%s)\n", feller,
               feller >= 1.0 ? "satisfied: variance stays positive"
                             : "violated: variance can hit zero; full truncation is used");
    }
    printf("\n");
    printf("Model-free variance swap strike (flat smile sigma=%.4f): %.6f  (sigma^2=%.6f)\n",
           sigma, swap_strike, sigma * sigma);
    printf("Heston expected annualized variance E[V]: %.6f\n", exp_var);
    printf("Deterministic-BS variance option (V=sigma^2): call %.6f, put %.6f\n",
           det_call, det_put);

    if (run_mc) {
        var_mc_result mc;
        if (heston_variance_mc(&mc, S, T, r, q, Kvar, v0, kappa, theta, xi, rho, paths) != 0) {
            fprintf(stderr, "Invalid Monte-Carlo input: require PATHS > 0 and "
                            "at most 1e7 time steps (252 per year)\n");
            return EXIT_FAILURE;
        }
        printf("\nMonte-Carlo variance options (%ld paths, antithetic, daily fixings):\n", paths);
        printf("Variance call: %.6f +/- %.6f\n", mc.call_price, mc.call_std_err);
        printf("Variance put:  %.6f +/- %.6f\n", mc.put_price, mc.put_std_err);
        printf("MC implied annualized variance: %.6f +/- %.6f (analytic %.6f)\n",
               mc.mean_variance, mc.mean_variance_std_err, exp_var);
    }

    return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "var") == 0)
        return main_variance(argv[0], argc - 1, argv + 1);

    int run_mc = 0;
    long paths = 0;

    if (argc == 9 && strcmp(argv[7], "mc") == 0) {
        run_mc = 1;
        if (parse_long(argv[8], &paths) != 0) {
            fprintf(stderr, "Invalid PATHS '%s'\n", argv[8]);
            return EXIT_FAILURE;
        }
    } else if (argc != 7) {
        fprintf(stderr,
                "Usage: %s S K T r q sigma [mc PATHS]\n"
                "       %s var S T r q sigma Kvar v0 kappa theta xi rho [mc PATHS]\n"
                "  S: spot price, K: strike, T: years to expiry,\n"
                "  r: risk-free rate, q: dividend yield, sigma: volatility,\n"
                "  optional 'mc PATHS': also run Monte-Carlo with PATHS paths\n"
                "  (set BS_SEED=<n> for a reproducible Monte-Carlo run)\n"
                "  Vega is per 1.00 of vol; theta and rho are per year.\n",
                argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    static const char *const names[6] = {"S", "K", "T", "r", "q", "sigma"};
    double a[6];
    if (parse_args(argv + 1, names, 6, a) != 0)
        return EXIT_FAILURE;
    double S = a[0], K = a[1], T = a[2], r = a[3], q = a[4], sigma = a[5];

    bs_result call, put;
    if (black_scholes(&call, S, K, T, r, q, sigma, 1) != 0 ||
        black_scholes(&put, S, K, T, r, q, sigma, 0) != 0) {
        fprintf(stderr, "Invalid input: require S>0, K>0, T>=0, sigma>=0\n");
        return EXIT_FAILURE;
    }

    printf("Black-Scholes  S=%.4f K=%.4f T=%.4f r=%.4f q=%.4f sigma=%.4f\n\n",
           S, K, T, r, q, sigma);
    printf("%-8s %10s %10s %10s %12s %10s\n",
           "Option", "Price", "Delta", "Gamma", "Vega", "Theta");
    printf("Call     %10.6f %10.6f %10.6f %12.6f %10.6f\n",
           call.price, call.delta, call.gamma, call.vega, call.theta);
    printf("Put      %10.6f %10.6f %10.6f %12.6f %10.6f\n",
           put.price, put.delta, put.gamma, put.vega, put.theta);
    printf("\nRho: call %.6f, put %.6f\n", call.rho, put.rho);

    if (run_mc) {
        mc_result mc;
        if (black_scholes_mc(&mc, S, K, T, r, q, sigma, paths) != 0) {
            fprintf(stderr, "Invalid Monte-Carlo input: require PATHS > 0\n");
            return EXIT_FAILURE;
        }
        printf("\nMonte-Carlo (%ld paths, antithetic):\n", paths);
        printf("Call: %.6f +/- %.6f (analytic %.6f)\n",
               mc.call_price, mc.call_std_err, call.price);
        printf("Put:  %.6f +/- %.6f (analytic %.6f)\n",
               mc.put_price, mc.put_std_err, put.price);
    }

    return EXIT_SUCCESS;
}
