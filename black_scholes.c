#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double norm_cdf(double x)
{
    return 0.5 * erfc(-x / sqrt(2.0));
}

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

    unsigned short seed[3] = {0x1234, 0x5678, (unsigned short)(time(NULL) & 0xffff)};
    seed48(seed);

    double drift = (r - q - 0.5 * sigma * sigma) * T;
    double vol = sigma * sqrt(T);
    double disc = exp(-r * T);

    double call_sum = 0.0, put_sum = 0.0;
    double call_sq_sum = 0.0, put_sq_sum = 0.0;
    long n = 0;

    for (long i = 0; i < paths; i++) {
        double u1 = drand48();
        double u2 = drand48();
        while (u1 <= 0.0)
            u1 = drand48();
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);

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

/* Standard normal via Box-Muller on drand48(). */
static double bm_rand(void)
{
    double u1 = drand48();
    double u2 = drand48();
    while (u1 <= 0.0)
        u1 = drand48();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
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
 * Integrates on a log-spaced strike grid, trapezoidal in ln K.
 * Returns 0 on success, -1 on invalid input.
 */
int variance_swap_strike(double *out, double S, double T, double r, double q,
                         double (*smile)(double, void *), void *smile_ctx,
                         int n_points)
{
    if (!out || !smile || S <= 0.0 || T <= 0.0 || n_points < 10)
        return -1;

    double F = S * exp((r - q) * T);
    double lmin = log(0.02 * F), lmax = log(5.0 * F);
    double dx = (lmax - lmin) / (n_points - 1);

    double integral = 0.0, prev = 0.0;
    for (int i = 0; i < n_points; i++) {
        double K = exp(lmin + i * dx);
        bs_result res;
        if (black_scholes(&res, S, K, T, r, q, smile(K, smile_ctx), K > F) != 0)
            return -1;
        double g = res.price / K; /* price/K^2 * dK, with dK = K dx */
        if (i > 0)
            integral += 0.5 * (prev + g) * dx;
        prev = g;
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
 * Antithetic variates. Returns 0 on success, -1 on invalid input.
 */
int heston_variance_mc(var_mc_result *out, double S, double T, double r, double q,
                       double Kvar, double v0, double kappa, double theta,
                       double xi, double rho, long paths)
{
    if (!out || S <= 0.0 || T <= 0.0 || Kvar < 0.0 || v0 < 0.0 ||
        kappa <= 0.0 || theta < 0.0 || xi < 0.0 || fabs(rho) > 1.0 || paths <= 0)
        return -1;

    unsigned short seed[3] = {0x1234, 0x5678, (unsigned short)(time(NULL) & 0xffff)};
    seed48(seed);

    int N = (int)(252.0 * T + 0.5);
    if (N < 1)
        N = 1;
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
 */
static int main_variance(int argc, char **argv)
{
    int run_mc = 0;
    long paths = 0;

    if (argc == 14 && strcmp(argv[12], "mc") == 0) {
        run_mc = 1;
        paths = atol(argv[13]);
    } else if (argc != 12) {
        fprintf(stderr,
                "Usage: %s S T r q sigma Kvar v0 kappa theta xi rho [mc PATHS]\n"
                "  S: spot, T: years to expiry, r: risk-free rate, q: dividend yield\n"
                "  sigma: flat vol for the model-free variance-swap replication\n"
                "  Kvar: variance strike (annualized variance, e.g. 0.04 = 20%% vol)\n"
                "  v0: initial variance, kappa: mean-reversion speed, theta: long-run variance\n"
                "  xi: vol of vol, rho: corr(S,v) in [-1,1]\n"
                "  optional 'mc PATHS': also price variance options by Monte-Carlo\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    double S = atof(argv[1]), T = atof(argv[2]), r = atof(argv[3]), q = atof(argv[4]);
    double sigma = atof(argv[5]), Kvar = atof(argv[6]);
    double v0 = atof(argv[7]), kappa = atof(argv[8]), theta = atof(argv[9]);
    double xi = atof(argv[10]), rho = atof(argv[11]);

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
    printf("Heston: v0=%.4f kappa=%.4f theta=%.4f xi=%.4f rho=%.4f\n\n",
           v0, kappa, theta, xi, rho);
    printf("Model-free variance swap strike (flat smile sigma=%.4f): %.6f  (sigma^2=%.6f)\n",
           sigma, swap_strike, sigma * sigma);
    printf("Heston expected annualized variance E[V]: %.6f\n", exp_var);
    printf("Deterministic-BS variance option (V=sigma^2): call %.6f, put %.6f\n",
           det_call, det_put);

    if (run_mc) {
        var_mc_result mc;
        if (heston_variance_mc(&mc, S, T, r, q, Kvar, v0, kappa, theta, xi, rho, paths) != 0) {
            fprintf(stderr, "Invalid Monte-Carlo input: require PATHS > 0\n");
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
        return main_variance(argc - 1, argv + 1);

    int run_mc = 0;
    long paths = 0;

    if (argc == 9 && strcmp(argv[7], "mc") == 0) {
        run_mc = 1;
        paths = atol(argv[8]);
    } else if (argc != 7) {
        fprintf(stderr,
                "Usage: %s S K T r q sigma [mc PATHS]\n"
                "  S: spot price, K: strike, T: years to expiry,\n"
                "  r: risk-free rate, q: dividend yield, sigma: volatility,\n"
                "  optional 'mc PATHS': also run Monte-Carlo with PATHS paths\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    bs_result call, put;
    double S = atof(argv[1]), K = atof(argv[2]), T = atof(argv[3]);
    double r = atof(argv[4]), q = atof(argv[5]), sigma = atof(argv[6]);

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
