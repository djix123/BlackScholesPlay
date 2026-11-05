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

    if (sigma == 0.0 || T == 0.0) {
        double intrinsic = is_call ? fmax(S - K, 0.0) : fmax(K - S, 0.0);
        double disc = exp(-r * T);
        out->price = is_call ? S * exp(-q * T) - K * disc : K * disc - S * exp(-q * T);
        if (out->price < intrinsic)
            out->price = intrinsic;
        out->delta = out->gamma = out->vega = out->theta = out->rho = 0.0;
        if (is_call && S > K)
            out->delta = exp(-q * T);
        if (!is_call && S < K)
            out->delta = -exp(-q * T);
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

int main(int argc, char **argv)
{
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
