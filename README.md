# BlackScholesPlay

A Black-Scholes option pricer in C with two pricing methods:

- **Analytic**: closed-form Black-Scholes-Merton prices and greeks
  (Delta, Gamma, Vega, Theta, Rho), supporting dividend yield.
- **Monte-Carlo**: risk-neutral GBM terminal simulation with
  Box-Muller normals and antithetic variates, reporting the estimate
  with its standard error.

It also prices **variance options** under Heston stochastic volatility
(`var` mode, see below) and ships a **calibration tool**
(`heston_calibrate`) that fits Heston parameters to an implied-vol
surface quoted by relative strike and tenor.

## Build

    gcc -O2 -o black_scholes black_scholes.c -lm

## Usage

    ./black_scholes S K T r q sigma [mc PATHS]

`sigma = 0` prices the deterministic-forward limit
`max(+-(S e^{-qT} - K e^{-rT}), 0)`; `T = 0` gives intrinsic value.
Vega is per 1.00 of vol; theta and rho are per year. Every numeric argument
is parsed strictly (`1x` or `abc` is an error, not silently `0`).

Monte-Carlo runs are seeded from the time, pid and clock; set `BS_SEED=<n>`
(an unsigned integer) for a reproducible run, e.g.
`BS_SEED=42 ./black_scholes 100 100 1 0.05 0 0.2 mc 1000000`.

Example:

    ./black_scholes 100 100 1 0.05 0 0.2 mc 1000000

## Variance mode

Prices annualized variance options (options on realized variance) and the
model-free variance-swap strike:

    ./black_scholes var S T r q sigma Kvar v0 kappa theta xi rho [mc PATHS]

- `sigma`: flat implied vol used for the model-free variance-swap
  replication, `Kvar = (2/T) e^{rT} [int P/K^2 + int C/K^2]` over a strip of
  vanilla options (with a flat smile this equals `sigma^2`).
- `v0 kappa theta xi rho`: Heston stochastic-vol parameters (initial
  variance, mean-reversion speed, long-run variance, vol of vol,
  correlation). Variance options are priced by Monte-Carlo under Heston
  (full-truncation Euler, N = round(252 T) steps, antithetic variates;
  realized variance is `V = (1/T) sum (ln S_i/S_{i-1})^2`) when
  `mc PATHS` is given; the run also reports the analytic Heston
  `E[V] = theta + (v0-theta)(1-e^{-kappa T})/(kappa T)` for comparison.
  The run also prints the Feller ratio `2 kappa theta / xi^2` (variance stays
  positive when it is >= 1; full truncation handles the other case).
- The replication integrates on an asymmetric log-strike grid running from
  `ln F - sigma^2 T/2 - 8 sigma sqrt(T)` to `ln F + 8 sigma sqrt(T)`. With
  `s = sigma sqrt(T)`, `ln K` under the risk-neutral measure is normal with
  mean `ln F - s^2/2` and standard deviation `s`, so that window covers the
  whole distribution however large `s` gets (the half-width is kept at 4 or
  more whenever the spacing makes it affordable). Spacing is
  `<= min(0.001, s/100)` with no floor, and the puts and the calls are
  integrated as two grids that both end at `K = F`, where the integrand has a
  kink. If the grid would need more than 5e7 points, or a strike limit would
  overflow `exp()`, the run reports `Failed to compute variance swap strike`
  instead of a wrong number.

Example:

    ./black_scholes var 100 1 0.05 0 0.2 0.04 0.04 2.0 0.04 0.3 -0.7 mc 1000000

## Calibration tool: heston_calibrate

Calibrates Heston parameters (v0, kappa, theta, xi, rho) to an implied-vol
surface quoted by relative strike and tenor, using a semi-analytic Heston
pricer (characteristic function, "little trap" formulation, Gauss-Legendre
quadrature) and Nelder-Mead on the mean squared implied-vol error.

Build:

    gcc -O2 -o heston_calibrate heston_calibrate.c -lm

Modes:

    ./heston_calibrate price S K T r q v0 kappa theta xi rho
    ./heston_calibrate mc S K T r q v0 kappa theta xi rho PATHS
    ./heston_calibrate gen OUTFILE S r q v0 kappa theta xi rho [ns nt kmin kmax tmin tmax]
    ./heston_calibrate fit SURFFILE S r q [v0 kappa theta xi rho] [maxiter]

Surface file format (`#` comment lines allowed), i.e. `impliedvols[strike][tenor]`:

    ns nt
    krel_1 ... krel_ns     # relative strikes k = K/F
    t_1 ... t_nt           # tenors in years
    vol row 1 (nt values) ... vol row ns

For `fit`, the five start parameters are all-or-nothing (default start
`0.04 2.0 0.04 0.3 -0.7`), and `maxiter` (default 2000) may follow either
form; any other argument count prints usage and exits non-zero.

`gen` writes a synthetic surface for given Heston params (useful for
testing); `fit` prints the calibrated parameters, the Feller ratio, a
per-point error table, and the ready-to-use
`black_scholes var ... v0 kappa theta xi rho` arguments. `mc` prices one
vanilla by Heston Monte-Carlo (full-truncation Euler, antithetic) and reports
the difference from the Fourier price in standard errors, as an independent
check of the semi-analytic pricer; it honours `BS_SEED`.

Implementation notes:

- The Fourier integral uses composite Gauss-Legendre panels and continues
  until the integrand has died out, so short tenors (which decay slowly) are
  priced accurately; the characteristic function is evaluated once per tenor
  and shared across all strikes.
- Implied vols are inverted from the out-of-the-money option (call for
  `K >= F`, otherwise the put via parity), which stays well conditioned in
  the wings. A model price outside the invertible range gives no vol: it is
  shown as `n/a`, penalised in the fit objective (not clamped), and makes
  `gen` refuse to write the surface.
- `fit` warns if a parameter ends on its search bound.
- `ns = 1` (`nt = 1`) in `gen` uses `kmin` (`tmin`) only, so an explicit
  `kmax` (`tmax`) must equal it.
- Surface files are validated: `1 <= ns, nt <= 10000`, `ns*nt <= 10^6`,
  strikes, tenors and vols must be finite and positive, and trailing data
  after the last vol is an error. Command-line numbers are parsed strictly.

Example:

    ./heston_calibrate gen surf.txt 100 0.05 0.01 0.04 2.0 0.05 0.3 -0.7
    ./heston_calibrate fit surf.txt 100 0.05 0.01 0.02 1.0 0.03 0.5 -0.3 3000

## Testing

Build both tools:

    gcc -O2 -o black_scholes black_scholes.c -lm
    gcc -O2 -o heston_calibrate heston_calibrate.c -lm

Analytic vs Monte-Carlo (expect agreement within about two standard errors;
a single run exceeds one SE roughly a third of the time):

    ./black_scholes 100 100 1 0.05 0 0.2 mc 1000000

Variance mode: with a flat smile the model-free swap strike must equal
`sigma^2`, and the MC implied variance must match the analytic Heston `E[V]`
(0.04 here):

    ./black_scholes var 100 1 0.05 0 0.2 0.04 0.04 2.0 0.04 0.3 -0.7 mc 1000000

BS limit (`xi=0`, `v0=theta=sigma^2`): the variance call must approach
`e^{-rT} max(sigma^2 - Kvar, 0)` = 0.009512:

    ./black_scholes var 100 1 0.05 0 0.2 0.03 0.04 2.0 0.04 0.0 0.0 mc 200000

Heston Fourier pricer in isolation (expect call ~10.2295, implied vol
~0.2107):

    ./heston_calibrate price 100 100 1 0.05 0.01 0.04 2.0 0.05 0.3 -0.7

Cross-check against Heston Monte-Carlo (expect a difference of a couple of
SE at most; the Euler time-stepping adds a small bias of a few 1e-3 in price):

    ./heston_calibrate mc 100 100 1 0.05 0.01 0.04 2.0 0.05 0.3 -0.7 1000000

Short tenors need the adaptive Fourier cutoff (expect call 0.442926, implied
vol 0.199964; a fixed cutoff of 200 gave 0.1943):

    ./heston_calibrate price 100 100 0.003 0.05 0.01 0.04 2.0 0.05 0.5 -0.7

Calibration round-trip with `surf.txt`, a synthetic 9x6 surface (relative
strikes 0.70-1.30, tenors 1m-3y) generated under Heston with
`v0=0.04 kappa=2 theta=0.05 xi=0.3 rho=-0.7` (S=100, r=0.05, q=0.02).
Regenerate it, then verify the fit recovers those parameters from a
perturbed start (expect RMSE ~5e-9 vol points, Feller ratio 2.22):

    ./heston_calibrate gen surf.txt 100 0.05 0.02 0.04 2.0 0.05 0.3 -0.7 9 6 0.7 1.3 0.0833333333 3.0
    ./heston_calibrate fit surf.txt 100 0.05 0.02 0.02 1.0 0.03 0.5 -0.3 3000

## Regression checks

    ./regress.sh

Builds both tools into a temp dir and runs all of the checks above plus
regressions for previously fixed bugs, exiting non-zero on any failure
(Monte-Carlo checks use 5-sigma tolerances, so they are not flaky despite the
time-based seed; set `BS_SEED` to reproduce a run):

- `sigma = 0` European prices have no intrinsic-value floor (put 100/105 at
  r=5% is 0, not 5; call with q>r is 0.4837, not 10), and the variance-swap
  strike with `sigma = 0` is 0. `T = 0` is still intrinsic.
- Variance MC `E[V]` matches the analytic Heston value when `252 T` is not an
  integer (T = 0.3, 0.55, 2.37); annualizing by `252/N` used to bias it low.
- `fit` does not depend on uninitialised heap memory: it is rebuilt with
  `malloc` pre-filled with `0x00`, `0x7f` and `0xff` and must converge to the
  same parameters each time (the Nelder-Mead initial simplex used to leave one
  vertex uninitialised and never perturb `v0`).
- `fit` honours `maxiter` when given without start parameters, and rejects
  partial start-parameter lists and trailing arguments instead of ignoring them.
- The variance-swap strike equals `sigma^2` for large `sigma sqrt(T)`
  (`sigma=1, T=5` used to give 0.893 because the strike grid was cut at fixed
  multiples of the forward, and the old grid symmetric about `ln F` missed the
  `-s^2/2` shift of the log-price, giving 3.216 for `sigma=2, T=100`) and for
  small ones (`sigma=0.05, T=0.1`, and `sigma=0.001, T=0.001` where the old
  `1e-4` floor on the spacing gave 0.000003 instead of 0.000001); grids that
  would need more than 5e7 points (`sigma=10, T=100` and `sigma=1e5, T=1e4`)
  are refused with an error instead of overflowing the int point count.
- The Heston Fourier price is accurate at short and long tenors (T = 0.003 and
  T = 10 against high-resolution references) and agrees with an independent
  Heston Monte-Carlo (`heston_calibrate mc`) within 5 SE.
- The committed `surf.txt` regenerates byte-for-byte and round-trips through
  `fit` with RMSE < 1e-6 vol points, with no `n/a` points or warnings.
- Surface files with truncated, non-positive, non-finite or trailing data are
  rejected; `gen` refuses to write a surface with non-invertible points (and
  writes no file); non-numeric command-line arguments are rejected; the `var`
  usage message names the program.
- `BS_SEED` makes `black_scholes mc`, `black_scholes var ... mc` and
  `heston_calibrate mc` reproducible, and different seeds give different runs.
