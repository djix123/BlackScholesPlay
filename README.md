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
Vega is per 1.00 of vol; theta and rho are per year.

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
testing); `fit` prints the calibrated parameters, a per-point error table,
and the ready-to-use `black_scholes var ... v0 kappa theta xi rho` arguments.

Example:

    ./heston_calibrate gen surf.txt 100 0.05 0.01 0.04 2.0 0.05 0.3 -0.7
    ./heston_calibrate fit surf.txt 100 0.05 0.01 0.02 1.0 0.03 0.5 -0.3 3000

## Testing

Build both tools:

    gcc -O2 -o black_scholes black_scholes.c -lm
    gcc -O2 -o heston_calibrate heston_calibrate.c -lm

Analytic vs Monte-Carlo (should agree within the reported standard error):

    ./black_scholes 100 100 1 0.05 0 0.2 mc 1000000

Variance mode: with a flat smile the model-free swap strike must equal
`sigma^2`, and the MC implied variance must match the analytic Heston `E[V]`
(0.04 here):

    ./black_scholes var 100 1 0.05 0 0.2 0.04 0.04 2.0 0.04 0.3 -0.7 mc 1000000

BS limit (`xi=0`, `v0=theta=sigma^2`): the variance call must approach
`e^{-rT} max(sigma^2 - Kvar, 0)` = 0.009512:

    ./black_scholes var 100 1 0.05 0 0.2 0.03 0.04 2.0 0.04 0.0 0.0 mc 200000

Heston Fourier pricer in isolation (expect call ~10.2295, implied vol
~0.2107; cross-checked against Heston Monte-Carlo to within its SE):

    ./heston_calibrate price 100 100 1 0.05 0.01 0.04 2.0 0.05 0.3 -0.7

Calibration round-trip with `surf.txt`, a synthetic 9x6 surface (relative
strikes 0.70-1.30, tenors 1m-3y) generated under Heston with
`v0=0.04 kappa=2 theta=0.05 xi=0.3 rho=-0.7` (S=100, r=0.05, q=0.02).
Regenerate it, then verify the fit recovers those parameters from a
perturbed start (expect RMSE ~1e-8 vol points):

    ./heston_calibrate gen surf.txt 100 0.05 0.02 0.04 2.0 0.05 0.3 -0.7 9 6 0.7 1.3 0.0833333333 3.0
    ./heston_calibrate fit surf.txt 100 0.05 0.02 0.02 1.0 0.03 0.5 -0.3 3000

## Regression checks

    ./regress.sh

Builds both tools into a temp dir and runs all of the checks above plus
regressions for previously fixed bugs, exiting non-zero on any failure
(Monte-Carlo checks use 5-sigma tolerances, so they are not flaky despite the
time-based seed):

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
