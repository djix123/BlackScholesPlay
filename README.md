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
  (full-truncation Euler, daily fixings, antithetic variates) when
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

`gen` writes a synthetic surface for given Heston params (useful for
testing); `fit` prints the calibrated parameters, a per-point error table,
and the ready-to-use `black_scholes var ... v0 kappa theta xi rho` arguments.

Example:

    ./heston_calibrate gen surf.txt 100 0.05 0.01 0.04 2.0 0.05 0.3 -0.7
    ./heston_calibrate fit surf.txt 100 0.05 0.01 0.02 1.0 0.03 0.5 -0.3 3000
