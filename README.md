# BlackScholesPlay

A Black-Scholes option pricer in C with two pricing methods:

- **Analytic**: closed-form Black-Scholes-Merton prices and greeks
  (Delta, Gamma, Vega, Theta, Rho), supporting dividend yield.
- **Monte-Carlo**: risk-neutral GBM terminal simulation with
  Box-Muller normals and antithetic variates, reporting the estimate
  with its standard error.

## Build

    gcc -O2 -o black_scholes black_scholes.c -lm

## Usage

    ./black_scholes S K T r q sigma [mc PATHS]

Example:

    ./black_scholes 100 100 1 0.05 0 0.2 mc 1000000
