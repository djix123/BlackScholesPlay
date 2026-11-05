#!/bin/sh
# Regression checks for black_scholes and heston_calibrate.
# Usage: ./regress.sh        (builds into a temp dir; exits non-zero on failure)
#
# Monte-Carlo checks use 5-sigma tolerances so they are not flaky despite the
# time-based seed.

set -u
cd "$(dirname "$0")" || exit 1
TMP=$(mktemp -d) || exit 1
trap 'rm -rf "$TMP"' EXIT
BS="$TMP/black_scholes"
HC="$TMP/heston_calibrate"
fails=0

gcc -O2 -Wall -Wextra -o "$BS" black_scholes.c -lm || exit 1
gcc -O2 -Wall -Wextra -o "$HC" heston_calibrate.c -lm || exit 1

# check NAME VALUE EXPECTED TOL   (|VALUE-EXPECTED| <= TOL)
check() {
    if awk -v v="$2" -v e="$3" -v t="$4" \
        'BEGIN { if (v !~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/) exit 1; d = v - e; if (d < 0) d = -d; exit !(d <= t) }'; then
        echo "PASS  $1 ($2)"
    else
        echo "FAIL  $1: got '$2', expected $3 +/- $4"
        fails=$((fails + 1))
    fi
}
# check_lt NAME VALUE LIMIT       (VALUE < LIMIT)
check_lt() {
    if awk -v v="$2" -v l="$3" 'BEGIN { if (v !~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/) exit 1; exit !(v < l) }'; then
        echo "PASS  $1 ($2 < $3)"
    else
        echo "FAIL  $1: got '$2', expected < $3"
        fails=$((fails + 1))
    fi
}
# check_gt NAME VALUE LIMIT       (VALUE > LIMIT)
check_gt() {
    if awk -v v="$2" -v l="$3" 'BEGIN { if (v !~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/) exit 1; exit !(v > l) }'; then
        echo "PASS  $1 ($2 > $3)"
    else
        echo "FAIL  $1: got '$2', expected > $3"
        fails=$((fails + 1))
    fi
}
check_fail() { # NAME CMD...   (command must exit non-zero)
    name=$1; shift
    if "$@" >/dev/null 2>&1; then
        echo "FAIL  $name: expected non-zero exit"
        fails=$((fails + 1))
    else
        echo "PASS  $name"
    fi
}

# ---- #2 sigma = 0 European prices (no intrinsic floor) ----------------------
# put K=105, r=5%: max(105 e^-0.05 - 100, 0) = 0 ; call: 100 - 105 e^-0.05
out=$("$BS" 100 105 1 0.05 0 0)
check "sigma=0 call (K=105)"  "$(echo "$out" | awk '$1=="Call"{print $2}')" 0.120910 1e-6
check "sigma=0 put  (K=105)"  "$(echo "$out" | awk '$1=="Put"{print $2}')"  0 1e-9
# call, r=0, q=0.1: max(100 e^-0.1 - 90, 0) = 0.483742 (intrinsic would be 10)
out=$("$BS" 100 90 1 0 0.1 0)
check "sigma=0 call (q>r)"    "$(echo "$out" | awk '$1=="Call"{print $2}')" 0.483742 1e-6
# T = 0 is still intrinsic
out=$("$BS" 100 105 0 0.05 0 0.2)
check "T=0 put intrinsic"     "$(echo "$out" | awk '$1=="Put"{print $2}')"  5 1e-9
# swap strike with sigma = 0 must be 0, not a replication artefact
v=$("$BS" var 100 1 0.05 0 0 0.04 0.04 2 0.04 0.3 -0.7 | awk '/swap strike/{print $8}')
check "swap strike sigma=0"   "$v" 0 1e-9

# ---- baseline: analytic vs MC, flat-smile swap strike -----------------------
out=$("$BS" 100 100 1 0.05 0 0.2 mc 1000000)
c=$(echo "$out" | awk '/^Call:/{print $2, $4, $6}' | tr -d ')')
set -- $c
check "BS call MC vs analytic (5 SE)" "$1" "$3" "$(awk -v s="$2" 'BEGIN{print 5*s}')"
v=$("$BS" var 100 1 0.05 0 0.2 0.04 0.04 2.0 0.04 0.3 -0.7 | awk '/swap strike/{print $8}')
check "swap strike = sigma^2 (flat smile)" "$v" 0.04 1e-6

# ---- #4 MC realized variance matches analytic E[V] for any T ---------------
# For T=0.3 and T=0.55, 252*T is non-integer; annualizing by 252/N instead of
# 1/T biased E[V] by about -2.0e-4 and -1.2e-4 (roughly 18 and 12 SE).
for T in 0.3 0.55 2.37; do
    line=$("$BS" var 100 $T 0.05 0 0.2 0.04 0.04 2.0 0.04 0.3 -0.7 mc 400000 |
           awk '/MC implied/{print $5, $7, $9}' | tr -d ')')
    set -- $line
    check "E[V] MC vs analytic, T=$T (5 SE + 5e-5 discretization)" "${1:-}" "${3:-}" \
          "$(awk -v s="${2:-0}" 'BEGIN{print 5*s + 5e-5}')"
done

# ---- Heston Fourier pricer --------------------------------------------------
out=$("$HC" price 100 100 1 0.05 0.01 0.04 2.0 0.05 0.3 -0.7)
check "Heston call"        "$(echo "$out" | awk '/^Call:/{print $2}')" 10.2295 5e-5
check "Heston implied vol" "$(echo "$out" | awk '/^Call:/{print $7}')" 0.2107  5e-5

# ---- calibration ------------------------------------------------------------
"$HC" gen "$TMP/surf.txt" 100 0.05 0.02 0.04 2.0 0.05 0.3 -0.7 9 6 0.7 1.3 0.0833333333 3.0 >/dev/null
if diff -q "$TMP/surf.txt" surf.txt >/dev/null; then
    echo "PASS  gen reproduces committed surf.txt"
else
    echo "FAIL  gen output differs from surf.txt"; fails=$((fails + 1))
fi

mse() { awk '/mean squared vol error/{print $5}'; }
param() { awk -v k="$1" '/^v0=/{for(i=1;i<=NF;i++){split($i,a,"="); if(a[1]==k) print a[2]}}'; }

# #1 the fit must not depend on uninitialised heap contents. Build variants
# whose malloc() returns memory pre-filled with 0x00, 0x7f and 0xff (NaN
# doubles) and require the same converged fit from each.
for fill in 00 7f ff; do
    cat > "$TMP/poison$fill.h" <<HDR
#include <stdlib.h>
#include <string.h>
static void *poison_malloc(size_t n) { void *p = malloc(n); if (p) memset(p, 0x$fill, n); return p; }
#define malloc poison_malloc
HDR
    gcc -O2 -include "$TMP/poison$fill.h" -o "$TMP/hc_$fill" heston_calibrate.c -lm || exit 1
    out=$("$TMP/hc_$fill" fit surf.txt 100 0.05 0.02 0.02 1.0 0.03 0.5 -0.3 3000)
    check_lt "fit round-trip MSE, malloc filled 0x$fill" "$(echo "$out" | mse)" 1e-12
done
check "fit recovers rho"   "$(echo "$out" | param rho)"   -0.7 1e-4
check "fit recovers theta" "$(echo "$out" | param theta)" 0.05 1e-4

# #1 v0 (coordinate 0) must be searched: start wrong in v0 only.
out=$("$HC" fit surf.txt 100 0.05 0.02 0.02 2.0 0.05 0.3 -0.7 3000)
check_lt "fit from v0-only perturbed start" "$(echo "$out" | mse)" 1e-12
check "  ... recovers v0" "$(echo "$out" | param v0)" 0.04 1e-4

# fit is deterministic
a=$("$HC" fit surf.txt 100 0.05 0.02 0.02 1.0 0.03 0.5 -0.3 500)
b=$("$HC" fit surf.txt 100 0.05 0.02 0.02 1.0 0.03 0.5 -0.3 500)
[ "$a" = "$b" ] && echo "PASS  fit deterministic" || { echo "FAIL  fit not deterministic"; fails=$((fails + 1)); }

# ---- #5 fit argument handling ------------------------------------------------
# maxiter alone (argc 7) is honoured: 1 iteration must NOT converge ...
check_gt "fit ... maxiter=1 honoured (no params)" \
    "$("$HC" fit surf.txt 100 0.05 0.02 1 | mse)" 1e-6
# ... and the default start with a large maxiter does converge.
check_lt "fit with defaults + maxiter=3000" \
    "$("$HC" fit surf.txt 100 0.05 0.02 3000 | mse)" 1e-6
# 5 params + maxiter=1 also honoured
check_gt "fit 5 params + maxiter=1 honoured" \
    "$("$HC" fit surf.txt 100 0.05 0.02 0.02 1.0 0.03 0.5 -0.3 1 | mse)" 1e-6
# partial parameter lists and trailing junk are rejected, not ignored
check_fail "fit rejects 2 start params" "$HC" fit surf.txt 100 0.05 0.02 0.04 2.0
check_fail "fit rejects 4 start params" "$HC" fit surf.txt 100 0.05 0.02 0.04 2.0 0.05 0.3
check_fail "fit rejects extra args"     "$HC" fit surf.txt 100 0.05 0.02 0.04 2.0 0.05 0.3 -0.7 100 junk

echo
if [ "$fails" -eq 0 ]; then echo "All checks passed."; else echo "$fails check(s) FAILED."; fi
[ "$fails" -eq 0 ]
