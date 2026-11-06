#!/bin/sh
# Regression checks for black_scholes and heston_calibrate.
# Usage: ./regress.sh        (builds into a temp dir; exits non-zero on failure)
#
# Monte-Carlo checks use 5-sigma tolerances so they are not flaky despite the
# time-based seed (set BS_SEED for reproducible runs).

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
        'BEGIN { if (v !~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/) exit 1; d = v - e; if (d < 0) d = -d; exit !(d <= t) }'; then
        echo "PASS  $1 ($2)"
    else
        echo "FAIL  $1: got '$2', expected $3 +/- $4"
        fails=$((fails + 1))
    fi
}
# check_lt NAME VALUE LIMIT       (VALUE < LIMIT)
check_lt() {
    if awk -v v="$2" -v l="$3" 'BEGIN { if (v !~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/) exit 1; exit !(v < l) }'; then
        echo "PASS  $1 ($2 < $3)"
    else
        echo "FAIL  $1: got '$2', expected < $3"
        fails=$((fails + 1))
    fi
}
# check_gt NAME VALUE LIMIT       (VALUE > LIMIT)
check_gt() {
    if awk -v v="$2" -v l="$3" 'BEGIN { if (v !~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/) exit 1; exit !(v > l) }'; then
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

# ---- swap-strike grid covers the wings and resolves the K=F kink ------------
# (the old fixed [0.02F, 5F] grid gave 0.893 for sigma=1, T=5)
for case in "1.0 5 1" "2.0 10 4" "0.2 1 0.04" "0.05 0.1 0.0025" "0.3 0.02 0.09" "0.5 0.003 0.25"; do
    set -- $case
    v=$("$BS" var 100 "$2" 0.05 0.02 "$1" 0.04 0.04 2 0.04 0.3 -0.7 | awk '/swap strike/{print $8}')
    check "swap strike sigma=$1 T=$2 equals sigma^2" "$v" "$3" "$(awk -v e="$3" 'BEGIN{print e*2e-4}')"
done

# ---- Heston quadrature adapts to the tenor ----------------------------------
# Fixed 128-node/u_max=200 quadrature gave IV 0.1943 (3% off) at T=0.003.
out=$("$HC" price 100 100 0.003 0.05 0.01 0.04 2.0 0.05 0.5 -0.7)
check "Heston call, T=0.003"   "$(echo "$out" | awk '/^Call:/{print $2}')" 0.442926 2e-6
check "Heston iv, T=0.003"     "$(echo "$out" | awk '/^Call:/{print $7}')" 0.199964 2e-6
out=$("$HC" price 100 100 10 0.05 0.01 0.04 2.0 0.05 0.5 -0.7)
check "Heston call, T=10"      "$(echo "$out" | awk '/^Call:/{print $2}')" 38.864187 2e-6

# ---- independent check: Heston Monte-Carlo vs the Fourier pricer ------------
for case in "100 100 1 0.04 2.0 0.05 0.3 -0.7" "100 120 0.25 0.04 2.0 0.05 0.5 -0.7"; do
    set -- $case
    out=$("$HC" mc "$1" "$2" "$3" 0.05 0.01 "$4" "$5" "$6" "$7" "$8" 300000)
    check "Heston MC call vs Fourier K=$2 T=$3 (|diff| < 5 SE)" \
          "$(echo "$out" | awk '/^Call:/{d=$8; if(d<0)d=-d; print d}')" 0 5
    check "Heston MC put  vs Fourier K=$2 T=$3 (|diff| < 5 SE)" \
          "$(echo "$out" | awk '/^Put:/{d=$8; if(d<0)d=-d; print d}')" 0 5
done

# ---- calibration quality (OTM inversion): round-trip RMSE ---------------------
out=$("$HC" fit surf.txt 100 0.05 0.02 0.02 1.0 0.03 0.5 -0.3 3000)
check_lt "fit RMSE (vol points)" "$(echo "$out" | awk '/^RMSE/{print $2}')" 1e-6
if echo "$out" | grep -q 'n/a\|warning'; then
    echo "FAIL  fit output has warnings / n/a on the committed surface"; fails=$((fails + 1))
else
    echo "PASS  fit output clean (no n/a, no warnings)"
fi
if echo "$out" | grep -q '^Feller ratio'; then echo "PASS  fit reports Feller ratio"
else echo "FAIL  fit does not report Feller ratio"; fails=$((fails + 1)); fi
if "$BS" var 100 1 0.05 0 0.2 0.04 0.04 2 0.04 0.3 -0.7 | grep -q '^Feller ratio'; then
    echo "PASS  var mode reports Feller ratio"
else echo "FAIL  var mode does not report Feller ratio"; fails=$((fails + 1)); fi

# ---- surface file validation -------------------------------------------------
printf '2 2\n0.9 1.1\n1 2\n0.2 0.2 0.2\n'      > "$TMP/trunc.txt"
printf '1 1\n1.0\n1\n-0.2\n'                   > "$TMP/negvol.txt"
printf '1 1\n1.0\n1\n0.2\n0.3\n'               > "$TMP/trail.txt"
printf '1 1\n1.0\n1\nnan\n'                    > "$TMP/nan.txt"
printf '1 1 # hdr\n1.0 # k\n1 # t\n0.2 # vol\n' > "$TMP/comments.txt"
check_fail "surface: truncated file rejected"  "$HC" fit "$TMP/trunc.txt"  100 0.05 0.02
check_fail "surface: negative vol rejected"    "$HC" fit "$TMP/negvol.txt" 100 0.05 0.02
check_fail "surface: trailing data rejected"   "$HC" fit "$TMP/trail.txt"   100 0.05 0.02
check_fail "surface: nan vol rejected"         "$HC" fit "$TMP/nan.txt"     100 0.05 0.02
check_fail "surface: missing file rejected"    "$HC" fit "$TMP/nope.txt"    100 0.05 0.02
if "$HC" fit "$TMP/comments.txt" 100 0.05 0.02 >/dev/null 2>&1; then
    echo "PASS  surface: trailing comments accepted"
else echo "FAIL  surface: trailing comments rejected"; fails=$((fails + 1)); fi

# gen refuses (and writes nothing) when a point cannot be inverted
rm -f "$TMP/bad.txt"
check_fail "gen rejects non-invertible surface" \
    "$HC" gen "$TMP/bad.txt" 100 0.05 0.02 0.04 2 0.05 0.3 -0.7 3 2 0.2 3.0 0.002 0.01
[ ! -e "$TMP/bad.txt" ] && echo "PASS  gen wrote no file on failure" ||
    { echo "FAIL  gen left a file behind"; fails=$((fails + 1)); }
check_fail "gen rejects ns=1 with kmax != kmin" \
    "$HC" gen "$TMP/bad.txt" 100 0.05 0.02 0.04 2 0.05 0.3 -0.7 1 3 0.9 1.1 0.5 1

# ---- strict numeric parsing and usage text -----------------------------------
check_fail "black_scholes rejects non-numeric r"      "$BS" 100 100 1 abc 0 0.2
check_fail "black_scholes rejects trailing junk"      "$BS" 100 100 1x 0.05 0 0.2
check_fail "black_scholes rejects non-numeric PATHS"  "$BS" 100 100 1 0.05 0 0.2 mc 10x
check_fail "var rejects non-numeric v0"               "$BS" var 100 1 0.05 0 0.2 0.04 oops 2 0.04 0.3 -0.7
check_fail "heston_calibrate price rejects non-numeric" "$HC" price 100 100 1 x 0.01 0.04 2 0.05 0.3 -0.7
check_fail "heston_calibrate fit rejects bad maxiter"   "$HC" fit surf.txt 100 0.05 0.02 abc
if "$BS" var 2>&1 | head -1 | grep -q "black_scholes var S T"; then
    echo "PASS  var usage shows the program name"
else echo "FAIL  var usage does not show the program name"; fails=$((fails + 1)); fi

# ---- reproducible Monte-Carlo via BS_SEED -----------------------------------
a=$(BS_SEED=42 "$BS" 100 100 1 0.05 0 0.2 mc 20000 | tail -2)
b=$(BS_SEED=42 "$BS" 100 100 1 0.05 0 0.2 mc 20000 | tail -2)
c=$(BS_SEED=43 "$BS" 100 100 1 0.05 0 0.2 mc 20000 | tail -2)
[ "$a" = "$b" ] && echo "PASS  BS_SEED reproduces black_scholes mc" ||
    { echo "FAIL  BS_SEED did not reproduce black_scholes mc"; fails=$((fails + 1)); }
[ "$a" != "$c" ] && echo "PASS  different BS_SEED gives a different run" ||
    { echo "FAIL  different BS_SEED gave the same run"; fails=$((fails + 1)); }
a=$(BS_SEED=7 "$BS" var 100 1 0.05 0 0.2 0.04 0.04 2 0.04 0.3 -0.7 mc 5000 | tail -3)
b=$(BS_SEED=7 "$BS" var 100 1 0.05 0 0.2 0.04 0.04 2 0.04 0.3 -0.7 mc 5000 | tail -3)
[ "$a" = "$b" ] && echo "PASS  BS_SEED reproduces var mc" ||
    { echo "FAIL  BS_SEED did not reproduce var mc"; fails=$((fails + 1)); }
a=$(BS_SEED=9 "$HC" mc 100 100 1 0.05 0.01 0.04 2 0.05 0.3 -0.7 2000)
b=$(BS_SEED=9 "$HC" mc 100 100 1 0.05 0.01 0.04 2 0.05 0.3 -0.7 2000)
[ "$a" = "$b" ] && echo "PASS  BS_SEED reproduces heston_calibrate mc" ||
    { echo "FAIL  BS_SEED did not reproduce heston_calibrate mc"; fails=$((fails + 1)); }

# ---- Fourier noise floor: noise is not an implied vol ------------------------
# At S=100 the Fourier call price carries absolute noise of about 1e-13, so an
# OTM price below PRICE_NOISE_FLOOR * S = 1e-11 is noise, not a price. At
# T=0.01 these strikes have prices at or below the noise level, yet used to
# invert to garbage (K=50 gave 0.938465, K=62.5 gave 0.647547, K=75 gave
# 0.384044) while K=70 and K=80 printed n/a, so gen accepted or rejected wing
# points at random. All of them must now report n/a.
for K in 50 62.5 70 75 80 120; do
    out=$("$HC" price 100 $K 0.01 0.05 0.01 0.04 2 0.05 0.3 -0.7)
    # 1 when a numeric vol was printed, 0 when the line says n/a
    printed=$(echo "$out" | awk '/^Call:/{print ($7 ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/) ? 1 : 0}')
    check_lt "wing K=$K (price is noise) implied vol is n/a" "$printed" 0.5
done
# The floor must not eat real points: these near-money strikes have OTM prices
# of 6.1e-7, 7.8e-1 and 1.0e-8, far above 1e-11, and keep their vols.
out=$("$HC" price 100 90 0.01 0.05 0.01 0.04 2 0.05 0.3 -0.7)
check "wing K=90 implied vol still reported"  "$(echo "$out" | awk '/^Call:/{print $7}')" 0.22607 1e-4
out=$("$HC" price 100 100 0.01 0.05 0.01 0.04 2 0.05 0.3 -0.7)
check "wing K=100 implied vol still reported" "$(echo "$out" | awk '/^Call:/{print $7}')" 0.20014 1e-4
out=$("$HC" price 100 110 0.01 0.05 0.01 0.04 2 0.05 0.3 -0.7)
check "wing K=110 implied vol still reported" "$(echo "$out" | awk '/^Call:/{print $7}')" 0.17546 1e-4
# A wing-heavy short-tenor grid is therefore unusable: gen refuses it and writes
# nothing. The floor widens what gen rejects on this grid from 11 points to 16;
# the newly caught ones are pure noise inverted to an arbitrary vol (krel=0.625
# at T=0.01, a true price of ~1e-13, used to come out as 0.648).
rm -f "$TMP/noise.txt"
check_fail "gen refuses noise-dominated wing surface" \
    "$HC" gen "$TMP/noise.txt" 100 0.05 0.01 0.04 2 0.05 0.3 -0.7 9 3 0.5 1.5 0.01 0.05
[ ! -e "$TMP/noise.txt" ] && echo "PASS  gen wrote no file for the wing surface" ||
    { echo "FAIL  gen left a file behind for the wing surface"; fails=$((fails + 1)); }
# krel=0.75 (K=75) at T=0.01 is such a point: the pricer returns ~1e-14 of
# noise there, which price mode used to print as vol 0.384044. gen must name it
# among the non-invertible points it refuses to write.
if "$HC" gen "$TMP/noise.txt" 100 0.05 0.01 0.04 2 0.05 0.3 -0.7 9 3 0.5 1.5 0.01 0.05 2>&1 >/dev/null |
        grep -q 'krel=0.75 T=0.01'; then
    echo "PASS  gen names the krel=0.75 T=0.01 noise point as non-invertible"
else echo "FAIL  gen does not report krel=0.75 T=0.01 as non-invertible"; fails=$((fails + 1)); fi

echo
if [ "$fails" -eq 0 ]; then echo "All checks passed."; else echo "$fails check(s) FAILED."; fi
[ "$fails" -eq 0 ]
