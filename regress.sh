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

# CC and REGRESS_CFLAGS may be set by the caller (see the Makefile).
CC=${CC:-gcc}
CFLAGS=${REGRESS_CFLAGS:--O2}
$CC $CFLAGS -Wall -Wextra -o "$BS" black_scholes.c -lm || exit 1
$CC $CFLAGS -Wall -Wextra -o "$HC" heston_calibrate.c -lm || exit 1

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
# Numeric, not byte-for-byte: the 1-month k=1.3 price is 2.5e-9, where the last
# digits of the vol are rounding noise that differs between math libraries.
surf_close() {
    awk -v tol=1e-6 '!/^#/ { for (i = 1; i <= NF; i++)
            if (NR == FNR) a[++n] = $i
            else { m++; d = $i - a[m]; if (d < 0) d = -d; if (d > tol) bad = 1 } }
        END { exit (bad || n != m || n == 0) }' "$1" "$2"
}
if surf_close "$TMP/surf.txt" surf.txt; then
    echo "PASS  gen reproduces committed surf.txt (within 1e-6)"
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
    $CC $CFLAGS -include "$TMP/poison$fill.h" -o "$TMP/hc_$fill" heston_calibrate.c -lm || exit 1
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

# ---- fit reports non-convergence --------------------------------------------
# Running out the maxiter budget must be visible in the output (a truncated
# Nelder-Mead run otherwise looks like a converged fit), and a fit that does
# meet the convergence test must stay quiet. Exit status is 0 either way.
out=$("$HC" fit surf.txt 100 0.05 0.02 0.02 1.0 0.03 0.5 -0.3 20); rc=$?
check "fit at maxiter=20 still exits 0" "$rc" 0 0
check "fit at maxiter=20 prints one such warning" \
    "$(echo "$out" | grep -c 'stopped at maxiter')" 1 0
if echo "$out" | grep -q 'warning: stopped at maxiter=20 without converging; increase maxiter'; then
    echo "PASS  fit warns when it stops at maxiter=20"
else
    echo "FAIL  fit is silent when it stops at maxiter=20"; fails=$((fails + 1))
fi
# the same start does meet the convergence test by maxiter=3000
out=$("$HC" fit surf.txt 100 0.05 0.02 0.02 1.0 0.03 0.5 -0.3 3000); rc=$?
check "converged fit exits 0" "$rc" 0 0
if echo "$out" | grep -q 'stopped at maxiter'; then
    echo "FAIL  converged fit (maxiter=3000) warns anyway"; fails=$((fails + 1))
else
    echo "PASS  converged fit (maxiter=3000) does not warn"
fi

# ---- swap-strike grid: mean shift, no spacing floor, size caps ---------------
# r = q = 0 (F = S): the strike is sigma^2 however large sigma sqrt(T) is. The
# old grid was symmetric about ln F, so the put wing was cut at ln F - 8 s
# instead of ln F - s^2/2 - 8 s (s = sigma sqrt(T)); sigma=2, T=100 gave 3.216.
for case in "1.0 100 1" "2.0 25 4" "2.0 100 4" "3.0 100 9" "5.0 25 25"; do
    set -- $case
    v=$("$BS" var 100 "$2" 0 0 "$1" 0.04 0.04 2 0.04 0.3 -0.7 | awk '/swap strike/{print $8}')
    check "swap strike sigma=$1 T=$2 (wing follows -s^2/2)" "$v" "$3" \
          "$(awk -v e="$3" 'BEGIN{print e*1e-4}')"
done
# Tiny sigma sqrt(T): the old 1e-4 floor on the spacing left the kink at K = F
# unresolved, so sigma=0.001, T=0.001 gave 0.000003 instead of 0.000001. The
# strike of sigma=0.0001, T=0.01 is 1e-8, which prints as 0.000000 at 6
# decimals; sigma=0.0316227766, T=1e-07 lands on the same grid (same
# sigma sqrt(T) = 1e-5) at a scale where the printed digits bound the relative
# error to 5e-4.
for case in "0.05 0.003 0.0025" "0.01 0.001 0.0001" "0.001 0.001 0.000001" "0.0001 0.01 0.00000001"; do
    set -- $case
    v=$("$BS" var 100 "$2" 0.05 0.02 "$1" 0.04 0.04 2 0.04 0.3 -0.7 | awk '/swap strike/{print $8}')
    check "swap strike sigma=$1 T=$2 (no spacing floor)" "$v" "$3" 1e-6
done
check "swap strike sigma=0.0316227766 T=1e-07 (same grid, scaled up)" \
    "$("$BS" var 100 1e-07 0.05 0.02 0.0316227766 0.04 0.04 2 0.04 0.3 -0.7 | awk '/swap strike/{print $8}')" \
    0.001 5e-7
# A grid too large to integrate is refused with the error message rather than a
# wrong number (the point count used to go through an int cast that overflowed
# for sigma=1e5, T=1e4).
check_fail "swap strike rejects sigma=10, T=100 (strike limit beyond exp range)" \
    "$BS" var 100 100 0 0 10 0.04 0.04 2 0.04 0.3 -0.7
check_fail "swap strike rejects sigma=100000, T=10000" \
    "$BS" var 100 10000 0 0 100000 0.04 0.04 2 0.04 0.3 -0.7
if "$BS" var 100 100 0 0 10 0.04 0.04 2 0.04 0.3 -0.7 2>&1 >/dev/null |
        grep -q "Failed to compute variance swap strike"; then
    echo "PASS  oversized swap grid reports the Failed message"
else
    echo "FAIL  oversized swap grid does not report the Failed message"
    fails=$((fails + 1))
fi

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
# krel=0.75 (K=75) at T=0.01 is such a point: the pricer returns ~1e-13 of
# noise there, which price mode used to print as vol 0.384044. gen must name it
# among the non-invertible points it refuses to write.
if "$HC" gen "$TMP/noise.txt" 100 0.05 0.01 0.04 2 0.05 0.3 -0.7 9 3 0.5 1.5 0.01 0.05 2>&1 >/dev/null |
        grep -q 'krel=0.75 T=0.01'; then
    echo "PASS  gen names the krel=0.75 T=0.01 noise point as non-invertible"
else echo "FAIL  gen does not report krel=0.75 T=0.01 as non-invertible"; fails=$((fails + 1)); fi

# ---- Heston pricer, Lewis form (sharp feature at u = 0) ---------------------
# The old P1/P2 form with uniform width-10 panels gave call 40.968498 and put
# -4.150339 on the first case; even width-0.01 panels gave put -0.29 there.
# Four of these five cases have kappa < rho*xi; the rho = -0.7 case was only
# mildly off before (51.899300 against 51.8879). All of them are now priced by
# the Lewis single-integral form with panels graded from PANEL_WIDTH0.
for case in "100 100 30 0.02 0 0.04 0.1 0.1 1.2 0.7|52.5295|7.4106" \
            "100 100 10 0.02 0 0.04 0.1 0.1 1.2 0.7|24.2685|6.1416" \
            "100 100 30 0.02 0 0.04 0.1 0.1 1.2 -0.7|51.8879|6.7691" \
            "100 100 10 0.02 0 0.04 0.5 0.06 1.5 0.7|27.4971|9.3701" \
            "100 100 5 0.02 0 0.25 0.5 0.25 2 0.8|38.6427|29.1265"; do
    args=${case%%|*}
    rest=${case#*|}
    ec=${rest%%|*}
    ep=${rest#*|}
    # shellcheck disable=SC2086
    out=$("$HC" price $args)
    c=$(echo "$out" | awk '/^Call:/{print $2}')
    p=$(echo "$out" | awk '/^Call:/{print $4}')
    check "Heston call $args" "$c" "$ec" 1e-3
    check "Heston put  $args" "$p" "$ep" 1e-3
    # the old quadrature produced a negative put on the first case
    check_gt "Heston put non-negative $args" "$p" -1e-9
done

# ---- small xi: Heston tends to Black-Scholes --------------------------------
# With v0 = theta, xi -> 0 is Black-Scholes at vol sqrt(v0). a-d and the log both
# used to cancel before being divided by xi^2: xi=1e-6 gave 9.827775, xi=1e-8
# gave 8.626045 and xi=1e-10 gave 3.890051.
for xi in 1e-6 1e-7 1e-8 1e-10 1e-100; do
    check "Heston call xi=$xi equals BS" \
        "$("$HC" price 100 100 1 0.05 0.01 0.04 2 0.04 "$xi" -0.7 | awk '/^Call:/{print $2}')" 9.826298 2e-6
done
check "Heston call xi=1e-8 K=120 T=2 rho=-0.9 equals BS" \
    "$("$HC" price 100 120 2 0.05 0.01 0.04 2 0.04 1e-8 -0.9 | awk '/^Call:/{print $2}')" 7.083326 2e-6
check "Heston call xi=1e-8 K=100 T=30 rho=0.9 equals BS" \
    "$("$HC" price 100 100 30 0.05 0.01 0.04 2 0.04 1e-8 0.9 | awk '/^Call:/{print $2}')" 54.565904 2e-6
check_fail "price rejects xi below 1e-100" "$HC" price 100 100 1 0.05 0.01 0.04 2 0.04 1e-160 -0.7

# ---- Fourier cutoff and Monte-Carlo step limits ------------------------------
# A cutoff of 40000 silently truncated the integral at T=1e-7 (iv 0.108).
for T in 1e-7 1e-8; do
    check "Heston iv, T=$T (vol 0.1)" \
        "$("$HC" price 100 100 "$T" 0 0 0.01 2 0.01 0.3 -0.7 | awk '/^Call:/{print $7}')" 0.1 1e-5
done
check_fail "price refuses a Fourier integral that cannot converge" \
    "$HC" price 100 100 1e-10 0 0 0.01 2 0.01 0.3 -0.7
# The step count used to overflow its int (undefined behaviour, then garbage).
check_fail "var mc rejects too many time steps" \
    "$BS" var 100 1e10 0 0 0 0.04 0.04 2 0.04 0.3 -0.7 mc 1
check_fail "heston mc rejects too many time steps" \
    "$HC" mc 100 100 1e8 0 0 0.04 2 0.05 0.3 -0.7 2

# ---- mc warns when put-call parity shows the sample missed the tail ----------
# Positive rho with high xi makes S_T so heavy-tailed that the call is far too
# low while its SE looks fine. E[S_T] = F holds for any step size, so a parity
# gap can only be sampling error. A single run need not warn, because one huge
# path also inflates the SE (this case fired in 29 of 30 seeds, and a fixed seed
# does not give the same paths on every math library), so any of three will do.
warned=0
for seed in 1 2 3; do
    if BS_SEED=$seed "$HC" mc 100 100 30 0.02 0 0.04 0.1 0.1 1.2 0.7 2000 |
            grep -q 'warning: call - put misses parity'; then warned=1; break; fi
done
if [ "$warned" -eq 1 ]; then echo "PASS  mc warns on a heavy-tailed sample"
else echo "FAIL  mc is silent on a heavy-tailed sample"; fails=$((fails + 1)); fi
if BS_SEED=1 "$HC" mc 100 100 1 0.05 0.01 0.04 2 0.05 0.3 -0.7 20000 | grep -q 'warning'; then
    echo "FAIL  mc warns on ordinary parameters"; fails=$((fails + 1))
else echo "PASS  mc does not warn on ordinary parameters"; fi

echo
if [ "$fails" -eq 0 ]; then echo "All checks passed."; else echo "$fails check(s) FAILED."; fi
[ "$fails" -eq 0 ]
