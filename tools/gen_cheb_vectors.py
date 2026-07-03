#!/usr/bin/env python3
"""Generate a large (input, expected-output) vector corpus for the new
range-reduced transcendentals, using libmath/tools/unum_cheb_check.py's
already mpmath-verified reference functions (exp_ref/log_g/log2_g/log10_g/
sin_g/cos_g/tan_g/atan_g/asin_g/acos_g) -- the SAME algorithm-of-record the
Hoon /lib/unum rewrite was checked against.  This is deliberately independent
of both the Hoon test vectors AND oracle.py's own transliteration: it re-runs
the ORIGINAL Python model that produced the Hoon's baked-in constants.

Output: a plain-text file, one line per case:  WIDTH FUNC IN_HEX OUT_HEX
Consumed by tools/cheb_check.c (linked directly against libsoftunum, no
ctypes -- exercises the real C ABI).

    /opt/anaconda3/bin/python tools/gen_cheb_vectors.py [outfile] [n_random]
"""
import os
import random
import sys
from fractions import Fraction

NUMERICS_TOOLS = os.path.expanduser("~/urbit/numerics/libmath/tools")
sys.path.insert(0, NUMERICS_TOOLS)

import mpmath as mp
mp.mp.dps = 60

import unum_cheb_check as ucc
from posit_check import ref_value_encode

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "cheb_vectors.txt")
N_RANDOM = int(sys.argv[2]) if len(sys.argv) > 2 else 5000

# ---- populate unum_cheb_check's module-global coefficient tables (mirrors
# what check_exp()/check_log()/check_trig()/check_atan() do, minus their
# verbose sweep printouts) -----------------------------------------------
#
# NB: use exp_g (the g-triple re-simulation), NOT exp_ref/exp_exact -- the
# latter materializes `Fraction(2) ** k` exactly, which is fine for
# check_exp()'s own deliberately bounded x in [-40,40] grid but is
# catastrophic (an integer with ~maxpos/ln2 bits, i.e. >10^16 bits at
# posit16/32 maxpos) for a full-range sweep that includes large-magnitude
# posits.  exp_g never materializes 2^k -- it tracks k as a plain int folded
# only into the final `encode()` exponent, exactly like the Hoon +exp and
# this project's own g_bit/ll_add_clamp -- and is verified below to agree
# with exp_ref everywhere the latter is actually tractable.
EXP_COEFFS = ucc.gen_exp_coeffs(7)
EXP_COEFFS_G = [(True, -128, ucc.raw128(c)) for c in reversed(EXP_COEFFS)]
ucc.LOG_COEFFS = ucc.gen_log_coeffs(16)
ucc.SIN_COEFFS, ucc.COS_COEFFS = ucc.gen_sincos_coeffs(6)
ucc.PI2_G = (True, -ucc.TRIG_WBITS, ucc.raw_w(ucc.PI2, ucc.TRIG_WBITS))
ucc.INVPI2_G = (True, -ucc.TRIG_WBITS, ucc.raw_w(ucc.INVPI2, ucc.TRIG_WBITS))
ucc.ATAN_COEFFS = ucc.gen_atan_coeffs(13)
ucc.BP_ANGLES = [(True, 0, 0),
                 (True, -128, ucc.raw128(ucc.quant(mp.atan(mp.mpf(1) / 2)))),
                 (True, -128, ucc.raw128(ucc.quant(mp.pi / 4))),
                 (True, -128, ucc.raw128(ucc.quant(mp.atan(mp.mpf(3) / 2)))),
                 (True, -128, ucc.raw128(ucc.quant(mp.pi / 2)))]

WIDTHS = [8, 16, 32]

def hx(v, n): return f"0x{v & ((1 << n) - 1):x}"

FUNCS = {
    "exp":    lambda p, n: ucc.exp_g(p, n, EXP_COEFFS_G),
    "log":    ucc.log_g,
    "log2":   ucc.log2_g,
    "log10":  ucc.log10_g,
    "sin":    ucc.sin_g,
    "cos":    ucc.cos_g,
    "tan":    ucc.tan_g,
    "atan":   ucc.atan_g,
    "asin":   ucc.asin_g,
    "acos":   ucc.acos_g,
}

def edge_patterns(n):
    nar = 1 << (n - 1)
    pats = [0, nar, ref_value_encode(Fraction(1), n), ref_value_encode(Fraction(-1), n),
            ref_value_encode(Fraction(1, 2), n), ref_value_encode(Fraction(-1, 2), n),
            ref_value_encode(Fraction(2), n), ref_value_encode(Fraction(-2), n),
            ref_value_encode(Fraction(10), n), ref_value_encode(Fraction(-10), n),
            ref_value_encode(Fraction(100), n), ucc.true_pattern(mp.pi, n),
            ucc.true_pattern(mp.pi / 2, n), ucc.true_pattern(mp.pi / 4, n),
            (1 << (n - 1)) - 1,   # maxpos
            1,                    # minpos
            nar + 1, nar - 1]
    return pats

def _sanity_check_exp_g():
    """exp_g (fast, used for the full sweep) must agree with exp_ref (slow,
    exact-Fraction, only tractable for bounded |x|) wherever both are cheap
    to evaluate -- guards against the two reductions silently diverging."""
    for xv in [Fraction(t, 7) for t in range(-280, 281, 3)]:
        for n in WIDTHS:
            p = ref_value_encode(xv, n)
            a = ucc.exp_ref(p, n, EXP_COEFFS)
            b = ucc.exp_g(p, n, EXP_COEFFS_G)
            assert a == b, f"exp_g/exp_ref disagree at n={n} x={xv}: {hex(a)} vs {hex(b)}"

def main():
    _sanity_check_exp_g()
    random.seed(20260702)
    lines = []
    for n in WIDTHS:
        pats = list(range(1 << n)) if n == 8 else \
               list(dict.fromkeys(edge_patterns(n) + [random.randrange(1 << n) for _ in range(N_RANDOM)]))
        for name, fn in FUNCS.items():
            for p in pats:
                try:
                    out = fn(p, n)
                except Exception as e:
                    print(f"!! {name} n={n} p={hx(p,n)} raised {e}", file=sys.stderr)
                    continue
                lines.append(f"{n} {name} {hx(p, n)} {hx(out, n)}")
    with open(OUT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {len(lines)} vectors to {OUT}")

if __name__ == "__main__":
    main()
