#!/usr/bin/env python3
"""Verify SoftUnum's posit<->IEEE-754 conversions against numpy (the IEEE
reference) and the exact posit value.

  to_r*  : posit value (exact) -> correctly-rounded binary{16,32,64}, vs numpy.
           posit8 exhaustive; posit16/32 sampled.  (binary128 via round-trip.)
  from_r*: binary16 EXHAUSTIVE (all 65,536, incl subnormals/inf/nan) -> posit,
           vs an independent reference (numpy value -> from-scratch posit encode);
           binary32/64 sampled.
  rq     : posit -> binary128 -> posit is the identity (binary128 holds every
           posit8/16/32 value exactly), checked over posit8 + samples.

    /opt/anaconda3/bin/python tools/ieee_check.py
"""
import ctypes, os, struct, sys, random
from fractions import Fraction
import numpy as np
np.seterr(over="ignore")   #  huge posit values cast to inf in float16 -- expected

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "libsoftunum.dylib")
if not os.path.exists(LIB): LIB = os.path.join(HERE, "..", "libsoftunum.so")
lib = ctypes.CDLL(LIB)
u32, u64 = ctypes.c_uint32, ctypes.c_uint64

#  --- independent posit decode (value) and encode (mirrors /lib/unum) ---------
def pvalue(p, n):                       #  exact value of posit bits -> Fraction | 'z' | 'n'
    p &= (1 << n) - 1
    if p == 0: return Fraction(0)
    if p == 1 << (n - 1): return "n"
    neg = (p >> (n - 1)) & 1
    mag = ((1 << n) - p) if neg else p
    pw = n - 1; r0 = (mag >> (pw - 1)) & 1; k = 1
    while k < pw and (mag >> (pw - 1 - k)) & 1 == r0: k += 1
    r = (k - 1) if r0 == 1 else -k
    remwid = pw - (k + 1); rem = mag & ((1 << max(remwid, 0)) - 1)
    if remwid >= 2: elo = rem >> (remwid - 2); fw = remwid - 2
    elif remwid == 1: elo = rem << 1; fw = 0
    else: elo = 0; fw = 0
    frac = rem & ((1 << fw) - 1)
    val = Fraction(2) ** (4 * r + elo - fw) * ((1 << fw) + frac)
    return -val if neg else val

def pencode(neg, e, a, n):              #  /lib/unum +bit (round-nearest-even)
    if a == 0: return 0
    M = (1 << n) - 1; maxpos = (1 << (n - 1)) - 1
    lead = a.bit_length() - 1; x = e + lead; frac = a & ((1 << lead) - 1)
    r = x >> 2; elo = x - 4 * r
    if r >= n - 2: return ((1 << n) - maxpos) & M if neg else maxpos
    if r <= -(n - 1): return ((1 << n) - 1) & M if neg else 1
    if r >= 0: regval = ((1 << (r + 1)) - 1) << 1; regwid = r + 2
    else: regval = 1; regwid = -r + 1
    totw = regwid + 2 + lead; pay = (regval << (2 + lead)) | (elo << lead) | frac; pw = n - 1
    if totw <= pw:
        mag = pay << (pw - totw)
    else:
        sh = totw - pw; keep = pay >> sh
        guard = (pay >> (sh - 1)) & 1; low = pay & ((1 << (sh - 1)) - 1)
        if guard and (low or (keep & 1)): keep += 1
        if keep > maxpos: keep = maxpos
        mag = keep
    return ((1 << n) - mag) & M if neg else mag

def frac_to_posit(v, n):               #  exact dyadic Fraction -> posit bits
    if v == 0: return 0
    neg = v < 0; m = -v if neg else v
    num, den = m.numerator, m.denominator   #  den is a power of two
    k = den.bit_length() - 1
    return pencode(neg, -k, num, n)

#  --- IEEE helpers via numpy / struct ----------------------------------------
def f_to_bits(d, kind):
    if kind == 16: return int(np.float16(d).view(np.uint16))
    if kind == 32: return struct.unpack("<I", struct.pack("<f", np.float32(d)))[0]
    return struct.unpack("<Q", struct.pack("<d", float(d)))[0]
def bits_to_np(bits, kind):
    if kind == 16: return np.uint16(bits).view(np.float16)
    if kind == 32: return struct.unpack("<f", struct.pack("<I", bits))[0]
    return struct.unpack("<d", struct.pack("<Q", bits))[0]
def is_nan_bits(bits, w, p):
    e = (bits >> p) & ((1 << w) - 1); f = bits & ((1 << p) - 1)
    return e == (1 << w) - 1 and f != 0

FMT = {16: (5, 10), 32: (8, 23), 64: (11, 52)}
fails = []
def fail(t): fails.append(t)

#  ---- to_r* : posit -> float -----------------------------------------------
def check_to(n):
    pfx = f"p{n}"
    sig = {16: (f"{pfx}_to_rh", u32, 16), 32: (f"{pfx}_to_rs", u32, 32), 64: (f"{pfx}_to_rd", u64, 64)}
    for nm, rt, _k in sig.values():
        getattr(lib, nm).argtypes = [u32]; getattr(lib, nm).restype = rt
    inputs = range(1 << n) if n == 8 else (
        [0, 1 << (n - 1), 1 << (n - 2), (1 << (n - 1)) - 1, 1] +
        [random.randrange(1 << n) for _ in range(20000)])
    for p in inputs:
        v = pvalue(p, n)
        for kind in (16, 32, 64):
            w, pp = FMT[kind]
            got = getattr(lib, sig[kind][0])(p) & ((1 << (w + pp + 1)) - 1)
            if v == "n":                       #  NaR -> some NaN
                if not is_nan_bits(got, w, pp): fail(("to", n, kind, p, hex(got), "NaN"))
                continue
            exp = f_to_bits(float(v), kind)
            if got != exp: fail(("to", n, kind, p, hex(got), hex(exp)))

#  ---- from_r* : float -> posit ---------------------------------------------
def check_from(n):
    pfx = f"p{n}"
    for nm, at in ((f"{pfx}_from_rh", u32), (f"{pfx}_from_rs", u32), (f"{pfx}_from_rd", u64)):
        getattr(lib, nm).argtypes = [at]; getattr(lib, nm).restype = u32
    def oracle(bits, kind):
        w, p = FMT[kind]
        e = (bits >> p) & ((1 << w) - 1); f = bits & ((1 << p) - 1)
        if e == (1 << w) - 1: return 1 << (n - 1)        #  inf / nan -> NaR
        d = bits_to_np(bits, kind)
        return frac_to_posit(Fraction(float(d)), n)
    #  binary16 exhaustive (subnormals, inf, nan all covered)
    for bits in range(1 << 16):
        got = lib.__getattr__(f"{pfx}_from_rh")(bits) & ((1 << n) - 1)
        exp = oracle(bits, 16)
        if got != exp: fail(("from", n, 16, hex(bits), hex(got), hex(exp)))
    #  binary32 / binary64 sampled
    rng = random.Random(5)
    for kind, fn in ((32, f"{pfx}_from_rs"), (64, f"{pfx}_from_rd")):
        w, p = FMT[kind]
        for _ in range(30000):
            bits = rng.randrange(1 << (w + p + 1))
            got = getattr(lib, fn)(bits) & ((1 << n) - 1)
            exp = oracle(bits, kind)
            if got != exp: fail(("from", n, kind, hex(bits), hex(got), hex(exp)))

#  ---- rq round-trip identity (binary128 holds every posit8/16/32 value) -----
def check_rq(n):
    pfx = f"p{n}"
    toq = getattr(lib, f"{pfx}_to_rq"); toq.argtypes = [u32, ctypes.POINTER(u64 * 2)]
    frq = getattr(lib, f"{pfx}_from_rq"); frq.argtypes = [ctypes.POINTER(u64 * 2)]; frq.restype = u32
    nar = 1 << (n - 1)
    inputs = range(1 << n) if n == 8 else [0, nar, 1 << (n - 2), (1 << (n - 1)) - 1, 1] + \
        [random.randrange(1 << n) for _ in range(20000)]
    for p in inputs:
        buf = (u64 * 2)(); toq(p, ctypes.byref(buf))
        rt = frq(ctypes.byref(buf)) & ((1 << n) - 1)
        exp = nar if p == nar else p           #  identity (NaR -> NaR)
        if rt != exp: fail(("rq-rt", n, p, hex(rt), hex(exp)))

random.seed(135)
for n in (8, 16, 32):
    print(f"posit{n}: to_r* (numpy) / from_r* (binary16 exhaustive) / rq round-trip ...", flush=True)
    check_to(n); check_from(n); check_rq(n)

if not fails:
    print("\nALL IEEE conversion checks PASS (vs numpy + exact value).")
    sys.exit(0)
from collections import Counter
print(f"\n{len(fails)} FAILURES: {dict(Counter(f[0] for f in fails))}")
for f in fails[:20]: print("  ", f)
sys.exit(1)
