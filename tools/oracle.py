#!/usr/bin/env python3
"""Exhaustive bit-exactness check: SoftUnum (C) vs SoftPosit pX2 (es=2).

SoftPosit's generic pX2 path is es=2 at any width -- the 2022-standard
convention SoftUnum targets, and the same oracle /lib/unum was validated
against.  posit8 is small enough to check EVERY value / pair.

    /opt/anaconda3/bin/python tools/oracle.py
"""
import ctypes, os, random, sys

try:
    import softposit as sp
except ImportError:
    sys.exit("need the `softposit` pip package (the pX2 oracle)")

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "libsoftunum.dylib")
if not os.path.exists(LIB):
    LIB = os.path.join(HERE, "..", "libsoftunum.so")
lib = ctypes.CDLL(LIB)

u32, u64, i64 = ctypes.c_uint32, ctypes.c_uint64, ctypes.c_int64
for nm in ("p8_neg", "p8_abs", "p8_sgn", "p8_nearest_int", "p8_floor", "p8_ceil"):
    getattr(lib, nm).argtypes = [u32]; getattr(lib, nm).restype = u32
for nm in ("p8_add", "p8_sub", "p8_mul", "p8_div"):
    getattr(lib, nm).argtypes = [u32, u32]; getattr(lib, nm).restype = u32
for nm in ("p8_eq", "p8_lt", "p8_le", "p8_gt", "p8_ge"):
    getattr(lib, nm).argtypes = [u32, u32]; getattr(lib, nm).restype = ctypes.c_int
lib.p8_sqrt.argtypes = [u32]; lib.p8_sqrt.restype = u32
lib.p8_fma.argtypes = [u32, u32, u32]; lib.p8_fma.restype = u32
lib.p8_from_u64.argtypes = [u64]; lib.p8_from_u64.restype = u32
lib.p8_from_i64.argtypes = [i64]; lib.p8_from_i64.restype = u32
lib.p8_to_i64.argtypes = [u32, ctypes.POINTER(i64)]; lib.p8_to_i64.restype = ctypes.c_int
lib.p8_fdp.argtypes = [ctypes.POINTER(u32), ctypes.POINTER(u32), i64]
lib.p8_fdp.restype = u32

N = 8
def mk(pat):  o = sp.convertDoubleToPX2(0.0, N); o.v = (pat << (32 - N)) & 0xffffffff; return o
def pt(p):    return p.v >> (32 - N)

fails = []
def check(name, got, exp, *args):
    if got != exp:
        fails.append((name, args, got, exp))

#  ---- arithmetic: all 65,536 pairs ----------------------------------------
print("posit8: add/sub/mul/div over all 65,536 pairs ...", flush=True)
for a in range(256):
    for b in range(256):
        check("add", lib.p8_add(a, b), pt(sp.pX2_add(mk(a), mk(b), N)), a, b)
        check("sub", lib.p8_sub(a, b), pt(sp.pX2_sub(mk(a), mk(b), N)), a, b)
        check("mul", lib.p8_mul(a, b), pt(sp.pX2_mul(mk(a), mk(b), N)), a, b)
        check("div", lib.p8_div(a, b), pt(sp.pX2_div(mk(a), mk(b), N)), a, b)
        check("lt", lib.p8_lt(a, b), 1 if sp.pX2_lt(mk(a), mk(b)) else 0, a, b)
        check("le", lib.p8_le(a, b), 1 if sp.pX2_le(mk(a), mk(b)) else 0, a, b)
        check("eq", lib.p8_eq(a, b), 1 if sp.pX2_eq(mk(a), mk(b)) else 0, a, b)

#  ---- unary over all 256 ---------------------------------------------------
print("posit8: sqrt / round / neg over all 256 ...", flush=True)
for a in range(256):
    check("sqrt", lib.p8_sqrt(a), pt(sp.pX2_sqrt(mk(a), N)), a)
    check("neg", lib.p8_neg(a), (256 - a) & 0xff, a)
    #  SoftPosit's roundToInt mishandles 0x81 (a known quirk); skip it.
    if a != 0x81:
        check("rint", lib.p8_nearest_int(a), pt(sp.pX2_roundToInt(mk(a), N)), a)

#  ---- integer conversion ---------------------------------------------------
print("posit8: i64->posit over [-300,300] ...", flush=True)
for v in range(-300, 301):
    check("from_i64", lib.p8_from_i64(v), pt(sp.i32_to_pX2(v, N)), v)

#  ---- fma: random sample ---------------------------------------------------
print("posit8: fma (200k random triples) ...", flush=True)
random.seed(135)
for _ in range(200_000):
    a, b, c = random.randrange(256), random.randrange(256), random.randrange(256)
    check("fma", lib.p8_fma(a, b, c), pt(sp.pX2_mulAdd(mk(a), mk(b), mk(c), N)), a, b, c)

#  ---- fdp: random vectors vs the quire -------------------------------------
print("posit8: fdp (5k random vectors) ...", flush=True)
for _ in range(5_000):
    L = random.randint(1, 12)
    av = [random.randrange(256) for _ in range(L)]
    bv = [random.randrange(256) for _ in range(L)]
    q = sp.qX2Clr()
    for x, y in zip(av, bv):
        q = sp.qX2_fdp_add(q, mk(x), mk(y))
    exp = pt(sp.qX2_to_pX2(q, N))
    ca = (u32 * L)(*av); cb = (u32 * L)(*bv)
    check("fdp", lib.p8_fdp(ca, cb, L), exp, tuple(av), tuple(bv))

#  ---- report ---------------------------------------------------------------
if not fails:
    print("\nALL posit8 checks PASS (bit-exact vs SoftPosit pX2).")
    sys.exit(0)
from collections import Counter
c = Counter(f[0] for f in fails)
print(f"\n{len(fails)} FAILURES: {dict(c)}")
for f in fails[:20]:
    print("  %-9s args=%s got=0x%x exp=0x%x" % (f[0], f[1], f[2], f[3]))
sys.exit(1)
