#!/usr/bin/env python3
"""Bit-exactness check: SoftUnum (C) vs SoftPosit pX2 (es=2).

SoftPosit's generic pX2 path is es=2 at any width -- the 2022-standard
convention SoftUnum targets, and the same oracle /lib/unum was validated
against.  posit8 is checked EXHAUSTIVELY (every value / every pair); posit16
and posit32 are checked on a large random sample plus a structured edge grid
(0, NaR, +-1, +-maxpos, +-minpos, small integers).

    /opt/anaconda3/bin/python tools/oracle.py
"""
import ctypes, os, random, sys
from collections import Counter

try:
    import softposit as sp
except ImportError:
    sys.exit("need the `softposit` pip package (the pX2 oracle)")

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "libsoftunum.dylib")
if not os.path.exists(LIB):
    LIB = os.path.join(HERE, "..", "libsoftunum.so")
lib = ctypes.CDLL(LIB)

u32, u64, i64, ci = ctypes.c_uint32, ctypes.c_uint64, ctypes.c_int64, ctypes.c_int

def bind(pfx):
    fns = {}
    for nm in ("neg", "abs", "sgn", "nearest_int", "floor", "ceil"):
        f = getattr(lib, f"{pfx}_{nm}"); f.argtypes = [u32]; f.restype = u32; fns[nm] = f
    for nm in ("add", "sub", "mul", "div"):
        f = getattr(lib, f"{pfx}_{nm}"); f.argtypes = [u32, u32]; f.restype = u32; fns[nm] = f
    for nm in ("eq", "lt", "le", "gt", "ge"):
        f = getattr(lib, f"{pfx}_{nm}"); f.argtypes = [u32, u32]; f.restype = ci; fns[nm] = f
    f = getattr(lib, f"{pfx}_sqrt"); f.argtypes = [u32]; f.restype = u32; fns["sqrt"] = f
    f = getattr(lib, f"{pfx}_fma"); f.argtypes = [u32, u32, u32]; f.restype = u32; fns["fma"] = f
    f = getattr(lib, f"{pfx}_from_i64"); f.argtypes = [i64]; f.restype = u32; fns["from_i64"] = f
    f = getattr(lib, f"{pfx}_fdp")
    f.argtypes = [ctypes.POINTER(u32), ctypes.POINTER(u32), i64]; f.restype = u32; fns["fdp"] = f
    return fns

def oracle(n):
    #  posit32 (es=2) coincides with the 2022 standard, so use SoftPosit's
    #  DEDICATED p32 path -- its generic pX2 misrounds tiny values at width 32
    #  (pX2 was only validated to X~=20).  posit8/16 use pX2 (es=2 at any width).
    if n == 32:
        def mk(p):  o = sp.convertDoubleToP32(0.0); o.v = p & 0xffffffff; return o
        def pt(p):  return p.v & 0xffffffff
        return dict(
            mk=mk, pt=pt,
            add=lambda a, b: sp.p32_add(a, b), sub=lambda a, b: sp.p32_sub(a, b),
            mul=lambda a, b: sp.p32_mul(a, b), div=lambda a, b: sp.p32_div(a, b),
            lt=lambda a, b: sp.p32_lt(a, b), le=lambda a, b: sp.p32_le(a, b),
            eq=lambda a, b: sp.p32_eq(a, b), sqrt=lambda a: sp.p32_sqrt(a),
            rint=lambda a: sp.p32_roundToInt(a),
            fma=lambda a, b, c: sp.p32_mulAdd(a, b, c), from_i=lambda v: sp.i32_to_p32(v),
            qclr=sp.q32Clr, qadd=sp.q32_fdp_add, qtop=lambda q: sp.q32_to_p32(q))
    def mk(p):  o = sp.convertDoubleToPX2(0.0, n); o.v = (p << (32 - n)) & 0xffffffff; return o
    def pt(p):  return p.v >> (32 - n)
    return dict(
        mk=mk, pt=pt,
        add=lambda a, b: sp.pX2_add(a, b, n), sub=lambda a, b: sp.pX2_sub(a, b, n),
        mul=lambda a, b: sp.pX2_mul(a, b, n), div=lambda a, b: sp.pX2_div(a, b, n),
        lt=lambda a, b: sp.pX2_lt(a, b), le=lambda a, b: sp.pX2_le(a, b),
        eq=lambda a, b: sp.pX2_eq(a, b), sqrt=lambda a: sp.pX2_sqrt(a, n),
        rint=lambda a: sp.pX2_roundToInt(a, n),
        fma=lambda a, b, c: sp.pX2_mulAdd(a, b, c, n), from_i=lambda v: sp.i32_to_pX2(v, n),
        qclr=sp.qX2Clr, qadd=sp.qX2_fdp_add, qtop=lambda q: sp.qX2_to_pX2(q, n))

def run(n, exhaustive):
    pfx = {8: "p8", 16: "p16", 32: "p32"}[n]
    fn = bind(pfx)
    M = (1 << n) - 1
    O = oracle(n)
    mk, pt = O["mk"], O["pt"]
    fails = []
    def chk(name, got, exp, *a):
        if got != exp: fails.append((name, a, got, exp))

    nar, one, maxp, minp = 1 << (n - 1), 1 << (n - 2), (1 << (n - 1)) - 1, 1
    edges = [0, nar, one, M ^ (one - 0) & M, maxp, minp,
             (M + 1 - one) & M, (M + 1 - maxp) & M, (M + 1 - minp) & M,
             2, 3, M, nar + 1, nar - 1]

    if exhaustive:
        pairs = [(a, b) for a in range(1 << n) for b in range(1 << n)]
        singles = list(range(1 << n))
    else:
        random.seed(135)
        rnd = [random.randrange(1 << n) for _ in range(4000)]
        pairs = [(a, b) for a in edges for b in edges]                  # edge x edge
        pairs += [(a, random.randrange(1 << n)) for a in edges for _ in range(200)]
        pairs += [(random.randrange(1 << n), random.randrange(1 << n)) for _ in range(600_000)]
        singles = edges + rnd

    print(f"posit{n}: add/sub/mul/div/compare ({len(pairs):,} pairs)%s ..."
          % (" [exhaustive]" if exhaustive else ""), flush=True)
    for a, b in pairs:
        chk("add", fn["add"](a, b), pt(O["add"](mk(a), mk(b))), a, b)
        chk("sub", fn["sub"](a, b), pt(O["sub"](mk(a), mk(b))), a, b)
        chk("mul", fn["mul"](a, b), pt(O["mul"](mk(a), mk(b))), a, b)
        chk("div", fn["div"](a, b), pt(O["div"](mk(a), mk(b))), a, b)
        chk("lt", fn["lt"](a, b), 1 if O["lt"](mk(a), mk(b)) else 0, a, b)
        chk("le", fn["le"](a, b), 1 if O["le"](mk(a), mk(b)) else 0, a, b)
        chk("eq", fn["eq"](a, b), 1 if O["eq"](mk(a), mk(b)) else 0, a, b)

    print(f"posit{n}: sqrt / nearest_int / neg ({len(singles):,}) ...", flush=True)
    for a in singles:
        chk("sqrt", fn["sqrt"](a), pt(O["sqrt"](mk(a))), a)
        chk("neg", fn["neg"](a), ((1 << n) - a) & M, a)
        if a != ((nar | 1) & M):        # SoftPosit roundToInt quirk at 0x..81
            chk("rint", fn["nearest_int"](a), pt(O["rint"](mk(a))), a)

    print(f"posit{n}: i64->posit over [-3000,3000] ...", flush=True)
    for v in range(-3000, 3001):
        chk("from_i64", fn["from_i64"](v), pt(O["from_i"](v)), v)

    print(f"posit{n}: fma (150k triples) ...", flush=True)
    random.seed(7)
    for _ in range(150_000):
        a, b, c = (random.randrange(1 << n) for _ in range(3))
        chk("fma", fn["fma"](a, b, c), pt(O["fma"](mk(a), mk(b), mk(c))), a, b, c)

    print(f"posit{n}: fdp (5k vectors) ...", flush=True)
    for _ in range(5_000):
        L = random.randint(1, 12)
        av = [random.randrange(1 << n) for _ in range(L)]
        bv = [random.randrange(1 << n) for _ in range(L)]
        q = O["qclr"]()
        for x, y in zip(av, bv): q = O["qadd"](q, mk(x), mk(y))
        exp = pt(O["qtop"](q))
        ca, cb = (u32 * L)(*av), (u32 * L)(*bv)
        chk("fdp", fn["fdp"](ca, cb, L), exp, tuple(av), tuple(bv))

    if not fails:
        print(f"  -> posit{n} PASS\n")
        return 0
    c = Counter(f[0] for f in fails)
    print(f"  -> posit{n} {len(fails)} FAILURES: {dict(c)}")
    for f in fails[:15]:
        print("     %-9s args=%s got=0x%x exp=0x%x" % (f[0], f[1], f[2], f[3]))
    return 1

rc = run(8, True) | run(16, False) | run(32, False)
print("ALL PASS (bit-exact vs SoftPosit pX2)." if rc == 0 else "FAILURES above.")
sys.exit(rc)
