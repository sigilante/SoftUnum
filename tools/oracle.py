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
    for nm in ("exp", "sin", "cos", "tan", "log", "log2", "log10",
               "cbrt", "factorial", "atan", "asin", "acos"):
        f = getattr(lib, f"{pfx}_{nm}"); f.argtypes = [u32]; f.restype = u32; fns[nm] = f
    f = getattr(lib, f"{pfx}_pow"); f.argtypes = [u32, u32]; f.restype = u32; fns["pow"] = f
    f = getattr(lib, f"{pfx}_pow_n"); f.argtypes = [u32, u64]; f.restype = u32; fns["pow_n"] = f
    f = getattr(lib, f"{pfx}_is_close"); f.argtypes = [u32, u32, u32]; f.restype = ci; fns["is_close"] = f
    for nm in ("one", "pi", "tau", "e", "phi", "sqrt2", "invsqrt2", "ln2", "invln2", "ln10"):
        f = getattr(lib, f"{pfx}_{nm}"); f.argtypes = []; f.restype = u32; fns[nm] = f
    return fns

#  Independent g-layer encoder (mirrors /lib/unum +bit) -- the constants oracle.
def encode(neg, e, a, n):
    if a == 0: return 0
    M = (1 << n) - 1; maxpos = (1 << (n - 1)) - 1
    lead = a.bit_length() - 1; x = e + lead; frac = a & ((1 << lead) - 1)
    r = x >> 2; elo = x - 4 * r
    if r >= n - 2: return ((1 << n) - maxpos) & M if neg else maxpos
    if r <= -(n - 1): return ((1 << n) - 1) & M if neg else 1
    if r >= 0: regval = ((1 << (r + 1)) - 1) << 1; regwid = r + 2
    else: regval = 1; regwid = -r + 1
    totw = regwid + 2 + lead
    pay = (regval << (2 + lead)) | (elo << lead) | frac
    pw = n - 1
    if totw <= pw:
        mag = pay << (pw - totw)
    else:
        sh = totw - pw; keep = pay >> sh
        guard = (pay >> (sh - 1)) & 1; low = pay & ((1 << (sh - 1)) - 1)
        if guard and ((1 if low else 0) or (keep & 1)): keep += 1
        if keep > maxpos: keep = maxpos
        mag = keep
    return ((1 << n) - mag) & M if neg else mag

#  Reference transcendentals: transliterate the Hoon /lib/unum arms over
#  SoftPosit posit arithmetic (the underlying ops are already verified above).
#  SoftPosit has NO transcendentals, so the Hoon structure IS the spec.
def run_trans(n):
    pfx = {8: "p8", 16: "p16", 32: "p32"}[n]
    fn = bind(pfx); O = oracle(n); mk, pt = O["mk"], O["pt"]
    M, nar = (1 << n) - 1, 1 << (n - 1)
    ADD = lambda a, b: pt(O["add"](mk(a), mk(b)))
    SUB = lambda a, b: pt(O["sub"](mk(a), mk(b)))
    MUL = lambda a, b: pt(O["mul"](mk(a), mk(b)))
    DIV = lambda a, b: pt(O["div"](mk(a), mk(b)))
    SQT = lambda a: pt(O["sqrt"](mk(a)))
    NEG = lambda a: ((1 << n) - a) & M
    ABS = lambda a: NEG(a) if (a & nar) else a
    LT = lambda a, b: bool(O["lt"](mk(a), mk(b)))
    LE = lambda a, b: bool(O["le"](mk(a), mk(b)))
    EQ = lambda a, b: bool(O["eq"](mk(a), mk(b)))
    SUN = lambda v: pt(O["from_i"](v))
    ONE = encode(False, 0, 1, n)
    PI = encode(False, -52, 0x3243f6a8885a31, n)
    LN2 = encode(False, -52, 0x0b17217f7d1cf8, n)
    LN10 = encode(False, -52, 0x24d763776aaa2b, n)
    consts = {"one": ONE, "pi": PI, "tau": encode(False, -52, 0x6487ed5110b461, n),
              "e": encode(False, -52, 0x2b7e151628aed3, n),
              "phi": encode(False, -52, 0x19e3779b97f4a8, n),
              "sqrt2": encode(False, -52, 0x16a09e667f3bcd, n),
              "invsqrt2": encode(False, -52, 0x0b504f333f9de6, n),
              "ln2": LN2, "invln2": encode(False, -52, 0x171547652b82fe, n), "ln10": LN10}

    def r_exp(x):
        s = t = ONE
        for k in range(1, 21): t = MUL(t, DIV(x, SUN(k))); s = ADD(s, t)
        return s
    def r_sin(x):
        t = s = x; x2 = MUL(x, x)
        for nn in range(1, 21):
            k = 2 * nn; t = NEG(MUL(t, DIV(x2, MUL(SUN(k), SUN(k + 1))))); s = ADD(s, t)
        return s
    def r_cos(x):
        t = s = ONE; x2 = MUL(x, x)
        for nn in range(1, 21):
            k = 2 * nn; t = NEG(MUL(t, DIV(x2, MUL(SUN(k - 1), SUN(k))))); s = ADD(s, t)
        return s
    def r_tan(x): return DIV(r_sin(x), r_cos(x))
    def r_log(x):
        if LE(x, 0): return nar
        y = DIV(SUB(x, ONE), ADD(x, ONE)); y2 = MUL(y, y); s = t = y
        for nn in range(1, 31):
            t = MUL(t, y2); s = ADD(s, MUL(DIV(ONE, SUN(2 * nn + 1)), t))
        return MUL(SUN(2), s)
    def r_log2(x): return DIV(r_log(x), LN2)
    def r_log10(x): return DIV(r_log(x), LN10)
    def r_pow(x, y): return r_exp(MUL(y, r_log(x)))
    def r_pow_n(x, p):
        if x == nar: return nar
        res = ONE
        while p: res = MUL(res, x); p -= 1
        return res
    def r_factorial(x):
        if x == nar or LT(x, 0): return nar
        t = ONE
        while not LE(x, ONE): t = MUL(t, x); x = SUB(x, ONE)
        return t
    def r_cbrt(x):
        if x == nar: return nar
        if x == 0: return 0
        if LT(x, 0): return nar
        return r_pow(x, DIV(ONE, SUN(3)))
    def r_atan(x):
        if x == nar: return nar
        rt = SQT(ADD(ONE, MUL(x, x))); a = DIV(ONE, rt); b = ONE
        for _ in range(41):
            ai = MUL(DIV(ONE, SUN(2)), ADD(a, b)); b = SQT(MUL(ai, b)); a = ai
        return DIV(x, MUL(rt, b))
    def r_asin(x):
        if x == nar: return nar
        if LT(ABS(x), ONE): return r_atan(DIV(x, SQT(SUB(ONE, MUL(x, x)))))
        if EQ(x, ONE): return MUL(PI, DIV(ONE, SUN(2)))
        if EQ(x, NEG(ONE)): return NEG(MUL(PI, DIV(ONE, SUN(2))))
        return nar
    def r_acos(x):
        if x == nar: return nar
        if LT(ABS(x), ONE):
            if EQ(x, 0): return MUL(PI, DIV(ONE, SUN(2)))
            return r_atan(DIV(SQT(SUB(ONE, MUL(x, x))), x))
        if EQ(x, ONE): return 0
        if EQ(x, NEG(ONE)): return PI
        return nar
    #  NB: factorial is excluded from the random sweep -- its naive loop
    #  (while x>1: x-=1) does not terminate for large x (huge-1 rounds back to
    #  huge in posit arithmetic), exactly as in /lib/unum.  It is checked below
    #  on its safe small-integer domain only.
    refs = {"exp": r_exp, "sin": r_sin, "cos": r_cos, "tan": r_tan, "log": r_log,
            "log2": r_log2, "log10": r_log10, "cbrt": r_cbrt,
            "atan": r_atan, "asin": r_asin, "acos": r_acos}

    import random as _r
    _r.seed(99)
    one_, pi_ = ONE, PI
    edges = [0, nar, one_, NEG(one_), pi_, encode(False, -1, 1, n), encode(False, 1, 1, n),
             2, 3, encode(False, -3, 5, n)]
    vals = edges + [_r.randrange(1 << n) for _ in range(200)]
    pvals = edges + [_r.randrange(1 << n) for _ in range(40)]   # cheaper for pow
    fails = []
    def chk(name, got, exp, *a):
        if got != exp: fails.append((name, a, got, exp))

    for nm, rf in refs.items():
        for x in vals: chk(nm, fn[nm](x), rf(x), x)
    #  factorial: safe domain only (small non-negative integers, plus the
    #  fast-exit cases: NaR and negatives -> NaR immediately).
    fvals = [encode(False, 0, k, n) for k in range(0, 11)] + [0, nar, NEG(ONE), encode(False, -1, 1, n)]
    for x in fvals: chk("factorial", fn["factorial"](x), r_factorial(x), x)
    for x in pvals:
        for y in pvals:
            chk("pow", fn["pow"](x, y), r_pow(x, y), x, y)
    for x in vals:
        for p in (0, 1, 2, 3, 7):
            chk("pow_n", fn["pow_n"](x, p), r_pow_n(x, p), x, p)
    for nm, want in consts.items():
        chk(nm, fn[nm](), want)

    if not fails:
        print(f"  -> posit{n} transcendentals PASS\n"); return 0
    c = Counter(f[0] for f in fails)
    print(f"  -> posit{n} transcendentals {len(fails)} FAILURES: {dict(c)}")
    for f in fails[:15]:
        print("     %-10s args=%s got=0x%x exp=0x%x" % (f[0], f[1], f[2], f[3]))
    return 1

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

rc = 0
for _n, _ex in [(8, True), (16, False), (32, False)]:
    rc |= run(_n, _ex)
    rc |= run_trans(_n)
print("ALL PASS (bit-exact vs SoftPosit + Hoon-faithful transcendental reference)."
      if rc == 0 else "FAILURES above.")
sys.exit(rc)
