#!/usr/bin/env python3
"""Verify SoftUnum's granular quire ops against an exact big-integer model.

The quire is a 16n-bit two's-complement fixed-point accumulator with scale
2^-qscale (qscale = 8n-16): a posit value v contributes the exact integer
v * 2^qscale.  That makes a from-scratch Python model trivial and exact -- an
oracle independent of SoftPosit.  We run random op sequences on both the C
quire (via its uint64-word buffer) and the model, comparing the raw words after
every step plus the final q_to_p, and also check the curated Hoon
test-quire-rpb vectors.

    /opt/anaconda3/bin/python tools/quire_check.py
"""
import ctypes, os, sys, random
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "libsoftunum.dylib")
if not os.path.exists(LIB): LIB = os.path.join(HERE, "..", "libsoftunum.so")
lib = ctypes.CDLL(LIB)
u32, u64 = ctypes.c_uint32, ctypes.c_uint64
P = ctypes.POINTER(u64)

#  --- exact posit value + /lib/unum encoder ---------------------------------
def pvalue(p, n):
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

def pencode(neg, e, a, n):
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

class QModel:
    """Exact quire model: a Python int in [0, 2^qbits) (two's complement),
    or the string 'n' for q-NaR."""
    def __init__(self, n):
        self.n = n; self.qbits = 16 * n; self.qscale = 8 * n - 16
        self.qmod = 1 << self.qbits; self.qnar = 1 << (self.qbits - 1)
        self.one = 1 << (n - 2); self.M = (1 << n) - 1
        self.q = 0
    def neg(self, p): return ((1 << self.n) - p) & self.M
    def contrib(self, a, b):                  #  exact integer for value(a)*value(b)*2^qscale
        va, vb = pvalue(a, self.n), pvalue(b, self.n)
        if va == "n" or vb == "n": return "n"
        v = va * vb * (Fraction(2) ** self.qscale)
        assert v.denominator == 1
        return int(v) % self.qmod
    def zero(self): self.q = 0
    def nar(self): self.q = "n"
    def p_to_q(self, p):
        v = pvalue(p, self.n)
        if v == "n": self.q = "n"; return
        self.q = int(v * (Fraction(2) ** self.qscale)) % self.qmod
    def mul_add(self, a, b):
        if self.q == "n": return
        c = self.contrib(a, b)
        if c == "n": self.q = "n"
        else: self.q = (self.q + c) % self.qmod
    def mul_sub(self, a, b): self.mul_add(a, self.neg(b))
    def add_p(self, p): self.mul_add(p, self.one)
    def sub_p(self, p): self.mul_add(self.neg(p), self.one)
    def negate(self):
        if self.q != "n": self.q = (self.qmod - self.q) % self.qmod
    def add_q(self, y):
        if self.q == "n" or y == "n": self.q = "n"
        else: self.q = (self.q + y) % self.qmod
    def words(self, qwords):
        q = self.qnar if self.q == "n" else self.q
        return [(q >> (64 * i)) & 0xffffffffffffffff for i in range(qwords)]
    def to_p(self):
        posit_nar = 1 << (self.n - 1)
        if self.q == "n": return posit_nar
        q = self.q & (self.qmod - 1)
        if q == self.qnar: return posit_nar
        neg = (q >> (self.qbits - 1)) & 1
        acc = (self.qmod - q) if neg else q
        if acc == 0: return 0
        return pencode(neg, -self.qscale, acc, self.n)

fails = []
def run(n, rounds):
    qwords = (16 * n) // 64
    pfx = f"p{n}"
    F = {nm: getattr(lib, f"{pfx}_q_{nm}") for nm in
         ("zero", "nar", "mul_add", "mul_sub", "add_p", "sub_p", "add_q", "sub_q", "negate", "to_p")}
    F["zero"].argtypes = [P]; F["nar"].argtypes = [P]; F["negate"].argtypes = [P]
    F["to_p"].argtypes = [P]; F["to_p"].restype = u32
    for nm in ("mul_add", "mul_sub"): F[nm].argtypes = [P, u32, u32]
    for nm in ("add_p", "sub_p"): F[nm].argtypes = [P, u32]
    for nm in ("add_q", "sub_q"): F[nm].argtypes = [P, P]
    ptoq = getattr(lib, f"{pfx}_p_to_q"); ptoq.argtypes = [u32, P]
    isn = getattr(lib, f"{pfx}_q_is_nar"); isn.argtypes = [P]; isn.restype = ctypes.c_int

    def buf(): return (u64 * qwords)()
    def cmp(tag, b, m, p=None):
        cw = list(b); mw = m.words(qwords)
        if cw != mw: fails.append((n, tag, p, [hex(x) for x in cw], [hex(x) for x in mw]))

    nar = 1 << (n - 1)
    pool = [0, nar, 1 << (n - 2), (1 << (n - 1)) - 1, 1, 2, 3, nar + 1]
    rng = random.Random(n * 7 + 1)

    #  p_to_q / q_to_p round-trip
    rt = range(1 << n) if n == 8 else pool + [rng.randrange(1 << n) for _ in range(3000)]
    for p in rt:
        b = buf(); ptoq(p, b)
        m = QModel(n); m.p_to_q(p)
        cmp("p_to_q", b, m, p)
        if F["to_p"](b) != m.to_p(): fails.append((n, "q_to_p", p, hex(F["to_p"](b)), hex(m.to_p())))

    #  random op sequences
    for _ in range(rounds):
        b = buf(); F["zero"](b); m = QModel(n); m.zero()
        for _ in range(rng.randint(1, 30)):
            op = rng.choice(["mul_add", "mul_sub", "add_p", "sub_p", "negate", "add_q", "sub_q", "nar"])
            if op in ("mul_add", "mul_sub"):
                a, c = rng.choice(pool + [rng.randrange(1 << n)]), rng.choice(pool + [rng.randrange(1 << n)])
                F[op](b, a, c); getattr(m, op)(a, c)
            elif op in ("add_p", "sub_p"):
                p = rng.choice(pool + [rng.randrange(1 << n)]); F[op](b, p); getattr(m, op)(p)
            elif op == "negate":
                F[op](b); m.negate()
            elif op in ("add_q", "sub_q"):
                y = QModel(n); y.zero()
                for _ in range(rng.randint(0, 4)):
                    y.mul_add(rng.choice(pool + [rng.randrange(1 << n)]), rng.choice(pool + [rng.randrange(1 << n)]))
                yb = buf()
                for i, wv in enumerate(y.words(qwords)): yb[i] = wv
                if op == "add_q": F[op](b, yb); m.add_q(y.q if y.q == "n" else y.q)
                else: F[op](b, yb); m.add_q("n" if y.q == "n" else (m.qmod - y.q) % m.qmod)
            else:
                F["nar"](b); m.nar()
            cmp(op, b, m)
        if F["to_p"](b) != m.to_p(): fails.append((n, "seq-to_p", None, hex(F["to_p"](b)), hex(m.to_p())))

#  curated Hoon test-quire-rpb (p8) -- composed ops, exact expected values
def hoon_p8():
    b = (u64 * 2)()
    lib.p8_p_to_q.argtypes = [u32, P]; lib.p8_q_to_p.argtypes = [P]; lib.p8_q_to_p.restype = u32
    lib.p8_q_add_p.argtypes = [P, u32]; lib.p8_q_mul_add.argtypes = [P, u32, u32]
    lib.p8_q_zero.argtypes = [P]
    cases = []
    lib.p8_p_to_q(0x42, b); cases.append(("rt42", lib.p8_q_to_p(b), 0x42))
    lib.p8_p_to_q(0x38, b); cases.append(("rt38", lib.p8_q_to_p(b), 0x38))
    lib.p8_p_to_q(0x40, b); lib.p8_q_add_p(b, 0x48); cases.append(("addp", lib.p8_q_to_p(b), 0x4c))
    lib.p8_q_zero(b); lib.p8_q_mul_add(b, 0x48, 0x4c); cases.append(("muladd", lib.p8_q_to_p(b), 0x54))
    for tag, got, exp in cases:
        if got != exp: fails.append((8, "hoon-" + tag, None, hex(got), hex(exp)))

for n, rounds in ((8, 4000), (16, 4000), (32, 4000)):
    print(f"posit{n}: quire round-trip + {rounds} random op sequences ...", flush=True)
    run(n, rounds)
hoon_p8()

if not fails:
    print("\nALL quire checks PASS (vs exact big-integer model + Hoon vectors).")
    sys.exit(0)
from collections import Counter
print(f"\n{len(fails)} FAILURES: {dict(Counter(f[1] for f in fails))}")
for f in fails[:20]: print("  ", f)
sys.exit(1)
