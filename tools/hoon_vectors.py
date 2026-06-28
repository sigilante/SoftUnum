#!/usr/bin/env python3
"""Check SoftUnum against the curated /lib/unum Hoon test vectors.

These are the exact expected values from numerics/libmath/desk/tests/lib/
unum-{fns,edge}.hoon -- authoritative output of the pure-Hoon /lib/unum
(cross-checked against mpmath in the Hoon comments).  This oracle is
INDEPENDENT of SoftPosit (which has no transcendentals) and independent of
oracle.py's Python transliteration, so agreement here is strong evidence the
C transcendentals are bit-exact to the Hoon spec.

    /opt/anaconda3/bin/python tools/hoon_vectors.py    # no deps but ctypes
"""
import ctypes, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "libsoftunum.dylib")
if not os.path.exists(LIB): LIB = os.path.join(HERE, "..", "libsoftunum.so")
lib = ctypes.CDLL(LIB)
u32, u64, ci = ctypes.c_uint32, ctypes.c_uint64, ctypes.c_int

def call(w, op, args):
    #  resolve symbolic args: ('sun',N)->from_u64, ('nsun',N)->neg(from_u64)
    raw = []
    for a in args:
        if isinstance(a, tuple) and a[0] == "sun":
            raw.append(getattr(lib, f"p{w}_from_u64")(u64(a[1])) & 0xffffffff)
        elif isinstance(a, tuple) and a[0] == "nsun":
            s = getattr(lib, f"p{w}_from_u64")(u64(a[1])) & 0xffffffff
            raw.append(getattr(lib, f"p{w}_neg")(s) & 0xffffffff)
        else:
            raw.append(a)
    f = getattr(lib, f"p{w}_{op}")
    if op == "pow_n":
        f.argtypes = [u32, u64]; f.restype = u32; return f(raw[0], raw[1]) & 0xffffffff
    f.argtypes = [u32] * len(raw); f.restype = u32
    return f(*raw) & 0xffffffff

#  (width, op, [args], expected)  -- transcribed from unum-fns / unum-edge.hoon
V = [
  # ---- posit8 (test-transcendental-rpb, -domain-rpb, -new-rpb, -domain-new) --
  (8,"exp",[0x0],0x40),(8,"exp",[0x40],0x4b),(8,"exp",[0x80],0x80),
  (8,"sin",[0x0],0x0),(8,"sin",[0x40],0x3d),
  (8,"cos",[0x0],0x40),(8,"cos",[0x40],0x39),(8,"tan",[0x40],0x44),
  (8,"log",[0x40],0x0),(8,"log",[("sun",2)],0x3b),
  (8,"log",[0x0],0x80),(8,"log",[0xc0],0x80),(8,"log",[0x80],0x80),(8,"log",[0xb4],0x80),
  (8,"pow_n",[("sun",2),3],0x58),(8,"pow_n",[0x80,0],0x80),(8,"pow_n",[0x40,0],0x40),
  (8,"pow",[("sun",2),0x0],0x40),(8,"pow",[0x0,("sun",2)],0x80),(8,"pow",[0xc0,("sun",2)],0x80),
  (8,"factorial",[("sun",3)],0x54),(8,"factorial",[("sun",4)],0x62),(8,"factorial",[0xb4],0x80),
  (8,"cbrt",[("sun",1)],0x40),(8,"cbrt",[0xc0],0x80),(8,"cbrt",[0x0],0x0),
  (8,"atan",[("sun",1)],0x3d),(8,"atan",[0x38],0x38),
  (8,"asin",[0x38],0x39),(8,"asin",[("sun",4)],0x80),(8,"asin",[("sun",1)],0x45),
  (8,"acos",[0x38],0x41),(8,"acos",[0x0],0x45),(8,"acos",[("sun",4)],0x80),(8,"acos",[("sun",1)],0x0),
  # ---- posit16 (test-transcendental-rph) ------------------------------------
  (16,"exp",[0x3800],0x4531),(16,"sin",[0x3800],0x3757),(16,"cos",[0x3800],0x3e0b),
  (16,"tan",[0x3800],0x38bd),(16,"log",[("sun",2)],0x3b18),
  (16,"factorial",[("sun",3)],0x5400),(16,"factorial",[("sun",4)],0x6200),
  (16,"cbrt",[("sun",1)],0x4000),(16,"atan",[("sun",1)],0x3c90),(16,"atan",[0x3800],0x36d5),
  (16,"asin",[0x3800],0x3861),(16,"acos",[0x3800],0x4060),(16,"acos",[0x0],0x4491),
  (16,"pow_n",[("sun",2),3],0x5800),
  # ---- posit32 (test-transcendental-rps + unum-edge) ------------------------
  (32,"exp",[0x38000000],0x453094c8),(32,"sin",[0x38000000],0x3757743a),
  (32,"cos",[0x38000000],0x3e0a9404),(32,"tan",[0x38000000],0x38bda7ad),
  (32,"log",[("sun",2)],0x3b172180),
  (32,"factorial",[("sun",3)],0x54000000),(32,"factorial",[("sun",4)],0x62000000),
  (32,"cbrt",[("sun",1)],0x40000000),(32,"atan",[("sun",1)],0x3c90fdab),
  (32,"atan",[0x38000000],0x36d63381),(32,"asin",[0x38000000],0x3860a91c),
  (32,"acos",[0x38000000],0x4060a91d),(32,"acos",[0x0],0x4490fdaa),
  (32,"pow_n",[("sun",2),3],0x58000000),
  (32,"exp",[0x0],0x40000000),(32,"cos",[0x0],0x40000000),(32,"sin",[0x0],0x0),
  (32,"tan",[0x0],0x0),(32,"log",[("sun",1)],0x0),
  (32,"factorial",[0x0],0x40000000),(32,"factorial",[("sun",5)],0x6b800000),
  (32,"atan",[0x0],0x0),(32,"asin",[0x0],0x0),(32,"acos",[("sun",1)],0x0),
  (32,"exp",[("nsun",1)],0x33c5ab1c),(32,"atan",[("nsun",1)],0xc36f0255),
  (32,"exp",[0x80000000],0x80000000),(32,"sin",[0x80000000],0x80000000),
  (32,"factorial",[0x80000000],0x80000000),
  (32,"exp",[("sun",10)],0x7a579ded),(32,"exp",[("sun",50)],0x7ffe1b02),
  (32,"exp",[("sun",100)],0x7ffff02b),(32,"sin",[("sun",10)],0xc74ba64a),
  (32,"log",[("sun",100)],0x50e9b7f7),
  # ---- sqrt edge identities (unum-edge) -------------------------------------
  (32,"sqrt",[0x0],0x0),(32,"sqrt",[("sun",1)],0x40000000),(32,"sqrt",[("sun",4)],0x48000000),
  (32,"sqrt",[0x80000000],0x80000000),(32,"sqrt",[("nsun",1)],0x80000000),
  (8,"sqrt",[0xc0],0x80),
]

fails = []
for w, op, args, exp in V:
    got = call(w, op, args)
    if got != exp: fails.append((w, op, args, got, exp))

print(f"checked {len(V)} Hoon /lib/unum vectors across posit8/16/32")
if not fails:
    print("ALL Hoon vectors PASS (SoftUnum bit-exact to /lib/unum).")
    sys.exit(0)
print(f"{len(fails)} FAILURES:")
for w, op, args, got, exp in fails:
    print(f"  p{w}_{op}{tuple(args)}: got=0x{got:x} exp=0x{exp:x}")
sys.exit(1)
