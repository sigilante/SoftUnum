# SoftUnum

A software implementation of the **2022 Posit Standard** (Type-III unums:
posits, quires, and eventually valids), in C and eventually Rust.

![](./img/hero.png)

SoftUnum is a companion to [SoftBLAS](https://github.com/urbit/SoftBLAS):
where SoftBLAS does software IEEE-754 linear algebra on top of Berkeley
SoftFloat, SoftUnum does software *posit* arithmetic.

> **Status (2026-06-27):** posit8 / posit16 / posit32 (`posit<{8,16,32},2>`)
> complete and **bit-exact** across the whole core surface — decode/encode,
> sign/compare, add/sub/mul/div/fma/sqrt, round (near/floor/ceil), integer
> conversion, and the 16n-bit quire + fused dot product. posit8 is verified
> *exhaustively* (all 65,536 pairs); posit16/posit32 over ~600k random pairs
> plus a structured edge grid. Oracle: SoftPosit `pX2` for posit8/16, and the
> **dedicated `p32`** path for posit32 (SoftPosit's generic `pX2` misrounds
> tiny values at width 32 — it was only validated to X≈20; our clean-room
> implementation does not have this bug). IEEE-754 conversions and the
> granular quire ops are next, then the vere vendoring + jets.

## Standard, not legacy

The 2022 Posit Standard fixes the exponent size at **`es = 2` for every width**
(the exponent field is a 2-bit unsigned integer, `0..3`), so `useed = 2^2^es =
16`. This differs from the 2017 draft and from SoftPosit's *fast* `p8`/`p16`
types (`es = 0` / `es = 1`). Only `posit32` coincides between the conventions;
our `posit8`/`posit16` layouts match SoftPosit's **generic `pX2`** path (es=2 at
any width) — which is exactly our verification oracle.

## Bit-exact twin of `/lib/unum`

SoftUnum is the C twin of the pure-Hoon `/lib/unum`
(`numerics/libmath/desk/lib/unum.hoon`). Each routine runs the **identical
algorithm** as its Hoon arm — same decode (`sea`), same single-rounding encode
(`bit`), same exact g-layer combine — so a jet built on SoftUnum produces
output *bit-identical* to the unjetted Hoon. (This mirrors how the `math.c`
jets relate to `math.hoon`.) The Hoon arm names are cited in the C source.

Posit bit patterns are passed **right-justified** in the low `n` bits of a
`uint32_t` (the Hoon `@`-atom convention), not left-justified like SoftPosit's
internal `posit_2_t`.

## Per-width native types

Rather than one arbitrary-precision core, each width uses the smallest *fixed*
native integer that covers its worst-case intermediate:

| width | aura | backing | quire | notes |
|---|---|---|---|---|
| `posit8`  | `@rpb` | `unsigned __int128` | 128-bit | everything fits one native int (`p8.c`) |
| `posit16` | `@rph` | 512-bit `wide_t` | 256-bit | generic core (`pcore.h`) |
| `posit32` | `@rps` | 512-bit `wide_t` | 512-bit | generic core; `add`/`sqrt`/quire need multi-word |

posit32's `add` can shift a significand to a common exponent by ~240 bits
(maxpos + minpos) and its quire is 512 bits, so a fixed-size multi-word integer
(`src/posit/pwide.h`, just shift/add/sub/compare + a digit-by-digit `isqt` — no
wide multiply or division) is needed there. posit16 and posit32 share one
generic algorithm (`src/posit/pcore.h`) instantiated per width; posit8 keeps its
own `__int128` copy (`p8.c`) as the readable reference. Still per-width native
fixed-size integers, no general/heap bignum.

## Build & test

```sh
make                       # libsoftunum.a + libsoftunum.{dylib,so}
make test                  # exhaustive bit-exactness vs SoftPosit pX2
```

`make test` needs the `softposit` pip package (the `pX2` es=2 oracle) — the same
oracle `/lib/unum` was validated against. Point `ORACLE_PY` at a Python that has
it (default `/opt/anaconda3/bin/python`).

## Layout

```
include/softunum.h     public API (raw-bit posit types + prototypes)
src/posit/p8.c         posit8 (posit<8,2>) -- __int128 reference copy
src/posit/pwide.h      fixed 512-bit integer (shift/add/sub/cmp/mask/isqt)
src/posit/pcore.h      generic posit<N,2> core, instantiated per width
src/posit/p16.c        posit16 (posit<16,2>)  =  #define PW_N 16 + pcore.h
src/posit/p32.c        posit32 (posit<32,2>)  =  #define PW_N 32 + pcore.h
tools/oracle.py        ctypes harness: SoftUnum vs SoftPosit (pX2 / p32)
```

## Planned: a parallel Rust implementation

A same-behavior Rust implementation is planned (for NockApp integration; not
urgent). Because C and Rust can't share code, "same behavior" will be guaranteed
by a **language-neutral test-vector corpus** generated from the oracle: the C
side and the Rust crate both reproduce `spec/vectors/*` bit-for-bit. The C lib
stays at the repo root; a top-level `rust/` crate will be added alongside it.

## License

MIT. See [LICENSE](LICENSE).
