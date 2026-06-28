# SoftUnum

A software implementation of the **2022 Posit Standard** (Type-III unums:
posits, quires, and — eventually — valids), in **pure-integer C**.

SoftUnum is the companion to [SoftBLAS](https://github.com/urbit/SoftBLAS):
where SoftBLAS does software IEEE-754 linear algebra on top of Berkeley
SoftFloat, SoftUnum does software *posit* arithmetic. Crucially, **SoftUnum has
no floating-point dependency at all** — a posit is an integer, and every value
here is a raw bit pattern (an unsigned integer), never a C `float`/`double`.
There is no SoftFloat submodule and nothing to link.

> **Status (2026-06-27):** posit8 (`posit<8,2>`) complete and **bit-exact**
> across its whole surface — verified exhaustively against SoftPosit `pX2`
> (all 65,536 pairs for add/sub/mul/div/compare, all 256 values for
> sqrt/round/neg, integer conversion, plus sampled fma and quire/fdp).
> posit16 and posit32 are next.

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
| `posit8`  | `@rpb` | `unsigned __int128` | 128-bit | everything fits one native int |
| `posit16` | `@rph` | `__int128` + a small `u256` | 256-bit | *planned* |
| `posit32` | `@rps` | `u256` / `u512` helpers | 512-bit | `add`/`sqrt`/quire need multi-word; *planned* |

posit32's `add` can shift a significand to a common exponent by ~240 bits
(maxpos + minpos) and its quire is 512 bits, so a couple of fixed-size
multi-word helpers are needed there — still per-width, still no general bignum.

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
src/posit/p8.c         posit8 (posit<8,2>)
tools/oracle.py        ctypes harness: SoftUnum vs SoftPosit pX2
```

## Planned: a parallel Rust implementation

A same-behavior Rust implementation is planned (for NockApp integration; not
urgent). Because C and Rust can't share code, "same behavior" will be guaranteed
by a **language-neutral test-vector corpus** generated from the oracle: the C
side and the Rust crate both reproduce `spec/vectors/*` bit-for-bit. The C lib
stays at the repo root; a top-level `rust/` crate will be added alongside it.

## License

MIT. See [LICENSE](LICENSE).
