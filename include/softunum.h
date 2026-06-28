//  SoftUnum -- a software implementation of the 2022 Posit Standard (unums).
//
//  Companion to SoftBLAS: where SoftBLAS does software IEEE-754 (via Berkeley
//  SoftFloat), SoftUnum does software *posits/quires/valids* (Type-III unums).
//  It has NO floating-point dependency -- posits are integers, and every value
//  here is a raw bit pattern (an unsigned integer), never a C float/double.
//
//  STANDARD, NOT LEGACY.  The 2022 Posit Standard fixes the exponent size at
//  es = 2 for EVERY width (the exponent field is a 2-bit unsigned integer,
//  0..3), so useed = 2^2^es = 16.  This differs from the 2017 draft and from
//  SoftPosit's *fast* p8/p16 types (es = 0 / es = 1).  Only posit32 coincides
//  between the two conventions; our posit8/posit16 layouts match SoftPosit's
//  GENERIC pX2 (es = 2 at any width), which is our verification oracle.
//
//  SoftUnum is the bit-exact C twin of the pure-Hoon `/lib/unum` (numerics
//  repo, libmath/desk/lib/unum.hoon): each routine runs the IDENTICAL algorithm
//  as its Hoon arm (same decode, same single-rounding encode), so a jet built
//  on SoftUnum produces output bit-identical to the unjetted Hoon.
//
//  Posit bit patterns are passed RIGHT-justified in the low n bits of a
//  uint32_t (the Hoon @-atom convention), NOT left-justified like SoftPosit's
//  internal `posit_2_t`.

#ifndef SOFTUNUM_H
#define SOFTUNUM_H

#include <stdint.h>

//  A posit is its raw n-bit pattern in the low bits of a uint32_t.
typedef uint32_t posit8_t;    //  posit<8,2>,  "byte"   (@rpb)
typedef uint32_t posit16_t;   //  posit<16,2>, "half"   (@rph)
typedef uint32_t posit32_t;   //  posit<32,2>, "single" (@rps)

//  Special bit patterns (per width): zero is all-zero; NaR (Not a Real) is the
//  most-negative two's-complement pattern (1000...0).
#define SU_P8_ZERO    0x00u
#define SU_P8_NAR     0x80u
#define SU_P8_ONE     0x40u
#define SU_P8_MAXPOS  0x7fu
#define SU_P8_MINPOS  0x01u

#define SU_P16_ZERO   0x0000u
#define SU_P16_NAR    0x8000u
#define SU_P16_ONE    0x4000u

#define SU_P32_ZERO   0x00000000u
#define SU_P32_NAR    0x80000000u
#define SU_P32_ONE    0x40000000u

//  =====================================================================
//  posit8  (posit<8,2>)
//  =====================================================================

//  Sign / comparison (raw two's-complement integer ordering of the bits, §5.3).
posit8_t p8_neg(posit8_t a);
posit8_t p8_abs(posit8_t a);
posit8_t p8_sgn(posit8_t a);
int      p8_eq(posit8_t a, posit8_t b);   //  a == b   (bitwise; NaR==NaR true)
int      p8_lt(posit8_t a, posit8_t b);   //  a <  b
int      p8_le(posit8_t a, posit8_t b);   //  a <= b
int      p8_gt(posit8_t a, posit8_t b);   //  a >  b
int      p8_ge(posit8_t a, posit8_t b);   //  a >= b

//  Arithmetic (§5.4): exact g-layer combine, single round.
posit8_t p8_add(posit8_t a, posit8_t b);
posit8_t p8_sub(posit8_t a, posit8_t b);
posit8_t p8_mul(posit8_t a, posit8_t b);
posit8_t p8_div(posit8_t a, posit8_t b);
posit8_t p8_fma(posit8_t a, posit8_t b, posit8_t c);   //  round(a*b + c) once
posit8_t p8_sqrt(posit8_t a);

//  Rounding to integral value.
posit8_t p8_nearest_int(posit8_t a);   //  round-nearest-even
posit8_t p8_floor(posit8_t a);
posit8_t p8_ceil(posit8_t a);

//  Integer conversion.
posit8_t p8_from_u64(uint64_t v);   //  unsigned -> posit
posit8_t p8_from_i64(int64_t v);    //  signed   -> posit
int      p8_to_i64(posit8_t a, int64_t *out);  //  posit -> int; 0 if NaR (none)

//  Quire (§3.4 / §5.11): the 16n-bit exact accumulator and the fused dot
//  product (single rounding).  Vectors are right-justified posit patterns.
posit8_t p8_fdp(const posit8_t *av, const posit8_t *bv, int64_t len);

//  Constants (correctly rounded at the width).
posit8_t p8_one(void);
posit8_t p8_pi(void);
posit8_t p8_tau(void);
posit8_t p8_e(void);
posit8_t p8_phi(void);
posit8_t p8_sqrt2(void);
posit8_t p8_invsqrt2(void);
posit8_t p8_ln2(void);
posit8_t p8_invln2(void);
posit8_t p8_ln10(void);

//  Elementary / transcendental functions (naive reproducible series; accurate
//  only near the expansion point -- not correctly rounded, matching /lib/unum).
posit8_t p8_exp(posit8_t x);
posit8_t p8_sin(posit8_t x);
posit8_t p8_cos(posit8_t x);
posit8_t p8_tan(posit8_t x);
posit8_t p8_log(posit8_t x);      //  natural log
posit8_t p8_log2(posit8_t x);     //  base-2 log
posit8_t p8_log10(posit8_t x);    //  base-10 log
posit8_t p8_pow(posit8_t x, posit8_t y);
posit8_t p8_pow_n(posit8_t x, uint64_t p);
posit8_t p8_factorial(posit8_t x);
posit8_t p8_cbrt(posit8_t x);
posit8_t p8_atan(posit8_t x);
posit8_t p8_asin(posit8_t x);
posit8_t p8_acos(posit8_t x);
int      p8_is_close(posit8_t a, posit8_t b, posit8_t tol);

//  IEEE-754 conversion (value-based: any posit width <-> any float width).
//  to_rd returns the 64-bit pattern; to_rq writes {lo, hi} to out[2].
posit8_t p8_to_rh(posit8_t p);              //  -> binary16 bits
posit8_t p8_to_rs(posit8_t p);              //  -> binary32 bits
uint64_t p8_to_rd(posit8_t p);              //  -> binary64 bits
void     p8_to_rq(posit8_t p, uint64_t out[2]);   //  -> binary128 {lo,hi}
posit8_t p8_from_rh(uint32_t r);            //  binary16  -> posit
posit8_t p8_from_rs(uint32_t r);            //  binary32  -> posit
posit8_t p8_from_rd(uint64_t r);            //  binary64  -> posit
posit8_t p8_from_rq(const uint64_t in[2]);  //  binary128 {lo,hi} -> posit

//  =====================================================================
//  posit16 (posit<16,2>)  and  posit32 (posit<32,2>)
//  =====================================================================
//  Identical surface to posit8, generated from the shared generic core.

posit16_t p16_neg(posit16_t a);
posit16_t p16_abs(posit16_t a);
posit16_t p16_sgn(posit16_t a);
int       p16_eq(posit16_t a, posit16_t b);
int       p16_lt(posit16_t a, posit16_t b);
int       p16_le(posit16_t a, posit16_t b);
int       p16_gt(posit16_t a, posit16_t b);
int       p16_ge(posit16_t a, posit16_t b);
posit16_t p16_add(posit16_t a, posit16_t b);
posit16_t p16_sub(posit16_t a, posit16_t b);
posit16_t p16_mul(posit16_t a, posit16_t b);
posit16_t p16_div(posit16_t a, posit16_t b);
posit16_t p16_fma(posit16_t a, posit16_t b, posit16_t c);
posit16_t p16_sqrt(posit16_t a);
posit16_t p16_nearest_int(posit16_t a);
posit16_t p16_floor(posit16_t a);
posit16_t p16_ceil(posit16_t a);
posit16_t p16_from_u64(uint64_t v);
posit16_t p16_from_i64(int64_t v);
int       p16_to_i64(posit16_t a, int64_t *out);
posit16_t p16_fdp(const posit16_t *av, const posit16_t *bv, int64_t len);
posit16_t p16_one(void);
posit16_t p16_pi(void);
posit16_t p16_tau(void);
posit16_t p16_e(void);
posit16_t p16_phi(void);
posit16_t p16_sqrt2(void);
posit16_t p16_invsqrt2(void);
posit16_t p16_ln2(void);
posit16_t p16_invln2(void);
posit16_t p16_ln10(void);
posit16_t p16_exp(posit16_t x);
posit16_t p16_sin(posit16_t x);
posit16_t p16_cos(posit16_t x);
posit16_t p16_tan(posit16_t x);
posit16_t p16_log(posit16_t x);
posit16_t p16_log2(posit16_t x);
posit16_t p16_log10(posit16_t x);
posit16_t p16_pow(posit16_t x, posit16_t y);
posit16_t p16_pow_n(posit16_t x, uint64_t p);
posit16_t p16_factorial(posit16_t x);
posit16_t p16_cbrt(posit16_t x);
posit16_t p16_atan(posit16_t x);
posit16_t p16_asin(posit16_t x);
posit16_t p16_acos(posit16_t x);
int       p16_is_close(posit16_t a, posit16_t b, posit16_t tol);
posit16_t p16_to_rh(posit16_t p);
posit16_t p16_to_rs(posit16_t p);
uint64_t  p16_to_rd(posit16_t p);
void      p16_to_rq(posit16_t p, uint64_t out[2]);
posit16_t p16_from_rh(uint32_t r);
posit16_t p16_from_rs(uint32_t r);
posit16_t p16_from_rd(uint64_t r);
posit16_t p16_from_rq(const uint64_t in[2]);

posit32_t p32_neg(posit32_t a);
posit32_t p32_abs(posit32_t a);
posit32_t p32_sgn(posit32_t a);
int       p32_eq(posit32_t a, posit32_t b);
int       p32_lt(posit32_t a, posit32_t b);
int       p32_le(posit32_t a, posit32_t b);
int       p32_gt(posit32_t a, posit32_t b);
int       p32_ge(posit32_t a, posit32_t b);
posit32_t p32_add(posit32_t a, posit32_t b);
posit32_t p32_sub(posit32_t a, posit32_t b);
posit32_t p32_mul(posit32_t a, posit32_t b);
posit32_t p32_div(posit32_t a, posit32_t b);
posit32_t p32_fma(posit32_t a, posit32_t b, posit32_t c);
posit32_t p32_sqrt(posit32_t a);
posit32_t p32_nearest_int(posit32_t a);
posit32_t p32_floor(posit32_t a);
posit32_t p32_ceil(posit32_t a);
posit32_t p32_from_u64(uint64_t v);
posit32_t p32_from_i64(int64_t v);
int       p32_to_i64(posit32_t a, int64_t *out);
posit32_t p32_fdp(const posit32_t *av, const posit32_t *bv, int64_t len);
posit32_t p32_one(void);
posit32_t p32_pi(void);
posit32_t p32_tau(void);
posit32_t p32_e(void);
posit32_t p32_phi(void);
posit32_t p32_sqrt2(void);
posit32_t p32_invsqrt2(void);
posit32_t p32_ln2(void);
posit32_t p32_invln2(void);
posit32_t p32_ln10(void);
posit32_t p32_exp(posit32_t x);
posit32_t p32_sin(posit32_t x);
posit32_t p32_cos(posit32_t x);
posit32_t p32_tan(posit32_t x);
posit32_t p32_log(posit32_t x);
posit32_t p32_log2(posit32_t x);
posit32_t p32_log10(posit32_t x);
posit32_t p32_pow(posit32_t x, posit32_t y);
posit32_t p32_pow_n(posit32_t x, uint64_t p);
posit32_t p32_factorial(posit32_t x);
posit32_t p32_cbrt(posit32_t x);
posit32_t p32_atan(posit32_t x);
posit32_t p32_asin(posit32_t x);
posit32_t p32_acos(posit32_t x);
int       p32_is_close(posit32_t a, posit32_t b, posit32_t tol);
posit32_t p32_to_rh(posit32_t p);
posit32_t p32_to_rs(posit32_t p);
uint64_t  p32_to_rd(posit32_t p);
void      p32_to_rq(posit32_t p, uint64_t out[2]);
posit32_t p32_from_rh(uint32_t r);
posit32_t p32_from_rs(uint32_t r);
posit32_t p32_from_rd(uint64_t r);
posit32_t p32_from_rq(const uint64_t in[2]);

#endif  //  SOFTUNUM_H
