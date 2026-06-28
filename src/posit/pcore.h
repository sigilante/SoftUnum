//  pcore.h -- the generic posit<N,2> algorithm, instantiated per width by
//  #include from p16.c / p32.c with PW_N (the bit width) and PW_PREFIX (the
//  function-name prefix, e.g. p16) defined.  It is the bit-exact twin of the
//  Hoon `++pp` core in /lib/unum: same decode (sea), same single-rounding
//  encode (bit), same exact g-layer combine.  posit8 has its own __int128
//  copy (p8.c) which doubles as the readable reference; this body widens that
//  same logic to the 512-bit wide_t for posit16/posit32.
//
//  Requires (before include): PW_N, PW_PREFIX.

#include "softunum.h"
#include "pwide.h"
#include <stdlib.h>   //  llabs

#define PW_CAT2(a, b) a ## b
#define PW_CAT(a, b)  PW_CAT2(a, b)
#define PFN(name)     PW_CAT(PW_PREFIX, name)

#define N        PW_N
#define MSK      ((PW_N >= 32) ? 0xffffffffu : ((1u << PW_N) - 1u))
#define NARBITS  (1u << (PW_N - 1))
#define MAXPOS   (NARBITS - 1u)
#define MINPOS   1u
#define QSCALE   (8 * PW_N - 16)
#define QBITS    (16 * PW_N)

enum { K_REAL = 0, K_ZERO = 1, K_NAR = 2 };
typedef struct { int kind; int sign; long long e; wide_t a; } up_t;

static long long fdiv4(long long x) { return (x >= 0) ? (x / 4) : -(((-x) + 3) / 4); }

//  +sea: decode an N-bit posit (right-justified) into g-layer form.
static up_t sea(uint32_t p) {
  up_t u = { K_REAL, 1, 0, { { 0 } } };
  p &= MSK;
  if ( p == 0 )       { u.kind = K_ZERO; return u; }
  if ( p == NARBITS ) { u.kind = K_NAR;  return u; }
  uint64_t base = (uint64_t)1 << N;
  int neg = (p >> (N - 1)) & 1;
  uint32_t mag = neg ? (uint32_t)(base - p) : p;
  int pw = N - 1;
  int r0 = (mag >> (pw - 1)) & 1;
  int k = 1, r;
  while ( 1 ) {
    if ( k == pw ) {
      r = r0 ? (k - 1) : -k;
      u.sign = !neg; u.e = 4 * r; u.a = w_from_u64(1);
      return u;
    }
    int nb = (mag >> (pw - 1 - k)) & 1;
    if ( nb == r0 ) { k++; continue; }
    break;
  }
  r = r0 ? (k - 1) : -k;
  int remwid = pw - (k + 1);
  uint32_t rem = mag & ((1u << remwid) - 1);
  int elo, fw;
  if ( remwid >= 2 ) { elo = rem >> (remwid - 2); fw = remwid - 2; }
  else if ( remwid == 1 ) { elo = rem << 1; fw = 0; }
  else { elo = 0; fw = 0; }
  uint32_t frac = rem & ((1u << fw) - 1);
  long long x = 4 * (long long)r + elo;
  u.sign = !neg; u.e = x - fw; u.a = w_from_u64(((uint64_t)1 << fw) + frac);
  return u;
}

static uint32_t smag(int neg, uint64_t mag) {
  uint64_t m = mag & MSK;
  uint64_t base = (uint64_t)1 << N;
  return neg ? (uint32_t)((base - m) & MSK) : (uint32_t)m;
}

//  +bit: encode g-layer form to an N-bit posit, round-nearest-even, saturating.
static uint32_t bit(up_t u) {
  if ( u.kind == K_ZERO ) return 0;
  if ( u.kind == K_NAR )  return NARBITS;
  if ( w_is_zero(u.a) )   return 0;
  int neg = !u.sign;
  int lead = w_bits(u.a) - 1;
  long long x = u.e + lead;
  wide_t frac = w_masklow(u.a, lead);
  long long r = fdiv4(x);
  int elo = (int)(x - 4 * r);
  if ( r >= N - 2 )    return smag(neg, MAXPOS);
  if ( r <= -(N - 1) ) return smag(neg, MINPOS);
  wide_t regval; int regwid;
  if ( r >= 0 ) { regval = w_shl(w_from_u64(((uint64_t)1 << (r + 1)) - 1), 1); regwid = (int)r + 2; }
  else          { regval = w_from_u64(1); regwid = (int)(-r) + 1; }
  int totw = regwid + 2 + lead;
  wide_t pay = w_or(w_or(w_shl(regval, 2 + lead),
                         w_shl(w_from_u64((uint64_t)elo), lead)), frac);
  int pw = N - 1;
  uint64_t mag;
  if ( totw <= pw ) {
    mag = (uint64_t)w_to_u128(w_shl(pay, pw - totw));
  } else {
    int sh = totw - pw;
    wide_t keepw = w_shr(pay, sh);
    int guard = w_testbit(pay, sh - 1);
    int sticky = !w_is_zero(w_masklow(pay, sh - 1));
    int lsbit = w_testbit(keepw, 0);
    int roundup = guard && (sticky || lsbit);
    uint64_t keep = (uint64_t)w_to_u128(keepw);
    mag = roundup ? keep + 1 : keep;
    if ( mag > MAXPOS ) mag = MAXPOS;
  }
  return smag(neg, mag);
}

//  ---- sign / comparison ----------------------------------------------------

static long long s_ext(uint32_t a) {
  a &= MSK;
  return (a & NARBITS) ? (long long)a - ((long long)1 << N) : (long long)a;
}

uint32_t PFN(_neg)(uint32_t a) { uint64_t base = (uint64_t)1 << N; return (uint32_t)((base - (a & MSK)) & MSK); }
uint32_t PFN(_abs)(uint32_t a) { return ((a & MSK) & NARBITS) ? PFN(_neg)(a) : (a & MSK); }
uint32_t PFN(_sgn)(uint32_t a) {
  a &= MSK;
  if ( a == 0 )       return 0;
  if ( a == NARBITS ) return NARBITS;
  uint32_t one = NARBITS >> 1;     //  posit 1.0 = 0x40.. = NARBITS/2
  return (a & NARBITS) ? PFN(_neg)(one) : one;
}
int PFN(_eq)(uint32_t a, uint32_t b) { return (a & MSK) == (b & MSK); }
int PFN(_lt)(uint32_t a, uint32_t b) { return s_ext(a) <  s_ext(b); }
int PFN(_le)(uint32_t a, uint32_t b) { return s_ext(a) <= s_ext(b); }
int PFN(_gt)(uint32_t a, uint32_t b) { return s_ext(a) >  s_ext(b); }
int PFN(_ge)(uint32_t a, uint32_t b) { return s_ext(a) >= s_ext(b); }

//  ---- arithmetic -----------------------------------------------------------

uint32_t PFN(_mul)(uint32_t a, uint32_t b) {
  up_t ua = sea(a), ub = sea(b);
  if ( ua.kind == K_NAR || ub.kind == K_NAR )   return NARBITS;
  if ( ua.kind == K_ZERO || ub.kind == K_ZERO ) return 0;
  up_t r = { K_REAL, (ua.sign == ub.sign), ua.e + ub.e,
             w_from_u128(w_to_u128(ua.a) * w_to_u128(ub.a)) };
  return bit(r);
}

uint32_t PFN(_add)(uint32_t a, uint32_t b) {
  up_t ua = sea(a), ub = sea(b);
  if ( ua.kind == K_NAR || ub.kind == K_NAR ) return NARBITS;
  if ( ua.kind == K_ZERO ) return b & MSK;
  if ( ub.kind == K_ZERO ) return a & MSK;
  long long emin = (ua.e < ub.e) ? ua.e : ub.e;
  wide_t s1 = w_shl(ua.a, (int)(ua.e - emin));
  wide_t s2 = w_shl(ub.a, (int)(ub.e - emin));
  up_t r = { K_REAL, 1, emin, { { 0 } } };
  if ( ua.sign == ub.sign ) { r.sign = ua.sign; r.a = w_add(s1, s2); return bit(r); }
  int c = w_cmp(s1, s2);
  if ( c > 0 ) { r.sign = ua.sign; r.a = w_sub(s1, s2); return bit(r); }
  if ( c < 0 ) { r.sign = ub.sign; r.a = w_sub(s2, s1); return bit(r); }
  return 0;
}

uint32_t PFN(_sub)(uint32_t a, uint32_t b) { return PFN(_add)(a, PFN(_neg)(b)); }

uint32_t PFN(_div)(uint32_t a, uint32_t b) {
  up_t ua = sea(a), ub = sea(b);
  if ( ua.kind == K_NAR || ub.kind == K_NAR ) return NARBITS;
  if ( ub.kind == K_ZERO ) return NARBITS;
  if ( ua.kind == K_ZERO ) return 0;
  int g = 2 * N;
  unsigned __int128 num = w_to_u128(ua.a) << g;     //  a<=~28 bits, g<=64 -> fits
  unsigned __int128 den = w_to_u128(ub.a);
  unsigned __int128 q = num / den;
  unsigned __int128 qs = (num % den == 0) ? q : (q | 1);
  up_t r = { K_REAL, (ua.sign == ub.sign), ua.e - ub.e - g, w_from_u128(qs) };
  return bit(r);
}

uint32_t PFN(_fma)(uint32_t a, uint32_t b, uint32_t c) {
  up_t ua = sea(a), ub = sea(b), uc = sea(c);
  if ( ua.kind == K_NAR || ub.kind == K_NAR || uc.kind == K_NAR ) return NARBITS;
  if ( ua.kind == K_ZERO || ub.kind == K_ZERO ) return c & MSK;
  int ps = (ua.sign == ub.sign);
  long long pe = ua.e + ub.e;
  wide_t pa = w_from_u128(w_to_u128(ua.a) * w_to_u128(ub.a));
  if ( uc.kind == K_ZERO ) { up_t r = { K_REAL, ps, pe, pa }; return bit(r); }
  long long emin = (pe < uc.e) ? pe : uc.e;
  wide_t s1 = w_shl(pa, (int)(pe - emin));
  wide_t s2 = w_shl(uc.a, (int)(uc.e - emin));
  up_t r = { K_REAL, 1, emin, { { 0 } } };
  if ( ps == uc.sign ) { r.sign = ps; r.a = w_add(s1, s2); return bit(r); }
  int cc = w_cmp(s1, s2);
  if ( cc > 0 ) { r.sign = ps;      r.a = w_sub(s1, s2); return bit(r); }
  if ( cc < 0 ) { r.sign = uc.sign; r.a = w_sub(s2, s1); return bit(r); }
  return 0;
}

uint32_t PFN(_sqrt)(uint32_t a) {
  up_t u = sea(a);
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.kind == K_ZERO ) return 0;
  if ( !u.sign )          return NARBITS;
  int odd = (int)(llabs(u.e) & 1);
  wide_t aa = odd ? w_shl(u.a, 1) : u.a;
  long long ee = odd ? (u.e - 1) : u.e;
  int g = 2 * N;
  wide_t m = w_shl(aa, 2 * g);
  int exact;
  wide_t s = w_isqt(m, &exact);
  wide_t sx = exact ? s : w_or(s, w_from_u64(1));
  up_t r = { K_REAL, 1, ee / 2 - g, sx };
  return bit(r);
}

//  ---- rounding to integral value ------------------------------------------

static uint32_t pround(uint32_t p, int mode) {   //  0 near, 1 down, 2 up
  up_t u = sea(p);
  if ( u.kind != K_REAL ) return p & MSK;
  if ( u.e >= 0 ) return p & MSK;
  int sh = (int)(-u.e);
  wide_t hi = w_shr(u.a, sh);
  wide_t rem = w_masklow(u.a, sh);
  wide_t half = w_shl(w_from_u64(1), sh - 1);
  int chalf = w_cmp(rem, half);
  if ( mode == 0 && (chalf > 0 || (chalf == 0 && w_testbit(hi, 0))) ) hi = w_add(hi, w_from_u64(1));
  if ( mode == 1 && !u.sign && !w_is_zero(rem) ) hi = w_add(hi, w_from_u64(1));
  if ( mode == 2 &&  u.sign && !w_is_zero(rem) ) hi = w_add(hi, w_from_u64(1));
  if ( w_is_zero(hi) ) return 0;
  up_t r = { K_REAL, u.sign, 0, hi };
  return bit(r);
}

uint32_t PFN(_nearest_int)(uint32_t a) { return pround(a, 0); }
uint32_t PFN(_floor)(uint32_t a)       { return pround(a, 1); }
uint32_t PFN(_ceil)(uint32_t a)        { return pround(a, 2); }

//  ---- integer conversion ---------------------------------------------------

uint32_t PFN(_from_u64)(uint64_t v) {
  if ( v == 0 ) return 0;
  up_t r = { K_REAL, 1, 0, w_from_u64(v) };
  return bit(r);
}
uint32_t PFN(_from_i64)(int64_t v) {
  if ( v == 0 ) return 0;
  int neg = v < 0;
  uint64_t mg = neg ? (uint64_t)(-(unsigned __int128)v) : (uint64_t)v;
  up_t r = { K_REAL, !neg, 0, w_from_u64(mg) };
  return bit(r);
}
int PFN(_to_i64)(uint32_t p, int64_t *out) {
  up_t u = sea(p);
  if ( u.kind == K_NAR ) return 0;
  up_t ur = sea(pround(p, 0));
  if ( ur.kind == K_ZERO ) { *out = 0; return 1; }
  int sh = (int)llabs(ur.e);
  wide_t mag = (ur.e >= 0) ? w_shl(ur.a, sh) : w_shr(ur.a, sh);
  long long m = (long long)w_to_u128(mag);
  *out = ur.sign ? m : -m;
  return 1;
}

//  ---- quire (QBITS-bit) + fused dot product --------------------------------

static wide_t q_nar(void) { return w_shl(w_from_u64(1), QBITS - 1); }

static wide_t q_mul_add(wide_t q, uint32_t a, uint32_t b) {
  if ( w_cmp(q, q_nar()) == 0 ) return q_nar();
  up_t ua = sea(a), ub = sea(b);
  if ( ua.kind == K_NAR || ub.kind == K_NAR )   return q_nar();
  if ( ua.kind == K_ZERO || ub.kind == K_ZERO ) return q;
  long long sh = ua.e + ub.e + QSCALE;
  wide_t m = w_shl(w_from_u128(w_to_u128(ua.a) * w_to_u128(ub.a)), (int)llabs(sh));
  wide_t qc = (ua.sign == ub.sign) ? m : w_sub(w_zero(), m);
  return w_masklow(w_add(q, qc), QBITS);
}

static uint32_t q_to_p(wide_t q) {
  q = w_masklow(q, QBITS);
  if ( w_cmp(q, q_nar()) == 0 ) return NARBITS;
  int neg = w_testbit(q, QBITS - 1);
  wide_t acc = neg ? w_masklow(w_sub(w_zero(), q), QBITS) : q;
  if ( w_is_zero(acc) ) return 0;
  up_t r = { K_REAL, !neg, -(long long)QSCALE, acc };
  return bit(r);
}

uint32_t PFN(_fdp)(const uint32_t *av, const uint32_t *bv, int64_t len) {
  wide_t q = w_zero();
  for ( int64_t i = 0; i < len; i++ ) q = q_mul_add(q, av[i], bv[i]);
  return q_to_p(q);
}

//  ---- elementary / transcendental functions (shared body) ------------------
//  +pconst: encode a constant given as a fixed-point significand `a` with
//  binary exponent `e` -- one rounding through +bit, matching the Hoon arms.
static uint32_t pconst(long long e, uint64_t a) {
  up_t u = { K_REAL, 1, e, w_from_u64(a) };
  return bit(u);
}

#include "ptrans.h"
