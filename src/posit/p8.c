//  posit8 (posit<8,2>) -- a bit-exact C twin of `++rpb` in /lib/unum.
//
//  Backed by `unsigned __int128`: posit8's widest intermediate is its 128-bit
//  quire, and `bit`'s pre-round payload stays well within 128 bits (the regime
//  saturates first), so the whole width lives in one native type -- no
//  multi-word helper.  posit16/posit32 widen this same code to u256 / u512.
//
//  Every routine transliterates its Hoon arm (libmath/desk/lib/unum.hoon)
//  line-for-line; the comments cite the Hoon names.

#include "softunum.h"
#include <stdlib.h>   //  llabs

typedef unsigned __int128 u128;

#define N        8
#define MSK      0xffu
#define NARBITS  0x80u
#define MAXPOS   0x7fu
#define MINPOS   0x01u
#define QSCALE   48          //  8n - 16
#define QNAR     ((u128)1 << 127)   //  most-negative 128-bit pattern

//  g-layer decoded posit, mirroring +$up: value = (sign?+1:-1) * a * 2^e.
enum { K_REAL = 0, K_ZERO = 1, K_NAR = 2 };
typedef struct { int kind; int sign; long long e; u128 a; } up_t;

static int u128_bits(u128 x) {           //  bit length (met 0)
  int n = 0;
  while ( x ) { n++; x >>= 1; }
  return n;
}

static long long fdiv4(long long x) {    //  floor(x / 4), matching +bit's rel
  return (x >= 0) ? (x / 4) : -(((-x) + 3) / 4);
}

//  +sea: decode an 8-bit posit into g-layer form.
static up_t sea(unsigned p) {
  up_t u = { K_REAL, 1, 0, 0 };
  p &= MSK;
  if ( p == 0 )       { u.kind = K_ZERO; return u; }
  if ( p == NARBITS ) { u.kind = K_NAR;  return u; }
  int neg = (p >> (N - 1)) & 1;
  unsigned mag = neg ? ((1u << N) - p) : p;
  int pw = N - 1;
  int r0 = (mag >> (pw - 1)) & 1;
  int k = 1, r;
  while ( 1 ) {
    if ( k == pw ) {                     //  all-regime, no exponent/fraction
      r = r0 ? (k - 1) : -k;
      u.sign = !neg; u.e = 4 * r; u.a = 1;
      return u;
    }
    int nb = (mag >> (pw - 1 - k)) & 1;
    if ( nb == r0 ) { k++; continue; }
    break;
  }
  r = r0 ? (k - 1) : -k;
  int remwid = pw - (k + 1);
  unsigned rem = mag & ((1u << remwid) - 1);
  int elo, fw;
  if ( remwid >= 2 ) { elo = rem >> (remwid - 2); fw = remwid - 2; }
  else if ( remwid == 1 ) { elo = rem << 1; fw = 0; }
  else { elo = 0; fw = 0; }
  unsigned frac = rem & ((1u << fw) - 1);
  long long x = 4 * (long long)r + elo;
  u.sign = !neg; u.e = x - fw; u.a = ((u128)1 << fw) + frac;
  return u;
}

//  +smag: apply sign to a magnitude pattern (two's complement if negative).
static unsigned smag(int neg, u128 mag) {
  unsigned m = (unsigned)mag & MSK;
  return neg ? (((1u << N) - m) & MSK) : m;
}

//  +bit: encode g-layer form to an 8-bit posit, round-nearest-even, saturating.
static unsigned bit(up_t u) {
  if ( u.kind == K_ZERO ) return 0;
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.a == 0 )         return 0;
  int neg = !u.sign;
  int lead = u128_bits(u.a) - 1;
  long long x = u.e + lead;
  u128 frac = u.a & (((u128)1 << lead) - 1);
  long long r = fdiv4(x);
  int elo = (int)(x - 4 * r);
  if ( r >= N - 2 )      return smag(neg, MAXPOS);
  if ( r <= -(N - 1) )   return smag(neg, MINPOS);
  u128 regval; int regwid;
  if ( r >= 0 ) { regval = (((u128)1 << (r + 1)) - 1) << 1; regwid = (int)r + 2; }
  else          { regval = 1; regwid = (int)(-r) + 1; }
  int totw = regwid + 2 + lead;
  u128 pay = (regval << (2 + lead)) | ((u128)elo << lead) | frac;
  int pw = N - 1;
  u128 mag;
  if ( totw <= pw ) {
    mag = pay << (pw - totw);
  } else {
    int sh = totw - pw;
    u128 keep = pay >> sh;
    int guard = (int)((pay >> (sh - 1)) & 1);
    int sticky = (pay & (((u128)1 << (sh - 1)) - 1)) != 0;
    int lsbit = (int)(keep & 1);
    int roundup = guard && (sticky || lsbit);
    mag = roundup ? keep + 1 : keep;
    if ( mag > MAXPOS ) mag = MAXPOS;
  }
  return smag(neg, mag);
}

//  ---- sign / comparison ----------------------------------------------------

static long long s_ext(unsigned a) {     //  sign-extend the n-bit pattern
  a &= MSK;
  return (a & NARBITS) ? (long long)a - (1 << N) : (long long)a;
}

posit8_t p8_neg(posit8_t a) { return ((1u << N) - (a & MSK)) & MSK; }

posit8_t p8_abs(posit8_t a) {
  return ((a & MSK) & NARBITS) ? p8_neg(a) : (a & MSK);
}

posit8_t p8_sgn(posit8_t a) {
  a &= MSK;
  if ( a == 0 )       return 0;
  if ( a == NARBITS ) return NARBITS;
  return (a & NARBITS) ? p8_neg(SU_P8_ONE) : SU_P8_ONE;
}

int p8_eq(posit8_t a, posit8_t b) { return (a & MSK) == (b & MSK); }
int p8_lt(posit8_t a, posit8_t b) { return s_ext(a) <  s_ext(b); }
int p8_le(posit8_t a, posit8_t b) { return s_ext(a) <= s_ext(b); }
int p8_gt(posit8_t a, posit8_t b) { return s_ext(a) >  s_ext(b); }
int p8_ge(posit8_t a, posit8_t b) { return s_ext(a) >= s_ext(b); }

//  ---- arithmetic -----------------------------------------------------------

posit8_t p8_mul(posit8_t a, posit8_t b) {
  up_t ua = sea(a), ub = sea(b);
  if ( ua.kind == K_NAR || ub.kind == K_NAR )   return NARBITS;
  if ( ua.kind == K_ZERO || ub.kind == K_ZERO ) return 0;
  up_t r = { K_REAL, (ua.sign == ub.sign), ua.e + ub.e, ua.a * ub.a };
  return bit(r);
}

posit8_t p8_add(posit8_t a, posit8_t b) {
  up_t ua = sea(a), ub = sea(b);
  if ( ua.kind == K_NAR || ub.kind == K_NAR ) return NARBITS;
  if ( ua.kind == K_ZERO ) return b & MSK;
  if ( ub.kind == K_ZERO ) return a & MSK;
  long long emin = (ua.e < ub.e) ? ua.e : ub.e;
  u128 s1 = ua.a << (ua.e - emin);
  u128 s2 = ub.a << (ub.e - emin);
  up_t r = { K_REAL, 1, emin, 0 };
  if ( ua.sign == ub.sign ) { r.sign = ua.sign; r.a = s1 + s2; return bit(r); }
  if ( s1 > s2 )            { r.sign = ua.sign; r.a = s1 - s2; return bit(r); }
  if ( s2 > s1 )            { r.sign = ub.sign; r.a = s2 - s1; return bit(r); }
  return 0;
}

posit8_t p8_sub(posit8_t a, posit8_t b) { return p8_add(a, p8_neg(b)); }

posit8_t p8_div(posit8_t a, posit8_t b) {
  up_t ua = sea(a), ub = sea(b);
  if ( ua.kind == K_NAR || ub.kind == K_NAR ) return NARBITS;
  if ( ub.kind == K_ZERO ) return NARBITS;
  if ( ua.kind == K_ZERO ) return 0;
  int g = 2 * N;
  u128 num = ua.a << g;
  u128 q = num / ub.a;
  u128 qs = (num % ub.a == 0) ? q : (q | 1);
  up_t r = { K_REAL, (ua.sign == ub.sign), ua.e - ub.e - g, qs };
  return bit(r);
}

posit8_t p8_fma(posit8_t a, posit8_t b, posit8_t c) {
  up_t ua = sea(a), ub = sea(b), uc = sea(c);
  if ( ua.kind == K_NAR || ub.kind == K_NAR || uc.kind == K_NAR ) return NARBITS;
  if ( ua.kind == K_ZERO || ub.kind == K_ZERO ) return c & MSK;
  int ps = (ua.sign == ub.sign);
  long long pe = ua.e + ub.e;
  u128 pa = ua.a * ub.a;
  if ( uc.kind == K_ZERO ) { up_t r = { K_REAL, ps, pe, pa }; return bit(r); }
  long long emin = (pe < uc.e) ? pe : uc.e;
  u128 s1 = pa << (pe - emin);
  u128 s2 = uc.a << (uc.e - emin);
  up_t r = { K_REAL, 1, emin, 0 };
  if ( ps == uc.sign ) { r.sign = ps;      r.a = s1 + s2; return bit(r); }
  if ( s1 > s2 )       { r.sign = ps;      r.a = s1 - s2; return bit(r); }
  if ( s2 > s1 )       { r.sign = uc.sign; r.a = s2 - s1; return bit(r); }
  return 0;
}

//  +isqt: integer floor square root (Hoon's Newton variant; floor sqrt is
//  unique, so any correct floor isqt is bit-identical).
static u128 isqt(u128 x) {
  if ( x == 0 ) return 0;
  u128 r = (u128)1 << ((u128_bits(x) + 1) / 2);
  while ( 1 ) {
    u128 nr = (r + x / r) / 2;
    if ( nr >= r ) {
      while ( r * r > x ) r--;
      return r;
    }
    r = nr;
  }
}

posit8_t p8_sqrt(posit8_t a) {
  up_t u = sea(a);
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.kind == K_ZERO ) return 0;
  if ( !u.sign )          return NARBITS;        //  sqrt of a negative is NaR
  int odd = (int)(llabs(u.e) & 1);
  u128 aa = odd ? (u.a << 1) : u.a;
  long long ee = odd ? (u.e - 1) : u.e;
  int g = 2 * N;
  u128 m = aa << (2 * g);
  u128 s = isqt(m);
  u128 sx = (m == s * s) ? s : (s | 1);
  up_t r = { K_REAL, 1, ee / 2 - g, sx };
  return bit(r);
}

//  ---- rounding to integral value ------------------------------------------

static posit8_t p8_round(posit8_t p, int mode) {  //  0 near, 1 down, 2 up
  up_t u = sea(p);
  if ( u.kind != K_REAL ) return p & MSK;
  if ( u.e >= 0 ) return p & MSK;
  int sh = (int)(-u.e);
  u128 hi = u.a >> sh;
  u128 rem = u.a & (((u128)1 << sh) - 1);
  u128 half = (u128)1 << (sh - 1);
  if ( mode == 0 && (rem > half || (rem == half && (hi & 1))) ) hi++;
  if ( mode == 1 && !u.sign && rem != 0 ) hi++;
  if ( mode == 2 &&  u.sign && rem != 0 ) hi++;
  if ( hi == 0 ) return 0;
  up_t r = { K_REAL, u.sign, 0, hi };
  return bit(r);
}

posit8_t p8_nearest_int(posit8_t a) { return p8_round(a, 0); }
posit8_t p8_floor(posit8_t a)       { return p8_round(a, 1); }
posit8_t p8_ceil(posit8_t a)        { return p8_round(a, 2); }

//  ---- integer conversion ---------------------------------------------------

posit8_t p8_from_u64(uint64_t v) {
  if ( v == 0 ) return 0;
  up_t r = { K_REAL, 1, 0, (u128)v };
  return bit(r);
}

posit8_t p8_from_i64(int64_t v) {
  if ( v == 0 ) return 0;
  int neg = v < 0;
  u128 mg = neg ? (u128)(-(unsigned __int128)v) : (u128)v;
  up_t r = { K_REAL, !neg, 0, mg };
  return bit(r);
}

int p8_to_i64(posit8_t p, int64_t *out) {
  up_t u = sea(p);
  if ( u.kind == K_NAR ) return 0;
  up_t ur = sea(p8_round(p, 0));
  if ( ur.kind == K_ZERO ) { *out = 0; return 1; }
  int sh = (int)llabs(ur.e);
  u128 mag = (ur.e >= 0) ? (ur.a << sh) : (ur.a >> sh);
  long long m = (long long)mag;
  *out = ur.sign ? m : -m;
  return 1;
}

//  ---- quire (128-bit) + fused dot product ----------------------------------

static u128 q_mul_add(u128 q, unsigned a, unsigned b) {
  if ( q == QNAR ) return QNAR;
  up_t ua = sea(a), ub = sea(b);
  if ( ua.kind == K_NAR || ub.kind == K_NAR )   return QNAR;
  if ( ua.kind == K_ZERO || ub.kind == K_ZERO ) return q;
  long long sh = ua.e + ub.e + QSCALE;
  u128 m = (ua.a * ub.a) << llabs(sh);
  u128 qc = (ua.sign == ub.sign) ? m : ((u128)0 - m);
  return q + qc;            //  mod 2^128 wraps naturally
}

static posit8_t q_to_p(u128 q) {
  if ( q == QNAR ) return NARBITS;
  int neg = (int)(q >> 127) & 1;
  u128 acc = neg ? ((u128)0 - q) : q;
  if ( acc == 0 ) return 0;
  up_t r = { K_REAL, !neg, -QSCALE, acc };
  return bit(r);
}

posit8_t p8_fdp(const posit8_t *av, const posit8_t *bv, int64_t len) {
  u128 q = 0;
  for ( int64_t i = 0; i < len; i++ ) q = q_mul_add(q, av[i], bv[i]);
  return q_to_p(q);
}

//  ---- granular quire ops ---------------------------------------------------
//  The quire is exposed as QWORDS=2 little-endian uint64 words (128 bits).
static u128 p_to_q(unsigned p) {
  up_t u = sea(p);
  if ( u.kind == K_NAR )  return QNAR;
  if ( u.kind == K_ZERO ) return 0;
  u128 m = u.a << llabs(u.e + QSCALE);
  return u.sign ? m : ((u128)0 - m);
}
static u128 q_neg_(u128 q) { return (q == QNAR) ? QNAR : ((u128)0 - q); }
static u128 q_addq_(u128 x, u128 y) { return (x == QNAR || y == QNAR) ? QNAR : (x + y); }
static inline u128 q_load(const uint64_t *w) { return (u128)w[0] | ((u128)w[1] << 64); }
static inline void q_store(uint64_t *w, u128 q) { w[0] = (uint64_t)q; w[1] = (uint64_t)(q >> 64); }

void p8_q_zero(uint64_t *q) { q[0] = 0; q[1] = 0; }
void p8_q_nar(uint64_t *q)  { q_store(q, QNAR); }
int  p8_q_is_nar(const uint64_t *q) { return q_load(q) == QNAR; }
void p8_p_to_q(posit8_t p, uint64_t *q) { q_store(q, p_to_q(p)); }
posit8_t p8_q_to_p(const uint64_t *q)   { return q_to_p(q_load(q)); }
void p8_q_mul_add(uint64_t *q, posit8_t a, posit8_t b) { q_store(q, q_mul_add(q_load(q), a, b)); }
void p8_q_mul_sub(uint64_t *q, posit8_t a, posit8_t b) { q_store(q, q_mul_add(q_load(q), a, p8_neg(b))); }
void p8_q_add_p(uint64_t *q, posit8_t p) { q_store(q, q_mul_add(q_load(q), p, NARBITS >> 1)); }
void p8_q_sub_p(uint64_t *q, posit8_t p) { q_store(q, q_mul_add(q_load(q), p8_neg(p), NARBITS >> 1)); }
void p8_q_add_q(uint64_t *x, const uint64_t *y) { q_store(x, q_addq_(q_load(x), q_load(y))); }
void p8_q_sub_q(uint64_t *x, const uint64_t *y) { q_store(x, q_addq_(q_load(x), q_neg_(q_load(y)))); }
void p8_q_negate(uint64_t *q) { q_store(q, q_neg_(q_load(q))); }

//  ---- IEEE-754 conversion (value-based, any posit width <-> any float) ------
#include "pieee.h"

static fnt up_to_fn(up_t u) {
  fnt f = { FN_FIN, u.sign, u.e, u.a };
  if ( u.kind == K_NAR )  { f.kind = FN_NAN; f.sign = 1; f.e = 0; f.a = 0; }
  else if ( u.kind == K_ZERO ) { f.sign = 1; f.e = 0; f.a = 0; }
  return f;
}
static up_t fn_to_up(fnt f) {
  up_t u = { K_REAL, f.sign, f.e, f.a };
  if ( f.kind != FN_FIN ) { u.kind = K_NAR; u.sign = 1; u.e = 0; u.a = 0; }
  else if ( f.a == 0 )    { u.kind = K_ZERO; u.sign = 1; u.e = 0; u.a = 0; }
  return u;
}
posit8_t p8_to_rh(posit8_t p) { return (posit8_t)f_bit(up_to_fn(sea(p)), IEEE_RH); }
posit8_t p8_to_rs(posit8_t p) { return (posit8_t)f_bit(up_to_fn(sea(p)), IEEE_RS); }
uint64_t p8_to_rd(posit8_t p) { return (uint64_t)f_bit(up_to_fn(sea(p)), IEEE_RD); }
void p8_to_rq(posit8_t p, uint64_t out[2]) {
  u128i b = f_bit(up_to_fn(sea(p)), IEEE_RQ); out[0] = (uint64_t)b; out[1] = (uint64_t)(b >> 64);
}
posit8_t p8_from_rh(uint32_t r) { return bit(fn_to_up(f_sea(r, IEEE_RH))); }
posit8_t p8_from_rs(uint32_t r) { return bit(fn_to_up(f_sea(r, IEEE_RS))); }
posit8_t p8_from_rd(uint64_t r) { return bit(fn_to_up(f_sea((u128i)r, IEEE_RD))); }
posit8_t p8_from_rq(const uint64_t in[2]) {
  return bit(fn_to_up(f_sea(((u128i)in[1] << 64) | in[0], IEEE_RQ)));
}

//  ---- elementary / transcendental functions (shared body) ------------------
//  +pconst: encode a constant given as Q(.)52 fixed-point significand `a` with
//  binary exponent `e` -- one rounding through +bit, matching the Hoon arms.
static uint32_t pconst(long long e, uint64_t a) {
  up_t u = { K_REAL, 1, e, (u128)a };
  return bit(u);
}

#define PW_PREFIX p8
#include "ptrans.h"
