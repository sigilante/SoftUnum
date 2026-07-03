//  pgmp.h -- exact arbitrary-precision g-layer arithmetic (GMP mpz_t), the
//  bit-exact C twin of unum.hoon's +gmul/+gadd/+gneg/+gsub/+gdiv/+gpoly/
//  +g-round/+glt.
//
//  WHY GMP: the new range-reduced transcendentals in ptrans.h decode a posit
//  ONCE (mirroring +sea), then do EXACT (never-rounded) arithmetic on a
//  [sign exponent significand] triple, rounding only once at the very end.
//  Hoon `@` atoms are arbitrary precision, so a Horner-chain-accumulated
//  significand (e.g. the degree-16 poly in +lr, used by +log) can grow to
//  thousands of bits -- +gmul/+gadd never normalize (see the comment on
//  +gdiv's `g + y-bit-length` shift in unum.hoon for why this is
//  intentional).  The existing fixed-512-bit `wide_t` (pwide.h) is NOT
//  enough for this; truncating mid-computation would be a DIFFERENT
//  algorithm needing its own precision-margin proof.  GMP's mpz_t gives us
//  the same "never truncate until +bit" property C already gets in Hoon.
//
//  gval_t mirrors Hoon's +$up (%p case only -- %z/%n are handled by the
//  caller before/after touching gval_t, exactly as the Hoon arms do):
//  value = (sign ? +1 : -1) * a * 2^e, `sign=1` meaning NON-NEGATIVE (the
//  Hoon +$up / +si convention: "s=%.y is non-negative").  A gval_t with
//  a==0 stands in for Hoon's up = [%p %.y --0 0] (the exact-zero case
//  +gadd's cancellation branch produces) -- g_bit/g_round/glt all treat
//  a==0 the same way the Hoon arms special-case `=(0 a.x)`.
//
//  OWNERSHIP CONVENTION: every gval_t the caller creates (via gval_new /
//  gval_from_* / gval_const, or receives as the return value of any g*
//  function below) owns its own `mpz_t a` and must eventually be cleared
//  with gval_clear.  Functions that take gval_t arguments BY VALUE borrow
//  them -- they read but never clear/mutate the caller's mpz internals, so
//  the caller must still separately clear every gval_t it created.  (This
//  bookkeeping is the price of exact bignum in C; Hoon's purely functional
//  +up values need none of it.)
//
//  DESIGN NOTE (deviation from a literal `g_bit`-via-existing-`bit()` port):
//  the existing per-width `bit()` in p8.c/pcore.h takes a FIXED-width
//  significand (u128 / 512-bit wide_t).  Bridging a possibly-thousands-of-
//  bits gval_t into that fixed width would require either (a) truncating
//  and losing exactness, needing a sticky-bit-preserving-compression
//  correctness proof of its own, or (b) growing wide_t to match every
//  algorithm's worst case, which changes with the polynomial degree.
//  Instead, g_bit() below is a direct GMP transliteration of Hoon's own
//  +bit arm (round-to-nearest-even, saturating) operating on the
//  arbitrary-precision mpz_t significand directly -- zero truncation risk,
//  and it *is* the "existing +bit arm", just implemented once generically
//  over `n` instead of once per fixed width.  (No `up_from_gval` bridge is
//  provided for the same reason: routing back through the fixed-width
//  `up_t`/`bit()` for an arbitrary-size gval_t would be the unsafe
//  truncation this note describes.)
//
//  DESIGN NOTE (g_round precision): Hoon's +g-round returns an arbitrary-
//  precision `@s`.  For posit32, |x| can reach ~2^120 (maxpos), so e.g.
//  `q = round(ax * 2/pi)` in +sc can itself need ~120 bits -- genuinely
//  more than a 64-bit `long long` can hold.  g_round() therefore writes its
//  *exact* result into a caller-supplied mpz_t, never truncating.  The one
//  place a plain scalar is needed (+exp's `(sum:si e.p k)` final-exponent
//  combine) uses ll_add_clamp() below, which is provably safe to saturate:
//  g_bit()'s own saturation thresholds only compare against +-n (tiny, <=32),
//  so ANY sufficiently-huge (but sign-correct) clamped magnitude drives the
//  identical maxpos/minpos branch the true, possibly-huge k would.

#ifndef SOFTUNUM_PGMP_H
#define SOFTUNUM_PGMP_H

#include <gmp.h>
#include <stdint.h>
#include "pwide.h"     //  wide_t, for gval_from_wide's bridge from pcore.h's up_t

//  sign: 1 = non-negative (Hoon's up.s convention), 0 = negative.
typedef struct { int sign; long long e; mpz_t a; } gval_t;

static inline gval_t gval_new(void) {
  gval_t g; g.sign = 1; g.e = 0; mpz_init(g.a); return g;
}
static inline void gval_clear(gval_t *g) { mpz_clear(g->a); }
static inline gval_t gval_copy(gval_t x) {
  gval_t r = gval_new(); r.sign = x.sign; r.e = x.e; mpz_set(r.a, x.a); return r;
}

//  ---- bridging from the existing fixed-width up_t ---------------------
//  p8.c's up_t.a is `unsigned __int128`; pcore.h's (p16/p32) up_t.a is the
//  512-bit `wide_t`.  Each includer defines UP_TO_GVAL(u) using whichever of
//  these matches its own up_t before #include "ptrans.h" (see p8.c / pcore.h).
static inline gval_t gval_from_u128(int sign, long long e, unsigned __int128 a) {
  gval_t g = gval_new(); g.sign = sign; g.e = e;
  uint64_t limbs[2] = { (uint64_t)a, (uint64_t)(a >> 64) };
  mpz_import(g.a, 2, -1, sizeof(uint64_t), 0, 0, limbs);
  return g;
}
static inline gval_t gval_from_wide(int sign, long long e, wide_t a) {
  gval_t g = gval_new(); g.sign = sign; g.e = e;
  mpz_import(g.a, WLIMBS, -1, sizeof(uint64_t), 0, 0, a.w);
  return g;
}

//  Build a gval_t constant from a plain (non-negative) hex significand
//  string (no "0x" prefix, no separators) and a binary exponent -- the C
//  equivalent of Hoon's `[%p %.y e 0x...]` constant literals.  `sign`
//  applies to the OVERALL constant (used for the alternating +/- terms in
//  the sin/cos/atan Taylor tables).
static inline gval_t gval_const(int sign, long long e, const char *hex) {
  gval_t g = gval_new(); g.sign = sign; g.e = e;
  mpz_set_str(g.a, hex, 16);
  return g;
}
static inline gval_t gval_small(int sign, long long e, unsigned long a) {
  gval_t g = gval_new(); g.sign = sign; g.e = e; mpz_set_ui(g.a, a);
  return g;
}

//  ---- +gmul / +gadd / +gneg / +gsub / +gdiv / +gpoly / +g-round / +glt ----

static inline gval_t gmul(gval_t x, gval_t y) {
  gval_t r = gval_new();
  r.sign = (x.sign == y.sign);
  r.e = x.e + y.e;
  mpz_mul(r.a, x.a, y.a);
  return r;
}

static inline gval_t gadd(gval_t x, gval_t y) {
  if ( mpz_sgn(x.a) == 0 ) return gval_copy(y);
  if ( mpz_sgn(y.a) == 0 ) return gval_copy(x);
  long long emin = (x.e < y.e) ? x.e : y.e;
  mpz_t s1, s2;
  mpz_init(s1); mpz_init(s2);
  mpz_mul_2exp(s1, x.a, (unsigned long)(x.e - emin));
  mpz_mul_2exp(s2, y.a, (unsigned long)(y.e - emin));
  gval_t r = gval_new(); r.e = emin;
  if ( x.sign == y.sign ) {
    r.sign = x.sign;
    mpz_add(r.a, s1, s2);
  } else {
    int c = mpz_cmp(s1, s2);
    if ( c > 0 )      { r.sign = x.sign; mpz_sub(r.a, s1, s2); }
    else if ( c < 0 ) { r.sign = y.sign; mpz_sub(r.a, s2, s1); }
    else              { r.sign = 1; mpz_set_ui(r.a, 0); }
  }
  mpz_clear(s1); mpz_clear(s2);
  return r;
}

static inline gval_t gneg(gval_t x) {
  gval_t r = gval_new(); r.sign = !x.sign; r.e = x.e; mpz_set(r.a, x.a); return r;
}

static inline gval_t gsub(gval_t x, gval_t y) {
  gval_t ny = gneg(y);
  gval_t r = gadd(x, ny);
  gval_clear(&ny);
  return r;
}

//  +gdiv:  x/y truncated to g bits of EXTRA precision beyond x's own (not
//  exact -- g is a target-width margin).  Shift is g+(bit-width of y), NOT
//  bare g: +gmul/+gadd never normalize, so a Horner-accumulated x/y pair
//  (e.g. tan's sin/cos) can have huge, unrelated raw bit-lengths -- a bare
//  (x<<g) would silently produce a zero-bit/garbage quotient.  See the
//  +gdiv comment in unum.hoon for the full story (this bug bit the Hoon
//  draft first).
static inline gval_t gdiv(gval_t x, gval_t y, long g) {
  gval_t r = gval_new();
  if ( mpz_sgn(x.a) == 0 ) { r.sign = 1; r.e = 0; return r; }
  unsigned long ybits = (mpz_sgn(y.a) == 0) ? 0 : mpz_sizeinbase(y.a, 2);
  unsigned long shift = (unsigned long)g + ybits;
  mpz_t num;
  mpz_init(num);
  mpz_mul_2exp(num, x.a, shift);
  mpz_tdiv_q(r.a, num, y.a);
  mpz_clear(num);
  r.sign = (x.sign == y.sign);
  r.e = x.e - y.e - (long long)shift;
  return r;
}

//  +gpoly:  Horner over `cs` (highest degree first, cs[0]..cs[n-1]) at `z`.
static inline gval_t gpoly(const gval_t *cs, int n, gval_t z) {
  gval_t acc = gval_copy(cs[0]);
  for ( int i = 1; i < n; i++ ) {
    gval_t m = gmul(acc, z);
    gval_t s = gadd(m, cs[i]);
    gval_clear(&acc); gval_clear(&m);
    acc = s;
  }
  return acc;
}

//  +g-round:  exact nearest integer of x, ties away from zero.  Writes the
//  exact SIGNED result into `out` (caller-owned, already mpz_init'd) -- see
//  the module comment for why this can't be a bounded `long long`.
static inline void g_round(mpz_t out, gval_t x) {
  if ( mpz_sgn(x.a) == 0 ) { mpz_set_ui(out, 0); return; }
  if ( x.e >= 0 ) {
    mpz_mul_2exp(out, x.a, (unsigned long)x.e);
  } else {
    unsigned long sh = (unsigned long)(-x.e);
    mpz_t rm, half;
    mpz_init(rm); mpz_init(half);
    mpz_tdiv_q_2exp(out, x.a, sh);
    mpz_tdiv_r_2exp(rm, x.a, sh);
    mpz_ui_pow_ui(half, 2, sh - 1);
    if ( mpz_cmp(rm, half) >= 0 ) mpz_add_ui(out, out, 1);
    mpz_clear(rm); mpz_clear(half);
  }
  if ( !x.sign ) mpz_neg(out, out);
}

//  Wrap a signed exact-integer result of g_round (or any signed mpz) into
//  the [sign, e=0, a=|k|] gval_t form every caller immediately builds
//  (Hoon's `kup`/`(True,0,qn)` pattern).
static inline gval_t gval_from_signed_mpz(const mpz_t k) {
  gval_t r = gval_new();
  r.sign = (mpz_sgn(k) >= 0);
  r.e = 0;
  mpz_abs(r.a, k);
  return r;
}

//  Add a (possibly astronomically large) exact signed integer `k` to a
//  small `long long` exponent `e`, saturating to +-lim.  Safe specifically
//  because g_bit()'s saturation test only compares against +-n (n<=32): any
//  clamped magnitude >= lim (chosen far larger than any real n) drives the
//  same maxpos/minpos branch the true sum would.  Used by +exp's final
//  `(sum:si e.p k)` combine, the one place a raw (unwrapped) big k meets a
//  scalar exponent.
static inline long long ll_add_clamp(long long e, const mpz_t k, long long lim) {
  mpz_t t;
  mpz_init_set_si(t, e);
  mpz_add(t, t, k);
  long long res;
  if ( mpz_cmp_si(t, lim) > 0 )       res = lim;
  else if ( mpz_cmp_si(t, -lim) < 0 ) res = -lim;
  else                                 res = mpz_get_si(t);
  mpz_clear(t);
  return res;
}

//  +glt:  exact g-layer x<y compare (shared reduction helper).
static inline int glt(gval_t x, gval_t y) {
  gval_t d = gsub(x, y);
  int r = ( mpz_sgn(d.a) == 0 ) ? 0 : !d.sign;
  gval_clear(&d);
  return r;
}

//  ---- g_bit:  arbitrary-precision +bit (round-to-nearest-even, saturating)
//  Direct GMP transliteration of unum.hoon's +bit, generic over the target
//  width `n` (posit8/16/32 all call this with their own N).  See the module
//  comment for why this exists instead of reusing the fixed-width bit().

static inline long long g_bit_fdiv4(long long x) {
  return (x >= 0) ? (x / 4) : -(((-x) + 3) / 4);
}

static inline uint32_t g_bit_smag(int neg, uint32_t mag, int n) {
  uint64_t msk = (n >= 32) ? 0xffffffffull : ((1ull << n) - 1ull);
  if ( !neg ) return (uint32_t)((uint64_t)mag & msk);
  uint64_t base = 1ull << n;
  return (uint32_t)((base - (uint64_t)mag) & msk);
}

static inline uint32_t g_bit(gval_t x, int n) {
  if ( mpz_sgn(x.a) == 0 ) return 0;
  int neg = !x.sign;
  unsigned long lead = (unsigned long)(mpz_sizeinbase(x.a, 2) - 1);
  long long xexp = x.e + (long long)lead;
  mpz_t frac;
  mpz_init(frac);
  mpz_tdiv_r_2exp(frac, x.a, lead);
  long long r = g_bit_fdiv4(xexp);
  int elo = (int)(xexp - 4 * r);
  uint32_t narbits = 1u << (n - 1);
  uint32_t maxpos = narbits - 1u;
  uint32_t minpos = 1u;
  if ( r >= n - 2 )    { mpz_clear(frac); return g_bit_smag(neg, maxpos, n); }
  if ( r <= -(n - 1) ) { mpz_clear(frac); return g_bit_smag(neg, minpos, n); }
  mpz_t regval;
  mpz_init(regval);
  long long regwid;
  if ( r >= 0 ) {
    mpz_ui_pow_ui(regval, 2, (unsigned long)(r + 1));
    mpz_sub_ui(regval, regval, 1);
    mpz_mul_2exp(regval, regval, 1);
    regwid = r + 2;
  } else {
    mpz_set_ui(regval, 1);
    regwid = -r + 1;
  }
  long long totw = regwid + 2 + (long long)lead;
  mpz_t pay, tmp;
  mpz_init(pay); mpz_init(tmp);
  mpz_mul_2exp(pay, regval, (unsigned long)(2 + lead));
  mpz_set_ui(tmp, (unsigned long)elo);
  mpz_mul_2exp(tmp, tmp, lead);
  mpz_ior(pay, pay, tmp);
  mpz_ior(pay, pay, frac);
  mpz_clear(tmp); mpz_clear(frac); mpz_clear(regval);
  long long pw = n - 1;
  uint32_t mag;
  if ( totw <= pw ) {
    mpz_mul_2exp(pay, pay, (unsigned long)(pw - totw));
    mag = (uint32_t)mpz_get_ui(pay);
  } else {
    long long sh = totw - pw;
    mpz_t keep, sticky_r;
    mpz_init(keep); mpz_init(sticky_r);
    mpz_tdiv_q_2exp(keep, pay, (unsigned long)sh);
    int guard = mpz_tstbit(pay, (unsigned long)(sh - 1));
    mpz_tdiv_r_2exp(sticky_r, pay, (unsigned long)(sh - 1));
    int sticky = mpz_sgn(sticky_r) != 0;
    int lsbit = mpz_tstbit(keep, 0);
    int roundup = guard && (sticky || lsbit);
    if ( roundup ) mpz_add_ui(keep, keep, 1);
    if ( mpz_cmp_ui(keep, maxpos) > 0 ) mpz_set_ui(keep, maxpos);
    mag = (uint32_t)mpz_get_ui(keep);
    mpz_clear(keep); mpz_clear(sticky_r);
  }
  mpz_clear(pay);
  return g_bit_smag(neg, mag, n);
}

#endif  //  SOFTUNUM_PGMP_H
