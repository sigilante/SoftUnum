//  pwide.h -- a fixed-width (512-bit) unsigned integer for the posit16/posit32
//  cores.  Just the operations the encode/decode/arithmetic/quire algorithm
//  needs: shift, add, sub, compare, bit-length, low-mask, test-bit, and a
//  bit-by-bit integer square root.  No wide multiply and no wide division --
//  products of significands stay within __int128, and the digit-by-digit isqt
//  needs only shifts/adds/subs.  512 bits covers posit32's 512-bit quire (and
//  hence everything narrower); posit16 masks down to 256 where it matters.

#ifndef SOFTUNUM_PWIDE_H
#define SOFTUNUM_PWIDE_H

#include <stdint.h>

#define WLIMBS 8                       //  8 * 64 = 512 bits
#define WBITS  (64 * WLIMBS)

typedef struct { uint64_t w[WLIMBS]; } wide_t;

static inline wide_t w_zero(void) {
  wide_t r; for ( int i = 0; i < WLIMBS; i++ ) r.w[i] = 0; return r;
}
static inline wide_t w_from_u64(uint64_t x) {
  wide_t r = w_zero(); r.w[0] = x; return r;
}
static inline wide_t w_from_u128(unsigned __int128 x) {
  wide_t r = w_zero(); r.w[0] = (uint64_t)x; r.w[1] = (uint64_t)(x >> 64); return r;
}
static inline unsigned __int128 w_to_u128(wide_t a) {
  return ((unsigned __int128)a.w[1] << 64) | a.w[0];
}
static inline int w_is_zero(wide_t a) {
  for ( int i = 0; i < WLIMBS; i++ ) if ( a.w[i] ) return 0;
  return 1;
}
static inline int w_bits(wide_t a) {            //  bit length (0 for zero)
  for ( int i = WLIMBS - 1; i >= 0; i-- )
    if ( a.w[i] ) return i * 64 + (64 - __builtin_clzll(a.w[i]));
  return 0;
}
static inline int w_cmp(wide_t a, wide_t b) {   //  -1 / 0 / 1
  for ( int i = WLIMBS - 1; i >= 0; i-- )
    if ( a.w[i] != b.w[i] ) return a.w[i] < b.w[i] ? -1 : 1;
  return 0;
}
static inline int w_testbit(wide_t a, int i) {
  if ( i < 0 || i >= WBITS ) return 0;
  return (a.w[i >> 6] >> (i & 63)) & 1;
}
static inline wide_t w_or(wide_t a, wide_t b) {
  wide_t r; for ( int i = 0; i < WLIMBS; i++ ) r.w[i] = a.w[i] | b.w[i]; return r;
}
static inline wide_t w_add(wide_t a, wide_t b) {     //  wrapping mod 2^512
  wide_t r; unsigned __int128 c = 0;
  for ( int i = 0; i < WLIMBS; i++ ) {
    unsigned __int128 s = (unsigned __int128)a.w[i] + b.w[i] + c;
    r.w[i] = (uint64_t)s; c = s >> 64;
  }
  return r;
}
static inline wide_t w_sub(wide_t a, wide_t b) {     //  wrapping mod 2^512
  wide_t r; unsigned __int128 brw = 0;
  for ( int i = 0; i < WLIMBS; i++ ) {
    unsigned __int128 d = (unsigned __int128)a.w[i] - b.w[i] - brw;
    r.w[i] = (uint64_t)d; brw = (d >> 64) & 1;
  }
  return r;
}
static inline wide_t w_shl(wide_t a, int k) {
  wide_t r = w_zero();
  if ( k <= 0 ) return k == 0 ? a : r;
  if ( k >= WBITS ) return r;
  int limb = k / 64, bit = k % 64;
  for ( int i = 0; i < WLIMBS; i++ ) {
    int dst = i + limb;
    if ( dst >= WLIMBS ) break;
    r.w[dst] |= a.w[i] << bit;
    if ( bit && dst + 1 < WLIMBS ) r.w[dst + 1] |= a.w[i] >> (64 - bit);
  }
  return r;
}
static inline wide_t w_shr(wide_t a, int k) {
  wide_t r = w_zero();
  if ( k <= 0 ) return k == 0 ? a : r;
  if ( k >= WBITS ) return r;
  int limb = k / 64, bit = k % 64;
  for ( int i = 0; i < WLIMBS; i++ ) {
    int src = i + limb;
    if ( src >= WLIMBS ) break;
    r.w[i] |= a.w[src] >> bit;
    if ( bit && src + 1 < WLIMBS ) r.w[i] |= a.w[src + 1] << (64 - bit);
  }
  return r;
}
static inline wide_t w_masklow(wide_t a, int k) {    //  keep the low k bits
  wide_t r = w_zero();
  if ( k >= WBITS ) return a;
  if ( k <= 0 ) return r;
  for ( int i = 0; i < WLIMBS; i++ ) {
    int lo = i * 64;
    if ( k <= lo ) break;
    if ( k >= lo + 64 ) r.w[i] = a.w[i];
    else { r.w[i] = a.w[i] & ((1ull << (k - lo)) - 1); break; }
  }
  return r;
}
//  floor(sqrt(x)), digit-by-digit; sets *exact iff the remainder is zero.
static inline wide_t w_isqt(wide_t x, int *exact) {
  wide_t res = w_zero();
  if ( w_is_zero(x) ) { *exact = 1; return res; }
  int p = (w_bits(x) - 1) & ~1;                 //  highest even <= bits-1
  wide_t bit = w_shl(w_from_u64(1), p);         //  largest 4^k <= x
  wide_t num = x;
  while ( !w_is_zero(bit) ) {
    wide_t rb = w_add(res, bit);
    if ( w_cmp(num, rb) >= 0 ) {
      num = w_sub(num, rb);
      res = w_add(w_shr(res, 1), bit);
    } else {
      res = w_shr(res, 1);
    }
    bit = w_shr(bit, 2);
  }
  *exact = w_is_zero(num);
  return res;
}

#endif  //  SOFTUNUM_PWIDE_H
