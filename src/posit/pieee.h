//  pieee.h -- IEEE-754 binary{16,32,64,128} codec for posit <-> float
//  conversion, shared by all posit widths.  Still no SoftFloat: a float is just
//  a bit pattern, and we encode/decode it with integer arithmetic.
//
//  This mirrors the Hoon stdlib float g-layer +$fn and its ++ff sea/bit:
//    +$fn = [%f s=? e=@s a=@u] (finite, value = +-a*2^e, a includes hidden bit)
//         | [%i s=?] (infinity) | [%n ~] (NaN)
//  decode (f_sea) is exact; encode (f_bit) rounds to nearest, ties to even --
//  the door's default mode (%n) -- which is the UNIQUE correctly-rounded result,
//  hence bit-identical to the Hoon ++fl rounding without transliterating it.
//
//  Per-width parameters (w=exponent bits, p=mantissa bits, bias):
//    binary16  w=5  p=10  bias=15      binary32  w=8  p=23  bias=127
//    binary64  w=11 p=52  bias=1023    binary128 w=15 p=112 bias=16383

#ifndef SOFTUNUM_PIEEE_H
#define SOFTUNUM_PIEEE_H

#include <stdint.h>

typedef unsigned __int128 u128i;

enum { FN_FIN = 0, FN_INF = 1, FN_NAN = 2 };
typedef struct { int kind; int sign; long long e; u128i a; } fnt;   // a = significand+hidden
typedef struct { int w; int p; long long bias; } ifmt;

static const ifmt IEEE_RH = { 5, 10, 15 };
static const ifmt IEEE_RS = { 8, 23, 127 };
static const ifmt IEEE_RD = { 11, 52, 1023 };
static const ifmt IEEE_RQ = { 15, 112, 16383 };

static inline int u128_clz(u128i x) {            //  x != 0
  uint64_t hi = (uint64_t)(x >> 64);
  return hi ? __builtin_clzll(hi) : 64 + __builtin_clzll((uint64_t)x);
}
static inline int u128_len(u128i x) { return x ? 128 - u128_clz(x) : 0; }

//  round(a * 2^-s) to nearest, ties to even (s<0 is an exact left shift).
static inline u128i rne_rshift(u128i a, int s) {
  if ( s <= 0 ) return a << (-s);
  if ( s >= 128 ) {
    if ( s > 128 ) return 0;
    u128i half = (u128i)1 << 127;
    return (a > half) ? 1 : 0;                  //  tie -> even (0)
  }
  u128i keep = a >> s;
  u128i rem = a & (((u128i)1 << s) - 1);
  u128i half = (u128i)1 << (s - 1);
  if ( rem > half || (rem == half && (keep & 1)) ) keep++;
  return keep;
}

//  decode an IEEE bit pattern (exact), mirroring ++ff sea.
static fnt f_sea(u128i bits, ifmt F) {
  int w = F.w, p = F.p; long long b = F.bias;
  u128i f = bits & (((u128i)1 << p) - 1);
  u128i e = (bits >> p) & (((u128i)1 << w) - 1);
  long long me = (1 - b) - p;
  fnt r = { FN_FIN, ((bits >> (p + w)) & 1) ? 0 : 1, 0, 0 };  //  sign 1 = positive
  if ( e == 0 ) {                                //  zero / subnormal
    if ( f != 0 ) { r.e = me; r.a = f; }
    return r;
  }
  if ( e == (((u128i)1 << w) - 1) ) {            //  inf / NaN
    r.kind = (f == 0) ? FN_INF : FN_NAN;
    return r;
  }
  r.e = (long long)e + me - 1;                   //  normal
  r.a = f + ((u128i)1 << p);
  return r;
}

//  encode to an IEEE bit pattern, round-nearest-ties-even, mirroring ++ff bit.
static u128i f_bit(fnt x, ifmt F) {
  int w = F.w, p = F.p; long long b = F.bias;
  u128i sb = (u128i)1 << (w + p);
  u128i allexp = ((u128i)1 << w) - 1;
  u128i inf = allexp << p;
  if ( x.kind == FN_NAN ) return ((u128i)((1u << (w + 1)) - 1)) << (p - 1);
  if ( x.kind == FN_INF ) return x.sign ? inf : (inf + sb);
  if ( x.a == 0 )         return x.sign ? 0 : sb;
  int L = u128_len(x.a) - 1;
  long long Ebias = x.e + (long long)L + b;
  long long efield; u128i mant;
  if ( Ebias >= (long long)allexp ) return x.sign ? inf : (inf + sb);   //  overflow
  if ( Ebias >= 1 ) {                            //  normal
    u128i sh = rne_rshift(x.a, L - p);
    if ( sh >= ((u128i)1 << (p + 1)) ) { sh >>= 1; Ebias++; }   //  rounding carry
    if ( Ebias >= (long long)allexp ) return x.sign ? inf : (inf + sb);
    efield = Ebias; mant = sh & (((u128i)1 << p) - 1);
  } else {                                       //  subnormal / underflow
    long long k = x.e - 1 + b + (long long)p;    //  m = round(a * 2^k)
    u128i m = rne_rshift(x.a, (int)(-k));
    if ( m == 0 ) return x.sign ? 0 : sb;
    if ( m >= ((u128i)1 << p) ) { efield = 1; mant = m & (((u128i)1 << p) - 1); }
    else { efield = 0; mant = m; }
  }
  u128i bits = ((u128i)efield << p) | mant;
  return x.sign ? bits : (bits + sb);
}

#endif  //  SOFTUNUM_PIEEE_H
