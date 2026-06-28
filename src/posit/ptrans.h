//  ptrans.h -- elementary / transcendental functions, shared across all widths.
//
//  These are the naive, reproducible series of /lib/unum (and /lib/math):
//  fixed term counts, posit arithmetic throughout, NOT correctly rounded and
//  accurate only near the expansion point (no range reduction).  They are pure
//  compositions of the public posit ops + the rounded constants, so this one
//  body serves every width: it is #included by p8.c, p16.c, p32.c after the
//  core ops and `pconst` are in scope.  Each function transliterates its Hoon
//  arm in /lib/unum line-for-line.
//
//  Requires (before include): PW_PREFIX, PFN, pconst(), MSK, NARBITS, and the
//  public PFN(_add/_sub/_mul/_div/_neg/_sqrt/_from_u64/_abs/_lt/_le/_eq) ops.

#ifndef PFN
#define PW_CAT2(a, b) a ## b
#define PW_CAT(a, b)  PW_CAT2(a, b)
#define PFN(name)     PW_CAT(PW_PREFIX, name)
#endif

#define ADD PFN(_add)
#define SUB PFN(_sub)
#define MUL PFN(_mul)
#define DIV PFN(_div)
#define NEG PFN(_neg)
#define SQT PFN(_sqrt)
#define SUN PFN(_from_u64)
#define ABS PFN(_abs)
#define LT  PFN(_lt)
#define LE  PFN(_le)
#define EQ  PFN(_eq)

//  Constants: Q2.52 fixed-point hex of the value, rounded by +bit at each
//  width (sign +, exponent -52), exactly as /lib/unum bakes them in.
uint32_t PFN(_one)(void)      { return pconst(0,   1); }
uint32_t PFN(_pi)(void)       { return pconst(-52, 0x3243f6a8885a31ull); }
uint32_t PFN(_tau)(void)      { return pconst(-52, 0x6487ed5110b461ull); }
uint32_t PFN(_e)(void)        { return pconst(-52, 0x2b7e151628aed3ull); }
uint32_t PFN(_phi)(void)      { return pconst(-52, 0x19e3779b97f4a8ull); }
uint32_t PFN(_sqrt2)(void)    { return pconst(-52, 0x16a09e667f3bcdull); }
uint32_t PFN(_invsqrt2)(void) { return pconst(-52, 0x0b504f333f9de6ull); }
uint32_t PFN(_ln2)(void)      { return pconst(-52, 0x0b17217f7d1cf8ull); }
uint32_t PFN(_invln2)(void)   { return pconst(-52, 0x171547652b82feull); }
uint32_t PFN(_ln10)(void)     { return pconst(-52, 0x24d763776aaa2bull); }

//  exp(x) = sum_{k>=0} x^k / k!   (20 terms)
uint32_t PFN(_exp)(uint32_t x) {
  uint32_t one = PFN(_one)(), sum = one, term = one;
  for ( uint64_t nn = 1; nn <= 20; nn++ ) {
    term = MUL(term, DIV(x, SUN(nn)));
    sum = ADD(sum, term);
  }
  return sum;
}

//  sin(x) = sum (-1)^k x^(2k+1) / (2k+1)!
uint32_t PFN(_sin)(uint32_t x) {
  uint32_t term = x, sum = x, x2 = MUL(x, x);
  for ( uint64_t nn = 1; nn <= 20; nn++ ) {
    uint64_t k = 2 * nn;
    term = NEG(MUL(term, DIV(x2, MUL(SUN(k), SUN(k + 1)))));
    sum = ADD(sum, term);
  }
  return sum;
}

//  cos(x) = sum (-1)^k x^2k / (2k)!
uint32_t PFN(_cos)(uint32_t x) {
  uint32_t one = PFN(_one)(), term = one, sum = one, x2 = MUL(x, x);
  for ( uint64_t nn = 1; nn <= 20; nn++ ) {
    uint64_t k = 2 * nn;
    term = NEG(MUL(term, DIV(x2, MUL(SUN(k - 1), SUN(k)))));
    sum = ADD(sum, term);
  }
  return sum;
}

uint32_t PFN(_tan)(uint32_t x) { return DIV(PFN(_sin)(x), PFN(_cos)(x)); }

//  pow-n: integer power by repeated multiplication (NaR propagates even at p=0).
uint32_t PFN(_pow_n)(uint32_t x, uint64_t p) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  uint32_t res = PFN(_one)();
  while ( p != 0 ) { res = MUL(res, x); p--; }
  return res;
}

//  log(x) = 2 * atanh((x-1)/(x+1))   (30 terms); NaR for x <= 0.
uint32_t PFN(_log)(uint32_t x) {
  if ( LE(x, 0) ) return NARBITS;
  uint32_t one = PFN(_one)();
  uint32_t y = DIV(SUB(x, one), ADD(x, one));
  uint32_t y2 = MUL(y, y), sum = y, term = y;
  for ( uint64_t nn = 1; nn <= 30; nn++ ) {
    term = MUL(term, y2);
    uint32_t coef = DIV(one, SUN(2 * nn + 1));
    sum = ADD(sum, MUL(coef, term));
  }
  return MUL(SUN(2), sum);
}

uint32_t PFN(_log2)(uint32_t x)  { return DIV(PFN(_log)(x), PFN(_ln2)()); }
uint32_t PFN(_log10)(uint32_t x) { return DIV(PFN(_log)(x), PFN(_ln10)()); }

uint32_t PFN(_pow)(uint32_t x, uint32_t y) { return PFN(_exp)(MUL(y, PFN(_log)(x))); }

//  factorial: x! by repeated multiplication; NaR for x < 0 (and NaR-propagates).
uint32_t PFN(_factorial)(uint32_t x) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  if ( LT(x, 0) ) return NARBITS;
  uint32_t one = PFN(_one)(), t = one;
  while ( !LE(x, one) ) { t = MUL(t, x); x = SUB(x, one); }
  return t;
}

//  cbrt(x) = exp(log(x)/3); domain x > 0, cbrt(0) = 0.
uint32_t PFN(_cbrt)(uint32_t x) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  if ( (x & MSK) == 0 ) return 0;
  if ( LT(x, 0) ) return NARBITS;
  uint32_t one = PFN(_one)();
  return PFN(_pow)(x, DIV(one, SUN(3)));
}

//  atan via Gauss/AGM iteration (41 steps).
uint32_t PFN(_atan)(uint32_t x) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  uint32_t one = PFN(_one)();
  uint32_t rt = SQT(ADD(one, MUL(x, x)));
  uint32_t a = DIV(one, rt), b = one;
  for ( uint64_t nn = 0; nn <= 40; nn++ ) {
    uint32_t ai = MUL(DIV(one, SUN(2)), ADD(a, b));
    uint32_t bi = SQT(MUL(ai, b));
    a = ai; b = bi;
  }
  return DIV(x, MUL(rt, b));
}

//  asin(x) = atan(x / sqrt(1-x^2)) for |x|<1; +-pi/2 at +-1; NaR outside [-1,1].
uint32_t PFN(_asin)(uint32_t x) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  uint32_t one = PFN(_one)();
  if ( LT(ABS(x), one) )
    return PFN(_atan)(DIV(x, SQT(SUB(one, MUL(x, x)))));
  if ( EQ(x, one) )      return MUL(PFN(_pi)(), DIV(one, SUN(2)));
  if ( EQ(x, NEG(one)) ) return NEG(MUL(PFN(_pi)(), DIV(one, SUN(2))));
  return NARBITS;
}

//  acos(x) = atan(sqrt(1-x^2)/x) for 0<|x|<1 (pi/2 at 0); 0 at 1, pi at -1.
uint32_t PFN(_acos)(uint32_t x) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  uint32_t one = PFN(_one)();
  if ( LT(ABS(x), one) ) {
    if ( EQ(x, 0) ) return MUL(PFN(_pi)(), DIV(one, SUN(2)));
    return PFN(_atan)(DIV(SQT(SUB(one, MUL(x, x))), x));
  }
  if ( EQ(x, one) )      return 0;
  if ( EQ(x, NEG(one)) ) return PFN(_pi)();
  return NARBITS;
}

int PFN(_is_close)(uint32_t a, uint32_t b, uint32_t tol) {
  return LE(ABS(SUB(a, b)), tol);
}

#undef ADD
#undef SUB
#undef MUL
#undef DIV
#undef NEG
#undef SQT
#undef SUN
#undef ABS
#undef LT
#undef LE
#undef EQ
