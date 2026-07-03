//  ptrans.h -- elementary / transcendental functions, shared across all widths.
//
//  These are BIT-EXACT range-reduced Chebyshev-minimax / exact-Taylor kernels,
//  the C twin of the rewritten unum.hoon transcendentals (numerics repo,
//  NEXT-STEPS.md item #4): correctly rounded (0 ULP vs mpmath) for exp/log/
//  log2/log10/sin/cos/atan at posit8/16/32, faithful (a few ULP) for
//  tan/asin/acos.  Each function transliterates its Hoon arm line-for-line;
//  the comments cite the Hoon arm names.  Decode a posit ONCE (+sea, i.e. the
//  existing per-width `sea()`), do EXACT (never-rounded) arithmetic on the
//  [sign exponent significand] triple via the shared GMP g-layer (pgmp.h --
//  +gmul/+gadd/+gneg/+gsub/+gdiv/+gpoly/+g-round/+glt/+bit), and round only
//  ONCE at the very end via g_bit.  This one body serves every width: it is
//  #included by p8.c, p16.c, p32.c after the core ops, `pconst()`, and
//  UP_TO_GVAL() are in scope.
//
//  Requires (before include): PW_PREFIX, PFN, pconst(), MSK, NARBITS, N, the
//  public PFN(_add/_sub/_mul/_div/_neg/_sqrt/_from_u64/_abs/_lt/_le/_eq) ops,
//  and UP_TO_GVAL(up_t) -- see p8.c / pcore.h for the two definitions (u128
//  vs wide_t significand).

#ifndef PFN
#define PW_CAT2(a, b) a ## b
#define PW_CAT(a, b)  PW_CAT2(a, b)
#define PFN(name)     PW_CAT(PW_PREFIX, name)
#endif

#include "pgmp.h"

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

//  ---- shared wide fixed-point constants + lazy one-time init ---------------
//  Hex significands transcribed VERBATIM from unum.hoon (dots stripped only;
//  no re-derivation -- cross-checked against libmath/tools/unum_cheb_check.py's
//  independently-computed hex dumps).  `static` (internal linkage): this
//  header is textually included once per translation unit (p8.c, and once
//  each for p16/p32 via pcore.h), so each TU gets its own private copy --
//  exactly like the existing per-width statics (sea/bit/pconst).

static int          PFN(_gci_ready) = 0;
static gval_t        PFN(_ONE_G);
static gval_t        PFN(_TWO_G);
static gval_t        PFN(_HALF_G);
static gval_t        PFN(_THREEHALF_G);
static gval_t        PFN(_LN2_WIDE);
static gval_t        PFN(_INVLN2_WIDE);
static gval_t        PFN(_LOG10_TWO_WIDE);
static gval_t        PFN(_INVLN10_WIDE);
static gval_t        PFN(_PI2_WIDE);
static gval_t        PFN(_INVPI2_WIDE);
static gval_t        PFN(_EXP_CS)[8];
static gval_t        PFN(_SIN_CS)[6];
static gval_t        PFN(_COS_CS)[6];
static gval_t        PFN(_LR_CS)[16];
static gval_t        PFN(_ATAN_CS)[13];
static gval_t        PFN(_ATAN_B1); static gval_t PFN(_ATAN_B2);
static gval_t        PFN(_ATAN_B3); static gval_t PFN(_ATAN_B4);
static gval_t        PFN(_ATAN_BP0); static gval_t PFN(_ATAN_BPHALF);
static gval_t        PFN(_ATAN_BPPI4); static gval_t PFN(_ATAN_BPTHREEHALF);
static gval_t        PFN(_ATAN_BPPI2);

static void PFN(_gconst_init)(void) {
  if ( PFN(_gci_ready) ) return;
  PFN(_gci_ready) = 1;
  PFN(_ONE_G)       = gval_small(1, 0, 1);
  PFN(_TWO_G)       = gval_small(1, 1, 1);
  PFN(_HALF_G)      = gval_small(1, -1, 1);
  PFN(_THREEHALF_G) = gval_small(1, -1, 3);
  //  +ln2-wide / +invln2-wide / +log10-two-wide / +invln10-wide (unum.hoon)
  PFN(_LN2_WIDE)        = gval_const(1, -128, "b17217f7d1cf79abc9e3b39803f2f6af");
  PFN(_INVLN2_WIDE)     = gval_const(1, -128, "171547652b82fe1777d0ffda0d23a7d12");
  PFN(_LOG10_TWO_WIDE)  = gval_const(1, -128, "4d104d427de7fbcc47c4acd605be48bc");
  PFN(_INVLN10_WIDE)    = gval_const(1, -128, "6f2dec549b9438ca9aadd557d699ee19");
  //  +pi2-wide / +invpi2-wide (TRIG_WBITS=200 -- see unum.hoon's comment: trig
  //  reduction error scales with the argument, unlike exp/log's flat margin).
  PFN(_PI2_WIDE)    = gval_const(1, -200, "1921fb54442d18469898cc51701b839a252049c1114cf98e804");
  PFN(_INVPI2_WIDE) = gval_const(1, -200, "a2f9836e4e441529fc2757d1f534ddc0db6295993c439041fe");
  //  +exp's `cs`: c7..c0, highest degree first (degree-7 Chebyshev-minimax).
  PFN(_EXP_CS)[0] = gval_const(1, -128, "d0bef97615e3036b3f095cd13a64d");
  PFN(_EXP_CS)[1] = gval_const(1, -128, "5b69d4c616b405aa1ffc644cd9fc9e");
  PFN(_EXP_CS)[2] = gval_const(1, -128, "222214c445d4691bf8c9a8af541fc45");
  PFN(_EXP_CS)[3] = gval_const(1, -128, "aaaa3250783cf22e306b08bd797037e");
  PFN(_EXP_CS)[4] = gval_const(1, -128, "2aaaaaafce24bb583aef94f8c614beae");
  PFN(_EXP_CS)[5] = gval_const(1, -128, "8000002e44f26dbc722719432ee4604c");
  PFN(_EXP_CS)[6] = gval_const(1, -128, "fffffffffb0fe8442a05c50f9cfaea45");
  PFN(_EXP_CS)[7] = gval_const(1, -128, "ffffffffd389a50d1cdb82b4c34e2f1d");
  //  +sc's `sin-cs` / `cos-cs`: k=6..1, highest degree first.
  PFN(_SIN_CS)[0] = gval_const(1, -128, "b092309d43684be51c198e92");
  PFN(_SIN_CS)[1] = gval_const(0, -128, "6b99159fd5138e3f9d1f92e0df");
  PFN(_SIN_CS)[2] = gval_const(1, -128, "2e3bc74aad8e671f5583911ca003");
  PFN(_SIN_CS)[3] = gval_const(0, -128, "d00d00d00d00d00d00d00d00d00d0");
  PFN(_SIN_CS)[4] = gval_const(1, -128, "2222222222222222222222222222222");
  PFN(_SIN_CS)[5] = gval_const(0, -128, "2aaaaaaaaaaaaaaaaaaaaaaaaaaaaaab");
  PFN(_COS_CS)[0] = gval_const(1, -128, "8f76c77fc6c4bdaa26d4c3d68");
  PFN(_COS_CS)[1] = gval_const(0, -128, "49f93edde27d71cbbc05b4fa99a");
  PFN(_COS_CS)[2] = gval_const(1, -128, "1a01a01a01a01a01a01a01a01a01a");
  PFN(_COS_CS)[3] = gval_const(0, -128, "5b05b05b05b05b05b05b05b05b05b0");
  PFN(_COS_CS)[4] = gval_const(1, -128, "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaab");
  PFN(_COS_CS)[5] = gval_const(0, -128, "80000000000000000000000000000000");
  //  +lr's `cs`: c_16..c_1 = 1/33..1/3, highest degree first.
  PFN(_LR_CS)[0]  = gval_const(1, -128, "7c1f07c1f07c1f07c1f07c1f07c1f08");
  PFN(_LR_CS)[1]  = gval_const(1, -128, "8421084210842108421084210842108");
  PFN(_LR_CS)[2]  = gval_const(1, -128, "8d3dcb08d3dcb08d3dcb08d3dcb08d4");
  PFN(_LR_CS)[3]  = gval_const(1, -128, "97b425ed097b425ed097b425ed097b4");
  PFN(_LR_CS)[4]  = gval_const(1, -128, "a3d70a3d70a3d70a3d70a3d70a3d70a");
  PFN(_LR_CS)[5]  = gval_const(1, -128, "b21642c8590b21642c8590b21642c86");
  PFN(_LR_CS)[6]  = gval_const(1, -128, "c30c30c30c30c30c30c30c30c30c30c");
  PFN(_LR_CS)[7]  = gval_const(1, -128, "d79435e50d79435e50d79435e50d794");
  PFN(_LR_CS)[8]  = gval_const(1, -128, "f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f");
  PFN(_LR_CS)[9]  = gval_const(1, -128, "11111111111111111111111111111111");
  PFN(_LR_CS)[10] = gval_const(1, -128, "13b13b13b13b13b13b13b13b13b13b14");
  PFN(_LR_CS)[11] = gval_const(1, -128, "1745d1745d1745d1745d1745d1745d17");
  PFN(_LR_CS)[12] = gval_const(1, -128, "1c71c71c71c71c71c71c71c71c71c71c");
  PFN(_LR_CS)[13] = gval_const(1, -128, "24924924924924924924924924924925");
  PFN(_LR_CS)[14] = gval_const(1, -128, "33333333333333333333333333333333");
  PFN(_LR_CS)[15] = gval_const(1, -128, "55555555555555555555555555555555");
  //  +atan-core's fdlibm breakpoint thresholds (exact dyadic rationals) and
  //  angle constants, plus its `cs`: k=13..1, highest degree first.
  PFN(_ATAN_B1) = gval_small(1, -4, 7);
  PFN(_ATAN_B2) = gval_small(1, -4, 11);
  PFN(_ATAN_B3) = gval_small(1, -4, 19);
  PFN(_ATAN_B4) = gval_small(1, -4, 39);
  PFN(_ATAN_BP0)         = gval_small(1, 0, 0);
  PFN(_ATAN_BPHALF)      = gval_const(1, -128, "76b19c1586ed3da2b7f222f65e1d4682");
  PFN(_ATAN_BPPI4)       = gval_const(1, -128, "c90fdaa22168c234c4c6628b80dc1cd1");
  PFN(_ATAN_BPTHREEHALF) = gval_const(1, -128, "fb985e940fb4d9007887af0cbbc9e142");
  PFN(_ATAN_BPPI2)       = gval_const(1, -128, "1921fb54442d18469898cc51701b839a2");
  PFN(_ATAN_CS)[0]  = gval_const(0, -128, "97b425ed097b425ed097b425ed097b4");
  PFN(_ATAN_CS)[1]  = gval_const(1, -128, "a3d70a3d70a3d70a3d70a3d70a3d70a");
  PFN(_ATAN_CS)[2]  = gval_const(0, -128, "b21642c8590b21642c8590b21642c86");
  PFN(_ATAN_CS)[3]  = gval_const(1, -128, "c30c30c30c30c30c30c30c30c30c30c");
  PFN(_ATAN_CS)[4]  = gval_const(0, -128, "d79435e50d79435e50d79435e50d794");
  PFN(_ATAN_CS)[5]  = gval_const(1, -128, "f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f");
  PFN(_ATAN_CS)[6]  = gval_const(0, -128, "11111111111111111111111111111111");
  PFN(_ATAN_CS)[7]  = gval_const(1, -128, "13b13b13b13b13b13b13b13b13b13b14");
  PFN(_ATAN_CS)[8]  = gval_const(0, -128, "1745d1745d1745d1745d1745d1745d17");
  PFN(_ATAN_CS)[9]  = gval_const(1, -128, "1c71c71c71c71c71c71c71c71c71c71c");
  PFN(_ATAN_CS)[10] = gval_const(0, -128, "24924924924924924924924924924925");
  PFN(_ATAN_CS)[11] = gval_const(1, -128, "33333333333333333333333333333333");
  PFN(_ATAN_CS)[12] = gval_const(0, -128, "55555555555555555555555555555555");
}

//  ---- +exp:  e^x = 2^k * poly(r), x = k*ln2 + r (Chebyshev-minimax degree-7,
//  correctly rounded at posit8/16/32).  +bit's own maxpos/minpos saturation
//  gives overflow/underflow for free.
uint32_t PFN(_exp)(uint32_t x) {
  up_t u = sea(x);
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.kind == K_ZERO ) return PFN(_one)();
  PFN(_gconst_init)();
  gval_t g = UP_TO_GVAL(u);
  gval_t gk = gmul(g, PFN(_INVLN2_WIDE));
  mpz_t kmpz; mpz_init(kmpz);
  g_round(kmpz, gk);
  gval_clear(&gk);
  gval_t kup = gval_from_signed_mpz(kmpz);
  gval_t kln2 = gmul(kup, PFN(_LN2_WIDE));
  gval_clear(&kup);
  gval_t r = gsub(g, kln2);
  gval_clear(&g); gval_clear(&kln2);
  gval_t p = gpoly(PFN(_EXP_CS), 8, r);
  gval_clear(&r);
  p.e = ll_add_clamp(p.e, kmpz, 1LL << 40);
  mpz_clear(kmpz);
  uint32_t result = g_bit(p, N);
  gval_clear(&p);
  return result;
}

//  ---- +sc:  shared quarter-turn reduction for sin/cos/tan.  ax must be a
//  non-negative gval_t (sign forced by the caller); returns the UNROUNDED
//  [sin(ax), cos(ax)] pair, still exact g-triples.
typedef struct { gval_t sn; gval_t cs; } PFN(_sc_t);

static PFN(_sc_t) PFN(_sc)(gval_t ax) {
  PFN(_gconst_init)();
  gval_t qg = gmul(ax, PFN(_INVPI2_WIDE));
  mpz_t qmpz; mpz_init(qmpz);
  g_round(qmpz, qg);
  gval_clear(&qg);
  mpz_abs(qmpz, qmpz);                  //  qn = abs(q); always >=0 here anyway
  gval_t qn_g = gval_from_signed_mpz(qmpz);
  gval_t qnpi2 = gmul(qn_g, PFN(_PI2_WIDE));
  gval_clear(&qn_g);
  gval_t r = gsub(ax, qnpi2);
  gval_clear(&qnpi2);
  gval_t z = gmul(r, r);
  gval_t sinpoly = gpoly(PFN(_SIN_CS), 6, z);
  gval_t sin_zp = gmul(z, sinpoly); gval_clear(&sinpoly);
  gval_t sin_1p = gadd(PFN(_ONE_G), sin_zp); gval_clear(&sin_zp);
  gval_t sink = gmul(r, sin_1p); gval_clear(&sin_1p); gval_clear(&r);
  gval_t cospoly = gpoly(PFN(_COS_CS), 6, z);
  gval_t cos_zp = gmul(z, cospoly); gval_clear(&cospoly);
  gval_t cosk = gadd(PFN(_ONE_G), cos_zp); gval_clear(&cos_zp);
  gval_clear(&z);
  unsigned long m = mpz_fdiv_ui(qmpz, 4);
  mpz_clear(qmpz);
  PFN(_sc_t) out;
  if ( m == 0 )      { out.sn = sink; out.cs = cosk; }
  else if ( m == 1 ) { out.sn = cosk; out.cs = gneg(sink); gval_clear(&sink); }
  else if ( m == 2 ) { out.sn = gneg(sink); out.cs = gneg(cosk); gval_clear(&sink); gval_clear(&cosk); }
  else               { out.sn = gneg(cosk); out.cs = sink; gval_clear(&cosk); }
  return out;
}

//  ---- +sin:  odd fn, reduce |x| via +sc, reapply x's sign.
uint32_t PFN(_sin)(uint32_t x) {
  up_t u = sea(x);
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.kind == K_ZERO ) return 0;
  PFN(_gconst_init)();
  gval_t ax = UP_TO_GVAL(u); ax.sign = 1;
  PFN(_sc_t) sc = PFN(_sc)(ax);
  gval_clear(&ax); gval_clear(&sc.cs);
  sc.sn.sign = u.sign ? sc.sn.sign : !sc.sn.sign;
  uint32_t result = g_bit(sc.sn, N);
  gval_clear(&sc.sn);
  return result;
}

//  ---- +cos:  even fn, reduce |x| via +sc, sign unaffected by x's sign.
uint32_t PFN(_cos)(uint32_t x) {
  up_t u = sea(x);
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.kind == K_ZERO ) return PFN(_one)();
  PFN(_gconst_init)();
  gval_t ax = UP_TO_GVAL(u); ax.sign = 1;
  PFN(_sc_t) sc = PFN(_sc)(ax);
  gval_clear(&ax); gval_clear(&sc.sn);
  uint32_t result = g_bit(sc.cs, N);
  gval_clear(&sc.cs);
  return result;
}

//  ---- +tan:  sin(ax)/cos(ax) via +gdiv on the UNROUNDED kernel outputs --
//  one final rounding; odd fn, same sign handling as +sin.
uint32_t PFN(_tan)(uint32_t x) {
  up_t u = sea(x);
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.kind == K_ZERO ) return 0;
  PFN(_gconst_init)();
  gval_t ax = UP_TO_GVAL(u); ax.sign = 1;
  PFN(_sc_t) sc = PFN(_sc)(ax);
  gval_clear(&ax);
  if ( mpz_sgn(sc.cs.a) == 0 ) {         //  cos(ax) exactly 0 (measure-zero, guard anyway)
    gval_clear(&sc.sn); gval_clear(&sc.cs);
    return NARBITS;
  }
  gval_t raw = gdiv(sc.sn, sc.cs, 160);
  gval_clear(&sc.sn); gval_clear(&sc.cs);
  raw.sign = u.sign ? raw.sign : !raw.sign;
  uint32_t result = g_bit(raw, N);
  gval_clear(&raw);
  return result;
}

//  pow-n: integer power by repeated multiplication (unchanged -- exact
//  compose of already-rounded posit ops, no exp/log involved).
uint32_t PFN(_pow_n)(uint32_t x, uint64_t p) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  uint32_t res = PFN(_one)();
  while ( p != 0 ) { res = MUL(res, x); p--; }
  return res;
}

//  ---- +lr:  shared mantissa/exponent split for +log/+log-2/+log-10.
//  x = m*2^E (free from +sea's own split), m in [1,2); log(m) = 2*atanh(s),
//  s = f/(m+1), f = m-1, via the EXACT atanh Taylor series (z=s*s converges
//  fast; degree-16 is the smallest correctly-rounded at p8/16/32).
typedef struct { long long e; gval_t m; } PFN(_lr_t);

static PFN(_lr_t) PFN(_lr)(gval_t g) {
  PFN(_gconst_init)();
  unsigned long lead = (unsigned long)(mpz_sizeinbase(g.a, 2) - 1);
  gval_t mup = gval_new(); mup.sign = 1; mup.e = -(long long)lead; mpz_set(mup.a, g.a);
  long long bige = (long long)lead + g.e;
  gval_t f = gsub(mup, PFN(_ONE_G));
  gval_t mp1 = gadd(mup, PFN(_ONE_G));
  gval_t s = gdiv(f, mp1, 160);
  gval_clear(&f); gval_clear(&mp1); gval_clear(&mup);
  gval_t z = gmul(s, s);
  gval_t poly = gpoly(PFN(_LR_CS), 16, z);
  gval_t zpoly = gmul(z, poly); gval_clear(&poly); gval_clear(&z);
  gval_t onep = gadd(PFN(_ONE_G), zpoly); gval_clear(&zpoly);
  gval_t twos = gmul(PFN(_TWO_G), s); gval_clear(&s);
  gval_t logm = gmul(twos, onep); gval_clear(&twos); gval_clear(&onep);
  PFN(_lr_t) out; out.e = bige; out.m = logm;
  return out;
}

//  ---- +log:  ln x = E*ln2 + log(m), x = m*2^E via +lr.  Domain x > 0.
uint32_t PFN(_log)(uint32_t x) {
  up_t u = sea(x);
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.kind == K_ZERO ) return NARBITS;
  if ( !u.sign )          return NARBITS;
  PFN(_gconst_init)();
  gval_t g = UP_TO_GVAL(u);
  PFN(_lr_t) em = PFN(_lr)(g);
  gval_clear(&g);
  gval_t eup = gval_small(em.e >= 0, 0, em.e >= 0 ? (unsigned long)em.e : (unsigned long)(-em.e));
  gval_t eln2 = gmul(eup, PFN(_LN2_WIDE)); gval_clear(&eup);
  gval_t sum = gadd(eln2, em.m); gval_clear(&eln2); gval_clear(&em.m);
  uint32_t result = g_bit(sum, N);
  gval_clear(&sum);
  return result;
}

//  ---- +log-2 / +log-10:  base-2 / base-10 log, via +lr directly (E +
//  log(m)/ln(b)) -- avoids a second rounding through a posit-rounded
//  log2/log10 constant.
uint32_t PFN(_log2)(uint32_t x) {
  up_t u = sea(x);
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.kind == K_ZERO ) return NARBITS;
  if ( !u.sign )          return NARBITS;
  PFN(_gconst_init)();
  gval_t g = UP_TO_GVAL(u);
  PFN(_lr_t) em = PFN(_lr)(g);
  gval_clear(&g);
  gval_t eup = gval_small(em.e >= 0, 0, em.e >= 0 ? (unsigned long)em.e : (unsigned long)(-em.e));
  gval_t mlog = gmul(em.m, PFN(_INVLN2_WIDE)); gval_clear(&em.m);
  gval_t sum = gadd(eup, mlog); gval_clear(&eup); gval_clear(&mlog);
  uint32_t result = g_bit(sum, N);
  gval_clear(&sum);
  return result;
}

uint32_t PFN(_log10)(uint32_t x) {
  up_t u = sea(x);
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.kind == K_ZERO ) return NARBITS;
  if ( !u.sign )          return NARBITS;
  PFN(_gconst_init)();
  gval_t g = UP_TO_GVAL(u);
  PFN(_lr_t) em = PFN(_lr)(g);
  gval_clear(&g);
  gval_t eup = gval_small(em.e >= 0, 0, em.e >= 0 ? (unsigned long)em.e : (unsigned long)(-em.e));
  gval_t t1 = gmul(eup, PFN(_LOG10_TWO_WIDE)); gval_clear(&eup);
  gval_t t2 = gmul(em.m, PFN(_INVLN10_WIDE)); gval_clear(&em.m);
  gval_t sum = gadd(t1, t2); gval_clear(&t1); gval_clear(&t2);
  uint32_t result = g_bit(sum, N);
  gval_clear(&sum);
  return result;
}

uint32_t PFN(_pow)(uint32_t x, uint32_t y) { return PFN(_exp)(MUL(y, PFN(_log)(x))); }

//  factorial: x! by repeated multiplication; NaR for x < 0 (unchanged).
uint32_t PFN(_factorial)(uint32_t x) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  if ( LT(x, 0) ) return NARBITS;
  uint32_t one = PFN(_one)(), t = one;
  while ( !LE(x, one) ) { t = MUL(t, x); x = SUB(x, one); }
  return t;
}

//  cbrt(x) = exp(log(x)/3); domain x > 0, cbrt(0) = 0 (unchanged).
uint32_t PFN(_cbrt)(uint32_t x) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  if ( (x & MSK) == 0 ) return 0;
  if ( LT(x, 0) ) return NARBITS;
  uint32_t one = PFN(_one)();
  return PFN(_pow)(x, DIV(one, SUN(3)));
}

//  ---- +atan-core:  fdlibm breakpoint reduction (7/16, 11/16, 19/16, 39/16 --
//  exact dyadic thresholds) + degree-13 EXACT Taylor kernel in z=xr*xr.
//  ax must already be a non-negative gval_t.
static gval_t PFN(_atan_core)(gval_t ax) {
  PFN(_gconst_init)();
  gval_t xr, bp;
  if ( glt(ax, PFN(_ATAN_B1)) ) {
    xr = gval_copy(ax);
    bp = gval_copy(PFN(_ATAN_BP0));
  } else if ( glt(ax, PFN(_ATAN_B2)) ) {
    gval_t num = gsub(ax, PFN(_HALF_G));
    gval_t axh = gmul(ax, PFN(_HALF_G));
    gval_t den = gadd(PFN(_ONE_G), axh); gval_clear(&axh);
    xr = gdiv(num, den, 160); gval_clear(&num); gval_clear(&den);
    bp = gval_copy(PFN(_ATAN_BPHALF));
  } else if ( glt(ax, PFN(_ATAN_B3)) ) {
    gval_t num = gsub(ax, PFN(_ONE_G));
    gval_t den = gadd(ax, PFN(_ONE_G));
    xr = gdiv(num, den, 160); gval_clear(&num); gval_clear(&den);
    bp = gval_copy(PFN(_ATAN_BPPI4));
  } else if ( glt(ax, PFN(_ATAN_B4)) ) {
    gval_t num = gsub(ax, PFN(_THREEHALF_G));
    gval_t axt = gmul(ax, PFN(_THREEHALF_G));
    gval_t den = gadd(PFN(_ONE_G), axt); gval_clear(&axt);
    xr = gdiv(num, den, 160); gval_clear(&num); gval_clear(&den);
    bp = gval_copy(PFN(_ATAN_BPTHREEHALF));
  } else {
    gval_t inv = gdiv(PFN(_ONE_G), ax, 160);
    xr = gneg(inv); gval_clear(&inv);
    bp = gval_copy(PFN(_ATAN_BPPI2));
  }
  gval_t z = gmul(xr, xr);
  gval_t poly = gpoly(PFN(_ATAN_CS), 13, z);
  gval_t zpoly = gmul(z, poly); gval_clear(&poly); gval_clear(&z);
  gval_t series = gadd(PFN(_ONE_G), zpoly); gval_clear(&zpoly);
  gval_t at_xr = gmul(xr, series); gval_clear(&xr); gval_clear(&series);
  gval_t result = gadd(bp, at_xr); gval_clear(&bp); gval_clear(&at_xr);
  return result;
}

//  ---- +atan:  odd fn, reduce |x| via +atan-core, reapply x's sign.
uint32_t PFN(_atan)(uint32_t x) {
  up_t u = sea(x);
  if ( u.kind == K_NAR )  return NARBITS;
  if ( u.kind == K_ZERO ) return 0;
  PFN(_gconst_init)();
  gval_t ax = UP_TO_GVAL(u); ax.sign = 1;
  gval_t raw = PFN(_atan_core)(ax);
  gval_clear(&ax);
  raw.sign = u.sign ? raw.sign : !raw.sign;
  uint32_t result = g_bit(raw, N);
  gval_clear(&raw);
  return result;
}

//  ---- +asin:  arcsin(x) = atan(x/sqrt(1-x^2)) for |x|<1; +-pi/2 at +-1; NaR
//  outside [-1,1].  Composes existing correctly-rounded posit ops with the
//  new +atan -- faithful (a few ULP near |x|~1), not a dedicated kernel.
uint32_t PFN(_asin)(uint32_t x) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  uint32_t one = PFN(_one)();
  if ( LT(ABS(x), one) )
    return PFN(_atan)(DIV(x, SQT(SUB(one, MUL(x, x)))));
  if ( EQ(x, one) )      return MUL(PFN(_pi)(), DIV(one, SUN(2)));
  if ( EQ(x, NEG(one)) ) return NEG(MUL(PFN(_pi)(), DIV(one, SUN(2))));
  return NARBITS;
}

//  ---- +acos:  arccos(x) = atan(sqrt(1-x^2)/x) for x>0; pi-atan(sqrt(1-x^2)/
//  |x|) for x<0 (pi/2 at 0); 0 at 1, pi at -1; NaR outside [-1,1].  The x<0
//  branch is the FIX vs. the old naive version, which called
//  `(atan (div (sqt ...) x))` unconditionally -- wrong for x<0 (returns a
//  negative angle instead of the correct (pi/2,pi]).
uint32_t PFN(_acos)(uint32_t x) {
  if ( (x & MSK) == NARBITS ) return NARBITS;
  uint32_t one = PFN(_one)();
  if ( LT(ABS(x), one) ) {
    if ( EQ(x, 0) ) return MUL(PFN(_pi)(), DIV(one, SUN(2)));
    if ( LT(x, 0) )
      return SUB(PFN(_pi)(), PFN(_atan)(DIV(SQT(SUB(one, MUL(x, x))), NEG(x))));
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
