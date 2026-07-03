#  SoftUnum -- 2022 Posit Standard (unums) in pure-integer C.  No SoftFloat,
#  no floating point: every value is a raw bit pattern.

CC ?= cc
#  The transcendentals (src/posit/pgmp.h + ptrans.h) do their g-layer combine
#  in arbitrary precision via GMP -- see pgmp.h's module comment for why a
#  fixed-width int can't stand in for Hoon's unbounded `@` atoms here.
GMP_PREFIX ?= /opt/homebrew
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -Iinclude -I$(GMP_PREFIX)/include
LDFLAGS ?= -L$(GMP_PREFIX)/lib
LDLIBS  ?= -lgmp

SRCS = src/posit/p8.c src/posit/p16.c src/posit/p32.c
OBJS = $(SRCS:.c=.o)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHLIB   = libsoftunum.dylib
  SHFLAGS = -dynamiclib
else
  SHLIB   = libsoftunum.so
  SHFLAGS = -shared
endif

ORACLE_PY ?= /opt/anaconda3/bin/python

.DEFAULT_GOAL := library

library: libsoftunum.a $(SHLIB)

libsoftunum.a: $(OBJS)
	ar rcs $@ $^

$(SHLIB): $(SRCS)
	$(CC) $(CFLAGS) -fPIC $(SHFLAGS) -o $@ $(SRCS) $(LDFLAGS) $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

#  Bit-exactness checks: (1) curated /lib/unum Hoon vectors (fast, independent
#  of SoftPosit -- covers the transcendentals), (2) the exhaustive/sampled
#  SoftPosit pX2 + p32 oracle for the core ops and the mpmath-verified
#  unum_cheb_check.py model for the transcendentals, (3) a large from-scratch
#  vector corpus (tools/gen_cheb_vectors.py + cheb_check.c) linked directly
#  against the C ABI (no ctypes) for the new range-reduced transcendentals.
test: $(SHLIB) test-cheb
	$(ORACLE_PY) tools/hoon_vectors.py
	$(ORACLE_PY) tools/ieee_check.py
	$(ORACLE_PY) tools/quire_check.py
	$(ORACLE_PY) tools/oracle.py

#  Just the fast, dependency-light Hoon-vector check.
test-vectors: $(SHLIB)
	$(ORACLE_PY) tools/hoon_vectors.py

tools/cheb_check: tools/cheb_check.c libsoftunum.a
	$(CC) $(CFLAGS) -o $@ tools/cheb_check.c libsoftunum.a $(LDFLAGS) $(LDLIBS)

tools/cheb_vectors.txt:
	$(ORACLE_PY) tools/gen_cheb_vectors.py tools/cheb_vectors.txt

#  Large from-scratch bit-exactness sweep for exp/log/log2/log10/sin/cos/tan/
#  atan/asin/acos against libmath/tools/unum_cheb_check.py (independent of
#  both the Hoon vectors and oracle.py's own transliteration).
test-cheb: tools/cheb_check tools/cheb_vectors.txt
	./tools/cheb_check tools/cheb_vectors.txt

clean:
	rm -f $(OBJS) libsoftunum.a $(SHLIB) tools/cheb_check tools/cheb_vectors.txt

.PHONY: library test test-vectors test-cheb clean
