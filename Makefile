#  SoftUnum -- 2022 Posit Standard (unums) in pure-integer C.  No SoftFloat,
#  no floating point: every value is a raw bit pattern.

CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Iinclude

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
	$(CC) $(CFLAGS) -fPIC $(SHFLAGS) -o $@ $(SRCS)

%.o: %.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

#  Bit-exactness checks: (1) curated /lib/unum Hoon vectors (fast, independent
#  of SoftPosit -- covers the transcendentals), (2) the exhaustive/sampled
#  SoftPosit pX2 + p32 oracle for the core ops and a Hoon-faithful reference
#  for the transcendentals.
test: $(SHLIB)
	$(ORACLE_PY) tools/hoon_vectors.py
	$(ORACLE_PY) tools/ieee_check.py
	$(ORACLE_PY) tools/quire_check.py
	$(ORACLE_PY) tools/oracle.py

#  Just the fast, dependency-light Hoon-vector check.
test-vectors: $(SHLIB)
	$(ORACLE_PY) tools/hoon_vectors.py

clean:
	rm -f $(OBJS) libsoftunum.a $(SHLIB)

.PHONY: library test test-vectors clean
