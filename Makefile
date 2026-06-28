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

#  Exhaustive bit-exactness check against the SoftPosit pX2 oracle.
test: $(SHLIB)
	$(ORACLE_PY) tools/oracle.py

clean:
	rm -f $(OBJS) libsoftunum.a $(SHLIB)

.PHONY: library test clean
