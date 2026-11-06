CFLAGS ?= -O2 -Wall -Wextra
LDLIBS = -lm
SANITIZE = -g -O1 -fsanitize=address,undefined -fno-sanitize-recover=undefined

all: black_scholes heston_calibrate

black_scholes: black_scholes.c common.h
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

heston_calibrate: heston_calibrate.c common.h
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

test:
	CC="$(CC)" ./regress.sh

sanitize:
	CC="$(CC)" REGRESS_CFLAGS="$(SANITIZE)" ./regress.sh

clean:
	rm -f black_scholes heston_calibrate

.PHONY: all test sanitize clean
