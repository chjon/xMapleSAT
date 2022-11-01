MROOT := $(shell pwd)
export MROOT

maplesat:
	make -C core rs
	mv core/maplesat_static core/maplesat

maplesat_random:
	make CMD_CFLAGS="-D RANDOM_DECISION=1" -C core rs
	mv core/maplesat_static core/maplesat_random

maplesat_reset:
	make CMD_CFLAGS="-D RANDOM_RESET=1" -C core rs
	mv core/maplesat_static core/maplesat_reset

maplesat_reset_weird:
	make CMD_CFLAGS="-D RANDOM_RESET=1 -D RANDOM_RESET_WEIRD=1" -C core rs
	mv core/maplesat_static core/maplesat_reset_weird

clean:
	make -C core clean
