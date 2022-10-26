MROOT := $(shell pwd)
export MROOT

maplesat:
	make -C core rs
	mv core/maplesat_static core/maplesat

maplesat_reset:
	make CMD_CFLAGS="-D RANDOM_RESET=1" -C core rs
	mv core/maplesat_static core/maplesat_reset

clean:
	make -C core clean
