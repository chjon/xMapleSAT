MROOT := $(shell pwd)
export MROOT

maplesat:
	make -C core rs

maplesat_reset:
	make CMD_CFLAGS="-D RANDOM_RESET=1" -C core rs

clean:
	make -C core clean
