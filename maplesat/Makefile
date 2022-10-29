MROOT := $(shell pwd)
export MROOT

maplesat:
	make -C core rs
	mv core/maplesat_static core/maplesat

maplesat_er_priority:
	make CMD_CFLAGS="-D PRIORITIZE_ER=1" -C core rs 
	mv core/maplesat_static core/maplesat_er_priority

clean:
	make -C core clean
