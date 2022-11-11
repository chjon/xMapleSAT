MROOT := $(shell pwd)
export MROOT

maplesat:
	make -C core rs
	mv core/maplesat_static core/maplesat

maplesat_bcp_priority:
	make CMD_CFLAGS="-D BCP_PRIORITY=1" -C core rs 
	mv core/maplesat_static core/maplesat_bcp_priority

maplesat_er_priority:
	make CMD_CFLAGS="-D PRIORITIZE_ER=1" -C core rs 
	mv core/maplesat_static core/maplesat_er_priority

maplesat_er_priority_low:
	make CMD_CFLAGS="-D PRIORITIZE_ER=1 -D PRIORITIZE_ER_LOW=1" -C core rs 
	mv core/maplesat_static core/maplesat_er_priority_low

maplesat_bcp_er_priority:
	make CMD_CFLAGS="-D PRIORITIZE_ER=1 -D BCP_PRIORITY=1" -C core rs 
	mv core/maplesat_static core/maplesat_bcp_er_priority

clean:
	make -C core clean
