MROOT := $(shell pwd)
export MROOT

maplesat:
	make -C core rs
	mv core/maplesat_static core/maplesat

maplesat_bcp_priority:
	make CMD_CFLAGS="-D BCP_PRIORITY_MODE=1" -C core rs 
	mv core/maplesat_static core/maplesat_bcp_priority

maplesat_bcp_priority2:
	make CMD_CFLAGS="-D BCP_PRIORITY_MODE=2" -C core rs 
	mv core/maplesat_static core/maplesat_bcp_priority2

maplesat_er_priority:
	make CMD_CFLAGS="-D PRIORITIZE_ER=1" -C core rs 
	mv core/maplesat_static core/maplesat_er_priority

maplesat_bcp_er_priority:
	make CMD_CFLAGS="-D PRIORITIZE_ER=1 -D BCP_PRIORITY_MODE=1" -C core rs 
	mv core/maplesat_static core/maplesat_bcp_er_priority

maplesat_bcp_er_priority2:
	make CMD_CFLAGS="-D PRIORITIZE_ER=1 -D BCP_PRIORITY_MODE=2" -C core rs 
	mv core/maplesat_static core/maplesat_bcp_er_priority2

clean:
	make -C core clean
