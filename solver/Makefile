MROOT := $(shell pwd)
export MROOT

curr_solver:
	make -C core rs

all_solvers:
	make -C core clean
	make CMD_CFLAGS="-D EXTENSION_HEURISTIC=1 -D EXTENSION_FORCE_BRANCHING=0" -C core rs
	mv core/maplesat_static core/xmaplesat_rnd_force0
	make -C core clean
	make CMD_CFLAGS="-D EXTENSION_HEURISTIC=1 -D EXTENSION_FORCE_BRANCHING=1" -C core rs
	mv core/maplesat_static core/xmaplesat_rnd_force1
	make -C core clean
	make CMD_CFLAGS="-D EXTENSION_HEURISTIC=2 -D EXTENSION_FORCE_BRANCHING=0" -C core rs
	mv core/maplesat_static core/xmaplesat_sub_force0
	make -C core clean
	make CMD_CFLAGS="-D EXTENSION_HEURISTIC=2 -D EXTENSION_FORCE_BRANCHING=1" -C core rs
	mv core/maplesat_static core/xmaplesat_sub_force1

clean:
	make -C core clean
