MROOT := $(shell pwd)
export MROOT

define make_er_solver
	make -C core clean
	make CMD_CFLAGS="\
		-D ER_USER_ADD_HEURISTIC=$(2) \
		-D ER_USER_FILTER_HEURISTIC=$(3) \
		-D ER_USER_SELECT_HEURISTIC=$(4) \
		-D ER_USER_SUBSTITUTE_HEURISTIC=$(5) \
		-D EXTENSION_SUBSTITUTION=1 \
		-D EXTENSION_FORCE_BRANCHING=0 \
		-D ER_USER_ADD_SUBEXPR_SET_INTERSECTION=1 \
	" -C core rs
	mv core/maplesat_static core/$(1)
endef

default:
	make -C core clean
	make -C core rs

xmaplesat:
	$(call make_er_solver,xmaplesat_rnd_rng,1,1,1,1)
	$(call make_er_solver,xmaplesat_rnd_lbd,1,3,1,3)
	$(call make_er_solver,xmaplesat_sub_rng,2,1,1,1)
	$(call make_er_solver,xmaplesat_sub_lbd,2,3,1,3)

glucoser:
	$(call make_er_solver,glucoser,3,4,3,0)

clean:
	make -C core clean
