MROOT := $(shell pwd)
export MROOT

define make_er_solver
	make -C core clean
	make CMD_CFLAGS="\
		-D EXTENSION_HEURISTIC=$(2) \
		-D ER_USER_FILTER_HEURISTIC=$(3) \
		-D ER_USER_SELECT_HEURISTIC=1 \
		-D ER_USER_SUBSTITUTE_HEURISTIC=3 \
		-D EXTENSION_SUBSTITUTION=1 \
		-D EXTENSION_FORCE_BRANCHING=0 \
		-D ER_USER_ADD_SUBEXPR_SET_INTERSECTION=1 \
		-D ER_USER_SELECT_CACHE_ACTIVE_CLAUSES=0 \
	" -C core rs
	mv core/maplesat_static core/$(1)
endef

default:
	$(call make_er_solver,maplesat,0,0,1,0,0,0,0)

all:
	$(call make_er_solver,xmaplesat_rnd_rng_lbd0,1,1)
	$(call make_er_solver,xmaplesat_rnd_lng_lbd0,1,2)
	$(call make_er_solver,xmaplesat_rnd_rng_lbd1,1,1)
	$(call make_er_solver,xmaplesat_rnd_lng_lbd1,1,2)
	$(call make_er_solver,xmaplesat_sub_rng_lbd0,2,1)
	$(call make_er_solver,xmaplesat_sub_lng_lbd0,2,2)
	$(call make_er_solver,xmaplesat_sub_rng_lbd1,2,1)
	$(call make_er_solver,xmaplesat_sub_lng_lbd1,2,2)

clean:
	make -C core clean
