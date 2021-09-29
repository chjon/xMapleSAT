MROOT := $(shell pwd)
export MROOT

define make_er_solver
	make -C core clean
	make CMD_CFLAGS="\
		-D EXTENSION_HEURISTIC=$(2) \
		-D ER_USER_FILTER_HEURISTIC=$(3) \
		-D ER_USER_SELECT_HEURISTIC=$(4) \
		-D EXTENSION_FORCE_BRANCHING=$(5) \
		-D EXTENSION_SUBSTITUTION=$(6) \
		-D ER_USER_ADD_SUBEXPR_SET_INTERSECTION=$(7) \
		-D ER_USER_SELECT_CACHE_ACTIVE_CLAUSES=$(8) \
	" -C core rs
	mv core/maplesat_static core/$(1)
endef

default:
	$(call make_er_solver,maplesat,0,0,1,0,0,0,0)

all:
	$(call make_er_solver,xmaplesat_rnd_rng_0,1,1,1,0,1,0,0)
	$(call make_er_solver,xmaplesat_sub_rng_0,2,1,1,0,1,0,0)
	$(call make_er_solver,xmaplesat_sub_rng_1,2,1,1,0,1,1,0)
	$(call make_er_solver,xmaplesat_rnd_lng_0,1,2,1,0,1,0,0)
	$(call make_er_solver,xmaplesat_sub_lng_0,2,2,1,0,1,0,0)
	$(call make_er_solver,xmaplesat_sub_lng_1,2,2,1,0,1,1,0)

clean:
	make -C core clean
