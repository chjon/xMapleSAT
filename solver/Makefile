MROOT := $(shell pwd)
export MROOT

define make_er_solver
	make -C core clean
	make CMD_CFLAGS="\
		-D EXTENSION_HEURISTIC=$(2) \
		-D EXTENSION_FORCE_BRANCHING=$(3) \
		-D EXTENSION_SUBSTITUTION=$(4) \
		-D ER_USER_ADD_SUBEXPR_SET_INTERSECTION=$(5) \
		-D ER_USER_ADD_SUBEXPR_MAX_CHECK_WIDTH=$(6) \
		-D ER_USER_SELECT_CACHE_ACTIVE_CLAUSES=$(7) \
	" -C core rs
	mv core/maplesat_static core/$(1)
endef

default:
	$(call make_er_solver,maplesat,0,0,0,0,0,0)

all:
	$(call make_er_solver,xmaplesat_rnd_force0,1,0,0,1,50,0)
	$(call make_er_solver,xmaplesat_rnd_force1,1,1,0,1,50,0)
	$(call make_er_solver,xmaplesat_sub_force0,2,0,0,1,50,0)
	$(call make_er_solver,xmaplesat_sub_force1,2,1,0,1,50,0)

clean:
	make -C core clean
