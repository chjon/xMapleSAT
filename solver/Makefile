MROOT := $(shell pwd)
export MROOT

define make_er_solver
	make -C core clean
	make CMD_CFLAGS="\
		-D ER_USER_ADD_HEURISTIC=$(2) \
		-D ER_USER_FILTER_HEURISTIC=$(3) \
		-D ER_USER_SELECT_HEURISTIC=$(4) \
		-D ER_USER_SUBSTITUTE_HEURISTIC=$(5) \
		-D EXTENSION_SUBSTITUTION=$(6) \
		-D ER_USER_DELETE_HEURISTIC=$(7) \
		-D EXTENSION_FORCE_BRANCHING=0 \
		-D ER_USER_ADD_SUBEXPR_SET_INTERSECTION=1 \
	" -C core rs
	mv core/maplesat_static core/$(1)
endef

all: maplesat xmaplesat glucoser

maplesat:
	$(call make_er_solver,maplesat,0,0,0,0,0,0)

xmaplesat:
	$(call make_er_solver,xmaplesat_rnd_rng,1,1,1,1,1,3)
	$(call make_er_solver,xmaplesat_rnd_lbd,1,3,1,3,1,3)
	$(call make_er_solver,xmaplesat_sub_rng,2,1,1,1,1,3)
	$(call make_er_solver,xmaplesat_sub_lbd,2,3,1,3,1,3)

glucoser:
	$(call make_er_solver,glucoser_dcnst,3,4,3,1,1,2)
	$(call make_er_solver,glucoser_dfrac,3,4,3,1,1,3)

clean:
	make -C core clean
