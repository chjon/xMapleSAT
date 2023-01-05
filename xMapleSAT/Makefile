MROOT := $(shell pwd)
export MROOT

###############################
#                             #
# Extended Resolution options #
#                             #
###############################

############################
# Filter heuristic options #
############################
# Consider all clauses
ER_FILTER_HEURISTIC_NONE     =0
# Consider clauses whose widths are in a certain range 
ER_FILTER_HEURISTIC_RANGE    =1
# Consider the longest clauses
ER_FILTER_HEURISTIC_LONGEST  =2
# Consider clauses whose LBDs are in a certain range
ER_FILTER_HEURISTIC_LBD      =3
# Consider the most recently learnt clauses
ER_FILTER_HEURISTIC_GLUCOSER =4

###############################
# Selection heuristic options #
###############################
# Consider all clauses
ER_SELECT_HEURISTIC_NONE      =0
# Select most active clauses
ER_SELECT_HEURISTIC_ACTIVITY  =1
# Select most active clauses using quickselect
ER_SELECT_HEURISTIC_ACTIVITY2 =2
# Only consider the previous two learnt clauses
ER_SELECT_HEURISTIC_GLUCOSER  =3

##################################
# Substitution heuristic options #
##################################
# Consider all clauses
ER_SUBSTITUTE_HEURISTIC_NONE  =0
# Consider clauses within a clause width range
ER_SUBSTITUTE_HEURISTIC_WIDTH =1
# Consider clauses within an LBD range
ER_SUBSTITUTE_HEURISTIC_LBD   =2
ER_SUBSTITUTE_HEURISTIC_BOTH  =3

#########################################
# Variable definition heuristic options #
#########################################
# Do not add extension variables
ER_ADD_HEURISTIC_NONE     =0
# Add extension variables by selecting random pairs of literals
ER_ADD_HEURISTIC_RANDOM   =1
# Add extension variables by selecting the most common pairs of literals
ER_ADD_HEURISTIC_SUBEXPR  =2
# Add extension variables according to the scheme prescribed by GlucosER
ER_ADD_HEURISTIC_GLUCOSER =3

#######################################
# Variable deletion heuristic options #
#######################################
# Do not delete extension variables
ER_DELETE_HEURISTIC_NONE      =0
# Delete all extension variables
ER_DELETE_HEURISTIC_ALL       =1
# Delete low-activity extension variables based on a constant activity threshold
ER_DELETE_HEURISTIC_ACTIVITY  =2
# Delete low-activity extension variables based on a proportional activity threshold
ER_DELETE_HEURISTIC_ACTIVITY2 =3

#########################################################
# New variable generation/introduction location options #
#########################################################
# Do not generate extension variables
NONE=0
# Generate extension variables after a restart
AFTER_RESTART=1
# Generate extension variables after a conflict
AFTER_CONFLICT=2

.PHONY: all clean test

all: maplesat xmaplesat

.prev_%:
	make -C core clean
	rm -f .prev_*
	touch $@

maplesat: .prev_maplesat
	make -C core rs || exit 1
	mkdir -p bin
	mv core/maplesat_static bin/$@

xmaplesat: .prev_xmaplesat
	make CMD_CFLAGS="\
		-D ER_USER_ADD_HEURISTIC=$(ER_ADD_HEURISTIC_SUBEXPR) \
		-D ER_USER_FILTER_HEURISTIC=$(ER_FILTER_HEURISTIC_LBD) \
		-D ER_USER_SELECT_HEURISTIC=$(ER_SELECT_HEURISTIC_ACTIVITY) \
		-D ER_USER_SUBSTITUTE_HEURISTIC=$(ER_SUBSTITUTE_HEURISTIC_BOTH) \
		-D ER_USER_DELETE_HEURISTIC=$(ER_DELETE_HEURISTIC_ACTIVITY2) \
		-D ER_USER_GEN_LOCATION=$(AFTER_RESTART) \
		-D ER_USER_ADD_LOCATION=$(AFTER_RESTART) \
		-D PRIORITIZE_ER=0 \
		-D BCP_PRIORITY_MODE=1 \
	" -C er rs || exit 1
	mkdir -p bin
	mv er/xmaplesat_static bin/$@

xmaplesat_rnd_rng: .prev_xmaplesat_rnd_rng
	make CMD_CFLAGS="\
		-D ER_USER_ADD_HEURISTIC=$(ER_ADD_HEURISTIC_RANDOM) \
		-D ER_USER_FILTER_HEURISTIC=$(ER_FILTER_HEURISTIC_RANGE) \
		-D ER_USER_SELECT_HEURISTIC=$(ER_SELECT_HEURISTIC_ACTIVITY) \
		-D ER_USER_SUBSTITUTE_HEURISTIC=$(ER_SUBSTITUTE_HEURISTIC_WIDTH) \
		-D ER_USER_DELETE_HEURISTIC=$(ER_DELETE_HEURISTIC_ACTIVITY2) \
		-D ER_USER_GEN_LOCATION=$(AFTER_RESTART) \
		-D ER_USER_ADD_LOCATION=$(AFTER_RESTART) \
	" -C er rs || exit 1
	mkdir -p bin
	mv er/xmaplesat_static bin/$@

xmaplesat_sub_lbd: .prev_xmaplesat_sub_lbd
	make CMD_CFLAGS="\
		-D ER_USER_ADD_HEURISTIC=$(ER_ADD_HEURISTIC_SUBEXPR) \
		-D ER_USER_FILTER_HEURISTIC=$(ER_FILTER_HEURISTIC_LBD) \
		-D ER_USER_SELECT_HEURISTIC=$(ER_SELECT_HEURISTIC_ACTIVITY) \
		-D ER_USER_SUBSTITUTE_HEURISTIC=$(ER_SUBSTITUTE_HEURISTIC_BOTH) \
		-D ER_USER_DELETE_HEURISTIC=$(ER_DELETE_HEURISTIC_ACTIVITY2) \
		-D ER_USER_GEN_LOCATION=$(AFTER_RESTART) \
		-D ER_USER_ADD_LOCATION=$(AFTER_RESTART) \
	" -C er rs || exit 1
	mkdir -p bin
	mv er/xmaplesat_static bin/$@

xmaplesat_glucoser: .prev_xmaplesat_glucoser
	make CMD_CFLAGS="-D ER_ENABLE_GLUCOSER=1" -C er rs || exit 1
	mkdir -p bin
	mv er/xmaplesat_static bin/$@

clean:
	make -C er   clean
	make -C core clean
	make -C test clean

test:
	make CMD_CFLAGS="-D TESTING" -C test