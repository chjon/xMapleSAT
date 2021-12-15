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

############################################
# New variable generation location options #
############################################
# Do not generate extension variables
ER_GEN_LOCATION_NONE           =0
# Generate extension variables after a restart
ER_GEN_LOCATION_AFTER_RESTART  =1
# Generate extension variables after a conflict
ER_GEN_LOCATION_AFTER_CONFLICT =2

##########################################
# New variable addition location options #
##########################################
# Do not add extension variables
ER_ADD_LOCATION_NONE           =0
# Add extension variables after a restart
ER_ADD_LOCATION_AFTER_RESTART  =1
# Add extension variables after a conflict
ER_ADD_LOCATION_AFTER_CONFLICT =2

define make_er_solver
	make -C core clean
	make CMD_CFLAGS="\
		-D ER_USER_ADD_HEURISTIC=$(2) \
		-D ER_USER_FILTER_HEURISTIC=$(3) \
		-D ER_USER_SELECT_HEURISTIC=$(4) \
		-D ER_USER_SUBSTITUTE_HEURISTIC=$(5) \
		-D EXTENSION_SUBSTITUTION=$(6) \
		-D ER_USER_DELETE_HEURISTIC=$(7) \
		-D ER_USER_GEN_LOCATION=$(8) \
		-D ER_USER_ADD_LOCATION=$(9) \
		-D EXTENSION_FORCE_BRANCHING=0 \
		-D ER_USER_ADD_SUBEXPR_SET_INTERSECTION=1 \
	" -C core rs
	mv core/maplesat_static core/$(1)
endef

all: maplesat xmaplesat

maplesat:
	$(call make_er_solver,maplesat,0,0,0,0,0,0,0,0)

xmaplesat: xmaplesat_rnd_rng xmaplesat_rnd_lbd xmaplesat_sub_rng xmaplesat_sub_lbd xmaplesat_glucoser

xmaplesat_rnd_rng:
	$(call make_er_solver,xmaplesat_rnd_rng,$(ER_ADD_HEURISTIC_RANDOM),$(ER_FILTER_HEURISTIC_RANGE),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_WIDTH),1,$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART))

xmaplesat_rnd_lbd:
	$(call make_er_solver,xmaplesat_rnd_lbd,$(ER_ADD_HEURISTIC_RANDOM),$(ER_FILTER_HEURISTIC_LBD),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_LBD),1,$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART))

xmaplesat_sub_rng:
	$(call make_er_solver,xmaplesat_sub_rng,$(ER_ADD_HEURISTIC_SUBEXPR),$(ER_FILTER_HEURISTIC_RANGE),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_WIDTH),1,$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART))

xmaplesat_sub_lbd:
	$(call make_er_solver,xmaplesat_sub_lbd,$(ER_ADD_HEURISTIC_SUBEXPR),$(ER_FILTER_HEURISTIC_LBD),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_LBD),1,$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART))

xmaplesat_glucoser:
	$(call make_er_solver,xmaplesat_glucoser,$(ER_ADD_HEURISTIC_GLUCOSER),$(ER_FILTER_HEURISTIC_GLUCOSER),$(ER_SELECT_HEURISTIC_GLUCOSER),$(ER_SUBSTITUTE_HEURISTIC_WIDTH),1,$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_CONFLICT),$(ER_ADD_LOCATION_AFTER_CONFLICT))

clean:
	make -C core clean
