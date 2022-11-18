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

###############################
# Selection heuristic options #
###############################
# Consider all clauses
ER_SELECT_HEURISTIC_ALL       =0
# Select most active clauses
ER_SELECT_HEURISTIC_ACTIVITY  =1
# Select most active clauses using quickselect
ER_SELECT_HEURISTIC_ACTIVITY2 =2

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

##########################################
# New variable addition location options #
##########################################
# Do not add extension variables
ER_ADD_LOCATION_NONE           =0
# Add extension variables after a restart
ER_ADD_LOCATION_AFTER_RESTART  =1

define make_er_solver
	make CMD_CFLAGS="\
		-D ER_USER_ADD_HEURISTIC=$(2) \
		-D ER_USER_FILTER_HEURISTIC=$(3) \
		-D ER_USER_SELECT_HEURISTIC=$(4) \
		-D ER_USER_SUBSTITUTE_HEURISTIC=$(5) \
		-D ER_USER_DELETE_HEURISTIC=$(6) \
		-D ER_USER_GEN_LOCATION=$(7) \
		-D ER_USER_ADD_LOCATION=$(8) \
		-D ER_ENABLE_GLUCOSER=$(9) \
		$(10) \
	" -C core rs
	mv core/glucose_static core/$(1)
endef

.PHONY: all clean test

all: maplelcm xmaplelcm

.prev_%:
	make -C core clean
	rm -f .prev_*
	touch $@

maplelcm: .prev_maplelcm
	$(call make_er_solver,maplelcm,0,0,0,0,0,0,0,0)

xmaplelcm: xmaplelcm_rnd_rng xmaplelcm_rnd_lbd xmaplelcm_sub_rng xmaplelcm_sub_lbd xmaplelcm_glucoser

xmaplelcm_rnd_rng: .prev_xmaplelcm_rnd_rng
	$(call make_er_solver,$@,$(ER_ADD_HEURISTIC_RANDOM),$(ER_FILTER_HEURISTIC_RANGE),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_WIDTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART),0)

xmaplelcm_sub_lbd: .prev_xmaplelcm_sub_lbd
	$(call make_er_solver,$@,$(ER_ADD_HEURISTIC_SUBEXPR),$(ER_FILTER_HEURISTIC_LBD),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_BOTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART),0)

xmaplelcm_glucoser: .prev_xmaplelcm_glucoser
	$(call make_er_solver,$@,0,0,0,$(ER_SUBSTITUTE_HEURISTIC_WIDTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_NONE),$(ER_ADD_LOCATION_NONE),1)

xmaplelcm_rnd_rng_bcp_priority: .prev_xmaplelcm_rnd_rng_bcp_priority
	$(call make_er_solver,$@,$(ER_ADD_HEURISTIC_RANDOM),$(ER_FILTER_HEURISTIC_RANGE),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_WIDTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART),0,-D BCP_PRIORITY_MODE=1)

xmaplelcm_rnd_rng_bcp_priority2: .prev_xmaplelcm_rnd_rng_bcp_priority2
	$(call make_er_solver,$@,$(ER_ADD_HEURISTIC_RANDOM),$(ER_FILTER_HEURISTIC_RANGE),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_WIDTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART),0,-D BCP_PRIORITY_MODE=2)

xmaplelcm_sub_lbd_bcp_priority: .prev_xmaplelcm_sub_lbd_bcp_priority
	$(call make_er_solver,$@,$(ER_ADD_HEURISTIC_SUBEXPR),$(ER_FILTER_HEURISTIC_LBD),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_BOTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART),0,-D BCP_PRIORITY_MODE=1)

xmaplelcm_sub_lbd_bcp_priority2: .prev_xmaplelcm_sub_lbd_bcp_priority2
	$(call make_er_solver,$@,$(ER_ADD_HEURISTIC_SUBEXPR),$(ER_FILTER_HEURISTIC_LBD),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_BOTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART),0,-D BCP_PRIORITY_MODE=2)

xmaplelcm_glucoser_bcp_priority: .prev_xmaplelcm_glucoser_bcp_priority
	$(call make_er_solver,$@,0,0,0,$(ER_SUBSTITUTE_HEURISTIC_WIDTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_NONE),$(ER_ADD_LOCATION_NONE),1,-D BCP_PRIORITY_MODE=1)

xmaplelcm_rnd_rng_er_priority: .prev_xmaplelcm_rnd_rng_er_priority
	$(call make_er_solver,$@,$(ER_ADD_HEURISTIC_RANDOM),$(ER_FILTER_HEURISTIC_RANGE),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_WIDTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART),0,-D PRIORITIZE_ER=1)

xmaplelcm_sub_lbd_er_priority: .prev_xmaplelcm_sub_lbd_er_priority
	$(call make_er_solver,$@,$(ER_ADD_HEURISTIC_SUBEXPR),$(ER_FILTER_HEURISTIC_LBD),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_BOTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART),0,-D PRIORITIZE_ER=1)

xmaplelcm_glucoser_er_priority: .prev_xmaplelcm_glucoser_er_priority
	$(call make_er_solver,$@,0,0,0,$(ER_SUBSTITUTE_HEURISTIC_WIDTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_NONE),$(ER_ADD_LOCATION_NONE),1,-D PRIORITIZE_ER=1)

xmaplelcm_rnd_rng_er_bump: .prev_xmaplelcm_rnd_rng_er_bump
	$(call make_er_solver,$@,$(ER_ADD_HEURISTIC_RANDOM),$(ER_FILTER_HEURISTIC_RANGE),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_WIDTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART),0,-D BUMP_ER=1)

xmaplelcm_sub_lbd_er_bump: .prev_xmaplelcm_sub_lbd_er_bump
	$(call make_er_solver,$@,$(ER_ADD_HEURISTIC_SUBEXPR),$(ER_FILTER_HEURISTIC_LBD),$(ER_SELECT_HEURISTIC_ACTIVITY),$(ER_SUBSTITUTE_HEURISTIC_BOTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_AFTER_RESTART),$(ER_ADD_LOCATION_AFTER_RESTART),0,-D BUMP_ER=1)

xmaplelcm_glucoser_er_bump: .prev_xmaplelcm_glucoser_er_bump
	$(call make_er_solver,$@,0,0,0,$(ER_SUBSTITUTE_HEURISTIC_WIDTH),$(ER_DELETE_HEURISTIC_ACTIVITY2),$(ER_GEN_LOCATION_NONE),$(ER_ADD_LOCATION_NONE),1,-D BUMP_ER=1)

clean:
	make -C core clean
	make -C test clean

test:
	make CMD_CFLAGS="-D TESTING" -C test