MROOT := $(shell pwd)
export MROOT

maplelcm:
	make -C core rs
	mv core/glucose_static core/maplelcm

maplelcm_bcp_priority:
	make CMD_CFLAGS="-D BCP_PRIORITY=1" -C core rs 
	mv core/glucose_static core/maplelcm_bcp_priority

maplelcm_er_priority:
	make CMD_CFLAGS="-D PRIORITIZE_ER=1" -C core rs 
	mv core/glucose_static core/maplelcm_er_priority

maplelcm_bcp_er_priority:
	make CMD_CFLAGS="-D PRIORITIZE_ER=1 -D BCP_PRIORITY=1" -C core rs 
	mv core/glucose_static core/maplelcm_bcp_er_priority

clean:
	make -C core clean
