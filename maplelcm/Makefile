MROOT := $(shell pwd)
export MROOT

maplelcm:
	make -C core rs
	mv core/glucose_static core/maplelcm

maplelcm_er_priority:
	make CMD_CFLAGS="-D PRIORITIZE_ER=1" -C core rs 
	mv core/glucose_static core/maplelcm_er_priority

clean:
	make -C core clean
