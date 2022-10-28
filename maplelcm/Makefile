MROOT := $(shell pwd)
export MROOT

maplelcm:
	make -C core rs
	mv core/glucose_static core/maplelcm

maplelcm_random:
	make CMD_CFLAGS="-D RANDOM_DECISION=1" -C core rs
	mv core/glucose_static core/maplelcm_random

maplelcm_reset:
	make CMD_CFLAGS="-D RANDOM_RESET=1" -C core rs
	mv core/glucose_static core/maplelcm_reset

clean:
	make -C core clean
