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

maplelcm_reset_weird:
	make CMD_CFLAGS="-D RANDOM_RESET=1 -D RANDOM_RESET_WEIRD=1" -C core rs
	mv core/glucose_static core/maplelcm_reset_weird

clean:
	make -C core clean
