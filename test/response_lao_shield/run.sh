#!/bin/sh

export MOL=LiH.lsym.dir.mol

export MENU=LiH.levyle.shield_lao_go=cm
../run.sh

export MENU=LiH.bss.shield_lao_go=cm
../run.sh

export MENU=LiH.bss+mfsso.shield_lao_go=cm
../run.sh

export MENU=LiH.bss_sfB.shield_lao_go=cm
../run.sh

export MENU=LiH.dc_rkb.shield_lao_go=cm
../run.sh

export MENU=LiH.dc_sf.shield_lao_go=cm
../run.sh

