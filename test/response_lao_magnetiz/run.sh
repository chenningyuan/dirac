#!/bin/sh

export MOL=LiH.lsym.dir.mol

export MENU=LiH.dc_rkb.suscp_lao_go=cm
../SCRIPTS/run.sh

export MENU=LiH.dc_ukb.suscp_lao_go=cm
../SCRIPTS/run.sh

export MENU=LiH.levyle.suscp_lao_go=cm
../SCRIPTS/run.sh



