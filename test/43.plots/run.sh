#!/bin/sh
#
# Set up calculation
#  
MOL=NeGh
export MOL
  MENU=1d.density.dc
  export MENU
  ../run.sh -get "PLOT*"
  for file in PLOT*; do
    mv $file $MENU.$MOL.$file
  done
  MENU=1d.density.ll
  export MENU
  ../run.sh -get "PLOT*"
  for file in PLOT*; do
    mv $file $MENU.$MOL.$file
  done

MOL=NeGh.C1
export MOL
  MENU=1d.density.dc
  export MENU
  ../run.sh -get "PLOT*"
  for file in PLOT*; do
    mv $file $MENU.$MOL.$file
  done

MOL=ArGh
export MOL
  MENU=1d.elf.ll
  export MENU
  ../run.sh -get "PLOT*"
  for file in PLOT*; do
    mv $file $MENU.$MOL.$file
  done

MENU=3d.1s2s2p.density.ll
export MENU
  MOL=Ar
  export MOL
  ../run.sh -get "PLOT*"
  for file in PLOT*; do
    mv $file $MENU.$MOL.$file
  done
  MOL=Ne
  export MOL
  ../run.sh -get "PLOT*"
  for file in PLOT*; do
    mv $file $MENU.$MOL.$file
  done

MOL=Ne.C1
export MOL
  MENU=3d.density.ll
  export MENU
  ../run.sh -get "PLOT*"
  for file in PLOT*; do
    mv $file $MENU.$MOL.$file
  done
  MENU=3d.1s2s.density.ll.C1
  export MENU
  ../run.sh -get "PLOT*"
  for file in PLOT*; do
    mv $file $MENU.$MOL.$file
  done

MOL=Ne.C1
export MOL
  MENU=polarizability.dc
  export MENU
  ../run.sh -get "PAMXVC"
  MENU=3d.density.reaxvec.dc
  export MENU
  ../run.sh -copy "PAMXVC" -get "PLOT*"
  for file in PLOT*; do
    mv $file $MENU.$MOL.$file
  done
