#!/bin/sh

# 1d

MOL=negh.auto; export MOL
MENU=1d.density.dc; export MENU
../run.sh -get "PLOT*"; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done

MOL=negh.auto; export MOL
MENU=1d.density.ll; export MENU
../run.sh -get "PLOT*"; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done

MOL=negh.c1; export MOL
MENU=1d.density.dc; export MENU
../run.sh -get "PLOT*"; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done

MOL=argh.auto; export MOL
MENU=1d.elf.ll; export MENU
../run.sh -get "PLOT*"; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done

# 3d

MOL=ar.auto; export MOL
MENU=3d.1s2s2p.density.ll; export MENU
../run.sh -get "PLOT*"; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done

MOL=ne.auto; export MOL
MENU=3d.1s2s2p.density.ll; export MENU
../run.sh -get "PLOT*"; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done

MOL=ne.c1; export MOL
MENU=3d.density.ll; export MENU
../run.sh -get "PLOT*"; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done

MOL=ne.c1; export MOL
MENU=3d.1s2s.density.ll.c1; export MENU
../run.sh -get "PLOT*"; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done

# 3d T+ perturbed densities

property=polarizability
for symmetry in auto c1 c2 c2h c2v ci cs d2 d2h; do
  MOL=ne.$symmetry; export MOL
  MENU=$property.dc; export MENU
  ../run.sh -get "PAMXVC" -outcmo

  quantity=density
  MENU=3d.$quantity.readxvec.dc.$symmetry; export MENU
  ../run.sh -copy "PAMXVC" -get "PLOT*" -incmo; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done
done

# 3d T- perturbed densities

property=magnetizability
for symmetry in auto c1 c2 c2h c2v ci cs d2 d2h; do
  MOL=ne.$symmetry; export MOL
  MENU=$property.dc; export MENU
  ../run.sh -get "PAMXVC" -outcmo

  quantity=spindensity
  MENU=3d.$quantity.readxvec.dc.$symmetry; export MENU
  ../run.sh -copy "PAMXVC" -get "PLOT*" -incmo; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done

  quantity=currentdensity_xy
  MENU=3d.$quantity.readxvec.dc.$symmetry; export MENU
  ../run.sh -copy "PAMXVC" -get "PLOT*" -incmo; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done

  quantity=currentdensity_xz
  MENU=3d.$quantity.readxvec.dc.$symmetry; export MENU
  ../run.sh -copy "PAMXVC" -get "PLOT*" -incmo; for file in PLOT*; do echo $file; mv $file $MENU\_$MOL.PLOT; done
done

exit 0
