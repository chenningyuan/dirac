#!/bin/sh

export MOL=h2o2.mol
export MENU=gamma5.inp

../run.sh

rm -f *tgz

exit 0
