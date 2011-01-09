#!/bin/sh

MOL=h2o2.mol MENU=point_charge.inp   ../SCRIPTS/run.sh
MOL=h2o2.mol MENU=gaussian_model.inp ../SCRIPTS/run.sh

rm -f *tgz

exit 0
