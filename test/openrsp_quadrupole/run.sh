#!/bin/bash

export MOL=fh.mol
export MENU=nr_hf.inp
../SCRIPTS/run.sh

# remove restart files
rm -f *tgz

exit 0
