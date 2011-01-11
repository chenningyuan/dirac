#!/bin/bash

for mol in *.mol; do
    for inp in *.inp; do
        MENU=$inp MOL=$mol ../SCRIPTS/run.sh
    done
done
rm -f *tgz

exit 0
