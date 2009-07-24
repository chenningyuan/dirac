#!/bin/bash

for mol in h2o2.mol; do
    for inp in *inp; do
        export MENU=$inp
        export MOL=$mol
        ../run.sh
    done
done

exit 0
