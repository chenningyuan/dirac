#!/bin/bash

for mol in *mol; do
    for inp in *inp; do
        export MENU=$inp
        export MOL=$mol
        ../run.sh -copy DALTON.HES
    done
done

exit 0
