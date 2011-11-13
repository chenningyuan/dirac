#!/bin/sh

for molecule in *spherical.mol; do
    for input in *spherical.inp; do
        export MOL=$molecule
        export MENU=$input
        ../run.sh
    done
done

for molecule in *cartesian.mol; do
    for input in *cartesian.inp; do
        export MOL=$molecule
        export MENU=$input
        ../run.sh
    done
done

exit 0
