#!/bin/sh

for molecule in *mol; do
    export MOL=$molecule

    for input in *inp; do
        export MENU=$input

        ../run.sh -get PLOT_3D
        mv PLOT_3D ${molecule%\.mol}_${input%\.inp}.3d

    done
done

exit 0
