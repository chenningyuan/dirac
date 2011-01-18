#!/bin/sh

for molecule in *mol; do
    export MOL=$molecule

    for input in *inp; do
        export MENU=$input

        ../SCRIPTS/run.sh --get=plot.3d.cube
        mv plot.3d.cube ${molecule%\.mol}_${input%\.inp}.3d

    done
done

exit 0
