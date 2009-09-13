#!/bin/sh

for molecule in ne*mol; do
    export MOL=${molecule%\.mol}

    export MENU=response_ez_ez.inp
    ../run.sh -outcmo -get  PAMXVC
    export MENU=visual_zo.inp
    ../run.sh -incmo  -copy PAMXVC
    export MENU=visual_ez_ez.inp
    ../run.sh -incmo  -copy PAMXVC
done

exit 0
