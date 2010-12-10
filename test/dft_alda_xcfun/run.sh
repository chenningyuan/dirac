#!/bin/bash

for mol in *mol; do
    for inp in *inp; do
        export MENU=$inp
        export MOL=$mol
        ../SCRIPTS/run.sh
    done
done

rm -f *tgz

exit 0
