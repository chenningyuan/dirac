#!/bin/bash

if [[ ! -z "$1" ]]; then
    echo "WARNING: Code sections related to $* will be removed IN PLACE"
    echo "This script should be run from the Dirac root."
    echo "Continue? y/N"
    read cont
    if [[ "$cont" == "y" ]]; then
	CPPARG=
	while [[ ! -z "$1" ]]; do
	    CPPARG="$CPPARG -U$1"
	    shift
	done
	echo "Running utils/cppp.py $CPPARG on the source tree (may take some time)"
	python utils/cppp.py "$CPPARG" `find . -name '*.c' -or -name '*.F' -or -name '*.F90' -or -name '*.f90'`
	echo "Done"
    else
	echo "Aborted"
    fi
else
    echo "Usage: remove_modules MOD_FOO [MOD_BAR ..]"
    echo "Remove sections enclosed by #ifdef MOD_FOO (etc)"
    echo "from the current source tree."
fi