#!/bin/bash

if [[ -z "$1" ]]; then
cat <<EOF
Script to execute another command with modified arguments.
Usage: replace_args_and_execute.sh SED_EXPR COMMAND [ARGS ..]
 Runs sed with SED_EXPR on all args, bulding a new 
 command line which is finally executed.
 Example, the xlf fortran compiler needs -WF,D instead
 of -D. Run it like
 replace_args_and_execute.sh s/^-D/-WF,D/g xlf [ARGS]
 This will automatically replace all -D arguments with
 -WF,D arguments.
By Ulf Ekstrom 2009
EOF
exit 0
fi

EXPR="$1"
shift
declare -a newargs
for a in "$@"; do
    newargs=( "${newargs[@]}" "`echo "$a" | sed $EXPR`" )
done
"${newargs[@]}"