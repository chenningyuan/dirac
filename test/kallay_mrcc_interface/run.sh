#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MOL=Be_D2h
         ;;
      "regular")
         MOL=Be_D2h
         ;;
      "small_benchmark")
         MOL=CO_C2v
         ;;
   esac
else
   MOL=Be_D2h
fi
MENU=ccsd;
export MENU;
export MOL;
../SCRIPTS/run.sh -mw=600 --noarch --get="MRCONEE MDCINT"
../../utils/dirac_mointegral_export.x >>${MENU}_$MOL.$TAG
$MRCCDIR/goldstone >>${MENU}_$MOL.$TAG
$MRCCDIR/mrcc      >>${MENU}_$MOL.$TAG
rm MRCONEE MDCINT fort.10 fort.14 fort.16 fort.55 fort.56 iface
