#!/bin/sh
#
# Set up calculation
#

if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
	 MOL=CH3F.xyz;
	 MENU=xyz.inp;
	 export MOL;
	 export MENU;
	 ../run.sh
	 ;;

      "medium")
         MOL=CH3F.xyz;
	 MENU=xyz.inp;
	 export MOL;
	 export MENU;
	 ../run.sh
         ;;

      "long")

         MOL=CH3F.xyz;
	 MENU=xyz.inp;
	 export MOL;
	 export MENU;
	 ../run.sh;

	 MOL=CH3F.xyz;
	 MENU=xyz-2.inp;
	 export MOL;
	 export MENU;
	 ../run.sh 
         ;;
   esac
else
    MOL=CH3F.xyz;
    MENU=xyz.inp;
    export MOL;
    export MENU;
    ../run.sh

fi
      


	 



