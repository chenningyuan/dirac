#!/bin/ksh
if [ -z "$PAM_PATH" ]; then
   PAM_PATH=`which pam`
   if [ -n "`echo $PAM_PATH | grep 'Command not found'`" ]; then
      echo "pam not found. Please either set the enviroment variable"
      echo "PAM_PATH or put pam in your path."
      exit 1
   fi
   PAM_PATH=${PAM_PATH%/*}
fi
echo $PAM_PATH/pamq
if [ -z "$DIRAC" ]; then
   DIRAC=`which dirac.x`
fi
echo $DIRAC
if [ -z "$TAG" ]; then
   TAG=out
fi
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 U.dyall.2p 1el.dirac
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 U.dyall.2p 1el.diracsf
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 U.dyall.2p 1el.zora
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 U.dyall.2p 1el.zorasf
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 U.dyall.2p 1el.zorasc
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 U.dyall.2p 1el.zorascsf
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas nonrel
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas zora
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas zora4
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas dirac
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas zorasf
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas zora4sf
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas diracsf
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas zorasc
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas zora4sc
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas zorascsf
pamq $MPIARG -dirac $DIRAC -run $TAG -t 3600 Xe.Faas zora4scsf
