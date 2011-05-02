#
# http://www.robvanderwoude.com/batchfiles.php
#
# for crontab&queue system you have to set up enviroment !

# the directory where your DIRAC sits
DIRAC=C:\Users\dv6-2145ec\Dirac-Windows\trunk

# go to the local DIRAC directory
cd $DIRAC

# radovan: bad idea, use svn co, not svn up
svn up

./pamadm -a --xcfun 

# configure step
cmake -G "MinGW Makefiles" "FC=gfortran"  -DENABLE_64BIT_INTEGERS=1 -DCMAKE_BUILD_TYPE=Debug  ..\

cd $DIRAC\build

make clean

make dirac.x

