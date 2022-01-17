#!/bin/bash

#Sometimes one needs to set additional environment variables, on MACOSX the following can be used:
#source /Users/visscher/.bash_profile

BRANCH='master'
COMPILER='GNU'
BUILDNAME='undefined'
INT64=0
TRACK='master'
DEPLOY_KEY='undefined'
#Example for mac, replace by proper path and uncomment
#DEPLOY_KEY='/Users/visscher/.ssh/nightly_rsa'
#Assuming that the machine is empty and you have 4 or more cores available, reduce or increase NUM_JOBS as desired.
NUM_JOBS=4

while getopts b:r:c:n:i:t:k:j: opt; do
  case $opt in
  b)
      BRANCH=$OPTARG
      ;;
  c)
      COMPILER=$OPTARG
      ;;
  n)
      BUILDNAME=$OPTARG
      ;;
  i)
      INT64=$OPTARG
      ;;
  t)
      TRACK=$OPTARG
      ;;
  k)
      DEPLOY_KEY=$OPTARG
      ;;
  j)
      NUM_JOBS=$OPTARG
      ;;
  esac
done

shift $((OPTIND - 1))

CTEST_TEMP_DIR=/tmp/cdash-scratch
mkdir -p $CTEST_TEMP_DIR

export DIRAC_TMPDIR=$CTEST_TEMP_DIR/run

ssh-agent bash -c "ssh-add ${DEPLOY_KEY}; git clone --recursive -b $BRANCH git@gitlab.com:dirac/dirac-private.git $CTEST_TEMP_DIR/compile"

cd $CTEST_TEMP_DIR/compile

if [ "$COMPILER" == 'Intel' ]; then
    source /opt/intel/bin/compilervars.sh intel64
    SETUP_FLAGS="--fc=ifort --cc=icc --cxx=icpc --mkl=parallel"
else
    export DIRAC_MPI_COMMAND="mpirun -np 4"
    export MATH_ROOT=/opt/intel/mkl
    SETUP_FLAGS="--fc=gfortran --cc=gcc --cxx=g++ --mpi"
fi

if [ $INT64 -eq 1 ]; then
    SETUP_FLAGS="$SETUP_FLAGS --int64"
fi

./setup $SETUP_FLAGS --cmake-options="-DBUILDNAME=$BUILDNAME"

cd build
make -j$NUM_JOBS
ctest -D Nightly -j$NUM_JOBS --track $TRACK

rm -rf $CTEST_TEMP_DIR

exit 0
