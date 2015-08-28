#!/bin/bash

BRANCH='master'
COMPILER='GNU'
PAMADM_FLAGS='--adm'
NUM_CORES=1
MAKE_RELEASE=0
BUILDNAME='undefined'
INT64=0
TRACK='Nightly'

while getopts b:p:m:r:c:n:i:t: opt; do
  case $opt in
  b)
      BRANCH=$OPTARG
      ;;
  p)
      PAMADM_FLAGS=$OPTARG
      ;;
  m)
      NUM_CORES=$OPTARG
      ;;
  r)
      MAKE_RELEASE=$OPTARG
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
  esac
done

shift $((OPTIND - 1))

CTEST_TEMP_DIR=/tmp/cdash-scratch
mkdir -p $CTEST_TEMP_DIR

export DIRAC_TMPDIR=$CTEST_TEMP_DIR/run

git clone --recursive -b $BRANCH git@gitlab.com:dirac/dirac.git $CTEST_TEMP_DIR/compile

cd $CTEST_TEMP_DIR/compile
./pamadm $PAMADM_FLAGS

if [ $MAKE_RELEASE -eq 1 ]; then
    mkdir -p build
    cd build
    cmake ..
    make release
    tar xvzf DIRAC-*-Source.tar.gz
    cd DIRAC-*-Source
fi

if [ "$COMPILER" == 'Intel' ]; then
    source /opt/intel/bin/compilervars.sh intel64
    SETUP_FLAGS="--fc=ifort --cc=icc --cxx=icpc --mkl=sequential"
else
    export MATH_ROOT=/opt/intel/mkl
    SETUP_FLAGS="--fc=gfortran --cc=gcc --cxx=g++"
fi

if [ $INT64 -eq 1 ]; then
    SETUP_FLAGS="$SETUP_FLAGS --int64"
fi

./setup.py --run-unit-tests $SETUP_FLAGS -D BUILDNAME=$BUILDNAME

cd build
make -j$NUM_CORES
ctest -D Nightly -j$NUM_CORES --track $TRACK

rm -rf $CTEST_TEMP_DIR

exit 0
