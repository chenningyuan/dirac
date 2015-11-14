#!/bin/bash

BRANCH='master'
COMPILER='GNU'
PAMADM_FLAGS='--adm'
MAKE_RELEASE=0
BUILDNAME='undefined'
INT64=0
TRACK='Nightly'
DEPLOY_KEY='undefined'
NUM_JOBS=1

while getopts b:p:m:r:c:n:i:t:k:j: opt; do
  case $opt in
  b)
      BRANCH=$OPTARG
      ;;
  p)
      PAMADM_FLAGS=$OPTARG
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

ssh-agent bash -c "ssh-add ${DEPLOY_KEY}; git clone --recursive -b $BRANCH git@gitlab.com:dirac/dirac.git $CTEST_TEMP_DIR/compile"

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

./setup $SETUP_FLAGS --cmake-options="-DBUILDNAME=$BUILDNAME"

cd build
make -j$NUM_JOBS
ctest -D Nightly -j$NUM_JOBS --track $TRACK

rm -rf $CTEST_TEMP_DIR

exit 0
