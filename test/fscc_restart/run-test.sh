#!/usr/bin/bash

export my_nproc=2
export my_scratchdir=/tmp/andre/test-dirac/cc_restart/fscc_restart
export my_bindir=/Users/gomes/devel/git/dirac_cc_restart/parallel_gfortran490_cc_restart


# first step for restarting, generate the sorted integrals and do some cc iterations
rm -rf $my_scratchdir/*
$my_bindir/pam --keep_scratch --scratchfull=$my_scratchdir --mol Mg.mol --inp  fsccsd_original.inp --mpi=$my_nproc

# first step for restarting, generate the sorted integrals and do some cc iterations
rm -rf $my_scratchdir/*
$my_bindir/pam --keep_scratch --scratchfull=$my_scratchdir --mol Mg.mol --inp  fsccsd_integrals.inp --mpi=$my_nproc

# now do the sectors
for sector in '00' '10' '01' '11' '20' '02'; do

   inpname="fsccsd_restart_s"$sector"_a.inp"
   $my_bindir/pam --keep_scratch --scratchfull=$my_scratchdir --mol Mg.mol --inp  $inpname --mpi=$my_nproc
   inpname="fsccsd_restart_s"$sector"_b.inp"
   $my_bindir/pam --keep_scratch --scratchfull=$my_scratchdir --mol Mg.mol --inp  $inpname --mpi=$my_nproc
done

rm *.tgz*
