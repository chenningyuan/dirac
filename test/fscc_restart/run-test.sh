#!/usr/bin/bash

export my_nproc=4
export my_scratchdir=/state/partition1/gomes/dirac/runtest/gomes/cc_restart/fscc_restart
export my_bindir=/home/gomes/devel/Dirac/git/cc-restart/impi-4.1.0.024_intel13.0.1_lp64_cc-restart


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
