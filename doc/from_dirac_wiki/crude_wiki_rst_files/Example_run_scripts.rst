Example run scripts
===================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Linux cluster "STALLO" at the University of Tromsø,
   Norway <#Linux_cluster_.22STALLO.22_at_the_University_of_Troms.C3.B8.2C_Norway>`_

Linux cluster "STALLO" at the University of Tromsø, Norway
----------------------------------------------------------

::

    #!/bin/bash

    # this is an example 8 processor job for 5 minutes
    # no commands before the PBS directives
    # please adapt [things in brackets]

    #PBS -S /bin/bash
    #PBS -A [your account in here]
    #PBS -l nodes=8:ppn=1:gige
    #PBS -l walltime=00:05:00
    #PBS -l pvmem=3800mb

    module load openmpi/1.3.2
    export LD_LIBRARY_PATH=/global/apps/intel/mkl/10.1.0.015/lib/em64t:$LD_LIBRARY_PATH
    export OMP_NUM_THREADS=1

    workdir=/global/work/$USER/$PBS_JOBID
    mkdir -p $workdir

    cd $PBS_O_WORKDIR
    /[edit path to pam]/pam -keepwrk -noarch -wrkdir $workdir -mpi 8 -mem 1700 -slavemem 1700 -aw 5 molecule.mol input.inp

    cd $PBS_O_WORKDIR
    rm -rf $workdir

    exit 0

