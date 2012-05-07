Standalone programs and scripts
===============================

DIRAC provides several standalone executable codes and scripts which are
either useful for us and perhaps also useful for you.

pamq.py 
-------

This is script for submitting DIRAC for a specified queue system.

The current version supports few batch systems, like PBS, LoadLeveler,
and Grid Engine.

For example, save PBS as the default batch system, 1 node with 2 CPU's
(these are necessary parameters):

::

    pamq.py --config --qsys=pbs --nodes=1 --pnodes=2 

Submit and run DIRAC on 2 nodes in the 'express' queue (for 12 minutes).
One assumes that the queue system and number of CPU's per node are saved
in the configuration file:

::

    pamq.py --nodes=2 --time=00:12:00 --queue=express --pam_args="-noarch Mg.mol fsccd_IH.inp"

Only prepare script, do not sumbit

::

    pamq.py --nosub

For help, type

::

    pamq.py -h

twofit.x 
---------

This standalone executable serves for calculating spectroscopic parameters of diatomic molecules.



gpu_test.x
-----------

You can compile it only on the computer which is provided at least with GPU and appropriate CUDA 
library. In addition, the CULA library is welcome to enjoy further functionalities.


