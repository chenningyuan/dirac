Standalone programs and scripts
===============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

DIRAC provides several standalone executable codes and scripts which are
either useful for us and perhaps also useful for you.

Contents
--------

-  `1 pamq.py script <#pamq.py_script>`_
-  `2 twofit.x executable <#twofit.x_executable>`_

pamq.py script
--------------

Script for submitting DIRAC for a specified queue system.

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

twofit.x executable
-------------------

Serves for calculating spectroscopic parameters of diatomic molecules.

