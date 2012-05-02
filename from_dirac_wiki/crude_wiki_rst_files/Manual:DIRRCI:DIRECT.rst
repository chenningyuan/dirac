Manual:DIRRCI:DIRECT
====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &DIRECT -- Convergence
   control <#.26DIRECT_--_Convergence_control>`_

   -  `1.1 CONVERE <#CONVERE>`_
   -  `1.2 MAXITER <#MAXITER>`_
   -  `1.3 CONVERR <#CONVERR>`_
   -  `1.4 RESTART <#RESTART>`_
   -  `1.5 CPUMAX <#CPUMAX>`_

&DIRECT -- Convergence control
==============================

CONVERE
-------

Convergence threshold for energy.

*Default:*

::

    CONVERR=1.0D-9

MAXITER
-------

Maximum number of direct CI iterations.

*Default:*

::

    MAXITER=10

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

CONVERR
-------

Convergence threshold for residual vector.

*Default:*

::

    CONVERR=1.0D-10

RESTART
-------

Restart on CI-vectors present in MRCFINV.

*Default:*

::

    RESTART=F

CPUMAX
------

Maximum amount of CPU-seconds to be used.

*Default:*

::

    CPUMAX=604800

