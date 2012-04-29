Manual:DIRRCI:RASORB
====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &RASORB -- Specify the type of CI and the active
   space <#.26RASORB_--_Specify_the_type_of_CI_and_the_active_space>`_

   -  `1.1 NELEC <#NELEC>`_
   -  `1.2 NRAS1 <#NRAS1>`_
   -  `1.3 NRAS2 <#NRAS2>`_
   -  `1.4 MAXH1 <#MAXH1>`_
   -  `1.5 MAXE3 <#MAXE3>`_

&RASORB -- Specify the type of CI and the active space
======================================================

NELEC
-----

Number of electrons (excluding frozen core electrons).

*Default:*

::

    NELEC=0

NRAS1
-----

Number of spinors in the RAS1 space for each abelian irrep.

*Default:*

::

    NRAS1=NSYMRP*0

NRAS2
-----

Number of spinors in the RAS2 space for each abelian irrep.

*Default:*

::

    NRAS2=NSYMRP*0

MAXH1
-----

Maximum number of holes in RAS1 spinors.

*Default:*

::

    MAXH1=0

MAXE3
-----

Maximum number of electrons in RAS3 spinors.

*Default:*

::

    MAXE3=0

