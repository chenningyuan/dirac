Manual:ANALYZE:PRIVEC
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*PRIVEC -- Print vectors <#.2APRIVEC_--_Print_vectors>`_

   -  `1.1 .VECPRI <#.VECPRI>`_
   -  `1.2 .PRICMP <#.PRICMP>`_
   -  `1.3 .AOLAB <#.AOLAB>`_
   -  `1.4 .PRINT <#.PRINT>`_

\*PRIVEC -- Print vectors
=========================

Print vectors read from DFCOEF.

.VECPRI
-------

For each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of orbitals to print.

*Default:* Print the occupied electronic solutions.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.PRICMP
-------

Separate control for printing of large and small components.

*Default:* Print large components.

::

    .PRICMP
     1

Print small components:

::

    .PRICMP
     2

.AOLAB
------

Print vectors in AO basis.

*Default:* Print vectors in SO basis.

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

