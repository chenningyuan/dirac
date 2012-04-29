Manual:INTEGRALS:READIN
=======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*READIN -- Reading basis set
   file <#.2AREADIN_--_Reading_basis_set_file>`_

   -  `1.1 .UNCONTRACT <#.UNCONTRACT>`_
   -  `1.2 .PRINT <#.PRINT>`_
   -  `1.3 .MAXPRI <#.MAXPRI>`_

\*READIN -- Reading basis set file
==================================

This subsection allows changes of defaults in the reading of the `basis
file </dirac/index.php/Manual:Basis:Title_and_symmetry>`_.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.UNCONTRACT
-----------

Decontract basis sets specified as contracted (when being read from the
library, for example). In the case of using non-relativistically
contracted basis set the decontraction is necessary for heavy elements.
Two-component quasirelativistic Hamiltonians (like IOTC) work only with
decontracted basis set.

*Default:* Decontract only the small component.

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level in the reading of the basis file.

*Default:* Print level from
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.MAXPRI
-------

Maximum number of primitive functions in a given block in basis file.

*Default:*

::

    .MAXPRI
     22

