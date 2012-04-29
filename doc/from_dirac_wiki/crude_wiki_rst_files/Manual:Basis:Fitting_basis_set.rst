Manual:Basis:Fitting basis set
==============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

`|image0| </dirac/index.php/Image:Dirac.png>`_

**This page describes unreleased functionality. The keywords may not be
available in your version of DIRAC.**

Fitting basis set
=================

In the development version of DIRAC it is now possible to make use of
density fitting. In such calculations, it is mandatory to include a
basis set used to expand the density.

This is done adding a line with the FTSET keyword to the MOLECULE file,
after the coordinates for a given atom type

::

    FTSET [options]

Basically all the options available in the definition of large or small
basis sets are applicable to the fit set (i.e. the use of a set from the
library, explicitly typed basis, even/well-tempered etc). Here are some
examples:

Basis sets from the basis set library
-------------------------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE BASIS cc-pVDZ
    FTSET BASIS Ahlrichs-Coulomb-Fit

