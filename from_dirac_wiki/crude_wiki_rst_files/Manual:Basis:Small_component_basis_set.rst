Manual:Basis:Small component basis set
======================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Small component basis set <#Small_component_basis_set>`_

   -  `1.1 By kinetic balance <#By_kinetic_balance>`_
   -  `1.2 Explicitly typed basis sets <#Explicitly_typed_basis_sets>`_
   -  `1.3 MOLFDIR-type basis sets <#MOLFDIR-type_basis_sets>`_

Small component basis set
=========================

After the large component basis set the small component basis set *can*
be specified. If nothing is specified it is equivalent to specifying
SMALL KINBAL (see below). There are three possibilities for giving the
basis set.

By kinetic balance
------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE BASIS aug-cc-pVDZ
    SMALL KINBAL

The small component basis is generated using kinetic balance. The small
component basis set will *always* be uncontracted, independent of
whether the large component basis set is contracted or not.

Explicitly typed basis sets
---------------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE BASIS aug-cc-pVDZ
    SMALL EXPLICIT  3    1    1    1
    ...

Syntax is analogous to `LARGE
EXPLICIT </dirac/index.php/Manual:Basis:Large_component_basis_set#Explicitly_typed_basis_sets>`_.

MOLFDIR-type basis sets
-----------------------

::

    DIRAC
    Water
    MOLFDIR basis
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE MOLFBAS Oxygen-xyz.bas
    SMALL MOLFBAS Oxygen-xyz.bas

The MOLFDIR basis set file(s) (here: Oxygen-xyz.bas) must be copied to
the scratch area, for example using the pam script:

::

    pam ... -copy "Oxygen-xyz.bas ..."

