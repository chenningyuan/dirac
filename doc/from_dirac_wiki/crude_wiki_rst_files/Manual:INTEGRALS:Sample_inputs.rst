Manual:INTEGRALS:Sample inputs
==============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Use the point nucleus model to compare with nonrelativistic programs
(most of them use point nucleus):

::

    **DIRAC
    .TITLE
      INTEGRALS sample input 1
    .WAVE FUNCTION
    **INTEGRALS
    .NUCMOD
     1
     [...]

 Read very large contracted basis set (from a basis set file, for
example) as uncontracted.

::

    **DIRAC
    .TITLE
      INTEGRALS sample input 2
    .WAVE FUNCTION
    **INTEGRALS
    *READIN
    .UNCONTRACT
    ! comment: for demonstration we have 30s25p18d/22s18p12d...
    ! ...consumes a lot of memory
    .MAXPRI
     30 
     [...]

