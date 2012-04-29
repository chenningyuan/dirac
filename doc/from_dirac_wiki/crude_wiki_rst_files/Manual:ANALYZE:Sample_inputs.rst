Manual:ANALYZE:Sample inputs
============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Print all orbitals in fermion irrep 1 and some in irrep 2, large
components only and using AO and not the SO basis:

::

    **ANALYZE
    .PRIVEC
    *PRIVEC
    .VECPRI     # what molecular orbitals to print
     all
    -4..2,4
    .PRICMP     # print only large components
     1 0
    .AOLAB      # use AO basis

Mulliken population analysis defining labels for each atomic shell and
specifying orbital to analyze:

::

    **ANALYZE
    .MULPOP
    *MULPOP
    .LABDEF     # redefinition of labels
    3
    H(s)        1
    H(p)        2..4
    H(d)        5..10
    .VECPOP     # what molecular orbitals to analyze
    1..10
    1,2

Projection analysis using two fragments, represented by coefficient
files DFFRG1 and DFFRG2, respectively, and each calculated in their own
basis:

::

    **ANALYZE
    ..PROJEC
    *PROJECTION
    .OWNBAS     # use the individual basis set of fragments
    .VECPRJ     # what molecular orbitals to analyze
    1,2
    1..3
    .VECREF     # definition of reference orbitals
     2
     DFFRG1
     all
     all
     DFFRG2
     1..4
     1..2

