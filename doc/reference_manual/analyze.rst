\*\*ANALYZE -- analyzing tools
==============================

This section is aimed to analyze the final Hartree-Fock wave function by
using one or more analysis modules. By default none of them is
activated.

.PRIVEC
-------

Print vectors. Activates the
`\*PRIVEC </dirac/index.php/Manual:ANALYZE:PRIVEC>`_ subsection.


\*PRIVEC -- Print vectors
=========================

Print vectors read from DFCOEF.

.VECPRI
-------

For each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of orbitals to print.

*Default:* Print the occupied electronic solutions.

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

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0


.MULPOP
-------

Perform Mulliken population analysis [Mulliken1955]_
Continues to
`\*MULPOP </dirac/index.php/Manual:ANALYZE:MULPOP>`_ subsection.

\*MULPOP -- Mulliken population analysis
========================================

Mulliken population analysis is performed in AO basis. The analysis is
based on the concept of *labels*. Each basis function is labeled by its
functional type and center (and boson irrep when in SO basis, the
default). The labels are given in the output. A set of primitive labels
may be collected to group labels as specified by the user.

.VECPOP
-------

For each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of orbitals to
analyze.

*Default:* Analyze the occupied electronic solutions.

**Advanced options**

.NETPOP
-------

Give Mulliken net/overlap populations in addition to gross populations.

.LABDEF
-------

Define labels for use in Mulliken population analysis.

Give first the number of labels to define, then the label name (12
characters) and `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ for each fermion
irrep.

*Example:*

::

    .LABDEF
    4 (total number of lines to follow)
    Re s        1
    Re p        2..4
    Re d        5..10
    Re f        11..20
     ...        21..35

*Pitfall:*

If you decide to redefine labels for the population analysis you have to
reassign **all** of them, otherwise Dirac is unable to continue. The
list of functions with the initial label definitions (*default:
SO-labels*) can be drawn from the DIRAC output starting at the section :

::

                                   GETLAB: SO-labels
                                                     
       Large components:   35
        1  L Ag AR s        2  L Ag AR dxx      3  L Ag AR dyy      4  L Ag AR dzz      etc.

.AOLAB
------

Base definition of labels on AO basis.

*Default:* Base definition of labels on SO basis.

.INDSML
-------

Use individual small component labels.

*Default:* Gather all small component functions belonging to a given
center and irrep, irrespective of function type.

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0


.PROJEC
-------

Perform projection analysis [Faegri2001]_
Activates the
`\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_ subsection.

\*PROJECTION -- Projection analysis
===================================

The current solution may be projected down onto another set of
coefficients generated from the basis, e.g. one may project molecular
solutions down onto atomic solutions in order to evaluate atomic
contributions. The coefficients of the fragments are read simultaneously
and the overlap between them is taken into account [Faegri2001]_ .

The projection analysis can be thought of as a Mulliken analysis based
on the atomic (fragment) orbitals; it is therefore very much less
sensitive to the basis set.

Normally, the fragments are calculated in the full molecular basis. This
is done by zeroing charges of the remaining atoms of the molecule in the
basis set input and adjusting occupation in the menu file. It is,
however, possible to calculate the fragments in the own basis (a subset
of the full molecular basis) using `.OWNBAS <#.OWNBAS>`_. The advantage
is faster calculations and conservation of atomic symmetry for atomic
fragments. To avoid working with symmetry-combinations of atomic centers
for the fragments, it may sometimes be advantageous to dump molecular
coefficients in *C*\ :sub:`1`\  symmetry using
`.ACMOUT </dirac/index.php/Manual:GENERAL#.ACMOUT>`_ 
under
`\*\*GENERAL </dirac/index.php/Manual:GENERAL>`_ 
and do the analysis
without symmetry. When 
`.OWNBAS <#.OWNBAS>`_ 
is used, one then needs to
calculate each atomic type only once in its own basis in order to do the
complete analysis.

.VECPRJ
-------

For each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of orbitals to
analyze.

*Default:* Analyze the occupied electronic solutions.

.VECREF
-------

First give number of fragments to project onto, and then for each
fragment give filename of MO coefficients and the number of
symmetry-independent nuclei in this fragment and for each fermion irrep,
give an `orbital string </dirac/index.php/Manual:Orbital_strings>`_ of
reference orbitals.

*Example:*

::

    .VECREF
    4
    H00001
    1
    1
    H00002
    1
    1
    X00003
    1
    1..43
    X00004
    1
    1..43

.OWNBAS
-------

Calculate fragments in their own basis.

This keyword must be used with some care as the list of fragments is
assumed to be identical to that of symmetry independent centers.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.PROTHR
-------

Set threshold for absolute value of projection coefficients to be
printed.

*Default:*

::

    .PROTHR
     0.01

.WGPOP
------

Split overlap densities according to weight of contributions.

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.DENSIT
-------

Write density to a formatted file in Gaussian cube format. Activates the
`\*DENSITY </dirac/index.php/Manual:ANALYZE:DENSITY>`_ subsection.

\*DENSITY -- Write the density to a formatted file
==================================================

The density can be dumped to formatted file in `CUBE file
format <http://local.wasp.uwa.edu.au/~pbourke/dataformats/cube/>`_ for
subsequent visualization (e.g. using
`Molekel <http://www.cscs.ch/molekel/>`_ or
`MOLDEN <http://www.cmbi.ru.nl/molden/molden.html>`_ or other programs;
see also `POV-Ray <http://www.povray.org/>`_).

The large and small component density is written to files rhol.cube and
rhoS.cube, respectively. They can be recovered from the work area by the
command

::

    pam -get ``rhoL.cube rhoS.cube ...

Note also that the coefficient file DFCOEF must be present in order to
generate the density.

::

    pam --incmo ...

Large component density of iodobenzene plotted using Molekel (cutoff:
0.01). 
</dirac/index.php/Image:Iodobenzene_l.jpg>

Large component density of iodobenzene plotted using
`Molekel <http://www.cscs.ch/molekel/>`_ (cutoff: 0.01).

Small component density of iodobenzene plotted using Molekel (cutoff:
0.0001). 
</dirac/index.php/Image:Iodobenzene_s.jpg>

Small component density of iodobenzene plotted using
`Molekel <http://www.cscs.ch/molekel/>`_ (cutoff: 0.0001).

It is also possible to plot the density of an orbital by subtracting
densities generated from a given coefficient file, varying the
occupation. The figure shows the large component density of the HOMO of
the uranyl cation obtained by subtracting from the density generated
from regular occupation the density generated by reducing the occupation
by two electrons. 

</dirac/index.php/Image:Uo2_homo.jpg>

It is also possible to plot the density of an orbital by subtracting
densities generated from a given coefficient file, varying the
occupation. The figure shows the large component density of the HOMO of
the uranyl cation obtained by subtracting from the density generated
from regular occupation the density generated by reducing the occupation
by two electrons.

.ORBITALS
---------

This keyword is only available in the Development version of Dirac. Make
individual cube files for each give orbital, instead of the total
density.

*Example:*

::

    .ORBITALS
     4..9
     1,3

This generates cube files for gerade orbitals 4-9 and ungerade orbitals
1 and 3. The numbering starts from the first occupied orbital, so
"positronic" orbitals have negative indices. Be aware that the cube
generation can be rather slow for large basis sets. Writing hundreds of
cube files can take much more time than a SCF calculation.

.DOCUBE
-------

Specify large/small component density generation (1 = on; 0 = off).

*Default:*

::

    .DOCUBE
     1 1

.NCUBE
------

Number of cells along the sides of the cube.

*Default:*

::

    .CUBE
     80 80 80

.CUBADJ
-------

Adjust the volume of the cube.

*Arguments:* Real CUBADJ(1), CUBADJ(2)

*Default:*

::

    .CUBADJ
     4.0 8.0

DIRAC will look at atomic positions and then try to place the cube such
that most density is included. The origin of the cube is determined by
finding the minimum atomic position in the x-, y-and z-direction
separately and then add CUBADJ(1). The length of the vectors spanning
the cube (parallelepiped) are analogously defined by maximal atomic
positions, with CUBADJ(2) added.

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0


.LOCALI
-------

Localize orbitals using the Pipek-Mezey criterion [Dubillard2006]_ .

Continues to the activated
`\*LOCALIZATION </dirac/index.php/Manual:ANALYZE:LOCALIZATION>`_
subsection. Note that in the present implementation this only works in
C\ :sub:`1`\  symmetry.


\*LOCALIZATION -- Molecular orbital localization
================================================

Performs molecular orbitals localization using the Pipek-Mezey
criterion. The implementation uses an exponential parametrization and a
quasi-Newton minimization method instead of the Edminston-Ruedenberg
algorithm. For the moment, it only works with the *C*\ :sub:`1`\ 
symmetry and for the occupied orbitals.

.THRESH
-------

Convergence threshold.

*Default:*

::

    .THRESH
     1.0D-6

.MAXITR
-------

Maximum number of iterations.

*Default:*

::

    .MAXITR
     100

.PRJLOC
-------

Localization from projection analysis instead of Mulliken analysis.

::

    .PRJLOC
     1

*Default:*

::

    .PRJLOC
     0

.OWNBAS
-------

See `.OWNBAS </dirac/index.php/Manual:ANALYZE:PROJECTION#.OWNBAS>`_
under `\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_.

.VECPRJ
-------

See `.VECPRJ </dirac/index.php/Manual:ANALYZE:PROJECTION#.VECPRJ>`_
under `\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_.

.VECREF
-------

See `.VECREF </dirac/index.php/Manual:ANALYZE:PROJECTION#.VECREF>`_
under `\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_.

.HESSUN
-------

Hessian approximated by the unity matrix.

*Default:* Use the diagonal elements of the (diagonally dominant)
Hessian.

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     1



Manual:ANALYZE:VISUAL
=====================

this page moved to...which is also empty


Manual:ANALYZE:RHO1
===================

\*RHO1 -- Print density along lines between atomic centers
==========================================================

For each unique pair (A,B) of centers a formatted file is created,
containing four columns of numbers: The origin is placed at atom A, and
the first column gives the coordinate in Angstroms along the line
between centers A and B. The following columns contain the total, large
and small component density, respectively (in atomic units). The file is
named 'RH',NAMN(ICENTA),NAMN(ICENTB),IDEGB where ICENTB is the name of
symmetry-independent center B and IDEGB the numbering of the dependent
center. Underscores are inserted wherever there are blanks. The
generation of the density requires the coefficient file DFCOEF. Note
that one can in principle obtain the density along any line simply by
introducing ghost centers.

For uranyl one may use the command

::

    pam -incmo -get ``RHU___O____1  RHU___O____2 ...

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.MESH
-----

Step length for grid.

*Default:*

::

    .MESH
     0.01


ANALYZE:Sample inputs
============================

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

