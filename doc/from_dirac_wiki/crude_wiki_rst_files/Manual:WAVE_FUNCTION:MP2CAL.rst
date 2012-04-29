Manual:WAVE FUNCTION:MP2CAL
===========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Direct relativistic parallelized closed-shell MP2-correlated
calculation.

Contents
--------

-  `1 \*MP2CAL -- MP2 calculation <#.2AMP2CAL_--_MP2_calculation>`_

   -  `1.1 .OCCUPIED <#.OCCUPIED>`_
   -  `1.2 .VIRTUAL <#.VIRTUAL>`_
   -  `1.3 .INTFLG <#.INTFLG>`_
   -  `1.4 .SCREEN <#.SCREEN>`_
   -  `1.5 .PRINT <#.PRINT>`_
   -  `1.6 .VIRTHR <#.VIRTHR>`_
   -  `1.7 .IJTSK <#.IJTSK>`_
   -  `1.8 .SCLMEM <#.SCLMEM>`_
   -  `1.9 .ORGALG <#.ORGALG>`_
   -  `1.10 Sample input <#Sample_input>`_
   -  `1.11 References <#References>`_

\*MP2CAL -- MP2 calculation
===========================

This section gives directives for the integral-direct closed-shell MP2
calculation. The default algorithm is described
`here <http://folk.uio.no/vebjornb/usit/dirac-mp2/dirac-mp2.php>`_ and
is a modification of the original algorithm published here.
\ :sup:``[1] <#_note-saue:TlAt>`_`\ 

.OCCUPIED
---------

Active, doubly occupied, electronic solutions.

For each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of active orbitals.

*Default:* The occupied electronic solutions.

.VIRTUAL
--------

Active virtual solutions.

For each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of active orbitals.

*Default:* All the unoccupied virtual solutions.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.INTFLG
-------

Specify what two-electron integrals to include in the direct MP2
calculation (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.SCREEN
-------

Screening threshold in the 4-index transformation.

A negative number deactivates screening.

*Default:*

::

    .SCREEN
     1.0D-14

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.VIRTHR
-------

Threshold for neglecting high-lying virtual orbitals. Neglect all
orbitals with energy above this threshold. This option is retained for
backwards compatibility, it is advised to use the more general energy
option in the `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ input.

*Arguments:* Real DMP2\_VIRTHR.

*Default:* No threshold.

.IJTSK
------

Max. number of I (or J) orbitals (occupied orbitals) in each calculation
batch. Usually one batch is sufficient but larger calculations may be
possible by reducing the batch size. Reducing the batch size to half
will approximately reduce the memory requirements to 1/4 and increase
the CPU time by a factor 4.

*Arguments:* Integer IJTSK.

*Default:* All active occupied electronic solutions.

.SCLMEM
-------

Internal memory available for scalar untransformed integrals. Larger
batches will be written to disk.

*Arguments:* Integer MAXSCL.

*Default:* Set by program and is usually sufficient.

.ORGALG
-------

Use original alorithm as described in
\ :sup:``[1] <#_note-saue:TlAt>`_`\ .

*Default:* Use `new
algorithm <http://folk.uio.no/vebjornb/usit/dirac-mp2/dirac-mp2.php>`_.

Sample input
------------

Valence + subvalence \ :sub:`2`\  calculation. The 19 active occupied
Kramers pairs are split in two and the calculation is done in three
batches (1..10+1..10, 11..19+1..10, and 11..19+11..19 for I+J). All
orbitals with energy above 100 hartree (a.u.) are neglected in the
virtual space. No SS integrals are included. Note that the orbitals in
the active space need not be consecutive.

::

    .NELECT
    80 78
    .....
    *MP2CAL
    .OCCUP
    24,32..40
    24,32..39
    .VIRTUAL
    all
    all
    .VIRTHR
    100.0
    .INTFLG
    1 1 0
    .IJTSK
    10

References
----------

#. ↑ \ :sup:``1.0 <#_ref-saue:TlAt_0>`_`\ 
   \ :sup:``1.1 <#_ref-saue:TlAt_1>`_`\  J. K. Lærdahl, T. Saue and K.
   Fægri, *Direct relativistic MP2: Properties of ground state CuF, AgF
   and AuF*, {Theor. Chem. Acc. **97**, *177* (1997) `electronic
   version <http://dx.doi.org/10.1007/s002140050251>`_.

