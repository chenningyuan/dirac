Manual:DIRAC
============

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*\*DIRAC -- Job
   specification <#.2A.2ADIRAC_--_Job_specification>`_

   -  `1.1 .WAVE FUNCTION <#.WAVE_FUNCTION>`_
   -  `1.2 .ANALYZE <#.ANALYZE>`_
   -  `1.3 .PROPERTIES <#.PROPERTIES>`_
   -  `1.4 .OPTIMIZE <#.OPTIMIZE>`_
   -  `1.5 .4INDEX <#.4INDEX>`_
   -  `1.6 .NO4INDEX <#.NO4INDEX>`_
   -  `1.7 .TITLE <#.TITLE>`_
   -  `1.8 .INPTEST <#.INPTEST>`_
   -  `1.9 .ONLY INTEGRALS <#.ONLY_INTEGRALS>`_

\*\*DIRAC -- Job specification
==============================

By default DIRAC only activates the generation of the basis set (e.g.
kinetic balance conditions for small component functions) and the
one-electron modules. Two-electron integrals over the atomic basis
functions will be calculated when needed by the job modules given below.

We recommend that you use the
`DALTON <http://www.kjemi.uio.no/software/dalton/dalton.html>`_ program
package if you want to e.g. save the two-electron integrals to disk for
another purpose, this is not possible with DIRAC.

.WAVE FUNCTION
--------------

Activates the wave function module(s).

This activates the reading of the `\*\*WAVE
FUNCTION </dirac/index.php/Manual:WAVE_FUNCTION>`_ section, where the
desired wave function type(s) must be specified.

.ANALYZE
--------

Activates the Hartree-Fock/Kohn-Sham analysis module.

This activates the reading of the \*\*ANALYZE section.

.PROPERTIES
-----------

Activates the property module (which will call the integral module for
property integrals).

This activates the reading of the
`\*\*PROPERTIES </dirac/index.php/Manual:PROPERTIES>`_ section.

.OPTIMIZE
---------

Activates the geometry optimization.

This activates the reading of the
`\*OPTIMIZE </dirac/index.php/Manual:DIRAC:OPTIMIZE>`_ subsection.

.4INDEX
-------

Explicitly activates transformation of integrals to molecular orbital
basis.

This activates the reading of the
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ section.

These transformed integrals are currently only used by the
`.RELCCSD </dirac/index.php/Manual:WAVE_FUNCTION#.RELCCSD>`_,
`.RELADC </dirac/index.php/Manual:WAVE_FUNCTION#.RELADC>`_,
`.DIRRCI </dirac/index.php/Manual:WAVE_FUNCTION#.DIRRCI>`_, and
`.LUCITA </dirac/index.php/Manual:WAVE_FUNCTION#.LUCITA>`_ modules, and
if one of these three modules are requested under `\*\*WAVE
FUNCTION </dirac/index.php/Manual:WAVE_FUNCTION>`_, then this flag is
automatically activated unless `.NO4INDEX <#.NO4INDEX>`_ is specified in
this input module.

By default, molecular orbitals with orbital energy between -10 and +20
hartree (a.u.) are included, this can be modified in the
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ section.

.NO4INDEX
---------

Do not automatically activate integral transformation to molecular
orbital basis if any of
`.RELCCSD </dirac/index.php/Manual:WAVE_FUNCTION#.RELCCSD>`_,
`.RELADC </dirac/index.php/Manual:WAVE_FUNCTION#.RELADC>`_,
`.DIRRCI </dirac/index.php/Manual:WAVE_FUNCTION#.DIRRCI>`_,
`.LUCITA </dirac/index.php/Manual:WAVE_FUNCTION#.LUCITA>`_, and `.KR
CI </dirac/index.php/Manual:WAVE_FUNCTION#.KR_CI>`_ are requested under
`\*\*WAVE FUNCTION </dirac/index.php/Manual:WAVE_FUNCTION>`_.

This keyword is utilized when repeating correlated CC or CI calculations
(with different parameters for instance) based on saved files after the
integral transformation.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.TITLE
------

Title line (max. 50 characters).

*Default:*

::

    .TITLE
    DIRAC: No title specified !!!                     

.INPTEST
--------

Input test - no job modules are called, only verification of DIRAC input
files.

It is often useful to start a new set of calculations with an input test
in order to check that input file processing is correct before
submitting your (long-term run) job.

.ONLY INTEGRALS
---------------

Stop after the calculation of the one-electron integrals for the
Hamiltonian and the one-electron integrals specified under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_. The integrals are
written to disk.

