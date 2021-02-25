:orphan:
 
.. _exacc:

=====================================================
Highly parallelised relativistic Coupled Cluster Code
=====================================================

starstar(EXACC)

Highly parallelised relativistic Coupled Cluster Code. 
In this release only computations using the X2C Hamiltonian are possible. 


**Mandatory keywords**
======================

keyword(OCCUPIED)

Defines occupied orbitals. Specification of list or energy range (see :ref:`orbital_strings`).

::

    .OCCUPIED
    energy -1.0 0.0 0.00001


keyword(VIRTUAL)

Defines virtual orbitals. Specification of list or energy range (see :ref:`orbital_strings`).

::

    .VIRTUAL
    20..30

**Optional keywords**
=====================

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

keyword(TCONVERG)

Set convergence criteria (CC iterations, Lambda equations) 

*Default:*

::

    .TCONVERG
     1.0D-9

keyword(NCYCLES)

Maximum number of allowed CC iterations to solve the CC and LAMBDA equations. 

*Default:*

::

    .NCYCLES
     30

keyword(EXATENSOR)

This keyword activates the full multinode EXATENSOR library, which is designed for massively
parallel supercomputers. The additional infrastructure needed for parallel
communication makes this implementation inefficient when used for single node runs.
For such purposes the use of only the TALSH library component is recommended, which is designed for one node but
will make use of GPUs (if available and suitable).

*Default:*

::

    Do not use EXATENSOR


keyword(LAMBDA)

Solve Lambda-equations, needs to be activated in order to compute the one particle density matrix and molecular properties.

*Default:*

::

    .LAMBDA
     FALSE

keyword(NOTRIPLES)

Deactivates computation of triples energy corrections (useful for ExaTENSOR as the current implementation is not efficient)

*Default:*

::

    Triples are done

keyword(MOINT_SCHEME)

Expert option to choose another AO to MO integral transformation scheme. Change at your own risk. 

In TALSH only schemes 3 (default) and 42 (using Cholesky decompostion) are available. 

In ExaTensor schemes 1-4 and 42 are available with 42 using Cholesky decompostion. 
Scheme 4 is default for ExaTensor as it reduces the memory footprint by only keeping part of the AO integrals in memory.
The other methds keep all AO integrals in memeory. 
Scheme 0 prints the memory requirements and attempts to allocate the memory 
without doing the calculation. 

*Default:*

::

    .MOINT_SCHEME
     3

keyword(OCC_BETA)

Can be used to specify a "high-spin" reference determinant with a different number of "barred" occupied orbitals,
than "unbarred" occupied spinors. If .OCC_BETA is specified .OCCUPIED is interpreted as a list of unbarred (alpha) spinors.
NB: alpha and beta are used in a loose sense in relativistic calculations to indicate the (un)barred spinors.

::

    .OCC_BETA
    energy -1.0 0.0 0.00001

keyword(VIR_BETA)

Can be used to specify a different number of "barred" virtual orbitals than "unbarred" occupied spinors. 
If .VIR_BETA is specified .VIRTUAL is interpreted as a list of unbarred (alpha) spinors.
NB: alpha and beta are used in a loose sense in relativistic calculations to indicate the (un)barred spinors.

::

    .VIR_BETA
    20..30

keyword(CCDOUBLES)

Perform CCD calculation (switch off the contributions of single excitations).

*Default:*

::

    .CCDOUBLES
     FALSE

keyword(EXA_BLOCKSIZE)

Expert option: Number to tune the parallel distribution (branching) of the spinor spaces.

*Default:*

::

    .EXA_BLOCKSIZE
     75

keyword(TALSH_BUFF)

Maximum memory (in gigabytes) used in TALSH, aim at about 80% of available memory on your machine.

*Default:*

::

    .TALSH_BUFF
     50

keyword(CHOLESKY)

Threshold to define the accuracy of the Cholesky decomposition (MOINT scheme 42), resulting
in inaccuracies of the computed energy of this order of magnitude (in Hartree units).

*Default:*

::

    .CHOLESKY
     1.0D-9

keyword(LSHIFT)

Expert option: Level shift of orbital energies, ignored for values smaller 0.

*Default:*

::

    .LSHIFT
     0.0D0



