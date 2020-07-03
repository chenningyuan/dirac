:orphan:
 

============================================
Highly parallelised relativistic Coupled Cluster Code
============================================

starstar(EXACC)

Highly parallelised relativistic Coupled Cluster Code. 
So far only computations using the X2C Hamiltonian are possible. 


**Mandatory keywords**
======================

keyword(OCCUPIED)

Defines occupied orbitals. List or energy range is possible.

::

    .OCCUPIED
    energy -1.0 0.0 0.00001


keyword(VIRTUAL)

Defines virtual orbitals. List or energy range is possible.

::

    .VIRTUAL
    20..30

**Optional keywords**
=====================

keyword(EXATENSOR)

This keyword activates the multinode EXATENSOR library. 
Otherwise TALSH isused, which is designed for one node and makes use of GPUs.

*Default:*

::

    Do not use EXATENSOR


keyword(LAMBDA)

Solve Lambda-equations, needs to be activated for properties

*Default:*

::

    .LAMBDA
     FALSE

keyword(NOTRIPLES)

keyword to deactivate triples (recommended for ExaTENSOR)

*Default:*

::

    Triples are done

keyword(MOINT_SCHEME)

AO to MO integral transformation scheme. In TALSH only 3 and 42 are available. 
In Exatensor 1-4 and 42 are available. 42 uses cholesky decompostion. 
Scheme 4 reduces the memory by only treating a part of the AO integrals. 
The other methds keep all AO integrals in memeory. 
Scheme 0 prints the memory requirements and attempts to allocate the memory 
without doing the calculation. 

*Default:*

::

    .MOINT_SCHEME
     3

keyword(OCC_BETA)

Defines occupied orbitals for beta spin, alpha spin defined by .OCCUPIED. List or energy range is possible.

::

    .OCC_BETA
    energy -1.0 0.0 0.00001

keyword(VIR_BETA)

Defines virtual orbitals for beta spin, alpha spin defined by .VIRTUAL. List or energy range is possible.

::

    .VIR_BETA
    20..30

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

keyword(CCDOUBLES)

perform CCD calculation

*Default:*

::

    .CCDOUBLES
     FALSE

keyword(EXA_BLOCKSIZE)

Number to branch the tensors, should be smaller than nocc

*Default:*

::

    .EXA_BLOCKSIZE
     75

keyword(TALSH_BUFF)

Maximum memory used in TALSH / for density matrix, should be about 80 % of avialable memory

*Default:*

::

    .TALSH_BUFF
     50

keyword(NCYCLES)

Maximum number of iterations in CC and LAMBDA. 

*Default:*

::

    .NCYCLES
     30

keyword(TCONVERG)

Set convergence criterria (CC iterations, Lambda equations) 

*Default:*

::

    .TCONVERG
     1.0D-9

keyword(CHOLESKY)

Threshold for cholesky decomposition (MOINT scheme 42)

*Default:*

::

    .CHOLESKY
     1.0D-9

keyword(LSHIFT)

Level shift, ignored for values smaller 0.

*Default:*

::

    .LSHIFT
     0.0D0



