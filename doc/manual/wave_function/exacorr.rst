:orphan:
 
.. _exacc:

=====================================================
Highly parallelised relativistic Coupled Cluster Code
=====================================================

starstar(EXACC)

Highly parallelised relativistic Coupled Cluster Code. 
In this release only computations using the X2C Hamiltonian (with either :ref:`HAMILTONIAN_.X2Cmmf` or :ref:`HAMILTONIAN_.X2C`) are possible. 

This code is based on the math libraries `TAL-SH <https://github.com/DmitryLyakh/TAL_SH>`_ and `ExaTENSOR <https://github.com/ORNL-QCI/ExaTENSOR>`_ by Dmitry Lyakh. The tensors are kept in working memeory, sufficent RAM needs to be available. In order to test memory requirements instructions can be found in the ``exacorr_talsh_memory`` and ``exacorr_exatensor_memory`` tests. TAL-SH runs on a single node which has to have enough memory (:ref:`EXACC_.TALSH_BUFF`). In ExaTENSOR the memory is distributed, so each additional node will contribute its memory to the memory pool accessible by the library. Currently, it is recommended to use enough nodes that the tensors fit, but not substantially more. 

In the current release, if the library runs out of memory the code will not stop but enter a blocked state and calculations will not advance. So carefully control the advancement of your calculations, stopping them if they appear to hang.

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

This calculation generates the file CCDENS, which contains the CC ground-state density
matrix in AO basis. In this release, CCDENS is used by the property module to calculate 
ground-state expectation values.

If saved, CCDENS can be used in a property calculation (see :ref:`PROPERTIES_.RDCCDENS`) 
without the need to invoke this module.


*Default:*

::

    Lambda equations are not solved

keyword(NOTRIPLES)

Deactivates computation of triples energy corrections (useful for ExaTENSOR as the current implementation is not efficient)

*Default:*

::

    Triples are done

keyword(CC2)

Performs a CC2 calculation instead of the default CCSD. Currently supported only for energies.

*Default:*

::

    CC2 is not activated


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

Performs a CCD calculation instead of the default CCSD (switch off the contributions of single excitations).

*Default:*

::

    CCDOUBLES is not activated

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

keyword(MP2NO)

Perform a closed shell MP2 calculation and generate the frozen natural orbitals(FNOs) by diagonalization of the virtual-virtual block of the MP2 density matrix, outputting the original occupied orbitals and a set of truncated virtuals in AO basis (the FNOs are transformed into AO basis and saved on file MP2NOs\_AO, which has the same structure as DFCOEF.)

The user must specify a threshold indicating a NO occupation number, and orbitals with occupation below that will not be retained in the truncated AO basis. 

::

    .MP2NO
    1.0d-3

If one would like to use the new set including the FNOs to do higher-level computation like CC, the MP2NOs\_AO should be retrieved from the work directory for the calculati, and used as a standard DFCOEF file in subsequent calculations (currently a run in which FNOs are generated and used in the same post-SCF calculation is not supported).  
    
Apart from MP2NOs\_AO, this option also ouputs the full, untruncated FNOs space (MO basis) in file NOs\_MO. The Fock matrix (FNO basis) is also outputted to file FM\_in\_NatOrb. 

The NOs\_MO is provided so that users wishing to change the truncation thresold above don't have to repeat the same MP2 calculation (If present in the scratch directory, the NOs\_MO will make the code skip the MP2 calculation). 

    
*Default:*

::

    MP2FNO is not activated.


keyword(DONATORB)

Calculate natural occupation numbers and orbitals in AO (quaternion) basis, from a density matrix in the same basis (such as the one generated by the first-order property code, which saved by default in file CCDENS), and store them in file DFNOSAO. 

This option assumes that CCDENS is in the work directory for the calculation.

*Default:*

::

    DONATORB is not activated.


keyword(EXEOM)

Activate the equation of motion (EOMCC) module. This option can be used to obtain excitation energies (EOM-EE), singly ionized (EOM-IP) and electron attached (EOM-EA) states starting from a closed shell reference and at the CCSD level, by solving for the right-hand eigenvalues and eigenvectors of the appropriate similarity-transformed Hamiltonian. Users should specify the number of eom runs.

::

    .EXEOM
    1

*Default:*

::

    EXEOM is not activated.


star(EXEOM)

This menu controls the parameters for the definition of the EOM model used and the number of roots.

Currently the implementation supports the excitation energy (EE), single electron attachment (EA) and detachment (IP) models for CCSD wavefunctions only.

Note that there is no default, if no options are selected no calculations will be performed.


::

   .NROOTS
   4

Specify the number of roots in the calculation.

::

   .FLAVOR
   EE

Specify the flavor of EOM: EE (excitation energy), EA (single electron attachment energies), and IP (single electron detachment energies)



keyword(CCCILR)

Perform coupled-cluster (CCSD) linear response calculation. The current code can run both static and frequency-dependent(complex frequency) response calculation. There are two wave-function models can be used in the calculation: CC-CC and CC-CI(EOM). CC-CC is the standard CC-LR, and CC-CI is the approximation version, which give the same results as EOM-LR.

To run this code, users need to first do AO to MO transformation for the required property by activating the PRPTRA in MOLTRA. In addition, one also need to solve Lamda-equations.

::

    .CCCILR
    1

1 indicates CC-CC type wave-function, and other numbers indicate CC-CI type wave-function. 

keyword(CCCIQR)

Perform coupled-cluster (CCSD) quadratic response calculation. The current code can run both static and frequency-dependent(complex frequency) response calculation. There are two wave-function models can be used in the calculation: CC-CC and CC-CI(EOM). CC-CC is the standard CC-QR, and CC-CI is the approximation version, which give the same results as EOM-QR.

::

    .CCCIQR
    1

1 indicates CC-CC type wave-function, and other numbers indicate CC-CI(EOM) type wave-function. Note: At the moment, CC-CC quadratic response is not supported. 


keyword(LRFRE)

Specify the frequencies (in atomic unit) in the coupled cluster linear response calculation

::

   .LRFRE
   3
   0.122,0.124,0.126
   0.010,0.010,0.010

The first line contains the number of frequencies to be used. The second line indicates the real part of each frequency. The third line indicates the imaginary part of each frequency.


keyword(QR1FRE)

Specify the first frequencies in the coupled cluster quadratic response calculation

::

   .QR1FRE
   1
   0.1
   0.0

The setup for the frequency is the same as LRFRE above.

keyword(QR2FRE)

Specify the second frequencies in the coupled cluster quadratic response calculation

::

   .QR2FRE
   1
   0.15
   0.0

The setup for the frequency is the same as LRFRE above.


keyword(AProp)

Specify the A operator in the response functions.

::

   .AProp
   3
   1,2,3

The first line contains the number of A operators. The second line indicates property number, which is consistent with the number PRPTRA section of MOLTRA.


keyword(BProp)

Specify the B operator in the response functions.

::

   .BProp
   3
   1,2,3

The first line contains the number of B operators. The second line indicates property number, which is consistent with the number PRPTRA section
of MOLTRA.

The code always calculates all possible combinations of A,B for linear response functions.


keyword(CProp)

Specify the C operator in the response functions.

::

   .CProp
   3
   1,2,3

The first line contains the number of C operators. The second line indicates property number, which is consistent with the number PRPTRA section
of MOLTRA.

The code always calculates all possible combinations of A,B,C for quadratic response functions.


keyword(CCTPA)

Evaluate the coupled-cluster two-photon scattering amplitudes (S_(AB,CD)), obtained by combining EOM-EE and quadratic response.

.. warning:: This code is experimental. At the moment, we only can do calculations for the target state, which is non-degenerate, and it may have problems using gcc complier. Only CC-CI (EOM) version works, CC-CC is in development.

*Example:* requesting CC-CI (EOM) two-photon scattering amplitudes between the reference state (ground state) and the target state (number 4 excited state in EOM-EE). The frequency used in the response equastions is a fixed value, which is half of the excitation energy of the target state ::

        .CCTPA
        20
        4
	.AProp
	1
	1
	.BProp
	1
	3
	.CProp
	1
	1
	.DProp
	1
	3 


