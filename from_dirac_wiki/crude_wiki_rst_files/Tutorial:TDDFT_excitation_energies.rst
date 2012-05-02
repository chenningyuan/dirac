Tutorial:TDDFT excitation energies
==================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

This examples shows how one can compute excitation energies with DFT and
the exact two component approximation (X2C). The input file
(`b3lyp.inp <#b3lyp.inp>`_) specifies the boson irreps of the
excitations to be considered. Since this depends on the molecule
symmetry it has to be adapted for each particular molecule.

The Br2 molecule (`Br2.xyz <#Br2.xyz>`_) has linear symmetry (D\_oo,h),
and the program will detect this automatically. In most parts of the
program the highest abelian subgroup (D2h) will be used instead of the
full linear symmetry.

In this example we use Dyall's augmented TZ basis dyall.av3z in
uncontracted form.

Run the example using:

::

    pam --mol=Br2.xyz --inp=b3lyp.inp

Output will be written to the file b3lyp\_Br2.out, which can be
monitored during the calculation. Restart information will be saved to
b3lyp\_Br2.tgz.

b3lyp.inp
'''''''''

::

    **DIRAC
    .WAVE FUNCTION
    .PROPERTIES
    **MOLECULE
    *BASIS
    ! Use the same basis for all atoms
    .DEFAULT
    dyall.av3z
    **HAMILTONIAN
    ! Use the exact two-component Hamiltonian
    ! with atomic mean field spin-orbit
    .X2C
    ! Use the B3LYP functional with default parameters
    .DFT
     B3LYP
    **INTEGRALS
    *READIN
    ! A limitation in the AMFI code requires an uncontracted basis
    .UNCONTRACT
    **WAVE FUNCTION
    .SCF
    **PROPERTIES
    ! Run excitation energies
    *EXCITATIONS
    ! Compute transition moments in the length gauge
    .A OPERATOR
     XDIPLEN
    .A OPERATOR
     YDIPLEN
    .A OPERATOR
     ZDIPLEN
    .EXCITATIONS
     2 8
    !  8 excitations in boson irrep 2 (here B3u)
    .EXCITATIONS
     5 8
    !  8 excitations in boson irrep 5 (here B1u)
    ! Analyze orbital contributions to each excitation
    .ANALYZE
    *END OF

Br2.xyz
'''''''

::

    2
    Br 0.0 0.0  1.244
    Br 0.0 0.0 -1.244

