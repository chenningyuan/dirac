TDDFT calculations with DIRAC
=============================

Here I would like to give a typical input for a TDDFT calculation of
excitation energies with DIRAC and discuss the important keywords
(colored).

Zn atom
-------

This is a TDDFT calculation of the Zn atom with the PBE functional.

To obtain good TDDFT excitation energies it is *essential* to include
the spin density contribution (spin magnetization contribution). This is
done with the .SDFT keyword. Without spin density, triplet excitation
energies correspond to KS orbital energy differences.

We suppress orbital rotations between positive and negative energy
orbitals with .SKIPEP.

Also include the .TRIPLET keyword. This keyword is without effect (but
also without harm) when spin-orbit coupling is turned on. However,
without spin-orbit coupling (nonrelativistic or spin-free calculations)
this keyword allows the spin density contribution in spin-forbidden
transitions. Without spin-orbit coupling and without .TRIPLET, triplet
energies again correspond to KS orbital energy differences. This keyword
does *not* mean that only triplet states are solved for.

::

    **DIRAC
    .WAVE FUNCTION
    .PROPERTIES
    .ANALYZE
    **HAMILTONIAN
    .LVCORR
    .DFT
     PBE
    *DFT
    .SDFT
    **INTEGRALS
    *READIN
    .UNCONTRACT
    **WAVE FUNCTION
    .DHF
    *DHFCAL
    .EVCCNV
     1.0D-9 1.0D-6
    .CLOSED
     18 12      !zinc atom
    **PROPERTIES
    *EXCITATION ENERGIES
    .TRIPLET
    .SKIPEP
    .ANALYZE
    .EXCITATION ENERGIES
     3 4
    .EXCITATION ENERGIES
     5 4
    .EXCITATION ENERGIES
     8 4
    .OPERATOR
     XDIPLEN
    .OPERATOR
     YDIPLEN
    .OPERATOR
     ZDIPLEN
    **ANALYZE
    .MULPOP
    *MULPOP
    .VECPOP
     1..oo
     1..oo      !remove this line if there is only one Fermion irreducible corepresentation
    *END OF

With

::

    .EXCITATION ENERGIES
     3 4
    .EXCITATION ENERGIES
     5 4
    .EXCITATION ENERGIES
     8 4

we ask for the four lowest excitation energies in the irreps 3, 5, and
8.

In this example (Zn, *D*\ :sub:`2*h*`\ ) these are the ungerade irreps
B2u, B1u, and Au.

Finally, it is a good idea to always include

::

    .OPERATOR
     XDIPLEN
    .OPERATOR
     YDIPLEN
    .OPERATOR
     ZDIPLEN

to get transition dipole moments and oscillator strengths (if symmetry
allowed).

Br2 molecule
-------------

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

