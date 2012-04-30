TDDFT calculations with DIRAC
=============================

Here I would like to give a typical input for a TDDFT calculation of
excitation energies with DIRAC and discuss the important keywords
(colored).

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
     18 12      !zinc atom
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
     1..oo      !remove this line if there is only one Fermion irreducible corepresentation
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
