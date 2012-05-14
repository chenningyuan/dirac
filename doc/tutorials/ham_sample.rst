HAMILTONIAN:Sample inputs
================================

Some examples of the Hamiltonian input.

Add dipole perturbation of strength 0.01 au in the z-direction to the
nonrelativistic Lévy-Leblond Hamiltonian:

::

    **HAMILTONIAN
    .LEVY-LEBLOND
    .OPERATOR
     ZDIPLEN
     COMFACTOR
     0.01

Do a 4-component B3LYP calculation with high precision in the radial
grid:

::

    **DIRAC
    .WAVE FUNCTIONS
    **WAVE FUNCTIONS
    .DHF
    **HAMILTONIAN
    .DFT
     B3LYP
    *DFT
    .RADINT
     1.0D-15

A functional can also be specified manually through the GGAKEY keyword.
In the next input example we use the PW86X and P86C functionals for
exchange and correlation respectively. We also tighten the angular grid
and make the pruning less pronounced:

::

    **DIRAC
    .WAVE FUNCTIONS
    **WAVE FUNCTIONS
    .DHF
    **HAMILTONIAN
    .DFT
     GGAKEY PW86X=1.0 P86C=1.0
    *DFT
    .ANGINT
     53
    .ANGMIN
     15


