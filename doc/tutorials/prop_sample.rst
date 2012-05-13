PROPERTIES:Sample inputs
========================

Calculation of NMR parameters:

::

    **PROPERTIES
    .NMR

Calculation of frequency-dependent polarizabilities, skipping rotations
between occupied positive-energy and virtual negative-energy orbitals:

::

    **PROPERTIES
    .POLARI
    *LINEAR RESPONSE
    .B FREQ
     3
     0.0 0.1 0.2
    .SKIPEP


