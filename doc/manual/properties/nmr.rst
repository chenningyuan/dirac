star(NMR)

This section gives directives for the calculation of NMR parameters.

**Advanced options**

keyword(LONDON)

Activate calculations of magnetic properties (NMR shielding constants
and magnetizabilities with London atomic orbitals.

*Default:* Use conventional atomic orbitals.

keyword(GAUGEO)

Reads in a user defined gauge origin (in bohr). Note that both the gauge
origin and dipole origin are changed in order to have consistency.

*Default for conventional atomic orbitals:*

::

    .GAUGEO
     0.0 0.0 0.0

keyword(GO ANG)

Same as :ref:`NMR_.GAUGEO` but
reads coordinates in angstrom.

keyword(USECM)

Use the center of mass as the gauge origin.

keyword(INTFLG)

Specify what two-electron integrals to include in the two-electron
London contributions to the magnetic field property gradient
(default: :ref:`HAMILTONIAN_.INTFLG` under :ref:`**HAMILTONIAN`).

keyword(NOTWO)

Do not calculate the two-electron London contributions for the magnetic
field property gradient when London atomic orbitals are used.

keyword(NOONEI)

Do not calculate the {H(0),T(B)} reorthonormalization terms for the
magnetic field property gradient when London atomic orbitals are used.

keyword(NOORTH)

Do not calculate the {T(B),h(mK)} reorthonormalization contributions for
the expectation value term when London atomic orbitals are used.

keyword(SYMCON)

Employ the symmetric connection for reorthonormalization terms when
using London atomic orbitals.

*Default:* Use the natural connection.


