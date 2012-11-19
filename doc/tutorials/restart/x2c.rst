

Restarting X2C calculations
===========================

Here we will show how one can run a Hartree--Fock (HF) calculation in
2-component mode with the exact two component (X2C) Hamiltonian and restart a
property calculation at the HF level in a second run without the need to re-do
the X2C transformation.  For the initial SCF step we use the following input
file (copy the content below into the ``scf_x2c.inp`` file)::

  **DIRAC
  .WAVE FUNCTION
  **MOLECULE
  *BASIS
  ! Use the same basis for all atoms
  .DEFAULT
  dyall.av3z
  **HAMILTONIAN
  ! Use the exact two-component Hamiltonian
  ! with atomic mean field spin-orbit
  .X2C
  **INTEGRALS
  *READIN
  ! A limitation in the AMFI code requires an uncontracted basis
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  *END OF

together with the following molecule file (copy it into the ``br2.xyz`` file)::

  2

  Br 0.0 0.0  1.244
  Br 0.0 0.0 -1.244

This molecule has linear symmetry (Dinf,h), and the program will detect this
automatically. In most parts of the program the highest abelian subgroup (D2h)
will be used instead of the full linear symmetry. In this example we use
Dyall's augmented TZ basis dyall.av3z in uncontracted form.

Prior to the HF calculation, the Hamiltonian transformation from 4- to
2-components will be performed and for restart purposes one should therefore
save the files ``X2CMAT`` and ``AOMOMAT`` containing the Hamiltonian
information and the file ``DFCOEF`` containing the (converged) MO-coefficients.

We run the initial SCF step using::

  pam --inp=scf_x2c.inp --mol=br2.xyz --get="X2CMAT AOMOMAT DFCOEF"

The output will be written to the file ``scf_x2c_br2.out``, which can be
monitored during the calculation.

In order to proceed in a second step with a property calculation of,
e.g. the contact density at the Br nucleus at the HF level, we use the following
input file (copy and paste it into the ``property_x2c.inp`` input file)::

  **DIRAC
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
  **INTEGRALS
  *READIN
  ! A limitation in the AMFI code requires an uncontracted basis
  .UNCONTRACT
  **PROPERTIES
  .RHONUC
  *END OF

Restarting from the (converged) MO-coefficients and reading the transformed X2C
Hamiltonian from file is done automatically once the files are present. We
therefore restart the calculation using::

  pam --inp=property_x2c.inp --mol=br2.xyz --put="X2CMAT AOMOMAT DFCOEF"
