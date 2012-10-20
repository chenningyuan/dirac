

Nb atom
=======

DIRAC presently does not perform an open-shell Hartree--Fock in a strict sense,
rather an average-of-configuration calculation, which amounts to the
optimization of the average energy of a set of configurations (or determinants)
generated from the specification of a given number of open shells and their
electron occupations. As an example of such a calculation, showing in
particular tricks to help convergence, we consider the niobium atom.

The ground state configuration of niobium is
[Kr].4d\ :sup:`4`\ .5s\ :sup:`1`\  and so a reasonable input for an
average Hartree--Fock calculation would be::

  **DIRAC
  .WAVE FUNCTION
  .ANALYZE
  **HAMILTONIAN
  .X2C
  **INTEGRALS
  *READIN
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  *SCF
  .CLOSED SHELL
   18 18
  .OPEN SHELL
   2
   4/10,0
   1/2,0
  **ANALYZE
  .MULPOP
  *MULPOP
  .VECPOP
   1..oo
   1..oo
  *END OF

Here we create separate open shells for the four *4d* electrons and the
single *5s* electron to reduce the number of states in the average.

We use the following mol file::

    DIRAC


    C   1
           41.    1
    Nb      0.0000000000        0.0000000000        0.0000000000
    LARGE BASIS dyall.v2z
    FINISH

Unfortunately this calculation does not converge::

  ERROR, SCF not converged, Dirac stops!

  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---

The reason is that the open-shell *4d* and *5s* orbitals mix.  We can see this
in the Mulliken population analysis::

  * Electronic eigenvalue no. 10: -0.2308050163326       (Occupation : f = 0.4000)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb s      L Ag Nb dxx    L Ag Nb dyy    L Ag Nb dzz
  -----------------------------------------------------------------------------------
   alpha    1.0000  |      0.8647        -0.0540        -0.0540         0.2432
   beta     0.0000  |      0.0000         0.0000         0.0000         0.0000

  * Electronic eigenvalue no. 11: -0.2230586181555       (Occupation : f = 0.4000)  m_j= -3/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb dxx    L Ag Nb dyy    L B1gNb dxy    L B2gNb dxz    L B3gNb dyz
  --------------------------------------------------------------------------------------------------
   alpha    0.2650  |      0.0663         0.0663         0.1325         0.0000         0.0000
   beta     0.7350  |      0.0000         0.0000         0.0000         0.3675         0.3675

  * Electronic eigenvalue no. 12: -0.2205588496606       (Occupation : f = 0.4000)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb s      L Ag Nb dxx    L Ag Nb dyy    L Ag Nb dzz    L B2gNb dxz    L B3gNb dyz
  -----------------------------------------------------------------------------------------------------------------
   alpha    0.2178  |      0.0285         0.0481         0.0481         0.0931         0.0000         0.0000
   beta     0.7822  |      0.0000         0.0000         0.0000         0.0000         0.3911         0.3911

  * Electronic eigenvalue no. 13: -0.2186062869816       (Occupation : f = 0.4000)  m_j= -3/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb dxx    L Ag Nb dyy    L B1gNb dxy    L B2gNb dxz    L B3gNb dyz
  --------------------------------------------------------------------------------------------------
   alpha    0.7350  |      0.1837         0.1837         0.3675         0.0000         0.0000
   beta     0.2650  |      0.0000         0.0000         0.0000         0.1325         0.1325

  * Electronic eigenvalue no. 14: -0.2178694489641       (Occupation : f = 0.4000)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb s      L Ag Nb dxx    L Ag Nb dyy    L Ag Nb dzz    L B2gNb dxz    L B3gNb dyz
  -----------------------------------------------------------------------------------------------------------------
   alpha    0.7822  |      0.1068         0.1728         0.1728         0.3298         0.0000         0.0000
   beta     0.2178  |      0.0000         0.0000         0.0000         0.0000         0.1089         0.1089

  * Electronic eigenvalue no. 15: -0.2161776019252       (Occupation : f = 0.5000)  m_j=  5/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb dxx    L Ag Nb dyy    L B1gNb dxy
  --------------------------------------------------------------------
   alpha    1.0000  |      0.2500         0.2500         0.5000
   beta     0.0000  |      0.0000         0.0000         0.0000

We see significant mixing of *4d* and *5s* orbitals for m\ :sub:`j`\  = 1/2.
This mixing is forbidden from atomic symmetry. The problem is that we run the
atom only in linear supersymmetry and so the mixing is introduced by numerical
noise due to the energetic proximity of these orbitals.

An efficient way to handle such a case is to strip off the open-shell
electrons and simply calculate first Nb\ :sup:`5+`\  with the electron
configuration of Kr and to save the converged DFCOEF file::

  **DIRAC
  .WAVE FUNCTION
  .ANALYZE
  **HAMILTONIAN
  .X2C
  **INTEGRALS
  *READIN
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  *SCF
  .CLOSED SHELL
   18 18
  #.OPEN SHELL
  # 2
  # 4/10,0
  # 1/2,0
  **ANALYZE
  .MULPOP
  *MULPOP
  .VECPOP
   1..oo
   1..oo
  *END OF

This calculation converges smoothly and
now the Mulliken population analysis tells us that
the gerade orbitals 10--14 are purely *d*, whereas
gerade orbital 15 is pure *s*. Please verify this.

We will now restart from the converged DFCOEF
and reinsert the open-shell electrons and impose
vector selection by overlap::

  **DIRAC
  .WAVE FUNCTION
  .ANALYZE
  **HAMILTONIAN
  .X2C
  **INTEGRALS
  *READIN
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  *SCF
  .CLOSED SHELL
   18 18
  .OPEN SHELL
   2
   4/10,0
   1/2,0
  .OVLSEL
  .NODYNSEL
  **ANALYZE
  .MULPOP
  *MULPOP
  .VECPOP
   1..oo
   1..oo
  *END OF

This makes sure that the orbitals obtained in the calculation on
Nb\ :sup:`5+`\  stay in place. Again, we keep the DFCOEF file
for further restart.

The final step is again to restart from DFCOEF and this time
we relax the SCF by removing overlap selection::

  **DIRAC
  .WAVE FUNCTION
  .ANALYZE
  **HAMILTONIAN
  .X2C
  **INTEGRALS
  *READIN
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  *SCF
  .CLOSED SHELL
   18 18
  .OPEN SHELL
   2
   4/10,0
   1/2,0
  #.OVLSEL
  #.NODYNSEL
  **ANALYZE
  .MULPOP
  *MULPOP
  .VECPOP
   1..oo
   1..oo
  *END OF

We obtain convergence and our open-shell orbitals are pure *d* and *s*::


  * Electronic eigenvalue no. 10: -0.2557501903224       (Occupation : f = 0.4000)  m_j= -3/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb dxx    L Ag Nb dyy    L B1gNb dxy    L B2gNb dxz    L B3gNb dyz
  --------------------------------------------------------------------------------------------------
   alpha    0.8000  |      0.2000         0.2000         0.4000         0.0000         0.0000
   beta     0.2000  |      0.0000         0.0000         0.0000         0.1000         0.1000

  * Electronic eigenvalue no. 11: -0.2557501876526       (Occupation : f = 0.4000)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb dxx    L Ag Nb dyy    L Ag Nb dzz    L B2gNb dxz    L B3gNb dyz
  --------------------------------------------------------------------------------------------------
   alpha    0.4000  |      0.0667         0.0667         0.2667         0.0000         0.0000
   beta     0.6000  |      0.0000         0.0000         0.0000         0.3000         0.3000

  * Electronic eigenvalue no. 12: -0.2505426815019       (Occupation : f = 0.4000)  m_j=  5/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb dxx    L Ag Nb dyy    L B1gNb dxy
  --------------------------------------------------------------------
   alpha    1.0000  |      0.2500         0.2500         0.5000
   beta     0.0000  |      0.0000         0.0000         0.0000

  * Electronic eigenvalue no. 13: -0.2505426814916       (Occupation : f = 0.4000)  m_j= -3/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb dxx    L Ag Nb dyy    L B1gNb dxy    L B2gNb dxz    L B3gNb dyz
  --------------------------------------------------------------------------------------------------
   alpha    0.2000  |      0.0500         0.0500         0.1000         0.0000         0.0000
   beta     0.8000  |      0.0000         0.0000         0.0000         0.4000         0.4000

  * Electronic eigenvalue no. 14: -0.2505426791823       (Occupation : f = 0.4000)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb dxx    L Ag Nb dyy    L Ag Nb dzz    L B2gNb dxz    L B3gNb dyz
  --------------------------------------------------------------------------------------------------
   alpha    0.6000  |      0.1000         0.1000         0.4000         0.0000         0.0000
   beta     0.4000  |      0.0000         0.0000         0.0000         0.2000         0.2000

  * Electronic eigenvalue no. 15: -0.2044678705701       (Occupation : f = 0.5000)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Nb s
  --------------------------------------
   alpha    1.0000  |      1.0000
   beta     0.0000  |      0.0000
