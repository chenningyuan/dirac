

12.8
----

- Bugfix for screening in MP2. Now you can use
  screening in conjunction with batched processing (for large systems) by
  specifying .IJTSK in input.


12.7
----

- Bugfix for permanent dipole moment calculations with KRCI. Missing nuclear contribution term 
  have been added and the correct sign for the dipole integrals (length gauge) has been introduced.
- Bugfix for MP2-NO calculations. correct writing of 2-electron file MDCINT.
- Workaround for memory allocator segmentation fault with gfortran can now be triggered with:
  ./setup -D CPP="-DGFORTRAN_ALLOCATOR_WORKAROUND"
- improved final state analysis for MCSCF wave functions: print natural orbital occupation numbers by default
- Bugfix for parallel CC/MP2 calculations using the RELCCSD module and the new input style (the old Namelist input was
  fine); thanks to Kirk Peterson for pointing it out.
- DSO integrals can be limited by the .SELECT keyword. Previously DSO integral code
  ignored .SELECT.
- Extended addressing of DSO integrals up to 333 atom centers.
- Compilation fix for 64bit integer Intel MKL.
- DIRAC compiles using Portland group compilers.
- Introduced block for any structure optimization runs using the xyz format.
- Removed --force option from setup script.


12.6
----

- Add basis_ecp to "make install" target.
- Make it possible to compile without XCFun.
- Enable creation of static executables with older linkers.
- Add --replace option to the pam script.
- Introduce keyword .THRECI in KR-CI input to allow the user to change the default CI energy convergence criterion,
  present default::

    .THRECI
     1.0d-10

- Bugfix for MPI interface on Cray/MPICH2 systems (thanks to Juho for pointing out the problem).


12.5
----

- IMPORTANT fix for two-component calculations based on the X2C/BSS Hamiltonian with two-electron spin-orbit corrections from AMFI.
  All DIRAC versions since DIRAC10 were affected.
  All two-component calculation based on a molecular input (.mol/.xyz) where atoms of the same nuclear charge were
  assigned different atom types in DIRAC were errorneous due to a multiple addition of two-electron spin-orbit corrections.
  Consider the following example h2.xyz input::

    2
    H2 molecule
    Ha 0.0 0.0 0.0
    Hb 0.0 0.0 0.8

  and the SCF input in C1 symmetry for DIRAC (scf_x2c.inp)::

    **DIRAC
    .WAVE FUNCTION
    **INTEGRALS
    *READIN
    .UNCONTRACTED
    **HAMILTONIAN
    .X2C
    **WAVE FUNCTION
    .SCF
    **MOLECULE
    *SYMMETRY
    .NOSYM
    *BASIS
    .DEFAULT
     6-31G-star-star
    *CENTERS
    .NUCLEUS
    Ha 1.0d0
    .NUCLEUS
    Hb 1.0d0
    **END OF

  In this case DIRAC will assign each hydrogen atom Ha/Hb a different atom type (cat from the output)::

    ...

    Atoms and basis sets
    --------------------

    Number of atom types:     2
    Total number of atoms:    2

    label    atoms   charge   prim    cont     basis
    ----------------------------------------------------------------------
    Ha          1       1       7       7      L  - [4s1p|4s1p]
                                0       0      S  - No small component basis functions attached
    Hb          1       1       7       7      L  - [4s1p|4s1p]
                                0       0      S  - No small component basis functions attached
    ----------------------------------------------------------------------
                               14      14      L  - large components
                                0       0      S  - small components
    ----------------------------------------------------------------------
    total:      2       2      14      14

  In this case 2-electron spin-orbit corrections would have been calculated for both Ha and Hb (which is correct as such)
  but since the two centers have the same nuclear charge the old code (DIRAC version < 12.5) would have added
  2-electron spin-orbit corrections, calculated for Ha, also to Hb and vice versa. As a result, spin-orbit corrections
  would have been added twice rather than once.
  Advice: if you are uncertain whether this applies in your case,
  please have a look at your output file(s) and check whether atoms of the same nuclear
  charge (e.g. charge == 10) had been assigned different atom types (due to symmetry or other reasons)
  which would list them as a separate line (X, Y) in the output table::

    label    atoms   charge   prim    cont     basis
    ----------------------------------------------------------------------
    X         ...       10     ...    ...       ...
    Y         ...       10     ...    ...       ...

- Allow pam --get and --put with rename (example: pam --get DFCOEF=dz.coef).
- Correct pam exit message when dirac.x crashes.


12.4
----

- Fix for sequential KRCI property runs - file label was not properly dimensioned.
- Remove unsupported options from CMakeLists.txt and setup.


12.3
----

- XCFun will not compile hidden cpp files.
- Code stops if user gives orbital occupation which is inconsistent with NFSYM.
- runtest --mpi=1 runs tests on one processor with MPI interface.
- Fixed wrong program stop in open-shell geometry optimization.


12.2
----

- Fixed memory allocations in dirgp.F which now allows pre- and post-SCF orbital reordering
  as well as simple magnetic balance for larger basis sets.


12.1
----

- Improved documentation.
- Added test "kallay_mrcc_interface".
- Made it possible to run MPI binary on one processor.
- Compilation fix for older ifort (10.1) together with gcc 4.1.
- Fix for too slow MCSCF convergence.
