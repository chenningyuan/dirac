15.1 (2016-xx-xx)
-----------------

- added OPENBLAS in search list for blas and lapack libraries
- made keyword .X2Cmmf case insensitive
- BUGfix for spinfree X2Cmmf and linear symmetry.
- Error corrections and updates in isotope properties for the following atoms:
      Br isotope 2: quadrupole moment  .2620 ->  .2615
      Ag isotope 2: magnetic moment    .130563 -> -.130691 (note sign change)
      In isotope 2: quadrupole moment   .790 ->  .799
      Nd magnetic moments of isotopes 4 and 5 were interchanged: -0.065 -> -1.065 and -1.065 -> -0.065
      Gd: quadrupole meoments of isotopes 4 and 5 updated: 1.36 -> 1.35 and 1.30 -> 1.27
      Ho isotope 1: quadrupole moment updated 3.49 -> 3.58
      Lu isotope 2: quadrupole moment updtaed 4.92 -> 4.97
      Hf isotope 1: mass was real*4, not real*8, thus 7 digits instead of 179.9465457D0 (i.e. approx 179.9465)
      Ta isotope 1: quadrupole moment added 0.00 -> 3,17
      Tl isotope 1: nuclear moment 1.63831461D0 -> 1.63821461D0 (typo, error 1.d-4)
      Pb isotope 3: nuclear moment 0.582583D0 -> 0.592583D0 (typo, error 1.d-2)
      Po isotope 1: nuclear moment added: 0.000 -> 0.688


15.0 (2015-12-16)
-----------------

- FanoADC-Stieltjes: Calculation of decay widths of electronic decay processes.
- Better convergence for open-shell calculations, using .OPENFAC 0.5 as default.
  Final orbital energies are recalculated with .OPENFAC 1.0, for IP interpretation.
- Performance improvement for MP2 NO: only transform (G O|G O) integrals instead of (G G|G G).
- Geometry optimization with xyz input.
- Performance improvements for determinant generation in GASCIP.
- Use Kramers conjugation on CI vectors (cuts time for CI in half for ESR doublets).
- ANO-RCC basis: Fixed Carbon basis set (wrong contraction coefficients, see [MOLCAS ANO-RCC](http://www.molcas.org/ANO/).
- ANO-RCC basis: Modified the 3 Th h-functions by replacing them with the 3 Ac h-functions to Th.
  (A mistake was made in the original work when the 3 Th h-functions were made,
  and this has never been redone. They are too diffuse, exponents
  (0.3140887600, 0.1256355100, 0.0502542000) for Th, Z=90, compared to
  (0.7947153600, 0.3149038200, 0.1259615200) for Ac, Z=89, and
  (0.8411791300, 0.3310795400, 0.1324318200) for Pa, Z=91.
  We have selected to just replace the 3 Th h-functions with those from the Ac basis set,
  because the Ac g-functions are quite close tot he Th g-functions, closer than Ac g-functions,
  and therefore differences in results compared to optimized Th h-functions should be minimal.)
  Thanks to Kirk Peterson for pointing out the Th problem on http://daltonforum.org.
- Fixed reading of ANO-RCC and ANO-DK3 basis sets from the included basis set library.
- Update PCMSolver module to its [1.0.4 version](https://github.com/PCMSolver/pcmsolver/releases).
- Do not write Hartree-Fock in output if it is a DFT calculation, write DFT.
- Fix a bunch of problems where the number of arguments of the call does not match the subroutine definition (submitted by Martin Siegert).
- Fixed plotting of the electrostatic potential (electronic part was factor 2 too large).
- Configuration framework uses [Autocmake](http://autocmake.org).


14.1 (2015-07-20)
-----------------

- Added OpenBLAS to default blas search list.
- Fixed misleading error message for .OPEN SHELL input.
- Workaround for CMake warning about nonexistent run_pamadm target.
- Intel compilers xHost flag automatic detection only upon request (-D ENABLE_XHOST_FLAG_DETECTION=ON).
- Fixed minor issues to allow compilation using IBM XL Fortran compiler (thanks to Bob Walkup, IBM).
- Update PCMSolver module to its [1.0.3 version](https://github.com/PCMSolver/pcmsolver/releases).


14.0 (2014-12-12)
-----------------

- DIRAC14 release, see doc/release/release-statement.txt
