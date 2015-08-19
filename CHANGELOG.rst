15.x
____

- Better convergence for open-shell calculations, using .OPENFAC 0.5 as default.
  Final orbital energies are recalculated with .OPENFAC 1.0, for IP interpretation.
- Performance improvement for MP2 NO: only transform (G O|G O) integrals instead of (G G|G G).
- Geometry optimization with xyz input.
- Performance improvements for determinant generation in GASCIP
- Use Kramers conjugation on CI vectors (cuts time for CI in half for ESR doublets)

14.2 (20XX-XX-XX)
-----------------

- Update PCMSolver module to its [1.0.4 version](https://github.com/PCMSolver/pcmsolver/releases).
- Do not write Hartree-Fock in output if it is a DFT calculation
- This will be the next patch.
- Write me ...


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
