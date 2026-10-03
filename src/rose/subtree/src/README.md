Programs to manipulate orbital spaces (occupied and virtual).

Usage 1: Generate intrinsic fragment orbitals and localized orbitals.

References:
B. Senjean, S. Sen, M. Repisky, G. Knizia, L. Visscher, J. Chem. Theory Comput. 17 (2021) 1337–1354.
S. Sen, B. Senjean, L. Visscher, J. Chem. Phys. 158 (2023) 054115.

Executable:
rose.x

Usage 2: Project virtual space on the space spanned by functions intended for polarization and correlation.
B. Senjean, K.G. Dyall, L. Visscher, to be published

Executable:
project_corrfuns.x

Prerequisite files:
HDF5 or KF files containing the supermolecular and fragment orbitals and goemtry information

Produces files:
Modified HDF5 or KF files containing the intrinsic fragment orbitals or the localized orbitals

Remarks:
- The fragment files can be obtained with a different basis set than used in the supermolecular basis (e.g. larger for enhanced precision).
- Orbitals produced by different GTO-type codes can be used together (e.g. fragments by PySCF, main molecule by Psi4).
- Fragments can contain additional atoms to saturate dangling bonds.

# Notes for working within AMS
- Do not introduce dependencies on SCM modules in any of source files except rose_inAMS.
- ROSE no longer reads or writes KF files itself. Data is moved between ADF and ROSE
  by python_scripts/rose_kf_bridge.py, run by hand under amspython and never called by
  the build or the test suite, so the standalone build needs no KF library.
- rose_noAMS.f90 stubs the handful of externals that only the in-AMS build resolves.
  Keep it in step with what rose_inAMS provides.
- HDF5 is needed for the GTO interfaces. This interface can be added from the public repo but does requires some work on getting includes right with foray.

