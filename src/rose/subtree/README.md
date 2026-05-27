# ROSE

Reduction of Orbital Space Extent

# INSTALL without the AMS file interface

```
export ROSEDIR=`pwd`
export PATH=$ROSEDIR/bin:$PATH
mkdir build
cd build
cmake ..
make install
ctest

```

# INSTALL with the AMS file interface

```
export ROSEDIR=`pwd`
export AMSBUILDNAME="auto"
export PATH=$AMSHOME/build/$AMSBUILDNAME/lib:$AMSHOME/bin.$AMSBUILDNAME/zlib/lib:$PATH
export PATH=$ROSEDIR/bin:$PATH
mkdir build_withams
cd build_withams
cmake .. -DCMAKE_Fortran_FLAGS="-I$AMSHOME/build/libscm_base.build/$AMSBUILDNAME"
make install
ctest

```

# REQUIREMENTS FOR the AMS file interface

- The ADF file interface works only with AMS2021 and later
- Use the same compiler to compile both ROSE and AMS
- The environment variables AMSHOME and AMSBUILDNAME should be set correctly. The former is standard,
  when you have an AMS installation, the latter currently requires developer access and should be
  consistent with the name you have given your build (default is "auto", change if you have defined something else).


# INSTALL inside AMS (remove files that conflict with AMS and then use regular AMS compilation)

```
cd src
rm rose_main.F90 ams_placeholders.F90 rose_noAMS.F90 hdf5_util/mh5_capi.c hdf5_util/mh5.F90
mv rose_ao_ADF.F90 rose_ao.F90
cd $AMSHOME
$AMSBIN/foray

```

# Notes for working within AMS (developers only)

- Do not introduce dependencies on SCM modules in any of source files listed in Cmakelist.txt
- Cmake is not used so the placeholder and hdf5 files should be removed before compilation (see above)
- HDF5 (needed for the GTO interfaces) should work, but requires some work on getting includes right with foray

# RUN

The automatic testing by ctest attempts to run all four supported codes. These test will of course only work if you have preinstalled these codes, please
consult the documentation of each supported code for the best install procedure.
Further examples are provided in the ```examples``` directory.
They are based on a python interface with DIRAC, PSI4, PYSCF and ADF to generate the MO coefficient files (and for ADF, the overlap matrices as well), needed to construct the Intrinsic Fragment (and Localized Molecular) Orbitals.
For the FCHK extensions, a python script can be used in the repository python_script/ to transform the FCHK file to another FCHK file with only the informations that are read by ROSE.
At the end, ROSE creates new MO coefficient files (.h5 for Gaussian-type orbitals, and .rkf files for ADF which is a Slater-type orbital program).
The .rkf files can be read by ADF for visualization or other purposes.

# CURRENT STATUS OF THE CODE THAT GENERATES IBOs

The code generates a minimal set of IAOs (spanning the occupied space exactly and a few valence virtuals) along with the localized IBOs for the occupied and valence virtuals (with the option to leave out the core orbitals from the localization). Additional hard virtuals are localized using an user-defined energy threshold. The code re-canonicalizes the set of occupied and the virtual(valence + hard) over each fragment after partitioning the orbitals according to a Mulliken population analysis (useful for post-SCF calculations on an embedded fragment).   

The standard file format for Gaussian-type codes is hdf5 (extension .h5) using a data formatting schema taken from DIRAC with slight modifications to allow interfacing PSI4 and PYSCF.

DIRAC INTERFACE (extension .h5)
- Works for real and quaternion spinors.
- Note that if one wants to generate CUBE files, the exact same MOLECULE.XYZ as the one used in ROSE has to be used. The ordering of the atoms in this XYZ file is important or your cube file will be wrong.
- Cannot restart a unrestricted calculation with complex algebra, as DIRAC uses quaternion algebra (Kramers-restricted formalism).

PSI4 INTERFACE (extension .h5)
- Works for real orbitals (non-relativistic and scalar-X2C) and using cartesian functions.
- Atomic fragments can be done with DIRAC or PYSCF to get spherical averaging of the orbital (using the average-of-configurations OpenShell SCF or the fractional occupation option of PYSCF).

PYSCF INTERFACE (extension .h5)
- Works for real orbitals (non-relativistic and scalar-X2C) and using cartesian functions.

GAUSSIAN INTERFACE (extension .fchk)
- Works for real orbitals (non-relativistic and scalar-X2C). Only uncontracted basis set tested.

ADF INTERFACE (extension .rkf)
- Works for non-relativistic and with the zeroth order regular approximation (ZORA) to the Dirac equation.
