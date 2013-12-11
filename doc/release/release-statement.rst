Dear DIRAC-users,

We are happy to announce the release of DIRAC13, the latest version
of the DIRAC code for 2- and 4-component relativistic  molecular
electronic structure calculations. The date has been chosen to mark the
80th anniversary of the `Nobel lecture <http://www.nobelprize.org/nobel_prizes/physics/laureates/1933/dirac-lecture.html>`_ 
given by P. A. M. Dirac Dec 12 1933.

See also our web site: http://diracprogram.org/


New features in DIRAC13
~~~~~~~~~~~~~~~~~~~~~~~

* Magnetizabilities at the Hartree-Fock level with London atomic orbitals (LAOs)
* `Improved output for TDDFT excitation energies <tutorials/dft/tddft.html>`_
* Hueckel start with atomic fragments
* XML output functionality
* Enhancements to frozen density embedding (`FDE <manual/fde.html>`_)
* Polarization propagator for 4-C excitations (ADC2 extended)
* Enhancements to X2C: local spin-free and spin-orbit X2C (`X2C and local X2C <tutorials/two_component_hamiltonians/x2c_mol_loc.html>`_)
* Dyall basis sets redefined to reduce linear dependence and conform to basis archive files, including fixes
* Basis sets for 2s, 2p, 3s, 3p, 3d added to Dyall 2z, Dyall 3z sets and for 2p, 3s, 3p for Dyall 4z basis sets
* Dyall aenz (all-electron) basis sets added, with correlating functions for all shells


Improved compilation and testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Support for Windows 7/8 with GNU MinGW32/64 suite and with native math libraries (`installation instructions <installation/windows.html>`_)
* New test script
* Simplified testing using MPI
* Updated math library detection
* Better support for MKL libraries
* Support for Cray
* Support for MPI runs which do not use mpirun


New defaults
~~~~~~~~~~~~

* The pam script sets::

   MKL_NUM_THREADS=1
   MKL_DYNAMIC="FALSE"
   OMP_NUM_THREADS=1
   OMP_DYNAMIC="FALSE" 

unless these variables are set by the user.


A full list of features of DIRAC13 is found here:
http://diracprogram.org/doku.php?id=features

On http://diracprogram.org/ under "News" you can also find the papers
describing the new developments and some examples of applications.
Our documentation pages
http://diracprogram.org/doc/release-13/
contains a number of new tutorials and is continually updated.

As usual, the code is freely distributed to all academic users;
please find detailed information about licensing and download here:
http://dirac.chem.sdu.dk/DIRAC13/

For proper citation of DIRAC13, please consult
http://www.diracprogram.org/doku.php?id=citation


On behalf of the DIRAC team::

   Lucas Visscher
   Hans Jørgen Aagaard Jensen
   Radovan Bast
   Trond Saue
