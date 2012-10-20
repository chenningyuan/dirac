

Setting up DIRAC
================

* **Installation:**
  :doc:`Supported platforms <installation/platforms>` |
  :doc:`Linux, Unix, Mac <installation/installation>` |
  :doc:`Windows <installation/windows>`

* **Testing:**
  :doc:`Running the test set <installation/tests>` |
  :doc:`How it works <programmers/testing>`

* **MPI:**
  :doc:`64bit integer support <installation/mpi_int64>`

* **Math libraries:**
  :doc:`64bit integer support <installation/math_int64>` |
  :doc:`MKL environment variables <installation/mkl>`


First steps in DIRAC territory
==============================

* **Getting started:**
  :doc:`First calculation <tutorials/getting_started>`

* **Molecule input formats:**
  :doc:`mol format <molecule_and_basis/molecule_using_mol>` |
  :doc:`xyz format <molecule_and_basis/molecule_using_xyz>`

* **Basis sets:**
  :doc:`General information <molecule_and_basis/basis>` |
  :doc:`Howto uncontract basis sets <molecule_and_basis/howto_uncontract>`


Tutorials and walkthrus
=======================

* **SCF convergence:**
  :doc:`Atomic start <tutorials/atomic_start>`

* **2- and 4-component Hamiltonians:**
  :doc:`Sample inputs <tutorials/ham_sample>` |
  :doc:`2-component Hamiltonians <tutorials/ham2c>` |
  :doc:`AMFI <tutorials/amfi>`

* **Restart:**
  :doc:`Restarting and multi-step jobs <tutorials/restart>`

* **Visualization:**
  :doc:`Vector fields <tutorials/visualize>` |
  :doc:`Orbital densities <tutorials/plotting_orbital_densities>`

* **Analysis:**
  :doc:`Mulliken and projection analysis <tutorials/analyze_sample>`

* **Modeling the solvent:**
  :doc:`Polarizable continuum model <tutorials/pcm>`

* **DFT:**
  :doc:`TDDFT <tutorials/tddft>` |
  :doc:`BSSE <tutorials/bsse>` |
  :doc:`CAM functional <tutorials/dft_cam>`

* **Molecular properties:**
  :doc:`Example calculation <tutorials/prop_sample>`

* **Correlated calculations:**
  :doc:`MP2 <tutorials/mp2cal_sample>` |
  :doc:`MOLTRA <tutorials/moltra_sample>` |
  :doc:`LUCITA <tutorials/lucita_sample>` |
  :doc:`RELCCSD <tutorials/relccsd_sample>` |
  :doc:`RELADC <tutorials/reladc_sample>`

* **Open-shell systems:**
  :doc:`DIRRCI <tutorials/dirrci_sample>` |
  :doc:`GOSCIP <tutorials/goscip_sample>` |
  :doc:`Open-shell examples <tutorials/open_shells>` |
  :doc:`Nb atom <tutorials/open_shell_atoms/nb>` |
  :doc:`Np atom <tutorials/open_shell_atoms/np>` |
  :doc:`SCF convergence <tutorials/scf_conv>`

* **Troubleshooting:**
  :doc:`Memory problems <getting_help/memory>` |
  :doc:`Known problems <getting_help/known_problems>` |
  :doc:`Known bugs <getting_help/known_bugs>` |
  :doc:`Frequently asked questions <tutorials/faq>`


Reference manual
================

* **Input sections:**
  :doc:`DIRAC <reference_manual/dirac>` |
  :doc:`GENERAL <reference_manual/general>` |
  :doc:`INTEGRALS <reference_manual/integrals>` |
  :doc:`HAMILTONIAN <reference_manual/hamiltonian>` |
  :doc:`WAVE FUNCTION <reference_manual/wavefunction>` |
  :doc:`GOSCIP <reference_manual/goscip>` |
  :doc:`ANALYZE <reference_manual/analyze>` |
  :doc:`PROPERTIES <reference_manual/properties>` |
  :doc:`MOLTRA <reference_manual/moltra>` |
  :doc:`DIRRCI <reference_manual/dirrci>` |
  :doc:`RELCCSD (old) <reference_manual/relccsd_old>` |
  :doc:`RELCCSD <reference_manual/relccsd_new>` |
  :doc:`RELADC <reference_manual/reladc>` |
  :doc:`GRID <reference_manual/grid>`

* **Operators and strings:**
  :doc:`One-electron operators <reference_manual/oneel_oper>` |
  :doc:`Orbital strings <reference_manual/orb_strings>`


Developers
==========

* **Releasing:**
  :doc:`Release preparation <programmers/release_preparation>` |
  :doc:`Beta testing <programmers/beta_testing>`

* **Programming:**
  :doc:`Programmers' manual <programmers/programmers_manual>` |
  :doc:`Input parsing <programmers/input_reading>` |
  :doc:`Git submodules <programmers/external_projects>` |
  :doc:`Nightly tests <programmers/cdash>` |
  :doc:`Profiling <programmers/profiling>` |
  :doc:`Debugging <programmers/debugging>`

* **Notes:**
  :doc:`Standalone programs <programmers/standalones>` |
  :doc:`Shell programming <programmers/shell_prog>` |
  :doc:`DFCOEF <programmers/dfcoef>` |
  :doc:`Screening <programmers/screening>` |
  :doc:`64bit integers <programmers/int8>` |
  :doc:`XML <programmers/xml>` |
  :doc:`Static linking <programmers/static_linking>` |
  :doc:`How this documentation works <programmers/documentation_howto>`


Bibliography
============

* :doc:`References <bib>`


Rest - move it to the right place
=================================

.. toctree::
   :maxdepth: 1

   tutorials/reproducing_nr_results.rst
   tutorials/general_scf_sample.rst


Does not belong here - will be moved elsewhere else
===================================================

.. toctree::
   :maxdepth: 1

   welcome/features.rst
   welcome/citing.rst
   installation/licensing.rst
   programmers/grasp.rst
