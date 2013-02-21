

Setting up DIRAC
================

* **Installation:**
  :doc:`Supported platforms <installation/platforms>` |
  :doc:`Linux, Unix, Mac <installation/installation>` |
  :doc:`Windows <installation/windows>`

* **Math libraries:**
  :doc:`Detection and linking <installation/math>` |
  :doc:`MKL environment variables <installation/mkl>`

* **MPI:**
  :doc:`Forwarding environment variables <installation/mpi>`

* **64-bit integer support:**
  :doc:`Do I need it? <installation/int64/faq>` |
  :doc:`MPI <installation/int64/mpi>` |
  :doc:`Math libraries <installation/int64/math>`

* **Testing:**
  :doc:`Running the test set <installation/testing>`

* **Example installations and run scripts:**
  :doc:`Hyperion (SGI Altix ICE 8200; Toulouse) <installation/examples/hyperion>`


Updates and patches
===================

* :doc:`ChangeLog <patches/CHANGELOG>`


First steps in DIRAC territory
==============================

* **Getting started:**
  :doc:`First calculation <tutorials/getting_started>`

* **Molecule input formats:**
  :doc:`mol format <molecule_and_basis/molecule_using_mol>` |
  :doc:`xyz format <molecule_and_basis/molecule_using_xyz>` |
  :doc:`ecp input <molecule_and_basis/molecule_with_ecp>`

* **Basis sets:**
  :doc:`General information <molecule_and_basis/basis>` |
  :doc:`Howto uncontract basis sets <molecule_and_basis/howto_uncontract>`

* **Troubleshooting:**
  :doc:`Known problems <getting_help/known_problems>` |
  :doc:`Memory problems <getting_help/memory>`


Tutorials and walkthrus
=======================

* **SCF convergence:**
  :doc:`Atomic start <tutorials/atomic_start/tutorial>`

* **Restarting and multi-step jobs:**
  :doc:`SCF <tutorials/restart/scf>` |
  :doc:`X2C <tutorials/restart/x2c>` |
  :doc:`Coupled cluster <tutorials/restart/coupled_cluster>` |
  :doc:`Troubleshooting <tutorials/restart/troubleshooting>`

* **2-component Hamiltonians:**
  :doc:`Menu of 2-component Hamiltonians <tutorials/two_component_hamiltonians/hamiltonians>` |
  :doc:`Molecular mean field <tutorials/two_component_hamiltonians/molecular_mean_field>`

* **Relativistic effective core potentials:**
  :doc:`Getting started <tutorials/ecp/first>`

* **Nonrelativistic limit:**
  :doc:`Reproducing nonrelativistic results <tutorials/reproducing_nr_results>`

* **DFT:**
  :doc:`TDDFT <tutorials/dft/tddft>` |
  :doc:`BSSE <tutorials/dft/bsse>` |
  :doc:`CAM functional <tutorials/dft/cam>` |
  :doc:`Troubleshooting <tutorials/dft/troubleshooting>`

* **Long-range WFT/short-range DFT:**
  :doc:`General <tutorials/srDFT/general>`

* **Property calculations:**
  :doc:`Calculation of NMR shieldings using simple magnetic balance <tutorials/simple_magnetic_balance/tutorial>` |
  :doc:`Case study: NMR shielding of UF3Cl3 <tutorials/nmr_uf3cl3/tutorial>` |
  :doc:`An introduction to complex reponse <tutorials/properties/complex_response>`

* **Visualization:**
  :doc:`General overview <tutorials/visual/general/tutorial>` |
  :doc:`Orbital densities <tutorials/visual/orbital_densities/tutorial>`

* **Open-shell SCF:**
  :doc:`Basics <tutorials/open_shell_scf/aoc>` |
  :doc:`Nb atom <tutorials/open_shell_scf/nb>` |
  :doc:`Np atom <tutorials/open_shell_scf/np>`

* **ADC:**
  :doc:`Triatomic molecule <tutorials/reladc/triatomic_molecule>` |
  :doc:`Atom <tutorials/reladc/atom>`

* **ECP:**
  :doc:`First calculation <tutorials/ecp/first>` |
  :doc:`Correlation calculations <tutorials/ecp/combine>`  

* **Outdated tutorials (need update):**
  :doc:`DIRRCI <tutorials/dusty_attic/dirrci>` |
  :doc:`GOSCIP <tutorials/dusty_attic/goscip>` |
  :doc:`LUCITA <tutorials/dusty_attic/lucita>` |
  :doc:`MOLTRA <tutorials/dusty_attic/moltra>` |
  :doc:`MP2 <tutorials/dusty_attic/mp2>` |
  :doc:`Open shells <tutorials/dusty_attic/open_shells>` |
  :doc:`RELCCSD <tutorials/dusty_attic/relccsd>`


Exercises
=========

.. toctree::
   :maxdepth: 1

   exercises/index.rst


Reference manual
================

.. toctree::
   :maxdepth: 1

   manual/index.rst


Developers
==========

* **Releasing:**
  :doc:`Release preparation <programmers/release_preparation>` |
  :doc:`Beta testing <programmers/beta_testing>` |
  :doc:`Where to commit after the release is out <programmers/commit_policy>`

* **Programming:**
  :doc:`Programmers' manual <programmers/programmers_manual>` |
  :doc:`Input parsing <programmers/input_reading>` |
  :doc:`Git submodules <programmers/external_projects>` |
  :doc:`Nightly tests <programmers/cdash>` |
  :doc:`Profiling <programmers/profiling>` |
  :doc:`Debugging <programmers/debugging>` |
  :doc:`Further development <programmers/further_development>` 

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

* :doc:`References <references>`
