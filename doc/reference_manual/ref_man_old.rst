Reference Manual
================


This page lists the available keywords for the DIRAC input and molecule
files. For information on how to get started and how to use the program
output see the `Tutorials </dirac/index.php/Tutorials>`_.

DIRAC input file
----------------

-  `\*\*DIRAC </dirac/index.php/Manual:DIRAC>`_ -- Job specification

   -  `\*OPTIMIZE </dirac/index.php/Manual:DIRAC:OPTIMIZE>`_ -- Geometry
      optimization directives

-  `\*\*MOLECULE </dirac/index.php/Manual:Basis:XYZ#.2A.2AMOLECULE>`_ --
   Specify molecule and basis set options

   -  `\*BASIS </dirac/index.php/Manual:Basis:XYZ#.2ABASIS>`_ -- Basis
      set
   -  `\*CENTERS </dirac/index.php/Manual:Basis:XYZ#.2ACENTERS>`_ --
      Physical properties of the nuclei
   -  `\*UNITS </dirac/index.php/Manual:Basis:XYZ#.2AUNITS>`_ -- Units
      of lengths used in the xyz file

-  `\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_ --
   Specification of the Hamiltonian

   -  `\*DFT </dirac/index.php/Manual:HAMILTONIAN:DFT>`_ -- DFT
      directives
   -  `\*AMFI </dirac/index.php/Manual:HAMILTONIAN:AMFI>`_ -- AMFI
      directives
   -  `\*FDE </dirac/index.php/Manual:HAMILTONIAN:FDE>`_ -- FDE
      directives
   -  `Sample
      inputs </dirac/index.php/Manual:HAMILTONIAN:Sample_inputs>`_

-  `\*\*WAVE FUNCTION </dirac/index.php/Manual:WAVE_FUNCTION>`_ -- Get
   the wave function

   -  `\*SCF </dirac/index.php/Manual:WAVE_FUNCTION:SCF>`_ --
      Hartree-Fock/Kohn-Sham calculation
   -  `\*MP2CAL </dirac/index.php/Manual:WAVE_FUNCTION:MP2CAL>`_ -- MP2
      calculation
   -  `\*RESOLVE </dirac/index.php/Manual:WAVE_FUNCTION:RESOLVE>`_ --
      Resolve open-shell states
   -  `\*MVOCAL </dirac/index.php/Manual:WAVE_FUNCTION:MVOCAL>`_ --
      Modified virtual orbitals
   -  `\*KRMCSCF </dirac/index.php/Manual:KR-MCSCF>`_ -- Two- and
      four-component (large-scale) KR-MCSCF module
   -  `\*LUCITA </dirac/index.php/Manual:LUCITA>`_ -- Direct GAS
      (large-scale) CI module

      -  `Sample inputs </dirac/index.php/Manual:LUCITA:Sample_inputs>`_

   -  `\*KRCICALC </dirac/index.php/Manual:KR-CI>`_ -- Direct two- and
      four-component GAS (large-scale) CI module

      -  `KR-CI properties </dirac/index.php/Manual:KR-CI_properties>`_
         -- CI property module

-  `\*\*GRID </dirac/index.php/Manual:GRID>`_ -- Numerical integration
   grid

-  `\*\*ANALYZE </dirac/index.php/Manual:ANALYZE>`_ -- Analyze the wave
   function

   -  `\*PRIVEC </dirac/index.php/Manual:ANALYZE:PRIVEC>`_ -- Print
      vectors
   -  `\*MULPOP </dirac/index.php/Manual:ANALYZE:MULPOP>`_ -- Mulliken
      population analysis
   -  `\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_ --
      Projection analysis
   -  `\*DENSITY </dirac/index.php/Manual:ANALYZE:DENSITY>`_ -- Write
      the density to a formatted file
   -  `\*RHO1 </dirac/index.php/Manual:ANALYZE:RHO1>`_ -- Print density
      along lines between atomic centers
   -  `\*LOCALIZATION </dirac/index.php/Manual:ANALYZE:LOCALIZATION>`_
      -- Molecular orbital localization
   -  `Sample inputs </dirac/index.php/Manual:ANALYZE:Sample_inputs>`_

-  `\*\*VISUAL </dirac/index.php/How_to_visualize_densities_and_properties_using_DIRAC11>`_
   -- Get various unperturbed and perturbed densities for subsequent
   visualization

-  `\*\*PROPERTIES </dirac/index.php/Manual:PROPERTIES>`_ -- Calculation
   of molecular properties

   -  `\*EXPECTATION
      VALUE </dirac/index.php/Manual:PROPERTIES:EXPECTATION_VALUE>`_ --
      Evaluate expectation values
   -  `\*LINEAR
      RESPONSE </dirac/index.php/Manual:PROPERTIES:LINEAR_RESPONSE>`_ --
      Linear response
   -  `\*QUADRATIC
      RESPONSE </dirac/index.php/Manual:PROPERTIES:QUADRATIC_RESPONSE>`_
      -- Quadratic response
   -  `\*EXCITATION
      ENERGIES </dirac/index.php/Manual:PROPERTIES:EXCITATION_ENERGIES>`_
      -- Evaluate excitation energies
   -  `\*NMR </dirac/index.php/Manual:PROPERTIES:NMR>`_ -- Controls
      evaluation of NMR parameters
   -  `\*MOLGRD </dirac/index.php/Manual:PROPERTIES:MOLGRD>`_ --
      Controls evaluation of the molecular gradient
   -  `\*STEX </dirac/index.php/Manual:PROPERTIES:STEX>`_ -- Core
      excitation spectra in the Static Exchange approximation
   -  `Sample
      inputs </dirac/index.php/Manual:PROPERTIES:Sample_inputs>`_

-  `\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ -- Integral
   transformation module

   -  `\*PRPTRA </dirac/index.php/Manual:MOLTRA:PRPTRA>`_ -- Property
      integrals transformation
   -  `Sample inputs </dirac/index.php/Manual:MOLTRA:Sample_inputs>`_

-  `RELCCSD </dirac/index.php/Manual:RELCCSD:New_input>`_ -- Coupled
   cluster module (new input)
-  `RELCCSD </dirac/index.php/Manual:RELCCSD>`_ -- Coupled cluster
   module (old input)

   -  `Sample inputs </dirac/index.php/Manual:RELCCSD:Sample_inputs>`_

-  `DIRRCI </dirac/index.php/Manual:DIRRCI>`_ -- Direct CI module

   -  `&RASORB </dirac/index.php/Manual:DIRRCI:RASORB>`_ -- Specify the
      type of CI and the active space
   -  `&CIROOT </dirac/index.php/Manual:DIRRCI:CIROOT>`_ -- Select the
      state to converge on
   -  `&DIRECT </dirac/index.php/Manual:DIRRCI:DIRECT>`_ -- Convergence
      control
   -  `&OPTIM </dirac/index.php/Manual:DIRRCI:OPTIM>`_ -- Fine tuning of
      the algorithm
   -  `&LEADDET </dirac/index.php/Manual:DIRRCI:LEADDET>`_ -- Analyze
      the CI wave function
   -  `Sample inputs </dirac/index.php/Manual:DIRRCI:Sample_inputs>`_

-  `GOSCIP </dirac/index.php/Manual:GOSCIP>`_ -- COSCI module

   -  `&GOSCIP </dirac/index.php/Manual:GOSCIP:GOSCIP>`_ -- Specify the
      CI space
   -  `&POPANA </dirac/index.php/Manual:GOSCIP:POPANA>`_ -- Analyze the
      wave function
   -  `Sample inputs </dirac/index.php/Manual:GOSCIP:Sample_inputs>`_

-  `RELADC </dirac/index.php/Manual:RELADC>`_ -- Propagator module
   (Green's functions) for ionizations

   -  `Sample inputs </dirac/index.php/Manual:RELADC:Sample_inputs>`_

-  `\*\*GENERAL </dirac/index.php/Manual:GENERAL>`_ -- General input

   -  `\*PARALLEL </dirac/index.php/Manual:GENERAL:PARALLEL>`_ --
      Parallelization directives
   -  `Sample inputs </dirac/index.php/Manual:GENERAL:Sample_inputs>`_

-  `\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_ -- Integral
   directives

   -  `\*READIN </dirac/index.php/Manual:INTEGRALS:READIN>`_ -- Reading
      basis set file
   -  `\*ONEINT </dirac/index.php/Manual:INTEGRALS:ONEINT>`_ --
      One-electron integrals
   -  `\*TWOINT </dirac/index.php/Manual:INTEGRALS:TWOINT>`_ --
      Two-electron integrals
   -  `Sample inputs </dirac/index.php/Manual:INTEGRALS:Sample_inputs>`_

DIRAC molecule file
-------------------

-  `General information </dirac/index.php/Manual:Basis:General>`_
-  `.MOL files </dirac/index.php/Manual:Basis:MOL>`_

   -  `Title and
      symmetry </dirac/index.php/Manual:Basis:Title_and_symmetry>`_
   -  `Atomic
      coordinates </dirac/index.php/Manual:Basis:Atomic_coordinates>`_
   -  `Large component basis
      set </dirac/index.php/Manual:Basis:Large_component_basis_set>`_
   -  `Small component basis
      set </dirac/index.php/Manual:Basis:Small_component_basis_set>`_
   -  `Fitting basis
      set </dirac/index.php/Manual:Basis:Fitting_basis_set>`_

-  `.XYZ files </dirac/index.php/Manual:Basis:XYZ>`_

Special features
----------------

-  `Specification of orbital
   strings </dirac/index.php/Manual:Orbital_strings>`_
-  `Specification of one-electron
   operators </dirac/index.php/Manual:One-electron_operators>`_

Restarting calculations
-----------------------

-  `Restarting DIRAC </dirac/index.php/Restarting>`_

