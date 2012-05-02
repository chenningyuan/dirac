Manual:Basis:MOL
================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

DIRAC Molecule File Structure
-----------------------------

In this section we provide more information about the MOLECULE format as
is used in Dirac. This file is divided in sections, with the following
being mandatory:

-  `Title and
   symmetry </dirac/index.php/Manual:Basis:Title_and_symmetry>`_
-  `Atomic
   coordinates </dirac/index.php/Manual:Basis:Atomic_coordinates>`_
-  `Large component basis
   set </dirac/index.php/Manual:Basis:Large_component_basis_set>`_

It is possible to define explicitly the small component basis set to be
use with an additional section:

-  `Small component basis
   set </dirac/index.php/Manual:Basis:Small_component_basis_set>`_

but this is however not generally advised, given that the small
component basis set is more conveniently generated via the kinetic
balance procedure.

In the development version(s) of DIRAC it is also possible to define the
so-called fit sets via an additional section:

-  `Fitting basis
   set </dirac/index.php/Manual:Basis:Fitting_basis_set>`_

