Manual:Basis:General
====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

General information
===================

Choosing input type
-------------------

Presently Dirac can handle two types of files with molecular
coordinates:

-  `.MOL </dirac/index.php/Manual:Basis:MOL>`_
-  `.XYZ </dirac/index.php/Manual:Basis:XYZ>`_

Which file is present, is determined on the basis of whether a dot is
present in the name of the submitted molecule file. If no dot is
present, it is assumed to be a .MOL file, if the name ends with .MOL or
.XYZ dirac will detect this and handle the molecule file accordingly. If
an other extension is detected Dirac will exit (as more filetypes may be
included later). In a .MOL file all the information of the molecule,
like basis set, coordinates, charges etc. is present in the molecule
file. In a .XYZ file only the coordinates of the atoms are given. The
setup of the calculation is given in the .INP file.

Choosing basis sets
-------------------

Apart from the practical aspects involved in setting up a MOLECULE file
for Dirac, one must keep in mind also the more fundamental issue of
choosing a suitable basis set for a given calculation.

The development of basis sets suitable to use in relativistic
calculations reflects the relatively lateness of the field's
development. Because the more consistent efforts in method development
started at about the mid 1980's, it wasn't until well into the late
1990's that the pioneering works of the early and mid 1990's were
substantially complemented and improved upon.

This situations has dramatically improved in recent years, notably with
the work of K. G. Dyall, and Dirac users are strongly advised to use
Dyall's basis sets whenever they are available. These sets follow
roughly the "correlation-consisten" philosophy introduced by Dunning and
coworkers \ :sup:``[1] <#_note-dunninc:cc:1>`_`\ , so they already
contain polarization functions.

At the time of writing, these are avaiable from Ga onwards, and a
detailed description can be found at the dirac website, under the "Basis
Set Repository"
link`[1] <http://dirac.chem.sdu.dk/basisarchives/dyall/index.html>`_.
For convenience, Dirac collects all of Dyall's available basis sets
(those published and areadly in distributable form, as the ones listed
in the repository below, as well as those in preparation) in basis set
files. One must keep in mind that currently these files support only the
use of fully uncontracted basis sets.

The following files are present in Dirac's distribution:

dyall.v2z contains the double-zeta basis dyall.v3z contains the
triple-zeta basis dyall.v4z contains the quadruple-zeta basis dyall.av2z
contains the double-zeta basis (as in dyall.v2z) plus diffuse functions
(p-block only) dyall.av3z contains the triple-zeta basis (as in
dyall.v3z) plus diffuse functions (p-block only) dyall.av4z contains the
quadruple-zeta basis (as in dyall.v4z) plus diffuse functions (p-block
only)

 For light elements (say up to Argon) it is advisable to continue using
standard non-relativistic basis set, such as the correlation-consistent
sets of Dunning and coworkers. It is advisable that, in order to have a
balanced description when light and heavy elements are present, that one
uses either contracted or uncontracted sets thoughout.

Apart from Dyall's sets, one can choose several different basis sets
based upon geometric progressions of exponents. One such set is that of
K. Faegri, also available in the basis set library, but with the
drawback that the user needs to extended it by adding polarization
functions.

The procedure used in Dirac is as follows: first, the basis set for the
large component will be read in from the the MOLECULE file; then the
small component basis set is generated automatically via the kinetic
balance prescription, but can (by experienced users) also be specified
explicitly. Try this only if you have studied the theory and know what
the criteria for small component basis sets are, there are many pitfalls
that one needs to be aware of!

--------------

#. `↑ <#_ref-dunninc:cc:1_0>`_ T. H. Dunning, *bla*, J. Chem. Phys.
   **106**, *0* (1989) `electronic
   version <http://dx.doi.org/10.1007/s002140100263>`_.

