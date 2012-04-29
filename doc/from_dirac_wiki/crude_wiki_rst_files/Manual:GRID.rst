Manual:GRID
===========

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

The numerical integration scheme uses the Becke partitioning
\ :sup:``[1] <#_note-Becke:1988:jchemphys>`_`\  of the molecular volume
into atomic ones, for which the quadrature is performed in spherical
coordinates. Radial integration is carried out using the scheme proposed
by Lindh *et al.* \ :sup:``[2] <#_note-Lindth:2001:theorchemacc>`_`\ ,
while the angular integration is handled by a set of highly accurate
Lebedev grids. Note that the Becke partitioning scheme generally uses
Slater-Bragg radii for atomic size adjustments. For the heavier
elements, with their more variable oxidation states, this can lead to
errors. If the user invokes the `.ATSIZE <#.ATSIZE>`_ keyword, DIRAC
tries to deduce relative atomic sizes from available densities, if it
can find coefficients.

Contents
--------

-  `1 .RADINT <#.RADINT>`_
-  `2 .ANGINT <#.ANGINT>`_
-  `3 .NOPRUN <#.NOPRUN>`_
-  `4 .ANGMIN <#.ANGMIN>`_
-  `5 .ATSIZE <#.ATSIZE>`_
-  `6 .IMPORT <#.IMPORT>`_
-  `7 .NOZIP <#.NOZIP>`_
-  `8 .4CGRID <#.4CGRID>`_
-  `9 .DEBUG <#.DEBUG>`_
-  `10 .COARSE <#.COARSE>`_
-  `11 .ULTRAFINE <#.ULTRAFINE>`_
-  `12 .INTCHK <#.INTCHK>`_
-  `13 References <#References>`_

.RADINT
-------

Specify the maximum error in the radial integration.

*Default:*

::

    .RADINT
     1.0D-13

.ANGINT
-------

Specify the precision of the Lebedev angular grid. The angular
integration of spherical harmonics will be exact to the L-value given by
the user. By default the grid will be pruned. The highest implemented
value is 64.

*Default:*

::

    .ANGINT
     41

.NOPRUN
-------

Turn off the pruning of the angular grid.

.ANGMIN
-------

Specify the minimum precision of the Lebedev angular grid after pruning.

::

    .ANGMIN
     LEBMIN

Close to the nucleus, the precision of the angular grid will be less
than L-value given by `.ANGINT <#.ANGINT>`_, but it will not be less
than the integer LEBMIN. Note that giving a LEBMIN value greater than or
equal to the L-value given by `.ANGINT <#.ANGINT>`_ is equivalent to
turning the pruning off.

*Default:*

::

    .ANGMIN
     15

.ATSIZE
-------

Generate new estimates for atomic size ratios for use in the Becke
partitioning scheme. Relative sizes for a pair of atoms A and B are
calculated from their contribution to the large component density along
the line connecting A and B.

.IMPORT
-------

::

    *DFT
    .IMPORT
    grid_file

Import previously `.EXPORTed <#.EXPORT>`_ DFT grid from grid\_file for
`BSSE calculations </dirac/index.php/BSSE_calculations_with_DIRAC>`_.

For debugging you can also create your own grid file. The grid file is
formatted - in free format. A grid file for three points could look like
this:

::

    3
    0.1  0.1 0.1    1.0
    0.01 0.2 0.4    0.9
    9.9  9.9 9.0    0.8
    -1

The first integer is the number of points, then come the x-, y-, and
z-coordinate, and the weight for each point. The last line is a negative
integer. Works also in parallel.

You can export the DIRAC grid simply by copying back the file
"numerical\_grid" from the scratch directory.

.NOZIP
------

DIRAC will try to remove redundant grid points when symmetry is present.
Each reflection plane reduces the number of grid points by a factor of
two. With this keyword you can turn this default symmetry grid
compression off.

.4CGRID
-------

Include the small component basis in the generation of the DFT grid even
if you run a 1- or 2-component calculation.

This can be useful if you want to compare to a 4-component calculation
using the identical integration grid.

By default the small component is ignored in the grid generation when
running 1- or 2-component DFT calculations.

.DEBUG
------

Very poor grid - corresponds to

::

    .RADINT
     1.0D-3
    .ANGINT
     10

.COARSE
-------

Coarse grid - corresponds to

::

    .RADINT
     1.0D-11
    .ANGINT
     35

.ULTRAFINE
----------

A better grid than default - corresponds to

::

    .RADINT
     2.0D-15
    .ANGINT
     64

.INTCHK
-------

Test the performance of the grid by computing the overlap matrix
numerically and analyzing the errors. In addition, the error matrix can
be printed.

*Default (no test):*

::

    .INTCHK
     0

*Error analysis:*

::

    .INTCHK
     1

*Print error matrix:*

::

    .INTCHK
     2

References
----------

#. `↑ <#_ref-Becke:1988:jchemphys_0>`_ A. D. Becke, *A multicenter
   numerical integration scheme for polyatomic molecules*, J. Chem.
   Phys. **88**, *2547* (1988) `electronic
   version <http://dx.doi.org/10.1063/1.454033>`_.
#. `↑ <#_ref-Lindth:2001:theorchemacc_0>`_ R. Lindh, P. A. Malmqvist,
   and L. Galgiardi, *Molecular integrals by numerical quadrature. I.
   Radial integration*, Theor. Chem. Acc. **106**, *178* (2001)
   `electronic version <http://dx.doi.org/10.1007/s002140100263>`_.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:GRID <http://wiki.chem.vu.nl/dirac/index.php/Manual:GRID>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:GRID>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:GRID&action=edit>`_
-  `View source </dirac/index.php?title=Manual:GRID&action=edit>`_
-  `History </dirac/index.php?title=Manual:GRID&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:GRID>`_

` </dirac/index.php/Dirac_Program>`_

Navigation
''''''''''

-  `Home </dirac/index.php/Dirac_Program>`_
-  `Documentation </dirac/index.php/Documentation>`_
-  `Reference manual </dirac/index.php/Manual>`_
-  `Discussion forum </dirac/index.php/DiscussionForum>`_
-  `News Archive </dirac/index.php/News_archive>`_
-  `Basis sets </dirac/index.php/Basis_Sets>`_

Search
''''''

Toolbox
'''''''

-  `What links
   here </dirac/index.php/Special:Whatlinkshere/Manual:GRID>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:GRID>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:GRID&printable=yes>`_
-  `Permanent link </dirac/index.php?title=Manual:GRID&oldid=32131>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 23:12, 9 November 2011.
-  This page has been accessed 444 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
