Manual:ANALYZE:MULPOP
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*MULPOP -- Mulliken population
   analysis <#.2AMULPOP_--_Mulliken_population_analysis>`_

   -  `1.1 .VECPOP <#.VECPOP>`_
   -  `1.2 .NETPOP <#.NETPOP>`_
   -  `1.3 .LABDEF <#.LABDEF>`_
   -  `1.4 .AOLAB <#.AOLAB>`_
   -  `1.5 .INDSML <#.INDSML>`_
   -  `1.6 .PRINT <#.PRINT>`_

\*MULPOP -- Mulliken population analysis
========================================

Mulliken population analysis is performed in AO basis. The analysis is
based on the concept of *labels*. Each basis function is labeled by its
functional type and center (and boson irrep when in SO basis, the
default). The labels are given in the output. A set of primitive labels
may be collected to group labels as specified by the user.

.VECPOP
-------

For each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of orbitals to
analyze.

*Default:* Analyze the occupied electronic solutions.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.NETPOP
-------

Give Mulliken net/overlap populations in addition to gross populations.

.LABDEF
-------

Define labels for use in Mulliken population analysis.

Give first the number of labels to define, then the label name (12
characters) and `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ for each fermion
irrep.

*Example:*

::

    .LABDEF
    4 (total number of lines to follow)
    Re s        1
    Re p        2..4
    Re d        5..10
    Re f        11..20
     ...        21..35

*Pitfall:*

If you decide to redefine labels for the population analysis you have to
reassign **all** of them, otherwise Dirac is unable to continue. The
list of functions with the initial label definitions (*default:
SO-labels*) can be drawn from the DIRAC output starting at the section :

::

                                   GETLAB: SO-labels
                                   -----------------
      * Large components:   35
        1  L Ag AR s        2  L Ag AR dxx      3  L Ag AR dyy      4  L Ag AR dzz      etc.

.AOLAB
------

Base definition of labels on AO basis.

*Default:* Base definition of labels on SO basis.

.INDSML
-------

Use individual small component labels.

*Default:* Gather all small component functions belonging to a given
center and irrep, irrespective of function type.

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:ANALYZE:MULPOP <http://wiki.chem.vu.nl/dirac/index.php/Manual:ANALYZE:MULPOP>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:ANALYZE:MULPOP>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:ANALYZE:MULPOP&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:ANALYZE:MULPOP&action=edit>`_
-  `History </dirac/index.php?title=Manual:ANALYZE:MULPOP&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:ANALYZE:MULPOP>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:ANALYZE:MULPOP>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:ANALYZE:MULPOP>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:ANALYZE:MULPOP&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:ANALYZE:MULPOP&oldid=29921>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 09:16, 1 July 2011.
-  This page has been accessed 1,440 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image1| image:: /dirac/images/thumb/6/65/HJAaJ.png/41px-HJAaJ.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
