Manual:ANALYZE:Sample inputs
============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Print all orbitals in fermion irrep 1 and some in irrep 2, large
components only and using AO and not the SO basis:

::

    **ANALYZE
    .PRIVEC
    *PRIVEC
    .VECPRI     # what molecular orbitals to print
     all
    -4..2,4
    .PRICMP     # print only large components
     1 0
    .AOLAB      # use AO basis

Mulliken population analysis defining labels for each atomic shell and
specifying orbital to analyze:

::

    **ANALYZE
    .MULPOP
    *MULPOP
    .LABDEF     # redefinition of labels
    3
    H(s)        1
    H(p)        2..4
    H(d)        5..10
    .VECPOP     # what molecular orbitals to analyze
    1..10
    1,2

Projection analysis using two fragments, represented by coefficient
files DFFRG1 and DFFRG2, respectively, and each calculated in their own
basis:

::

    **ANALYZE
    ..PROJEC
    *PROJECTION
    .OWNBAS     # use the individual basis set of fragments
    .VECPRJ     # what molecular orbitals to analyze
    1,2
    1..3
    .VECREF     # definition of reference orbitals
     2
     DFFRG1
     all
     all
     DFFRG2
     1..4
     1..2

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:ANALYZE:Sample\_inputs <http://wiki.chem.vu.nl/dirac/index.php/Manual:ANALYZE:Sample_inputs>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:ANALYZE:Sample_inputs>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:ANALYZE:Sample_inputs&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:ANALYZE:Sample_inputs&action=edit>`_
-  `History </dirac/index.php?title=Manual:ANALYZE:Sample_inputs&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:ANALYZE:Sample_inputs>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:ANALYZE:Sample_inputs>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:ANALYZE:Sample_inputs>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:ANALYZE:Sample_inputs&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:ANALYZE:Sample_inputs&oldid=3071>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 16:16, 9 March 2008.
-  This page has been accessed 948 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
