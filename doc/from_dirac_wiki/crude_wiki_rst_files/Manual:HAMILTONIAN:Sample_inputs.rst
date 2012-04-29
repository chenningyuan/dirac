Manual:HAMILTONIAN:Sample inputs
================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Add dipole perturbation of strength 0.01 au in the z-direction to the
nonrelativistic Lévy-Leblond Hamiltonian:

::

    **HAMILTONIAN
    .LEVY-LEBLOND
    .OPERATOR
     ZDIPLEN
     COMFACTOR
     0.01

Do a 4-component B3LYP calculation with high precision in the radial
grid:

::

    **DIRAC
    .WAVE FUNCTIONS
    **WAVE FUNCTIONS
    .DHF
    **HAMILTONIAN
    .DFT
     B3LYP
    *DFT
    .RADINT
     1.0D-15

A functional can also be specified manually through the GGAKEY keyword.
In the next input example we use the PW86X and P86C functionals for
exchange and correlation respectively. We also tighten the angular grid
and make the pruning less pronounced:

::

    **DIRAC
    .WAVE FUNCTIONS
    **WAVE FUNCTIONS
    .DHF
    **HAMILTONIAN
    .DFT
     GGAKEY PW86X=1.0 P86C=1.0
    *DFT
    .ANGINT
     53
    .ANGMIN
     15

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:HAMILTONIAN:Sample\_inputs <http://wiki.chem.vu.nl/dirac/index.php/Manual:HAMILTONIAN:Sample_inputs>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:HAMILTONIAN:Sample_inputs>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:HAMILTONIAN:Sample_inputs&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:HAMILTONIAN:Sample_inputs&action=edit>`_
-  `History </dirac/index.php?title=Manual:HAMILTONIAN:Sample_inputs&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:HAMILTONIAN:Sample_inputs>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:HAMILTONIAN:Sample_inputs>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:HAMILTONIAN:Sample_inputs>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:HAMILTONIAN:Sample_inputs&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:HAMILTONIAN:Sample_inputs&oldid=11951>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 19:28, 3 September 2008.
-  This page has been accessed 1,228 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
