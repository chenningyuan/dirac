Manual:Basis:Small component basis set
======================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Small component basis set <#Small_component_basis_set>`_

   -  `1.1 By kinetic balance <#By_kinetic_balance>`_
   -  `1.2 Explicitly typed basis sets <#Explicitly_typed_basis_sets>`_
   -  `1.3 MOLFDIR-type basis sets <#MOLFDIR-type_basis_sets>`_

Small component basis set
=========================

After the large component basis set the small component basis set *can*
be specified. If nothing is specified it is equivalent to specifying
SMALL KINBAL (see below). There are three possibilities for giving the
basis set.

By kinetic balance
------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE BASIS aug-cc-pVDZ
    SMALL KINBAL

The small component basis is generated using kinetic balance. The small
component basis set will *always* be uncontracted, independent of
whether the large component basis set is contracted or not.

Explicitly typed basis sets
---------------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE BASIS aug-cc-pVDZ
    SMALL EXPLICIT  3    1    1    1
    ...

Syntax is analogous to `LARGE
EXPLICIT </dirac/index.php/Manual:Basis:Large_component_basis_set#Explicitly_typed_basis_sets>`_.

MOLFDIR-type basis sets
-----------------------

::

    DIRAC
    Water
    MOLFDIR basis
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE MOLFBAS Oxygen-xyz.bas
    SMALL MOLFBAS Oxygen-xyz.bas

The MOLFDIR basis set file(s) (here: Oxygen-xyz.bas) must be copied to
the scratch area, for example using the pam script:

::

    pam ... -copy "Oxygen-xyz.bas ..."

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:Basis:Small\_component\_basis\_set <http://wiki.chem.vu.nl/dirac/index.php/Manual:Basis:Small_component_basis_set>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:Basis:Small_component_basis_set>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:Basis:Small_component_basis_set&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:Basis:Small_component_basis_set&action=edit>`_
-  `History </dirac/index.php?title=Manual:Basis:Small_component_basis_set&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:Basis:Small_component_basis_set>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:Basis:Small_component_basis_set>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:Basis:Small_component_basis_set>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:Basis:Small_component_basis_set&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:Basis:Small_component_basis_set&oldid=510>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 18:03, 7 February 2008.
-  This page has been accessed 1,136 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
