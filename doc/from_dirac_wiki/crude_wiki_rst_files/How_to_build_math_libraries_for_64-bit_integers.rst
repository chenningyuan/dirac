How to build math libraries for 64-bit integers
===============================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Before you continue please read the
`FAQ </dirac/index.php/FAQ#I_am_about_to_install_DIRAC_on_a_64-bit_machine.2C_should_I_compile_using_64-bit_integers_or_32-bit_integers.3F>`_
and verify whether you really need 64-bit integers.

It is extra work and perhaps not needed.

If you decide to build DIRAC for 64-bit integers, you have to link to
64-bit integer aware math library/libraries providing both BLAS and
LAPACK.

The following things can happen:

-  Your math library is compiled for 64-bit integers and DIRAC will link
   and run properly.
-  Your math library is not compiled for 64-bit integers and DIRAC will
   not link.
-  Your math library is not compiled for 64-bit integers and DIRAC will
   link but stop at runtime.
-  Your math library is not compiled for 64-bit integers and DIRAC will
   link and not stop at runtime because you deactivated the self-test
   and may produce wrong numbers.

For cases 2-4 you may have to either:

-  Verify linking.
-  Compile your own math library.
-  Use DIRAC's internal math implementation (slow).
-  Go back to 32-bit integers.

MKL
---

MKL provides bindings for 64-bit integers and normally DIRAC should
correctly link if you have specified MATH\_ROOT correctly (see the
installation guide).

Atlas
-----

It is not completely trivial but doable to compile your own Atlas
library with 64-bit integer support.

The difficulty is that you need to compile both LAPACK and Atlas.

To do this we recommend to follow the nice example: `Installing ATLAS
with full LAPACK on
Linux/AMD64 <http://math-atlas.sourceforge.net/atlas_install/atlas_install.html#SECTION00090000000000000000>`_.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/How\_to\_build\_math\_libraries\_for\_64-bit\_integers <http://wiki.chem.vu.nl/dirac/index.php/How_to_build_math_libraries_for_64-bit_integers>`_"

Views
'''''

-  `Article </dirac/index.php/How_to_build_math_libraries_for_64-bit_integers>`_
-  `Discussion </dirac/index.php?title=Talk:How_to_build_math_libraries_for_64-bit_integers&action=edit>`_
-  `View
   source </dirac/index.php?title=How_to_build_math_libraries_for_64-bit_integers&action=edit>`_
-  `History </dirac/index.php?title=How_to_build_math_libraries_for_64-bit_integers&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=How_to_build_math_libraries_for_64-bit_integers>`_

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
   here </dirac/index.php/Special:Whatlinkshere/How_to_build_math_libraries_for_64-bit_integers>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/How_to_build_math_libraries_for_64-bit_integers>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=How_to_build_math_libraries_for_64-bit_integers&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=How_to_build_math_libraries_for_64-bit_integers&oldid=31661>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 15:36, 3 November 2011.
-  This page has been accessed 248 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
