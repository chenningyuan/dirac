Manual:RELCCSD:CCFSPC
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCFSPC -- Fock space module <#.26CCFSPC_--_Fock_space_module>`_

   -  `1.1 DOEXC <#DOEXC>`_
   -  `1.2 DOEA <#DOEA>`_
   -  `1.3 DOIE <#DOIE>`_
   -  `1.4 DOEA2 <#DOEA2>`_
   -  `1.5 DOIE2 <#DOIE2>`_
   -  `1.6 NACTP <#NACTP>`_
   -  `1.7 NACTH <#NACTH>`_
   -  `1.8 NTOL <#NTOL>`_
   -  `1.9 MAXIT <#MAXIT>`_

      -  `1.9.1 MAXITmn <#MAXITmn>`_

   -  `1.10 MAXDIM <#MAXDIM>`_
   -  `1.11 FSSECT <#FSSECT>`_

&CCFSPC -- Fock space module
============================

DOEXC
-----

Calculate excitation energies.

*Default:*

::

    DOEXC=F

DOEA
----

Calculate electron affinities.

*Default:*

::

    DOEA=F

DOIE
----

Calculate ionization energies.

*Default:*

::

    DOIE=F

DOEA2
-----

Calculate electron affinities (add two electrons).

*Default:*

::

    DOEA2=F

DOIE2
-----

Calculate ionization energies (remove two electrons).

*Default:*

::

    DOIE2=F

NACTP
-----

Number of active spinors in which electrons may be put. This variable
determines the active space for the electron affinities (EA)
calculation.

*Arguments:* Integer (NACTP(I),I=1,NFSYM\*2).

*Default:* No default; should always be set when using `DOEA <#DOEA>`_
or `DOEA2 <#DOEA2>`_!

NACTH
-----

Number of active spinors in which holes can be created. This variable
determines the active space for the ionization energies (IE)
calculation.

*Arguments:* Integer (NACTH(I),I=1,NFSYM\*2).

*Default:* No default; should always be set when using `DOIE <#DOIE>`_
or `DOIE2 <#DOIE2>`_!

--------------

An example of the use of the Fock space method concerns the calculation
of the lowest three states of molecular oxygen. These are obtained by
distributing two electrons in the degenerate π\ :sub:`g`\  orbitals. The
X ground state is a triplet that is split by second-order spin-orbit
coupling and spin-spin interactions into a lowest M\ :sub:`S`\ =0
component and two higher M\ :sub:`S`\ =±2 components. This case is also
discussed under the RELCCSD input as an example of the use of the
`NELEC </dirac/index.php/Manual:RELCCSD:RELCCSD#NELEC>`_ keyword. In the
Fock space approach we first perform a calculation on
O\ :sub:`2`\ \ :sup:`+2`\  and calculate the three X states, the
degenerate open shell singlet a and the singlet b state in one single
Fock space run. The FSCCSD input reads then

&RELCCSD NELEC=2,2,3,3, DOFSPC=T, DOENER=F &END

&CCFSPC MAXIT=100, NACTP=2,2,0,0, DOEA2=T, &END

 This example can be found as test #52 in the Dirac test set.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

NTOL
----

Convergence tolerance of the DIIS error vector (10\ :sup:`--NTOL`\ ).

*Default:*

::

    NTOL=12

MAXIT
-----

Maximum number of iterations to solve the CCSD equations.

*Default:*

::

    MAXIT=30

MAXITmn
~~~~~~~

Maximum number of iterations to solve the CCSD equations on each of the
Fock-space sectors.

*Default:*

::

    MAXIT00 = MAXIT01 = MAXIT10 = MAXIT11 = MAXIT20 = MAXIT02 = MAXIT

MAXDIM
------

Maximum dimension of the reduced space.

*Default:*

::

    MAXDIM=8

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

FSSECT
------

Active Fock space sectors (m,n) specifying the number of valence holes
(m) and valence particles (n) with respect to the reference determinant
(ordered 00,01,10,11,02,20).

*Default:*

::

    FSSECT=0,0,0,0,0,0

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:CCFSPC <http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:CCFSPC>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:RELCCSD:CCFSPC>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:RELCCSD:CCFSPC&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:RELCCSD:CCFSPC&action=edit>`_
-  `History </dirac/index.php?title=Manual:RELCCSD:CCFSPC&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:RELCCSD:CCFSPC>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:RELCCSD:CCFSPC>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:RELCCSD:CCFSPC>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:RELCCSD:CCFSPC&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:RELCCSD:CCFSPC&oldid=30201>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 10:11, 7 August 2011.
-  This page has been accessed 1,467 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image1| image:: /dirac/images/thumb/6/65/HJAaJ.png/41px-HJAaJ.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
