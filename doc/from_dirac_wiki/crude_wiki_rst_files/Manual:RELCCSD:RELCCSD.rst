Manual:RELCCSD:RELCCSD
======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &RELCCSD -- Generic and control
   input <#.26RELCCSD_--_Generic_and_control_input>`_

   -  `1.1 DOENER <#DOENER>`_
   -  `1.2 DOFSPC <#DOFSPC>`_
   -  `1.3 DOFOPR <#DOFOPR>`_
   -  `1.4 NELEC <#NELEC>`_
   -  `1.5 NELEC\_F1 <#NELEC_F1>`_
   -  `1.6 NELEC\_F2 <#NELEC_F2>`_
   -  `1.7 IPRNT <#IPRNT>`_
   -  `1.8 TIMING <#TIMING>`_
   -  `1.9 DEBUG <#DEBUG>`_
   -  `1.10 DOSOPR <#DOSOPR>`_
   -  `1.11 DOSORT <#DOSORT>`_
   -  `1.12 NFROZ <#NFROZ>`_
   -  `1.13 INTERFACE <#INTERFACE>`_

&RELCCSD -- Generic and control input
=====================================

DOENER
------

Calculate the CCSD(T) energy based on the standard single determinant
ansatz.

*Default:*

::

    DOENER=T

DOFSPC
------

Calculate CCSD energies using Fock space CC ansatz. In this method one
chooses a closed shell configuration from which the open shell
configurations of interest can be obtained by adding or subtracting one
or two electrons. The theory allows full coupling between the states
that belong to the same "sector" of Fock space (the set of
configurations that differ by the same amount of electrons from the
closed shell reference) and thus enables a multireference treatment of
open shell states. Energies are calculated and reported as electron
affinities (sectors [0,1] and [0,2]), ionization energies(sectors [1,0]
and [2,0]), or excitation energies (sector [1,1]) relative to the
reference closed shell state.

*Default:*

::

    DOFSPC=F

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

DOFOPR
------

Calculate the (relaxed) density matrix and first-order properties (as
specified in the section
`Manual:PROPERTIES </dirac/index.php/Manual:PROPERTIES>`_). These are
only available for the MP2 wave function. The MP2 density matrix is
**not** available for spinfree (or nonrelatvistic) calculations and the
program will stop when this combination of options is chosen.

*Default:*

::

    DOFOPR=F

NELEC
-----

Number of active electrons. This variable determines the reference
determinant to be used in the exponential expansion of the wave
function. Since the default values correspond to the information passed
on by the MOLTRA code on basis of the Hartree-Fock occupations and
chosen range of active orbitals in
`MOLTRA </dirac/index.php/Manual:MOLTRA#.ACTIVEMOLTRA>`_, there is
usually no need to specify this variable manually. Notable exception is
the treatment of a high-spin open shell state with the single reference
CCSD(T) ansatz. In such a case one first needs to determine the number
of closed shell electrons in each irrep and then add to these
occupations the open shell electrons (see below for an example).

*Arguments:* Integer (NELEC(I),I=1,NFSYM\*2).

*Default:* Active electrons in these irreps (written by
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_).

 A simple example is the oxygen molecule in which the 1s electrons are
kept frozen. The active space then consists of 12 electrons, of which
there are 10 closed shell electrons. The closed shell electrons go in
the nonbonding 2s σ\ :sub:`g`\  and σ\ :sub:`u`\  orbitals and the
bonding 2p σ\ :sub:`g`\  and π\ :sub:`u`\  orbitals, while the 2 open
shell electrons are distributed over the two antibonding 2p
π\ :sub:`g`\  orbitals. Let us suppose that we want to take the
M\ :sub:`S`\ =1 state as our reference and assign alpha spin to both our
open shell electrons. In the double group (D\ :sub:`2h`\ \ :sup:`\*`\ )
we do not distinguish between σ and π, so we need to add all alpha
electrons in gerade orbitals (2s σ\ :sub:`g`\  + 2p σ\ :sub:`g`\  + 2p
π\ :sub:`g`\  = 4), beta electrons in gerade orbitals (2s σ\ :sub:`g`\ 
+ 2p σ\ :sub:`g`\  = 2), alpha electrons in ungerade orbitals (2s
σ\ :sub:`u`\  + 2p π\ :sub:`u`\ pi\_u = 3) and beta electrons in
ungerade orbitals (2s σ\ :sub:`u`\  + 2p π\ :sub:`u`\  = 3). This gives
the following occupation:

 NELEC=4,2,3,3

 Note that this determinant does not represent the exact ground state of
the oxygen molecule as this triplet is split by a few wave numbers due
to spin-spin and (second order) spin-orbit interactions. The lowest
state is the Omega=0 component that cannot be represented by a single
determinant. This state can be calculated using the Fock space method
(see under `&CCFSPC </dirac/index.php/Manual:RELCCSD:CCFSPC#DOEA2>`_).

NELEC\_F1
---------

Number of electrons in the gerade irreps of the Abelian symmetry group.

*Arguments:* Integer (NELEC\_F1(I),I=1,number of gerade irreps).

*Default:* Values set by `\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_
or `NELEC <#NELEC>`_.

NELEC\_F2
---------

Number of electrons in the ungerade irreps of the Abelian symmetry
group.

*Arguments:* Integer (NELEC\_F2(I),I=1,number of ungerade irreps).

*Default:* Values set by `\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_
or `NELEC <#NELEC>`_.

 Coming back to the example given above, the oxygen molecule, we now
show how this is done with the keywords given above. We again want to
take the M\_S=1 state as our reference. The irreps of Dinfh\* are
ordered as 1/2,-1/2,3/2,-3/2,.... so we need to consider the omega value
(giving the projection on the molecular axis of both spin and orbital
angular momentum) of the occupied oribitals. The sigma orbitals go in
the irreps 1/2 and -1/2 while the pi-orbitals span the four irreps
(1/2,-1/2,3/2,-3/2). Putting an alpha electron in a sigma orbital will
give an omega value of 1/2, while putting it in a pi-orbital can either
give -1/2 (when put in the orbitals with orbital momentum -1) or 3/2
(when put in the +1 orbital). Similarly the beta electrons go in irreps
-1/2 (for the sigma), -3/2 and 1/2 (for the pi). This makes the input
for our example :

 NELEC\_F1=2,3,1,0 NELEC\_F2=2,2,1,1

IPRNT
-----

Print level. Level 1 gives information on the convergence of the CCSD
iterations and may be useful.

*Default:*

::

    IPRNT=0

TIMING
------

Write timing information.

*Default:*

::

    TIMING=F

The detailed `NELEC\_F1 <#NELEC_F1>`_ and `NELEC\_F2 <#NELEC_F2>`_ input
is mainly of interest in atomic calculations where the program may use
the |D\_{\\infty h}| symmetry group. In that case this option permits
precise specification of the reference determinant.
`Here </dirac/index.php/Manual:RELCCSD:Sample_inputs>`_ we give an
example for the high spin state of the nitrogen atom.

`|image2| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

DEBUG
-----

Write debug information.

*Default:*

::

    DEBUG=F

DOSOPR
------

Calculate second-order properties (only for testing and based on the HF
wave function).

*Default:*

::

    DOSOPR=F

DOSORT
------

Sort the incoming integrals over molecular spinors. This is a necessary
step that can only be skipped in special (restart) cases.

*Default:*

::

    DOSORT=T

NFROZ
-----

Number of frozen electrons in each irrep. Makes T1 and T2 amplitudes
that concern frozen spinors zero. Only for test purposes.

*Arguments:* Integer (NFROZ(I),I=1,NFSYM\*2).

*Default:* 0

INTERFACE
---------

This keyword controls how RELCCSD should handle the (transformed)
integral files coming from `MOLTRA </dirac/index.php/Manual:MOLTRA>`_.
Used together with scheme 6 this keyword should be specified as

::

    INTERFACE='DIRAC6'

Since scheme 6 is the default for MOLTRA and RELCCSD, this is normally
not necessary. However, if one uses scheme 4 in MOLTRA, the keyword
**MUST** be specified as shown below.

::

    INTERFACE='DIRAC'

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:RELCCSD <http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:RELCCSD>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:RELCCSD:RELCCSD>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:RELCCSD:RELCCSD&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:RELCCSD:RELCCSD&action=edit>`_
-  `History </dirac/index.php?title=Manual:RELCCSD:RELCCSD&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:RELCCSD:RELCCSD>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:RELCCSD:RELCCSD>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:RELCCSD:RELCCSD>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:RELCCSD:RELCCSD&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:RELCCSD:RELCCSD&oldid=24121>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 19:48, 21 October 2010.
-  This page has been accessed 3,421 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |D\_{\\infty
h}| image:: /dirac/images/math/1/2/f/12fb65a005c771345f0bb591f3e8e9bd.png
.. |image2| image:: /dirac/images/thumb/6/65/HJAaJ.png/41px-HJAaJ.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
