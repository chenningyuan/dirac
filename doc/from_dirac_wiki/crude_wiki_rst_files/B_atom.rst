B atom
======

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Open-shell SCF and small CI relativistic calculations on the Boron atom
-----------------------------------------------------------------------

We seek to get the first excited state, 2P\ :sub:`3/2`\ , of B (X
2P\ :sub:`1/2`\ ), and find out the 2P\ :sub:`1/2`\ -2P\ :sub:`3/2`\ 
fine structure splitting.

Contents
--------

-  `1 Open-shell SCF and small CI relativistic calculations on the Boron
   atom <#Open-shell_SCF_and_small_CI_relativistic_calculations_on_the_Boron_atom>`_

   -  `1.1 Complete open-shell CI <#Complete_open-shell_CI>`_
   -  `1.2 MO reordering <#MO_reordering>`_

-  `2 Correlated calculations <#Correlated_calculations>`_
-  `3 Conclusion <#Conclusion>`_

 We select the Dirac-Coulomb Hamiltonian and neglect the (SS\|SS) class
of integrals.

Complete open-shell CI
^^^^^^^^^^^^^^^^^^^^^^

Possible way to get the value is performing small CI over the active
open-shell space, containing six 2p orbitals, which are split due to the
spin-orbital interaction. Note, that for pedagogical purposes we have
employed several computational symmetries.

Starting from the average SCF occupation (one electron distributed over
2p\ :sub:`1/2`\  and 2p\ :sub:`3/2`\  shells),

::

    pam  --mol=B.sto-3g.lsym.mol -inp=B.2Pav.dc_rkb.scf_2fs_cosci.inp  

::

    pam  --mol=B.sto-3g.D2h.mol --inp=B.2Pav.dc_rkb.scf_2fs_cosci.inp

for one-fermion symmetries,

::

    pam  --mol=B.sto-3g.C2v.mol --inp=B.2Pav.dc_rkb.scf_1fs_cosci.inp

we get the value of 17.179488 cm-1 and in the outputs we see the proper
degeneracy of resulting total CI energies - 2 and 4.

By averaging the COSCI energies of the ground and the states we get the
total SCF energy of the average occupation;

::

     ((2*(-24.249116730099) + 4*(-24.249038454591)))/6 = -24.249064546427
                                         SCF energy:     -24.249064546425263

MO reordering
^^^^^^^^^^^^^

To apply the MO reordering scheme, first, one has to obtain the MO
coefficients file of the ground state:

::

    pam  --outcmo --mol=B.sto-3g.D2h.mol --inp=B.2P12.dc_rkb.scf_2fs.inp
    pam  --outcmo --mol=B.sto-3g.C1.mol  --inp=B.2P12.dc_rkb.scf_1fs.inp

and use them for calculation of the excited state:

::

    pam --incmo --mol=B.sto-3g.D2h.mol  --inp=B.2P32.dc_rkb.scf_2fs_reord_ovlsel.inp
    pam --incmo --mol=B.sto-3g.C1.mol     --inp=B.2P32.dc_rkb.scf_1fs_reord_ovlsel.inp

For the linear symmetry one can employ the M\_J specified occupation of
shells both for the ground and the first excited state. Only the preDHF
reordering works in this case, not the overlap selection.

::

     pam --outcmo --mol=B.sto-3g.lsym.mol --inp=B.2P12.dc_rkb.scf_2fs_mj.inp
     pam --incmo --mol=B.sto-3g.lsym.mol --inp=B.2P32.dc_rkb.scf_2fs_reord_mj.inp

Finally, in the output one can see proper order of spinors for the
excited state - first 2p\ :sub:`3/2`\  spinors (four-fold degeneracy),
and then the 2p\ :sub:`1/2`\  spinors (two-fold degeneracy):

::

    * Fermion symmetry E1u
      * Open shell #1, f = 0.2500
     -0.28691011961454  ( 4)
      * Virtual eigenvalues, f = 0.0000
      0.12724113329935  ( 2)

Regarding total energies we get

::

    SCF energy 2_P1/2: -24.249116730391307
    SCF energy 2_P3/2: -24.249038454664991
    FSS : (-24.249116730391307+24.249038454664991)*219474.631280634 cm-1 = 17.179536 cm-1

what is very close to the above given COSCI value of 17.179488 cm-1.

Correlated calculations
-----------------------

... to add ...

Conclusion
----------

Further quality improvement of the theoretical data would be in
enlarging basis set. The current one, decontracted STO-3G, is too small
and is chosen for demonstrative purposes only.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/B\_atom <http://wiki.chem.vu.nl/dirac/index.php/B_atom>`_"

Views
'''''

-  `Article </dirac/index.php/B_atom>`_
-  `Discussion </dirac/index.php?title=Talk:B_atom&action=edit>`_
-  `View source </dirac/index.php?title=B_atom&action=edit>`_
-  `History </dirac/index.php?title=B_atom&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=B_atom>`_

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

-  `What links here </dirac/index.php/Special:Whatlinkshere/B_atom>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/B_atom>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable version </dirac/index.php?title=B_atom&printable=yes>`_
-  `Permanent link </dirac/index.php?title=B_atom&oldid=32641>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 09:59, 11 November 2011.
-  This page has been accessed 532 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
