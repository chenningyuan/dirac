Tutorial:Getting started DIRAC11
================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 References <#References>`_

 DIRAC uses two files as inputs: the input file **\*.inp**, which
determines what should be calculated, and either a geometry file
**\*.xyz**, that provides the nuclear coordinates, or a molecule file
**\*.mol**, which defines the molecular geometry and also the basis set.
These files are described in detail in the
`Manual </dirac/index.php/Manual>`_. When DIRAC runs it writes its
output to a text file, named after the molecule and input file names. If
this file already exists, for example from a previous calculation, the
old file is renamed to make place for the new output. The **pam** python
script that launches the main DIRAC executable also writes some general
information to standard output.

Let's try it out and run our first DIRAC calculation with the following
minimal input (**hf.inp**) to calculate the Dirac-Coulomb Hartree-Fock
ground state (with the default simple Coulombic correction
(`LVCORR </dirac/index.php/Manual:HAMILTONIAN#.LVCORR>`_) to approximate
the contribution from SS\|SS integrals
\ :sup:``[1] <#_note-visscher-97>`_`\ . ):

::

    **DIRAC
    .WAVE FUNCTION
    **WAVE FUNCTION
    .SCF
    **MOLECULE
    *BASIS
    .DEFAULT
    cc-pVDZ
    *END OF INPUT

together with the following example geometry file (**methanol.xyz**):

::

    6
    Methanol. This is a comment line.
    C       0.000000000000       0.138569980000       0.355570700000   
    O       0.000000000000       0.187935770000      -1.074466460000  
    H       0.882876920000      -0.383123830000       0.697839450000  
    H      -0.882876940000      -0.383123830000       0.697839450000  
    H       0.000000000000       1.145042790000       0.750208830000  
    H       0.000000000000      -0.705300580000      -1.426986340000

Now start the calculation by executing:

::

    pam  --mol=methanol.xyz --inp=hf.inp

If everything works fine the results of the calculation will be written
to **hf\_methanol.out**. In addition an archive file *hf\_methanol.tgz*
will be created. Depending on the type of calculation the archive file
may include:

-  ``DFCOEF`` containing molecular orbital coefficients and possibly
   energies.
-  ``*.cube`` containing plots of molecular properties, typically
   electron densities, in the gridded Gaussian Cube format.
-  Additional files for restarting particular types of calculations.

References
----------

#. `↑ <#_ref-visscher-97_0>`_ L. Visscher, *Approximate molecular
   Dirac-Coulomb calculations using a simple Coulombic correction*,
   Theor. Chem. Acc. **98**, *68* (1997) `electronic
   version <http://dx.doi.org/10.1007/s002140050280>`_.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Tutorial:Getting\_started\_DIRAC11 <http://wiki.chem.vu.nl/dirac/index.php/Tutorial:Getting_started_DIRAC11>`_"

Views
'''''

-  `Article </dirac/index.php/Tutorial:Getting_started_DIRAC11>`_
-  `Discussion </dirac/index.php?title=Talk:Tutorial:Getting_started_DIRAC11&action=edit>`_
-  `View
   source </dirac/index.php?title=Tutorial:Getting_started_DIRAC11&action=edit>`_
-  `History </dirac/index.php?title=Tutorial:Getting_started_DIRAC11&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Tutorial:Getting_started_DIRAC11>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Tutorial:Getting_started_DIRAC11>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Tutorial:Getting_started_DIRAC11>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Tutorial:Getting_started_DIRAC11&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Tutorial:Getting_started_DIRAC11&oldid=32561>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 09:40, 11 November 2011.
-  This page has been accessed 730 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
