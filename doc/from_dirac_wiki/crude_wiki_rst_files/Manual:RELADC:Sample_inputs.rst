Manual:RELADC:Sample inputs
===========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Example 1:

A bent triatomic molecule in *C*\ :sub:`*s*`\  symmetry. We have only
one fermionic IRREP and complex MOs. This is indicated in the ADC output
reporting complex ADC matrices. We want to do single ionizations
(DOSIPS=T) at the ADC(3) level (ADCLEVEL=3) including constant diagrams
(DOCONST=default) and compute all final states in symmetry 1, the only
available one (XREPS=1). If the basis set is large, we restrict the
iterations for the constant diagrams by releasing convergence a bit
(VCONV=1.0E-05) beacause this step takes longest. We request 1000
Lanczos iterations and consider ionization energies up to 100 eV on
screen. The complete spectrum is written to SSPEC.1. You get as many
entries as you have Lanczos iterations. A characteristic of the
iterative diagonalizer is that the eigenvalues at the edge of the
spectrum converge very fast and are reproduced for high iteration
numbers. In order to get higher eigenvalues equally tightly converged
the number of iterations has to be increased accordingly. In this case
the spurious (already converged) eigenvalues are projected out. For the
*C*\ :sub:`*s*`\  molecule the ADC input would then look like this:

::

    $RELADC DOSIPS=T, ADCLEVEL=3, XREPS=1, VCONV=1.0D-05 $END
    $LANINP SIPITER=1000, SIPPRNT=100.0 $END

Since we do not activate DIPs nothing else needs to be specified in the
&LANINP Namelist.

Example 2:

The cadmium atom in *D*\ :sub:`2*h*`\  with .SPINFREE Hamiltonian. The
*D*\ :sub:`2*h*`\  double group has two fermionic IRREPS (g and u)
leading to XREPS=1,2. By using the .SPINFREE option a real calculation
in the ordinary *D*\ :sub:`2*h*`\  point group is performed in DIRAC
leading to eight symmetries as known from the character table. Therefore
we set XREPS=1,2,3,4,5,6,7,8 if we need all symmetries. The relation of
the irrep name to the number can be obtained from the Dirac output and
is reproduced in the ADC code as well. After the calculation the
*A*\ :sub:`1*g*`\  final states are the in SSPEC.01, the
*B*\ :sub:`2*g*`\  in SSPEC.02 asf. The complete spectrum is then the
merge of SSPEC.01 ... SSPEC.08. Remember that for the plot we only need
the IP and the corresponding pole strength. We therefore do a grep '@'
on the merged SSPEC.X files (cat SSPEC.\* > SPEC.all) and obtain all
lines in each symmetry. If we need only a specific range we do "sort -n
SPEC.all > SPEC.range" and edit according to our needs. Then we can use
gnuplot with "plot "SPEC.range" u 1:2 w i". For ADC(3) with constant
diagrams and 600 Lanczos iterations we would have an input like this:

::

    **DIRAC
    .TITLE
     input for Cd atom
    .WAVE F
    .4INDEX
    **GENERAL
    .DIRECT
     1 1 1
    **INTEGRALS
    *READIN
    *TWOINT
    .SOFOCK
    .SCREEN
    1.E-16
    **HAMILTONIAN
    .SPINFREE
    **WAVE FUNCTIONS
    .DHF
    .RELADC
    *DHFCAL
    .CLOSED
     30 18
    .FCKCNV
    5.0E-09
    .INTFLG
    1 1 1
    **MOLTRA
    .INTFLG
    1 1 1
    .CORE
    1..8
    1..6
    .ACTIVE
    9..25
    7..22
    *END OF
    $RELADC DOSIPS=T, ADCLEVEL=3, XREPS=1,2,3,4,5,6,7,8, $END
    $LANINP SIPITER=600 $END

 Back to `RELADC </dirac/index.php/Manual:RELADC>`_ manual.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:RELADC:Sample\_inputs <http://wiki.chem.vu.nl/dirac/index.php/Manual:RELADC:Sample_inputs>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:RELADC:Sample_inputs>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:RELADC:Sample_inputs&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:RELADC:Sample_inputs&action=edit>`_
-  `History </dirac/index.php?title=Manual:RELADC:Sample_inputs&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:RELADC:Sample_inputs>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:RELADC:Sample_inputs>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:RELADC:Sample_inputs>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:RELADC:Sample_inputs&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:RELADC:Sample_inputs&oldid=34221>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 17:01, 2 March 2012.
-  This page has been accessed 1,078 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
