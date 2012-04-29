Manual:One-electron operators
=============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 One-electron operators <#One-electron_operators>`_

   -  `1.1 Syntax for the specification of one-electron
      operators <#Syntax_for_the_specification_of_one-electron_operators>`_

      -  `1.1.1 Operator types <#Operator_types>`_
      -  `1.1.2 Operator specification <#Operator_specification>`_

   -  `1.2 List of one-electron
      operators <#List_of_one-electron_operators>`_
   -  `1.3 Examples of using various
      operators <#Examples_of_using_various_operators>`_

      -  `1.3.1 XAVECTOR <#XAVECTOR>`_
      -  `1.3.2 Kinetic part of the Dirac
         Hamiltonian <#Kinetic_part_of_the_Dirac_Hamiltonian>`_
      -  `1.3.3 Dipole moment as finite field
         parturbation <#Dipole_moment_as_finite_field_parturbation>`_
      -  `1.3.4 Fermi-contact integrals <#Fermi-contact_integrals>`_
      -  `1.3.5 Cartesian moment expectation
         value <#Cartesian_moment_expectation_value>`_

One-electron operators
======================

Syntax for the specification of one-electron operators
------------------------------------------------------

A general property operator in 4-component calculations is generated
from linear combinations of the basic form

.. figure:: /dirac/images/math/4/4/d/44d4fd16bd6198c74c77c33ee41cd9a3.png
   :align: center
   :alt:  \\hat{P} = f M\_{4 \\times 4} \\Omega

    \\hat{P} = f M\_{4 \\times 4} \\Omega
with the scalar factor *f* and the scalar operator Ω, and where |M\_{4
\\times 4}| is one of the following 4 |\\times| 4 matrices:

.. figure:: /dirac/images/math/c/0/c/c0ce829776217dff5580bfeee7b6ab94.png
   :align: center
   :alt:  I\_{4 \\times 4}, \\gamma\_5,

    I\_{4 \\times 4}, \\gamma\_5,
α\ :sub:`*x*`\ ,α\ :sub:`*y*`\ ,α\ :sub:`*z*`\ ,

Σ\ :sub:`*x*`\ ,Σ\ :sub:`*y*`\ ,Σ\ :sub:`*z*`\ ,

βΣ\ :sub:`*x*`\ ,βΣ\ :sub:`*y*`\ ,βΣ\ :sub:`*z*`\ ,

*i*βα\ :sub:`*x*`\ ,*i*βα\ :sub:`*y*`\ ,*i*βα\ :sub:`*y*`\ 

Operator types
~~~~~~~~~~~~~~

There is 19 basic operator types used in DIRAC.

*Operator types*

**Keyword**

**Operator form**

**Nr. Factors**

DIAGONAL

.. figure:: /dirac/images/math/c/b/6/cb65c6e8a7fe29140e2b4e7c9c3e7fab.png
   :align: center
   :alt:  f I\_{4 \\times 4} \\Omega

    f I\_{4 \\times 4} \\Omega
1

GAMMA5

*f*γ\ :sub:`5`\ Ω

1

XALPHA

*f*α\ :sub:`*x*`\ Ω

1

YALPHA

*f*α\ :sub:`*y*`\ Ω

1

ZALPHA

*f*α\ :sub:`*z*`\ Ω

1

XSIGMA

*f*Σ\ :sub:`*x*`\ Ω

1

YSIGMA

*f*Σ\ :sub:`*y*`\ Ω

1

ZSIGMA

*f*Σ\ :sub:`*z*`\ Ω

1

XAVECTOR

*f*\ :sub:`1`\ α\ :sub:`*y*`\ Ω\ :sub:`*z*`\  −
*f*\ :sub:`2`\ α\ :sub:`*z*`\ Ω\ :sub:`*y*`\ 

2

YAVECTOR

*f*\ :sub:`1`\ α\ :sub:`*z*`\ Ω\ :sub:`*x*`\  −
*f*\ :sub:`2`\ α\ :sub:`*x*`\ Ω\ :sub:`*z*`\ 

2

ZAVECTOR

*f*\ :sub:`1`\ α\ :sub:`*x*`\ Ω\ :sub:`*y*`\  −
*f*\ :sub:`2`\ α\ :sub:`*y*`\ Ω\ :sub:`*x*`\ 

2

ALPHADOT

*f*\ :sub:`1`\ α\ :sub:`*x*`\ Ω\ :sub:`*x*`\  +
*f*\ :sub:`2`\ α\ :sub:`*y*`\ Ω\ :sub:`*y*`\  +
*f*\ :sub:`3`\ α\ :sub:`*z*`\ Ω\ :sub:`*z*`\ 

3

XBETASIG

*f*βΣ\ :sub:`*x*`\ Ω

1

YBETASIG

*f*βΣ\ :sub:`*y*`\ Ω

1

ZBETASIG

*f*βΣ\ :sub:`*z*`\ Ω

1

XiBETAAL

*f**i*βα\ :sub:`*x*`\ Ω

1

YiBETAAL

*f**i*βα\ :sub:`*y*`\ Ω

1

ZiBETAAL

*f**i*βα\ :sub:`*z*`\ Ω

1

BETAGAMMA5

*f**i*βγ\ :sub:`5`\ Ω

1

Operator specification
~~~~~~~~~~~~~~~~~~~~~~

Operators are specified by the keyword .OPERATOR with the following
arguments:

::

    .OPERATOR
     'operator name'
     operator type
     labels for each component
     FACTORS
     factors for each component
     COMFACTOR
     common factor for all components

Note that the arguments following the keyword .OPERATOR must start with
a blank. The arguments are optional, except for the operator label.

List of one-electron operators
------------------------------

+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| **Keyword**   |                                                                                                                                     |                      | **Components**                                                          |                                                                 |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| MOLFIELD      | Nuclear attraction integrals                                                                                                        | *(Symmetric)*        | MOLFIELD                                                                | |image16|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| OVERLAP       | Overlap integrals                                                                                                                   | *(Symmetric)*        | OVERLAP                                                                 | Ω\ :sub:`1`\  = 1                                               |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| BETAMAT       | Overlap integrals, only SS-block                                                                                                    | *(Symmetric)*        | BETAMAT                                                                 | Ω\ :sub:`1`\  = 1                                               |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| DIPLEN        | Dipole length integrals                                                                                                             | *(Symmetric)*        | XDIPLEN                                                                 | Ω\ :sub:`1`\  = *x*                                             |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | YDIPLEN                                                                 | Ω\ :sub:`2`\  = *y*                                             |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | ZDIPLEN                                                                 | Ω\ :sub:`3`\  = *z*                                             |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| DIPVEL        | Dipole velocity integrals                                                                                                           | *(Anti-symmetric)*   | XDIPVEL                                                                 | |image17|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | YDIPVEL                                                                 | |image18|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | ZDIPVEL                                                                 | |image19|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| QUADRUP       | Quadrupole moments integrals                                                                                                        | *(Symmetric)*        | XXQUADRU                                                                | |image20|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | XYQUADRU                                                                | |image21|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | XZQUADRU                                                                | |image22|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | YYQUADRU                                                                | |image23|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | YZQUADRU                                                                | |image24|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | ZZQUADRU                                                                | |image25|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| SPNORB        | Spatial spin-orbit integrals                                                                                                        | *(Anti-symmetric)*   | X1SPNORB                                                                | |image26|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | Y1SPNORB                                                                | |image27|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | Z1SPNORB                                                                | |image28|                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| SECMOM        | Second moments integrals                                                                                                            | *(Symmetric)*        | XXSECMOM                                                                | Ω\ :sub:`1`\  = *x**x*                                          |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | XYSECMOM                                                                | Ω\ :sub:`1`\  = *x**y*                                          |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | XZSECMOM                                                                | Ω\ :sub:`1`\  = *x**z*                                          |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | YYSECMOM                                                                | Ω\ :sub:`1`\  = *y**y*                                          |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | YZSECMOM                                                                | Ω\ :sub:`1`\  = *y**z*                                          |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
|               |                                                                                                                                     |                      | ZZSECMOM                                                                | Ω\ :sub:`1`\  = *z**z*                                          |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| THETA         | Traceless theta quadrupole integrals                                                                                                |                      |                                                                         |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| CARMOM        | Cartesian moments integrals,Symmetric integrals,(l + 1)(l + 2)/2 components ( l = i + j + k) ; see also the example further below   |                      | Ω\ :sub:`*l*`\  = *x*\ :sup:`*i*`\ *y*\ :sup:`*j*`\ *z*\ :sup:`*k*`\    |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| SPHMOM        | Spherical moments integrals (real combinations), symmetric integrals, (2l+1) components, ( m = + 0, - 1, + 1,..., + l)              |                      | |image29|                                                               |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| SOLVENT       | Electronic solvent integrals                                                                                                        |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| FERMI C       | One-electron Fermi contact integrals                                                                                                |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| PSO           | Paramagnetic spin-orbit integrals                                                                                                   |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| SPIN-DI       | Spin-dipole integrals                                                                                                               |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| DSO           | Diamagnetic spin-orbit integrals                                                                                                    |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| SDFC          | Spin-dipole + Fermi contact integrals                                                                                               |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| HDO           | Half-derivative overlap integrals                                                                                                   |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| S1MAG         | Second order contribution from overlap matrix to magnetic properties                                                                |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| ANGLON        | Angular momentum around the nuclei                                                                                                  |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| ANGMOM        | Electronic angular momentum around the origin                                                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| LONMOM        | London orbital contribution to angular momentum                                                                                     |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| MAGMOM        | One-electron contributions to magnetic moment                                                                                       |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| KINENER       | Electronic kinetic energy                                                                                                           |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+
| DSUSNOL       | Diamagnetic susceptibility without London contribution                                                                              |
+---------------+-------------------------------------------------------------------------------------------------------------------------------------+----------------------+-------------------------------------------------------------------------+-----------------------------------------------------------------+

add more from
`http://dirac.chem.sdu.dk/dirlab/Dirlabdoc/Master.html <http://dirac.chem.sdu.dk/dirlab/Dirlabdoc/Master.html>`_
...

Examples of using various operators
-----------------------------------

We give here several concrete examples on how to construct operators for
various properties.

XAVECTOR
~~~~~~~~

An example:

::

    .OPERATOR
     'B_x'
     XAVECTOR
     ZDIPLEN
     YDIPLEN
     COMFACTOR
     -68.517999904721

Kinetic part of the Dirac Hamiltonian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As another example the kinetic part of the Dirac Hamiltonian may be
specified by:

::

    .OPERATOR
     'Kin energy'
     ALPHADOT
     XDIPVEL
     YDIPVEL
     ZDIPVEL
     COMFACTOR
     -68.51799475

where -68.51799475 is -c/2

The program will assume all operators to be Hermitian and will therefore
insert an imaginary phase *i* if necessary (applies to antisymmetric
scalar operators).

If no other arguments are given, the program assumes the operator to be
a diagonal operator and expects the operator name to be the component
label, for instance:

::

    .OPERATOR
     OVERLAP

Dipole moment as finite field parturbation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another example is the finite perturbation calculation with the
|\\hat{z}| dipole length operator added to the Hamiltonian (don't forget
to decrease the symmetry of your system):

.. figure:: /dirac/images/math/d/3/7/d3798e00fb55a6afa6f89f1bf96d80fc.png
   :align: center
   :alt: \\hat{H} = \\hat{H}\_0 + 0.01 \\cdot \\hat{z}

   \\hat{H} = \\hat{H}\_0 + 0.01 \\cdot \\hat{z}
::

    .OPERATOR
     ZDIPLEN
     COMFACTOR
     0.01

Fermi-contact integrals
~~~~~~~~~~~~~~~~~~~~~~~

Here is an example where the Fermi-contact (FC) integrals for a certain
nucleus are added to the Hamiltonian in a finite-field calculation.
Let's assume you are looking at a PbX dimer (order in the .mol file: 1.
Pb, 2. X) and you want to add to the Dirac-Coulomb
`Hamiltonian </dirac/index.php/Manual:HAMILTONIAN>`_ the FC integrals
for the Pb nucleus as a perturbation with a given field-strength
(FACTORS).

**Important note:** The raw density values obtained after the fit of
your finite-field energies need to be scaled by
|\\frac{1}{\\frac{4\*\\pi\*g\_{e}}{3}} = \\frac{1}{8.3872954891254192}|,
a factor that originates from the definition of the operator for
calculating the density at the nucleus.

::

    **HAMILTONIAN
    .OPERATOR
    'Density at nucleus'
    DIAGONAL
    'FC Pb 01'
    FACTORS
    -0.000000001

... and here is an example of how-to calculate the electron density at
the nucleus as an expectation value |\\langle 0 \\vert \\delta(r-R)
\\vert 0 \\rangle| for a Dirac-Coulomb HF wave function including a
decomposition of the molecular orbital contributions to the density:

::

    **DIRAC
    .WAVE FUNCTION
    .PROPERTIES
    **HAMILTONIAN
    **WAVE FUNCTION
    .DHF
    **PROPERTIES
    .RHONUC
    *EXPECTATION
    .ORBANA
    *END OF

Cartesian moment expectation value
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the following example I am calculating a cartesian moment expectation
value |\\langle 0 \\vert x^1 y^2 z^3 \\vert 0 \\rangle| for a
Levy-Leblond HF wave function:

::

    **DIRAC
    .WAVE FUNCTION
    .PROPERTIES
    **HAMILTONIAN
    .LEVY-LEBLOND
    **WAVE FUNCTION
    .DHF
    **PROPERTIES
    *EXPECTATION
    .OPERATOR
     CM010203
    *END OF

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:One-electron\_operators <http://wiki.chem.vu.nl/dirac/index.php/Manual:One-electron_operators>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:One-electron_operators>`_
-  `Discussion </dirac/index.php/Talk:Manual:One-electron_operators>`_
-  `View
   source </dirac/index.php?title=Manual:One-electron_operators&action=edit>`_
-  `History </dirac/index.php?title=Manual:One-electron_operators&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:One-electron_operators>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:One-electron_operators>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:One-electron_operators>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:One-electron_operators&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:One-electron_operators&oldid=33751>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 13:18, 16 February 2012.
-  This page has been accessed 3,925 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |M\_{4 \\times
4}| image:: /dirac/images/math/3/c/2/3c285f6a74ed8bcbe8fa513c1051dbc9.png
.. |\\times| image:: /dirac/images/math/9/e/e/9eedd61e32f7a8e70e171028a7e5dc08.png
.. | \\Omega\_1 = \\sum\_K \\hat{V}\_{iK}
| image:: /dirac/images/math/c/d/9/cd9010c50a5657b81faf33172ba938ee.png
.. | \\Omega\_1 = \\frac{\\partial}{\\partial x}
| image:: /dirac/images/math/7/4/b/74bf4a24f225395908b1f595004f0698.png
.. | \\Omega\_2 = \\frac{\\partial}{\\partial y}
| image:: /dirac/images/math/c/9/4/c949136d82a78936b5a5d0d42bea30ea.png
.. | \\Omega\_3 = \\frac{\\partial}{\\partial z}
| image:: /dirac/images/math/c/f/a/cfad7b2b0b1d3a7602e5160f08312b1c.png
.. | \\Omega\_1 = \\frac{\\partial^2}{\\partial x \\partial x}
| image:: /dirac/images/math/3/1/1/311110523bfe59709fe22cbc81d3df34.png
.. | \\Omega\_2 = \\frac{\\partial^2}{\\partial x \\partial y}
| image:: /dirac/images/math/7/1/e/71e2cb28a73d6296b583631124ed85e2.png
.. | \\Omega\_3 = \\frac{\\partial^2}{\\partial x \\partial z}
| image:: /dirac/images/math/4/c/e/4cec989091b60eac53988d93427df274.png
.. | \\Omega\_4 = \\frac{\\partial^2}{\\partial y \\partial y}
| image:: /dirac/images/math/a/7/5/a755d873c65b4f45d1210dea1ac0b956.png
.. | \\Omega\_5 = \\frac{\\partial^2}{\\partial y \\partial z}
| image:: /dirac/images/math/0/b/1/0b17260e58563d0a757fc2a3ee29e235.png
.. | \\Omega\_6 = \\frac{\\partial^2}{\\partial z \\partial z}
| image:: /dirac/images/math/3/d/2/3d2af67afc6283ea01a735299faca68f.png
.. | \\Omega\_1 = \\sum\_K \\frac{l\_x}{r\_{iK}^3}
| image:: /dirac/images/math/9/e/5/9e5d7e15ed6ce440250cacf050863fc6.png
.. | \\Omega\_2 = \\sum\_K \\frac{l\_y}{r\_{iK}^3}
| image:: /dirac/images/math/9/e/f/9ef128ed676ad54d70fe3870d5705358.png
.. | \\Omega\_3 = \\sum\_K \\frac{l\_z}{r\_{iK}^3}
| image:: /dirac/images/math/2/8/3/28326ab19cde00117469d456daeaf12f.png
.. | \\Omega\_l = R\_{l \\pm \_m}
| image:: /dirac/images/math/f/f/8/ff821df085ec7a5da60b396be2d35248.png
.. |image16| image:: /dirac/images/math/c/d/9/cd9010c50a5657b81faf33172ba938ee.png
.. |image17| image:: /dirac/images/math/7/4/b/74bf4a24f225395908b1f595004f0698.png
.. |image18| image:: /dirac/images/math/c/9/4/c949136d82a78936b5a5d0d42bea30ea.png
.. |image19| image:: /dirac/images/math/c/f/a/cfad7b2b0b1d3a7602e5160f08312b1c.png
.. |image20| image:: /dirac/images/math/3/1/1/311110523bfe59709fe22cbc81d3df34.png
.. |image21| image:: /dirac/images/math/7/1/e/71e2cb28a73d6296b583631124ed85e2.png
.. |image22| image:: /dirac/images/math/4/c/e/4cec989091b60eac53988d93427df274.png
.. |image23| image:: /dirac/images/math/a/7/5/a755d873c65b4f45d1210dea1ac0b956.png
.. |image24| image:: /dirac/images/math/0/b/1/0b17260e58563d0a757fc2a3ee29e235.png
.. |image25| image:: /dirac/images/math/3/d/2/3d2af67afc6283ea01a735299faca68f.png
.. |image26| image:: /dirac/images/math/9/e/5/9e5d7e15ed6ce440250cacf050863fc6.png
.. |image27| image:: /dirac/images/math/9/e/f/9ef128ed676ad54d70fe3870d5705358.png
.. |image28| image:: /dirac/images/math/2/8/3/28326ab19cde00117469d456daeaf12f.png
.. |image29| image:: /dirac/images/math/f/f/8/ff821df085ec7a5da60b396be2d35248.png
.. |\\hat{z}| image:: /dirac/images/math/8/7/5/875326710761d5ed42bdc6e30c4bf962.png
.. |\\frac{1}{\\frac{4\*\\pi\*g\_{e}}{3}} =
\\frac{1}{8.3872954891254192}| image:: /dirac/images/math/4/5/a/45a0f4745029158594b0d6785db4a0ba.png
.. |\\langle 0 \\vert \\delta(r-R) \\vert 0
\\rangle| image:: /dirac/images/math/c/1/9/c19cd81b1c90dd7c6e4f5c78e4368891.png
.. |\\langle 0 \\vert x^1 y^2 z^3 \\vert 0
\\rangle| image:: /dirac/images/math/7/9/b/79bc308bd48e8cd8e79b832d34b640f5.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
