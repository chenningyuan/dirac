Manual:HAMILTONIAN:DFT
======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Integration algorithm and
   thresholds <#Integration_algorithm_and_thresholds>`_

   -  `1.1 .SCREENING <#.SCREENING>`_
   -  `1.2 .POINTWISE <#.POINTWISE>`_
   -  `1.3 .TINYDENS <#.TINYDENS>`_

-  `2 Asymptotic correction <#Asymptotic_correction>`_

   -  `2.1 .GRAC <#.GRAC>`_
   -  `2.2 .SAOP! <#.SAOP.21>`_
   -  `2.3 .SAOP <#.SAOP>`_

-  `3 Spin magnetization TDDFT <#Spin_magnetization_TDDFT>`_

   -  `3.1 .NOSDFT <#.NOSDFT>`_
   -  `3.2 .COLLINEAR <#.COLLINEAR>`_
   -  `3.3 .BETASIGMA <#.BETASIGMA>`_

-  `4 Adiabatic local density approximation
   (ALDA) <#Adiabatic_local_density_approximation_.28ALDA.29>`_

   -  `4.1 .ALDA <#.ALDA>`_
   -  `4.2 .XALDA <#.XALDA>`_
   -  `4.3 .ALDA+ <#.ALDA.2B>`_
   -  `4.4 .ALDA- <#.ALDA->`_
   -  `4.5 .XALDA+ <#.XALDA.2B>`_
   -  `4.6 .XALDA- <#.XALDA->`_

-  `5 Other functionality <#Other_functionality>`_

   -  `5.1 .GAUNTSCALE <#.GAUNTSCALE>`_
   -  `5.2 References <#References>`_

 The specification of the functional is described in connection with the
`.DFT option in the
\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN#.DFT>`_ section

The DFT module is described in references
\ :sup:``[1] <#_note-Schipper:2000:jchemphys>`_`\ 
\ :sup:``[2] <#_note-Saue:2002:jcomputchem>`_`\ 
\ :sup:``[3] <#_note-Fossgaard:2003:jchemphys>`_`\ .

It uses standard non-relativistic functionals. This is justified in many
cases since various studies
\ :sup:``[4] <#_note-Varga:2000:jchemphys>`_`\ 
\ :sup:``[5] <#_note-Varga:1999:physreva>`_`\ 
\ :sup:``[6] <#_note-Mayer:1996:physreva>`_`\  indicate that
relativistic corrections to the exchange-correlation functionals have a
negligible effect on spectroscopic constants; however, for example for
core properties further studies are necessary.

Integration algorithm and thresholds
====================================

.SCREENING
----------

Specify the screening threshold used in the BLAS-3 integration scheme.

*Default:*

::

    .SCREENING
     1.0D-18

.POINTWISE
----------

Switch to old (unscreened) BLAS-2 integration scheme (default in
DIRAC10) instead of the new BLAS-3 scheme (default in DIRAC11).

.TINYDENS
---------

Specify the tiny density threshold. Grid points with a density smaller
than this threshold do not contribute to XC matrix elements.

*Default:*

::

    .TINYDENS
     1.0D-14

Asymptotic correction
=====================

.GRAC
-----

Activate the gradient regulated asymptotic correction scheme (GRAC)
\ :sup:``[7] <#_note-Gruning:2001:jchemphys>`_`\ , followed by two lines
of input.

First line: Asymptotic potential (currently LB94 or LBalpha).

Second line (free format): Parameters α and β (see reference
\ :sup:``[7] <#_note-Gruning:2001:jchemphys>`_`\ ), the ionization
potential, and the threshold for difference in HOMO eigenvalue below
which asymptotic correction is switched on. A sufficiently converged
density is needed before the correction is activated (before the bulk
potential is shifted in order to reproduce the desired IP (see reference
\ :sup:``[7] <#_note-Gruning:2001:jchemphys>`_`\ ).

*Typical input (here correcting PBE0):*

::

    **HAMILTONIAN
    .DFT
    PBE0
    *DFT
    .GRAC
    LB94
    0.5 40.0 0.79248 1.0D-6

.SAOP!
------

Activate the statistical averaging of (model) orbital potentials (SAOP)
as defined in reference
\ :sup:``[8] <#_note-Gritsenko:2000:jchemphys>`_`\ . This implies (and
activates) the ALDA kernel. The functional under
`.DFT </dirac/index.php/Manual:HAMILTONIAN#.DFT>`_ is expected to be
GLLBhole. The asymptotic potential is LBalpha.

*Recommended input:*

::

    **HAMILTONIAN
    .DFT
    GLLBhole
    *DFT
    .SAOP!

.SAOP
-----

The more general version of `.SAOP! <#.SAOP.21>`_

Followed by a single line of input: asymptotic potential (currently LB94
or LBalpha).

*Typical input (here correcting LDA):*

::

    **HAMILTONIAN
    .DFT
    LDA
    *DFT
    .SAOP
    LBalpha

Spin magnetization TDDFT
========================

.NOSDFT
-------

Turn off spin density contribution to XC response.

.COLLINEAR
----------

Use the collinear approximation as a definition of the spin density.

.. figure:: /dirac/images/math/7/b/b/7bba6e6e6ec5f4b9a6a6150321271ce6.png
   :align: center
   :alt:  s = m\_z = \\sum\_i \\phi\_i^\\dagger \\Sigma\_z \\phi\_i

    s = m\_z = \\sum\_i \\phi\_i^\\dagger \\Sigma\_z \\phi\_i
instead of the default noncollinear definition

.. figure:: /dirac/images/math/b/f/d/bfdef498df70433af2aca178b387f798.png
   :align: center
   :alt:  s = \\vert \\mathbf{m} \\vert

    s = \\vert \\mathbf{m} \\vert
.BETASIGMA
----------

use

.. figure:: /dirac/images/math/3/d/d/3ddf4ef5733771050baff36330159ef2.png
   :align: center
   :alt:  \\mathbf{m} = \\sum\_i \\phi\_i^\\dagger \\beta
   \\boldsymbol{\\Sigma} \\phi\_i

    \\mathbf{m} = \\sum\_i \\phi\_i^\\dagger \\beta
   \\boldsymbol{\\Sigma} \\phi\_i
instead of the default

.. figure:: /dirac/images/math/e/f/1/ef1e2146ae1ff966d11b2418e61a5346.png
   :align: center
   :alt:  \\mathbf{m} = \\sum\_i \\phi\_i^\\dagger \\boldsymbol{\\Sigma}
   \\phi\_i

    \\mathbf{m} = \\sum\_i \\phi\_i^\\dagger \\boldsymbol{\\Sigma}
   \\phi\_i
Adiabatic local density approximation (ALDA)
============================================

.ALDA
-----

Approximate all functional derivatives beyond the xc potential by SVWN
derivatives. For hybrid functionals exact exchange is switched off in
the solution of the response equation.

.XALDA
------

Same as `.ALDA <#.ALDA>`_ but keep the fraction of exact exchange of the
xc functional under `.DFT </dirac/index.php/Manual:HAMILTONIAN#.DFT>`_.

.ALDA+
------

Use `.ALDA <#.ALDA>`_ only for the Hermitian contribution (density
contribution) and use the proper xc kernel for the anti-Hermitian part
(spin density contribution).

.ALDA-
------

Use `.ALDA <#.ALDA>`_ only for the anti-Hermitian contribution (spin
density contribution) and use the proper xc kernel for the Hermitian
part (density contribution).

.XALDA+
-------

Use `.XALDA <#.XALDA>`_ only for the Hermitian contribution (density
contribution) and use the proper xc kernel for the anti-Hermitian part
(spin density contribution).

.XALDA-
-------

Use `.XALDA <#.XALDA>`_ only for the anti-Hermitian contribution (spin
density contribution) and use the proper xc kernel for the Hermitian
part (density contribution).

Other functionality
===================

`|image0| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.GAUNTSCALE
-----------

Scale Gaunt integrals (if included) with the same factor as for
Hartree-Fock exchange. This means that hybrid functionals include
fractional HF Gaunt interaction, and pure functionals no HF Gaunt
interaction at all. If this option is not given the Gaunt integrals will
be included to 100%, meaning that even an LDA calculation will include
full Hartree-Fock Gaunt interaction when the .GAUNT keyword is given in
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

References
----------

#. `↑ <#_ref-Schipper:2000:jchemphys_0>`_ P. R. T. Schipper, O. V.
   Gritsenko, S. J. A. van Gisbergen, and E. J. Baerends, *Molecular
   calculations of excitation energies and (hyper)polarizabilities with
   a statistical average of orbital model exchange-correlation
   potentials*, J. Chem. Phys. **112**, *1344* (2000) `electronic
   version <http://dx.doi.org/10.1063/1.480688>`_.
#. `↑ <#_ref-Saue:2002:jcomputchem_0>`_ T. Saue and T. Helgaker,
   *Four-component relativistic Kohn-Sham theory*, J. Comput. Chem.
   **23**, *814* (2002) `electronic
   version <http://dx.doi.org/10.1002/jcc.10066>`_.
#. `↑ <#_ref-Fossgaard:2003:jchemphys_0>`_ O. Fossgaard, O. Gropen, M.
   Corral Valero, and T. Saue, *On the performance of four-component
   relativistic density functional theory: Spectroscopic constants and
   dipole moments of the diatomics HX and XY (X, Y = F, Cl, Br and I)*,
   J. Chem. Phys. **118**, *10418* (2003) `electronic
   version <http://dx.doi.org/10.1063/1.1574317>`_.
#. `↑ <#_ref-Varga:2000:jchemphys_0>`_ S. Varga, B. Fricke, H.
   Nakaamatsu, T. Mukoyama, J. Anton, D. Geschke, A. Heitmann, E. Engel,
   and T. Bastug, *Four-component relativistic density functional
   calculations of heavy diatomic molecules*, J. Chem. Phys. **112**,
   *3499* (2000) `electronic
   version <http://dx.doi.org/10.1063/1.480934>`_.
#. `↑ <#_ref-Varga:1999:physreva_0>`_ S. Varga, E. Engel, W.-D. Sepp,
   and B. Fricke, *Systematic study of the Ib diatomic molecules
   Cu\ :sub:`2`\ , Ag\ :sub:`2`\ , and Au\ :sub:`2`\  using advanced
   relativistic density functionals*, Phys. Rev. A **59**, *4288* (1999)
   `electronic version <http://dx.doi.org/10.1103/PhysRevA.59.4288>`_.
#. `↑ <#_ref-Mayer:1996:physreva_0>`_ M. Mayer and O. D. Häberlen and N.
   Rösch, *Relevance of relativistic exchange-correlation functionals
   and of finite nuclei in molecular density-functional calulations*,
   Phys. Rev. A **54**, *4775* (54) `electronic
   version <http://dx.doi.org/10.1103/PhysRevA.54.4775>`_.
#. ↑ \ :sup:``7.0 <#_ref-Gruning:2001:jchemphys_0>`_`\ 
   \ :sup:``7.1 <#_ref-Gruning:2001:jchemphys_1>`_`\ 
   \ :sup:``7.2 <#_ref-Gruning:2001:jchemphys_2>`_`\  M. Grüning, O. V.
   Gritsenko, S. J. A. van Gisbergen, and E. J. Baerends, *Shape
   corrections to exchange-correlation potentials by gradient-regulated
   seamless connection of model potentials for inner and outer region*,
   J. Chem. Phys. **114**, *652* (2001) `electronic
   version <http://dx.doi.org/10.1063/1.1327260>`_.
#. `↑ <#_ref-Gritsenko:2000:jchemphys_0>`_ P. R. T. Schipper, O. V.
   Gritsenko, S. J. A. van Gisbergen, and E. J. Baerends, *Molecular
   calculations of excitation energies and (hyper)polarizabilities with
   a statistical average of orbital model exchange-correlation
   potentials*, J. Chem. Phys. **112**, *1344* (2000) `electronic
   version <http://dx.doi.org/10.1063/1.480688>`_.

