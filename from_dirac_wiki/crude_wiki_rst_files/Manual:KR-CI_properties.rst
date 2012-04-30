Manual:KR-CI properties
=======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

*2- and 4-component relativistic KR-CI property module* written by
Stefan Knecht and Hans Joergen Aa. Jensen, parallelization by Stefan
Knecht

The KR-CI property module \ :sup:``[1] <#_note-knecht:thesis>`_`\  takes
advantage of the 2- and 4-component KR-CI module. It can be used for the
computation of:

::

    - permanent dipole moments in ground and excited states, 

::

    - transition dipole moments 

::

    - computation of ,  and  (Ω values) expectation values (only for linear molecules).

Upon request the analysis of other one-electron properties may also be
implemented. The module could in principle be used for each one-electron
operator that is specified in
`[1] <http://wiki.chem.vu.nl/dirac/index.php/Manual:One-electron_operators>`_.

**NOTE: the following keywords must be placed under the input deck**

::

    *KRCICALC

together with the appropriate CI calculation input (GASSH, GASSPC, etc.)

Contents
--------

-  `1 implemented keywords <#implemented_keywords>`_

   -  `1.1 .DIPMOM <#.DIPMOM>`_
   -  `1.2 .TRDM <#.TRDM>`_
   -  `1.3 .OMEGAQ <#.OMEGAQ>`_

-  `2 References <#References>`_

implemented keywords
====================

.DIPMOM
-------

compute the permanent dipole moments in electronic ground and excited
states.

.TRDM
-----

compute the transition dipole moments between electronic states.

.OMEGAQ
-------

compute the expectation values of the spin- and angular momentum
operator in z-direction and print the total |<\\hat{j}\_z> =
<\\hat{s}\_z> + <\\hat{l}\_z>| expectation value (**Ω** value) for each
electronic state.

::

    Note: this applies only to linear molecules.

References
==========

#. `↑ <#_ref-knecht:thesis_0>`_ S. R. Knecht, *PhD thesis*, University
   of Duesseldorf (Germany) '***, '** (2009)* .

