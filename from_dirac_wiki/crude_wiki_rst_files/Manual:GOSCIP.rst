Manual:GOSCIP
=============

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

*General Open Shell CI Program* written by Olivier Visser

This module allows for small full configuration interaction
calculations. It is invoked from within DIRRCI if the namelist
`GOSCIP </dirac/index.php/Manual:GOSCIP:GOSCIP>`_ is present in the
input, or if
`.RESOLVE </dirac/index.php/Manual:WAVE_FUNCTION#.RESOLVE>`_ is
specified when doing open shell HF calculations.

The input should be given in namelist form.

*Historical note:* This program was originally written for MOLFDIR
\ :sup:``[1] <#_note-prog:molfdir>`_`\  and was later also included in
DIRAC.

References
----------

#. `↑ <#_ref-prog:molfdir_0>`_ L. Visscher, O. Visser, H. Aerts, H.
   Merenga, and W. C Nieuwpoort, *Relativistic quantum chemistry: the
   MOLFDIR program package*, Comput. Phys. Comm. **81**, *120* (1994)
   `electronic
   version <http://dx.doi.org/10.1016/0010-4655(94)90115-5>`_.

