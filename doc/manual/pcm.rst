star(PCM)

Polarizable continuum model (PCM) directives. 

General information on the PCM can be found in [Tomasi2005]_, whereas details 
concerning the Dirac implementation are found in [DiRemigio2014]_

**General**
===========

keyword(SKIPSS)

Performs a PCM-SCF calculation with the PCM-SCC approximation. The small-small block
of the electrostatic potential is neglected, in the spirit of the SCC approximation 
described in [Visscher1997a]_. It is not employed by default.
Currently not working for response calculations.

keyword(PRINT )

Print level in the PCM subroutines. Default::

  .PRINT
   0

**Advanced/Debug**
==================

keyword(SEPARA)

Performs PCM-SCF separating the nuclear and electronic electrostatic potential and apparent surface
charge. Results are unaffected. *For debug purpose only*

keyword(DOSPF )

Remove spin-orbit dependent part from PCM potential. *For debug purpose only*

keyword(SKIPOI)

Skips the calculation of the one-index transformed apparent surface charge in the
linear response function. *For debug purpose only* 

