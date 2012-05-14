GOSCIP:Sample inputs
====================

Calculate the spin-orbit components of the 2P ground state of the
fluorine atom using a COSCI wavefunction. This will represent the
"uncorrelated" result since the CI will merely give the energies of the
individual determinants in this simple case. The same calculation can
also be done automatically by specifying
`.RESOLVE </dirac/index.php/Manual:WAVE_FUNCTION:RESOLVE>`_ in the `WAVE
FUNCTION </dirac/index.php/Manual:WAVE_FUNCTION>`_ section of the input.

The full calculation is given as input example 6.cosci.energy REFERENCE in the
test directory:

::

     &GOSCIP NELEC=5, IPRNT=1 &END
     &POPANA THRESH=1.0D-4, DEGEN=1.0D-12, SELPOP=50.0 &END

**NOTE: In the developer's version the first line reads as**

::

     &GOSCIP NELACT=5, IPRNT=1 &END
     &POPANA THRESH=1.0D-4, DEGEN=1.0D-12, SELPOP=50.0 &END


