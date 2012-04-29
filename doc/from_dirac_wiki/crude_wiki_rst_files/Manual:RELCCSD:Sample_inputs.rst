Manual:RELCCSD:Sample inputs
============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Calculate the CCSD(T) energy and MP2 first order properties. The full
calculation is given as input example 1.energy+dipole in the test
directory:

::

     &RELCCSD TIMING=T, DOENER=T, DOFOPR=T &END
     &CCFOPR  DOMP2G=T &END

