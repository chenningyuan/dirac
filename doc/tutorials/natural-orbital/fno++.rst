:orphan:
 

================================
FNO++ Natural Orbital tutorial
================================

Introduction
============

The reader is advised to read the MP2NO tutorial before reading this one. MP2NO tend to perform pretty badly for recovering response properties such as the polarizability. FNO++ :cite:`Chakraborty2025` is another variant of NO that have been benchmarked on linear response properties.

FNO++ in Dirac
===============

With the water molecule as system of study:

.. literalinclude:: h2o.xyz

The input for generation the NO is:

.. literalinclude:: fno++_generation_nr_mean.inp

The new keywords are `.MP2NO++` with its truncation threshold below, `.NO_PROP` with first the number of operators to involve and followed by the operators index defined in `MOLTRA`. Finally, the `.NO_WEIGHTED_SUM` keyword specifies to take the averaged sum of the previously defined operators. It is recommended to not experiment with `.NO_WEIGHTED_SUM` equals 0. Then u can reuse the generated orbitals in the same fashion as the MP2NO tutorial.

FNO++ with frequency in Dirac
==============================

For computing frequency dependent polarizability, you can specify the keyword `.NO_OMEGA` followed by real and imaginary part of frequency of the FNO++:

.. code-block:: bash

   .NO_OMEGA
    1.0E-2 0.0

