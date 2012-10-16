\*\*GENERAL -- Sample inputs
============================

Use a reduced value for the speed of light to amplify the effect of
relativity (for testing):

::

    **GENERAL
    .CVALUE
     50.0

Contrary you may want to extrapolate to the nonrelativistic limit by setting
speed of light to very large value. You have to experiment a little to find out appropriate number
as this depends on the character of the studied system.

Very often one needs to export or to read formatted file with orbital
coefficients - DFPCMO. Therefore the user should provide

::

    **GENERAL
    .PCMOUT

for export of the molecular orbitals (spinors).


\*SCF -- Sample inputs
========================

*Example:* Minimum wave function input for water (the :math:`C_{2v}` point group):

::

    **WAVE FUNCTION
    .SCF
    *SCF
    .CLOSED SHELL
     10

For symmetry groups *with* inversion symmetry, you need to specify
the distribution of the electrons in the two fermion irreps, see Ref. [Saue2000]_.

*Example:* Minimum wave function input for the Neon atom:

::

    **WAVE FUNCTION
    .SCF
    *SCF
    .CLOSED SHELL
     4 6


Simple SCF input with 10 electrons in *Eg* irrep, with printing of positive eigenvalues only.

::

    *SCF
    .PRINT
     1
    .NELECT
     10 0
    .MAXITR
     30
    .EVCCNV
     1.0D-6
    .MXDIIS
     10
    .CNVINT
     1.0D-2 1.0D-4
    .EIGPRI
     1 0

Open shell example - Boron atom
-------------------------------

The simplest, and in most cases converging open-shell configuration is through the average-of-occupation means.
As an example:

::

    .OPEN SHELL
     1
     5/0,6

There is 1 open shell filled with 5 electrons distributed over 6 *p* spinors (= 3 Kramers pairs) in the
irrep 2. Thus, the fractional occupation is 5/6.

The open shell module in DIRAC is based on
average-of-configurations, see Ref. [Thyssen1998]_ .

The simplest case is one electron in two spinors (= one Kramers
pair). For this special case the average-of-configuration
calculation gives the same result as the usual restricted
open-shell Hartree-Fock. For all other cases the calculation gives
the average energy of many states.

Let us try to build up occupation of the Boron atom.
This is  example of one electron in two spinors:

::

 **WAVE FUNCTION
    .SCF
    *SCF
    .CLOSED SHELL
     4 0
    .OPEN SHELL
     1
     1/0,2

The input is interpreted in the following way: We have 4 inactive
electrons (in the *1s* and *2s* orbitals) and one open shell, where
1 electron is in two spinors (the :math:`2p_{1/2}`  and :math:`2\bar{p}_{1/2}` spinors) in irrep 2.
This calculation would give us the energy of the :math:`X P_{1/2}` state of the Boron atom.

Alternatively, we could also do a calculation with occupation of

::

    **WAVE FUNCTION
    .SCF
    *SCF
    .CLOSED SHELL
     4 0
    .OPEN SHELL
     1 
     1/0,6

where we have the active electrons distributed in all six
*p*-orbitals. This would give us the average energy of six states:

.. math::

 \frac{2 \cdot E_{P_{1/2}} + 4 \cdot E_{P_{3/2}}}{6}

Open shell example -- Platinum atom
------------------------------------

It is also possible to specify more than one open shell. Let us
take a look at two different inputs for the open shells of the
Platinum atom (:math:`6s^{1} 5d^{9}`):

In the first example:

::

    .OPEN SHELL
     2
     9/10,0
     1/2,0

we get the average energy of 20 possible electronic states:

.. math::
 \left(\begin{array}{c}
   10  \\
   9 
 \end{array}\right)
   \cdot
 \left(\begin{array}{c} 
   2 \\ 
   1 
 \end{array}\right) 
 \cdot = 10 \cdot 2 = 20.

All 20 states are arising from the :math:`s^{1} d^{9}` configuration.

In the second input:

::

    .OPEN SHELL
     1
     10/12,0

we distribute the same ten electrons  over all 12 spinors. This gives
a total of :math:`\left(\begin{array}{c}{12}\\{10}\end{array}\right)=66` 
states, which are a mixture of one *s*\ :sup:`0`\ *d*\ :sup:`10`\
state, 20 *s*\ :sup:`1`\ *d*\ :sup:`9`\  and 45
*s*\ :sup:`2`\ *d*\ :sup:`8`\  states.
