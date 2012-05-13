SCF sample input
==================

*Example:* Minimum wave function input for water (the
*C*\ :sub:`2*v*`\  point group):

::

    **WAVE FUNCTION
    .SCF
    *SCF
    .CLOSED SHELL
     10

For symmetry groups *with* inversion symmetry, you need to specify
the distribution of the electrons in the two fermion irreps [Saue2000]_.

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

*Example:*

::

    .OPEN SHELL
     1
     5/0,6

1 open shell with 5 electrons in 6 spinors (= 3 Kramers pairs) in
irrep 2. Thus, the fractional occupation is 5/6.

The open shell module in DIRAC is based on
average-of-configurations [Thyssen1998]_ .
The simplest case is one electron in two spinors (= one Kramers
pair). For this special case the average-of-configuration
calculation gives the same result as the usual restricted
open-shell Hartree-Fock. For all other cases the calculation gives
the average energy of many states.

An example of one electron in two spinors is the Boron atom:

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
electrons (in the 1*s* and 2*s* orbitals) and one open shell, where
1 electron is in two spinors (the 2*p*\ :sub:`1 / 2`\  and
2\\bar{p}\_{1/2}| spinors) in irrep 2. This calculation would give
you the *P*\ :sub:`1 / 2`\  energy of Boron. We could also do a
calculation with

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

   :align: center
   :alt:

\\frac{2\\cdot E\_{P\_{1/2}} + 4\\cdot E\_{P\_{3/2}}}{6}

\\frac{2\\cdot E\_{P\_{1/2}} + 4\\cdot E\_{P\_{3/2}}}{6}

It is also possible to specify more than one open shell. Let us
take a look at two different inputs for the open shells of the
Platinum atom (6*s*\ :sup:`1`\ 5*d*\ :sup:`9`\ ):

In the first example:

::

    .OPEN SHELL
     2
     9/10,0
     1/2,0

we get the average energy of 20 states:

\\left(\\begin{array}{c} 10 \\\\ 9 \\end{array}\\right)\\cdot
\\left(\\begin{array}{c} 2 \\\\ 1 \\end{array}\\right)\\cdot=10\\cdot 2=20|.

All 20 states are *s*\ :sup:`1`\ *d*\ :sup:`9`\  states.

In the second input:

::

    .OPEN SHELL
     1
     10/12,0

we distribute the same ten electrons in all 12 spinors. This gives
a total of
\\left(\\begin{array}{c} 12 \\\\ 10 \\end{array}\\right) = 66|
states, which are a mixture of one *s*\ :sup:`0`\ *d*\ :sup:`10`\
state, 20 *s*\ :sup:`1`\ *d*\ :sup:`9`\  and 45
*s*\ :sup:`2`\ *d*\ :sup:`8`\  states.

