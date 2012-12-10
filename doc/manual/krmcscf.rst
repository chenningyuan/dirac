
star(KRMCSCF)

*2- and 4-component relativistic KR-MCSCF module*

Original implementation written
by Joern Thyssen, Timo Fleig and Hans Joergen Aagaard Jensen,
parallelization and continuous improvement by Stefan Knecht and Hans
Joergen Aagaard Jensen. Linear symmetry adaptation by Stefan Knecht and Hans Joergen Aagaard Jensen.

See works for more details on the theory and actual implementation [Jensen1996]_, [Thyssen2004]_ , [Thyssen2008]_, [Knecht2009]_.


**Mandatory keywords**
======================

keyword(CI PROGRAM)

specifies the CI module behind KRMCSCF, choices are LUCIAREL and GASCIP (*note*: GASCIP does not work with linear
symmetry). *default*:

::

    .CI PROGRAM
    GASCIP

keyword(CIROOTS)

number of states (roots) in symmetry X to optimize on. This keyword may be
repeated several times for multi-root multi-symmetry calculations
step-by-step. *default:* none.

example for 5 roots in (boson/fermion) irrep 1:

::

    .CIROOTS
    1 5

If you run the calculation in **linear symmetry** you have to specify the **2 x OMEGA** value
of the state to optimize on. The doubling (*2 x*) stems from the fact
that we want to avoid non-integer input, e.g. in case of an odd number
of electrons we might have :math:`\Omega` = 1/2, 3/2, 5/2, etc. values and the
corresponding input for one root with the :math:`\Omega` =1/2 would then read as

::

    .CIROOTS
    1 1

If we have a system with an even number of electrons and inversion
symmetry the input for two :math:`\Omega` =2g states read as

::

    .CIROOTS
    4g 2

keyword(INACTI)

Inactive orbitals per boson symmetry, separated by commas. This keyword is not allowed in connection with the citype
choice :ref:`CITYPE` GASCI. *default*: all orbitals active. 
example for D:math:`_{2h}`:

::

    .INACTI
     1,1,0,0,1,0,0,2


keyword(GAS SHELLS)

Specification of CI calculation with electron distribution in orbital (GAS) spaces. 
The first line contains the number of GA spaces to be used (1-7), 
followed by one line per GAS with a separation by a "/" of the min/max number of electrons in each GAS and the 
number of orbitals per fermion correp (either one (no inversion symmetry) or two (inversion
symmetry: *gerade* *ungerade*) entries per line). 

The first entry before the "/" gives the minimal number of accumulated (!) electrons
after consideration of this GAS, the second the corresponding maximum
number, separated by blanks. The minimum and maximum accumulated
occupations allow for a very flexible parameterization of the wave
function. All determinants fulfilling the occupation constraints will be
constructed. The second entry after the "/" gives  the 
number of orbitals per fermion correp. See also the open-shell input in :ref:`SCF` which is similar to the syntax used
here. The design of a GAS scheme is non-trivial and should be motivated by the electronic structure of the system (e.g.
inner core, outer core, valence, virtual space). Sometimes it is useful to subdivide the valence space, for scientific reasons, or/and the
virtual space, for technical reasons (save core memory). See reference [Fleig2006a]_,
pp. 27 for more details. example for 4 GAS spaces with an advanced excitation pattern:

::

    .GAS SHELLS
     4
      1  2 /  1  0
      6  8 /  1  2
     14 16 /  7  7
     16 16 / 30 24

keyword(INACTIVE)

Inactive orbitals per fermion correp. *default*: all orbitals active. example for 4 *gerade*  and 2 *ungerade* orbitals frozen in a molecule with inversion center:

::

    .INACTIVE
     4 2

**Optional keywords**
=====================

keyword(NOOCCN)

compute natural orbital occupation numbers for each electronic state. *default*: do not compute natural orbital occupation numbers.


keyword(ANALYZ)

analyze the final CI wave function printing the coefficients for each determinant above a given threshold
:math:`10^{-2}`. *default*: do not analyze the final CI wave function.


keyword(MAX CI)

maximum number of CI iterations. *default*:

::

    .MAX CI
     5

keyword(MXCIVE)

maximum size of Davidson subspace. *default*: 3 times the number of eigenstates (see :ref:`CIROOTS`) to optimize on.
example:

::

    .MXCIVE
     24

