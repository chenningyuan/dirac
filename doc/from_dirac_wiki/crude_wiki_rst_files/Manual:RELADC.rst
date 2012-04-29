Manual:RELADC
=============

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

*Relativistic Green's function (propagator) module RELADC* written by
Markus Pernpointner.

Propagators (Green's functions) are a very useful tool for the
calculation of ionization energies, electron affinities, excitation
energies, double ionization energies and so forth. The mathematical form
of the propagator reveals the type of physical property which is
described but does not necessarily point to an efficient solution
strategy. One way to achieve this is the method of (A)lgebraic
(D)iagrammatic (C)onstruction originally formulated by Schirmer
\ :sup:``[1] <#_note-schirmer:adc>`_`\ . Hereby the diagonal
representation of the propagator is transformed to a nondiagonal form by
formally inserting a completeness relation in an arbitrary basis which
does not have to be specified further because this is just a formal
step. But now the nondiagonal form is expanded into orders of
perturbation theory and compared to the known diagrammatic expansion of
the propagator order by order. By this procedure the matrix entries of
the nondiagonal expansion can be determined and calculated resulting in
the ADC matrix. The final diagonalization of the ADC matrix then yields
the desired eigenvalues (poles) and eigenvectors (spectral intensities).

In the current version of DIRAC a four-component realization of ADC
providing single \ :sup:``[2] <#_note-pernpointner:dhfadc2>`_`\  and
double \ :sup:``[3] <#_note-pernpointner:dhfadc3>`_`\  ionization
potentials (SIPs and DIPs) is available starting from a \*closed shell\*
reference determinant. Open-shell ADC based on a single reference (e.g.
high-spin states) can be done in principle but the theoretical basis is
not very sound for this. The four-component ADC is based on the
Dirac-Coulomb operator and officially termed as DC-ADC. Be aware that
for the moment RELADC has to be compiled with the -i8 compiler flag that
means 64bit Integers are required. This will soon be abandoned.

The following perturbational orders are implemented in DC-ADC:

**SIPs:** strict second order (ADC-2), extended second order (ADC-2e)
and third order (ADC-3); **DIPs:** extended second order (ADC-2e).

Only little input is needed for program control and is given via a
&RELADC &END Namelist directive. It should be mentioned that as in any
method treating electron correlation the active orbital space also in
this case determines the quality of the results. Since after
construction of the ADC matrix a final diagonalization is required for
obtaining the total spectrum (in the space of active occupied orbitals,
of course) a Lanczos diagonalizer is invoked as well. Controlling the
diagonalizer also occurs via Namelist starting with &LANINP &END.

RELADC similar to many post-DHF programs sets out from the transformed
totally antisymmetric two-electron integrals |\\langle
ij\|\|kl\\rangle|. Therefore the
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ step is required and
automatically invoked if the .4INDEX keyword is not given.

Contents
--------

-  `1 &RELADC &END -- Control propagator
   calculation <#.26RELADC_.26END_--_Control_propagator_calculation>`_

   -  `1.1 DOSIPS,DODIPS <#DOSIPS.2CDODIPS>`_
   -  `1.2 ADCLEVEL <#ADCLEVEL>`_
   -  `1.3 XREPS <#XREPS>`_
   -  `1.4 READQKL <#READQKL>`_
   -  `1.5 DOCONST <#DOCONST>`_
   -  `1.6 VCONV <#VCONV>`_

-  `2 &LANINP &END -- Control Lanczos
   diagonalizer <#.26LANINP_.26END_--_Control_Lanczos_diagonalizer>`_

   -  `2.1 SIPITER,DIPITER <#SIPITER.2CDIPITER>`_
   -  `2.2 SIPPRNT,DIPPRNT <#SIPPRNT.2CDIPPRNT>`_
   -  `2.3 SIPEIGV <#SIPEIGV>`_

-  `3 Examples <#Examples>`_
-  `4 References <#References>`_

&RELADC &END -- Control propagator calculation
==============================================

DOSIPS,DODIPS
-------------

Logical variables. These invoke the main tasks of calculating Single
Ionization Potentials via DOSIPS=T and/or Double ionization potentials
via DODIPS=T. Either of them can be specified individually or
simultaneously.

*Default:*

::

    DOSIPS = F, DODIPS = F

ADCLEVEL
--------

Integer variable. Specifies the perturbational order for the SIP
calculation only. Remember: the order for DIPs is always extended second
order.

::

     ADCLEVEL = 1  strict second order
     ADCLEVEL = 2  extended second order 
     ADCLEVEL = 3  third order + constant diagrams

*Default:*

::

    ADCLEVEL = 3

XREPS
-----

Integer array of length 32. Specifies the symmetries of the one-hole
final states (SIP) to be calculated. Individual specification is useful
if one is interested in a few symmetries only. However, for the double
ionizations (DIPs) the determination of final state symmetries is more
complicated and done automatically by the program. Therefore the XREPS
directive has no influence on DIP calculations and all DIP final states
are calculated in whatever symmetry they may occur.

*Default:*

::

    XREPS(*) = 0  no symmetries preselected.

READQKL
-------

Logical variable. If the constant diagrams were calculated in a previous
ADC(3) run (time consuming step) they can be read back in via the file
QKLVAL. If this file does not exist an automatic recalculation is
performed.

*Default:*

::

    READQKL = F

DOCONST
-------

Logical variable. Turns off the calculation of constant diagrams in an
ADC(3) run (if DOCONST = F). All other perturbational orders of the
matrix blocks are retained. This option is be useful for the precise
analysis of individual contributions.

*Default:*

::

    DOCONST = T

VCONV
-----

Double precision variable. Determines convergence of the inverse
iteration in the calculation of the constant diagrams. The threshold is
set to 1.0E-06 being sufficent regarding the overall methodological
accuracy.

*Default:*

::

    VCONV=1.0E-06

&LANINP &END -- Control Lanczos diagonalizer
============================================

This is a second Namelist controlling the Lanczos diagonalizer. The
spectral information is written to the files SSPEC.#irrep (SIPs) and
DSPEC.#irrep (DIPs). Hereby the ionization potential, the pole strength
and the error estimate are written in a line terminated by the '@' for
grep purposes. Immediately after this line follows the (indented)
configuration information belonging to this final state. This is
imaginable as the one-hole or two-hole Slater-determinant forming this
state in zeroth order. The suffix #irrep corresponds to the final state
symmetry. For each symmetry calculated one therefore gets one (DS)SPEC
file.

SIPITER,DIPITER
---------------

Integer variables. Determines the number of Lanczos iterations in a SIP
or DIP calculation. There is no need to specify the number of iterations
per symmetry because the convergence behaviour is similar within a
specific ionization class. However, DIP calculations can require
substantially more iterations for a comparable accuracy. Due to the
iterative nature of the Lanczos diagonalizer the edge values converge
very fast and some may be reproduced if SIPITER or DIPITER are set to
high values. These reproduced eigenvalues are spurious and will be
projected out from the final result. If one observes very many spurious
solutions (mainly in the SIP case) it is recommended to reduce SIPITER
accordingly.

*Default:*

::

    SIPITER = 500, DIPITER = 500.

SIPPRNT,DIPPRNT
---------------

Real variables. These two variables only control *Screen output* of the
calculated eigenvalues and have no influence on the results in the
(DS)SPEC files! You can enter the threshold in eV up to which computed
IPs (SIPs or DIPs) will be printed on screen. Sometimes one is only
interested in a few lowest IPs and the screen output suffices.

*Default:*

::

    SIPPRNT = 50.0, DIPPRNT = 50.0.

SIPEIGV
-------

Integer array of length 32 with the same structure as XREPS. For each
chosen symmetry in XREPS you specify here the number of requested *long*
eigenvectors corresponding to the ADC matrix.

The following files are created by activating long eigenvector
calculation: a) LONGEVC.#irrep: it contains the number of requested
eigenvectors in binary format (real or complex) for a specific irrep and
is used in subsequent lifetime calculations. b) the file EVCANAL.#irrep:
it is generated in human readable form and contains the same information
as LONGEVC.#irrep but with all the contributing configurations (main and
satellite contributions to a specific final state). The creation of long
eigenvectors invokes a second set of Lanczos iterations. The usual hole
analysis and spectrum generation is not affected and long eigenvectors
are only needed if one wants to analyze the satellite structure of a
specific singly ionized final state or for lifetime calculations of the
ionized state (coming soon).

*Default:*

::

    SIPEIGV(*)=0

`Examples </dirac/index.php/Manual:RELADC:Sample_inputs>`_
==========================================================

References
==========

#. `↑ <#_ref-schirmer:adc_0>`_ J. Schirmer and L. S. Cederbaum and O.
   Walter, *New approach to the one-particle Green's function for finite
   Fermi systems*, Phys. Rev. A **28**, *1237* (1983) .
#. `↑ <#_ref-pernpointner:dhfadc2_0>`_ M. Pernpointner, *The
   one-particle Green's function method in the Dirac–Hartree–Fock
   framework. II. Third-order valence ionization energies of the noble
   gases, CO and ICN*, J. Chem. Phys. **121**, *8782* (2004) .
#. `↑ <#_ref-pernpointner:dhfadc3_0>`_ M. Pernpointner, *The
   four-component two-particle propagator for the calculation of double
   ionization spectra of heavy-element compounds I. Method.*, J. Phys. B
   **43**, *205102* (2010) .

