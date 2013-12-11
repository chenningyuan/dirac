.. _reladc:

RELADC & LANCZOS
================

This section lists the available keywords for RELADC and the closely connected
iterative LANCZOS diagonalizer in the new input format maintained in DIRAC.

If you want to perform RELADC/LANCZOS calculations for the one- and two-particle
propagator you invoke both by setting the .RELADC keyword in the \*\*WAVE FUNCTION
Section.

The RELADC/LANCZOS calculation is then individually controlled in the
\*\*RELADC and \*\*LANCZOS input sections.

star(RELADC)

keyword(DOSIPS)

Do Single Ionization Potentials.
If you intend to do single ionization spectra calculations then set this keyword.
Closely related to SIP calculations is the keyword

keyword(ADCLEVEL)

Here you determine the perturbational order for a SIP calculation.

::

     ADCLEVEL = 1  strict second order
     ADCLEVEL = 2  extended second order
     ADCLEVEL = 3  third order + constant diagrams

*Default:*

::

    ADCLEVEL = 3 (including constant diagrams)

*Input example:*

::

   .ADCLEVEL
   2


keyword(SIPREPS)

Integer array of length 32. Specifies the symmetries of the one-hole
final states (SIP) to be calculated. If you do not enter this array
all symmetries are checked for possible final states and those are
calculated in the respective perturbational order. You can also specify
symmetries individually if you are interested only in a few by listing the
number of requested symmetries followed by the individual irrep numbers.

*Default:*

::

    SIPREPS(1:32) = 0  (no symmetries preselected, all are calculated)

*Input example:*

::

   .SIPREPS
   8
   1,3,5,7,17,19,21,23


keyword(READQKL)

Read back previously calculated constant diagrams in SIP runs. This is a restart
option to avoid the most time consuming step in SIP runs. The constant diagrams
for all symmetries are stored in the file QKLVAL.

keyword(NOCONS)

Block calculation of constant diagrams in third-order SIP calculations.
ADC is still executed up to third order but in the hole/hole block the (time consuming)
constant diagrams are omitted. Skipping constant diagrams reduces accuracy but
increases speed considerably. Not recommended for production runs.

keyword(VCONV)

Determines convergence of the inverse iteration in the calculation of
the constant diagrams. If these iterations take much time you can reduce
tightness by setting VCONV to 1.0E-05, 1.0E-04 asf. Check accuracy
of results when you activate this option.

*Default:*

::

    VCONV=1.0E-06

*Input example:*

::

    .VCONV
    1.0E-04

keyword(DOFULL)

Lanczos is automatically invoked in a single and double ionization calculation. If you want to 
invoke an additional full diagonalizer set this keyword. Attention: ADC matrices can become very large. A full
diagonalization is therefore supported only for matrix dimensions up to 5000 x 5000 and is useful
for test purposes only.

keyword(DOLANC)

Additionally invoke a Lanczos diagonalizer in the excitation calculation
(not recommended because Davidson is superior in this case). Lanczos is automatically invoked in SIP/DIP runs.

keyword(DODAVI)

Activate the Davidson diagonalizer in the excitation run. At the moment Davidson is only effective in combination 
with DOEXCI. If not specified, Davidson will not start in the excitation run.
In the next release the user has to activate the desired diagonalizer explicitly because the then available lifetime 
calculations do not require a complete diagonalization of the ADC matrix.

keyword(DODIPS)

Do Double Ionization Potentials.
If you intend to do double ionization spectra calculations then set this keyword.
The keyword .ADCLEVEL does not apply to DIP runs because the perturbational order
is set fixed to 'extended second order'.

*Default:*

::

    DODIPS = F

keyword(DIPREPS)

Integer array of length 32. Specifies the symmetries of the two-hole
final states (DIP) to be calculated. If you do not enter this array
all symmetries are checked for possible final states and those are
calculated in the respective perturbational order. You can also specify
symmetries individually if you are interested only in a few. Input is
analogous to SIPREPS.

*Default:*

::

    DIPREPS(1:32) = 0  (no symmetries preselected, all are calculated)

keyword(ADCTHR)

In DIP calculations only matrix elements whose amount is larger than ADCTHR
will be written to disk. This is due to the large matrices occurring in DIP
runs. Tests have shown that you can reduce matrix
size by a factor of two setting ADCTHR to 1.0E-05 with an accuracy loss
in the meV region. Perform accuracy tests for production runs since
the behavior is system-dependent.

*Default:*

::

    ADCTHR = 0.0 (all nonzero matrix elements are written to disk)

*Input example:*

::

    .ADCTHR
    1.0E-05

keyword(DOEXCI)

Activate calculation of excited energies and states using the four-component 
polarization propagator. An excited state has even number of electrons and is 
therefore classified according to the bosonic irreducible representations of the 
corresponding molecular point group. If only one or a few final states should be 
calculated use the following keyword:

keyword(EXCREPS)

Integer array of length 32. Specifies the symmetries of the excited
final states to be calculated. If you do not enter this array
all symmetries are checked for the occurrence of possible final states and those are
calculated in the respective perturbational order.

*Default:*

::

    EXCREPS(1:32) = 0  (no symmetries preselected, all are calculated)

*Input example:*

::

   .EXCREPS
   4
   1,3,17,19

keyword(DOXEXT)

Logical variable. If set true extended ADC(2) calculations will be performed (default = .false.).

star(LANCZOS)

Once the ADC matrices are stored in packed form on disk they are diagonalized
by the iterative Lanczos algorithm. The spectral information is written to
the files SSPEC.#irrep (SIPs) and DSPEC.#irrep (DIPs).
Hereby the ionization potential, the pole strength
and the error estimate are written in a line terminated by the '@' for
grep purposes. Immediately after this line follows the (indented)
configuration information belonging to this final state. This is
imaginable as a one-hole or two-hole Slater-determinant forming this
state in zeroth order. Strictly speaking, the configuration coefficients
refer to the intermediate state basis.

.SIPITER, .DIPITER, .EXCITER
----------------------------

Determines the number of Lanczos iterations in a SIP,
DIP or EXC calculation for all symmetries.
There is no need to specify the number of iterations
per symmetry because the convergence behaviour is similar within a
specific ionization class. However, DIP calculations can require
substantially more iterations for a comparable accuracy. Due to the
iterative nature of the Lanczos diagonalizer the edge values converge
very fast and some may be reproduced if SIPITER or DIPITER are set to
high values. These reproduced eigenvalues are spurious and will be
projected out from the final result. If one observes very many spurious
solutions (mainly in the SIP case) it is recommended to reduce SIPITER
accordingly. Note that Lanczos is not activated automatically in EXC runs.

*Default:*

::

    SIPITER = 500, DIPITER = 500, EXCITER = 500.


*Input example:*

::

    .SIPITER
    1000
    .DIPITER
    2500
    .EXCITER
    1500

.SIPPRNT, .DIPPRNT, .EXCPRNT
----------------------------

Real values. These two variables only control screen output of the
calculated eigenvalues and have no influence on the results in the
(SDX)SPEC files. You can enter the threshold in eV up to which computed
IPs (SIPs, DIPs or EXCs) will be printed on screen. Sometimes one is only
interested in a few lowest IPs and the screen output suffices.

*Default:*

::

    SIPPRNT = 50.0, DIPPRNT = 50.0, EXCPRNT = 20.0.


*Input example:*

::

    .SIPPRNT
    20.0
    .DIPPRNT
    100.0
    .EXCPRNT
    50.0

.SIPEIGV, .DIPEIGV, .EXCEIGV
----------------------------

For each chosen symmetry selected by the XXXREPS array a lower and upper energy boundary
value for the corresponding method are specified. 
Only within this energy range the long eigenvectors are calculated. This is more
user-friendly since one can not anticipate the number of eigenvectors to be expected in a certain energy range.


*Default:*

::

    0.0 0.0   (no eigenvectors calculated)

*Input example:*

::

    .SIPEIGV   #(same for .DIPEIGV and EXCEIGV)
    4 # number of lines of ranges to follow
    10.0 20.0
    20.0 30.0
    0.0 0.0
    10.0 15.0

keyword(DOINCORE)

This keyword activates an incore Lanczos diagonalization. If you run jobs on machines with
large core memory you can speed up the diagonalization considerably by transferring large parts
of the ADC matrix to memory. This is especially noticeable in DIP jobs because the matrices
are much larger than in the single ionization case.
Attention: Be aware that the operating system allows you to allocate as much memory as you want.
If there is not enough physical memory the OS starts to swap portions of the memory to disk.
You have to avoid this situation since your job will not terminate in time. Check carefully that
the memory you request is physically available! In later realizations the amount of physically
available memory will be checked by the program itself.

keyword(DVXXXXX)

The following keywords all start with DV and control the iterative Davidson diagonalizer effective in the
excitation module.

.DVROOTS
--------

This is the main keyword that tells the diagonalizer how many of the lowest roots you wish to calculate. 
Be aware that for higher roots convergence is harder to achieve and therefore requires more Davidson
iterations.

.DVMAXSP
--------

Here you define the maximum size of the Krylov space that is generated. In general, this number is ten to twenty 
times the number of requested roots in order to achieve convergence of approximately 1.0E-06 which is perfectly
sufficient. The algortihm performs so-called microiterations up to the size of DVMAXSP. Then one macroiteration
has terminated. The resulting vectors then serve as start vectors for the next set of microiterations leading
to considerably improved excited eigenstates. (Default = 50). The number of macroiterations is fixed by the keyword

.DVMAXIT
--------

This number normally ranges between 5 and 20 for harder cases. Default value is 5.

.DVCONV
-------

Here you determine the convergence of the eigenvalues (default = 1.0E-05). Sometimes a bit tighter
values are desired. You can choose that but be aware that very tight thresholds lead to many
iterations and sometimes to instabilities of the algorithm.
