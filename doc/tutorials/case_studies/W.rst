================
The Wolfram atom
================

Here we show how we get to converged ground state of the Wolfram atom (Z=74).
It is open-shell system with unpaired electrons in the 5d shell::

 W: [Xe] 4f(14) 5d(4) 6s(2)

Based on experience with d open-shell systems, the first step is to obtain converged state of the W(4+) cation. 
For that we have to ensure proper occupation of spinors in each fermion irreducible representation.::

 $PAM --noarch --mw=300 --mol=W.dyall_v2z.lsym.mol --inp=W_4+.x2c.scf_30_40.2fs.inp  --get "DFCOEF=DFCOEF.W_4+"

(Files to download ...)

We get the total SCF energy of -16146.075957236933 a.u., and, as a sign of correctness of SCF calculation, 
properly degenerated p,d eigenvalues in the analysis output. 

In the next step we reuse saved molecular spinors and keep overlap selection ::

 $PAM --noarch --mw=300 --mol=W.dyall_v2z.lsym.mol --inp=W.x2c.scf_ovlsel.2fs.inp --put "DFCOEF.W_4+=DFCOEF"

(Files to download ...)

One has to employ the dynamic overlap selection as the mere M_{J} selection does not work.
For the W atom, we get the total SCF energy of -16149.950403420618 a.u. and correct occupation
of d open-shell
::
                            1/2   3/2  5/2
 open shell #1 (f=0.4000):    2    2    1

This can be verified as the 
d3/2 shell gives M_J=1/2,3/2 and the d5/2 has M_J=(1/2, 3/2,5/2), what makes the above shown open-shell occupation.

--------------------------
How about lower symmetry ?
--------------------------

We descend to the C2v symmetry with 1 ferion irrep. Again, we start from the closed shell W(4+) cation and use DFCOEF for the subsequent run.

First, the W(4+) cation does not have proper degeneracy of spinors, though it converges. For the subsequent, "spinors reordering" run one has to place
restrictions of the order of MO and get DFCOEF saved.

Unfortunately, calculations in C2v symmetry do not lead to  properly ordered and proprely degenerated orbitals of  W(4+) ,so no success in the neutral W atom.
Maybe more work with pre-SCF reordering of starting W(4+) spinors would fix this issue.

(Files to download ...)
