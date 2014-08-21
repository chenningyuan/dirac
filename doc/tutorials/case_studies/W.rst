================
The Wolfram atom
================

Here we show how we get to the converged ground state of the Wolfram atom (Z=74).
It is open-shell system with four unpaired electrons in the *5d* shell::

 W: [Xe] 4f(14) 5d(4) 6s(2)

----------------------------
Starting from :math:`W^{4+}`
----------------------------

Based on experience with *d* open-shell systems, the first step is to obtain converged state of the (closest) closed-shell cation, :math:`W^{4+}`.
For that we have to ensure proper occupation of spinors in each of two fermion irreducible representation:

::

 $PAM --noarch --mw=300 --mol=W.dyall_v2z.lsym.mol --inp=W_4+.x2c.scf_30_40.2fs.inp  --get "DFCOEF=DFCOEF.W_4+"

(Files to download are :download:`W.dyall_v2z.lsym.mol <../../../../test/tutorial_W/W.dyall_v2z.lsym.mol>` and
:download:`W_4+.x2c.scf_30_40.2fs.inp <../../../../test/tutorial_W/W_4+.x2c.scf_30_40.2fs.inp>` )

The total SCF energy of this cation is -16146.075957236933 a.u., and, as a sign of correctness of the SCF calculation, 
there are properly degenerated *p*,*d* etc eigenvalues in the analysis output 
(see the full output file :download:`W_4+.x2c.scf_30_40.2fs_W.dyall_v2z.lsym.out <../../../../test/tutorial_W/result/W_4+.x2c.scf_30_40.2fs_W.dyall_v2z.lsym.out>`). 
In the case of ill SCF process, though converged, atomic shells are not degenerated as they should be.

------------------------------
The W atom from :math:`W^{4+}`
------------------------------

In the next step we reuse saved molecular spinors from file *DFCOEF.W_4+* and keep the overlap selection ::

 $PAM --noarch --mw=300 --mol=W.dyall_v2z.lsym.mol --inp=W.x2c.scf_ovlsel.2fs.inp --put "DFCOEF.W_4+=DFCOEF"

(File to download :download:`W.x2c.scf_ovlsel.2fs.inp <../../../../../test/tutorial_W/W.x2c.scf_ovlsel.2fs.inp>`)

One has to take the dynamic overlap selection as the mere :math:`M_{J}` selection does not work 
(you can try with input file :download:`W.x2c.scf_Mjsel.inp <../../../../test/tutorial_W/W.x2c.scf_Mjsel.inp>`).

For the W atom, we get the total SCF energy of -16149.950403420618 a.u. and, for example, correct occupation
of the *5d* shell, which lies above the *6s* doubly occupied level ::

                             1/2   3/2  5/2  
  open shell #1 (f=0.4000):    2    2    1

This can be easily verified as the 
:math:`d_{3/2}` shell gives values of :math:`M_{J}=\frac{1}{2},\frac{3}{2}` 
and the :math:`d_{5/2}` spans values of :math:`M_{J}=\frac{1}{2},\frac{3}{2},\frac{5}{2}`, 
what makes the  *5d* open-shell occupation, 2 2 1, as shown above.
See also the complete output file, :download:`W.x2c.scf_ovlsel.2fs_W.dyall_v2z.lsym.out <../../../../test/tutorial_W/result/W.x2c.scf_ovlsel.2fs_W.dyall_v2z.lsym.out>`


--------------------------
How about lower symmetry ?
--------------------------

We descend to the C2v symmetry having 1 fermion irreducible representation. Again, we start from the closed shell cation, :math:`W^{4+}`, 
and use its DFCOEF for the next run.
First we discover, that the :math:`W^{4+}` system does not have proper degeneracy of (l>0) spinors in the C2v symmetry, though it is converging. 
For the subsequent, spinors reordering run, one has to place restrictions on the order of MO and get its DFCOEF saved.

Unfortunately, subsequent calculations in the C2v symmetry do not lead to properly ordered and proprely degenerated orbitals of :math:`W^{4+}`,
therefore (at least we assume), no success can be achieved with the neutral W atom.
Maybe some more test with various pre-SCF reordering of starting :math:`W^{4+}` spinors would fix this issue ?

(Input files of the :math:`C_{2v}` symmetry to download are :download:`W.dyall_v2z.C2v.mol <../../../../test/tutorial_W/W.dyall_v2z.C2v.mol>`, 
:download:`W_4+.x2c.scf.1fs.inp <../../../../test/tutorial_W/W_4+.x2c.scf.1fs.inp>`,
:download:`W.x2c.scf_ovlsel.1fs.inp <../../../../test/tutorial_W/W.x2c.scf_ovlsel.1fs.inp>`,
:download:`W_4+.x2c.scf_reorder_ovlsel.1fs.inp <../../../../test/tutorial_W/W_4+.x2c.scf_reorder_ovlsel.1fs.inp>`)

