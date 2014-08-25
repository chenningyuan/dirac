==========================================
Getting excited states of :math:`Ir^{16+}`
==========================================

In the following we are interested in getting the ground (:math:`^{2}S_{1/2}`) and two closest excited states 
(:math:`^{2}P_{1/2}`, :math:`^{2}P_{3/2}`) of the :math:`Ir^{16+}` cation and their correlation energies.
The elctronic configurations of these states are 
::

 (^2S)    : [Xe] 4f(14) 5s(1) 5p1/2(0)  5p3/2(0)  
 (^2P1/2) : [Xe] 4f(14) 5s(0) 5p1/2(1)  5p3/2(0)  
 (^2P3/2) : [Xe] 4f(14) 5s(0) 5p1/2(0)  5p3/2(1) 


For simplicity, we will work with the wo-component Hamiltonian (:ref:`.X2C`) and employ the smallest *v2z* decontracted basis set by K.Dyall.
Due to the convergence problem of the standalone  AMFI atomic SCF code we keep the +2 charge (:ref:`.AMFICH`) for mean-field orbitals.


There are two ways to obtain excited states -  (i) from the converged SCF state, and, (ii) only at the correlated level from the 2P_aver SCF state.

For subsequent Coupled Cluster (CC) correlated calculations please soften the DHOLU value in the subroutine DENOM (file *src/relccsd/cceqns.F*) to 5.0D-4 and
recompile DIRAC. However, this is not recommended approach as the *p32* states are in general not well described at the CC level.


The :math:`^{2}S_{1/2}` ground state
------------------------------------

Thanks to the linear symmetry, having two irreps, one can place the unpaired electron the 1st irrep to get the :math:`^{2}S_{1/2}` ground state.
SCF calculations are followed by two CC calculations, where in the second one uses orbital energies for denominators.
::

 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2S12.scf_cc33e.2fs.inp
 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2S12.scf_cc33e_oe.2fs.inp

Input files to download are
:download:`Ir.dyall_v2z.lsym.mol <../../../../test/tutorial_Ir_16plus/Ir.dyall_v2z.lsym.mol>`, 
:download:`Z61.x2c.2S12.scf_cc33e.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2S12.scf_cc33e.2fs.inp>`,
:download:`Z61.x2c.2S12.scf_cc33e_oe.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2S12.scf_cc33e_oe.2fs.inp>`.
Corresponding output files are 
:download:`Z61.x2c.2S12.scf_cc33e.2fs_Ir.dyall_v2z.lsym.out <../../../../test/tutorial_Ir_16plus/result/Z61.x2c.2S12.scf_cc33e.2fs_Ir.dyall_v2z.lsym.out>`, 
:download:`Z61.x2c.2S12.scf_cc33e_oe.2fs_Ir.dyall_v2z.lsym.out <../../../../test/tutorial_Ir_16plus/result/Z61.x2c.2S12.scf_cc33e_oe.2fs_Ir.dyall_v2z.lsym.out>`.


The SCF :math:`^{2}P_{1/2}`, :math:`^{2}P_{3/2}` excited states
----------------------------------------------------------------

By placing the unpaired electron into 2nd irrep one gets the :math:`^{2}P_{1/2}` first excited state:
::

 pam  --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2P12.scf_cc33e.2fs.inp --get "DFCOEF=DFCOEF.v2z.2P12"
 pam  --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2P12.scf_cc33e_oe.2fs.inp 

Input files to download
:download:`Z61.x2c.2P12.scf_cc33e.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2P12.scf_cc33e.2fs.inp>`,
:download:`Z61.x2c.2P12.scf_cc33e_oe.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2P12.scf_cc33e_oe.2fs.inp>`.
Corresponding output files are
:download:`Z61.x2c.2P12.scf_cc33e.2fs_Ir.dyall_v2z.lsym.out <../../../../test/tutorial_Ir_16plus/result/Z61.x2c.2P12.scf_cc33e.2fs_Ir.dyall_v2z.lsym.out>`, 
:download:`Z61.x2c.2P12.scf_cc33e_oe.2fs_Ir.dyall_v2z.lsym.out <../../../../test/tutorial_Ir_16plus/result/Z61.x2c.2P12.scf_cc33e_oe.2fs_Ir.dyall_v2z.lsym.out>`.

How to obtain the :math:`^{2}P_{3/2}` second excited state at the SCF level, and, consequently, at the RelCC level ? 
For that, we utilize the :ref:`.REORDER MO` keyword with readin DFCOEF.v2z.2P12 from the previous step.
Also one has use overlap selection in the SCF step:
::

  pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2P32.scf_cc33e.2fs.inp --put "DFCOEF.v2z.2P12=DFCOEF"
  pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2P32.scf_cc33e_oe.2fs.inp --put "DFCOEF.v2z.2P12=DFCOEF"

Corresponding input files to download are
:download:`Z61.x2c.2P32.scf_cc33e.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2P32.scf_cc33e.2fs.inp>`,
:download:`Z61.x2c.2P32.scf_cc33e_oe.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2P32.scf_cc33e_oe.2fs.inp>`.
Output files are 
:download:`Z61.x2c.2P12.scf_cc33e.2fs_Ir.dyall_v2z.lsym.out <../../../../test/tutorial_Ir_16plus/result/Z61.x2c.2P12.scf_cc33e.2fs_Ir.dyall_v2z.lsym.out>`, 


The CCSD(T) :math:`^{2}P_{1/2}`, :math:`^{2}P_{3/2}` excited states
-------------------------------------------------------------------

The other option is to start from the 2^P averaged single determinant state and distinguish 
between individual  ^2P1/2,^2P3/2 states at the Coupled Cluster correlated level thanks to the linear symmetry.

First we test the pure 2Paveraged state:
::

 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e.2fs.inp
 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e_oe.2fs.inp

Files to download are
:download:`Z61.x2c.2Paver.scf_cc33e.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e.2fs.inp>`,
:download:`Z61.x2c.2Paver.scf_cc33e_oe.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e_oe.2fs.inp>`.

And afterwards we proceed to the individual states, based on M_J splitted occupation of each fermion irrep.
First the 2P12 state:
::

 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e_2P12.2fs.inp
 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e_oe_2P12.2fs.inp

Input files to download are
:download:`Z61.x2c.2Paver.scf_cc33e_2P12.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e_2P12.2fs.inp>`,
:download:`Z61.x2c.2Paver.scf_cc33e_oe_2P12.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e_oe_2P12.2fs.inp>`.
Corresponding output files
:download:`Z61.x2c.2Paver.scf_cc33e_2P12.2fs_Ir.dyall_v2z.lsym.out <../../../../test/tutorial_Ir_16plus/result/Z61.x2c.2Paver.scf_cc33e_2P12.2fs_Ir.dyall_v2z.lsym.out`, 
:download:`Z61.x2c.2Paver.scf_cc33e_oe_2P12.2fs_Ir.dyall_v2z.lsym.out <../../../../test/tutorial_Ir_16plus/result/Z61.x2c.2Paver.scf_cc33e_oe_2P12.2fs_Ir.dyall_v2z.lsym.out`.


Then we proceed to the second excited 2P32 state:
::

 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e_2P32.2fs.inp
 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e_oe_2P32.2fs.inp

Files to download are
:download:`Z61.x2c.2Paver.scf_cc33e_2P32.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e_2P32.2fs.inp>`,
:download:`Z61.x2c.2Paver.scf_cc33e_oe_2P32.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e_oe_2P32.2fs.inp>`.
Corresponding output files are
:download:`Z61.x2c.2Paver.scf_cc33e_2P32.2fs_Ir.dyall_v2z.lsym.out <../../../../test/tutorial_Ir_16plus/result/Z61.x2c.2Paver.scf_cc33e_2P32.2fs_Ir.dyall_v2z.lsym.out`, 
:download:`Z61.x2c.2Paver.scf_cc33e_oe_2P32.2fs_Ir.dyall_v2z.lsym.out <../../../../test/tutorial_Ir_16plus/result/Z61.x2c.2Paver.scf_cc33e_oe_2P32.2fs_Ir.dyall_v2z.lsym.out`.


The :math:`^{2}S_{1/2}`,:math:`^{2}P_{1/2}` and :math:`^{2}P_{3/2}` FSCCSD states
---------------------------------------------------------------------------------

Simple and very stable approach to obtain ground and multiple excited states in one step is the Fock-space Coupled Cluster method.
Starting from the closed-shell system, Ir(+17), one gets - by solving (01) sector - 
all three states of interest, :math:`^{2}S_{1/2}`,:math:`^{2}P_{1/2}` and :math:`^{2}P_{3/2}`:
::

 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.scf_fscc01_33ce_5s5p.2fs.inp

Input file to download is
:download:`Z61.x2c.scf_fscc01_33ce_5s5p.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.scf_fscc01_33ce_5s5p.2fs.inp>`.
Corresponding output file is
:download:`Z61.x2c.scf_fscc01_33ce_5s5p.2fs_Ir.dyall_v2z.lsym.out <../../../../test/tutorial_Ir_16plus/result/Z61.x2c.scf_fscc01_33ce_5s5p.2fs_Ir.dyall_v2z.lsym.out`. 



Overview of excitation energies
===============================

In the following table we summarize excitation energies. All values are in a.u. Energies in the Table are not rounded, the are cut at 8
decimal places. (oe= used orbital energies)

============   ===============           ===============          ================          ===========   ===========
Method              2S1/2                     2P1/2                    2P3/2                  S-2P12         S-P32
============   ===============           ===============          =================         ===========   ===========
(SCF ref)
SCF            -17751.10181462           -17749.67796221          -17748.96107014            1.42385       2.140744
CCSD           -17751.90589433           -17750.48885480          -17749.77167089            1.41703953    2.13422344
CCSDoe         -17751.90589433           -17750.48885479          -17749.77167088            1.4170395     2.13422
CCSD(T)        -17751.90998637           -17750.49777405          -17749.78015403            1.41221232    2.12983234
CCSD(T)oe      -17751.90999205           -17750.49779342          -17749.78020119            1.412198      2.129790
(CC ref)
CCSD           -17751.90589433           -17750.48895732          -17749.77154045            1.416937      2.13435
CCSDoe         -17751.90589433           -17750.48895728          -17749.77154043            1.416937      2.134353
CCSD(T)        -17751.90998637           -17750.49779338          -17749.78012458            1.41219       2.12986
CCSD(T)oe      -17751.90999205           -17750.49784867          -17749.78024342            1.41214       2.12974
FSCCSD         -17751.90324662           -17750.48431436          -17749.76770431            1.4189        2.1355
============   ===============           ===============          =================         ===========   ===========




