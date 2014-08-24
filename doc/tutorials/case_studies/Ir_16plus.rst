==============================
Getting excited states Ir(16+)
==============================

In the following we are interested in getting two closest excited states (^2P1/2,^2P3/2) of the Ir(16+) cation and their CCSD(T) correlation energies.
::

 Ir(16+) : [Xe] 4f(14) 5s(1) 5p1/2(0)  5p3/2(0)  (^2S)
 Ir(16+) : [Xe] 4f(14) 5s(0) 5p1/2(1)  5p3/2(0)  (^2P1/2)
 Ir(16+) : [Xe] 4f(14) 5s(0) 5p1/2(0)  5p3/2(1)  (^2P3/2)


For simplicity, we will work with the .X2C Hamiltonian and with the smallest Dyall.v2z decontracted basis set.

There are two ways to achieve, (i) from the converged SCF state, and, (ii) only at the correlated level from the 2P_aver SCF state.

The ^2S1/2 ground state
-----------------------

Thanks to the symmetry with two irreps, one can place the unpaired electron either the 1st irrep (2^S1/2) or into 2nd irrep (2P1/2).
::

 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2S12.scf_cc33e.2fs.inp
 pam --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2S12.scf_cc33e_oe.2fs.inp

Files to download  
   :download:`Ir.dyall_v2z.lsym.mol <../../../../test/tutorial_Ir_16plus/Ir.dyall_v2z.lsym.mol>`
                   :download:`Z61.x2c.2S12.scf_cc33e.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2S12.scf_cc33e.2fs.inp>`
                   :download:`Z61.x2c.2S12.scf_cc33e_oe.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2S12.scf_cc33e_oe.2fs.inp>`

^2P1/2,^2P3/2 SCF states
------------------------

Thanks to the symmetry with two irreps, one can place the unpaired electron either the 1st irrep (2^S1/2) or into 2nd irrep (2P1/2).

How to obtain the ^2P3/2 excited state ? For that, we use the REORDER MO keyword and read MO from the previous step.
::

 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2P12.scf_cc33e.2fs.inp --get "DFCOEF=DFCOEF.v2z.2P12"
 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2P12.scf_cc33e_oe.2fs.inp 

Files to downlaod
   :download:`Z61.x2c.2P12.scf_cc33e.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2P12.scf_cc33e.2fs.inp>`
   :download:`Z61.x2c.2P12.scf_cc33e_oe.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2P12.scf_cc33e_oe.2fs.inp>`


For the 2P3/2 excited state one has use overlap selection
::

 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2P32.scf_cc33e.2fs.inp --put "DFCOEF.v2z.2P12=DFCOEF"
 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2P32.scf_cc33e_oe.2fs.inp --put "DFCOEF.v2z.2P12=DFCOEF"

Files to download ...
   :download:`Z61.x2c.2P32.scf_cc33e.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2P32.scf_cc33e.2fs.inp>`
   :download:`Z61.x2c.2P32.scf_cc33e_oe.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2P32.scf_cc33e_oe.2fs.inp>`


^2P1/2,^2P3/2 CCSD(T) states
----------------------------

The other option is to start from the 2P average single determinant state in the linear symmetry, and distinguish 
between ^2P1/2,^2P3/2 states at the Coupled Cluster correlated scheme. 

First we test the pure 2Paver state:
::

 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e.2fs.inp
 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e_oe.2fs.inp

Files to download ...
   :download:`Z61.x2c.2Paver.scf_cc33e.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e.2fs.inp>`
   :download:`Z61.x2c.2Paver.scf_cc33e_oe.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e_oe.2fs.inp>`

And afterwards we proceed to the individual states:
::

 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e_2P12.2fs.inp
 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e_oe_2P12.2fs.inp

 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e_2P32.2fs.inp
 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.2Paver.scf_cc33e_oe_2P32.2fs.inp

Files to download ...
   :download:`Z61.x2c.2Paver.scf_cc33e_2P12.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e_2P12.2fs.inp>`
   :download:`Z61.x2c.2Paver.scf_cc33e_oe_2P12.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e_oe_2P12.2fs.inp>`

   :download:`Z61.x2c.2Paver.scf_cc33e_2P32.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e_2P32.2fs.inp>`
   :download:`Z61.x2c.2Paver.scf_cc33e_oe_2P32.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.2Paver.scf_cc33e_oe_2P32.2fs.inp>`


^2P1/2,^2P3/2 FSCCSD states
---------------------------

Simple approach to obtain multiple excited states is the Fock-space Coupled Cluster method.
Starting from the closed-shell Ir(+17) one gets in (01) sector all three states, 2S, 2P12, 2P32.
::

 $PAM --noarch --mw=120  --mol=Ir.dyall_v2z.lsym.mol --inp=Z61.x2c.scf_fscc01_33ce_5s5p.2fs.inp

File to download ...
   :download:`Z61.x2c.scf_fscc01_33ce_5s5p.2fs.inp <../../../../test/tutorial_Ir_16plus/Z61.x2c.scf_fscc01_33ce_5s5p.2fs.inp>

Overview of excitation energies
===============================

In the following table we summarize excitation energies. All values are in a.u.


                    2S1/2                     2P1/2                    2P3/2                   S-2P12       S-P32
                --------------------------------------------------------------------------------------------------
SCF            -17751.101814622431       -17749.677962219408      -17748.961070144513         1.42385    2.140744
CCSD           -17751.905894332194293    -17750.488854800489207   xx
CCSDoe         -17751.905894338924554    -17750.488854797091335   -17749.771670889032976     1.4170395   2.13422
CCSD(T)        -17751.909986376944289    -17750.497774055980699  xx
CCSD(T)oe      -17751.909992051667359    -17750.497793424234260   -17749.780201190642401     1.412198    2.129790

(CC reference)
CCSD          -17751.905894332194293     -17750.488957325971569   -17749.771540450528846      1.416937   2.13435
CCSDoe        -17751.905894338924554     -17750.488957285222568   -17749.771540437406657      1.416937   2.134353
CCSD(T)       -17751.909986376944289     -17750.497793389531580   -17749.780124582986900      1.41219    2.12986
CCSD(T)oe     -17751.909992051667359     -17750.497848677994625   -17749.780243420427723      1.41214    2.12974


FSCCSD        -17751.903246623278        -17750.484314362893       -17749.767704318059         1.4189    2.1355

