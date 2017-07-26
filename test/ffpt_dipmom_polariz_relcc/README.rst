BeH electric properties
=======================

This toy molecule in the very small STO-2G decontracted basis
demonstrates the use of finite-field perturbation theory for 
calculating the dipole moment and the polarizability at the Hartree-Fock
and Coupled Cluster levels.

Because BeH is open-shell system, we can not use MP2 and CCSD gradients, which are working 
only for closed-shell systems. We have to resort to the finite-field perturbation scheme, 
where to the X2C relativistic Hamiltonian we add z-dipole moment operator as small perturbation.

Let us first compute the dependence of the total energies on the field strength. 
This can be done through replacing the input string by filed value, like

::
 pam --noarch  --replace zff=+0.0005 --mol=BeH.sto-2g.C2v.mol --inp=BeH.x2c_scf_relcc_ffz.inp --put "DFPCMO.BeH.x2c_scf_sto-2g.C2v=DFPCMO"

(Input files for download are 
:download:`BeH.x2c_scf_relcc_ffz.inp  <../../../test/ffpt_dipmom_polariz_relcc/BeH.x2c_scf_relcc_ffz.inp>` ,
:download:`BeH.sto-2g.C2v.mol <../../../test/ffpt_dipmom_polariz_relcc/BeH.sto-2g.C2v.mol>`)

Energy dependent on field strength are sorted in the following Table:

=======    ===================    ===================
zfield          E(SCF)                 E(CCSD(T))
=======    ===================    ===================
+0.0010    -14.424976681856371    -14.461324673372751
+0.0005    -14.425860788235255    -14.462193434302819
+0.0000    -14.426746572346620    -14.463063980838802
-0.0005    -14.427634036448282    -14.463936314225093
-0.0010    -14.428523182834038    -14.464810435757153
=======    ===================    ===================

(NOTE: In this test, only two output files are given for Coupled Cluster method,
:download:`BeH.x2c_scf_relcc_-0.0005_BeH.sto-2g.C2v.out  <../../../test/ffpt_dipmom_polariz_relcc/result/BeH.x2c_scf_relcc_-0.0005_BeH.sto-2g.C2v.out>` ,
:download:`BeH.x2c_scf_relcc_+0.0005_BeH.sto-2g.C2v.out  <../../../test/ffpt_dipmom_polariz_relcc/result/BeH.x2c_scf_relcc_+0.0005_BeH.sto-2g.C2v.out>`.)


Spreadsheet
~~~~~~~~~~~
In the spreadsheet, we interpolate the 4th-degree polynomial over this data - see the 
:download:`data.odf file <../../../test/ffpt_dipmom_polariz_relcc/data.ods>` file. For that, we employ function LINEST,
where we also prepare columns with 1st,2nd and 3rd finite-fields.

Dipole moment
-------------
The SCF expectation value dipole moment - electronic contribution - is -1.77324745 au
(see input file, :download:`BeH.x2c_scf_dipmom.inp  <../../../test/ffpt_dipmom_polariz_relcc/BeH.x2c_scf_dipmom.inp>` and the
the corresponding output file, :download:`BeH.x2c_scf_dipmom_BeH.sto-2g.C2v.out  <../../../test/ffpt_dipmom_polariz_relcc/result/BeH.x2c_scf_dipmom_BeH.sto-2g.C2v.out>`.)

From the LINEST function table, we obtain the value of 1.7732474 a.u., which is the first order expansion coefficient.

Simple first numerical derivation of the SCF perturbed energy according to one applied electric field, dm(SCF;F)=E(+F)-E(-F)/2*F,
gives 1.771568 au for F=0.0005.  For the field strength of F=0.001 it makes 1.76989 a.u.

Clearly, the zero field (F=0) derivative of the fourth order polynomial is more precise
when comparing against the SCF expectation value.

Polarizability
--------------
Simple second derivate calculation for given F=0.0005 at SCF level, 
polariz(SCF;F)=-(E(+F)+E(-F)-2*E(0))/(F*F), 
gives 6.71996 a.u., while
two times the second order expansion coefficient gives (-2)*(-3.3599745995852954)=6.71994919 a.u.

Similar calculation for CCSD(T), F=0.0005,
polariz(CCSD(T);F)=-(E(+F)+E(-F)-2*E(0))/(F*F),
gives
7.147401232 a.u.. 
The corresponding expansion coefficient (multiplied by factor two) makes 7.14738420 a.u. 
One can conclude that
the F=0.0005 field is sufficient enough to obtain polarizability accurate to 3 decimal places by simple numerical second derivation.
