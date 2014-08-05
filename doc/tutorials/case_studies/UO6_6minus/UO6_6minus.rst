 ==========================
 The :math:`UF_6(6-)` anion
 ==========================

Another (difficult) case study example is the UO6^{6-} anion.
The molecule is in squashed geometry (not fully octahedreal). Axial and equatorial bond distances
are not equal.


Working strategy:
=================

The ground state of :math:`UF_6(6-)` converges fine, but not 
with the default DIIS and starting from screened bare nucleus. Rather one has to use the robust damping scheme.
 Try first DAMPFC 0.9: ::

 pam  --noarch  --get "DFCOEF=DFCOEF.UO6_6m_ground_1pass" --replace dampfactor=0.90 --inp=UO6_6minus.x2c.scf_damp.inp  --mol=UO6.mol

(Files to download are :download:`UO6_6minus.x2c.scf_damp.inp` and :download:`UO6.mol`.)

Then, after we get DFCOEF, decrease the damping level to 0.65. One could try also some  lower value: ::

 pam  --noarch --put "DFCOEF.UO6_6m_ground_1pass=DFCOEF" --get "DFCOEF=DFCOEF.UO6_6m_ground_2pass" --replace dampfactor=0.50 --inp=UO6_6minus.x2c.scf_damp.inp  --mol=UO6.mol

And we get converge molecular spinors saved in the "DFCOEF=DFCOEF.UO6_6m_ground_2pass" file.

The :math:`UF_6(5-)` anion
===========================

Now to the next :math:`UF_6(5-)` state.
With 2pi_3/2 electron out of the U atom there is a difficult convergence. There are five runs for that. ::

 ./pam  --noarch --put "DFCOEF.UO6_6m_ground_2pass=DFCOEF" --get "DFCOEF=DFCOEF.UO6_5m_2p3_1pass" --replace dampfactor=0.90 --inp=UO6_5m.x2c.scf_2p3_damp.inp --mol=UO6.mol

The damping scheme does not converge, it slightly oscilates. Try OPENFAC=0.0.

./pam  --noarch --put "DFCOEF.UO6_5m_2p3_1pass=DFCOEF" --get "DFCOEF=DFCOEF.UO6_5m_2p3_2pass" --replace dampfactor=0.60 --inp=UO6_5m.x2c.scf_2p3_damp.inp --mol=UO6.mol

... does not converge...

Another ways is first run  cca 500 iterations with the damping factor of 0.65 and shift of virtuals of +0.65. This goes very slowly.
The pass two continues from pass 1. Pass 3 has ... Pas 4 has another 500 iteration. In pass 5 the convergence was reached with DIIS.
Damping factors were decreasing from 0.65 to 0.25.
