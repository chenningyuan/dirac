 =================
 The UO6_6- anion
 ================

Another case example is the UO6^{6-} anion.

The molecule is in squashed geometry (not fully octahedreal). Axial and equatorial bond distances
are not equal.

The ground state converges fine.  We 2pi_3/2 electron out of U atom and this has difficult convergence. There are five runs for that.

The first run is cca 500 iterations with the damping factor of 0.65 and shift of virtuals of +0.65. This goes very slowly.

The pass two continues from pass 1. Pass 3 has ... Pas 4 has another 500 iteration. In pass 5 the convergence was reached with DIIS.

Damping factors were decreasing from 0.65 to 0.25.


Working strategy:

The UO6_6m ground state does not converge with the default DIIS and starting from screened bare nucleus. Try first DAMPFC 0.9

/pam  --noarch  --get "DFCOEF=DFCOEF.UO6_6m_ground_1pass" --replace dampfactor=0.90 --inp=UO6_6minus.x2c.scf_damp.inp  --mol=UO6.mol

Then, after we get DFCOEF, decrease the damping level to 0.65. One could try also some  lower value

/pam  --noarch --put "DFCOEF.UO6_6m_ground_1pass=DFCOEF" --get "DFCOEF=DFCOEF.UO6_6m_ground_2pass" --replace dampfactor=0.50 --inp=UO6_6minus.x2c.scf_damp.inp  --mol=UO6.mol

And we get converge molecular spinors written in the "DFCOEF=DFCOEF.UO6_6m_ground_2pass" file.


No to the UO6_5m state:

 ./pam  --noarch --put "DFCOEF.UO6_6m_ground_2pass=DFCOEF" --get "DFCOEF=DFCOEF.UO6_5m_2p3_1pass" --replace dampfactor=0.90 --inp=UO6_5m.x2c.scf_2p3_damp.inp --mol=UO6.mol


./pam  --noarch --put "DFCOEF.UO6_5m_2p3_1pass=DFCOEF" --get "DFCOEF=DFCOEF.UO6_5m_2p3_2pass" --replace dampfactor=0.60 --inp=UO6_5m.x2c.scf_2p3_damp.inp --mol=UO6.mol
