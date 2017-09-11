W_A operators in DIRAC
======================

The operators are added as finite-field perturbation at the KU-RelCC level.

This operator preferably contains PVC representation of the nuclear density.

Test system: MgF, STO-2G uncontracted basis, all (13) correlated electrons


field  symm        SCF                        CCSD(T)
------------------------------------------------------------
0.0    C2v    -287.76439849222800      -288.008945929251695
0.0    C2     -287.76439849222265      -288.008945929251695
              -287.76439849222265      -288.008945929251695
0.0    Cs     -287.76439849222726      -288.008928821452059
0.0    C1     -287.76439849220532      -288.008944099824760 <--- !


1e-4   C2v    -287.76439849222800      -288.007312476243214 (zalpha)  <-----
       Cs                              -288.008928822665382 (zalpha) ??
       C1                              -288.008945953081763 (zalpha) ??
       C2                              -288.008945930465870 (not active)

1e-4   Cs     -287.76439849222879      -288.008928821451661 (xalpha)  <-----
       C2                              -288.008945929252320 (not active)
       C1                              not converged

1e-4   Cs     -287.76439849222879      -288.008928278187511 (yalpha)  <-----
       C2                              -288.008945929252320 (not active)
       C1                               not converged

