MgF molecule (2Sigma)
=====================

P-odd interaction constants
---------------------------

Demonstration of DC-SCF and DC-DFT calculations of the :math:`W_{A}` parameter of the P-odd spin-rotational Hamiltonian
for small diatomic molecule.

The operator is T-odd, and has the cector form of  "alpha_x*rho, alpha_y*rho, alpha_z*rho".

For more details see Ref. [Borschevsky2012]_ .

Checking the sole electron density at the nucleus
=================================================

.RHONUC expectation value is confirmed by individual components of the Fermi contact operator,
 multiplied with the factor (3/4*pi*g_e)= 0.11922794437D0
( DIAGONAL 'FC F  01' )

These results match with PVC operator MO-matrix elements.

Values for the valence (11th) spinor contributions (in a.u.)

GNuc - Gaussian (finite) nucleus
PNuc - point like nucleus

            F                   Mg
GNuc:  0.124281062802    0.110710184980 (Fermi-contact multiplied with the factor)
       0.1242810255      0.1107101196   (PVC integrals-prptra elements, factor of 1.0D0)

PVC integrals (factor 1.0D0),        varrying nuclear exponents (identical for both atoms)
       0.1242811947      0.11071005            1.0D+20
       0.1242811947      0.1107100466          1.0D+50
       0.1242811947      0.1107100466          1.0D+90
        NaN                NaN                 1.0D+120

FC integrals, exp. value,                    GNuc exp.value
       0.124281194663    0.110710046559        1.0D+90
       0.124281194663    0.110710046559       1.0D+120

PNuc:  0.124281194663    0.110710046559 (Fermi-contact multiplied with the factor)
        NaN                  NaN        (PVC integrals-prptra elements, factor of 1.0D0)

Checking alpha*rho operator matrix elements
===========================================

FC operator matrix elements are mutiplied with the above given factor.
For PVC operator we keep factors of 1.0D0.

MgF molecule
------------

Table: matrix elements +/- omega states of the valence open-shell spinor, GNuc model

                 F                   Mg
     PVC:  -0.0020102470       -0.0044718635
      FC:  -0.0020102465       -0.0044718611

RnF(65+) MgF-like molecule (21 electrons)
-----------------------------------------

We simulated diatomic open-shell system with heavy nucleus (Z=86 - Rn atom), RnF,
but with 21 electrons only, having thus high charge of 65+.

MO-transformed matrix elements of alpha*rho operators:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Table: matrix elements between +/- omega states of the valence open-shell spinor, GNuc model

                F                  Rn
      PVC   0.0000799991       -0.0000165875
       FC   0.0000799991       -0.0000166104   

Density at the nucleus (scalar) operators :
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GNuc model, valence spinor (#11) exp.val. (in a.u.):
                F                  Rn
      PVC   0.29382975D-03   0.26175495D-03
       FC   0.29383985D-03   0.26190598D-03

Total exp.val, GNuc model (in a.u.):
                F                  Rn
      PVC   190.2191404628     1564150.328748
       FC   190.2190907103     1562379.047052

