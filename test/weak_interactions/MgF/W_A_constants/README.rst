MgF molecule (2Sigma)
=====================

P-odd interaction constants
---------------------------

Demonstration of DC-SCF and DC-DFT calculations of the :math:`W_{A}` parameter of the P-odd spin-rotational Hamiltonian
for small diatomic molecule.

The operator is T-odd, and has the form of dot-product,  "alpha_x*rho + alpha_y*rho + alpha_z*rho".

For more details see Ref. [Borschevsky2012]_ .

Checking the sole electron density at the nucleus
=================================================

.RHONUC expectation value confirmed by individual components of the Fermi contact operator, multiplied with the factor (3/4*pi*g_e)= 0.11922794437D0
( DIAGONAL 'FC F  01' )

Values for the valence (11th) spinor contributions (in a.u.)

GNuc - Gaussian (finite) nucleus
PNuc - point like nucleus

            F                   Mg
GNuc:  0.124281062802    0.110710184980 (Fermi-contact multiplied with the factor)
       0.1242810255      0.1107101196   (PVC integrals-prptra elements, factor of 1.0D0)

PVC integrals (factor 1.0D0),        varrying nuclear exponents (the same for both atoms)
       0.1242811947      0.11071005            1.0D+20
       0.1242811947      0.1107100466          1.0D+50
       0.1242811947      0.1107100466          1.0D+90
        NaN                NaN                 1.0D+120

FC integrals, exp. value,                    GNuc exp.value
       0.124281194663    0.110710046559        1.0D+90
       0.124281194663    0.110710046559       1.0D+120

PNuc:  0.124281194663    0.110710046559 (Fermi-contact multiplied with the factor)
        NaN                  NaN        (PVC integrals-prptra elements, factor of 1.0D0)


