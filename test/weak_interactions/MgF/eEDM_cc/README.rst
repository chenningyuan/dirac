MgF molecule (2Sigma)
=====================

Testing finite field calculations with the (T,P)-odd eEDM operator.

The symmetry must be the lowest, C_{1}, due to symmetry properties of this operator.

The Kramers symmetry breaking eEDM operator is beeing added to the unperturbed Dirac Hamiltonian at the correlated KU-CCSD(T)/FSCCSD level.

CCSD(T) energies:
------------------
E(+0.0000)=-287.049504317674746
E(+0.0001)=-287.049504319253003
E(+0.0005)=-287.049504357131298 
E(+0.0010)=-287.049504475501067

1st derivations: 

(-287.049504319253003+287.049504317674746)/0.0001=-.00001578257
(-287.049504357131298+287.049504317674746)/0.0005=-.000078913104 
(-287.049504475501067+287.049504317674746)/0.0010=-.00015782632

extrapolation for the zero field strength....cca 10^{-11} - property is not stable numerically 

FSCC energies:
--------------
all 21 electrons
 
    +0.0000               +0.0005            +0.0010
-287.036567976389   -287.036568015927    -287.036568134540
-286.901251728430   -286.901251768261    -286.901251887753
-286.901130397672   -286.901130437453    -286.901130556795

1st derivations:
                        -.0000790760         -.0001581510
                        -.0000796620         -.0001593230 
                        -.0000795620         -.0001591230

17 active electrons (1s^2 of F and 1s^2 of Mg frozen)

    +0.0000                +0.0005                +0.0010
-287.011517211004      -287.011517235274     -287.011517308078
-286.876207905122      -286.876207929658     -286.876208003240
-286.876086634744      -286.876086659182     -286.876086732474

1st derivations:
                        -.000048540           -.00009707400
                        -.000049072           -.00009811800
                        -.000048876           -.00009773000

