Polarizability of the open-shell system - sodium atom 
======================================================

Neither RPA (or TDHF), nor the linear response method are working for 
open-shell systems for the moment.

For the polarizability of the sodium atom (STO-2G uncontracted basis set) the finite field method gives:

SCF,z-finite field:
polariz=-2*(E(F)-E(0))/(F*F)=(2*(-155,86637808348794+155,86636965856189))/(0,0005*0,0005)=+67,3994084
polariz=-(2*(-155,86640335028159+155,86636965856189))/(0,001*0,001)=+67,3834394

Contrary, the RPA method (from RELCCSD) gives different numbers:

x=y=z: 52.71783 au (RPA,orbital energies)
RPA,recomp. fock diagonal elements:
x=60.6522490299745
y=60.8187426643006 
z=54.4170346000506

For this case the linear response module quits upon missing the DV active matrix.
