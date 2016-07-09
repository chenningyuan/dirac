DFT-COSCI tests
===============

Test suite of the new feature DFT-COSCI.

Single determinant open-shell DFT with fractional occupation is followed by
Kramers unrestricted COSCI with exchange-correlation and spin-polarized potential.


Development notes
-----------------

subroutine xcint_potential_rks - whey is fmat intent(in) ? is it modified inside and returned back


dirscf.F:
call xcint_potential_rks(ntbas(0),dmat,fock)


-after DFT cycles - does the converged Fock contain the Vxc contribution, and the Fock is written into DFFCK1 ?

Fock F2 mtx in ERGCAL contains the Vxc contribution, which is not in MOLTRA-for-RESOLVE !

traone.F/GTCFC1 - injection of xcint to core two-el. Fock matrix, gives warning concerning number of electrons

Problem: xcint is added to the core transformed matrix ! we need active (open-shell) part !

src/xcint: does not depend on gp ! gp depends on xcint module ! keep this dependencies in local CMakeFiles.txt

moltra/dcbtra_h.F90 - encapsulating module for moltra/dcbtra.h
moltra/dcbtra.h - refurbished to be included in Fortran90 source files

dirone/dirscf.F ERGCAL - printouts of xcint contributions from KS potential and for the total energy



Preliminary Results
-------------------

Experim. FSS(F)=404.072 +/- 0.025 cm-1 ( see http://scitation.aip.org/content/aip/journal/jcp/72/12/10.1063/1.439122 )

     COSCI,x2c, sto-2g   : 230.95 cm-1
DFT-COSCI,x2c,bp86,sto-2g: 267.92 cm-1 (core dens. only)
DFT-COSCI,x2c,bp86,sto-2g: 283.97 cm-1 (full dens for Vxc)
DFT-COSCI,x2c,bp86,sto-2g: 212.298cm-1 (no xc for Fij elem)

