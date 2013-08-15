Basis set superposition error (BSSE)
====================================

Example for the HF molecule in the small (STO-2G) decontracted basis set.

For the 'ghost'(Gh) atoms the user has to type (copy) corresponding basis sets from the file.

For the F atom in the molecular H(_Gh)F basis you have to give the SCF occupation explicitly as the
program's autooccupation feature gives wrong (nonconverging) open-shell state.
Instead we define the averaged ^2P state (5 electrons distributed over 6 2p shells), which is converging well.

For the H atom in the HF(_Gh) molecular basis autoccupation choice is sufficient as there is only one open-shell electron.

SCF energies: 
-------------
HF      : -95.96439679866987
H(_Gh)F : -0.48328856732273590
HF(_Gh) : -95.419915918676935 (2Paver)
