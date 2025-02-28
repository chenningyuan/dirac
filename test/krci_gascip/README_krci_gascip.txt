hjaaj 11. Feb. 2025

The test/krci_gascip is a copy of test/krci_energy with luciarel,
only difference is that GASCIP is used instead of LUCIAREL.

The F atom gives slightly higher energies compared to krci_energy:

  LUCIAREL        Final CI energies  =         -99.49492593        -99.49492593        -99.49299297
  (GASCIP_RCISTD) Final CI energies  =         -99.49443766        -99.49423784        -99.49233285

but no obvious reason. Maybe some symmetry breaking in LUCIAREL? The number of determinants is also
greater in LUCIAREL.
(That is also the case for Be test, but there the final energies agree.)

The two H2O tests agree both in number of determinants and in final energy.

