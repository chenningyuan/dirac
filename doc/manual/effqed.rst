:orphan:
 

star(EFFQED)

Options for the numerical integration of effective QED potentials described in :cite:`Sunaga2022`.
An example of the input is shown in `tutorial <..//tutorials/qed/effqed.html>`__.

keyword(GRASP)
Activate the radial grid of GRASP code.

*Default:*

::

     Lindh's radial grid

keyword(ANGQED)
Angular momentum of Lebedev quadrature for the numerical intebration of effective QED potential.

*Default:*

::

    .ANGQED
     15

keyword(UEHCUT)
Due to the very local nature of the effective QED potentials, the upper limit of
the radial integration over a potential associated with atomic center A is cut off.
The radial integration is cut off based on the Uehling potential and is performed in a shorter region. 

*Default:*

::

    Cutoff based on Flambaum-Ginges low-frequency term


keyword(MINIMUM)
Decision of the smallest nuclear charge Z to include effective QED potentials. 

*Default:*

::

    .MINIMUM
     19

The low-frequency (LF) term of the Flambaum-Ginges potential has the range of the 1s-orbital of a hydrogen-like atom, so the LF term of light atoms are rather widely distributed and may overlap with the neighboring atoms in the molecule. This is the reason why light elements are skipped in the default.

keyword(ALLATOM)
Incorporating effective QED potentials for all atoms in the molecule.
*Default:*

::

    The element with Z < 19 is not calculated.

keyword(SCREEN)
The AO basis set less than the threshold in a grid point is omitted in the numerical integration,
which is done in the region close the nucleus.
For example, in AB molecule, tight exponent of atom B would not be considered
in the numerical integration of atom A's effQED potentials. 
*Default:*

::

    .SCREEN                                                  
     1.0E-10
