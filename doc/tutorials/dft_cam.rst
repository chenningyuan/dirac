CAM functional
==============

In DIRAC, you can calculate with the CAMB3LYP functional using

::

    .DFT
     CAMB3LYP

It is also possible to change the default parameters. The following line
will reproduce the above functional:

::

    .DFT
     CAM p:alpha=0.19 p:beta=0.46 p:mu=0.33 x:slater=1 x:becke=1 c:lyp=0.81 c:vwn5=0.19

but it will give you complete freedom over the parameters. Note that for
this to work correctly you have to apply the following
`patch </dirac/index.php/Known_bugs#User-defined_CAM_functional_line_is_cut_after_character_80>`_.
But as always, please verify your results carefully when changing these
parameters.

