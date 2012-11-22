

CAM functional
==============

In DIRAC, you can calculate with the CAMB3LYP functional using::

  .DFT
   CAMB3LYP

It is also possible to change the default parameters. The following line
will reproduce the above functional::

  .DFT
   CAM p:alpha=0.19 p:beta=0.46 p:mu=0.33 x:slater=1 x:becke=1 c:lyp=0.81 c:vwn5=0.19

but it will give you complete freedom over the parameters.  As always, please
verify your results carefully when changing these parameters.
