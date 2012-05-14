\*\*GENERAL -- Sample inputs
============================

Use a reduced value for the speed of light to amplify the effect of
relativity (for testing):

::

    **GENERAL
    .CVALUE
     50.0

Contrary you may want to extrapolate to the nonrelativistic limit by setting
speed of light to very large value. You have to experiment a little to find out appropriate number
as this depends on the character of the studied system.

Very often one needs to export or to read formatted file with orbital
coefficients - DFPCMO. Therefore the user should provide

::
    **GENERAL
    .PCMOUT

for the export.


