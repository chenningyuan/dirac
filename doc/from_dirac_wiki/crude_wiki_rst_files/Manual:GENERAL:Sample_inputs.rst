Manual:GENERAL:Sample inputs
============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Use a reduced value for the speed of light to amplify the effect of
relativity (for testing):

::

    **GENERAL
    .CVALUE
     50.0

Very often one needs to export or to read formatted file with orbital
coefficients - DFPCMO. Therefore the user should provide for the reading

::

    **GENERAL
    .PCMOIN

or

::

    **GENERAL
    .PCMOUT

for the export.

