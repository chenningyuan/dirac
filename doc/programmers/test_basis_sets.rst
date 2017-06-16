:orphan:
 
Testing basis set files in DIRAC
================================

The DIRAC program suite hosts more than different 300 basis set input files. These are either relativistic basis sets (all-electron, ECP)
or nonrelativistic basis sets from other sources (Dalton, EMSL etc). Therefore it is necessary to have (automated) tools
for testing basis set files.

Currently, there are two tests focused at reading basis set files:

.. include :: ../../test/basis_input_scripted/README.rst


and test/basis_input/README.rst.
