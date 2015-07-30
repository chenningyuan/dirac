==================
Input reader tests
==================

Collection of various tests for the DIRAC new inpur reader.

The .NEL_F1 occupation of \**RELCC
-----------------------------------
Motivation for this test is that the reading of the .NEL_F1 occupation
of the AsCl(+) molecule in the C(oo,v) is going wrong.

The proper order given in the input is:
::
 NEL_F1: 13 13  5  6  1  1  0  0  0  0  0  0  0  0  0  0

while the relcc_input_reader reads it uncorrectly as:
::
 NEL_F1: 13  3 13  3  5  6  1  1  0  0  0  0  0  0  0  0
