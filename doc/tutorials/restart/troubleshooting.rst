:orphan:


Troubleshooting
===============

Can I move the CHECKPOINT.h5 file from one machine to another and restart from it?
-----------------------------------------------------------------------------------

Yes. Since DIRAC 22 the default checkpoint file is CHECKPOINT.h5, which is in
HDF5 format and is machine-independent. You can move it freely between
different architectures and DIRAC will restart correctly.

If you have an old binary DFCOEF file from a pre-DIRAC22 calculation you
cannot in general use it on a different architecture. In that case you can
convert it to the portable DFPCMO text format using the .PCMOUT keyword
first, and restart from DFPCMO on any machine.


I tried to restart but get the error "Incompatible number of basis functions" - what does it mean?
--------------------------------------------------------------------------------------------------

This probably means that you have changed the basis set. Presently DIRAC cannot
restart from CHECKPOINT.h5 or DFPCMO if the basis set changes (more correctly if
matrix dimensions change).
