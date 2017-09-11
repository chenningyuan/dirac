:orphan:

Counting the RelCC memory consumption
=====================================

It is of great importance to know the exact number of
the memory needed for a high memory and time demanding  RelCC job.

SCF
---
The first step is the  SCF run with saving  molecular orbitals, preferably as DFPCMO.
Afterward, prepare the active space for MOLTRA and RelCC.

MOLTRA & RelCC
--------------
Equip the MOLTRA input with the ".NO4INDX" keyword and RelCC with the .COUNTMEMORY keyword.
Run SCF with saved molecular orbitals. The job passes quickly through MOLTRA as it does only
two-index transformations.

In the RelCC module, the program count the memory and leave. 

Example
-------

Demonstration of this feature is in the DIRAC *test/count_cc_memory*, see
:download:`test  <../../../../test/count_cc_memory/test>`

We found that minimal memory setting (via the *--ag* flag) for the memory counting run:

::

 pam --noarch --gb=0.20 --ag=0.366 --inp=H2O.ae4z.x2c.scf_countmem-relcc.inp  --mol=H2O.xyz --put "DFPCMO.H2O.x2c.ae4z=DFPCMO" 

Input files are :download:`H2O.ae4z.x2c.scf_countmem-relcc.inp  <../../../../test/count_cc_memory/H2O.ae4z.x2c.scf_countmem-relcc.inp>`, :download:`H2O.xyz  <../../../../test/count_cc_memory/H2O.xyz>`, producing the output 
:download:`file  <../../../../test/count_cc_memory/result/H2O.ae4z.x2c.scf_countmem-relcc_H2O.out>`.

We are to look the "Peak memory usage" to get reasonable estimate of dynamic memory:

::

 Predicted RelCC memory demand:           0.335 GB
 Peak memory usage    (Gb) :   0.363

Next is the sharp run based on previous relcc memory conuting.
First, the serial run 
needs at least 0.366 GB , rounding to 0.4 GB of the total memory:

::

 pam --noarch --gb=0.20 --ag=0.366 --inp=H2O.ae4z.x2c.scf_relcc.inp  --mol=H2O.xyz --put "DFPCMO.H2O.x2c.ae4z=DFPCMO"

The following parallel run
must have the  memory of at least 16*0.366 = 5.856 GB, rounding to 5.9 GB, 
which must be free at the running machine with at least 16 CPUs:

::

 pam --noarch --mpi=16 --gb=0.20 --ag=0.366 --inp=H2O.ae4z.x2c.scf_relcc.inp  --mol=H2O.xyz --put "DFPCMO.H2O.x2c.ae4z=DFPCMO"  

giving the output 
:download:`file  <../../../../test/count_cc_memory/result/H2O.ae4z.x2c.scf_relcc_H2O.out>`.
We check the peak memory usage which is cca 91.5% of the *--ag* parameter.

::

 Predicted RelCC memory demand:           0.335 GB
                Peak memory usage:        0.335 GB

