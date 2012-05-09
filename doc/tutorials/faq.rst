FAQ
===

Some frequently asked questions.

Building
--------

I am about to install DIRAC on a 64-bit machine, should I compile using 64-bit integers or 32-bit integers?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DIRAC runs smoothly (with few exceptions; see below) on 64-bit platforms
using either 32-bit or 64-bit integers.

It is easier to install for 32-bit integers because MPI and math
libraries are typically built for 32-bit integers. Therefore DIRAC is by
default built for 32-bit integers.

If you decide to build for 32-bit integers, you will **not** be able to:

-  Safely allocate more than 16 GB of memory **per cpu-core**
-  Run LUCITA
-  Run really, really large Coupled-Cluster calculations

If you want to be able to run the above calculations, you have to build
for 64-bit integers.

If you decide to build for 64-bit integers, you **have to** make sure
that:

-  The math library you link to can handle 64-bit integers
-  MPI library is built for 64-bit integers

Otherwise DIRAC will stop.

Running Jobs with DIRAC
-----------------------

Parallel run information
^^^^^^^^^^^^^^^^^^^^^^^^

My shell environment variables are not forwarded to the compute nodes. What do i have to do?
                                                                                            

Assuming you are using our submit script pam to submit your parallel
calculation you can forward any environment variable, e.g.,
LD\_LIBRARY\_PATH, to the slave-node shell by adding to you pam command
line:

-  ./pam --mpiarg="-x LD\_LIBRARY\_PATH" (for **OpenMPI**)
-  ./pam --mpiarg="-envall" (for **IntelMPI** and **MPICH(2)**)

In the latter case ALL environment variables will be forwarded in one
shot whereas for **OpenMPI** you may need to provide the most essential
environment variables individually. These may be for a typical parallel
Dirac run:

-  LD\_LIBRARY\_PATH
-  BASDIR
-  PATH

General run information
^^^^^^^^^^^^^^^^^^^^^^^

How can I restart from the MO coefficients file DFCOEF?
                                                       

Simply copy DFCOEF to the scratch directory using the pam script. If
DIRAC finds DFCOEF in the scratch directory it will automatically try to
restart from it.

Can I move the DFCOEF file from one machine to another and restart from it?
                                                                           

In general no, unless it is the same architecture.

What you can do instead is to export the DFCOEF in a machine independent
format using the .PCMOUT keyword. This creates a file called DFPCMO.
Copy DFPCMO to the other machine, make it avaliable in the scratch
directory using the pam script, then restart from it using .PCMOIN.

**Note:** if you are running a calculation in linear symmetry you will
find an indication in the DIRAC output (.PCMOUT run) of whether you need
the keyword .PCMOIN or .PCMOI3 for the restart. This is connected to the
formatted writing of DFCOEF and the space needed to store the mj-value
of each spinor.

I tried to restart but get the error "Incompatible number of basis functions" - what does it mean?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This probably means that you have changed the basis set. Presently DIRAC
cannot restart from DFCOEF if the basis set changes (more correctly if
matrix dimensions change).

I would like to use the 2c-molecular-mean-field Hamiltonian (X2Cmmf) for my post-HF (Coupled-Cluster [CC], Fock-Space CC, KR-MCSCF or KR-CI) calculations - what do i need to consider?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The 2c-molecular-mean-field Hamiltonian (X2Cmmf) is at present (Dirac11
release Nov 2011) **available** for the following methods:

-  Coupled-Cluster
-  Fock-space Coupled Cluster (including the intermediate Hamiltonian
   version)

The remaining modules are **not supported** at the Dirac11 release date
but patches are planned in the course of 2012:

-  KR-MCSCF
-  KR-CI
-  TD-DFT

The X2Cmmf Hamiltonian scheme can be used for the Dirac-Coulomb as well
as the Dirac-Coulomb-Gaunt Hamiltonian (add .GAUNT to \*\*HAMILTONIAN).

As follows we will consider the necessary inputs for an X2Cmmf
calculation where we perform the SCF + post-SCF calculation in **one**
shot or **two separate** runs:

-  **one-step** procedure:

   -  add .X2Cmmf to the \*\*HAMILTONIAN keyword section:

::

    **HAMILTONIAN
    .X2Cmmf

that's it. After the SCF calculation the program will proceed with the
transformation to 2-components prior to the post-HF steps.

-  **two-step** procedure (**EXPERT mode**!):

1. ***SCF step***. Add .X2Cmmf to the \*\*HAMILTONIAN keyword section:

::

    **HAMILTONIAN
    .X2Cmmf

and save the file DFCOEF (keyword *--outcmo* with our default submit
script *pam*)

2. ***post-HF step***. Remove .X2Cmmf from the \*\*HAMILTONIAN section
(**non-intuitive (!)** step which will be solved in Dirac12) and add the
following keywords to the \*\*MOLTRA and namelist RELCCSSD input
section.

::

    **MOLTRA
    .INTFL2
    1 0 0 0 
    .INTFL4
    1 0 0 0   
    ...
    ...
    &CCSORT USEOE=T, NORECMP=T &END

You are now ready to restart your post-HF calculation but do not forget
to copy your MO-coefficient file DFCOEF to your scratch directory, i.e.
add *--incmo* to your *pam* command line.

