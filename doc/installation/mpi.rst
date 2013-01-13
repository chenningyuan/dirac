

Parallel run information
------------------------

My shell environment variables are not forwarded to the compute nodes. What do
i have to do?

Assuming you are using our submit script pam to submit your parallel
calculation you can forward any environment variable, e.g., LD\_LIBRARY\_PATH,
to the slave-node shell by adding to you pam command line:

-  ./pam --mpiarg="-x LD\_LIBRARY\_PATH" (for **OpenMPI**)
-  ./pam --mpiarg="-envall" (for **IntelMPI** and **MPICH(2)**)

In the latter case ALL environment variables will be forwarded in one shot
whereas for **OpenMPI** you may need to provide the most essential environment
variables individually. These may be for a typical parallel Dirac run:

-  LD\_LIBRARY\_PATH
-  BASDIR
-  PATH
