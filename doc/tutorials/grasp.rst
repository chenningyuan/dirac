Grasp:Orbpri
============

I assume that you have updated Grasp and performed the following steps
in your Grasp-directory

::

    admin make
    ./configure
    make
    make utils

so that you have both grasp.x and orbpri.x in your installation
directory. Now let's get started:

At the 4-component level atomic orbitals are given as

::
 ADDME

Consider the following Grasp-input (on file Cl.ginp):

::

    Chlorine atom
      1 5 1                                                                
    1S                                                                      
    2S                                                                      
    2P                                                                      
    3S                                                                      
    3P 5
    ANG                                                                   
    -1
    MCP 1 2                                                                 
    MCDF 23
    2 0 3 0 4                                                                  
    17 35.453
    NUCLEUS GAUSSIAN                                                        
    ORBOUT 1 - 5
    AL                                                                       
    END   

It is pretty standard except that we ask the radial functions *P*(*r*)
and *Q*(*r*) to be dumped to unformatted file on unit 4. We run this
example:

::

    grasp.x < Cl.ginp > Cl.grasp   

Next, we run orbpri.x to get selected radial functions dumped to
formatted file ready for processing by a 2D visualization program.
orbpri.x is pretty self-explanatory. You may follow this exchange with
the program

::

    saue-laptop 178>orbpri.x
    *** OUTPUT from ORBPRI ***
    Give unit number of orbital file to read
    4
    *** OUTPUT from ORBPRI ***
    * Atomic number:               17.0
    * Number of grid points:       1995
    * First grid point:       0.588E-07
    * Grid step size :        0.100E-01
      
    * Number of orbitals read:        7
      
    Do you want info on orbitals (y/n) ?
    y
     1  1S       -105.291114936211
     2  2S        -10.675917467977
     3  2P-        -8.117105154258
     4  2P         -8.052414693447
     5  3S         -1.079593288072
     6  3P-        -0.509246554903
     7  3P         -0.504056491673
    Do you want some orbitals printed (y/n) ?
    y
    Name atom (A3)
    Cl
    RMAX = :   26.877
    Give index of orbital to print(-1 to quit):
    3
       PMIN :0.0000E+00   PMAX :0.1520E+01
       QMIN :-.3800E-02   QMAX :0.5684E-01
    Give index of orbital to print(-1 to quit):
    -1
    *** EXITING ORBPRI ***

There is now a formatted file *Cl\_2P-* containing three columns, first
the grid points in atomic units, followed by the numerical values of
*P*(*r*) and *Q*(*r*).

