
Non-relativistic results from Dalton:
=====================================


         Isotropic[au]  Diamagnetic[au]  Paramagnetic[au]  Parallel[au]  Perpendicular[au]  Anisotropy[au]

LDA:      -2.1051         -6.6772           4.5721           -2.0301        -2.1426            0.1125         (sc, go=(0,0,0))
          -2.1051         -6.6772           4.5721           -2.0301        -2.1426            0.1125         (nc, go=(0,0,0))

PBE:      -2.0884         -6.6564           4.5680           -2.0178        -2.1237            0.1060         (sc, go=(0,0,0))
          -2.0884         -6.6564           4.5680           -2.0178        -2.1237            0.1060         (nc, go=(0,0,0))



sc = symmetric connection
nc = natural connection
go = gauge origin




Dalton calculations - details
=============================

1. LDA, sc, go=(0,0,0) 

* input:

```
  **DALTON INPUT
  .RUN PROPERTIES
  **WAVE FUNCTIONS
  .DFT
  LDA
  *SCF INPUT
  .DOUBLY OCCUPIED
  5
  .THRESH
  1.0D-7
  .MAX MACRO ITERATIONS
  30
  **PROPERTIES
  .MAGNET
  .NODIFC
  .GAUGEO
   0.0 0.0 0.0
  .PRINT
  10
  *LINRES
  .PRINT
  10
  .THRESH
  1.0D-9
  .MAX IT
  100
  *EXPECT
  .PRINT
  10
  *END OF INPUT
```

* output:

```
  Magnetizabilities                    au          SI         cgs
  ---------------------------------------------------------------

  Isotropic magnetizability:      -2.1051   -166.1124    -10.0035
  Diamagnetic contribution:       -6.6772   -526.8975    -31.7305
  Paramagnetic contribution:       4.5721    360.7851     21.7270
  Parallel component:             -2.0301   -160.1926     -9.6470
  Perpendicular component:        -2.1426   -169.0722    -10.1818
  Anisotropy:                      0.1125      8.8796      0.5347

  Magnetizability is cylindrical.


  Principal values (au, SI, and cgs) and axes:
  --------------------------------------------

  khi_aa      -2.0301   -160.1926     -9.6470      0.0000    0.0000    1.0000
  khi_bb      -2.1426   -169.0722    -10.1818      1.0000    0.0000    0.0000
  khi_cc      -2.1426   -169.0722    -10.1818      0.0000    1.0000    0.0000


                        Total magnetizability tensor (au)
                        ---------------------------------

                       Bx                  By                  Bz

  Bx        -2.142586054864     -0.000000000000      0.000000000000
  By        -0.000000000000     -2.142586054864      0.000000000000
  Bz         0.000000000000      0.000000000000     -2.030058045886



                     Diamagnetic magnetizability tensor (au)
                     ---------------------------------------

                       Bx                  By                  Bz

  Bx        -9.000717658320      0.000000000000     -0.000000000000
  By         0.000000000000     -9.000717658320     -0.000000000000
  Bz        -0.000000000000     -0.000000000000     -2.030058045886



                     Paramagnetic magnetizability tensor (au)
                     ----------------------------------------

                       Bx                  By                  Bz

  Bx         6.858131603456     -0.000000000000      0.000000000000
  By        -0.000000000000      6.858131603456      0.000000000000
  Bz         0.000000000000      0.000000000000     -0.000000000000
```


2. LDA, nc, go=(0,0,0) 

* input:

```
  **DALTON INPUT
  .RUN PROPERTIES
  **WAVE FUNCTIONS
  .DFT
  LDA
  *SCF INPUT
  .DOUBLY OCCUPIED
  5
  .THRESH
  1.0D-7
  .MAX MACRO ITERATIONS
  30
  **PROPERTIES
  .MAGNET
  .GAUGEO
   0.0 0.0 0.0
  .PRINT
  10
  *LINRES
  .PRINT
  10
  .THRESH
  1.0D-9
  .MAX IT
  100
  *EXPECT
  .PRINT
  10
  *END OF INPUT
```

* output:

```
  Magnetizabilities                    au          SI         cgs
  ---------------------------------------------------------------

  Isotropic magnetizability:      -2.1051   -166.1124    -10.0035
  Diamagnetic contribution:       -6.6772   -526.8975    -31.7305
  Paramagnetic contribution:       4.5721    360.7851     21.7270
  Parallel component:             -2.0301   -160.1926     -9.6470
  Perpendicular component:        -2.1426   -169.0722    -10.1818
  Anisotropy:                      0.1125      8.8796      0.5347

  Magnetizability is cylindrical.


  Principal values (au, SI, and cgs) and axes:
  --------------------------------------------

  khi_aa      -2.0301   -160.1926     -9.6470      0.0000    0.0000    1.0000
  khi_bb      -2.1426   -169.0722    -10.1818      0.9938    0.1111    0.0000
  khi_cc      -2.1426   -169.0722    -10.1818     -0.1111    0.9938    0.0000


                        Total magnetizability tensor (au)
                        ---------------------------------

                       Bx                  By                  Bz

  Bx        -2.142586056154      0.000000000000     -0.000000000000
  By         0.000000000000     -2.142586056154     -0.000000000000
  Bz        -0.000000000000     -0.000000000000     -2.030058045886



                     Diamagnetic magnetizability tensor (au)
                     ---------------------------------------

                       Bx                  By                  Bz

  Bx        -9.000717658320      0.000000000000     -0.000000000000
  By         0.000000000000     -9.000717658320     -0.000000000000
  Bz        -0.000000000000     -0.000000000000     -2.030058045886



                     Paramagnetic magnetizability tensor (au)
                     ----------------------------------------

                       Bx                  By                  Bz

  Bx         6.858131602166      0.000000000000      0.000000000000
  By         0.000000000000      6.858131602166      0.000000000000
  Bz         0.000000000000      0.000000000000     -0.000000000000
```


3. PBE, sc, go=(0,0,0) 

* input:

```
  **DALTON INPUT
  .RUN PROPERTIES
  **WAVE FUNCTIONS
  .DFT
  PBE
  *SCF INPUT
  .DOUBLY OCCUPIED
  5
  .THRESH
  1.0D-7
  .MAX MACRO ITERATIONS
  30
  **PROPERTIES
  .MAGNET
  .NODIFC
  .GAUGEO
   0.0 0.0 0.0
  .PRINT
  10
  *LINRES
  .PRINT
  10
  .THRESH
  1.0D-9
  .MAX IT
  100
  *EXPECT
  .PRINT
  10
  *END OF INPUT
```

* output:

```
  Magnetizabilities                    au          SI         cgs
  ---------------------------------------------------------------

  Isotropic magnetizability:      -2.0884   -164.7973     -9.9243
  Diamagnetic contribution:       -6.6564   -525.2585    -31.6318
  Paramagnetic contribution:       4.5680    360.4613     21.7075
  Parallel component:             -2.0178   -159.2224     -9.5886
  Perpendicular component:        -2.1237   -167.5847    -10.0922
  Anisotropy:                      0.1060      8.3623      0.5036

  Magnetizability is cylindrical.


  Principal values (au, SI, and cgs) and axes:
  --------------------------------------------

  khi_aa      -2.0178   -159.2224     -9.5886      0.0000    0.0000    1.0000
  khi_bb      -2.1237   -167.5847    -10.0922      0.0000    1.0000    0.0000
  khi_cc      -2.1237   -167.5847    -10.0922      1.0000    0.0000    0.0000


                        Total magnetizability tensor (au)
                        ---------------------------------

                       Bx                  By                  Bz

  Bx        -2.123734646838     -0.000000000000      0.000000000000
  By        -0.000000000000     -2.123734646838      0.000000000000
  Bz         0.000000000000      0.000000000000     -2.017762828074



                     Diamagnetic magnetizability tensor (au)
                     ---------------------------------------

                       Bx                  By                  Bz

  Bx        -8.975710005976     -0.000000000000     -0.000000000000
  By        -0.000000000000     -8.975710005976     -0.000000000000
  Bz        -0.000000000000     -0.000000000000     -2.017762828074



                     Paramagnetic magnetizability tensor (au)
                     ----------------------------------------

                       Bx                  By                  Bz

  Bx         6.851975359138     -0.000000000000      0.000000000000
  By        -0.000000000000      6.851975359138      0.000000000000
  Bz         0.000000000000      0.000000000000     -0.000000000000
```

4. PBE, nc, go=(0,0,0) 

* input

```
  **DALTON INPUT
  .RUN PROPERTIES
  **WAVE FUNCTIONS
  .DFT
  PBE
  *SCF INPUT
  .DOUBLY OCCUPIED
  5
  .THRESH
  1.0D-7
  .MAX MACRO ITERATIONS
  30
  **PROPERTIES
  .MAGNET
  .GAUGEO
   0.0 0.0 0.0
  .PRINT
  10
  *LINRES
  .PRINT
  10
  .THRESH
  1.0D-9
  .MAX IT
  100
  *EXPECT
  .PRINT
  10
  *END OF INPUT
```

* output

```
  Magnetizabilities                    au          SI         cgs
  ---------------------------------------------------------------

  Isotropic magnetizability:      -2.0884   -164.7973     -9.9243
  Diamagnetic contribution:       -6.6564   -525.2585    -31.6318
  Paramagnetic contribution:       4.5680    360.4613     21.7075
  Parallel component:             -2.0178   -159.2224     -9.5886
  Perpendicular component:        -2.1237   -167.5847    -10.0922
  Anisotropy:                      0.1060      8.3623      0.5036

  Magnetizability is cylindrical.


  Principal values (au, SI, and cgs) and axes:
  --------------------------------------------

  khi_aa      -2.0178   -159.2224     -9.5886      0.0000    0.0000    1.0000
  khi_bb      -2.1237   -167.5847    -10.0922      1.0000    0.0000    0.0000
  khi_cc      -2.1237   -167.5847    -10.0922      0.0000    1.0000    0.0000


                        Total magnetizability tensor (au)
                        ---------------------------------

                       Bx                  By                  Bz

  Bx        -2.123734648029     -0.000000000000     -0.000000000000
  By        -0.000000000000     -2.123734648029     -0.000000000000
  Bz        -0.000000000000     -0.000000000000     -2.017762828074



                     Diamagnetic magnetizability tensor (au)
                     ---------------------------------------

                       Bx                  By                  Bz

  Bx        -8.975710005976     -0.000000000000     -0.000000000000
  By        -0.000000000000     -8.975710005976     -0.000000000000
  Bz        -0.000000000000     -0.000000000000     -2.017762828074



                     Paramagnetic magnetizability tensor (au)
                     ----------------------------------------

                       Bx                  By                  Bz

  Bx         6.851975357947     -0.000000000000      0.000000000000
  By        -0.000000000000      6.851975357947      0.000000000000
  Bz         0.000000000000      0.000000000000     -0.000000000000
```



