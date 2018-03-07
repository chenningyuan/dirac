
Levy-Leblond results from Dirac:
================================

                                                                LDA,sc,go=(0,0,0)       LDA,nc,go=(0,0,0)    PBE,sc,go=(0,0,0)    PBE,nc,go=(0,0,0)

Isotropic                                                          -2.1051                   -2.1051            -2.0884               -2.0884                                       
Diamagnetic                                                        -6.6772                   -6.6772            -6.6564               -6.6564                                       
Paramagnetic                                                        4.5721                    4.5721             4.5680                4.5680                                       
Parallel                                                           -2.0301                   -2.0301            -2.0178               -2.0178                                       
Perpendicular                                                      -2.1426                   -2.1426            -2.1237               -2.1237                                       
Anisotropy                                                          0.1125                    0.1125             0.1060                0.1060                                       
                                                                                                                                               
One-electron expectation values (average)                           2.707990                  2.707990           2.667206              2.667206                                                  
   2-el + dft                                                      -0.718297                 -0.718297          -0.698589             -0.698589
Highest-order reorthonormalization (average)                       -0.109949                 -0.109949          -0.106654             -0.106654
Lowest-order reorthonormalization (average)                         5.325702                  5.482021           5.300339              5.466566
Relaxation (average)                                               -5.100370                 -5.256689          -5.073891             -5.240118
   lowest order + relaxation                                        0.225332                  0.225332           0.226448              0.226448


all values in [au]

sc = symmetric connection
nc = natural connection
go = gauge origin




Non-relativistic results from Dalton:
=====================================

                                                                LDA,sc,go=(0,0,0)       LDA,nc,go=(0,0,0)    PBE,sc,go=(0,0,0)    PBE,nc,go=(0,0,0)

Isotropic[au]                                                      -2.1051                   -2.1051            -2.0884               -2.0884
Diamagnetic[au]                                                    -6.6772                   -6.6772            -6.6564               -6.6564
Paramagnetic[au]                                                    4.5721                    4.5721             4.5680                4.5680
Parallel[au]                                                       -2.0301                   -2.0301            -2.0178               -2.0178
Perpendicular[au]                                                  -2.1426                   -2.1426            -2.1237               -2.1237
Anisotropy[au]                                                      0.1125                    0.1125             0.1060                0.1060

One-electron expectation values (average)                           2.707990                  2.707990           2.667206              2.667206
Two-electron expectation values (average)                          -0.780471                 -0.780471          -0.757909             -0.757909
DFT static contribution (average)                                   0.062175                  0.062175           0.059320              0.059320
    2-el + dft (to compare with DIRAC)                             -0.718296                 -0.718296          -0.698589             -0.698589                                               
Highest-order reorthonormalization (average)                       -0.109949                 -0.109949          -0.106654             -0.106654
Lowest-order reorthonormalization (average)                         5.325703                  5.482021           5.300339              5.466566
Relaxation (average)                                               -5.100370                 -5.256688          -5.073891             -5.240118
    lowest order + relaxation (to compare with DIRAC)               0.225333                  0.225333           0.226448              0.226448


all values in [au]

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
                         One-electron expectation values
                         -------------------------------

                       Bx                  By                  Bz

  Bx         3.046955661064      0.000000000000     -0.000000000000
  By         0.000000000000      3.046955661064     -0.000000000000
  Bz        -0.000000000000     -0.000000000000      2.030058045886


       Average value:    2.707990
  
  
                         Two-electron expectation values
                         -------------------------------

                       Bx                  By                  Bz

  Bx        -1.170707114458     -0.000000000000      0.000000000000
  By        -0.000000000000     -1.170707114458      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:   -0.780471
  
  
                             DFT static contribution
                             -----------------------

                       Bx                  By                  Bz

  Bx         0.093262489397      0.000000000000      0.000000000000
  By         0.000000000000      0.093262489397      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:    0.062175


                        Highest-order reorthonormalization
                        ----------------------------------

                       Bx                  By                  Bz

  Bx        -0.164923710296      0.000000000000      0.000000000000
  By         0.000000000000     -0.164923710296      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:   -0.109949


                        Lowest-order reorthonormalization
                        ---------------------------------

                       Bx                  By                  Bz

  Bx         7.988553800087      0.000000000000      0.000000000000
  By         0.000000000000      7.988553800087     -0.000000000000
  Bz         0.000000000000     -0.000000000000      0.000000000000


       Average value:    5.325703


                                    Relaxation
                                    ----------

                       Bx                  By                  Bz

  Bx        -7.650555070931      0.000000000000     -0.000000000000
  By         0.000000000000     -7.650555070931     -0.000000000000
  Bz        -0.000000000000     -0.000000000000      0.000000000000


       Average value:   -5.100370

  Units:  JT**(-2)/(10**(-30))   (SI)
          ppm cm**(3) mol**(-1)  (cgs)

 London orbitals (GIAOs) used.


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
                         One-electron expectation values
                         -------------------------------

                       Bx                  By                  Bz

  Bx         3.046955661064      0.000000000000     -0.000000000000
  By         0.000000000000      3.046955661064     -0.000000000000
  Bz        -0.000000000000     -0.000000000000      2.030058045886


       Average value:    2.707990


                         Two-electron expectation values
                         -------------------------------

                       Bx                  By                  Bz

  Bx        -1.170707114458     -0.000000000000      0.000000000000
  By        -0.000000000000     -1.170707114458      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:   -0.780471


                             DFT static contribution
                             -----------------------

                       Bx                  By                  Bz

  Bx         0.093262489397      0.000000000000      0.000000000000
  By         0.000000000000      0.093262489397      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:    0.062175


                        Highest-order reorthonormalization
                        ----------------------------------

                       Bx                  By                  Bz

  Bx        -0.164923710296      0.000000000000      0.000000000000
  By         0.000000000000     -0.164923710296      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:   -0.109949


                        Lowest-order reorthonormalization
                        ---------------------------------

                       Bx                  By                  Bz

  Bx         8.223031344582     -0.000000000000      0.000000000000
  By        -0.000000000000      8.223031344582      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:    5.482021


                                    Relaxation
                                    ----------

                       Bx                  By                  Bz

  Bx        -7.885032614136     -0.000000000000     -0.000000000000
  By        -0.000000000000     -7.885032614136     -0.000000000000
  Bz        -0.000000000000     -0.000000000000      0.000000000000


       Average value:   -5.256688

  Units:  JT**(-2)/(10**(-30))   (SI)
          ppm cm**(3) mol**(-1)  (cgs)

 London orbitals (GIAOs) used.


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
                         One-electron expectation values
                         -------------------------------

                       Bx                  By                  Bz

  Bx         2.991927636075      0.000000000000     -0.000000000000
  By         0.000000000000      2.991927636075     -0.000000000000
  Bz        -0.000000000000     -0.000000000000      2.017762828074


       Average value:    2.667206


                         Two-electron expectation values
                         -------------------------------

                       Bx                  By                  Bz

  Bx        -1.136863936237     -0.000000000000      0.000000000000
  By        -0.000000000000     -1.136863936237      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:   -0.757909


                             DFT static contribution
                             -----------------------

                       Bx                  By                  Bz

  Bx         0.088980074676      0.000000000000      0.000000000000
  By         0.000000000000      0.088980074676      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:    0.059320


                        Highest-order reorthonormalization
                        ----------------------------------

                       Bx                  By                  Bz

  Bx        -0.159980320227      0.000000000000      0.000000000000
  By         0.000000000000     -0.159980320227      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:   -0.106654


                        Lowest-order reorthonormalization
                        ---------------------------------

                       Bx                  By                  Bz

  Bx         7.950508339024      0.000000000000      0.000000000000
  By         0.000000000000      7.950508339024      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:    5.300339


                                    Relaxation
                                    ----------

                       Bx                  By                  Bz

  Bx        -7.610837146474     -0.000000000000     -0.000000000000
  By        -0.000000000000     -7.610837146474     -0.000000000000
  Bz        -0.000000000000     -0.000000000000      0.000000000000


       Average value:   -5.073891

  Units:  JT**(-2)/(10**(-30))   (SI)
          ppm cm**(3) mol**(-1)  (cgs)

 London orbitals (GIAOs) used.


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
                         One-electron expectation values
                         -------------------------------

                       Bx                  By                  Bz

  Bx         2.991927636075      0.000000000000     -0.000000000000
  By         0.000000000000      2.991927636075     -0.000000000000
  Bz        -0.000000000000     -0.000000000000      2.017762828074


       Average value:    2.667206


                         Two-electron expectation values
                         -------------------------------

                       Bx                  By                  Bz

  Bx        -1.136863936237     -0.000000000000      0.000000000000
  By        -0.000000000000     -1.136863936237      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:   -0.757909


                             DFT static contribution
                             -----------------------

                       Bx                  By                  Bz

  Bx         0.088980074676      0.000000000000      0.000000000000
  By         0.000000000000      0.088980074676      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:    0.059320


                        Highest-order reorthonormalization
                        ----------------------------------

                       Bx                  By                  Bz

  Bx        -0.159980320227      0.000000000000      0.000000000000
  By         0.000000000000     -0.159980320227      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:   -0.106654


                        Lowest-order reorthonormalization
                        ---------------------------------

                       Bx                  By                  Bz

  Bx         8.199848445275     -0.000000000000      0.000000000000
  By        -0.000000000000      8.199848445275      0.000000000000
  Bz         0.000000000000      0.000000000000      0.000000000000


       Average value:    5.466566


                                    Relaxation
                                    ----------

                       Bx                  By                  Bz

  Bx        -7.860177251533      0.000000000000     -0.000000000000
  By         0.000000000000     -7.860177251533     -0.000000000000
  Bz        -0.000000000000     -0.000000000000      0.000000000000


       Average value:   -5.240118

  Units:  JT**(-2)/(10**(-30))   (SI)
          ppm cm**(3) mol**(-1)  (cgs)

 London orbitals (GIAOs) used.



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



