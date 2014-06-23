B+ system, 2p shells active
===========================

what works: 
-----------

         01 sector
sto-2g uncontracted;  -0.220068052835 (deg 6) - RELCCSD,01
(4s2p)                -0.220068080000 (   6) - hsfsccsd; D2h
                      -0.220068136469 (4)
                      -0.220062425815 (2) - hsfsccsd,C2v

does not work:
---------------

sto-2g contracted (2s1p)
Active space 2p, active electrons 1s(2)2s(0-1)2p(0)

             -0.160875456076 (6) - RELCCSD, C2v

 broken:   energy  1 -0.229806545146E+00. degen  2. 
           energy  2 -0.161171173231E+00. degen  2.   C2v
           energy  3 -0.136130558509E+00. degen  2.

 energy  1 -0.220388266182E+00. degen  2. 
 energy  2 -0.158274757401E+00. degen  2.      D2h
 energy  3 -0.148228693156E+00. degen  2

sto-2g contracted (2s1p)
------------------------
1s(2) frozen, active space 2p, active electrons 2s(0-1)2p(0)
-------------------------------------------------------------------

C1:

  (0,1) state     -0.160989772304  (   6 * )
        energy  1 -0.255082386176E+00. degen  2. 
        energy  2 -0.164213047207E+00. degen  2. 
        energy  3 -0.160989773020E+00. degen  2.

C2v:
    CCSD correlation energy : -0.062614099864741
    correlation energy(final)=-0.06261409771666
   (0,1) state    -0.160989772304 (   6 * )

 energy  1 -0.255082386176E+00. degen  2.
 energy  2 -0.164213047207E+00. degen  2.
 energy  3 -0.160989773020E+00. degen  2.

D2h:
       CCSD correlation energy :   -0.062614099864742
      correlation energy(final)=   -0.06261409771666

 energy  1 -0.255082386176E+00. degen  2
 energy  2 -0.164213047207E+00. degen  2
 energy  3 -0.160989773020E+00. degen  2

