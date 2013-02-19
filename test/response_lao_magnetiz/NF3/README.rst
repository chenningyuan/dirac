Magnetizabilities of the NF3 molecule
=======================================

With London atomic orbitals (LAO) and without - common gauge origin (CGO).

Used small STO-2G decontracted basis sets for N and F, respectivelly. 
Point nucleus model.

All values in a.u.

                        isotropic         anisotropic
-------------------------------------------------------------------

                            nonrelativistic results

LevyLe-LAO(go=cm)        -5.2733             0.5431 (nonrelativ. RDSUSLL form)
LevyLe-CGO(go=cm)       -26.1709           -14.6105
LevyLe-CGO(go=1)        -260.3916           -183.0540(1st),357.1151(2nd)

Dalton-LAO(go=cm)        -5.2733             0.5431
Dalton-LAO(go=1 )        -5.2733             0.5431
Dalton-LAO(go=1 )        -5.2733             0.5431 (symm.connection)

Dalton-CGO(go=cm)       -26.1709            -14.6105
Dalton-CGO(go=1)        -260.3916           -183.0540(1st),357.1151(2nd)

LevyLe-LAO(go=cm)        -5.8363             0.4126 (full RDSUSLL form)

                           relativistic results

RKB-ee-LAO+dia(go=cm)    -5.8300             0.4137
RKB-ee-LAO+dia(go=cm), increasing speed of light to get nonrelativ.limit 
   cvalue=2000           -5.8363             0.4126 (full RDSUSLL form)
                         -5.2733             0.5431 (nonrelativ. RDSUSLL form with factor 2)
RKB-ee-LAO+dia(go=1)
  cvalue=2000            -5.2733             0.5431 (nonrelativ. RDSUSLL form with factor 2)
RKB-ee-LAO+dia(go=1),SymCon
  cvalue=2000            -5.2732a            0.5431  (nonrelativ. RDSUSLL form with factor 2)

RKB-ee-CGO+dia(go=cm)   -26.1706            -14.6131

RKB-ee-CGO+dia(go=cm), increasing speed of light to get nonrelativ.results
   cvalue=300           -26.1708            -14.6110
   cvalue=1000          -26.1709            -14.6105
   cvalue=2000          -26.1709            -14.6105 

RKB-ep-LAO(go=cm)        
   LL+SL                -3.1852              0.4392
   LL+SL+SS             -3.1852              0.4392

RKB-ep-CGO(go=cm)        
   LL+SL                -10.4296             -6.0383

RKBIMP-ep-LAO(go=cm)
   LL+SL                -4.8265              0.4864 (C1 symmetry)
   LL+SL                -4.8264              0.4863 (automatic symmetry detection)

RKBIMP-ep-LAO(go=cm),LL+SL
   cvalue=2000          -4.8323              0.4856
   cvalue=2000          -4.2693              0.6161 (nonrel.for of RDSUSLL)
RKBIMP-ep-LAO(go=cm)+diamag,LL+SL <-- wrong, ep terms are confused
   cvalue=2000          -10.4380             0.8024 (nonrel.for of RDSUSLL)

UKB-ep-LAO(go=cm)
   LL+SL                -5.0942              0.3767
   LL+SL+SS             -5.0942              0.3767
   LL+SL,cval=2000      -5.1003              0.3758
   LL+SL,cval=5000      -5.1003              0.3758

UKB-ee-LAO(gp=cm)+diam.
   LL+SL                -5.9832              0.3019       
   LL+SL,cval=2000      -5.9906              0.3003

RKBIMP-ep-CGO(go=cm)  
   LL+SL                -22.3284             -12.7181 (automatic symm.detection)
   LL+SL                -22.3284             -12.7181 (C1 symmetry)                

UKB-ep-CGO(go=cm)
   LL+SL                -25.7754             -14.6864
   LL+SL+SS             -25.7754             -14.6864

         OpenRSP-ee-lao(Go=cm)+diamag; C1 symm. :
London Magnetizability = Ebb w:  0.000  0.000
 -5.844574004         0.160554026         0.100120776
  0.160554026        -5.758773468         0.130385923
  0.100120776         0.130385923        -5.886553460
         RKB-ee-lao(go=cm); C1 symmetry :
 -5.844570248197      0.160558500304      0.100123566256
  0.160558500240     -5.758767322001      0.130389556087
  0.100123566384      0.130389556333     -5.886550873651

----------------------------------------------------------------------------------
go=cm  : gauge orgin (go) is the center of mass
go=1   : go is  (-3.00,4.00,-5.00)
LL+SL, LL+SL+SS : selected classes of integrals in all methods
"ep" label: LAO reorthonormalization sums are over both negative and positive states
