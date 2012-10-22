NMR shielding constants
=======================

Relativistic calculation of NMR shielding constants using gauge-including atomic orbitals (GIAO),
preferably labeled as London atomic orbitals (LAO) in our cirles.

Testing the effect of moving the gauge origin (GO), effects of relativity by modifying the speed of light
and other parameters.

RnH(77+) molecule
-----------------

We have chosen heavy diatomic, RnH(77+), the HF-like molecule, with Gaussian type of nucleus.
Small, 'sucked-out-of-finger' decontracted basis set to have MOs of this system converged.
All values in ppm.

             
Method/Atom               Rn                H
---------------------------------------------------------------
ukb-ep-lao-go_cm 
    ee+ep-rn         10322.7416          20.7567
ukb-ep-lao-go_1  
    ee+ep-rn         10322.7416          20.7567
    ee+ep-rn,symcon  10322.7416          20.7567
   
rkbimp-ep-lao-go_cm, ee+ep-rn:
  2-Fock LL+SL+SS   10322.7416          20.7567
  2-Fock LL+SL

rkb-ep-lao-go_cm
    ee-rn            9515.5269           14.5035
    ee+ep-rn         9515.5266           14.5036
rkb-ep-lao-go_1    
    ee-rn            9515.5269           14.5035
    ee+ep-rn         9515.5266           14.5036

rkb-ee-cgo_go_cm     14060.9990           2.7921
rkb-ee-cgo_go_1      14315.9857         -45.0504

rkb-ee-lao-go_1      14026.1843          24.0643
rkb-ee-lao-go_cm     14026.1843          24.0643
 cval=2000            5795.7876          24.2080
 cval=4000            5769.9413          24.2080
 cval=10000           5762.5390          24.2080
 cval=15000 (linresp very slowly converging)
 cval=30000 (linresp very slowly converging)

levyle-lao-go_cm      5761.1176          24.2080
levyle-lao-go_1       5761.1176          24.2080

lel-pnuc-lao-go_cm    5761.5661          24.2080

dal-pnuc-lao-go_cm    5761.5661          24.2080
dal-pnuc-lao-go_1
-------------------------------------------------------------------------------------
go_cm : center of mass used as gauge origin
go_1  : (3.0 -5.0 4.0), shifted go
ee-rn : reorthonormalization summations over electronic shells only
ee+ep-rn : reorthonormalization summations over both electronic and positronic shells
-------------------------------------------------------------------------------------

Some observations: when the system is too heavy, it's quite difficult to get 
nonrelativistic data by extrapolatig speed of light.

In Dalton there is no nuclear radius for finite nucleus of Rn, we switched to point nucleus.

Rn2(154+) molecule
------------------

Another diatomic, Rn2(154+), the F2-like molecule, was chosen to simulate properties
of relativistic systems with two heavy elements, especially effect of SS-class of integrals,
which is numerically not observable for a molecule composed of heavy-atom -- light-atom.

Method/Atom                            Rn                
---------------------------------------------------------------------------
rkb-ep-lao-go_cm, ee+ep-rn:         8830.5190  (all methods with LL+SL) 
                                    8831.0659  (all methods with LL+SL+SS)

rkbimp-ep-lao-go_cm, ee+ep-rn: 
2-Fock regener. LL+SL+SS,MO-LL+SL   9635.3832  (dg(1,2)/dB LL+SL) 
                                    9635.3928  (dg(1,2)/dB LL+SL+SS)

2-Fock regen. LL+SL+SS,MO-LL+SL+SS  9637.0731  (dg(1,2)/dB LL+SL+SS)
                                    9637.0636  (dg(1,2)/dB LL+SL)

2-Fock regener. LL+SL,MO-LL+SL      9636.6726  (dg(1,2)/dB LL+SL)
                                    9636.6822  (dg(1,2)/dB LL+SL+SS)

2-Fock regener. LL+SL,MO-LL+SL+SS   9633.8154  (dg(1,2)/dB LL+SL)


ukb-ep-lao-go_cm, ee+ep-rn:         9637.3930  (all methods with LL+SL)
                                    9637.7873  (all methods with LL+SL+SS)


---------------------------------------------------------------

