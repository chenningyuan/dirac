Magnetizabilities (magnetic susceptibilities)
=============================================

Relativistic calculation of magnetizabilities using gauge-including atomic orbitals (GIAO),
preferably labeled as London atomic orbitals (LAO) in our cirles.

Testing the effect of moving the gauge origin (GO), effects of relativity by modifying the speed of light
and other parameters.

RnH(77+) molecule
-----------------

We have chosen heavy diatomic, RnH(77+), the HF-like molecule, with Gaussian type of nucleus.
Small, 'sucked-out-of-finger' decontracted basis set to have MOs of this system converged.
All values in ppm.

             
Method/Atom                                    
---------------------------------------------------------------
ukb-ep-lao-go_cm 
 ee+ep-rn, LL+SL         0.032
 ee+ep-rn, LL+SL+SS 
       
ukb-ep-lao-go_1  
    ee+ep-rn    
    ee+ep-rn,symcon  
   
rkbimp-ep-lao-go_cm, ee+ep-rn:
   LL+SL                            0.0137    -0.0253
   LL+SL+SS   

rkb-ep-lao-go_cm
    ee-rn           
    ee+ep-rn        
rkb-ep-lao-go_1    
    ee-rn           
    ee+ep-rn        

rkb-ee-cgo_go_cm   
rkb-ee-cgo_go_1    

rkb-ee-lao-go_1     
rkb-ee-lao-go_cm             -0.0047         
 cval=2000           
 cval=4000           
 cval=10000          
 cval=15000 (linresp very slowly converging)
 cval=30000 (linresp very slowly converging)

levyle-lao-go_cm   
levyle-lao-go_1     

lel-pnuc-lao-go_cm   

dal-pnuc-lao-go_cm    
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

Method/Atom                  isotropic   anisotropic                
---------------------------------------------------------------------------
rkb-ep-lao-go_cm, ee+ep-rn:                    (all methods with LL+SL) 
                                               (all methods with LL+SL+SS)

rkbimp-ep-lao-go_cm, ee+ep-rn: 
2-Fock regener. LL+SL+SS,MO-LL+SL              (dg(1,2)/dB LL+SL) 
                                               (dg(1,2)/dB LL+SL+SS)

2-Fock regen. LL+SL+SS,MO-LL+SL+SS             (dg(1,2)/dB LL+SL+SS)
                                               (dg(1,2)/dB LL+SL)

2-Fock regener. LL+SL,MO-LL+SL                 (dg(1,2)/dB LL+SL)
                                               (dg(1,2)/dB LL+SL+SS)

2-Fock regener. LL+SL,MO-LL+SL+SS              (dg(1,2)/dB LL+SL)


ukb-ep-lao-go_cm, ee+ep-rn:                    (all methods with LL+SL)
                                               (all methods with LL+SL+SS)


---------------------------------------------------------------

