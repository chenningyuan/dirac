      SUBROUTINE cceqn_driver_lambda(MAXIT,EPS,FOO,FVO,FVV,T1,S1,T2,S2,
     &                 AZ,CZ,H,HO,HV,HOV,GOO,GVV,
     &                 BUF1,BUF2,BUF3,CC1,BB1,BB2,CONV,
     &                 NBUF2,NBUF3                              )

      use interface_to_mpi
      use xmlout
      use intermediates_1b_2b
      IMPLICIT INTEGER (A-Z)

C---------------Description--------------------------------------------
C
C     Solves CCSD amplitude or lambda equations.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C     Intermediaries etc.

      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*),T1(*),S1(*),T2(*),S2(*)
      REAL*8 AZ(*),H(*),HO(*),HV(*),HOV(*),BUF1(*),BUF2(*),BUF3(*)
      REAL*8 CZ(*) ! is only allocated for Lambda equations
      REAL*8 GOO(*),GVV(*)

C     for diis ...

      REAL*8 CC1(*)
      REAL*8 BB1(*),BB2(*)

      REAL*8 RMS,CONV
      REAL*8 DIAG

C     for Fock space ...

                                                     

C---------------Common Blocks--------------------------------------

#include "files.inc"
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "ihm.inc"
#include "results.inc"
#include "complex.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
      DIMENSION ISTAT(df_MPI_STATUS_SIZE)
#endif

C---------------Local variables--------------------------------------

      REAL*8 CPU,CPU4V,CPUH1,CPUH2
      LOGICAL     FOCKSP,PRENE
      CHARACTER*2 DELS1
      CHARACTER*4 DELS2
      REAL*8 DDOT

C---------------Executable code--------------------------------------
C
C
C     Get initial values of timers for the most important contractions
C     This is important in Fock space calculations where we enter this
C     routine more than once. We do it in this manner rather than resetting
C     the timer because we want the overall timing at the end.

      CALL XTIMRD (CPU4V,'-- BINTM                       ')
      CALL XTIMRD (CPUH1,'--- HINTM: VVOO contribution   ')
      CALL XTIMRD (CPUH2,'--- T2EQNS: HINTM*T2           ')

!     CALL cceqn_driver_init(MAXIT,FVO,T1,S1,T2,S2,
!    &                 BUF1,BUF2,BUF3,CC1,BB1,BB2,
!    &                 NBUF2,NBUF3,RMS,DELS1,DELS2)

! this is the start of what the init above should replace
      IF (EQNS.EQ.'LAMBDA') THEN
         FOCKSP = .FALSE.
      ELSE
         CALL QUIT('Wrong keyword in call to CCEQNS')
      ENDIF

C     Zero amplitudes if the NOCCX options are specified

                                                      
                                                      

C     Initialize the work arrays to zero.

      CALL XCOPY (NDIMT2,A0,0,BUF1,1)
      CALL XCOPY (NBUF2,A0,0,BUF2,1)
      CALL XCOPY (NBUF3,A0,0,BUF3,1)

C     Initialize DIIS.

      IT = -1
      CALL DIISX(NDIMT1,NDIMT2,BUF1,T1,BUF2,T2,
     &           IT,CC1,BB1,ITC,BB2,RMS)

      ECCSD = 0.D0
      NIT = 0

      CALL T1DIAG(DIAG,T1,BUF1)
      CALL CCENG (T1,T2,FVO,BUF1,BUF2,ECCSD,ECCSDIM)
      IF (IPRNT.GE.3) CALL PAIREN (T1,T2,FVO,BUF1,BUF2)

      IF (IPRNT.GE.1) THEN
         WRITE(IW,6003)
         WRITE(IW,6001) NIT,ECCSD,RMS,DIAG
      ENDIF
! this is the start of what the init above should replace


!     Calculate fixed intermediaries for Lambda equations

                        
         CALL GETAMPT (S1,S2)
!         CALL HOINTM (FOO,EPS,S1,S2,BUF1,BUF2,BUF3,HO)

         call f_bar_mi(foo,fvo,s1,s2,ho)

         write(*,*)'f_bar_mi intermediate'
         write(*,*)'---------------------'
         write(*,*)(ho(i),i=1,nfoo)


         CALL HVINTM (FVV,EPS,S1,S2,BUF1,BUF2,BUF3,HV)
         CALL HOVINTM (FVO,S1,BUF1,BUF2,HOV)
         CALL GOINTM (FVO,S1,BUF1,BUF2,HO)
         CALL GVINTM (FVO,S1,BUF3,NBUF3,BUF2,NBUF2,HV)
         CALL AINTM (S1,S2,BUF1,BUF2,AZ)
         CALL CINTM (T1,T2,FVO,AZ,BUF1,BUF2,BUF3,NBUF3,CZ)
!avijit         CALL DINTM (S1,S2,HOV,BUF1,BUF2,NBUF2,BUF3,NBUF3,H)
!avijit         CALL HINTM (S1,S2,BUF1,BUF2,BUF3,NBUF3,H)
           

C
C  START ITERATIVE PROCEDURE
C
 1000 CONTINUE
      NIT = NIT + 1
      IF (NIT.GT.MAXIT) THEN
         NIT = NIT - 1
         WRITE(IW,6004) RMS
!        GO TO 2000
      ENDIF
C
                        
C
C        Make Gpq intermediates for lambda equations
C
         CALL GETAMPT (S1,S2)
         CALL MKGPP (S2,T2,BUF1,BUF2,GOO,GVV)
           

C  EVALUATE T1 EQUATIONS
C  CHECK FOR CASES WITH NO NON_ZERO T1'S TO AVOID CRASHES
C
      CALL XTIME(5,1,'- L1 equations                  ')
                   
                                                        
           
                   
                                                        
           
      IF (NDIMT1.GT.0) THEN
                      
                                             
                                                            
                                                     
             
          CALL lambda_equation_T1(HOV,EPS,T1,T2,HO,HV,HOV,H,CZ,GOO,GVV,
     &                 BUF1,BUF2,NBUF2,BUF3,NBUF3,S1)
              
      ENDIF
      CALL XTIME(5,2,'- L1 equations                  ')

C
C     In case of NOCCD we can skip the calculation of T2 amplitudes
C
                                 
      CALL XTIME(5,1,'- L2 equations                  ')


                             
                     
                                                           
                                                           
                                                               
                                                             
           
 
                            
 
                   
                                                            
                                          
                                                            
                                                            
                                                      
                                                            
                                                            
                                        
                                                            
                                                            
                                                  
                                                            

                     
                     
                                                               
                                                         
                                                                
                                                                    
                                                                   
 
                             

                                  
                                                             

                            

                                                              
                                      
                                                 
                
                                                        
                               
                                            
                          
                                                              
             
                                            
                                                              

           

      

           

 
                     
                                                           
                                                           
                              
                                                               
                            
                                                        
           
C
C  INTERMEDIATES FORMED --- START T2 EQUATIONS
C
      CALL XTIME(5,1,'-- L2 EQNS                     ')
      call lambda_equation_T2 (EPS,T1,T2,HO,HV,HOV,H,AZ,CZ,GOO,GVV,
     &             BUF1,BUF2,BUF3,NBUF3,S2)
      CALL XTIME(5,2,'-- T2 EQNS                     ')
 
 
                                              
 
                     
                                                             
           
C
      CALL XTIME(5,1,'-- BINTM                       ')
      CALL BINTM (T1,T2,BUF1,BUF2,NBUF2,BUF3,NBUF3,S2,IW)
      CALL XTIME(5,2,'-- BINTM                       ')

                     
                                                                        
           
C
      CALL XTIME(5,2,'- L2 equations                  ')
C
              
 
                                                                     
                                                 
                                        
           
 
                                                    
 
                                           
 
                                        
 
                                                          
                                        
                                                          
 
                     
                                                                       
                                                                       
           
 
                      
                                                               
                                                             
                           
                       
                                                                        
                                                                        
             
           
 
                                                            
 
                                   
                                   
 
                      
                                       
                                                                      
                                                                   
                                       
          
                                    
           
 
                                                        
 
                                                      
                                                      
 
                        
 
                                                          
                           
                                                           
                                                  
 
                             
 
                                                        
                                 
                                         
                                               
                                        
                                                        
 
                                        
 
                                                          
                                        
                                                          
 
                                                    
 
                   
                                                          
                                                        
           
 
                                
                                                       
                     
 
 
                                
                
 
              
                                
                    
 
                                                               
                    
 
                           
                     
          
                               
           
 
                     
                            
                            
                                        
                                
             
                                           
                                           
              
                            
                                                          
           
 
                         
                        
                        
                     
                                                    
                                                    
                                                    
            
C
                      
                                               
                                                               
                                                  
           
C
C     Read out timers for the most important contractions
C
      CALL XTIMRD (CPU,'-- BINTM                       ')
      CPU4V = CPU - CPU4V
      CALL XTIMRD (CPU,'--- HINTM: VVOO contribution   ')
      CPUH1 = CPU - CPUH1
      CALL XTIMRD (CPU,'--- T2EQNS: HINTM*T2           ')
      CPUH2 = CPU - CPUH2
C
C     Compute & print the performance in Mega/Giga/TeraFlops
C
      IF (TIMING) CALL BFLOPS (CPU4V,CPUH1,CPUH2,NIT)
C
C     Close and delete DIIS scratch files
C
      IT = -2
      CALL DIISX(NDIMT1,NDIMT2,BUF1,T1,BUF2,T2,
     &           IT,CC1,BB1,ITC,BB2,RMS)
C
 111  FORMAT(A20,I4,2f25.15)
 6001 FORMAT(2X,I3,5X,2F20.15,F10.5)
 6002 FORMAT(3X,'MP2   =',F20.15/3X,'EMP2  =',F20.15/)
 6003 FORMAT(//3X,'NIT',10X,'ENERGY',20X,'RMS',7X,'T1-DIAGN')
 6004 FORMAT(//' WARNING: Requested convergence not reached in CCSD :',
     &         'RMS = ',F20.15)
 6071 FORMAT(//,'  CCSD results',/)
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6074 FORMAT(' CCSD correlation energy :',T40,F25.15)
 6084 FORMAT(' CCD correlation energy :',T40,F25.15)
 6085 FORMAT(' CCS correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6076 FORMAT(' Total CCSD energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 6078 FORMAT(' Convergence :',T40,F25.15)
 6079 FORMAT(' Number or iterations used :',T60,I5)
      RETURN
      END
