      SUBROUTINE lambda_equation_T1 (FVO,EPS,T1,T2,HO,HV,HOV,H,CZ,GOO,
     &                   GVV,BUF1,BUF2,NBUF2,BUF3,NBUF3,S1)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Evaluate       Lambda1 equations
C     Here we have to take care of the fact that certain integral classes
C     are now split over individual nodes. The corresponding contributions
C     which are linear in the integrals are added at the end of the T1
C     routines. Therefore we have to take care of calculating the
C     non-distributed contributions to the T1 amplitudes only on the
C     master node in order to avoid double-counting in the PARTS routine !
C     all direct VOVV/VVVV contributions and VOVV/VVVV contributions stemming
C     from intermediates have to be calculated on EACH node.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER NBUF2,NBUF3
      REAL*8 EPS(*),FVO(*)
      REAL*8 T1(*),T2(*)
      REAL*8 S1(*)
      REAL*8 HO(*),HV(*),HOV(*),H(*)
      REAL*8 GOO(*),GVV(*),CZ(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "ihm.inc"
#include "complex.inc"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
      LOGICAL DONE,TEQ,USEDZ,RIGHT
      REAL*8 SDOT,DDOT
C
C---------------Executable code--------------------------------------
C
                                
                     
                                    
                      
                                    
                     
          
                                                     
           
C-------------------
C Initialize to zero
C-------------------
      CALL XCOPY (NDIMT1,A0,0,S1,1)
                                 
                                 
                                 
                   
                      
                      
                      
                                  
                                                             
                                                             
                                                             
                                                         
                                       
                                               
                                           
                                     
                 
                         
                       
                       
                        
                        
                                                                     
                                    
                                                             
                                       
                                     
              
                                                         
             
                                                      
                                                      
                                                      
                                                                  
                                         
                                                                  
                                                           
                                                      
                                                   

                            
                            
                            
                      
                                                                 
                                                           
           

C ********************
C ** TASK SPLITTING **
C ********************

      IF (MYPROC.EQ.MASTER) THEN

      CALL XTIME(0,1,'--- T1EQNS   HV*T / T*HO              ')
C-----------------------------------------------
C S(A,I) = S(A,I) + HV(A,C) * T(C,I)
C-----------------------------------------------
      CALL CNTRCT ('N','N',NV,NO,NV,A1,HV,T1,A1,S1,NREP)
C-----------------------------------------------
C S(A,I) = S(A,I) - T(A,K) * HO(K,I)
C-----------------------------------------------
      CALL CNTRCT ('N','N',NV,NO,NO,-A1,T1,HO,A1,S1,NREP)
      CALL XTIME(0,2,'--- T1EQNS   HV*T / T*HO            ')

C--------------------------------------------
C S(A,I) = S(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C--------------------------------------------
      CALL XTIME(0,1,'--- T1EQNS VOOO*TAU                ')
                   
                                 
          
         CALL XCOPY (NDIMT2,T2,1,BUF3,1)
           
      CALL T1EQN1 (BUF1,BUF2,BUF3,NBUF3,S1)
      CALL XTIME(0,2,'--- T1EQNS VOOO*TAU               ')

C ****************************************************
C ** END OF TASK SPLITTING. THE COMMON STUFF STARTS **
C ****************************************************

      ENDIF

                                             
                                                                      
                                             
                   
                                 
          
C -------------------------------------------
C TAU'(CD,IK) = L(CD,IK) +  G(C,B) * T(BD,IK)   <<<<<<------- LAMBDA EQNS
C -------------------------------------------
         CALL GETAMPT (BUF1,BUF2)
         CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &                BUF2,BUF1)
         CALL CNTRCT ('N','N',NV,NVOOT,NV,A1,GVV,BUF1,A0,BUF2,NREP)
         CALL SRT1T2 (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &                BUF2,BUF1)
         CALL XAXPY (NDIMT2,A1,T2,1,BUF1,1)
           
      CALL XTIME(0,1,'--- T1EQNS VOVV contribution  ')
      CALL T1EQN2 (BUF1,BUF2,BUF3,NBUF3,S1)   ! VOVV contractions, term 6 !
      CALL XTIME(0,2,'--- T1EQNS VOVV contribution  ')

C--------------------------------------
C S(A,I) = S(A,I) - W(A,K,C,I) * T(C,K)    ! VOVO contributions
C--------------------------------------
      CALL XTIME(0,1,'--- T1EQNS VOVO * T(C,K)         ')
                   
                            
                                               
                                                                        
                              
          
         CALL XCOPY (JVOVO(2),H,1,BUF2,1)
           

      M    = MVO(1)
      K    = IDIST(1,3,1)
      OFF1 = IDIST(2,3,1)
      OFF2 = IDIST(3,3,1)
      CALL XGEMV ('N',M,K,-A1,BUF2(OFF1),M,T1(OFF2),1,A1,S1,1)
      CALL XTIME(0,2,'--- T1EQNS VOVO * T(C,K)        ')
C------------------------------------
C Terms specific for Lambda1 equations
C------------------------------------
                        
C----------------------------------------------------
C S(A,I) = S(A,I) - T(A,L) * V(K,L,C,D) * L(C,D,K,I)
C        = S(A,I) - T(A,L) * V*(C,D,K,L) * L(C,D,K,I)
C----------------------------------------------------
      CALL GETVVOO(BUF2)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             BUF2,BUF1)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             T2,BUF2)
      CALL CNTRCT ('C','N',NO,NO,NVVOT,A1,BUF1,BUF2,A0,BUF3,NREP)
      CALL GETAMPT (BUF1,BUF2)
      CALL CNTRCT ('N','N',NV,NO,NO,-A1,BUF1,BUF3,A1,S1,NREP)
C----------------------------------------------------
C S(A,I) = S(A,I) - L(A,C,K,L) * V(K,L,D,C) * T(D,I)
C        = S(A,I) - L(A,C,K,L) * V*(D,C,K,L) * T(D,I)
C----------------------------------------------------
      CALL GETVVOO(BUF2)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             BUF2,BUF1)
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             T2,BUF2)
      CALL CNTRCT ('N','C',NV,NV,NVOOT,A1,BUF2,BUF1,A0,BUF3,NREP)
      CALL GETAMPT (BUF1,BUF2)
      OFF1 = 1
      OFF2 = 1
      OFF3 = 1
      DO IRP =1 , NREP
         M = NV(IRP)
         N = NO(IRP)
         K = NV(IRP)
         CALL XGEMM ('N','N',M,N,K,-A1,BUF3(OFF1),M,BUF1(OFF2),K,
     &               A1,S1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF2 = OFF2 + K * N * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
C--------------------------------------
C S(A,I) = S(A,I) - W(C,A,D,I) * G(C,D)
C--------------------------------------
      M = NFVV
      N = NFVO
      IF (NBUF3.LT.JVVVO(2)) THEN
         CALL QUIT(" OUT OF CORE STEP IN GV NOT IMPLEMENTED")
      ENDIF
      USEDZ = .TRUE.
      RIGHT = .FALSE.
clv : changed to SRT20D, untested as the rest of the lambda code
clv : addition should be supplied with a minus sign
      CALL SRT20D (NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
     &             JJVO,BUF2,NBUF2,S1,GVV,USEDZ,RIGHT)
C--------------------------------------
C S(A,I) = S(A,I) + W(A,K,L,I) * G(K,L)
C--------------------------------------
      CALL SRT1ST4 (NREP,MULTB,LTR,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
     &              CZ,BUF2)
      M = NFVO
      K = NFOO
      CALL XGEMV ('N',M,K,A1,BUF2,M,GOO,1,A1,S1,1)
C---------------------------------------------------
C S(A,I) = S(A,I) + G(C,B) * W(B,K,C,D) * T(A,D,K,I)
C---------------------------------------------------
      USEDZ = .FALSE.
      RIGHT = .FALSE.
      CALL SRT20D (NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
     &            JJVO,BUF2,NBUF2,BUF1,GVV,USEDZ,RIGHT)
      CALL GETAMPT (BUF3,BUF2)
      CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,BUF2,BUF3)
      M = MVO(1)
      K = MVO(1)
      CALL XGEMV ('N',M,K,A1,BUF3,M,BUF1,1,A1,S1,1)
           
C
            
      END
