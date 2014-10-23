      SUBROUTINE cceqn_driver_amplitudes(MAXIT,EPS,FOO,FVO,FVV,
     &                 T1,S1,T2,S2,AZ,CZ,H,HO,HV,HOV,GOO,GVV,
     &                 BUF1,BUF2,BUF3,CC1,BB1,BB2,CONV,
     &                 NBUF2,NBUF3,HEFF,EIGENVALUES,EIGENVECTORS)

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

      REAL*8 HEFF(*),EIGENVECTORS(*),EIGENVALUES(0:*)

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
      IF (EQNS.EQ.'AMPLTD') THEN
         FOCKSP = .FALSE.
      ELSEIF (EQNS.EQ.'FOCKSP') THEN
         FOCKSP = .TRUE.
         DELS1 = 'DD'
         IF (ICURFSS.EQ. 2) DELS1 = 'DK'
         IF (ICURFSS.EQ. 3) DELS1 = 'KD'
         IF (ICURFSS.EQ. 4) DELS1 = 'KK'
         IF (ICURFSS.EQ. 5) DELS1 = 'DK'
         IF (ICURFSS.EQ. 6) DELS1 = 'KD'
         DELS2 = 'DDDD'
         IF (ICURFSS.EQ. 2) DELS2 = 'DDDK'
         IF (ICURFSS.EQ. 3) DELS2 = 'KDDD'
         IF (ICURFSS.EQ. 4) DELS2 = 'KDDK'
         IF (ICURFSS.EQ. 5) DELS2 = 'DDKK'
         IF (ICURFSS.EQ. 6) DELS2 = 'KKDD'

C        Initialize S1 and S2 in case we want to do zero
C        iterations (tests). In normal runs S1 and S2 will
C        be initialized in T1EQNS and T2EQNS.

         CALL XCOPY(NDIMT1,A0,0,S1,1)
         CALL XCOPY(NDIMT2,A0,0,S2,1)
         CALL S1S2_ADDFO (FVO,BUF1,BUF2,S1,S2)

      ELSE
         CALL QUIT('Wrong keyword in call to CCEQNS')
      ENDIF

C     Zero amplitudes if the NOCCX options are specified

      IF (        NOCCS) CALL XCOPY (NDIMT1,A0,0,T1,1)
      IF (        NOCCD) CALL XCOPY (NDIMT2,A0,0,T2,1)

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
! this is the end   of what the init above should replace
                     
                                          
                                   

                                                      
                                            
                                          
                                                      
                                        
                                                          
                                                                   
                                                         
           

C
C  START ITERATIVE PROCEDURE
C
 1000 CONTINUE
      NIT = NIT + 1
      IF (NIT.GT.MAXIT) THEN
         NIT = NIT - 1
         WRITE(IW,6004) RMS
         GO TO 2000
      ENDIF
C
                        
 
                                                    
 
                             
                                             
           

C  EVALUATE T1 EQUATIONS
C  CHECK FOR CASES WITH NO NON_ZERO T1'S TO AVOID CRASHES
C
      CALL XTIME(5,1,'- T1 equations                  ')
                   
           CALL HOINTM (FOO,EPS,T1,T2,BUF1,BUF2,BUF3,HO)
           
                   
           CALL HVINTM (FVV,EPS,T1,T2,BUF1,BUF2,BUF3,HV)
           
      IF (NDIMT1.GT.0) THEN
                      
          CALL HOVINTM (FVO,T1,BUF1,BUF2,HOV)
          CALL amplitude_equation_T1 (FVO,EPS,T1,T2,HO,HV,HOV,H,CZ,GOO,
     &                 GVV,BUF1,BUF2,NBUF2,BUF3,NBUF3,S1)
             
                                                            
                                                     
      ENDIF
      CALL XTIME(5,2,'- T1 equations                  ')

C
C     In case of NOCCD we can skip the calculation of T2 amplitudes
C
      IF (        NOCCD) GOTO 100
      CALL XTIME(5,1,'- T2 equations                  ')


                             
                     
                                                           
                                                           
                                                               
                                                             
           
C
C  EVALUATE T2 INTERMEDIATES
C
                   
         CALL XTIME(0,1,'-- GOINTM                        ')
         CALL GOINTM (FVO,T1,BUF1,BUF2,HO)
         CALL XTIME(0,2,'-- GOINTM                        ')
         CALL XTIME(0,1,'-- GVINTM                        ')
         CALL GVINTM (FVO,T1,BUF3,NBUF3,BUF2,NBUF2,HV)
         CALL XTIME(0,2,'-- GVINTM                        ')
         CALL XTIME(0,1,'-- AINTM                         ')
         CALL AINTM (T1,T2,BUF1,BUF2,AZ)
         CALL XTIME(0,2,'-- AINTM                         ')
         CALL XTIME(0,1,'-- HINTM                         ')
         CALL HINTM (T1,T2,BUF1,BUF2,BUF3,NBUF3,H)
         CALL XTIME(0,2,'-- HINTM                         ')

#if defined (VAR_MPI)
C *******************
C ** parallel case: the H-intermediate was partially calculated
C ** on the separate slaves. In order to perform the H*T2
C ** contraction we need the complete H array ==> we complete it
C ** here and synchronize all the nodes. They then have the complete
C ** H-array and this lengthy contraction can also be parallelized.
C
      IF (NMPROC .GT. 1) THEN

        IF(NBUF3.LT.IVOVO(NREP+1))
     &    CALL QUIT('BUF3 can not keep the VOVO integrals !')

C ... we collect the H array

        CALL XTIME(4,1,'-- combining HINTM via MPI_ALLREDUCE')
        CALL XCOPY (NBUF3,A0,0,BUF3,1)
CSK        WRITE(IW,*) 'THIS IS IERR before',IERR
        IERR = 0
        call interface_mpi_allreduce_r1_work_f77(H,BUF3,
     &       RCW*IVOVO(NREP+1),
     &       op_MPI_SUM,global_communicator)
        IF(IERR.GT.0) THEN
           CALL QUIT('MPI_REDUCE ERROR in collecting HINTM !')
        ENDIF
        CALL XCOPY(IVOVO(NREP+1),BUF3,1,H,1)
        CALL XTIME(4,2,'-- combining HINTM via MPI_ALLREDUCE')

      ENDIF

#endif


C
                     
                                                           
                                                           
                              
                                                               
                            
                                                        
           
C
C  INTERMEDIATES FORMED --- START T2 EQUATIONS
C
      CALL XTIME(5,1,'-- T2 EQNS                     ')
      CALL amplitude_equation_T2 (EPS,T1,T2,HO,HV,HOV,H,AZ,CZ,GOO,GVV,
     &             BUF1,BUF2,BUF3,NBUF3,S2)
      CALL XTIME(5,2,'-- T2 EQNS                     ')
C
C
C  MUST NOW ADD IN B INTERMEDIATE CONTRIBUTION
C
                     
                                                             
           
C
      CALL XTIME(5,1,'-- BINTM                       ')
      CALL BINTM (T1,T2,BUF1,BUF2,NBUF2,BUF3,NBUF3,S2,IW)
      CALL XTIME(5,2,'-- BINTM                       ')

                     
                                                                        
           
C
      CALL XTIME(5,2,'- T2 equations                  ')
C
 100  CONTINUE
C
C     Delete the Pi->Q excitations in case of Int. Ham. CCFS scheme 2
      IF (FOCKSP.AND.DOIH.AND.IHSCHEME.EQ.2) THEN
         CALL DENOMF (EPS,S1,S2,S1,S2,3)
      ENDIF
C
C     Add the first order wave function contribution
C
      CALL S1S2_ADDFO (FVO,BUF1,BUF2,S1,S2)
C
C  ADD PARTIAL S1 AND S2 (PARALLEL RUNS)
C
      CALL XTIME(4,1,'- adding partial T1/T2 amplitudes ')
      CALL PARTS (1,RMS,S1,S2,BUF1,BUF2)
      CALL XTIME(4,2,'- adding partial T1/T2 amplitudes ')
C
                     
                                                                       
                                                                       
           
C
      IF (FOCKSP) THEN
C       Adding folded diagrams, use HOV and H as scratch arrays
        CALL FOLDED (ICURFSS,IHEQ,T1,T2,BUF1,BUF2,BUF3,HOV,H,
     &               S1,S2)
                       
                                                                        
                                                                        
             
      ENDIF
C
C  Make new T1s and T2 by dividing S1s and S2 by eigenvalues
C
      CALL XCOPY (NDIMT1,S1,1,T1,1)
      CALL XCOPY (NDIMT2,S2,1,T2,1)
C
      IF (FOCKSP) THEN
        CALL DENOMF (EPS,T1,T2,T1,T2,1)
C       Remove the Hamiltonian elements, T1 and T2 should only contain
C       the wave function (amplitudes) appropriate for this sector.
        CALL DENOMF (EPS,T1,T2,T1,T2,2)
      ELSE
        CALL DENOM (EPS,T1,T2,T1,T2)
      ENDIF
C
C     Zero amplitudes if the NOCCX options are specified
C
      IF (        NOCCS) CALL XCOPY (NDIMT1,A0,0,T1,1)
      IF (        NOCCD) CALL XCOPY (NDIMT2,A0,0,T2,1)
C
      CALL ZCORE (T1,T2)
C
                                                          
                           
                                                           
                                                  
C
C  EXTRAPOLATE THE AMPLITUDES
C
      CALL XTIME(5,1,'- DIIS extrapolation            ')
      CALL DELFCK ('VO',DELS1,T1)
      CALL DELINT ('VVOO',DELS2,T2,0,0,0)
      CALL DIISX(NDIMT1,NDIMT2,BUF1,T1,BUF2,T2,
     &           IT,CC1,BB1,ITC,BB2,RMS)
      CALL XTIME(5,2,'- DIIS extrapolation            ')
C
C  SYNCHRONIZE T1 AND T2 (PARALLEL RUNS)
C
      CALL XTIME(4,1,'- synchronizing T1 & T2 amplitudes')
      CALL PARTS (2,RMS,T1,T2,BUF1,BUF2)
      CALL XTIME(4,2,'- synchronizing T1 & T2 amplitudes')
C
C  EVALUATE CORRELATION ENERGY AND WRITE OUT RESULTS
C
                   
         IF (IPRNT.GE.3) CALL PAIREN (T1,T2,FVO,BUF1,BUF2)
         CALL CCENG  (T1,T2,FVO,BUF1,BUF2,ECCSD,ECCSDIM)
           
C
      CALL T1DIAG (DIAG,T1,BUF1)
      IF (IPRNT.GE.1) WRITE(IW,6001) NIT,ECCSD,RMS,DIAG
      CALL FLSHFO(IW)
C
C
      IF(RMS.LT.CONV) GO TO 2000
      GO TO 1000
C
 2000 CONTINUE
      CALL T1DIAG (DIAG,T1,BUF1)
      WRITE(IW,6071)
C
C     Only print energy when doing amplitudes for the reference
C     wave function.
C
      IF (.NOT.FOCKSP) THEN
         PRENE =  .TRUE.
      ELSE
         PRENE = ICURFSS .EQ. 1
      ENDIF
C
      IF (PRENE) THEN
         ETOT = ESCF + ECCSD
         WRITE(IW,6072) ESCF
         IF (.NOT.(NOCCD.OR.NOCCS)) THEN
            WRITE(IW,6074) ECCSD
         ELSE
            IF (NOCCS) WRITE(IW,6084) ECCSD
            IF (NOCCD) WRITE(IW,6085) ECCSD
         ENDIF
         WRITE(IW,6076) ETOT
         IF (IPRNT.GE.2) CALL PAIREN (T1,T2,FVO,BUF1,BUF2)
      ENDIF
C
      WRITE(IW,6077) DIAG
      WRITE(IW,6078) RMS
      WRITE(IW,6079) NIT
      if (doxml) then
         call xml_quantity('CC iterations used',NIT)
         call xml_quantity('CC convergence',RMS,' ')
         call xml_quantity('T1 diagnostic',DIAG,' ')
      end if
C
      IF (FOCKSP) THEN
C        Calculate the energies for this sector
         CALL HEFF1(ICURFSS,IPRNT,S1,S2,FVO,EPS,BUF1,BUF2,BUF3,
     &              HEFF,EIGENVALUES,EIGENVECTORS)
      ENDIF
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
