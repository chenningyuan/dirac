      SUBROUTINE cceqn_driver_lambda(MAXIT,EPS,FOO,FVO,FVV,L1,S1,L2,S2,
     &                 W_OOOO,CZ,W_VOOV,W_OVVO,W_OVOO,HO,HV,HOV,GOO,GVV,
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


      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*),L1(*),S1(*),L2(*),S2(*)
      REAL*8 W_OOOO(*),W_VOOV(*),W_OVVO(*),W_OVOO(*),HO(*),HV(*),
     & HOV(*),BUF1(*),BUF2(*),BUF3(*)
      REAL*8 CZ(*) ! is only allocated for Lambda equations
      REAL*8 GOO(*),GVV(*)
      INTEGER NBUF2,NBUF3 

C     for diis ...

      REAL*8 CC1(*)
      REAL*8 BB1(*),BB2(*)

      REAL*8 RMS,CONV,RMS1,ddot
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

      REAL*8     :: CPU,CPU4V,CPUH1,CPUH2
      LOGICAL    :: FOCKSP,PRENE
      real*8,allocatable :: resL1(:),resL2(:)
      real*8,allocatable :: G_oo(:),G_vv(:)
!      integer            :: i !remove
      real*8             :: val  

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

      IF (EQNS.EQ.'LAMBDA') THEN
         FOCKSP = .FALSE.
      ELSE
         CALL QUIT('Wrong keyword in call to CCEQNS')
      ENDIF


C     Initialize the work arrays to zero.

      CALL XCOPY (NDIMT2,A0,0,BUF1,1)
      CALL XCOPY (NBUF2,A0,0,BUF2,1)
      CALL XCOPY (NBUF3,A0,0,BUF3,1)

!     Calculate fixed intermediaries for Lambda equations

    5   format(10f10.6)  
                        
         call getampt (s1,s2)

         call f_bar_mi(foo,fvo,s1,s2,ho)

       write(*,*)'whh done',ddot(nfoo,ho,1,ho,1)

        write(*,*)'whh done'

        write(77,*)'f_bar_mi intermediate'
        write(77,*)'---------------------'
!       write(*,5)(ho(i),i=1,nfoo)

         call w_ijmn(S1,S2,w_oooo)

!       write(*,*)'whhhh done'

        write(*,*)'whhhh done',ddot(nv1,w_oooo,1,w_oooo,1)
        write(78,*)'w_ijmn intermediate'
        write(78,*)'---------------------'
        write(78,5)(w_oooo(i),i=1,nv1)

         call f_bar_me(fvo,S1,hov)

         write(*,*)'f_me intermediate'

         write(*,*)'whp done',ddot(nfvo,hov,1,hov,1)

         call f_bar_ae(fvv,fvo,s1,s2,hv)

         write(*,*)'wpp done',ddot(nfvv,hv,1,hv,1)
         write(79,*)'f_bar_ae intermediate'
         write(79,*)'---------------------'
         write(79,5)(hv(i),i=1,nfvv)

        call w_ejmb(S1,S2,w_voov)

        write(*,*)'w_ejmb done',ddot(nv4,w_voov,1,w_voov,1)
        write(80,*)'w_bar_ejmb intermediate'
        write(80,*)'---------------------'
        write(80,5)(w_voov(i),i=1,nv4)

         call wbar_mbej(S2,w_ovvo)

        write(*,*)'w_mbej done',ddot(nv4,w_ovvo,1,w_ovvo,1)
        write(81,*)'w_bar_mbej intermediate'
        write(81,*)'---------------------'
        write(81,5)(w_ovvo(i),i=1,nv4)
      
!       write(*,*)'starting w_iemn intermediate'

         call w_iemn(s1,s2,w_ovoo,w_ovvo,w_oooo,hov)

         write(*,*)'w_iemn done',ddot(nv2,w_ovoo,1,w_ovoo,1)
         write(*,*)'w_bar_iemn intermediate'
         write(*,*)'---------------------'
         write(82,5)(w_ovoo(i),i=1,nv2)

         call w_efam(s1,s2,w_ovvo,hov,cz)

         call w_efam_diagram2 (s1,s2,buf1,buf2,buf3,nbuf3,cz)
         write(*,*)'w_bar_efam intermediate'
         write(*,*)'---------------------'
         write(83,5)(cz(i),i=1,nv5)

! out-of-core section of the intemediates

!       call f_bar_ae_out_of_core (s1,buf2,nbuf2,hv)

!       call w_ejmb_out_of_core (s1,s2,buf1,buf2,buf3,nbuf3,w_voov)

!      write(*,*)'w_ejmb done',ddot(nv4,w_voov,1,w_voov,1)
!      write(80,*)'w_bar_ejmb intermediate'
!      write(80,*)'---------------------'
!      write(80,5)(w_voov(i),i=1,nv4)

!        call w_iemn_out_of_core (s1,s2,buf1,buf2,nbuf3,w_ovoo)

!        write(*,*)'w_bar_iemn intermediate'
!        write(*,*)'---------------------'
!        write(82,5)(w_ovoo(i),i=1,nv2)

!       call w_efam_out_of_core (s1,s2,hov,w_voov,w_ovvo,buf1,buf2,
!    &              buf3,nbuf3,nbuf2,CZ)

!         write(*,*)'w_bar_efam intermediate'
!         write(*,*)'---------------------'
!         write(83,5)(CZ(i),i=1,nv5)


!-------Initialization Step for the Left Eigenvectors-----------

      call srt1c1 (nrep,nv,no,s1,l1) 
      call srt1c1 (nrep,nvvt,noot,s2,l2) 

!     Initialize the work arrays to zero.

      call xcopy (ndimt2,A0,0,buf1,1)
      call xcopy (nbuf2,A0,0,buf2,1)
      call xcopy (nbuf3,A0,0,buf3,1)

!---------before starting the iterative procedure, we need to initialize the DIIS-------------
      rms = 0.0d0
      itc = 0
      it = -1
      call diisx(ndimt1,ndimt2,L1,L2,it,cc1,bb1,itc,bb2,rms)

      nit = 0

      if (iprnt.ge.1) then
         write(iw,6003)
         write(iw,6001) nit,rms
      endif

!--------------------------------------
!  Start Iterative Procedure
!--------------------------------------

       allocate(resL1(ndimt1*rcw)) 
       allocate(resL2(ndimt2*rcw)) 

      do nit = 1, maxit 

!     Making non-sequential G_ea, G_im intermediates for lambda equations

!      ResL1 = 0.0d0
!      ResL2 = 0.0d0

      allocate (G_vv(ivv(nrep+1)*rcw))
      allocate (G_oo(ioo(nrep+1)*rcw))

      call G_intm (s2,l2,G_vv,G_oo)

!     evaluate L1 equation  

      CALL XTIME(5,1,'- L1 equations                  ')
           
      if (ndimt1 > 0) then

       ResL1 = 0.0d0

       CALL lambda_equation_L1(s1,s2,l1,l2,ho,hv,hov,w_voov,w_ovoo,
     &                         cz,G_oo,G_vv,buf1,buf2,buf3,
     &                         nbuf2,nbuf3,resL1)

      endif

      CALL XTIME(5,2,'- L1 equations                  ')
                                 
      CALL XTIME(5,1,'- L2 equations                  ')

C  INTERMEDIATES FORMED --- START T2 EQUATIONS
C
      CALL XTIME(5,1,'-- L2 EQNS                     ')

       CALL lambda_equation_L2(s1,s2,l1,l2,ho,hv,hov,w_voov,w_oooo,
     &                            G_oo,G_vv,buf1,buf2,buf3,
     &                            nbuf3,resL2)

       write(*,*)'RESL2 afetr all diagrams',ddot(ndimt2,resL2,1,resL2,1)

      CALL XTIME(5,2,'-- T2 EQNS                     ')

      CALL XTIME(5,2,'- L2 equations                  ')

C     Read out timers for the most important contractions

      CALL XTIMRD (CPU,'-- BINTM                       ')
      CPU4V = CPU - CPU4V
      CALL XTIMRD (CPU,'--- HINTM: VVOO contribution   ')
      CPUH1 = CPU - CPUH1
      CALL XTIMRD (CPU,'--- T2EQNS: HINTM*T2           ')
      CPUH2 = CPU - CPUH2

       deallocate (G_oo)
       deallocate (G_vv)

C     Compute & print the performance in Mega/Giga/TeraFlops

      IF (TIMING) CALL BFLOPS (CPU4V,CPUH1,CPUH2,NIT)

!------------------------------------
!     calculation of RMS value 
!------------------------------------

!  Make new L1s and L2s by dividing ResL1s and ResL2 by eigenvalues


!      rms = ddot(ndimt1,resL1,1,resL1,1)+ddot(ndimt2,resL2,1,resL2,1)
!      rms = sqrt((rms)/(ndimt1+ndimt2))

       rms = 0.0d0
       do i =1,ndimt1
         rms = rms + resl1(i)*resl1(i)
       enddo

       do i =1,ndimt2
         rms = rms + resl2(i)*resl2(i)
       enddo


       call denom (eps,resL1,resL2,resL1,resL2)

!       write(*,*)'dot product of RESL1',ddot(ndimt1,resL1,1,resL1,1)
!       write(*,*)'dot product of RESL2',ddot(ndimt2,resL2,1,resL2,1)

!------------------------------------------------------------------------------------------------------------
!       here i have summed the residue and old L amplitude value with a negative sign because the denominatior
!       has opposite sign w.r.t. the T amplitudes, which has been used in the denom routine 
!-------------------------------------------------------------------------------------------------------------
      call xaxpy (ndimt1,-A1,resL1,1,l1,1)
      call xaxpy (ndimt2,-A1,resL2,1,l2,1)


      CALL XTIME(5,1,'- DIIS extrapolation            ')
      call diisx(ndimt1,ndimt2,l1,l2,   
     &           it,cc1,bb1,itc,bb2,rms1)
      CALL XTIME(5,2,'- DIIS extrapolation            ')

!  SYNCHRONIZE T1 AND T2 (PARALLEL RUNS)

      CALL XTIME(4,1,'- synchronizing T1 & T2 amplitudes')
      CALL PARTS (2,RMS,L1,L2,BUF1,BUF2)
      CALL XTIME(4,2,'- synchronizing T1 & T2 amplitudes')


      IF (IPRNT.GE.1) WRITE(IW,6001) NIT,RMS
      CALL FLSHFO(IW)

      if(rms.lt.conv) exit

      enddo  ! end of iteration loop

      deallocate(resL1) 
      deallocate(resL2) 

!     Close and delete DIIS scratch files

      IT = -2
      call diisx(ndimt1,ndimt2,L1,L2,
     &           it,cc1,bb1,itc,bb2,rms)
c
 111  FORMAT(A20,I4,2f25.15)
 6001 FORMAT(2X,I3,5X,2F20.15,F10.5)
 6002 FORMAT(3X,'MP2   =',F20.15/3X,'EMP2  =',F20.15/)
 6003 FORMAT(//3X,'NIT',10X,'RMS')
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
