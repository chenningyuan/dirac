      SUBROUTINE lambda_equation_mp2 (EPS,FVO,T1,S1,T2,S2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Solves MP2 amplitude or lambda equations.
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
C
      REAL*8 EPS(*),FVO(*),T1(*),S1(*),T2(*),S2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "results.inc"
#include "complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 DIAG
C
C---------------Executable code--------------------------------------
C
C  INITIALIZE T2 AND T1 AMPLITUDES
C
      CALL XCOPY (NDIMT1,FVO,1,T1,1)
      CALL GETVVOO (T2)
      CALL DENOM (EPS,T1,T2,T1,T2)
      CALL ZCORE (T1,T2)
C
C
C  INITIALIZE LAMDA AMPLITUDES
C
         CALL XCOPY (NDIMT1,T1,1,S1,1)
         CALL XCOPY (NDIMT2,T2,1,S2,1)
 6071 FORMAT(//,'  MP2 results',/)
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 7000 FORMAT(//' Timing of routine MP2EQN :'
     &/' Before MP2EQN :',T30,F12.3,' seconds'
     &/' Energy, T1-Diag :',T30,F12.3,' seconds'
     &/' Untimed parts :',T30,F12.3,' seconds'
     &/' Total time in MP2EQN :',T30,F12.3,' seconds')
      END
      SUBROUTINE lambda_equation_L2 (t1,t2,l1,l2,ho,hv,hov, 
     &             w_voov,w_oooo,G_oo,G_vv,buf1,buf2,
     &             buf3,nbuf3,resL2)

      use contraction
      use intermediates_1b_2b

      implicit none

!---------------Description--------------------------------------------
!     Evaluate residue for lambda equation.
!----------------------------------------------------------------------

!---intermediates required-------------
!f_bar_ae,f_bar_mi,f_bar_me,w_ejmb,w_ijmn,G_ea,G_im
!--------------------------------------

!---------------Common Blocks--------------------------------------
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "ccpar.inc"


!       calling variables
!       -----------------

      integer,intent(in) :: nbuf3
      real(8),intent(in) :: ho(rcw*nfoo),hv(rcw*nfvv),hov(nfvo*rcw)     ! one body intermediates
      real(8),intent(in) :: w_voov(nv4*rcw),w_oooo(nv1*rcw)! two body intermediates 
      real(8),intent(in) :: G_vv((ivv(nrep+1))*rcw),
     &                         G_oo((ioo(nrep+1))*rcw)              ! G-intermediates  
      real(8),intent(in) :: t1(ndimt1*rcw),t2(ndimt2*rcw)! converged T1 amplitude
      real(8),intent(in) :: L1(ndimt1*rcw),L2(ndimt2*rcw)! L1 and L2 amplitudes from last iteration.
      real(8),intent(in) :: buf1(*),buf2(*),buf3(*)      ! buffer arrays 
      real*8,intent(out) :: ResL2(ndimt2*rcw)          ! residue for L2 amplitude  


!     Local variables                                      
!     ---------------

      logical          :: done,teq
      real*8           :: t1val(2)
      real*8,allocatable           :: buf4(:)
      complex*16       :: zdotc
      real*8           :: ddot,val
      real*8,allocatable :: int_vv(:)
      real*8,allocatable :: int_oo(:),w_vvvv(:)
      real*8,allocatable :: L2_temp(:),ResL2_temp(:),tau(:)
      integer     :: jklrep,m2,irep,jrep,abrep,jkloff,ij,ijkl,ijkl1,mint
      integer     :: off2,i,j,istart,m,k,n,irp,ntot,nint 
      integer     :: jabrep,x(nrep,nrep),y(nrep) 
      

!---------------Executable code--------------------------------------

    5   format(10f10.6)  

       ResL2 = 0.0d0 

       call getoovv (ResL2) 

!-------------------------------------------------
!   term1 : ResL2 = ResL2 + P(ab) L2 (ij,ae)*f_bar_ae(e,b) 
!-------------------------------------------------

       call contraction_424 ((/"o1","o2","p1","p3"/),(/"p3","p2"/),
     &    (/"o1","o2","p1","p2"/),ResL2,1.0d0,0.0d0,nrep,L2,hv) 

       write(*,*)'RESL2 after term1',ddot(ndimt2,resL2,1,resL2,1)
!----------------------------------------------
!   term2 : ResL2 = ResL2 - P(ij) L2 (im,ab)*f_bar_im (j,m)   
!----------------------------------------------

       call contraction_424 ((/"o1","o3","p1","p2"/),(/"o2","o3"/),
     &    (/"o1","o2","p1","p2"/),ResL2,-1.0d0,0.0d0,nrep,L2,ho) 

       write(*,*)'RESL2 after term2',ddot(ndimt2,resL2,1,resL2,1)
!--------------------------------------------------------
!   term3 : ResL2 = ResL2 + 0.5 L2(mn,ab) * W_ijmn(ij,mn)
!--------------------------------------------------------

      call contraction_444((/"o3","o4","p1","p2"/),
     &    (/"o1","o2","o3","o4"/),(/"o1","o2","p1","p2"/),W_oooo,ResL2,
     &    0.5d0,0.0d0,nrep,L2)

       write(*,*)'RESL2 after term3',ddot(ndimt2,resL2,1,resL2,1)
!--------------------------------------------------------
!   term4 : ResL2 = ResL2 + P(ij)P(ab) L2(im,ae) * W_ejmb(ej,mb)
!--------------------------------------------------------

      call contraction_444((/"o1","o3","p1","p3"/),
     &    (/"p3","o2","o3","p2"/), (/"o1","o2","p1","p2"/),W_voov,ResL2,
     &    1.0d0,0.0d0,nrep,L2)

       write(*,*)'RESL2 after term4',ddot(ndimt2,resL2,1,resL2,1)
!-----------------------------------------------------------
!   term5 : ResL2 = ResL2 + P(ab) V(ij,ae)*(G_vv(e,b)-L1(m,b)*t1(e,m))
!-----------------------------------------------------------

      allocate(int_vv((ivv(nrep+1))*rcw))

      int_vv = 0.0d0

      call xcopy ((ivv(nrep+1))*rcw,G_vv,1,int_vv,1)

      call contraction_222((/"o3","p2"/),(/"p3","o3"/),(/"p3","p2"/),L1,
     &        t1,int_vv,-1.0d0,1.0d0,nrep)

      call contraction_424 ((/"o1","o2","p1","p3"/),(/"p3","p2"/),
     &        (/"o1","o2","p1","p2"/),ResL2,1.0d0, 
     &        0.0d0,nrep,RightTensor=int_vv)

      deallocate(int_vv) 

       write(*,*)'RESL2 after term5',ddot(ndimt2,resL2,1,resL2,1)
!--------------------------------------------------------------------
!    term6 : ResL2 = ResL2 + P(ij) V(im,ab)*(G_oo(j,m)+L1(j,e)*t(e,m)) 
!--------------------------------------------------------------------

      allocate(int_oo((ioo(nrep+1))*rcw))

      int_oo = 0.0d0

      call xcopy ((ioo(nrep+1))*rcw,G_oo,1,int_oo,1)

      call contraction_222((/"o2","p3"/),(/"p3","o3"/),(/"o2","o3"/),L1,
     &        t1,int_oo,1.0d0,1.0d0,nrep)

      call contraction_424 ((/"o1","o3","p1","p2"/),(/"o2","o3"/),
     &        (/"o1","o2","p1","p2"/),ResL2,1.0d0, 
     &        0.0d0,nrep,RightTensor=int_oo)

      deallocate(int_oo) 

!-----------------------------------------------------------------------
!   term 7 : ResL2 = ResL2 - L1(m,a)*V(ij,mb)
!-----------------------------------------------------------------------

      call contraction_244 ((/"o3","p1"/),(/"o1","o2","o3","p2"/),
     &     (/"o1","o2","p1","p2"/),ResL2,-1.0d0,1.0d0,nrep,L1)

       write(*,*)'RESL2 after term7',ddot(ndimt2,resL2,1,resL2,1)
!-------------------------------------------------------------------------
!   term 8 : ResL2 = ResL2 + P(ij)P(ab) L1(i,a)*f_bar_me(j,b) 
!-------------------------------------------------------------------------

      call lambda2_disconnected(l1,hov,ResL2)

!     out-of-core terms

      call contraction_244 ((/"o1","p3"/),(/"p3","o2","p1","p2"/),
     &     (/"o1","o2","p1","p2"/),ResL2,1.0d0,1.0d0,nrep,L1)

      write(*,*)'RESL2 after term8',ddot(ndimt2,resL2,1,resL2,1)

!--------------------------------------------------------
!    term 10 : ResL2 = ResL2 + L2(ij,ef) * W_vvvv(ef,ab) 
!--------------------------------------------------------

      done = .false.
      allocate(tau(ndimt2*rcw))
      call gettau(t1,t2,tau)
      off2 = 1
      do irp = 1, nrep
         istart = 0
         if (nvvt(irp).eq.0) cycle
 10      mint = nbuf3/nvvt(irp)         
         if (idist(1,1,irp).gt.istart) istart=idist(1,1,irp)
         ntot = idist(2,1,irp)
         nint = min0(mint,ntot-istart)
 
         allocate(w_vvvv(nvvt(irp)*nint))

!        call getvvvv (irp,istart,nint,done,w_vvvv,mint)
         call w_efab_out_of_core (irp,istart,done,tau,buf1,nbuf3,w_vvvv)
         n = nvvt(irp)
         m = noot(irp)
         k = nint
         call xgemm ('n','n',m,k,n,AP5,l2,m,
     &           w_vvvv,n,a1,ResL2(off2+istart*m*rcw),m)

         deallocate(w_vvvv) 
         if (.not.done) then
            istart = istart + nint
            goto 10
         endif
         off2 = off2 + m * n * rcw
       enddo

       deallocate(tau)

!     call xaxpy (ndimt2,a1,ResL2_temp,1,ResL2,1)

       write(*,*)'RESL2 after term10',ddot(ndimt2,resL2,1,resL2,1)

!-----------------------------------------------------------------------
!     term9 : ResL2 = ResL2 + P(ij) L1(i,e) * V(ej,ab)  
!-----------------------------------------------------------------------
       
!     done = .false.
!     allocate (L2_temp(loovvt(nrep+1)*rcw))
!     off2 = 1
!     do abrep = 1, nrep
!        if (nov(abrep).eq.0.or.nvvt(abrep).eq.0) cycle
!        mint = nbuf3/nvvt(abrep)
!        istart = 0
!110     call getvovv (abrep,istart,nint,done,buf1,mint)

!!   sort V(ej,ab) -> V'(e,jab)
!         ij = 1
!         ijkl = 0
!         do jrep = 1, nrep
!           
!          irep = multb(jrep,abrep+nrep,2)           

!          jklrep = irep

!           call symm_index_outofcore_sorter(istart,x)
!           jkloff = (x(jrep,abrep)*nv(irep)) * rcw

!           write(*,*)'jkloff',jkloff,nint,nv(irep) 
!           jkloff = 0
!          m2 = nv(irep) * no(jrep)
!          do j = 1, no(jrep)
!           do i = 1, nv(irep)
!            ijkl1 = jkloff + ((j-1) * nv(irep) + i - 1) * rcw + 1
!               call xcopy(nint,buf1(ijkl+ij),nvo(abrep),buf2(ijkl1),m2)
!               ij = ij + rcw
!           enddo
!          enddo


!!   do the contraction: L1(i,e)*V'(e,jab) -> L2'(i,jab)

!        call xgemm ('n','n',no(irep),nint*no(jrep),nv(irep),a1,L1,
!    &                no(irep),buf2,nv(irep),a1,
!    &                L2_temp(off2+istart*m2*rcw),no(irep)) 

!         enddo
!        if (.not.done) then
!        istart = istart + nint
!        goto 110
!        endif
!        off2 = off2 + no(irep)*novvt(jklrep)*rcw ! doubtful 
!     enddo

!!   sorting L2'(i,jab) -> L2(i>j,a>b) 

!     allocate(ResL2_temp(ndimt2)) 

!     call srt1t2 (nrep,multb,ltr,noot,no,no,nvvt,loovvt,llovvt,
!    &             L2_temp,ResL2_temp)

!     deallocate(L2_temp)

!     call xaxpy (ndimt2,a1,ResL2_temp,1,ResL2,1)

!     deallocate(ResL2_temp)

!     write(*,*)'RESL2 after term9',ddot(ndimt2,resL2,1,resL2,1)


!     term9 : ResL2 = ResL2 + P(ij) L1(i,e) * V(ej,ab)  


      END

      subroutine lambda2_disconnected(l1,f_bar,ResL2)
 
!---------------Description--------------------------------------------

!     Calculating the contribution from the disconnected piece in lambda2 residue

!---------------Common Blocks----------------------------------------

      implicit none
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------

      real*8, intent(in)    :: l1(ndimt1*rcw),f_bar(nfvo*rcw) 
      real*8, intent(inout) :: ResL2(ndimt2*rcw)
 
!---------------Local variables--------------------------------------

      integer :: ijaboff,abrep,arep,brep,a,b,amin,imin,i,j,ij,jb,ia,ja,
     &           ib,ni,irep,jrep  
      real*8  :: T1MIN(2)
!---------------Executable code--------------------------------------

      t1min=(/0.0d0,0.0d0/) 

      ijaboff = 0
      DO ABREP = 1, NREP
      DO BREP = 1, NREP
      AREP = MULTB(BREP,ABREP+NREP,2)
      IF (AREP.LT.BREP) CYCLE
      DO B = 1, NV(BREP)
         AMIN = 1
         IF (AREP.EQ.BREP) AMIN = B + 1
         DO A = AMIN, NV(AREP)
            DO JREP = 1, NREP
            IREP = MULTB(JREP,ABREP+NREP,2)
            IF (IREP.LT.JREP) CYCLE
!--------------------------------------------------------------
! Resl2(ij,ab) = Resl2(ij,ab) + l1(i,a) * f_bar(j,b) + f_bar(i,a) * l1(j,b)
!--------------------------------------------------------------
            IF (JREP.EQ.BREP) THEN
               IJ = 1
               JB = (IVO(BREP) + (B-1) * NO(BREP)) * RCW + 1
               DO J = 1, NO(JREP)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = J + 1
                  IA = (IVO(AREP) + (A-1) * NO(AREP) + IMIN - 1)*RCW + 1
                  NI = NO(IREP) - IMIN + 1
                  CALL XAXPY (NI,f_bar(jb),L1(ia),1,RESL2(IJABOFF+IJ),1)
                  CALL XAXPY (NI,L1(JB),f_bar(IA),1,RESL2(IJABOFF+IJ),1)
                  IJ = IJ + NI * RCW
                  JB = JB + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------------
C Resl2(ij,ab) = Resl2(ij,ab) - L1(j,a) * f_bar(i,b) - f_bar(j,a) * L1(i,b)
C--------------------------------------------------------------
            IF (JREP.EQ.AREP) THEN
               IJ = 1
               JA = (IVO(AREP) + (A-1) * NO(AREP)) * RCW + 1
               DO J = 1, NO(JREP)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = J + 1
                  IB = (IVO(BREP) + (B-1) * NO(BREP) + IMIN - 1)*RCW + 1
                  NI = NO(IREP) - IMIN + 1
                  T1MIN(1) = -f_bar(IB)
                  IF (CARITH) T1MIN(2) = -f_bar(IB+1)
                  CALL XAXPY (NI,T1MIN,L1(JA),1,ResL2(IJABOFF+IJ),1)
                  T1MIN(1) = -L1(IB)
                  IF (CARITH) T1MIN(2) = -L1(IB+1)
                  CALL XAXPY (NI,T1MIN,f_bar(JA),1,ResL2(IJABOFF+IJ),1)
                  IJ = IJ + NI * RCW
                  JA = JA + RCW
               ENDDO
            ENDIF
C--------------------------------------------------------
C UPDATE OFFSET AND GO TO NEXT IRREP PAIR
C--------------------------------------------------------
            IF (IREP.NE.JREP) THEN
               IJABOFF = IJABOFF + NO(IREP) * NO(JREP) * RCW
            ELSE
               IJABOFF = IJABOFF + NO(IREP) * (NO(IREP)-1) * RCW / 2
            ENDIF
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

C---------------------------------------------
C S(AB,IJ) = S(AB,IJ) + TAU(AB,KL) * AZ(KL,IJ)
C---------------------------------------------
!     CALL XTIME(0,1,'--- T2EQNS: TAU*AINTM contraction')
!        CALL XCOPY (NDIMT2,T2,1,BUF1,1)
!     OFF3 = 1
!     DO IREP = 1, NREP
!       M    = NVVT(IREP)
!       N    = NOOT(IREP)
!       KLD  = NOOT(IREP)
!       K    = IDIST(1,5,IREP)
!       OFF1 = IDIST(2,5,IREP)
!       OFF2 = IDIST(3,5,IREP)
!       CALL XGEMM ('N','N',M,N,K,A1,BUF1(OFF1),M,AZ(OFF2),KLD,
!    &              A1,S2(OFF3),M)
!       OFF3 = OFF3 + M * N * RCW
!     ENDDO
! ***
! *** original CNTRCT call in serial version
! ***
! *** CALL CNTRCT ('N','N',NVVT,NOOT,NOOT,A1,BUF1,AZ,A1,S2,NREP)
! ***
!     CALL XTIME(0,2,'--- T2EQNS: TAU*AINTM contraction')
! ** we will not parallelize a summation over just one occupied index
! ** if the system is small this can lead to zero batches. Only master
! ** does the T2(AB,IK) * GO(K,J) contraction
! ********************
! ** TASK SPLITTING **
! ********************
!     IF(MYPROC.EQ.MASTER) THEN
!----------------------------------------------------------------
! S(AB,IJ) = S(AB,IJ) - PIJ [ T2(AB,IK) * GO(K,J) ]
!----------------------------------------------------------------
!     CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
!    &             T2,BUF1)
!     CALL CNTRCT ('N','N',NVVOT,NO,NO,A1,BUF1,GO,A0,BUF2,NREP)
!     CALL SRT1T3 (NREP,MULTB,LTR,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
!    &             BUF2,BUF1)
!     CALL XAXPY (NDIMT2,-A1,BUF1,1,S2,1)
! ***************************
! ** END OF TASK SPLITTING **
! ***************************
!     ENDIF
!-----------------------------------------------------------------------
! S(AB,IJ) = S(AB,IJ) + PIJ [ W(AB,CJ) * T1(C,I) ]
!-----------------------------------------------------------------------
!     CALL XTIME(0,1,'--- T2EQNS: VOVV*T1              ')
!     DO 30 IRP = 1, NREP
!        IF (NVO(IRP).EQ.0) GOTO 30
!        DONE = .FALSE.
!        ISTART = 0
!        MINT = NBUF3/NVO(IRP)  ! max. # of NVO batches that fit in mem.
! 10     CONTINUE
!           CALL GETDZ   (IRP,ISTART,NINT,DONE,BUF3)
!        IF (CARITH) CALL CONJUGA (NVO(IRP)*NINT,BUF3,1)
!        IJ = 0
!        DO 20 JREP = 1, NREP
!          CREP = MULTB(JREP,IRP+NREP,2)
!          IREP = CREP
!          IF (IREP.LT.JREP) GOTO 20
!          DO J = 1, NO(JREP)
!             IMIN = 1
!             IF (JREP.EQ.IREP) IMIN = J + 1
!             DO I = IMIN, NO(IREP)
!                IJ = IJ + 1
!                ABIJ = (IVVOOTT(IRP)+(IJ-1)*NVVT(IRP)+ISTART)*RCW+1
!                DO C = 1, NV(CREP)
!                   CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
!                   T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
!                   T1VAL(1) = T1(T1CI)
!                   IF (CARITH) T1VAL(2) = T1(T1CI+1)
!                   CALL XAXPY(NINT,T1VAL,BUF3(CJ),NVO(IRP),S2(ABIJ),1)
!                ENDDO
!             ENDDO
!          ENDDO
!20      CONTINUE
!        JI = 0
!        DO 21 IREP = 1, NREP
!          JREP = MULTB(IREP,IRP+NREP,2)
!          CREP = IREP
!          IF (JREP.LT.IREP) GOTO 21
!          DO I = 1, NO(IREP)
!             JMIN = 1
!             IF (JREP.EQ.IREP) JMIN = I + 1
!             DO J = JMIN, NO(JREP)
!                JI = JI + 1
!                ABJI = (IVVOOTT(IRP)+(JI-1)*NVVT(IRP)+ISTART)*RCW+1
!                DO C = 1, NV(CREP)
!                   CJ = (IIVO(CREP,JREP)+(J-1)*NV(CREP)+C-1) * RCW + 1
!                   T1CI = (IVO(CREP)+(I-1)*NV(CREP)+C-1) * RCW + 1
!                   T1VAL(1) = - T1(T1CI)
!                   IF (CARITH) T1VAL(2) = - T1(T1CI+1)
!                   CALL XAXPY(NINT,T1VAL,BUF3(CJ),NVO(IRP),S2(ABJI),1)
!                ENDDO
!             ENDDO
!          ENDDO
!21      CONTINUE
!        IF (.NOT.DONE) THEN
!           ISTART = ISTART + NINT
!           GOTO 10
!        ENDIF
!30    CONTINUE
!     CALL XTIME(0,2,'--- T2EQNS: VOVV*T1              ')
! ****************************
! ** END OF TASK SPLITTING. **
! ****************************
! ** this contribution coming from the intermediate G(A,C) has to be
! ** calculated on each node. The fact that G(A,C) contains LOCAL +
! ** DISTRIBUTED parts as well is accounted for by the original calculation
! ** of the G(A,C) array before entering T2EQNS
!-------------------------------------------------
! S(AB,IJ) = S(AB,IJ) + PAB [ GV(A,C) * T2(CB,IJ)]
!-------------------------------------------------
!     CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,T2,BUF2)
!     CALL CNTRCT ('N','N',NV,NVOOT,NV,A1,GV,BUF2,A0,BUF1,NREP)
! ** the contraction of G(A,C) with the T2 amplitudes now is in BUF1 !
! ** these are the partial contributions available on EACH node
! ***
! ****** task splitting
! ***
!     IF(MYPROC.EQ.MASTER) THEN
!        CALL XCOPY (NV2,CZ,1,BUF2,1)
!     CALL SRT1L1 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,IVOOOT,IIOV,BUF2,BUF3)
!     CALL SRT1S2 (NREP,MULTB,LFA,NVO,NO,NV,NOOT,LOVOOT,LLVOOT,
!    &             BUF3,BUF2)
!     CALL CNTRCT ('N','N',NV,NVOOT,NO,-A1,T1,BUF2,A1,BUF1,NREP)
! ***
! *** end of task splitting
! ***
!     ENDIF
!     CALL SRT1T2 (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
!    &             BUF1,BUF2)
!     CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)
! ** now the final H*T2 contraction. H is complete since we communicated
! ** it already.
!-----------------------------------------------------------------------
! S(AB,IJ) = S(AB,IJ) + PAB.PIJ [H(AI,CK)*T2(CK,BJ)]
!-----------------------------------------------------------------------
!     CALL XTIME(0,1,'--- T2EQNS: HINTM*T2           ')
!     CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
!    &              T2,BUF1)
!     CALL XCOPY (J2VOVO(NREP+1),A0,0,BUF2,1)
!     OFF3 = 1
!     DO IREP = 1, NREP
!        JREP = MULTB(IREP+NREP,1+NREP,2)
!        M    = MVO(IREP)
!        N    = MVO(JREP)
!        KLD  = MVO(IREP)
!        K    = IDIST(1,10,JREP)
!        OFF1 = IDIST(2,10,JREP)
!        OFF2 = IDIST(3,10,JREP)
!        CALL XGEMM ('N','N',M,N,K,A1,H(OFF1),M,BUF1(OFF2),KLD,
!    &                A0,BUF2(OFF3),M)
!        OFF3 = OFF3 + M * N * RCW
!     ENDDO
!     CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
!    &              BUF2,BUF1)
!     CALL XAXPY (NDIMT2,A1,BUF1,1,S2,1)
!     CALL XTIME(0,2,'--- T2EQNS: HINTM*T2           ')
!==========================================================
!============= L A M B D A   E Q U A T I O N S ============
!==========================================================
!-------------------------------------------------------
! S(AB,IJ) = S(AB,IJ) + TAU(AB,KL) * V(KL,CD) * L(CD,IJ)
!-------------------------------------------------------
!     CALL GETAMPT (BUF1,BUF2)
!     CALL GETTAU (BUF1,BUF2,BUF3)
!     CALL GETVVOO (BUF1)
!     CALL CNTRCT ('C','N',NOOT,NOOT,NVVT,A1,BUF1,T2,A0,BUF2,NREP)
!     CALL CNTRCT ('N','N',NOOT,NOOT,NVVT,A1,BUF3,BUF2,A1,S2,NREP)
!-----------------------------------------------
! S(AB,IJ) = S(AB,IJ) + PAB.PIJ. [T(A,I)*H(B,J)]
!-----------------------------------------------
!     CALL LAMBD21 (T1,HOV,S2)
!---------------------------------------------
! S(AB,IJ) = S(AB,IJ) + PAB. [G(A,C)*W(CB,IJ)]
!---------------------------------------------
!     CALL GETVVOO (BUF1)
!     CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
!    &             BUF1,BUF2)
!     CALL CNTRCT ('N','N',NV,NVOOT,NV,A1,GVV,BUF2,A0,BUF1,NREP)
!     CALL SRT1T2 (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
!    &             BUF1,BUF2)
!     CALL XAXPY (NDIMT2,A1,BUF2,1,S2,1)
!---------------------------------------------
! S(AB,IJ) = S(AB,IJ) - PIJ. [W(AB,IK)*G(K,J)]
!---------------------------------------------
!     CALL GETVVOO(BUF2)
!     CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
!    &             BUF2,BUF1)
!     CALL CNTRCT ('N','N',NVVOT,NO,NO,A1,BUF1,GOO,A0,BUF2,NREP)
!     CALL SRT1T3 (NREP,MULTB,LTR,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
!    &             BUF2,BUF1)
!     CALL XAXPY (NDIMT2,-A1,BUF1,1,S2,1)
!-----------------------------------------------------------------------
! S(AB,IJ) = S(AB,IJ) + PAB.PIJ [L(D,I)*V(AK,CD)*T2(BC,KJ)]
!-----------------------------------------------------------------------
!     CALL GETAMPT (BUF1,BUF3)
!     CALL DINTM2 (T1,BUF1,BUF3,NBUF3,BUF2)
!     CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
!    &              BUF3,BUF1)
!     CALL XCOPY (J2VOVO(NREP+1),A0,0,BUF3,1)
!     OFF1 = 1
!     OFF2 = 1
!     OFF3 = 1
!     DO IRP =1 , NREP
!        JRP = MULTB(IRP+NREP,1+NREP,2)
!        M = MVO(IRP)
!        N = MVO(JRP)
!        K = MVO(IRP)
!        CALL XGEMM ('N','N',M,N,K,A1,BUF2(OFF1),M,BUF1(OFF2),K,
!    &               A0,BUF3(OFF3),M)
!        OFF1 = OFF1 + M * K * RCW
!        OFF2 = OFF2 + K * N * RCW
!        OFF3 = OFF3 + M * N * RCW
!     ENDDO
!     CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,J2VOVO,JJVO,JJVO,
!    &              BUF3,BUF1)
!     CALL XAXPY (NDIMT2,A1,BUF1,1,S2,1)
!

      END


      subroutine symm_index_outofcore_sorter(istart,x)

#include "symm.inc"
!------------------------------------------------
        integer, intent(in) :: istart
        integer, intent(out) :: x(nrep,nrep)
!------------------------------------------------

!------------------------------------------------
        integer :: y(nrep), krep, ijrep, ijkrep
!------------------------------------------------

            y = 0

          do ijrep = 1,nrep 
            do krep = 1,nrep

            ijkrep = multb(krep,ijrep+nrep,2) 
            x(krep,ijrep) = y(ijkrep)
            y(ijkrep) =  y(ijkrep) + (istart+1)*no(krep)

            enddo
         enddo

      end
