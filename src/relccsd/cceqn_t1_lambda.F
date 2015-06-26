      SUBROUTINE lambda_equation_L1(t1,t2,l1,l2,ho,hv,hov,w_voov,w_ovoo,
     &                               w_vvvo,G_oo,G_vv,buf1,buf2,buf3,
     &                               nbuf2,nbuf3,ResL1,iprnt)

       use contraction

       implicit none

!       Description                                            
!       -----------
!     Evaluate       Lambda1 equations
!     Here we have to take care of the fact that certain integral classes
!     are now split over individual nodes. The corresponding contributions
!     which are linear in the integrals are added at the end of the T1
!     routines. Therefore we have to take care of calculating the
!     non-distributed contributions to the T1 amplitudes only on the
!     master node in order to avoid double-counting in the PARTS routine !
!     all direct VOVV/VVVV contributions and VOVV/VVVV contributions stemming
!     from intermediates have to be calculated on EACH node.

!     intermediates required for L1 equation             
!     ---------------------------------------
!     f_bar_mi,f_bar_me,f_bar_ae,w_iemn,w_ejmb,G_ea,G_im        


!   common blocks
!   -------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
#include "ccpar.inc"

!       calling variables
!       -----------------

      integer,intent(in) :: nbuf2,nbuf3,iprnt
      real(8),intent(in) :: ho(nfoo*rcw),hv(nfvv*rcw),hov(nfvo*rcw)    ! one body intermediates
      real(8),intent(in) :: w_voov(ivoov(nrep+1)*rcw),
     &                      w_ovoo(iovoot(nrep+1)*rcw),w_vvvo(nv5*rcw) ! two body intermediates 
      real(8),intent(in) :: G_vv((ivv(nrep+1))*rcw),
     &                         G_oo((ioo(nrep+1))*rcw)                 ! G-intermediates  
      real(8),intent(in) :: t1(ndimt1*rcw),t2(ndimt2*rcw)              ! converged T1 amplitude
      real(8),intent(in) :: L1(ndimt1*rcw),L2(ndimt2*rcw)              ! L1 and L2 amplitudes from last iteration.
      real(8),intent(in) :: buf1(*),buf2(*),buf3(*)                    ! buffer arrays 
      real*8,intent(inout) :: ResL1(ndimt1*rcw)                        ! residue for L1 amplitude  


!     Local variables                                      
!     ---------------
      logical            :: done,teq,usedz,right
      real*8,allocatable :: int_temp(:), l2_sorted1(:)
      real*8,allocatable :: w_v_temp(:), w_o_temp(:)
      real*8             :: ddot 
      integer            :: i,j,ij
!-------------------
! Initialize residue
!-------------------

    5   format(10f10.6)  

!      call xcopy (ndimt1,hov,1,ResL1,1)


       if (iprnt.ge.1)
     &   write(*,*)'RESL1: init',ddot(ndimt1,resL1,1,resL1,1)
!---------------------------------
!    ResL1 = ResL1 + L(i,e) * f_bar_ae(e,a)
!---------------------------------

       call contraction_222 ((/"o1","p2"/),(/"p2","p1"/),(/"o1","p1"/),
     &    L1,hv,ResL1,1.0d0,1.0d0,nrep)


       if (iprnt.ge.1)
     & write(*,*)'RESL1: dia1',ddot(ndimt1,resL1,1,resL1,1)
!----------------------------------------
!    ResL1 = ResL1 - L(m,a)*f_bar_mi(i,m)
!----------------------------------------

        call contraction_222 ((/"o2","p1"/),(/"o1","o2"/),(/"o1","p1"/),
     &    L1,ho,ResL1,-1.0d0,1.0d0,nrep)

       if (iprnt.ge.1)
     & write(*,*)'RESL1: dia2',ddot(ndimt1,resL1,1,resL1,1)
!----------------------------------------------
!      ResL1 = ResL1 - G(m,n) * V(mi,na)  
!----------------------------------------------

! second way:

!   w_mnie(mn,ie) = w_mnie(mn,ie) + V(mn,fe)*t(f,i)

      allocate(w_o_temp(iovoot(nrep+1))) 

      call getooov(w_o_temp)

      call contraction_424 ((/"o2","o3","p1","p2"/),(/"p2","o1"/),
     &        (/"o2","o3","o1","p1"/),w_o_temp,-1.0d0, 
     &        1.0d0,nrep,RightTensor=t1)

      call contraction_242 ((/"o2","o3"/),(/"o2","o1","o3","p1"/), 
     & (/"o1","p1"/),G_oo,ResL1,-1.0d0,1.0d0,nrep,RightTensor=w_o_temp)

!     call contraction_242 ((/"o2","o3"/),(/"o2","o1","o3","p1"/), 
!    & (/"o1","p1"/),G_oo,ResL1,-1.0d0,1.0d0,nrep)

      deallocate(w_o_temp) 


       if (iprnt.ge.1)
     & write(*,*)'RESL1: dia3',ddot(ndimt1,resL1,1,resL1,1)
!-------------------------------------------------
!      ResL1 = ResL1 + L1(m,e) * W(ei,ma) 
!second      ResL1 = ResL1 + L1(m,e) * W(ie,am) 
!-------------------------------------------------

!     call contraction_242 ((/"o2","p2"/),(/"p2","o1","o2","p1"/),
!    &   (/"o1","p1"/),L1,ResL1,1.0d0,1.0d0,nrep,RightTensor=w_voov)

      call contraction_242 ((/"o2","p2"/),(/"o1","p2","p1","o2"/),
     &   (/"o1","p1"/),L1,ResL1,1.0d0,1.0d0,nrep,RightTensor=w_voov)

       if (iprnt.ge.1)
     & write(*,*)'RESL1 after term4',ddot(ndimt1,resL1,1,resL1,1)
!-----------------------------------------------------
!       ResL1 = ResL1 - 0.5d0* L2(mn,ae)*W(ie,mn) 
!-----------------------------------------------------

      call contraction_442((/"o2","o3","p1","p2"/),(/"o1","p2","o2","o3"
     &          /),(/"o1","p1"/),ResL1,-1.0d0,1.0d0,nrep,L2,w_ovoo)


       if (iprnt.ge.1)
     & write(*,*)'RESL1 after term5',ddot(ndimt1,resL1,1,resL1,1)

!-------------------------------------------------------
!      ResL1    = ResL1 + L2(im,ef) * W(ef,am) 
!-------------------------------------------------------

      allocate (l2_sorted1(ndimt2*rcw))

       l2_sorted1 = 0.0d0
      call srt1c1(nrep,noot,nvvt,L2,l2_sorted1)

      call srt1t3 (nrep,multb,.false.,nvvt,no,no,nvvot,kvvoot,kkvvot,
     &           l2_sorted1,buf2)

      deallocate (l2_sorted1)

      call srt19 (nrep,multb,.false.,nvvt,nv,no,nvvot,kvvovt,kkvvot,
     &               w_vvvo,buf3)

      call cntrct ('c','n',no,nv,nvvot,-a1,buf2,buf3,a1,
     &    ResL1,nrep)


!     call contraction_442((/"o1","o2","p2","p3"/),(/"p2","p3","p1","o2"
!    &          /),(/"o1","p1"/),ResL1,1.0d0,1.0d0,nrep,LeftTensor=L2)!,
!    &           RightTensor=w_vvvo)


       if (iprnt.ge.1)
     & write(*,*)'RESL1 after term6',ddot(ndimt1,resL1,1,resL1,1)

!-----------------------------------------------------------
!      int_temp(e,m) = (G(e,f) * t(f,m) - G(n,m)*t(e,n)) ! corrected from stanton and Gauss 
!      ResL1    = ResL1 + int_temp(e,m) * V(im,ae) 
!------------------------------------------------------------

      allocate(int_temp(ndimt1*rcw))

      int_temp = 0.0d0

!     call contraction_222 ((/"p2","p3"/),(/"p3","o2"/),(/"p2","o2"/),
!    &    G_vv,t1,int_temp,1.0d0,1.0d0,nrep)

!     call contraction_222 ((/"o3","o2"/),(/"p2","o3"/),(/"p2","o2"/),
!    &    G_oo,t1,int_temp,-1.0d0,1.0d0,nrep)

!------------------------------------------------------------

!     call contraction_242 ((/"p2","o2"/),(/"o1","o2","p1","p2"/),
!    &   (/"o1","p1"/),int_temp,ResL1,1.0d0,1.0d0,nrep)

      deallocate(int_temp)


       if (iprnt.ge.1)
     & write(*,*)'RESL1 after term7',ddot(ndimt1,resL1,1,resL1,1)

!---------------------------------------------------------------------------

!second way.

! w_amef(am,ef) = w_amef(am,ef) - V(nm,ef) * t(a,n) 

       allocate (w_v_temp(nv5*rcw))

        call getvovv_incore(w_v_temp)
       
      call contraction_424 ((/"o2","o1","p3","p1"/),(/"p2","o2"/),
     &        (/"p2","o1","p3","p1"/),w_v_temp,-1.0d0, 
     &        1.0d0,nrep,RightTensor=t1)

!-------------------------------------------------------------------------------



!-------------------------------------------------------
!  out-of-core term....
!    ResL1 = ResL1 + G(f,e) * V(ei,fa)
!-------------------------------------------------------

      call contraction_242 ((/"p2","p3"/),(/"p2","o1","p3","p1"/),
     &  (/"o1","p1"/),G_vv,ResL1,-1.0d0,1.0d0,nrep,RightTensor=w_v_temp)

!     call contraction_242 ((/"p3","p2"/),(/"p2","o1","p3","p1"/),
!    &  (/"o1","p1"/),G_vv,ResL1,-1.0d0,1.0d0,nrep)
       deallocate(w_v_temp)


       if (iprnt.ge.1)
     & write(*,*)'RESL1 after term8',ddot(ndimt1,resL1,1,resL1,1)
!     m = nfvv
!     n = nfvo
!     if (nbuf3.lt.jvvvo(2)) then
!        call quit(" out of core step in gv not implemented")
!     endif
!     usedz = .true.
!     right = .false.
!lv : changed to srt20d, untested as the rest of the lambda code
!lv : addition should be supplied with a minus sign
!     call srt20d (nrep,multb,nvo,nv,no,nv,nv,mvv,jvvvo,jjvv,
!    &             jjvo,buf2,nbuf2,ResL1,g_vv,usedz,right)


!        CALL QUIT('force quit. just before the old codes')

C!********************
C!** TASK SPLITTING **
C!********************

!!    IF (MYPROC.EQ.MASTER) THEN

!!    CALL XTIME(0,1,'--- T1EQNS   HV*T / T*HO              ')
C!----------------------------------------------
C!S(A,I) = S(A,I) + HV(A,C) * T(C,I)
C!----------------------------------------------
!!    CALL CNTRCT ('N','N',NV,NO,NV,A1,HV,L1,A1,resL1,NREP)
C!----------------------------------------------
C!S(A,I) = S(A,I) - T(A,K) * HO(K,I)
C!----------------------------------------------
!!    CALL CNTRCT ('N','N',NV,NO,NO,-A1,L1,HO,A1,resL1,NREP)
!!    CALL XTIME(0,2,'--- T1EQNS   HV*T / T*HO            ')

C!-------------------------------------------
C!S(A,I) = S(A,I) + V(K,L,C,I) * TAU(A,C,K,L)
C!-------------------------------------------
!!    CALL XTIME(0,1,'--- T1EQNS VOOO*TAU                ')
!!                 
!!        
!!       CALL XCOPY (NDIMT2,T2,1,BUF3,1)
!!         
!!    CALL T1EQN1 (BUF1,BUF2,BUF3,NBUF3,resL1)
!!    CALL XTIME(0,2,'--- T1EQNS VOOO*TAU               ')

C!****************************************************
C!** END OF TASK SPLITTING. THE COMMON STUFF STARTS **
C!****************************************************

!!    ENDIF

C!-------------------------------------------
C!TAU'(CD,IK) = L(CD,IK) +  G(C,B) * T(BD,IK)   <<<<<<------- LAMBDA EQNS
C!-------------------------------------------
!!       CALL GETAMPT (BUF1,BUF2)
!!       CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
!!   &                BUF2,BUF1)
!!       CALL CNTRCT ('N','N',NV,NVOOT,NV,A1,G_VV,BUF1,A0,BUF2,NREP)
!!       CALL SRT1T2 (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
!!   &                BUF2,BUF1)
!!       CALL XAXPY (NDIMT2,A1,T2,1,BUF1,1)
!!         
!!    CALL XTIME(0,1,'--- T1EQNS VOVV contribution  ')
!!    CALL T1EQN2 (BUF1,BUF2,BUF3,NBUF3,resL1)   ! VOVV contractions, term 6 !
!!    CALL XTIME(0,2,'--- T1EQNS VOVV contribution  ')

C!-------------------------------------
C!S(A,I) = S(A,I) - W(A,K,C,I) * T(C,K)    ! VOVO contributions
C!-------------------------------------
!!    CALL XTIME(0,1,'--- T1EQNS VOVO * T(C,K)         ')
!!        
!!       CALL XCOPY (JVOVO(2),H,1,BUF2,1)

!!    M    = MVO(1)
!!    K    = IDIST(1,3,1)
!!    OFF1 = IDIST(2,3,1)
!!    OFF2 = IDIST(3,3,1)
!!    CALL XGEMV ('N',M,K,-A1,BUF2(OFF1),M,T1(OFF2),1,A1,resL1,1)
!!    CALL XTIME(0,2,'--- T1EQNS VOVO * T(C,K)        ')
C!-----------------------------------
C!Terms specific for Lambda1 equations
C!-----------------------------------
!!                      
C!---------------------------------------------------
C!S(A,I) = S(A,I) - T(A,L) * V(K,L,C,D) * L(C,D,K,I)
C!       = S(A,I) - T(A,L) * V*(C,D,K,L) * L(C,D,K,I)
C!---------------------------------------------------
!!    CALL GETVVOO(BUF2)
!!    CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
!!   &             BUF2,BUF1)
!!    CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
!!   &             T2,BUF2)
!!    CALL CNTRCT ('C','N',NO,NO,NVVOT,A1,BUF1,BUF2,A0,BUF3,NREP)
!!    CALL GETAMPT (BUF1,BUF2)
!!    CALL CNTRCT ('N','N',NV,NO,NO,-A1,BUF1,BUF3,A1,resL1,NREP)
C!---------------------------------------------------
C!S(A,I) = S(A,I) - L(A,C,K,L) * V(K,L,D,C) * T(D,I)
C!       = S(A,I) - L(A,C,K,L) * V*(D,C,K,L) * T(D,I)
C!---------------------------------------------------
!!    CALL GETVVOO(BUF2)
!!    CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
!!   &             BUF2,BUF1)
!!    CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
!!   &             T2,BUF2)
!!    CALL CNTRCT ('N','C',NV,NV,NVOOT,A1,BUF2,BUF1,A0,BUF3,NREP)
!!    CALL GETAMPT (BUF1,BUF2)
!!    OFF1 = 1
!!    OFF2 = 1
!!    OFF3 = 1
!!    DO IRP =1 , NREP
!!       M = NV(IRP)
!!       N = NO(IRP)
!!       K = NV(IRP)
!!       CALL XGEMM ('N','N',M,N,K,-A1,BUF3(OFF1),M,BUF1(OFF2),K,
!!   &               A1,resL1(OFF3),M)
!!       OFF1 = OFF1 + M * K * RCW
!!       OFF2 = OFF2 + K * N * RCW
!!       OFF3 = OFF3 + M * N * RCW
!!    ENDDO
C!-------------------------------------
C!S(A,I) = S(A,I) - W(C,A,D,I) * G(C,D)
C!-------------------------------------
!!    M = NFVV
!!    N = NFVO
!!    IF (NBUF3.LT.JVVVO(2)) THEN
!!       CALL QUIT(" OUT OF CORE STEP IN GV NOT IMPLEMENTED")
!!    ENDIF
!!    USEDZ = .TRUE.
!!    RIGHT = .FALSE.
c!v : changed to SRT20D, untested as the rest of the lambda code
c!v : addition should be supplied with a minus sign
!!    CALL SRT20D (NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
!!   &             JJVO,BUF2,NBUF2,resL1,G_VV,USEDZ,RIGHT)
C!-------------------------------------
C!S(A,I) = S(A,I) + W(A,K,L,I) * G(K,L)
C!-------------------------------------
!!    CALL SRT1ST4 (NREP,MULTB,LTR,NV,NO,NO,NO,MVO,JVOOO,JJVO,JJOO,
!!   &              CZ,BUF2)
!!    M = NFVO
!!    K = NFOO
!!    CALL XGEMV ('N',M,K,A1,BUF2,M,G_OO,1,A1,resL1,1)
C!--------------------------------------------------
C!S(A,I) = S(A,I) + G(C,B) * W(B,K,C,D) * T(A,D,K,I)
C!--------------------------------------------------
!!    USEDZ = .FALSE.
!!    RIGHT = .FALSE.
!!    CALL SRT20D (NREP,MULTB,NVO,NV,NO,NV,NV,MVV,JVVVO,JJVV,
!!   &            JJVO,BUF2,NBUF2,BUF1,G_VV,USEDZ,RIGHT)
!!    CALL GETAMPT (BUF3,BUF2)
!!    CALL SRT1TT4 (NREP,MULTB,LTR,LFA,NV,NV,NO,NO,
!!   &              MVO,JVOVO,JJVO,JJVO,BUF2,BUF3)
!!    M = MVO(1)
!!    K = MVO(1)
!!    CALL XGEMV ('N',M,K,A1,BUF3,M,BUF1,1,A1,resL1,1)
           
      RETURN
      END
