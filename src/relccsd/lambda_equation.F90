  module lambda_equation

!this module consists of subroutines of lambda-equations at different level of theory.

  use contraction
  use intermediates_1b_2b

    public lambda_equation_mp2 
    public lambda_equation_L2 
    public lambda_equation_L1 


     private

      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface


     contains


      subroutine lambda_equation_mp2 (eps,fvo,t1,s1,t2,s2)
 
      IMPLICIT INTEGER (A-Z)
 
!---------------Description--------------------------------------------
!
!     Solves MP2 amplitude or lambda equations.
!
!---------------Routines called----------------------------------------
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Calling variables--------------------------------------

!     Intermediaries etc.

      REAL*8 EPS(*),FVO(*),T1(*),S1(*),T2(*),S2(*)

!---------------Common Blocks--------------------------------------

#include "files.inc"
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "results.inc"
#include "complex.inc"

!---------------Local variables--------------------------------------

      REAL*8 DIAG

!---------------Executable code--------------------------------------

!  INITIALIZE T2 AND T1 AMPLITUDES

      CALL XCOPY (NDIMT1,FVO,1,T1,1)
      CALL GETVVOO (T2)
      CALL DENOM (EPS,T1,T2,T1,T2)
      CALL ZCORE (T1,T2)


!  INITIALIZE LAMDA AMPLITUDES

         CALL XCOPY (NDIMT1,T1,1,S1,1)
         CALL XCOPY (NDIMT2,T2,1,S2,1)
 6071 FORMAT(//,'  MP2 results',/)                      
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 7000 FORMAT(//' Timing of routine MP2EQN :'      &      
     &/' Before MP2EQN :',T30,F12.3,' seconds'    &
     &/' Energy, T1-Diag :',T30,F12.3,' seconds'  &
     &/' Untimed parts :',T30,F12.3,' seconds'    &
     &/' Total time in MP2EQN :',T30,F12.3,' seconds')
      end subroutine lambda_equation_mp2


      subroutine lambda_equation_L2 (t1,t2,l1,l2,ho,hv,hov, &
     &             w_voov,w_oooo,G_oo,G_vv,buf1,buf2, &
     &             buf3,nbuf3,resL2,iprnt,niter)

      use contraction

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
#include "files.inc"

!       calling variables
!       -----------------

      integer,intent(in) :: nbuf3,iprnt,niter
      real(8),intent(in) :: ho(rcw*nfoo),hv(rcw*nfvv),hov(nfvo*rcw)     ! one body intermediates
      real(8),intent(in) :: w_voov(ivoov(nrep+1)*rcw),w_oooo(nv1*rcw)! two body intermediates 
      real(8),intent(in) :: G_vv(nfvv*rcw), &
     &                         G_oo(nfoo*rcw)            ! G-intermediates  
      real(8),intent(in) :: t1(ndimt1*rcw),t2(ndimt2*rcw)! converged T1 amplitude
      real(8),intent(in) :: L1(ndimt1*rcw),L2(ndimt2*rcw)! L1 and L2 amplitudes from last iteration.
      real(8),intent(in) :: buf1(*),buf2(*),buf3(*)      ! buffer arrays 
      real*8,intent(out) :: ResL2(*)                     ! residue for L2 amplitude  


!     Local variables                                      
!     ---------------

      logical          :: done,teq
      real*8           :: t1val(2)
      complex*16       :: zdotc
      real*8           :: ddot,val
      real*8,allocatable :: int_vv(:)
      real*8,allocatable :: int_oo(:),w_vvvv(:),w_intm(:),int_wooop(:)
      real*8,allocatable :: L2_temp(:),ResL2_temp(:),tau(:)
      integer     :: jklrep,m2,irep,jrep,abrep,jkloff,ij,ijkl,ijkl1,mint
      integer     :: off2,i,j,istart,m,k,n,irp,ntot,nint,kstart,off1 
      integer     :: jabrep,x(nrep,nrep),y(nrep) 
      
      interface
         subroutine gettau(t1,t2,tau)
         real*8,intent(in) :: t1(*),t2(*)
         real*8,intent(out) :: tau(*)
         end subroutine gettau
      end interface

!---------------Executable code--------------------------------------

    5   format(10f10.6)  

       call getoovv (ResL2) 

       if (iprnt.ge.1) &
     & write(*,*)'RESL2 after initialization',(ResL2(i),i=1,5)
!-------------------------------------------------
!   term1 : ResL2 = ResL2 + P(ab) L2 (ij,ae)*f_bar_ae(e,b) 
!-------------------------------------------------

       call contraction_424 ((/"o1","o2","p1","p3"/),(/"p3","p2"/), &
     &    (/"o1","o2","p1","p2"/),ResL2,1.0d0,1.0d0,nrep,LeftTensor=L2, &
     &    RightTensor=hv) 

       if (iprnt.ge.1) &

     & write(*,*)'RESL2 after term1',(ResL2(i),i=1,5)
!----------------------------------------------
!   term2 : ResL2 = ResL2 - P(ij) L2 (im,ab)*f_bar_im (j,m)   
!----------------------------------------------

       call contraction_424 ((/"o1","o3","p1","p2"/),(/"o2","o3"/), &
     &    (/"o1","o2","p1","p2"/),ResL2,-1.0d0,1.0d0,nrep, &
     &    LeftTensor=L2,RightTensor=ho) 

       if (iprnt.ge.1) &

     & write(*,*)'RESL2 after term2',(ResL2(i),i=1,5)
!--------------------------------------------------------
!   term3 : ResL2 = ResL2 + 0.5 L2(mn,ab) * W_ijmn(ij,mn)
!--------------------------------------------------------

      call contraction_444((/"o3","o4","p1","p2"/), &
     &    (/"o1","o2","o3","o4"/),(/"o1","o2","p1","p2"/),ResL2, &
     &    1.0d0,0.0d0,nrep,RightTensor=w_oooo,LeftTensor=L2)

       if (iprnt.ge.1) &
     & write(*,*)'RESL2 after term3',(ResL2(i),i=1,5)
!--------------------------------------------------------
!   term4 : ResL2 = ResL2 + P(ij)P(ab) L2(im,ae) * W_ejmb(ej,mb)
!   second : ResL2 = ResL2 + P(ij)P(ab) L2(im,ae) * W_bmej(je,bm)

!     allocate(w_intm(JOVOV(NREP+1)*rcw))
!     w_intm = 0.0d0 

!     CALL SRT1TT5 (NREP,MULTB,.FALSE.,.FALSE.,NO,NO,NV,NV,MOV,J2OVVO, 
!    & JJOV,JJVO,L2,BUF2)

!     CALL SRT1SS4 (NREP,MULTB,.FALSE.,.FALSE.,NO,NV,NV,NO,MOV,J2OVVO,
!    &  JJOV,JJVO,w_voov,BUF3)

!     CALL CNTRCT ('N','C',MOV,MOV,MVO,A1,BUF2,BUF3,A0,w_intm,NREP)
!     do i = 1, JOVOV(NREP+1)

!      if ((abs(w_intm(i)) <= 1.0d-14)) then

!       w_intm(i) = 0.0d0

!      endif

!      enddo


!     CALL SRT1TT4 (NREP,MULTB,.FALSE.,.TRUE.,NO,NO,NV,NV,MOV,JOVOV,
!    &    JJOV,JJOV,w_intm,BUF2)

!     CALL XAXPY (ndimt2,A1,BUF2,1,ResL2,1)

!     
!     deallocate(w_intm)

!     call contraction_444((/"o1","o3","p1","p3"/),
!    &    (/"p3","o2","o3","p2"/), (/"o1","o2","p1","p2"/),ResL2,
!    &    1.0d0,0.0d0,nrep,righttensor=w_voov,LeftTensor=L2)

      call contraction_444((/"o1","o3","p1","p3"/), &
     &    (/"o2","p3","p2","o3"/), (/"o1","o2","p1","p2"/),ResL2, &
     &    1.0d0,0.0d0,nrep,righttensor=w_voov,LeftTensor=L2)

       if (iprnt.ge.1) &
     & write(*,*)'RESL2 after term4',(ResL2(i),i=1,5)
!-----------------------------------------------------------
!   term5 : ResL2 = ResL2 + P(ab) V(ij,ae)*(G_vv(b,e)-L1(m,b)*t1(e,m))
!-----------------------------------------------------------

      call contraction_424 ((/"o1","o2","p1","p3"/),(/"p2","p3"/), &
     &        (/"o1","o2","p1","p2"/),ResL2,1.0d0,1.0d0,nrep,RightTensor=G_vv)

       if (iprnt.ge.1) &
     & write(*,*)'RESL2 after term5',(ResL2(i),i=1,5)

!      allocate(int_vv((ivv(nrep+1))*rcw))
!     allocate(int_vv(nfvv*rcw))

!      int_vv = 0.0d0

!     call xcopy ((ivv(nrep+1))*rcw,G_vv,1,int_vv,1)

!     call contraction_222((/"o3","p1"/),(/"p3","o3"/),(/"p3","p1"/),L1,
!    &        t1,int_vv,-1.0d0,1.0d0,nrep)

!     call contraction_424 ((/"o1","o2","p3","p2"/),(/"p3","p1"/),
!    &        (/"o1","o2","p1","p2"/),ResL2,1.0d0, 
!    &        1.0d0,nrep,RightTensor=int_vv)

!     deallocate(int_vv) 

       if (iprnt.ge.1) &
     & write(*,*)'RESL2 after term6',(ResL2(i),i=1,5)

!----------------------------------------------------------------------------
!     w(mn,ie) = V(mn,ie) + t(f,i) * V(mn,fe)
!----------------------------------------------------------------------------

      allocate (int_wooop(iovoot(nrep+1)*rcw))
     
       int_wooop = 0.0d0

      call getooov(int_wooop)

      call contraction_424 ((/"o1","o2","p2","p3"/),(/"p3","o3"/), &
     &        (/"o1","o2","o3","p2"/),int_wooop,-1.0d0, &
     &        1.0d0,nrep,RightTensor=t1)

!---------------------------------------------------------------------------
!     L2(ij,ab) = L2(ij,ab) + L1(m,a) * W(ij,mb) 
!---------------------------------------------------------------------------
      call contraction_244 ((/"o3","p1"/),(/"o1","o2","o3","p2"/), &
     &   (/"o1","o2","p2","p1"/),ResL2,-1.0d0,-1.0d0,nrep,LeftTensor=L1, &
     &   RightTensor=int_wooop)

      deallocate (int_wooop)

       if (iprnt.ge.1) &

     & write(*,*)'RESL2 after term7',(ResL2(i),i=1,5)
!     & write(*,*)'RESL2 after term5',ddot(ndimt2,resL2,1,resL2,1)
!--------------------------------------------------------------------
!    term6 : ResL2 = ResL2 - P(ij) V(im,ab)*(G_oo(j,m)+L1(j,e)*t(e,m)) 
!--------------------------------------------------------------------

      call contraction_424 ((/"o1","o3","p1","p2"/),(/"o3","o2"/), &
     &        (/"o1","o2","p1","p2"/),ResL2,-1.0d0, &
     &        1.0d0,nrep,RightTensor=G_oo)    ! be careful, recheck.

    if (iprnt.ge.1) &

    & write(*,*)'RESL2 after term8',(ResL2(i),i=1,5)

      allocate(int_oo(nfoo*rcw))

      int_oo = 0.0d0

!      call xcopy ((ioo(nrep+1))*rcw,G_oo,1,int_oo,1)

      call contraction_222((/"o1","p3"/),(/"p3","o3"/),(/"o1","o3"/),L1, &
     &        t1,int_oo,1.0d0,1.0d0,nrep)

      call contraction_424 ((/"o3","o2","p1","p2"/),(/"o1","o3"/), &
     &        (/"o1","o2","p1","p2"/),ResL2,-1.0d0, &
     &        1.0d0,nrep,RightTensor=int_oo)

      deallocate(int_oo) 

     call contraction_244 ((/"o1","p3"/),(/"p3","o2","p1","p2"/), &
    &  (/"o1","o2","p1","p2"/),ResL2,1.0d0,1.0d0,nrep,LeftTensor=L1)

       if (iprnt.ge.1) &

     & write(*,*)'RESL2 after term9',(ResL2(i),i=1,5)
!    & write(*,*)'RESL2 after term6',ddot(ndimt2,resL2,1,resL2,1)

!-----------------------------------------------------------------------
!   term 7 : ResL2 = ResL2 - L1(m,a)*V(ij,mb)
!-----------------------------------------------------------------------

!!    call contraction_244 ((/"o3","p1"/),(/"o1","o2","o3","p2"/),
!!   &   (/"o1","o2","p2","p1"/),ResL2,-1.0d0,-1.0d0,nrep,LeftTensor=L1)

!!     if (iprnt.ge.1)

!!   & write(*,*)'RESL2 after term7',(ResL2(i),i=1,5)
!    & write(*,*)'RESL2 after term7',ddot(ndimt2,resL2,1,resL2,1)
!-------------------------------------------------------------------------
!   term 8 : ResL2 = ResL2 + P(ij)P(ab) L1(i,a)*f_bar_me(j,b) 
!-------------------------------------------------------------------------

      call lambda2_disconnected(l1,hov,ResL2)

       if (iprnt.ge.1) &

     & write(*,*)'RESL2 after term10',(ResL2(i),i=1,5)
!     write(*,*)'RESL2 after term8',ddot(ndimt2,resL2,1,resL2,1)

!     out-of-core terms
!-----------------------------------------------------------------
!     term9 : ResL2 = ResL2 + P(ij) L1(i,e) * V(ej,ab)  
!-----------------------------------------------------------------
!     call contraction_244 ((/"o1","p3"/),(/"p3","o2","p1","p2"/), &
!    &     (/"o1","o2","p1","p2"/),ResL2,1.0d0,1.0d0,nrep,LeftTensor=L1)

!!     if (iprnt.ge.1)

!!   & write(*,*)'RESL2 after term9',(ResL2(i),i=1,5)
!    & write(*,*)'RESL2 after term9',ddot(ndimt2,resL2,1,resL2,1)

!--------------------------------------------------------
!    term 10 : ResL2 = ResL2 + L2(ij,ef) * W_vvvv(ef,ab) 
!--------------------------------------------------------

      
!        if (niter<=1) then 
!        call waio_open(itapt+7)
!        endif

      done = .false.
      allocate(tau(ndimt2*rcw))
      call gettau(t1,t2,tau)
      off2 = 1
      off1 = 1
      do irp = 1, nrep
         istart = 0
         if (nvvt(irp).eq.0) cycle
         mint = nbuf3/nvvt(irp)         
 10      kstart = istart
         if (idist(1,1,irp).gt.kstart) kstart=idist(1,1,irp)
         ntot = idist(2,1,irp)
         nint = min0(mint,ntot-kstart)
 
         allocate(w_vvvv(nvvt(irp)*nint*rcw))

!!         if (niter<=1) then 
!!         call w_efab_out_of_core (irp,istart,done,tau,buf1,t1,nbuf3, &
!!     &           w_vvvv,nint)

!!         call putdz (irp,istart,nint,w_vvvv)

!!         else
         
         call getdz (irp,istart,nint,done,w_vvvv,nbuf3)

!!         endif

         n = nvvt(irp)
         m = noot(irp)
         k = nint
         call xgemm ('n','n',m,k,n,A1,l2(off1),m, &
     &           w_vvvv,n,a1,ResL2(off2+istart*m*rcw),m)

         deallocate(w_vvvv) 
         if (.not.done) then
            istart = istart + nint
            goto 10
         endif
         off2 = off2 + m * n * rcw
         off1 = off1 + m * n * rcw
       enddo

       deallocate(tau)

       if (iprnt.ge.1) &

     & write(*,*)'RESL2 after term11',(ResL2(i),i=1,5)

!--------------------------------------------------------------------------------
!     missing term from w_vvvv
!     w_intm (ij,mf) = w_intm(ij,mf) + L2(ij,ef) * t(e,m)
!     ResL2(ij,ab) = ResL2(ij,ab) + w_intm(ij,mf) * v(mf,ab) 
!--------------------------------------------------------------------------------

!!    allocate(w_intm(iovoot(nrep+1)*rcw))

!!    w_intm = 0.0d0
!!    call contraction_424 ((/"o1","o2","p3","p4"/),(/"p3","o3"/),
!!   &        (/"o1","o2","o3","p4"/),w_intm,1.0d0, 
!!   &        1.0d0,nrep,LeftTensor=L2,RightTensor=t1)


!!!    call srt1t3 (nrep,multb,.false.,noot,nv,nv,noovt,koovvt,kkoovt, &
!!!  &             L2,buf2)

!!!    call xcopy(koovot(nrep+1),a0,0,buf3,1)

!!!    call cntrct ('n','n',noovt,no,nv,-a1,buf2,t1,a0,buf3,nrep)


!!!   call srt19 (nrep,multb,.true.,noot,no,nv,noovt,koovot,kkoovt,buf3, &
!!!  &  buf2)

!!!   call xaxpy (iovoot(nrep+1),a1,buf2,1,w_intm,1)

!!!   call contraction_444((/"o1","o2","o3","p4"/), &
!!!  &    (/"o3","p4","p1","p2"/), (/"o1","o2","p1","p2"/),ResL2, &
!!!  &    -1.0d0,1.0d0,nrep,LeftTensor=w_intm)

!!!    deallocate(w_intm)      

!!!    if (iprnt.ge.1) &

!!!  & write(*,*)'RESL2 after term12',(ResL2(i),i=1,5)

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

      end subroutine

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

      integer :: ijaboff,abrep,arep,brep,a,b,amin,imin,i,j,ij,jb,ia,ja, &
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
!--------------------------------------------------------------
! Resl2(ij,ab) = Resl2(ij,ab) - L1(j,a) * f_bar(i,b) - f_bar(j,a) * L1(i,b)
!--------------------------------------------------------------
            IF (JREP.EQ.AREP) THEN
               IJ = 1
               JA = (IVO(AREP) + (A-1) * NO(AREP)) * RCW + 1
               DO J = 1, NO(JREP)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = J + 1
                  IB = (IVO(BREP) + (B-1) * NO(BREP) + IMIN - 1)*RCW + 1
                  NI = NO(IREP) - IMIN + 1
                   T1MIN(1) = -f_bar(JA)
                   IF (CARITH) T1MIN(2) = -f_bar(JA+1)
                  CALL XAXPY (NI,T1MIN,L1(IB),1,ResL2(IJABOFF+IJ),1)
                   T1MIN(1) = -L1(JA)
                   IF (CARITH) T1MIN(2) = -L1(JA+1)
                  CALL XAXPY (NI,T1MIN,f_bar(IB),1,ResL2(IJABOFF+IJ),1)
                  IJ = IJ + NI * RCW
                  JA = JA + RCW
               ENDDO
            ENDIF
!--------------------------------------------------------
! UPDATE OFFSET AND GO TO NEXT IRREP PAIR
!--------------------------------------------------------
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

      END SUBROUTINE


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

      end subroutine


      subroutine lambda_equation_L1(t1,t2,l1,l2,ho,hv,hov,w_voov,w_ovoo, &
     &                               w_vvvo,G_oo,G_vv,buf1,buf2,buf3, &
     &                               nbuf2,nbuf3,ResL1,iprnt)

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
      real(8),intent(in) :: w_voov(ivoov(nrep+1)*rcw), &
     &                      w_ovoo(iovoot(nrep+1)*rcw),w_vvvo(nv5*rcw) ! two body intermediates 
      real(8),intent(in) :: G_vv(nfvv*rcw),  &
     &                         G_oo(nfoo*rcw)                 ! G-intermediates  
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

      call xcopy (ndimt1,hov,1,ResL1,1)

       if (iprnt.ge.1) &
     &   write(*,*)'RESL1: init',dot_product(resL1,resL1)
!---------------------------------
!    ResL1 = ResL1 + L(i,e) * f_bar_ae(e,a)
!---------------------------------

       call contraction_222 ((/"o1","p2"/),(/"p2","p1"/),(/"o1","p1"/), &
     &    L1,hv,ResL1,1.0d0,1.0d0,nrep)

       if (iprnt.ge.1) &
     & write(*,*)'RESL1: dia1',dot_product(resL1,resL1)
!----------------------------------------
!    ResL1 = ResL1 - L(m,a)*f_bar_mi(i,m)
!----------------------------------------

        call contraction_222 ((/"o2","p1"/),(/"o1","o2"/),(/"o1","p1"/), &
     &    L1,ho,ResL1,-1.0d0,1.0d0,nrep)

       if (iprnt.ge.1) &
     & write(*,*)'RESL1: dia2',dot_product(resL1,resL1)
!----------------------------------------------
!      ResL1 = ResL1 - G(m,n) * V(mi,na)  
!----------------------------------------------

! second way:

!   w_mnie(mn,ie) = w_mnie(mn,ie) + V(mn,fe)*t(f,i)

      allocate(w_o_temp(iovoot(nrep+1)*rcw)) 

      call getooov(w_o_temp)

      call contraction_424 ((/"o2","o3","p1","p2"/),(/"p2","o1"/), &
     &        (/"o2","o3","o1","p1"/),w_o_temp,-1.0d0, &
     &        1.0d0,nrep,RightTensor=t1)


      call contraction_242 ((/"o2","o3"/),(/"o2","o1","o3","p1"/), &
     & (/"o1","p1"/),G_oo,ResL1,-1.0d0,1.0d0,nrep,RightTensor=w_o_temp)


!!!!! call contraction_422 ((/"o2","o1","o3","p1"/),(/"o2","o3"/), &
!!!!!& (/"o1","p1"/),G_oo,ResL1,-1.0d0,1.0d0,nrep,LeftTensor=w_o_temp)

      deallocate(w_o_temp) 

       if (iprnt.ge.1) &
     & write(*,*)'RESL1: dia3',dot_product(resL1,resL1)
!-------------------------------------------------
!      ResL1 = ResL1 + L1(m,e) * W(ei,ma) 
!second      ResL1 = ResL1 + L1(m,e) * W(ie,am) 
!-------------------------------------------------

!     call contraction_242 ((/"o2","p2"/),(/"p2","o1","o2","p1"/),
!    &   (/"o1","p1"/),L1,ResL1,1.0d0,1.0d0,nrep,RightTensor=w_voov)

!!!!! call contraction_242 ((/"o2","p2"/),(/"o1","p2","p1","o2"/), &
!!!!!&   (/"o1","p1"/),L1,ResL1,1.0d0,1.0d0,nrep,RightTensor=w_voov)


      call contraction_422 ((/"o1","p2","p1","o2"/),(/"o2","p2"/), &
     &   (/"o1","p1"/),L1,ResL1,1.0d0,1.0d0,nrep,LeftTensor=w_voov)


       if (iprnt.ge.1) &
     & write(*,*)'RESL1 after term4',dot_product(resL1,resL1)
!-----------------------------------------------------
!       ResL1 = ResL1 -  L2(mn,ae)*W(ie,mn) 
!-----------------------------------------------------

      call contraction_442((/"o2","o3","p1","p2"/),(/"o1","p2","o2","o3" &
     &          /),(/"o1","p1"/),ResL1,-1.0d0,1.0d0,nrep,LeftTensor=L2,  &
     &          RightTensor=w_ovoo)


       if (iprnt.ge.1) &
     & write(*,*)'RESL1 after term5',dot_product(resL1,resL1)

!-------------------------------------------------------
!      ResL1    = ResL1 + L2(im,ef) * W(ef,am) 
!-------------------------------------------------------

   allocate (l2_sorted1(ndimt2*rcw))

   call srt1c1n(nrep,noot,nvvt,L2,l2_sorted1)

      call srt1t3 (nrep,multb,.false.,nvvt,no,no,nvvot,kvvoot,kkvvot, &
     &           l2_sorted1,buf2)

      deallocate (l2_sorted1)

      call srt19 (nrep,multb,.false.,nvvt,nv,no,nvvot,kvvovt,kkvvot, &
     &               w_vvvo,buf3)

      call cntrct ('t','n',no,nv,nvvot,-a1,buf2,buf3,a1, &
     &    ResL1,nrep)


!!!!! call contraction_442((/"p2","p3","p1","o2"/),(/"o1","o2","p2","p3"/),(/"o1","p1"/),ResL1,1.0d0,1.0d0,nrep, &
!!!!!  &         LeftTensor=w_vvvo,RightTensor=L2)

       if (iprnt.ge.1) &
     & write(*,*)'RESL1 after term6',dot_product(resL1,resL1)

!-----------------------------------------------------------
!      int_temp(e,m) = (G(e,f) * t(f,m) - G(n,m)*t(e,n)) ! corrected from stanton and Gauss 
!      ResL1    = ResL1 + int_temp(e,m) * V(im,ae) 
!------------------------------------------------------------

      allocate(int_temp(nfvo*rcw))

      int_temp = 0.0d0

!     call contraction_222 ((/"p2","p3"/),(/"p3","o2"/),(/"p2","o2"/),
!    &    G_vv,t1,int_temp,1.0d0,1.0d0,nrep)

!     call contraction_222 ((/"o3","o2"/),(/"p2","o3"/),(/"p2","o2"/),
!    &    G_oo,t1,int_temp,-1.0d0,1.0d0,nrep)

!------------------------------------------------------------

!     call contraction_242 ((/"p2","o2"/),(/"o1","o2","p1","p2"/),
!    &   (/"o1","p1"/),int_temp,ResL1,1.0d0,1.0d0,nrep)

      deallocate(int_temp)

       if (iprnt.ge.1) &
     & write(*,*)'RESL1 after term7',dot_product(resL1,resL1)

!---------------------------------------------------------------------------

!second way.

! w_amef(am,ef) = w_amef(am,ef) - V(nm,ef) * t(a,n) 

      allocate (w_v_temp(nv5*rcw))

      call getvovv_incore(w_v_temp)
       
      call contraction_424 ((/"o2","o1","p3","p1"/),(/"p2","o2"/), &
     &    (/"p2","o1","p3","p1"/),w_v_temp,-1.0d0, &
     &    1.0d0,nrep,RightTensor=t1)

!-------------------------------------------------------------------------------



!-------------------------------------------------------
!  out-of-core term....
!    ResL1 = ResL1 - G(e,f) * V(ei,fa) !!wrong
!-------------------------------------------------------

!!!!! call contraction_422 ((/"p2","o1","p3","p1"/),(/"p2","p3"/), &
!!!!!&  (/"o1","p1"/),G_vv,ResL1,-1.0d0,1.0d0,nrep,LeftTensor=w_v_temp)

!right---      
      call contraction_242 ((/"p2","p3"/),(/"p2","o1","p3","p1"/), &
     &  (/"o1","p1"/),G_vv,ResL1,-1.0d0,1.0d0,nrep,RightTensor=w_v_temp)
!right----
!     call contraction_242 ((/"p3","p2"/),(/"p2","o1","p3","p1"/), &
!    &  (/"o1","p1"/),G_vv,ResL1,-1.0d0,1.0d0,nrep,RightTensor=w_v_temp)
       deallocate(w_v_temp)


       if (iprnt.ge.1) &
     & write(*,*)'RESL1 after term8',dot_product(resL1,resL1)
           
      end subroutine 
     
     end module
