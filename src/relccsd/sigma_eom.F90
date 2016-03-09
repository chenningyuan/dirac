  module sigma_eom

!this module calculates the sigma vectors for IP(2h-1p), EA(1h-2p) and EE(2h-2p)  

  use contraction
  use intermediates_1b_2b

!    public eom_ip  
!    public eom_ea
    public eom_ee

    private

   contains

   subroutine sigma_ee_1h1p()

!---------------Description--------------------------------------------
!     Calculates 1h1p type sigma vectors for excitation energy
!-----------------------------------------------------------------------
#include "param.inc"
#include "symm.inc"
#include "complex.inc"
!---------------Calling variables--------------------------------------
  real(8), intent(in)  ::  ! Amplitudes.
  real(8), intent(out) :: w_oooo(rcw*nv1)                ! Intermediate.
!---------------Local Variables--------------------
  real(8), allocatable :: tau(:)
!--------------------------------------------------

!*****
! class : (HBar_SS*R)(a,i)
!*****

!S_ai = S_ai +  F^a_e*R^e_i   
!---------------------------
   
    call contraction_222 ((/"p1","p2"/),(/"p2","o1"/),(/"p1","o1"/), &
   &    hv,r1,sigma1,1.0d0,1.0d0,nrep)

!S_ai = S_ai - f_bar_mi(m,i)*r1(a,m) 
!------------------------------------

    call contraction_222 ((/"o2","o1"/),(/"p1","o2"/),(/"p1","o1"/), &
   &    ho,r1,sigma1,-1.0d0,1.0d0,nrep)

!S_ai = S_ai + w_mbej(ma,ei)*r1(e,m)
!------------------------------------ 

    call contraction_422 ((/"o2","p1","p2","o1"/),(/"p2","o2"/),(/"p1","o1"/),w_ovvo,r1,sigma1,1.0d0,1.0d0,nrep)


!******
! class : [HBar_SD*R](a,i)
!******

!S_ai = S_ai + F_bar_me(m,e)*r2(ae,im)
!--------------------------------------

    call contraction_242 ((/"o2","p2"/),(/"p1","p2","o1","o2"/), &
   &   (/"p1","o1"/),hov,sigma1,1.0d0,1.0d0,nrep,RightTensor=r2)

!S_ai = S_ai + w_amef(am,ef) * r2(ef,im)
!----------------------------------------
 
    call contraction_442((/"p1","o2","p2","p3"/),(/"p2","p3","o1","o2"/),(/"p1","o1"/),sigma1,1.0d0,1.0d0,nrep, &
   &         LeftTensor=w_v_temp,RightTensor=r2)


!S_ai = S_ai - W_mnie(mn,ie)*r2(ae,mn)
!---------------------------------------

    call contraction_442((/"o2","o3","o1","p2"/),(/"p1","p2","o2","o3"/),(/"p1","o1"/),sigma1,-1.0d0,1.0d0,nrep, &
   &         LeftTensor=int_wooop,RightTensor=r2)


!********
! class : [HBar_DS*R](ab,ij)
!********

!S(ab,ij) = S(ab,ij) + P(ab) w_iemn(ma,ij)*r1(b,m)
!--------------------------------------------------

    call contraction_424 ((/"o2","p1","o1","o2"/),(/"p2","o2"/), &
   &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=w_ovoo, &
   &    RightTensor=r1) 

!S(ab,ij) = S(ab,ij) + P(ij) w_efam(ab,ej)*r1(e,i)
!--------------------------------------------------

    call contraction_424 ((/"p1","p2","p3","o2"/),(/"p3","o1"/), &
   &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=w_vvvo, &
   &    RightTensor=r1) 

!int_1b(b,f) = int_1b(b,f) + w_amef (bm,fe)*r1(e,m)
!S(ab,ij) = S(ab,ij) + P(ab) int_1b(b,f)*t2(af,ij)
!---------------------------------------------------

    int_1b = 0.0d0

    call contraction_422 ((/"p2","o3","p4","p3"/),(/"p3","o3"/),(/"p2","p4"/),w_v_temp,r1,int_1b,1.0d0,1.0d0,nrep)

    call contraction_244 ((/"p2","p4"/),(/"p1","p4","o1","o2"/), &
    &
(/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=t2)

!int_1b(n,j) = int_1b(n,j) + w_mnie(nm,je)*r1(e,m)
!S(ab,ij) = S(ab,ij) + int_1b(n,j) * t2(ab,in)
!---------------------------------------------------

    int_1b = 0.0d0

    call contraction_422 ((/"o4","o3","o2","p3"/),(/"p3","o3"/),(/"o4","o2"/),int_woop,r1,int_1b,1.0d0,1.0d0,nrep)

    call contraction_244 ((/"o4","o2"/),(/"p1","p2","o1","o4"/), &
   & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=t2)

!*******************
! class : [HBar_DD*R](ab,ij) 
!*******************

! S(ab,ij) = S(ab,ij) + P(ab) f_bar_ae(b,e) * r2(ae,ij)
!-------------------------------------------------------

    call contraction_244 ((/"p2","p3"/),(/"p1","p3","o1","o2"/), &
   & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=hv,RightTensor=r2)
     
! S(ab,ij) = S(ab,ij) - P(ij) f_bar_mi(m,j)*r2(ab,im) 
!-----------------------------------------------------

    call contraction_244 ((/"o3","o2"/),(/"p1","p2","o1","o3"/), &
   & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=ho,RightTensor=r2)

! S(ab,ij) = S(ab,ij) + w_vvvv(ab,ef)*r2(ef,ij)
!------------------------------------------------

      off1 = 1
      off2 = 1
      do irp = 1, nrep
         istart = 0
         if (nvvt(irp).eq.0) cycle
         mint = nbuf3/nvvt(irp)
  10     call getdz (irp,istart,nint,done,buf3,mint)
         m = nvvt(irp)
         n = noot(irp)
         k = nint
         call xgemm ('n','n',m,n,k,a1,buf3,m,r2(off1+istart*rcw),m, &
     &               a1,sigma2(off2),m)
         if (.not.done) then
            istart = istart + nint
            goto 10
         endif
         off1 = off1 + m * n * rcw
         off2 = off2 + m * n * rcw
      enddo

! S(ab,ij) = S(ab,ij) + w_ijmn(mn,ij) * r2(ab,mn)
!--------------------------------------------------

      call contraction_444((/"o3","o4","o1","o2"/),(/"p1","p2","o1","o2"/), &
     &    (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,RightTensor=r2,LeftTensor=w_oooo)

! S(ab,ij) = S(ab,ij) + P(ab)P(ij) w_mbej(mb,ej)*r2(ae,im)
!----------------------------------------------------------

      call srt1tt5 (nrep,multb,.false.,.false.,nv,nv,no,no,mvo,j2voov,jjvo,jjov, &
     &              r2,buf2)


      call srt1ss4 (nrep,multb,.false.,.false.,no,nv,nv,no,mov,j2ovvo,jjov,jjvo, &
     &              w_voov,buf3)


      allocate (sigma2_temp(j2vovo(nrep+1)*rcw))

      resl2_temp = 0.0d0


      do i = 1, nrep
         j = multb(i+nrep,1+nrep,2)

          m = mvo(i)
          n = mvo(i)
          k = mov(j)

          nd = mvo(i)
          off1 = 1 + j2voov(i)*rcw  
          off2 = 1 + j2ovvo(i)*rcw  
          off3 = 1 + j2vovo(i)*rcw  


            call xgemm ('n','t',m,n,k,a1,buf2(off1),m, &
      &                  buf3(off2),nd,a0,sigma2_temp(off3),m)

      enddo

      call srt1tt4 (nrep,multb,.false.,.true.,nv,nv,no,no,mvo,j2vovo,jjvo,jjvo, &
     &              sigma2_temp,buf2)

      deallocate(sigma2_temp)

      call xaxpy (ndimt2,a1,buf2,1,sigma2,1)


!int_1b(a,f) = int_1b(a,f) - V(nm,fe)*r2(ea,mn)  
!S(ab,ij) = S(ab,ij) + int_1b(a,f)*t2(fb,ij)  
!-----------------------------------------------
    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"p3","p1","o3","o4"/),(/"p1","p4"/),int_1b,-1.0d0,1.0d0,nrep,RightTensor=r2)

    call contraction_244 ((/"p1","p4"/),(/"p4","p2","o1","o2"/), &
    & (/"p1","p2","o1","o2"/),sigma2,1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=t2)

!int_1b(i,n) = int_1b(i,n) + V(nm,fe)*r2(fe,im)
!S(ab,ij) = S(ab,ij) - int_1b(i,n)*t(ab,jn)
!------------------------------------------------

    int_1b = 0.0d0

    call contraction_442 ((/"o4","o3","p4","p3"/),(/"p4","p3","o1","o3"/),(/"o1","o4"/),int_1b,1.0d0,1.0d0,nrep,RightTensor=r2)

    call contraction_244 ((/"o1","o4"/),(/"p1","p2","o2","o4"/), &
    & (/"p1","p2","o1","o2"/),sigma2,-1.0d0,1.0d0,nrep,LeftTensor=int_1b,RightTensor=t2)


  end subroutine

   end module 
