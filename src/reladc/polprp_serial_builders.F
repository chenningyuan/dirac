C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_X_MAIN (daibj,eps,dimphph)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Computes HP/HP block (main block) of the Excitation ADC Matrix
!     (polarization propagator). 
!     Works for real and complex case.
!
!---------------Calling variables--------------------------------------
!
!
      REAL*8,dimension(:)         :: daibj
      REAL*8,dimension(:)         :: eps
      INTEGER                     :: dimphph
!
!---------------Common Blocks, symmetry arrays from CCSETI --------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
!
!---------------Local variables--------------------------------------
!
      interface
        SUBROUTINE denomvvoo_pp(ra1,ra2)
          real*8, dimension(:)  :: ra1,ra2
        END SUBROUTINE
      end interface
!
!---------------Local variables--------------------------------------
!
      real*8, allocatable        :: sa_vvoott(:)
      real*8, allocatable        :: sb_vvoott(:)
      real*8, allocatable        :: buf1(:)
      real*8, allocatable        :: buf2(:)
      real*8, allocatable        :: buf3(:)
      real*8, allocatable        :: buf4(:)
      real*8, allocatable        :: buf5(:)
      real*8, allocatable        :: bab(:)
      real*8, allocatable        :: dij(:)
      real*8, allocatable        :: vovo1(:)

      integer                    :: ialloc
      real*8                     :: qmasymh2
!
!---------------Executable code--------------------------------------
!
      IF(CARITH) THEN
        CALL PST('Building complex HP/HP block+')
      ELSE
        CALL PST('Building real HP/HP block+')
      ENDIF
      if(carith.and.(rcw.ne.2)) then
        call quit('Internal error in BUILD_X_MAIN!')
      endif
      if(.not.carith.and.(rcw.ne.1)) then
        call quit('Internal error in BUILD_X_MAIN!')
      endif
!
!******************* preparations
!
!  often needed v_ac,kl array which is V_ac,kl (VVOO)
!  divided by the energy denominator
!  e_a + e_c - e_k - e_l and the original V_ac,kl
!  we therefore calculate it once and store it in this subroutine.
!  the variables icX count the number of integrals. For the actual
!  storage requirements one has to multiply with RCW!
!
      ic1 = ivvoott(nrep+1)
      allocate(sa_vvoott(ic1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'allocation problem 1 in BUILD_X_MAIN'
      sa_vvoott = 0.0d0
      allocate(sb_vvoott(ic1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'allocation problem 2 in BUILD_X_MAIN'
      sb_vvoott = 0.0d0

      call getvvoo(sa_vvoott)
      call xcopy(ic1,sa_vvoott,1,sb_vvoott,1)
      call denomvvoo_pp(eps,sa_vvoott)
!
! we have v_ac,kl in sa_vvoott and
! V_ac,kl         in sb_vvoott
!
!  allocate all remaining arrays here in this section
!  and not throughout the rest of the code. Otherwise it's error-prone!
!  for BAB: since the delta_ij enforces gamma_i = gamma_j only the totally
!  symmetric entries in BAB will survive!
!
      ic2 = NFVV
      allocate(bab(ic2*rcw))
      bab = 0.0d0

      ic3 = max(lvvoot(nrep+1),kvvoot(nrep+1))
      allocate(buf1(ic3*rcw))  ! array has to accomodate A,BIJ and ABI,J storage
      allocate(buf2(ic3*rcw))
      buf1 = 0.0d0; buf2 = 0.0d0

      ic4 = NFOO
      allocate(dij(ic4*rcw)); 
      dij = 0.0d0

      ic5 = max(j2vovo(nrep+1),jvovo(nrep+1))
      allocate(buf3(ic5*rcw))  ! all arrays have to hold VOVO integrals
      allocate(buf4(ic5*rcw))  ! in MVO*MVO instead of NVO*NVO manner!
      allocate(buf5(ic5*rcw))
      buf3 = 0.0d0; buf4 = 0.0d0; buf5 = 0.0d0

      ic6 = ivovo(nrep+1)      ! size of the generic VOVO array in generic BKC
      allocate(vovo1(ic6*rcw))  ! in MVO*MVO instead of NVO*NVO manner!
      vovo1 = 0.0d0
!
!
!_____________________________________________________
!|
!| Calculate B_ab intermediate   TERM C(2)A
!| result in bab  <a||b> BKC
!
!  sort v_ac,kl  -->  v_a,ckl    (buf1)
!  sort V*_bc,kl  -->  V*_b,ckl  (buf2), c.c. will bw accounted for
!                                in the corresponding contrct call!
!  v_a,ckl  *  V*_b,ckl   -->  bab   <||>
!  'hermitize'  the result in bab (add the transposed cc)
!  **att** the array bab is packed as nv1*nv1(beta1),...
!
      call srt1t2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             sa_vvoott,buf1)
      call srt1t2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             sb_vvoott,buf2)
      call cntrct ('N','C',nv,nv,nvoot,a1,buf1,buf2,a1,bab,nrep)
      do krep=1,nrep
        ioff=jjvv(krep,krep)*rcw + 1
        if(carith) then
          call cmatsym(bab(ioff),nv(krep),nv(krep))
        else
          call matsym(bab(ioff),nv(krep),nv(krep))
        endif
      enddo
!
!_____________________________________________________
!|
!| Calculate D_ij intermediate   TERM C(2)B
!| result in dij  |i><j|
!
!  sort v_cd,ki   -->  v_cdk,i  <c|<d||k>,|i>  (buf1)
!  V_cd,kj        -->  V*_cd,kj
!  sort V*_cd,kj   -->  V*_cdk,j  |c>|d><k|,<j|  (buf2)
!  v_cdk,i  *  V*_cdk,j  -->  Dij   -->  'hermitize'
!  Dij   has |i><j| 
!
      buf1 = 0.0d0; buf2 = 0.0d0
      call srt1t3 (nrep,multb,lfa,nvvt,no,no,nvvot,kvvoot,kkvvot,
     &             sa_vvoott,buf1)
      if(carith) call conjuga(ic1,sb_vvoott,1)
      call srt1t3 (nrep,multb,lfa,nvvt,no,no,nvvot,kvvoot,kkvvot,
     &             sb_vvoott,buf2)
      if(carith) call conjuga(ic1,sb_vvoott,1)
      call cntrct ('T','N',no,no,nvvot,a1,buf1,buf2,a0,dij,nrep)
      do krep=1,nrep
        ioff=jjoo(krep,krep)*rcw + 1
        if(carith) then
          call cmatsym(dij(ioff),no(krep),no(krep))
        else
          call matsym(dij(ioff),no(krep),no(krep))
        endif
      enddo
! ****************** why this conjugation ???? ***************
      if(carith) call conjuga(ic4,dij,1)
!
!_____________________________________________________
!|
!|  Calculate TERM C(2)C, part[A] herm.    |
!  conjugation required, result in buf3   |
!
!  sort v_ac,ik    -->  v_ai,ck    (buf3)
!  V*_bc,jk = V*_cb,kj
!  sort V*_cb,kj   -->  V*_ck,bj   (buf4)
!
      call srt1tt4 (nrep,multb,lfa,lfa,nv,nv,no,no,
     &              mvo,j2vovo,jjvo,jjvo,sa_vvoott,buf3)
      if(carith) call conjuga(ic1,sb_vvoott,1)
      call srt1tt4 (nrep,multb,lfa,lfa,nv,nv,no,no,
     &              mvo,j2vovo,jjvo,jjvo,sb_vvoott,buf4)
      if(carith) call conjuga(ic1,sb_vvoott,1)
!
!      CALL CNTRCT('N','T',MVO,MVO,MVO,A1,buf3,buf4,A0,buf5,NREP)
!
!
!  here again special contraction of v_ai,ck  X  V*_ck,bj
!  buf3  *  buf4. Result has BKC:  <a||i>,|b><j| as the storage mode.
!
      ioff1 = 1
      ioff3 = 1
      do irep = 1, nrep
         jrep = multb(irep+nrep,1+nrep,2)
         if(multb(jrep+nrep,irep+nrep,1).ne.1)
     &      stop 'BRUTAL FUCK'
         m = mvo(irep)
         n = mvo(irep)
         k = mvo(jrep)
         ioff2 = j2vovo(jrep) * rcw + 1
         call xgemm ('N','N',m,n,k,a1,buf3(ioff1),m,buf4(ioff2),n,
     &               a0,buf5(ioff3),m)
         ioff1 = ioff1 + m * k * rcw
         ioff3 = ioff3 + m * n * rcw
      enddo
!
! This is now the first half contribution to C(2)C. Part[B] is the hermitian
! conjugated one. Multiply result by -0.5 according to the equations.
!
      do krep=1,nrep
!       write(*,*) 'JVOVO/J2VOVO of krep:',jvovo(krep),j2vovo(krep)
        ioff=jvovo(krep)*rcw + 1
        if(carith) then
          call cmatsym(buf5(ioff),mvo(krep),mvo(krep))
        else
          call matsym(buf5(ioff),mvo(krep),mvo(krep))
        endif
      enddo
      call xscal(jvovo(nrep+1),-ap5,buf5,1)
!
!  buf5 now contains C(2)C in <a||i>,|b><j| storage mode (Y_ai,bj).
!  Do intermediate resort of Y_ai,bj (buf5) ---->  Y_aj,bi  (buf3) with *generic* VOVO BKC.
!  in SRT 16 the *second* buffer is the source buffer when it comes to inverse !  sorting.
!
      buf3 = 0.0d0
      call srt16 (nrep,multb,lfa,ltr,nv,no,nv,no,mvo,jvovo,jjvo,jjvo,
     &            buf3,buf5)
!
!  ----------------------------------
!  Then we can add C(2)C and C(1)   |
!  due to identical BKC structure   |
!  ----------------------------------
!
      call getvovo(vovo1)
      call xscal(ic6,-A1,vovo1,1)
      call xaxpy (ic6,A1,buf3,1,vovo1,1)
!
!--------------------------------------------------------
!  'weave' into the C(0), C(2A) and C(2B) contributions |
!  as diagonal energies, 0.5*Bab and 0.5*Dij            |
!--------------------------------------------------------
!
      if(carith) then
        call bweaverc(vovo1,eps,bab,dij,ic6,ic2,ic4)
      else
        call bweaverr(vovo1,eps,bab,dij,ic6,ic2,ic4)
      endif
!
!--------------------------------------------------------
!  vovo1 now contains all the C(2) ai,bj contributions  |
!  we sort it back now to the PHPH storage convention.  |
!  This is now a forward sort <<>>  -->  <>><           |
!--------------------------------------------------------
!
      CALL SRT16(NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,
     &           JJVO,JJVO,vovo1,daibj)

      deallocate(vovo1)
      deallocate(buf5)
      deallocate(buf4)
      deallocate(buf3)
      deallocate(dij)
      deallocate(buf2)
      deallocate(buf1)
      deallocate(bab)
      deallocate(sb_vvoott)
      deallocate(sa_vvoott)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BWEAVERR(bphph,eps,bab,dij,ic6,ic2,ic4)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Crawl through the VOVO structure and distribute the
!     C(0), C(2A) and C(2B) contributions to the phph block
!
!---------------Calling variables--------------------------------------
!
      integer, intent(in)         :: ic6, ic2, ic4
      real*8, intent(inout)       :: bphph(ic6)
      real*8, intent(in)          :: eps(*)
      real*8, intent(in)          :: bab(ic2)
      real*8, intent(in)          :: dij(ic4)
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      real*8                           ::  esum
!
!---------------Executable code--------------------------------------
!
!  the loop structure follows the generic VO||VO storage with <|<|  |>|>
!
      AIBJ = 0
      DO BJRP = 1, NREP
      DO 10 JRP = 1, NREP
       BRP = MULTB(JRP,BJRP+NREP,2)
       DO J = 1, NO(JRP)
       DO B = 1, NV(BRP)
         DO 20 IRP = 1, NREP
          ARP = MULTB(IRP,BJRP+NREP,2)
          DO I = 1, NO(IRP)
          DO A = 1, NV(ARP)
            AIBJ = AIBJ + 1  ! count elements
!
!  from here we know arep,irep,brep,jrep, a,i,b,j
!
!
!  treat the diagonal element which is C(0)=(e_a-e_i)*d_ab*d_ij
!
            IF(ARP.eq.BRP.and.IRP.eq.JRP.and.a.eq.b.and.i.eq.j) THEN
              AOFF = IO(NREP+1) + IV(ARP) + A
              IOFF = IO(IRP) + I
              esum = EPS(AOFF) - EPS(IOFF)
              bphph(aibj) = bphph(aibj) + esum
            ENDIF
!
!  treat the C(2A) contribution: C(2A) = +0.5*Delta_ij*B_ab
!
            IF(IRP.eq.JRP) THEN
              IF(I.eq.J) THEN
                if(ARP.NE.BRP) call quit('bweaverr_error C(2A)')
                aboff = IVV(BRP) + (A-1)*NV(BRP) + B
                bphph(aibj) = bphph(aibj) + 0.5d0*bab(aboff)
              ENDIF
            ENDIF
!
!  treat the C(2B) contribution: C(2B) = +0.5*Delta_ab*D_ij
!
            IF(ARP.eq.BRP) THEN
              IF(A.eq.B) THEN
                if(IRP.NE.JRP) call quit('bweaverr_error C(2B)')
                ijoff = IOO(JRP) + (I-1)*NO(JRP) + J
                bphph(aibj) = bphph(aibj) + 0.5d0*dij(ijoff)
              ENDIF
            ENDIF
!
!  all contributions "woven into"
!
            ENDDO
            ENDDO
 20      CONTINUE
       ENDDO
       ENDDO
 10   CONTINUE
      ENDDO  !BJRP
!
      if(aibj.ne.ivovo(nrep+1)) then
        call quit('VOVO counting error in BWEAVERR!')
      endif

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BWEAVERC(bphph,eps,bab,dij,ic6,ic2,ic4)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Count through the VOVO integrals and distribute the
!     C(0), C(2A) and C(2B) contributions to the phph block
!
!---------------Calling variables--------------------------------------
!
      integer, intent(in)         :: ic6, ic2, ic4
      complex*16, intent(inout)   :: bphph(ic6)
      real*8, intent(in)          :: eps(*)
      complex*16, intent(in)      :: bab(ic2)
      complex*16, intent(in)      :: dij(ic4)
      complex*16,parameter        :: coha=(0.5d0,0.0d0)
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      real*8                           ::  esum
!
!---------------Executable code--------------------------------------
!
!  the loop structure follows the generic VO||VO storage with <|<|  |>|>
!
      AIBJ = 0
      DO BJRP = 1, NREP
      DO 10 JRP = 1, NREP
       BRP = MULTB(JRP,BJRP+NREP,2)
       DO J = 1, NO(JRP)
       DO B = 1, NV(BRP)
         DO 20 IRP = 1, NREP
          ARP = MULTB(IRP,BJRP+NREP,2)
          DO I = 1, NO(IRP)
          DO A = 1, NV(ARP)
            AIBJ = AIBJ + 1  ! count elements
!
!  from here we know arep,irep,brep,jrep, a,i,b,j
!
!
!  treat the diagonal element which is C(0)=(e_a-e_i)*d_ab*d_ij
!
            IF(ARP.eq.BRP.and.IRP.eq.JRP.and.a.eq.b.and.i.eq.j) THEN
              AOFF = IO(NREP+1) + IV(ARP) + A
              IOFF = IO(IRP) + I
              esum = EPS(AOFF) - EPS(IOFF)
              bphph(aibj) = bphph(aibj) + dcmplx(esum,0.0d0)
            ENDIF
!
!  treat the C(2A) contribution: C(2A) = +0.5*Delta_ij*B_ab
!
            IF(IRP.eq.JRP) THEN
              IF(I.eq.J) THEN
                if(ARP.NE.BRP) call quit('bweaverc_error C(2A)')
                aboff = IVV(BRP) + (A-1)*NV(BRP) + B
                bphph(aibj) = bphph(aibj) + coha*bab(aboff)
              ENDIF
            ENDIF
!
!  treat the C(2B) contribution: C(2B) = +0.5*Delta_ab*D_ij
!
            IF(ARP.eq.BRP) THEN
              IF(A.eq.B) THEN
                if(IRP.NE.JRP) call quit('bweaverc_error C(2B)')
                ijoff = IOO(JRP) + (I-1)*NO(JRP) + J
                bphph(aibj) = bphph(aibj) + coha*dij(ijoff)
              ENDIF
            ENDIF
!
!  all contributions "woven into"
!
            ENDDO
            ENDDO
 20      CONTINUE
       ENDDO
       ENDDO
 10   CONTINUE
      ENDDO  !BJRP
!
      if(aibj.ne.ivovo(nrep+1)) then
        call quit('VOVO counting error in BWEAVERC!')
      endif

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SOLVE_PHPH_R(daibj,n,iounit)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Diagonalizes the phph block in each symmetry separately.
!     For the singles amplitudes a local configuration analysis is
!     provided for convenience.
!
!---------------Calling variables--------------------------------------
!
      real*8,dimension(:)         :: daibj
      integer                     :: n,iounit
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      real*8, allocatable         :: amat(:,:),xmat(:)
      integer                     :: desrep,jrep,brep,trep
      integer                     :: ioff,nsiz,ialloc,icount
      integer                     :: nx,lfun,kfun,jabs,babs
      integer                     :: resh_array(2)
      real*8, allocatable         :: evl(:),work(:)
      integer                     :: lwork,info
      real*8,parameter            :: AUTOEV = 27.2113957D0
      real*8                      :: x2
      character*26,allocatable    :: cfield(:)
      character*11                :: davstartname
      character*8                 :: dsroot='DAVSTART'
!
!---------------Interface region-------------------------------------
!
      interface

        function qmasymr(ra1,i1)
          real*8                 :: qmasymr
          real*8,dimension(*)    :: ra1
          integer                :: i1
        end function qmasymr

      end interface
!
!---------------Executable code--------------------------------------
!
      call PST('Generating (real) Davidson start vectors+')
!
!  create local single configuration descriptor array
!  remember: first index is virtual^*, second is occupied
!  this corresponds to the storage mode!
!
      nx = jvoxx(nrep+1)
      allocate(cfield(nx))

      icount = 0
      do desrep=1,nrep   !desrep is a bosonic loop
        do jrep = 1,nrep   !jrep,brep are fermionic
          brep = multb(desrep+nrep,jrep,2)
          if(multb(brep,jrep,2).ne.desrep) then
            write(*,*) 'HP/HP ordering problem!'
            stop 'solve_phph_r'
          endif
          do jfun = 1,NO(jrep)
          do bfun = 1,NV(brep)
            jabs=IO(jrep) + jfun
            babs=IO(nrep+1) + IV(brep) + bfun
            icount = icount + 1
            write(cfield(icount),777) jabs,repna(jrep),
     &                                babs,repna(brep)
          enddo
          enddo
        enddo
      enddo
 777  Format(I3,' (',A4,') --->',1X,I3,' (',A4,')')
      if(icount.ne.jvoxx(nrep+1)) then
         write(*,*) icount,jvoxx(nrep+1)
         stop 'conf array error!'
      endif
!     write(*,*) 'HP configuration counting successful'
!
!  diagonalize phph submatrices
!  and write corresponding orthogonal eigenvectors
!  as starting vectors for the Davidson algorithm.
!
!___________________________________________________________
!|
!|
!|
      do desrep = 1,nrep

        ioff = JVOVO(desrep)   ! calculate total offset in corresponding irrep
        nsiz = MVO(desrep)
        if(nsiz.eq.0) then
          write(*,'(1X,A27,6X,I3,3X,A4)') 'No states in symmetry',
     &          desrep,repna(desrep)
          cycle
        endif
!       write(*,*) 'irep/name/offset/# of jvovo/size of main block:',
!    &        desrep,repna(desrep),
!    &        ioff,JVOVO(desrep+1) - JVOVO(desrep),MVO(desrep)
        write(*,*)
!
!  open corresponding file
!
        IF(desrep.GT.9) THEN
          WRITE(davstartname,'(A8,A1,I2)') dsroot,'.',desrep
        ELSE
          WRITE(davstartname,'(A8,A2,I1)') dsroot,'.0',desrep
        ENDIF
        OPEN(iounit,file=davstartname,access='sequential',
     &       status='unknown',form='unformatted')

        allocate(amat(nsiz,nsiz),stat=ialloc)
        if(ialloc.ne.0) STOP 'Allocation problem 1 in SOLVE_PHPH_R'
        allocate(xmat(nsiz*nsiz),stat=ialloc)
        if(ialloc.ne.0) STOP 'Allocation problem 2 in SOLVE_PHPH_R'
        amat = 0.0d0; xmat = 0.0d0
!
!  we convert the linear MVO*MVO array to quadratic array
!  hereby we need xmat, compiler does not understand it otherwise
!
        resh_array(1)=nsiz; resh_array(2)=nsiz
        xmat(1:nsiz*nsiz) = daibj(ioff+1:ioff + nsiz*nsiz)
        amat = reshape(xmat,resh_array)
!
!  check asymmetry of obtained PHPH matrix
!
        x2 = dabs(qmasymr(amat,nsiz))
        if(x2.gt.1.0E-10) then
          write(*,*) 'Warning: Asymmetry of HP/HP matrix',
     &         ' in symmetry',desrep,' amounts to',x2
!       else
!         write(*,*) 'Matrix symmetric with',x2
        endif
!
!  amat array has ph/ph structure and is transferred to LAPACK
!
        lwork=5*nsiz
        allocate(evl(nsiz));allocate(work(lwork))
        CALL DSYEV('V','L',nsiz,amat,nsiz,evl,work,lwork,info)
        if(info.ne.0) then
          write(*,*) 'DSYEV problem!  Info code:',info
          stop
        endif
!
!  write out eigenvalues and largest eigenvectors for the current symmetry
!  on screen.
!
!       write(*,*)
        write(*,'(1X,A27,6X,I3,3X,A4,A)') 'Vectors available for symm.',
     &            desrep,repna(desrep),'  (max. 10 are printed)'
        write(*,*)
        do i=1,min(nsiz,10)
          write(*,'(I5,2F22.12,2X,A4,X,A2)') i,evl(i),evl(i)*autoev,
     &             repna(desrep),'@X'
!         write(*,'(5X,A)') '-----------------------------------------'
!         do j=1,nsiz
!           if(dabs(amat(j,i)).gt.0.1d0) then
!             ioff=jvoxx(desrep)
!             write(*,*) j,':  ',cfield(ioff+j),amat(j,i)
!           endif
!         enddo
!         write(*,*)
        enddo
!
!  write out start vectors to specific file name.
!
        Write(iounit) nsiz
        Do i=1,nsiz
          Write(iounit) (amat(j,i),j=1,nsiz)
        Enddo
        Close(iounit)
        
!
!  deallocate all auxiliary arrays and go to next symmetry
!
        deallocate(work)
        deallocate(evl)
        deallocate(xmat)
        deallocate(amat)

      enddo  !desrep
      deallocate(cfield)
      write(*,*) '...success.'
!|
!|
!|
!|__________________________________________________________


      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SOLVE_PHPH_C(daibj,nele,iounit)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Diagonalizes the hermitian PHPH block in each symmetry separately.
!     For the singles amplitudes a local coniguration analysis is
!     provided for convenience.
!     Due to the more restrictive F90 compiler with function prototypes
!     we can not reinterpret the real array as a complex one. We introduce
!     an *intermediate complex* array. The main block size is small enough
!     not to run into trouble.
!     nele is the number of elements in the array and *not* the physical
!     extension of the real array. The real array has an allocation of
!     2*nele !!
!
!---------------Calling variables--------------------------------------
!
      real*8,dimension(:)         :: daibj
      integer                     :: nele,iounit
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      complex*16, allocatable, dimension(:,:)      :: amat
      complex*16, allocatable, dimension(:)        :: xmat,daibjc
      integer                             :: desrep,jrep,brep,trep
      integer                             :: ioff,nsiz,ialloc,icount
      integer                             :: nx,lfun,kfun,jabs,babs
      integer                             :: resh_array(2)
      integer                             :: i,ix
      real*8, allocatable, dimension(:)        :: evl,rwork
      complex*16, allocatable, dimension(:)    :: work
      integer                             :: lwork,info
      real*8,parameter                    :: AUTOEV = 27.2113957D0
      real*8                              :: x2
      character*26,allocatable, dimension(:)   :: cfield
      character*11                        :: davstartname
      character*8                         :: dsroot='DAVSTART'

      real*8,parameter        :: rzero=(0.0d0)
      complex*16,parameter    :: czero=(rzero,rzero)
!
!---------------Interface region-------------------------------------
!
      interface

        function qmasymh(ca1,i1)
          real*8                       :: qmasymh
          complex*16,dimension(:,:)    :: ca1
          integer                      :: i1
        end function qmasymh

        subroutine ZHEEV(c1,c2,i1,ca1,i2,ra1,ca2,i3,ra2,i4)
           character*1                         :: c1,c2
           integer                             :: i1,i2,i3,i4
           double complex, dimension(i2,*)     :: ca1
           double precision, dimension(*)      :: ra1,ra2
           double complex, dimension(*)        :: ca2
        end subroutine

      end interface
!
!---------------Executable code--------------------------------------
!
      call PST('Generating (complex) Davidson start vectors+')

      if(rcw.ne.2) then
        write(*,*) 'Variable RCW is not 2 despite complex mode!'
        call quit('*** internal consistency violated ***')
      endif
!
!  create local single configuration descriptor array
!  remember: first index is virtual^*, second is occupied
!  this corresponds to the storage mode!
!
      nx = jvoxx(nrep+1)
      allocate(cfield(nx))  !character field
!___________________________________________
!|
!|  configuration generator
!|
      icount = 0
      do desrep=1,nrep   !desrep is a bosonic loop
        do jrep = 1,nrep   !jrep,brep are fermionic
          brep = multb(desrep+nrep,jrep,2)
          if(multb(brep,jrep,2).ne.desrep) then
            write(*,*) 'HP/HP ordering problem!'
            stop 'SOLVE_PHPH_C'
          endif
          do jfun = 1,NO(jrep)
          do bfun = 1,NV(brep)
            jabs=IO(jrep) + jfun
            babs=IO(nrep+1) + IV(brep) + bfun
            icount = icount + 1
            write(cfield(icount),777) jabs,repna(jrep),
     &                                babs,repna(brep)
          enddo
          enddo
        enddo
      enddo
 777  Format(I3,' (',A4,') --->',1X,I3,' (',A4,')')
      if(icount.ne.jvoxx(nrep+1)) then
         write(*,*) icount,jvoxx(nrep+1)
         stop 'conf array error!'
      endif
!      write(*,*) 'HP configuration counting successful'
!|
!|  end configuration generator
!|__________________________________________
!
!  diagonalize PHPH submatrices
!  and write corresponding complex orthogonal eigenvectors
!  as starting vectors for the Davidson algorithm.
!
!___________________________________________________________
!|
!|  generate complex intermediate array from original real one.

      allocate(daibjc(nele))
      daibjc = czero
      do i=1,nele
        ix = 2*(i-1) + 1
        daibjc(i)=dcmplx(daibj(ix),daibj(ix+1))
      enddo
!|
!|__________________________________________
!
!___________________________________________________________
!|
!| form quadratic matrix in specific symmetry
!| from total linear array and diagonalize
!|
      do desrep = 1,nrep

        ioff = JVOVO(desrep)   ! calculate total offset in corresponding irrep
        nsiz = MVO(desrep)
        if(nsiz.eq.0) then
          write(*,'(1X,A27,6X,I3,3X,A4)') 'No states in symmetry',
     &          desrep,repna(desrep)
          cycle
        endif
!
!  open corresponding file
!
        IF(desrep.GT.9) THEN
          WRITE(davstartname,'(A8,A1,I2)') dsroot,'.',desrep
        ELSE
          WRITE(davstartname,'(A8,A2,I1)') dsroot,'.0',desrep
        ENDIF
        OPEN(iounit,file=davstartname,access='sequential',
     &       status='unknown',form='unformatted')

        allocate(amat(nsiz,nsiz),stat=ialloc)
        if(ialloc.ne.0) STOP 'Allocation problem 1 in SOLVE_PHPH_C'
        allocate(xmat(nsiz*nsiz),stat=ialloc)
        if(ialloc.ne.0) STOP 'Allocation problem 2 in SOLVE_PHPH_C'
        amat = czero
        xmat = czero
!
!  we convert the linear MVO*MVO array to quadratic array
!  hereby we need xmat, compiler does not understand it otherwise
!
        resh_array(1)=nsiz; resh_array(2)=nsiz
        xmat(1:nsiz*nsiz) = daibjc(ioff+1:ioff + nsiz*nsiz)
        amat = reshape(xmat,resh_array)
!
!  check deviation from hermiticity of obtained PHPH matrix
!
        x2 = dabs(qmasymh(amat,nsiz))
        if(x2.gt.1.0E-10) then
          write(*,*) 'Warning: Deviation from hermiticity of HP/HP',
     &       ' matrix in symmetry',desrep,' amounts to',x2
        else
          write(*,*) 'Hermiticity of HP/HP block:',x2
        endif
        write(*,*)
!
!  amat array has ph/ph structure and is transferred to LAPACK
!
        lwork=5*nsiz
        allocate(evl(nsiz))
        allocate(work(lwork))
        allocate(rwork(lwork))
        CALL ZHEEV('V','L',nsiz,amat,nsiz,evl,work,lwork,rwork,info)
        if(info.ne.0) then
          write(*,*) 'ZHEEV problem in SOLVE_PHPH_C! code:',info
          call quit('numerical problem')
        endif
!
!  write out eigenvalues and largest eigenvectors for the current symmetry
!  on screen.
!
!       write(*,*)
        write(*,'(1X,A27,6X,I3,3X,A4,A)') 'Vectors available for symm.',
     &            desrep,repna(desrep),'  (max. 10 are printed)'
        write(*,*)
        do i=1,min(nsiz,10)
          write(*,'(I5,2F22.12,2X,A4,X,A2)') i,evl(i),evl(i)*autoev,
     &             repna(desrep),'@X'
        enddo
!
!  write out start vectors to specific file name.
!
        Write(iounit) nsiz
        Do i=1,nsiz
          Write(iounit) (amat(j,i),j=1,nsiz)
        Enddo
        Close(iounit)
        
!
!  deallocate all auxiliary arrays and go to next symmetry
!
        deallocate(work)
        deallocate(rwork)
        deallocate(evl)
        deallocate(xmat)
        deallocate(amat)

      enddo
      deallocate(cfield)
      deallocate(daibjc)
      write(*,*) '...success.'
!|
!|
!|
!|__________________________________________________________

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_X_COUSAT_R (PHPH_BLOCK,EPS,DESREP,LENCOUPLBL,
     &                             INTBUF,NBUFS,DOXEXT,WTHR,IOCH)
!
      IMPLICIT INTEGER (A-Z)
!     IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!  Computes (real) 2H2P/HP and SAT block of the exc. ADC Matrix
!  In the corresponding symmetry DESREP.
!
!---------------Calling variables--------------------------------------
!
      REAL*8                           :: PHPH_BLOCK(:)
      REAL*8                           :: EPS(:)
      INTEGER                          :: DESREP
      INTEGER                          :: LENCOUPLBL,INTBUF,NBUFS
      LOGICAL                          :: DOXEXT
      REAL*8                           :: WTHR
      INTEGER                          :: IOCH
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/files.inc"
#include  "../relccsd/param.inc"
!
!---------------Local variables--------------------------------------
!
      INTEGER              :: n1,irep,k
      INTEGER              :: jrep,brep,jfun,bfun
      INTEGER              :: ialloc,mjcol,nmvod

      INTEGER              :: mxno,mxnv,ioff
      INTEGER, dimension (:,:,:), allocatable  :: oolt
      INTEGER, dimension (:,:,:), allocatable  :: vvlt

      REAL*8, dimension(:), allocatable :: vooo1
      REAL*8, dimension(:), allocatable :: vovv1
      REAL*8, dimension(:), allocatable :: bufcol,bufcolex
!
!  variables for writing out ADC matrix
!
      REAL*8,  dimension(:), allocatable :: wbufr
      INTEGER, dimension(:), allocatable :: wbufi, wbufj
      INTEGER              :: IBFP
      CHARACTER*6          :: NAMEX = 'ADCXPM'
!
!---------------Interface area --------------------------------------
!
      interface

        INTEGER FUNCTION XCOLDET(i1)
          INTEGER                         :: i1
        END FUNCTION XCOLDET

        SUBROUTINE XCOLMAKE(i1,i2,ra1,
     &                      i3,i4,i5,i6,
     &                      i7,i8,
     &                      ra2,ra3,ia1,ia2)
          INTEGER                         :: i1,i2
          REAL*8, dimension(:)            :: ra1
          INTEGER                         :: i3,i4,i5,i6
          INTEGER                         :: i7,i8
          REAL*8, dimension(:)            :: ra2,ra3
          INTEGER, dimension(:,:,:)       :: ia1,ia2
        END SUBROUTINE XCOLMAKE

        SUBROUTINE XCOLWRIT(i1,r1,ra1,ra2,ia1,ia2,i2,i3,i4,i5,i6)
          INTEGER                         :: i1
          REAL*8                          :: r1
          REAL*8, dimension(:)            :: ra1,ra2
          INTEGER, dimension(:)           :: ia1,ia2
          INTEGER                         :: i2,i3,i4,i5,i6
        END SUBROUTINE XCOLWRIT

        SUBROUTINE XSATMAKE(i1,i2,i3,i4,i5,r1,
     &              ra1,ia1,ia2,ra2,ia3,ia4,
     &              l1)
          INTEGER                   :: i1,i2,i3,i4,i5
          REAL*8                    :: r1
          REAL*8, dimension(:)      :: ra1
          INTEGER, dimension(:)     :: ia1,ia2
          REAL*8, dimension(*)      :: ra2
          INTEGER,dimension(:,:,:)  :: ia3,ia4
          LOGICAL                   :: l1
        END SUBROUTINE XSATMAKE

        SUBROUTINE XCOLFLSH(i1,i2,i3,i4,ra1,ia1,ia2)
          INTEGER                         :: i1,i2,i3,i4
          REAL*8, dimension(:)            :: ra1
          INTEGER, dimension(:)           :: ia1,ia2
        END SUBROUTINE XCOLFLSH

      end interface
!
!---------------Executable code--------------------------------------
!
!________________________________________________________
!|     allocate ADC MATRIX WRITE BUFFERS and counters
!|
!|
      allocate(wbufr(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 1 could not be allocated.'
      allocate(wbufi(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 2 could not be allocated.'
      allocate(wbufj(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 3 could not be allocated.'
      wbufr=0.0d0; wbufi = 0; wbufj = 0
      IBFP = 0   ! pointer to actual buffer level
      NBUFS = 0  ! counter for full buffers
!|
!|
!|_______________________________________________________
!________________________________________________________
!|     create OOT and VVT lookup tables
!|     tables are small, computation is too clumsy
!|
      mxno = 0; mxnv = 0
      do irep=1,nrep  
        if(no(irep).gt.mxno) mxno=no(irep)
        if(nv(irep).gt.mxnv) mxnv=nv(irep)
      enddo
      allocate(oolt(mxno,mxno,nrep),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #a in BUILD_X_COUSAT_R'
      allocate(vvlt(mxnv,mxnv,nrep),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #b in BUILD_X_COUSAT_R'
      OOLT = 0
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NO(IREP)
          DO I=J+1,NO(IREP)
            OOLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
      VVLT = 0
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NV(IREP)
          DO I=J+1,NV(IREP)
            VVLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
!|
!|
!|_______________________________________________________
!
!
!  Attention*** Some arrays that have to store *huge* numbers 
!  Are defined as Integer*8 in symm.inc ! The generic sorters
!  operate with unspecified INTEGER declarations. Transferring such
!  an array as an argument to the SRT_xxx results in an erroneous
!  behaviour of the sorting routine.
!
! fetch <VO||OO> integrals needed for coupling block 
! 
      n1=ivooot(nrep+1)
      allocate(vooo1(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error vooo1!'
      call getvooo(vooo1)
!
! fetch <VO||VV> integrals needed for coupling block 
! 
      n1=ivovvt(nrep+1)
      allocate(vovv1(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error vovv1!'
      call rdvovv (vovv1)
!
! length of main block:
!
      NMVOD = MVO(desrep)
!
! determine actual length of coupling block
! and allocate column buffer
!
      lencouplbl =  XCOLDET(desrep)
      WRITE(IW,*) '   *)) Current symmetry:          ',desrep
      WRITE(IW,*) '   *)) Block length of main block:',NMVOD
      WRITE(IW,*) '   *)) Block length of coupling block:',lencouplbl
      WRITE(IW,'(A,I8,A)') '    *)) Matrix dimension: [',
     &        lencouplbl + NMVOD,' ]**2'

      ALLOCATE(BUFCOL(lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #3 in BUILD_X_COUSAT_R'
      ALLOCATE(BUFCOLEX(NMVOD + lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #4 in BUILD_X_COUSAT_R'
!
!  according to the SRT16 convention, 
!  the columns also have DESREP symmetry. 
!  We loop over the columns in the same manner and construct
!  the columns of the 2H2P/HP coupling block. 
!  after construction we write out the column together with the
!  components of the main block.
!
!     write(*,*) '   ======= cousat desrep:',desrep,' ============='
!     write(*,*)
      mjcol = 1
      do jrep = 1,nrep

        brep = multb(desrep+nrep,jrep,2)

!       write(*,*) 'brep,jrep:',brep,jrep

        do jfun = 1,NO(jrep)
        do bfun = 1,NV(brep)
!       write(*,*) '         B,J:',bfun,jfun

          bufcol   = 0.0d0
          bufcolex = 0.0d0

          CALL XCOLMAKE(desrep,lencouplbl,bufcol,
     &                  jrep,brep,jfun,bfun,
     &                  mxno,mxnv,
     &                  vooo1,vovv1,oolt,vvlt)
!
! combine main block column with coupling block column
! istart is the start address for picking the right main block
! data to corresponding symetry DESREP.
!
          istart = JVOVO(desrep) + (mjcol-1)*nmvod
          bufcolex(1:nmvod)=phph_block(istart+1:istart+nmvod)
          bufcolex(nmvod+1:nmvod+lencouplbl)=
     &        bufcol(1:lencouplbl)

!         do k=1,nmvod+lencouplbl
!           write(*,*) 'coupling out:',k,mjcol,bufcolex(k)
!         enddo

!
! write out full column.
!
          CALL XCOLWRIT(IOCH,wthr,bufcolex,wbufr,wbufi,wbufj,
     &                  IBFP,NBUFS,LENCOUPLBL+NMVOD,mjcol,INTBUF)
          mjcol = mjcol + 1

        enddo   !bfun
        enddo   !jfun

      enddo  ! jrep

      IF((mjcol-1).ne.mvo(desrep))
     &  stop 'error 3 in build_x_cousat_r!'
!________________________________________________________
!|
!|
!| from here the main and coupling block are written to
!| the external file. The SAT entries are also written to
!| the same file.
!| ATT**  we can release, and should release *all*
!| memory which we do not need anymore.
!| Due to the large size of the matrices we can not afford
!| a four-indexed lookup table as it was used in the two-particle
!| propagator!
!|
!|_______________________________________________________


      CALL XSATMAKE(IOCH,IBFP,NBUFS,INTBUF,DESREP,WTHR,
     &              WBUFR,WBUFI,WBUFJ,EPS,OOLT,VVLT,
     &              DOXEXT)
!
!  the SAT block in DESREP is now also on the file
!  and we can return in this symmetry and diagonalize...
!
!
!  flush buffer
!
      IF(IBFP.gt.0) THEN
        CALL XCOLFLSH(IOCH,IBFP,NBUFS,INTBUF,WBUFR,WBUFI,WBUFJ)
      ENDIF

      write(*,*) '--------------------------------------------'
      write(*,*) '---- Coupling and Satellite blocks finished.'
      write(*,*) '---- Number of required buffers:',NBUFS
      write(*,*) '--------------------------------------------'

      deallocate(bufcolex)
      deallocate(bufcol)
      deallocate(vovv1)
      deallocate(vooo1)
      deallocate(vvlt)
      deallocate(oolt)
!
! release storage for matrix writing
!
      deallocate(wbufj)
      deallocate(wbufi)
      deallocate(wbufr)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_X_COUSAT_C (PHPH_BLOCK,EPS,DESREP,LENCOUPLBL,
     &                             INTBUF,NBUFS,DOXEXT,WTHR,IOCH)
!
      IMPLICIT INTEGER (A-Z)
!     IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!  Computes (complex) 2H2P/HP and SAT block of the exc. ADC Matrix
!  In the corresponding symmetry DESREP.
!  The matrix is built columnwise. Desrep determines column symmetry
!  according to |b><j| spinors.
!
!---------------Calling variables--------------------------------------
!
      real*8                           :: PHPH_BLOCK(:)
      REAL*8                           :: EPS(*)
      INTEGER                          :: DESREP
      INTEGER                          :: LENCOUPLBL
      INTEGER                          :: INTBUF
      INTEGER                          :: NBUFS
      LOGICAL                          :: DOXEXT
      REAL*8                           :: WTHR
      INTEGER                          :: IOCH
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/files.inc"
#include  "../relccsd/param.inc"
!
!---------------Local variables--------------------------------------
!
      INTEGER              :: n1,irep,k
      INTEGER*8            :: i8,n8  !needed for the ivovvt variable!
      INTEGER              :: jrep,brep,jfun,bfun
      INTEGER              :: ialloc,mjcol,nmvod
      integer              :: n2

      INTEGER              :: mxno,mxnv,ioff
      INTEGER, dimension (:,:,:), allocatable  :: oolt
      INTEGER, dimension (:,:,:), allocatable  :: vvlt

      REAL*8, dimension(:), allocatable     :: vooo1
      REAL*8, dimension(:), allocatable     :: vovv1
      complex*16, dimension(:), allocatable :: bufcol,bufcolex
!
!  variables for writing out ADC matrix
!
      complex*16,dimension(:), allocatable   :: wbufr
      INTEGER, dimension(:), allocatable     :: wbufi, wbufj
      INTEGER              :: IBFP
      CHARACTER*6          :: NAMEX = 'ADCXPM'
!
!  auxiliary array for complex data type
!
      complex*16, allocatable, dimension(:)    :: phph_blc

      real*8,parameter        :: rzero=(0.0d0)
      complex*16,parameter    :: czero=(rzero,rzero)
!
!---------------Interface area --------------------------------------
!
      interface

        INTEGER FUNCTION XCOLDET(i1)
          INTEGER                         :: i1
        END FUNCTION XCOLDET

        SUBROUTINE XCOLMAKE_C(i1,i2,ca1,
     &                        i3,i4,i5,i6,
     &                        i7,i8,
     &                        ra2,ra3,ia1,ia2)
          INTEGER                         :: i1,i2
          complex*16, dimension(:)        :: ca1
          INTEGER                         :: i3,i4,i5,i6
          INTEGER                         :: i7,i8
          REAL*8, dimension(:)            :: ra2,ra3
          INTEGER, dimension(:,:,:)       :: ia1,ia2
        END SUBROUTINE XCOLMAKE_C

        SUBROUTINE XCOLWRIT_C(i1,r1,ca1,ca2,ia1,ia2,i2,i3,i4,i5,i6)
          INTEGER                         :: i1
          REAL*8                          :: r1
          complex*16, dimension(:)        :: ca1,ca2
          INTEGER, dimension(:)           :: ia1,ia2
          INTEGER                         :: i2,i3,i4,i5,i6
        END SUBROUTINE XCOLWRIT_C

        SUBROUTINE XSATMAKE_C (i1,i2,i3,i4,i5,r1,
     &              ca1,ia1,ia2,ra2,ia3,ia4,l1)
          INTEGER                     :: i1,i2,i3,i4,i5
          REAL*8                      :: r1
          complex*16, dimension(:)    :: ca1
          INTEGER, dimension(:)       :: ia1,ia2
          REAL*8, dimension(*)        :: ra2
          INTEGER,dimension(:,:,:)    :: ia3,ia4
          LOGICAL                     :: l1
        END SUBROUTINE XSATMAKE_C

        SUBROUTINE XCOLFLSH_C(i1,i2,i3,i4,ca1,ia1,ia2)
          INTEGER                         :: i1,i2,i3,i4
          complex*16, dimension(:)        :: ca1
          INTEGER, dimension(:)           :: ia1,ia2
        END SUBROUTINE XCOLFLSH_C

      end interface
!
!---------------Executable code--------------------------------------
!
      if(rcw.ne.2) then
        call quit('Internal error in BUILD_X_COUSAT_C !')
      endif

!___________________________________________________________
!|
!|  generate complex auxiliary array from original real one.

      n2 = jvovo(nrep+1)
      allocate(phph_blc(n2))
      phph_blc = czero
      do i=1,n2
        ix = 2*(i-1) + 1
        phph_blc(i)=dcmplx(phph_block(ix),phph_block(ix+1))
      enddo

!|  configuration generator
!|__________________________________________

!________________________________________________________
!|     allocate ADC MATRIX WRITE BUFFERS and counters
!|
!|
      allocate(wbufr(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 1 could not be allocated.'
      allocate(wbufi(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 2 could not be allocated.'
      allocate(wbufj(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 3 could not be allocated.'

      wbufr=czero
      wbufi = 0
      wbufj = 0
      IBFP = 0   ! pointer to actual buffer level
      NBUFS = 0  ! counter for full buffers
!|
!|
!|_______________________________________________________
!________________________________________________________
!|     create OOT and VVT lookup tables
!|     tables are small, computation is too clumsy
!|
      mxno = 0; mxnv = 0
      do irep=1,nrep  
        if(no(irep).gt.mxno) mxno=no(irep)
        if(nv(irep).gt.mxnv) mxnv=nv(irep)
      enddo
      allocate(oolt(mxno,mxno,nrep),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #a in BUILD_X_COUSAT_R'
      allocate(vvlt(mxnv,mxnv,nrep),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #b in BUILD_X_COUSAT_R'
      OOLT = 0
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NO(IREP)
          DO I=J+1,NO(IREP)
            OOLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
      VVLT = 0
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NV(IREP)
          DO I=J+1,NV(IREP)
            VVLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
!|
!|
!|_______________________________________________________
!
!
!  Attention*** Some arrays that have to store *huge* numbers 
!  Are defined as Integer*8 in symm.inc ! The generic sorters
!  operate with unspecified INTEGER declarations. Transferring such
!  an array as an argument to the SRT_xxx results in an erroneous
!  behaviour of the sorting routine.

      if(rcw.ne.2) then
        call quit('*** RCW ne 2! Internal error (build_x_cousat_c) ***')
      endif
!
! fetch <VO||OO> integrals needed for coupling block 
! 
      n1=ivooot(nrep+1)
      allocate(vooo1(n1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error vooo1!'
      call getvooo(vooo1)
!
! fetch <VO||VV> integrals needed for coupling block 
! 
      n8=ivovvt(nrep+1)
      allocate(vovv1(n8*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error vovv1!'
      call rdvovv (vovv1)
!
! length of main block:
!
      NMVOD = MVO(desrep)
!
! determine actual length of coupling block
! and allocate column buffer
!
      lencouplbl = XCOLDET(desrep)
      WRITE(IW,*) '   *)) Current symmetry:          ',desrep
      WRITE(IW,*) '   *)) Block length of main block:',NMVOD
      WRITE(IW,*) '   *)) Block length of coupling block:',lencouplbl
      WRITE(IW,'(A,I8,A)') '    *)) Matrix dimension: [',
     &        lencouplbl + NMVOD,' ]**2'

      ALLOCATE(BUFCOL(lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #3 in BUILD_X_COUSAT_R'
      ALLOCATE(BUFCOLEX(NMVOD + lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #4 in BUILD_X_COUSAT_R'
!
!  according to the SRT16 convention, 
!  the columns also have DESREP symmetry. 
!  We loop over the columns in the same manner and construct
!  the columns of the 2P2H/PH coupling block. 
!  after construction we write out the column together with the
!  components of the main block.
!
      mjcol = 1
      do jrep = 1,nrep

        brep = multb(desrep+nrep,jrep,2)

        do jfun = 1,NO(jrep)
        do bfun = 1,NV(brep)

          bufcol   = czero
          bufcolex = czero

          CALL XCOLMAKE_C(desrep,lencouplbl,bufcol,
     &                    jrep,brep,jfun,bfun,
     &                    mxno,mxnv,
     &                    vooo1,vovv1,oolt,vvlt)
!
! combine main block column with coupling block column
! istart is the start address for picking the right main block
! data to corresponding symetry DESREP.
!
          istart = JVOVO(desrep) + (mjcol-1)*nmvod
          bufcolex(1:nmvod)=phph_blc(istart+1:istart+nmvod)
          bufcolex(nmvod+1:nmvod+lencouplbl)=
     &        bufcol(1:lencouplbl)
!
! write out full column.
!
          CALL XCOLWRIT_C(IOCH,wthr,bufcolex,wbufr,wbufi,wbufj,
     &                    IBFP,NBUFS,LENCOUPLBL+NMVOD,mjcol,INTBUF)
          mjcol = mjcol + 1
        enddo
        enddo
      enddo

      IF((mjcol-1).ne.mvo(desrep))
     &  stop 'error 3 in BUILD_X_COUSAT_C!'
!________________________________________________________
!|
!|
!| from here the main and coupling block are written to
!| the external file. The SAT entries are also written to
!| the same file.
!| ATT**  we can release, and should release *all*
!| memory which we do not need anymore.
!| Due to the large size of the matrices we can not afford
!| a four-indexed lookup table as it was used in the two-particle
!| propagator!
!|
!|_______________________________________________________


      CALL XSATMAKE_C(IOCH,IBFP,NBUFS,INTBUF,DESREP,WTHR,
     &                WBUFR,WBUFI,WBUFJ,EPS,OOLT,VVLT,
     &                DOXEXT)
!
!  the SAT block in DESREP is now also on the file
!  and we can return in this symmetry and diagonalize...
!
!
!  flush buffer
!
      IF(IBFP.gt.0) THEN
        CALL XCOLFLSH_C(IOCH,IBFP,NBUFS,INTBUF,WBUFR,WBUFI,WBUFJ)
      ENDIF

      write(*,*) '--------------------------------------------'
      write(*,*) '---- Coupling and Satellite blocks finished.'
      write(*,*) '---- Number of required buffers:',NBUFS
      write(*,*) '--------------------------------------------'

      deallocate(bufcolex)
      deallocate(bufcol)
      deallocate(vovv1)
      deallocate(vooo1)
      deallocate(vvlt)
      deallocate(oolt)
!
! release storage for matrix writing
!
      deallocate(wbufj)
      deallocate(wbufi)
      deallocate(wbufr)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION XCOLDET(desrep)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!  determines length of excitation coupling block (V>V O>O: desrep)
!  att. Length is not related to entries in ivvoott(:) !
!
!---------------Calling variables--------------------------------------
!
      INTEGER                       :: XCOLDET
      INTEGER, INTENT(IN)           :: desrep
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/symm.inc"
!
!---------------Executable code--------------------------------------
!
      ICNT = 0

      DO KLREP = 1,NREP     ! klrep, cdrep and desrep are bosonic
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                  ICNT = ICNT + 1
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
 
      XCOLDET = ICNT
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLMAKE (desrep,lencouplbl,bufcol,
     &                     jrep,brep,J,B,
     &                     mxno,mxnv,
     &                     voooz,vovvz,
     &                     oolt,vvlt)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
! Build the coupling block column for which we provide
! jrep,brep,J and B. BUFCOL was allocated by the caller
! and is filled in this subroutine.
!              *** real case ***
! all arrays are transferred as ** assumed shape **
!
      INTEGER                          :: desrep,lencouplbl
      REAL*8, dimension(:)             :: bufcol
      INTEGER                          :: jrep,brep,J,B
      INTEGER                          :: mxno,mxnv
      REAL*8, dimension(:)             :: voooz,vovvz
      INTEGER, dimension(:,:,:)        :: oolt,vvlt
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Executable code--------------------------------------
!

!     write(iw,*) 'Xcolmake:',desrep,lencouplbl,jrep,brep,j,b
      CDKL = 0
!
!  do the |>|> side...
!
      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      if(MULTB(CDREP+NREP,KLREP+NREP,2).ne.desrep) then
         Stop 'xcolmake error #aa1'
      endif
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
!
!  do the <|<| side...
!
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                   CDKL = CDKL + 1
!
! from here we have all ireps/spinor numbers.
! Now treat the four delta cases.
! CDKL is the local column index in the current irrep!
! The offsets into the integral fields are *calculated*
! since a lookup table is impossible.
!____________________________________________________________
!|
!|
!|  Case A:  -< dj||kl > delta_bc, we calc a VOOO-offset.
!|
      if(BREP.eq.CREP.and.B.eq.C) then
        DJREP = MULTB(DREP,JREP,1)
        IF(MULTB(DJREP+NREP,KLREP+NREP,2).ne.1) STOP 'IE #44'
        IF(KREP.EQ.LREP) THEN
          IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
        ELSE
          IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
        ENDIF
        IOFF=IOFF * NVO(KLREP) + IVOOOT(KLREP)
        IOFF=IOFF + IIVO(DREP,JREP) + (J-1)*NV(DREP) + D
        if(ioff.gt.ivooot(nrep+1).or.ioff.lt.1) 
     &     stop 'intern. error #66!'
        bufcol(cdkl) = bufcol(cdkl) - voooz(ioff)
      endif
!|
!|  Case B:  +< cj||kl > delta_bd, we calc a VOOO-offset.
!|
      if(BREP.eq.DREP.and.B.eq.D) then
        CJREP = MULTB(CREP,JREP,1)
        IF(MULTB(CJREP+NREP,KLREP+NREP,2).ne.1) STOP 'IE #45'
        IF(KREP.EQ.LREP) THEN
          IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
        ELSE
          IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
        ENDIF
        IOFF=IOFF * NVO(KLREP) + IVOOOT(KLREP)
        IOFF=IOFF + IIVO(CREP,JREP) + (J-1)*NV(CREP) + C
        if(ioff.gt.ivooot(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #67!'
        bufcol(cdkl) = bufcol(cdkl) + voooz(ioff)
      endif
!|
!|  Case C:  -< bl||cd >* delta_jk, we calc a VOVV offset
!|
      if(JREP.eq.KREP.and.J.eq.K) then
        BLREP = MULTB(BREP,LREP,1)
        IF(MULTB(BLREP+NREP,CDREP+NREP,2).ne.1) STOP 'IE #46'
        IF(CREP.EQ.DREP) THEN
          IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
        ELSE
          IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
        ENDIF
        IOFF=IOFF * NVO(CDREP) + IVOVVT(CDREP)
        IOFF=IOFF + IIVO(BREP,LREP) + (L-1)*NV(BREP) + B
        if(ioff.gt.ivovvt(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #68!'
        bufcol(cdkl) = bufcol(cdkl) - vovvz(ioff)
      endif
!|
!|  Case D:  +< bk||cd >* delta_jl, we calc a VOVV offset
!|
      if(JREP.eq.LREP.and.J.eq.L) then
        BKREP = MULTB(BREP,KREP,1)
        IF(MULTB(BKREP+NREP,CDREP+NREP,2).ne.1) STOP 'IE #47'
        IF(CREP.EQ.DREP) THEN
          IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
        ELSE
          IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
        ENDIF
        IOFF=IOFF * NVO(CDREP) + IVOVVT(CDREP)
        IOFF=IOFF + IIVO(BREP,KREP) + (K-1)*NV(BREP) + B
        if(ioff.gt.ivovvt(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #69!'
        bufcol(cdkl) = bufcol(cdkl) + vovvz(ioff)
      endif
!|
!|            END   Case A-D
!|___________________________________________________________
!
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      IF(CDKL.ne.LENCOUPLBL) STOP 'Internal error #77!'

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLMAKE_C (desrep,lencouplbl,bufcol,
     &                       jrep,brep,J,B,
     &                       mxno,mxnv,
     &                       voooz,vovvz,oolt,vvlt)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
! Build the coupling block column for which we provide
! jrep,brep,J and B. BUFCOL was allocated by the caller
! and is filled in this subroutine.
!              *** real case ***
! all arrays are transferred as ** assumed shape **
! the size in each dimension stems from the caller
!
      INTEGER                          :: desrep
      INTEGER                          :: lencouplbl
      complex*16, dimension(:)         :: bufcol
      INTEGER                          :: jrep,brep,J,B
      INTEGER                          :: mxno,mxnv
      REAL*8, dimension(:)             :: voooz
      REAL*8, dimension(:)             :: vovvz
      INTEGER, dimension(:,:,:)        :: oolt
      INTEGER, dimension(:,:,:)        :: vvlt
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
!
!---------------Executable code--------------------------------------
!
!
! check consistency
!
      if(rcw.ne.2) then
        call quit('Internal error in XCOLMAKE_C !')
      endif

! start

      CDKL = 0
!
!  do the |>|> side...
!
      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
!
!  do the <|<| side...
!
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                   CDKL = CDKL + 1
!
! from here we have all ireps and spinor numbers for filling the
! corresponding matrix element at the correct position.
! we have to abandon indentation here in order to fit some code
! on the line !
!
! treat the four delta cases.
! CDKL can be seen as the local column index here in this irrep!
! the offsets are *calculated* and not yielded via a lookup table
! again.
!____________________________________________________________
!|
!|
!|
!|  Case A:  < dj||kl > delta_bc, we calc a VOOO-offset.
!|
      if(BREP.eq.CREP.and.B.eq.C) then
        DJREP = MULTB(DREP,JREP,1)
        IF(MULTB(DJREP+NREP,KLREP+NREP,2).ne.1) STOP 'IE #44'
        IF(KREP.EQ.LREP) THEN
          IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
        ELSE
          IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
        ENDIF
        IOFF=IOFF * NVO(KLREP) + IVOOOT(KLREP)
        IOFF=IOFF + IIVO(DREP,JREP) + (J-1)*NV(DREP) + D
        if(ioff.gt.ivooot(nrep+1).or.ioff.lt.1) 
     &     stop 'intern. error #66!'
        iir = (ioff - 1)*rcw + 1
        bufcol(cdkl) = bufcol(cdkl) - dcmplx(voooz(iir),voooz(iir+1))
      endif
!|
!|  Case B:  < cj||kl > delta_bd, we calc a VOOO-offset.
!|
      if(BREP.eq.DREP.and.B.eq.D) then
        CJREP = MULTB(CREP,JREP,1)
        IF(MULTB(CJREP+NREP,KLREP+NREP,2).ne.1) STOP 'IE #45'
        IF(KREP.EQ.LREP) THEN
          IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
        ELSE
          IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
        ENDIF
        IOFF=IOFF * NVO(KLREP) + IVOOOT(KLREP)
        IOFF=IOFF + IIVO(CREP,JREP) + (J-1)*NV(CREP) + C
        if(ioff.gt.ivooot(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #67!'
        iir = (ioff - 1)*rcw + 1
        bufcol(cdkl) = bufcol(cdkl) + dcmplx(voooz(iir),voooz(iir+1))
      endif
!|
!|  Case C:  < bl||cd >* delta_jk, we calc a VOVV offset
!|
      if(JREP.eq.KREP.and.J.eq.K) then
        BLREP = MULTB(BREP,LREP,1)
        IF(MULTB(BLREP+NREP,CDREP+NREP,2).ne.1) STOP 'IE #46'
        IF(CREP.EQ.DREP) THEN
          IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
        ELSE
          IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
        ENDIF
        IOFF=IOFF * NVO(CDREP) + IVOVVT(CDREP)
        IOFF=IOFF + IIVO(BREP,LREP) + (L-1)*NV(BREP) + B
        if(ioff.gt.ivovvt(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #68!'
        iir = (ioff - 1)*rcw + 1
        bufcol(cdkl) = bufcol(cdkl) - 
     &    conjg(dcmplx(vovvz(iir),vovvz(iir+1)))
      endif
!|
!|  Case D:  < bk||cd >* delta_jl, we calc a VOVV offset
!|
      if(JREP.eq.LREP.and.J.eq.L) then
        BKREP = MULTB(BREP,KREP,1)
        IF(MULTB(BKREP+NREP,CDREP+NREP,2).ne.1) STOP 'IE #47'
        IF(CREP.EQ.DREP) THEN
          IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
        ELSE
          IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
        ENDIF
        IOFF=IOFF * NVO(CDREP) + IVOVVT(CDREP)
        IOFF=IOFF + IIVO(BREP,KREP) + (K-1)*NV(BREP) + B
        if(ioff.gt.ivovvt(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #69!'
        iir = (ioff - 1)*rcw + 1
        bufcol(cdkl) = bufcol(cdkl) + 
     &    conjg(dcmplx(vovvz(iir),vovvz(iir+1)))
      endif
!|
!|
!|___________________________________________________________
!
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      IF(CDKL.ne.LENCOUPLBL) STOP 'Internal error #77!'
!
!     write(iw,*) 'Xcolmake: current run finished.'
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLWRIT (IOCH,wthr,bufferin,buf,ioi,ioj,
     &                     IB,NBUFS,LENBUFFER,mjcol,INTBUF)
!
      IMPLICIT NONE

      integer                       ::  ioch
      real*8                        ::  wthr
      real*8, dimension(:)          ::  bufferin,buf
      integer, dimension(:)         ::  ioi,ioj
      integer                       ::  ib,nbufs,lenbuffer,mjcol,intbuf
!
!---------------Description--------------------------------------------
!
!  writes out a complete column of the ADC matrix held in bufferin
!  to the buckets provided by the caller. If the buckets are full,
!  they are emptied to the file in IOCH.
!
!---------------Local variables ---------------------------------------
!
      INTEGER            ::  i,jdummy,ixx
      REAL*8             ::  s
!
!---------------Execution ---------------------------------------------
!
      jdummy = 0
      DO i=mjcol,lenbuffer  !i hereby serves as the row index !
        s=bufferin(i)
        IF(abs(s).gt.wthr) THEN
          IB=IB+1
          BUF(IB)=S
          IOI(IB)=i
          IOJ(IB)=mjcol
!         write(*,*) 'coupling out:',i,mjcol,S
          IF(IB.EQ.INTBUF) THEN
            NBUFS = NBUFS + 1
            WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &                  (IOI(IXX),IXX=1,INTBUF),
     &                  (IOJ(IXX),IXX=1,INTBUF),
     &                   INTBUF,JDUMMY
            IB = 0
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLWRIT_C (ioch,wthr,bufferin,buf,ioi,ioj,
     &                       ib,nbufs,lenbuffer,mjcol,intbuf)
!
      IMPLICIT NONE

      integer                       ::  ioch
      real*8                        ::  wthr
      complex*16, dimension(:)      ::  bufferin,buf
      integer, dimension(:)         ::  ioi,ioj
      integer                       ::  ib,nbufs,lenbuffer,mjcol,intbuf
!
!---------------Description--------------------------------------------
!
!  writes out a complete column of the ADC matrix held in bufferin
!  to the buckets provided by the caller. If the buckets are full,
!  they are emptied to the file in IOCH.
!
!---------------Local variables ---------------------------------------
!
      INTEGER            ::  i,jdummy,ixx
      complex*16         ::  s
!
!---------------Execution ---------------------------------------------
!
      jdummy = 0
      DO i=mjcol,lenbuffer  !i hereby serves as the row index starting from mjcol 
        s=bufferin(i)
        IF(abs(s).gt.wthr) THEN
          IB=IB+1
          BUF(IB)=S
          IOI(IB)=i
          IOJ(IB)=mjcol
          IF(IB.EQ.INTBUF) THEN
            NBUFS = NBUFS + 1
            WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &                  (IOI(IXX),IXX=1,INTBUF),
     &                  (IOJ(IXX),IXX=1,INTBUF),
     &                   INTBUF,JDUMMY
            IB = 0
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLFLSH (IOCH,IB,NBUFS,INTBUF,buf,ioi,ioj)
!
      IMPLICIT NONE

      integer                       :: IOCH,IB,NBUFS,INTBUF
      real*8, dimension(:)          :: buf
      integer, dimension(:)         :: ioi,ioj

      integer                       :: ixx,jdummy = 0

      NBUFS = NBUFS + 1
      WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &            (IOI(IXX),IXX=1,INTBUF),
     &            (IOJ(IXX),IXX=1,INTBUF),
     &            IB,JDUMMY
      IB = 0

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLFLSH_C (IOCH,IB,NBUFS,INTBUF,buf,ioi,ioj)
!
      IMPLICIT NONE

      integer                       :: IOCH,IB,NBUFS,INTBUF
      complex*16, dimension(:)      :: buf
      integer, dimension(:)         :: ioi,ioj

      integer                       :: ixx,jdummy = 0

      NBUFS = NBUFS + 1
      WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &            (IOI(IXX),IXX=1,INTBUF),
     &            (IOJ(IXX),IXX=1,INTBUF),
     &            IB,JDUMMY
      IB = 0

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XSATMAKE(ioch,ibfp,nbufs,intbuf,desrep,wthr,
     &                    buf,ioi,ioj,eps,oolt,vvlt,dox)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!  generates the ADC-2x SAT block.
!  the column counter MJCOL carries the actual value at this stage!
!
!
!---------------Calling variables--------------------------------------
!
      INTEGER                   :: ioch,ibfp,nbufs,intbuf
      REAL*8                    :: wthr
      INTEGER                   :: mjcol,desrep
      REAL*8, dimension(:)      :: buf
      INTEGER, dimension(:)     :: ioi,ioj
      REAL*8, dimension(*)      :: eps
      INTEGER,dimension(:,:,:)  :: oolt, vvlt
      LOGICAL                   :: dox  !tells the routine ADC2-x
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8, allocatable, dimension(:)    :: vvvvz, ooooz, vovoz
      INTEGER                              :: n1,ialloc
      INTEGER*8                            :: i8,n8
      REAL*8                               :: EA,EB,EI,EJ
      REAL*8                               :: S
      INTEGER                              :: JDUMMY
!
!---------------Interface area --------------------------------------
!
      interface

        SUBROUTINE mpi_vvvv_complete(ra1,i1)
          real*8, dimension(:)            :: ra1
          integer*8                       :: i1
        END SUBROUTINE mpi_vvvv_complete

      end interface
!
!---------------Start execution--------------------------------------
!

! allocate space and read in integrals

      n8 = ivvvvtt(nrep+1)
      allocate(vvvvz(n8),stat=ialloc)
      if(ialloc.ne.0) STOP 'VVVV buffer could not be allocated.'

#if defined (VAR_MPI)
      do i8=1,n8
        vvvvz(i8) = -9999.0d0
      enddo
      call rdvvvv(vvvvz)  ! master reads its own chunks
      call mpi_vvvv_complete(vvvvz,n8)  !.... and completes the stream.
      write(*,*) 'VVVV check...'
      do i8=1,n8
        if(vvvvz(i8).eq. -9999.0d0)
     &    call quit('Internal error: gap in the VVVV stream!')
      enddo
      write(*,*) '... complete.'
#else
      call rdvvvv(vvvvz)  ! master reads complete stream in serial mode
#endif

      n1 = ioooott(nrep+1)
      allocate(ooooz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'OOOO buffer could not be allocated.'

      n1 = ivovo(nrep+1)
      allocate(vovoz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOVO buffer could not be allocated.'

      CALL GETOOOO(ooooz)
      CALL GETVOVO(vovoz)
!
! SAT elements are calculated only for the lower triangle which
! means that ABIJ >= CDKL !
! no lookup table (too big!)
!
      CDKL=0    !column counter

      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                   CDKL = CDKL + 1
!_______________________________________________________________
!
!  here we are at the innermost position of the cdkl loop.
!  we abandon indentation here and start the abij loop.
!  This part of code unfortunately scales as N**8 but we have
!  no remedy at the moment.

      ABIJ=0     !row counter

      DO IJREP = 1,NREP
      ABREP = MULTB(DESREP+NREP,IJREP+NREP,2)
      DO 44 JREP = 1, NREP
      IREP = MULTB(JREP,IJREP+NREP,2)
      IF (IREP.LT.JREP) GOTO 44
      DO J = 1, NO(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, NO(IREP)
            DO 33 BREP = 1, NREP
              AREP = MULTB(BREP,ABREP+NREP,2)
              IF (AREP.LT.BREP) GOTO 33
              DO B = 1, NV(BREP)
                AMIN = 1
                IF (AREP.EQ.BREP) AMIN = B + 1
                DO A = AMIN, NV(AREP)
                   ABIJ = ABIJ + 1
!_______________________________________________________________
!
!  here we are at the innermost position of the abij loop.
!  Again, we abandon indentation here and start the delta cases
!
      IF(ABIJ.GE.CDKL) THEN
!
! we are on or below the diagonal and start selection
!
! otherwise: DO NOTHING!
!
!  a: check diagonal
!
       S = 0.0d0

       if(arep.eq.crep.and.brep.eq.drep.and.
     &    irep.eq.krep.and.jrep.eq.lrep) then
         if(a.eq.c.and.b.eq.d.and.i.eq.k.and.j.eq.l) then
           aoff=io(nrep+1) + iv(arep) + a
           boff=io(nrep+1) + iv(brep) + b
           ioff=io(irep) + i
           joff=io(jrep) + j
           ea=eps(aoff); eb=eps(boff); ei=eps(ioff); ej=eps(joff)
           S = S + (ea + eb - ei - ej)
!          write(*,*) 'qaq diagonal:',S,ABIJ,CDKL
         endif
       endif
!
!  if strict adc2, we only check diagonal in the SAT block
!  and leave the rest.
!
       IF(DOX) THEN    !yes, we need more ... (ADC2X)
!
!
!  b: check VVVV contributions  <AB||CD> delta_ik delta_jl
!
       if(irep.eq.krep.and.jrep.eq.lrep) then
         if(i.eq.k.and.j.eq.l) then
           ZCDREP=MULTB(CREP,DREP,1)
           IF(ZCDREP.ne.CDREP) STOP 'internal error 1 XSATMAKE'
           IF(CREP.EQ.DREP) THEN
             IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
           ELSE
             IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
           ENDIF
           IOFF=IOFF * NVVT(ZCDREP) + IVVVVTT(ZCDREP)
           IF(AREP.EQ.BREP) THEN
             IOFF=IOFF + IIVVT(AREP,BREP) + VVLT(A,B,AREP)
           ELSE
             IOFF=IOFF + IIVVT(AREP,BREP) +
     &         (B-1)*NV(AREP) + A
           ENDIF
           if(ioff.gt.ivvvvtt(nrep+1)) stop 'vvvv offset error !'
           S = S + vvvvz(IOFF)
         endif
       endif
!
!  c: check OOOO contributions  <IJ||KL>* delta_ac delta_bd
!
       if(arep.eq.crep.and.brep.eq.drep) then
         if(a.eq.c.and.b.eq.d) then
           ZKLREP=MULTB(KREP,LREP,1)
           IF(ZKLREP.ne.KLREP) STOP 'internal error 2 XSATMAKE'
           IF(KREP.EQ.LREP) THEN
             IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
           ELSE
             IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
           ENDIF
           IOFF=IOFF * NOOT(ZKLREP) + IOOOOTT(ZKLREP)
           IF(IREP.EQ.JREP) THEN
             IOFF=IOFF + IIOOT(IREP,JREP) + OOLT(I,J,IREP)
           ELSE
             IOFF=IOFF + IIOOT(IREP,JREP) +
     &         (J-1)*NO(IREP) + I
           ENDIF
           if(ioff.gt.ioooott(nrep+1)) stop 'oooo offset error !'
           S = S + ooooz(IOFF)
         endif
       endif
!
! Block 1: check the four delta_bd  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(brep.eq.drep.and.b.eq.d)then

! 1a
         if(jrep.eq.lrep.and.j.eq.l) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 1b
         if(irep.eq.krep.and.i.eq.k) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 1c
         if(jrep.eq.krep.and.j.eq.k) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 1d
         if(irep.eq.lrep.and.i.eq.l) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 

       endif  ! delta_bd
!
! Block 2: check the four delta_ac  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(arep.eq.crep.and.a.eq.c)then
! 2a
         if(jrep.eq.lrep.and.j.eq.l) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 2b
         if(irep.eq.krep.and.i.eq.k) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 2c
         if(jrep.eq.krep.and.j.eq.k) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 2d
         if(irep.eq.lrep.and.i.eq.l) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 

       endif  !delta_ac
!
! Block 3: check the four delta_bc  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(brep.eq.crep.and.b.eq.c)then
! 3a
         if(jrep.eq.lrep.and.j.eq.l) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 3b
         if(irep.eq.krep.and.i.eq.k) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 3c
         if(jrep.eq.krep.and.j.eq.k) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 3d
         if(irep.eq.lrep.and.i.eq.l) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 

       endif  !delta_bc
!
! Block 4: check the four delta_ad  contributions
!
       if(arep.eq.drep.and.a.eq.d)then
! 4a
         if(jrep.eq.lrep.and.j.eq.l) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 4b
         if(irep.eq.krep.and.i.eq.k) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 4c
         if(jrep.eq.krep.and.j.eq.k) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 4d
         if(irep.eq.lrep.and.i.eq.l) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 

       endif  !delta_ad

       ENDIF     ! DOX
!
!  write out corresponding matrix entry with row and column number
!  note: local row and column counters both start at value 1.
!  row and column counter has to be offset by size of main block
!
       IF(ABS(S).gt.WTHR) THEN
         IBFP=IBFP+1
         BUF(IBFP)=S
         IOI(IBFP)=ABIJ + MVO(DESREP)
         IOJ(IBFP)=CDKL + MVO(DESREP)
!        write(*,*) 'Satmatout:',IOI(IBFP),IOJ(IBFP),S
         IF(IBFP.EQ.INTBUF) THEN
           NBUFS = NBUFS + 1
           WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &                 (IOI(IXX),IXX=1,INTBUF),
     &                 (IOJ(IXX),IXX=1,INTBUF),
     &                  INTBUF,JDUMMY
           IBFP = 0
         ENDIF
       ENDIF

      ENDIF     !(ABIJ.GE.CDKL)
!  
!
!_______________________________________________________________
                ENDDO
              ENDDO
 33         CONTINUE
         ENDDO
      ENDDO
 44   CONTINUE
      ENDDO
!_______________________________________________________________
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
!

      deallocate(vovoz)
      deallocate(ooooz)
      deallocate(vvvvz)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XSATMAKE_C (ioch,ibfp,nbufs,intbuf,desrep,wthr,
     &                       buf,ioi,ioj,eps,oolt,vvlt,dox)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!  generates the complex ADC-2x SAT block.
!  the column counter MJCOL carries the actual value at this stage!
!
!
!---------------Calling variables--------------------------------------
!
      INTEGER                    :: ioch,ibfp,nbufs,intbuf
      REAL*8                     :: wthr
      INTEGER                    :: mjcol,desrep
      complex*16, dimension(:)   :: buf
      INTEGER, dimension(:)      :: ioi,ioj
      REAL*8, dimension(*)       :: eps
      INTEGER,dimension(:,:,:)   :: oolt, vvlt
      LOGICAL                    :: dox  !tells the routine ADC2-x
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
!  att: the integral arrays (not the integrals!) remain real due to legacy code.
!  in F90 we have better type safety and need to convert to complex type explicitly
!
      REAL*8, allocatable, dimension(:)    :: vvvvz, ooooz, vovoz
      INTEGER                              :: n1,ialloc
      INTEGER*8                            :: i8,n8
      REAL*8                               :: EA,EB,EI,EJ
      INTEGER                              :: JDUMMY

      complex*16                           :: S
!
!---------------Interface area --------------------------------------
!
      INTERFACE

        SUBROUTINE mpi_vvvv_complete(ra1,i1)
          real*8, dimension(:)            :: ra1
          integer*8                       :: i1
        END SUBROUTINE mpi_vvvv_complete

      END INTERFACE
!
!---------------Start execution--------------------------------------
!
! check consistency
!
      if(rcw.ne.2) then
        call quit('Internal error in XSATMAKE_C !')
      endif

! allocate space and read in integrals

      n8 = ivvvvtt(nrep+1)*rcw
      allocate(vvvvz(n8),stat=ialloc)
      if(ialloc.ne.0) STOP 'VVVV buffer could not be allocated.'

#if defined (VAR_MPI)
      do i8=1,n8
        vvvvz(i8) = -9999.0d0
      enddo
      call rdvvvv(vvvvz)  ! master reads its own chunks
      call mpi_vvvv_complete(vvvvz,n8)  !.... and completes the stream.
      write(*,*) 'VVVV check...'
      do i8=1,n8
        if(vvvvz(i8).eq. -9999.0d0)
     &    call quit('Internal error: gap in the VVVV stream!')
      enddo
      write(*,*) '... complete.'
#else
      call rdvvvv(vvvvz)  ! master reads complete stream in serial mode
#endif

      n1 = ioooott(nrep+1)*rcw
      allocate(ooooz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'OOOO buffer could not be allocated.'

      n1 = ivovo(nrep+1)*rcw
      allocate(vovoz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOVO buffer could not be allocated.'

      CALL GETOOOO(ooooz)
      CALL GETVOVO(vovoz)
!________________________________________________________________________
!|
!|
!|               START     MAIN      LOOP
!|
!| SAT elements are calculated only for the lower triangle which
!| means that ABIJ >= CDKL !
!| no lookup table (too big!)
!|
!|
!|
      CDKL=0    !column counter

      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                   CDKL = CDKL + 1
!_______________________________________________________________
!
!  here we are at the innermost position of the cdkl loop.
!  we abandon indentation here and start the abij loop.
!  This part of code unfortunately scales as N**8 but we have
!  no remedy at the moment.

      ABIJ=0     !row counter

      DO IJREP = 1,NREP
      ABREP = MULTB(DESREP+NREP,IJREP+NREP,2)
      DO 44 JREP = 1, NREP
      IREP = MULTB(JREP,IJREP+NREP,2)
      IF (IREP.LT.JREP) GOTO 44
      DO J = 1, NO(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, NO(IREP)
            DO 33 BREP = 1, NREP
              AREP = MULTB(BREP,ABREP+NREP,2)
              IF (AREP.LT.BREP) GOTO 33
              DO B = 1, NV(BREP)
                AMIN = 1
                IF (AREP.EQ.BREP) AMIN = B + 1
                DO A = AMIN, NV(AREP)
                   ABIJ = ABIJ + 1
!_______________________________________________________________
!
!  here we are at the innermost position of the abij loop.
!  Again, we abandon indentation here and start the delta cases
!
      IF(ABIJ.GE.CDKL) THEN
!
! we are on or below the diagonal and start selection
!
! otherwise: DO NOTHING!
!
!  a: check diagonal, that's actually the K-contribution.
!
       S = (0.0d0,0.0d0)

       if(arep.eq.crep.and.brep.eq.drep.and.
     &    irep.eq.krep.and.jrep.eq.lrep) then
         if(a.eq.c.and.b.eq.d.and.i.eq.k.and.j.eq.l) then
           aoff=io(nrep+1) + iv(arep) + a
           boff=io(nrep+1) + iv(brep) + b
           ioff=io(irep) + i
           joff=io(jrep) + j
           ea=eps(aoff); eb=eps(boff); ei=eps(ioff); ej=eps(joff)
           S = S + dcmplx((ea + eb - ei - ej),0.0d0)
         endif
       endif
!
!  if strict adc2, we only check diagonal in the SAT block
!  and leave the rest.
!
       IF(DOX) THEN    !yes, we need more ... (ADC2X)
!
!
!  b: check VVVV contributions  <AB||CD> delta_ik delta_jl
!
       if(irep.eq.krep.and.jrep.eq.lrep) then
         if(i.eq.k.and.j.eq.l) then
           ZCDREP=MULTB(CREP,DREP,1)
           IF(ZCDREP.ne.CDREP) STOP 'internal error 1 XSATMAKE'
           IF(CREP.EQ.DREP) THEN
             IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
           ELSE
             IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
           ENDIF
           IOFF=IOFF * NVVT(ZCDREP) + IVVVVTT(ZCDREP)
           IF(AREP.EQ.BREP) THEN
             IOFF=IOFF + IIVVT(AREP,BREP) + VVLT(A,B,AREP)
           ELSE
             IOFF=IOFF + IIVVT(AREP,BREP) +
     &         (B-1)*NV(AREP) + A
           ENDIF
           if(ioff.gt.ivvvvtt(nrep+1)) stop 'vvvv offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vvvvz(iir),vvvvz(iir+1))
         endif
       endif
!
!  c: check OOOO contributions  <IJ||KL>* delta_ac delta_bd
!
       if(arep.eq.crep.and.brep.eq.drep) then
         if(a.eq.c.and.b.eq.d) then
           ZKLREP=MULTB(KREP,LREP,1)
           IF(ZKLREP.ne.KLREP) STOP 'internal error 2 XSATMAKE'
           IF(KREP.EQ.LREP) THEN
             IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
           ELSE
             IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
           ENDIF
           IOFF=IOFF * NOOT(ZKLREP) + IOOOOTT(ZKLREP)
           IF(IREP.EQ.JREP) THEN
             IOFF=IOFF + IIOOT(IREP,JREP) + OOLT(I,J,IREP)
           ELSE
             IOFF=IOFF + IIOOT(IREP,JREP) +
     &         (J-1)*NO(IREP) + I
           ENDIF
           if(ioff.gt.ioooott(nrep+1)) stop 'oooo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + conjg(dcmplx(ooooz(iir),ooooz(iir+1)))
         endif
       endif
!
! Block 1: check the four delta_bd  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(brep.eq.drep.and.b.eq.d)then

! 1a
         if(jrep.eq.lrep.and.j.eq.l) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 1b
         if(irep.eq.krep.and.i.eq.k) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 1c
         if(jrep.eq.krep.and.j.eq.k) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 1d
         if(irep.eq.lrep.and.i.eq.l) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 

       endif  ! delta_bd
!
! Block 2: check the four delta_ac  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(arep.eq.crep.and.a.eq.c)then
! 2a
         if(jrep.eq.lrep.and.j.eq.l) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 2b
         if(irep.eq.krep.and.i.eq.k) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 2c
         if(jrep.eq.krep.and.j.eq.k) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 2d
         if(irep.eq.lrep.and.i.eq.l) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 

       endif  !delta_ac
!
! Block 3: check the four delta_bc  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(brep.eq.crep.and.b.eq.c)then
! 3a
         if(jrep.eq.lrep.and.j.eq.l) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 3b
         if(irep.eq.krep.and.i.eq.k) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 3c
         if(jrep.eq.krep.and.j.eq.k) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 3d
         if(irep.eq.lrep.and.i.eq.l) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 

       endif  !delta_bc
!
! Block 4: check the four delta_ad  contributions
!
       if(arep.eq.drep.and.a.eq.d)then
! 4a
         if(jrep.eq.lrep.and.j.eq.l) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 4b
         if(irep.eq.krep.and.i.eq.k) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 4c
         if(jrep.eq.krep.and.j.eq.k) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 4d
         if(irep.eq.lrep.and.i.eq.l) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 

       endif  !delta_ad

       ENDIF     ! DOX
!
!  write out corresponding matrix entry with row and column number
!  note: local row and column counters both start at value 1.
!  row and column counter has to be offset by size of main block
!
       IF(ABS(S).gt.WTHR) THEN
         IBFP=IBFP+1
         BUF(IBFP)=S
         IOI(IBFP)=ABIJ + MVO(DESREP)
         IOJ(IBFP)=CDKL + MVO(DESREP)
!        write(*,*) 'Satmatout:',IOI(IBFP),IOJ(IBFP),S
         IF(IBFP.EQ.INTBUF) THEN
           NBUFS = NBUFS + 1
           WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &                 (IOI(IXX),IXX=1,INTBUF),
     &                 (IOJ(IXX),IXX=1,INTBUF),
     &                  INTBUF,JDUMMY
           IBFP = 0
         ENDIF
       ENDIF

      ENDIF     !(ABIJ.GE.CDKL)
!  
!
!_______________________________________________________________
                ENDDO
              ENDDO
 33         CONTINUE
         ENDDO
      ENDDO
 44   CONTINUE
      ENDDO
!_______________________________________________________________
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
!

      deallocate(vovoz)
      deallocate(ooooz)
      deallocate(vvvvz)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE print_adc_matrix(iw,iobase,intbuf,nbufs)
!
!     For small matrices only: print the whole matrix.

      IMPLICIT NONE

!-------------- formal parameters

      integer        ::  iw
      integer        ::  iobase
      integer        ::  intbuf
      integer        ::  nbufs

!-------------- common variables

!-------------- local variables

      integer                             :: i,ixx,nact,jdummy
      real*8, allocatable, dimension (:)  :: buf
      integer, allocatable, dimension (:) :: ioi,ioj
      integer                             :: icol,irow
      real*8                              :: a

!-------------- execution

      call PST('Printing the complete ADC matrix+')
      write(iw,*)
      write(iw,*)
!
!  allocate read/write buffers
!
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))

      do i=1,nbufs
        read(iobase,ERR=889) (buf(ixx),ixx=1,intbuf),
     &                       (ioi(ixx),ixx=1,intbuf),
     &                       (ioj(ixx),ixx=1,intbuf),
     &                       nact,jdummy
          do ixx = 1,nact
            irow = ioi(ixx)
            icol = ioj(ixx)
            a    = buf(ixx)
            write(iw,'(2I8,E22.14)') irow,icol,a
          enddo
      enddo
      write(iw,*)
      write(iw,*) '---------------------------------------------'

      deallocate(ioi,ioj,buf)
      return

 889  call quit('Error in reading ADC matrix')

      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function qmasymh(e,lde)
!
!     computes deviation from Hermiticity of a quadratic matrix
!
      implicit none
      real*8                       :: qmasymh
      complex*16, dimension(:,:)   :: e
      integer                      :: lde

! local variables

      real*8                       :: asym
      integer                      :: i,j
!
!---------------executable code--------------------------------------
!
      asym=0.0D0

      if(lde.eq.1) then
       qmasymh = asym
       return
      endif

      do i=2,lde
        do j=1,i-1
          asym = asym + abs(conjg(e(i,j)) - e(j,i))
        enddo
      enddo

      qmasymh = asym

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function qmasymh2(e,lde)
!
!     computes deviation from Hermiticity of a quadratic matrix
!     but here, no assumed-shape array argument is given.
!     Assumed-shape array arguments are possible in fortran90 in combination
!     with an explicit interface. However, we have a real array that is to be reinterprerted
!     as a quadratic complex array! We therefore have to give in addition the
!     ** leading dimension ** for a correct offset calculation by the compiler
!
      implicit none
      real*8                       :: qmasymh2
      integer                      :: lde
      complex*16                   :: e(lde,lde)

! local variables

      real*8                       :: asym
      integer                      :: i,j
!
!---------------executable code--------------------------------------
!
      asym=0.0D0

      if(lde.eq.1) then
       qmasymh2 = asym
       return
      endif

      do i=2,lde
        do j=1,i-1
!         write(*,*) '   i/j:',i,j,e(i,j),e(j,i)
          asym = asym + abs(conjg(e(i,j)) - e(j,i))
        enddo
      enddo

      qmasymh2 = asym

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE xmvmul_sp_r(iobase,intbuf,ladc,PO,PN)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      INTEGER                        :: iobase,intbuf,ladc
      REAL*8                         :: PO(LADC), PN(LADC)
!
!------------------------ local variables -------------------------
!
      REAL*8,  allocatable, dimension(:)      :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      REAL*8                                  :: A
      integer                                 :: nbufsloc
      integer                                 :: ioresult

!
!---------------Description--------------------------------------------
!
!  Perform the real symmetric multiplication of the partial or complete
!  ADC matrix with PO and store the result in the vector PN
!
!---------------Executable code--------------------------------------
!
!  for the matrix elements we note that only the lower triangular
!  part *AND* the diagonal of the ADC-matrix is stored. 
!  the if takes care of not double counting the diagonal.

      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
!
!  nbufsloc is an indicator if the ADC matrix is complete or partial
!  In serial mode it is complete and has -1 at the beginning. Here we need the
!  iostat mechanism to determine the file end
!
!  In parallel mode nbufsloc contains the actual number of buffers. This is known due
!  to the matrix distribution procedure.
!
!  By this we avoid an awkward handling of the number of buffers as arguments, or common blocks...
!
      pn = 0.0d0

      rewind(iobase)
      read(iobase) nbufsloc

      If(nbufsloc.lt.0) then  !any negative number indicates serial
!___________________________________________________________
!|
!|   serial case, exact buffer number is not known
!|   we read stream blocks until end of file is reached.

      Do

        READ(iobase,iostat=ioresult) (BUF(IXX),IXX=1,INTBUF),
     &                               (IOI(IXX),IXX=1,INTBUF),
     &                               (IOJ(IXX),IXX=1,INTBUF),
     &                               NACT,JDUMMY
        if (ioresult > 0) then
           call quit('ADC matrix read error (real). Terminating')
        endif
        if (ioresult < 0) then  ! EOF detected
           exit    !leave infinite do loop
        endif

! useful data in buffer

        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
        ENDDO

      Enddo
!|
!|   
!|__________________________________________________________
      else
!___________________________________________________________
!|
!|   parallel case, exact buffer number is known
!|
        DO I = 1,nbufsloc
          READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                         (IOI(IXX),IXX=1,INTBUF),
     &                         (IOJ(IXX),IXX=1,INTBUF),
     &                          NACT,JDUMMY
          DO K = 1, NACT
            IROW = IOI(K)
            ICOL = IOJ(K)
            A    = BUF(K)

            PN(IROW) = PN(IROW) + PO(ICOL) * A
            IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
          ENDDO
        ENDDO
!|
!|   
!|__________________________________________________________

      endif

      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 889  CALL QUIT('Matrix read error in xmvmul_sp_r')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE xmvmul_sp_c(iobase,intbuf,ladc,PO,PN)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase,intbuf,ladc
      complex*16,dimension(:)        :: PO,PN
!
!------------------------ local variables -------------------------
!
      complex*16, allocatable, dimension(:)   :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      complex*16                              :: A
      integer                                 :: nbufsloc
      integer                                 :: ioresult

!
!---------------Description--------------------------------------------
!
!  Perform the hermitian multiplication of the external partial 
!  ADC matrix with PO and store the result in the vector PN
!
!---------------Executable code--------------------------------------
!
!  for the matrix elements we note that only the lower triangular
!  part *AND* the diagonal of the ADC-matrix is stored. 
!  the if takes care of not double counting the diagonal.

      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
!
!  nbufsloc is an indicator if the ADC matrix is complete or partial
!  In serial mode it is complete and has -1 at the beginning. Here we need the
!  iostat mechanism to determine the file end
!
!  In parallel mode nbufsloc contains the actual number of buffers. This is known due
!  to the matrix distribution procedure.
!
!  By this we avoid an awkward handling of the number of buffers as arguments, or common blocks...
!
      pn = (0.0d0,0.0d0)

      rewind(iobase)
      read(iobase) nbufsloc

      If(nbufsloc.lt.0) then  !any negative number indicates serial
!___________________________________________________________
!|
!|   serial case, exact buffer number is not known
!|   we read stream blocks until end of file is reached.

        Do

          READ(iobase,iostat=ioresult) (BUF(IXX),IXX=1,INTBUF),
     &                                 (IOI(IXX),IXX=1,INTBUF),
     &                                 (IOJ(IXX),IXX=1,INTBUF),
     &                                 NACT,JDUMMY
          if (ioresult > 0) then
             call quit('ADC matrix read error (complex). Terminating')
          endif
          if (ioresult < 0) then  ! EOF detected
             exit    !leave infinite do loop
          endif

! useful data in buffer

          DO K = 1, NACT
            IROW = IOI(K)
            ICOL = IOJ(K)
            A    = BUF(K)

             PN(IROW) = PN(IROW) + PO(ICOL) * A
            IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW)*dconjg(A)
          ENDDO

        Enddo
!|
!|   
!|__________________________________________________________
      else
!___________________________________________________________
!|
!|   parallel case, exact buffer number is known
!|
        DO I = 1,nbufsloc
          READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                         (IOI(IXX),IXX=1,INTBUF),
     &                         (IOJ(IXX),IXX=1,INTBUF),
     &                          NACT,JDUMMY
          DO K = 1, NACT
            IROW = IOI(K)
            ICOL = IOJ(K)
            A    = BUF(K)

            PN(IROW) = PN(IROW) + PO(ICOL) * A
            IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW)*dconjg(A)
          ENDDO
        ENDDO
!|
!|   
!|__________________________________________________________

      endif

      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 889  CALL QUIT('Matrix read error in xmvmul_sp_c')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XDIAG_INI_R(iobase,io_out,intbuf,nbufs,ladc,ascout)

!  Fetch the diagonal entries of the sparse matrix A
!  This routine is only called in the serial case and obtains NBUFS
!  furthermore it allows writing out the ADC matrix in ASCII format
!  for a brute force comparison.

      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase, io_out
      integer                        :: intbuf,nbufs,ladc
      logical                        :: ascout
!
!------------------------ local variables -------------------------
!
      Character*8,parameter                   :: diagfile='ADCDGTMP'
      integer                                 :: serialindicator
      REAL*8,  allocatable, dimension(:)      :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      REAL*8,  allocatable, dimension(:)      :: adiag
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      REAL*8                                  :: A
      logical                                 :: isthere

!
!---------------Description--------------------------------------------
!
!---------------Executable code--------------------------------------
!
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      allocate(adiag(ladc))

! open ascii file if activated

      if(ascout) then
         open(unit=io_out, file="ADCMAT.ASC", status='unknown',
     &        form="formatted")
         rewind(io_out)
      endif

c  .. read elements

      REWIND(iobase)
      read(iobase) serialindicator   !discarded, but should be less than zero!

      DO I = 1,NBUFS
        READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          IF(IROW.EQ.ICOL) adiag(irow) = A

          IF(ascout) write(io_out,'(2I6,2X,F20.6)') IROW,ICOL,A

        ENDDO
      ENDDO

! close file if it was addressed, file handle io_out becomes available again

      if(ascout) then
         close(io_out)
      endif

!      Inquire(file=diagfile,exist=isthere)
!      write(*,*) 'Found an ADC diagonal file. Deleting...'
!      If(isthere) then
!        open(unit=io_out, file=diagfile, status='old')
!        close(io_out, status='delete')
!      Endif

       open(unit=io_out,file=diagfile, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='NEW')
       write(io_out,ERR=887) (adiag(ixx),ixx=1,ladc)
       close(io_out)

      deallocate(adiag)
      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 887  CALL QUIT('Matrix write error in XDIAG_INI_R')
 889  CALL QUIT('Matrix read error in XDIAG_INI_R')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XDIAG_INI_C(iobase,io_out,intbuf,nbufs,ladc)

!  Fetch the diagonal entries of the sparse matrix A
!  This routine is only called in the serial case and obtains NBUFS

      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase, io_out
      integer                        :: intbuf,nbufs,ladc
!
!------------------------ local variables -------------------------
!
      Character*8,parameter                     :: diagfile='ADCDGTMP'
      integer                                   :: serialindicator
      complex*16, allocatable, dimension(:)     :: BUF
      INTEGER, allocatable, dimension(:)        :: IOI,IOJ
      complex*16, allocatable, dimension(:)     :: adiag
      INTEGER                                   :: I,K,IXX
      INTEGER                                   :: NACT,JDUMMY,IROW,ICOL
      complex*16                                :: A

!
!---------------Description--------------------------------------------
!
!---------------Executable code--------------------------------------
!
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      allocate(adiag(ladc))

c  .. read elements

      REWIND(iobase)
      read(iobase) serialindicator   !discarded, but should be less than zero!

      DO I = 1,NBUFS
        READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          IF(IROW.EQ.ICOL) adiag(irow) = A

        ENDDO
      ENDDO

       open(unit=io_out,file=diagfile, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='NEW')
       write(io_out,ERR=887) (adiag(ixx),ixx=1,ladc)
       close(io_out)

      deallocate(adiag)
      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 887  CALL QUIT('Matrix write error in XDIAG_INI_R')
 889  CALL QUIT('Matrix read error in XDIAG_INI_R')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE fulldia_r3(io,fname,intbuf,ladc,desrep)

      implicit none

      INTEGER                       :: io,intbuf,ladc,desrep
      CHARACTER(len=6)              :: fname
!
!---------------Description--------------------------------------------
!
!  For not too large symmetric matrices the full diagonalization is
!  performed for debugging purposes. This adapted version calculates
!  eigenvalues and eigenvectors and writes them to file "FULLEVC.XX"
!
!---------------Local variables--------------------------------------
!
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      REAL*8, allocatable             :: A(:,:),WORK(:),EVL(:)
      REAL*8, allocatable             :: BUF(:)
      INTEGER, allocatable            :: IOI(:),IOJ(:)
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO,LWORK
      integer                         :: ioresult,serialindicator
      character*10                    :: fullevcname
      Real*8                          :: D0
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!
      IF(LADC.gt.2000) THEN
        write(*,*) 'Full diagonalization skipped'
        return
      ENDIF

      CALL PST('Entering POLPRP REAL FULL diagonalizer+')

      D0 = 0.0d0

      open(unit=io,file=fname, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='unknown')

      REWIND(io)
      read(io) serialindicator   !discarded, but necessary.


      allocate(a(ladc,ladc))
      allocate(evl(ladc))
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      lwork=10*ladc
      allocate(work(lwork))
      
      A=0.0d0
!___________________________________________________________
!|
!|   serial case, exact buffer number is not known
!|   we read stream blocks until end of file is reached.

      Do

        READ(io,iostat=ioresult) (BUF(IXX),IXX=1,INTBUF),
     &                           (IOI(IXX),IXX=1,INTBUF),
     &                           (IOJ(IXX),IXX=1,INTBUF),
     &                            NACT,JDUMMY
        if (ioresult > 0) then
           call quit('ADC matrix read error (real). Terminating')
        endif
        if (ioresult < 0) then  ! EOF detected
           exit    !leave infinite do loop
        endif

! useful data in buffer

        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A(irow,icol) = BUF(K)
        ENDDO

      Enddo

      close(io)
!|
!|
!|__________________________________________________________

!
!  call symmetric diagonalizer
!  calculate only eigenvalues being placed in A
!
      JOBZ='V'
      UPLO='L'

      CALL DSYEV( JOBZ, UPLO, LADC, A, LADC, EVL, WORK, LWORK, INFO )

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in DSYEV. Info:',INFO
           CALL QUIT ('*** Error in fulldia_r3  ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 20):'
      DO I=1,min(20,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO
!
!  write full set of eigenvectors to file
!  Keep hereby the *same* format as in the DAVILONG file in order
!  to make comparison more readable
!
      IF(desrep.GT.9) THEN
        WRITE(fullevcname,'(A8,I2)') 'FULLEVC.',desrep
      ELSE
        WRITE(fullevcname,'(A9,I1)') 'FULLEVC.0',desrep
      ENDIF
      open(io,file=fullevcname,access='sequential',
     &     status='unknown',form='unformatted')
      do k=1,ladc
!       write(*,*) 'Wrinting full eigenvector #',k
        write(io) k
        write(io) evl(k)
        write(io) D0
        write(io) (A(i,k),i=1,ladc)
      enddo
      close(io)

!
!  release all memory
!
      deallocate(a)
      deallocate(evl)
      deallocate(buf)
      deallocate(ioi)
      deallocate(ioj)
      deallocate(work)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE fulldia_c3(io,fname,intbuf,ladc)

      implicit none

      INTEGER                       :: io,intbuf,ladc
      CHARACTER(len=6)              :: fname
!
!---------------Description--------------------------------------------
!
!  For not too large Hermitian matrices the full diagonalization is
!  performed for debugging purposes. The required arrays are allocated
!  dynamically in order to achieve maximal independence.
!
!  The only specific entities are two available file handles,
!  the length of the write chunks and the length of the ADC matrix.  
!
!---------------Local variables--------------------------------------
!
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      COMPLEX*16,   allocatable       :: A(:,:),CWORK(:),BUF(:)
      REAL*8,       allocatable       :: EVL(:),RWORK(:)
      INTEGER,      allocatable       :: IOI(:),IOJ(:)
      INTEGER                         :: LWORK
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO
      INTEGER                         :: ioresult,serialindicator
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!


      IF(LADC.gt.5000) THEN                                
        write(*,*) 'Full diagonalization skipped'
        return                                             
      ENDIF

      CALL PST('Entering POLPRP FULL complex diagonalizer+')

      open(unit=io,file=fname, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='unknown')

      REWIND(io)
      read(io) serialindicator   !discarded, but necessary.


      allocate(a(ladc,ladc))
      allocate(evl(ladc))
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      lwork=10*ladc
      allocate(rwork(lwork))
      allocate(cwork(lwork))
      
      A=(0.0d0,0.0d0)
!_________________________ reading section
!|
      Do

        READ(io,iostat=ioresult) (BUF(IXX),IXX=1,INTBUF),
     &                           (IOI(IXX),IXX=1,INTBUF),
     &                           (IOJ(IXX),IXX=1,INTBUF),
     &                            NACT,JDUMMY
        if (ioresult > 0) then
           call quit('ADC matrix read error (complex). Terminating')
        endif
        if (ioresult < 0) then  ! EOF detected
           exit    !leave infinite do loop
        endif

! useful data in buffer

        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A(irow,icol) = BUF(K)
        ENDDO

      Enddo

      close(io)

!|________________________ end reading section

!
!  call hermitian diagonalizer
!  calculate only eigenvalues being placed in A
!
      JOBZ='N'
      UPLO='L'

      CALL ZHEEV( JOBZ, UPLO, LADC, A, LADC, EVL, CWORK,
     &            LWORK, RWORK, INFO )

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in ZHEEV. Info:',INFO
           CALL QUIT ('*** Error in fulldia_c3 ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 20):'
      DO I=1,min(20,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO
!
!  release all memory
!
      deallocate(a)
      deallocate(evl)
      deallocate(buf)
      deallocate(ioi)
      deallocate(ioj)
      deallocate(rwork)
      deallocate(cwork)

      RETURN

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QMOOCCALC_R(ioda,qmu,aij,ladc,maxdavsp,locdim,iroot)

      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      INTEGER             :: ioda,ladc,maxdavsp,locdim,iroot
      REAL*8              :: qmu(ladc),aij(maxdavsp,maxdavsp)
!
!------------------------ local variables -------------------------
!
      INTEGER                                 :: i,j
      real*8, allocatable, dimension(:)       :: aux
!
!---------------Description--------------------------------------------
!
!   This routine calculates the qm vector using the new b-vectors
!   from disk
!
      allocate(aux(ladc))
      qmu = 0.0d0
      do i=1,locdim
        read(ioda,rec=i,err=888) aux
        qmu(:) = qmu(:) + aij(i,iroot)*aux(:)
      enddo
      deallocate(aux)
      return

 888  call quit('Vector read error in QMOOCCALC_R!')

      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QMOOCCALC_C(ioda,qmu,aij,ladc,maxdavsp,locdim,iroot)
!
!  the call is:
!
!  Call QMOOCCALC_C(iodachannel,qm,aij,ladc,maxdavsp,locdim,k)
!

      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      INTEGER                      :: ioda,ladc,maxdavsp,locdim,iroot
      complex*16, dimension(:)     :: qmu
      complex*16, dimension(:,:)   :: aij
!
!------------------------ local variables -------------------------
!
      INTEGER                                     :: i,j
      complex*16, allocatable, dimension(:)       :: aux
!
!---------------Description--------------------------------------------
!
!   This routine calculates the qm vector using the new b-vectors
!   from disk
!
      allocate(aux(ladc))
      qmu = (0.0d0,0.0d0)
      do i=1,locdim
        read(ioda,rec=i,err=888) aux
        qmu(:) = qmu(:) + aij(i,iroot)*aux(:)
      enddo
      deallocate(aux)
      return

 888  call quit('Vector read error in QMOOCCALC_C!')

      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DISORT(n,arr,brr)
!
! adapted from numerical recipes in fortran 77
! to simultaneously shuffle the index array brr
! with the sorted array arr
!
      implicit none
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      INTEGER brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,atemp
      INTEGER b,btemp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        atemp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=atemp
        btemp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=btemp
        if(arr(l+1).gt.arr(ir))then
          atemp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=atemp
          btemp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=btemp
        endif
        if(arr(l).gt.arr(ir))then
          atemp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=atemp
          btemp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=btemp
        endif
        if(arr(l+1).gt.arr(l))then
          atemp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=atemp
          btemp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=btemp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        atemp=arr(i)
        arr(i)=arr(j)
        arr(j)=atemp
        btemp=brr(i)
        brr(i)=brr(j)
        brr(j)=btemp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)then
          stop 'increase NSTACK in DISORT'
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WCONDAT_X(LUN,FILENAME,IRECL,DESREP,RECCOUNT)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Writes explicit configuration information for the excited
C     states in symmetry DESREP.
C     Important: The absolute spinor numbers calculated here correspond
C     to the list in the CCSETI output and not the irepspi numbering!
C     This makes retrieval of orbitals easier for the user
C
C---------------Calling variables--------------------------------------
C
      INTEGER LUN
      CHARACTER*6 FILENAME
      INTEGER IRECL
      INTEGER DESREP
      INTEGER RECCOUNT
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*70 FIELD
      INTEGER bl_conf
      CHARACTER*4 bl_irna
      
C
C---------------Executable code--------------------------------------
C
C open config file
C
      IRECL = 70
      OPEN(LUN,FILE=FILENAME,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
C
C First do the P--->H part. Column counter is the decisive
C IRREP.
C
      BL_CONF = 0
      BL_IRNA = '    '
      RECCOUNT = 1

C____________ START WITH H ---> P __________________________
C|
C|
C|

!     write(*,*) '  **** wcondat output ****'
      do jrep = 1,nrep   !jrep,brep are fermionic, desrep is bosonic
        brep = multb(desrep+nrep,jrep,2)
!       write(*,*) '  B:',desrep,'  G_a:',brep,'  G_i:',jrep
        do jfun = 1,NO(jrep)
        do bfun = 1,NV(brep)
          jabs=IO(jrep) + jfun
          babs=IO(nrep+1) + IV(brep) + bfun

          WRITE(FIELD,333)
     &         babs,brep,REPNA(brep),
     &         jabs,jrep,REPNA(jrep),
     &         BL_CONF,BL_CONF,BL_IRNA,
     &         BL_CONF,BL_CONF,BL_IRNA 
 333  FORMAT(2(I3,' (',I2,',',A4,')  <---- ',I3,' (',I2,',',A4,') '))

          WRITE(LUN,REC=RECCOUNT) FIELD
          RECCOUNT = RECCOUNT + 1

        enddo
        enddo
      enddo
      if( (reccount-1).ne.MVO(desrep))
     &   stop 'Inconsistent rec counter in WCONDAT_X!'
      write(*,*)
C|
C|
C|
C|__________________________________________________________
C
C
C____________ NOW DO  2H ---> 2P ___________________________
C|
C| important: exactly the same order as the column counter in
C| XSATMAKE! Otherwise config. information will be erroneous
C|
C|
C|
      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)

          kabs =IO(krep) + K
          labs =IO(lrep) + L
          ccabs=IO(nrep+1) + IV(crep) + C
          ddabs=IO(nrep+1) + IV(drep) + D

          WRITE(FIELD,333)
     &         ccabs,crep,REPNA(crep),
     &         kabs,krep,REPNA(krep),
     &         ddabs,drep,REPNA(drep),
     &         labs,lrep,REPNA(lrep)

          WRITE(LUN,REC=RECCOUNT) FIELD
          RECCOUNT = RECCOUNT + 1

                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C|
C|
C|
C|__________________________________________________________

      CLOSE(LUN)

      RECCOUNT = RECCOUNT - 1
      write(*,*) 'Wcondat info:',reccount,' conf recs written.'

 
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMVVOO_PP (EPS,T2)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Divide VVOO integrals by denominators ea + eb - ej - ei
!     If we have complex numbers, the program accounts for this by
!     treating real/imaginary part in the *real* array.
!
!---------------Calling variables--------------------------------------
!
      REAL*8,dimension(:)                 ::  EPS
      REAL*8,dimension(:)                 ::  T2
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8 FAC,FAC1,FAC2,FAC3
!
!---------------Executable code--------------------------------------
!
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = - EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = - EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 + EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 + EPS(AA)
                  ABIJ = ABIJ + 1
                  IF (CARITH) THEN
                     T2(2*ABIJ-1) = T2(2*ABIJ-1)/FAC
                     T2(2*ABIJ)   = T2(2*ABIJ)/FAC
                  ELSE
                     T2(ABIJ) = T2(ABIJ)/FAC
                  ENDIF
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      IF(ABIJ.ne.IVVOOTT(nrep+1)) THEN
        WRITE(*,*) 'ABIJ/IVVOOTT:',ABIJ,IVVOOTT(NREP+1)
        STOP 'COUNTING ERROR IN DENOMVVOO_PP!'
      ENDIF
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
