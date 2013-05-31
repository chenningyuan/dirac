!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!      *****************************************************
!      ****        D H F - A D C (2x)  M O D U L E      ****
!      ****           POLARIZATION PROPAGATOR           ****
!      ****                   for                       ****
!      ****           ELECTRONIC EXCITATIONS            ****
!      *****************************************************
!
!                       AUTHOR: M.Pernpointner
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADCEXC_MAIN (IW,IOBASE,EPS)
      use adc_cfg
      use adc_mat
 
      IMPLICIT INTEGER(A-Z)

      integer iw       ! I/O file handle for screen output
      integer IOBASE   ! I/O base file handle, all purposes
      real*8  eps(*)   ! array of spinor energies, assumed size array!
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
!
!---------------Local variables---------------------------------
!
! file names for the Lanczos
      CHARACTER*6                            :: NAMEMATFIL
      CHARACTER*6                            :: NAMECNFFIL
      CHARACTER*5                            :: NMSPECF
      PARAMETER(NMSPECF='XSPEC')

! variables for the lanczos

      integer                 :: IONIZLEVEL, IRECL, INTBUF
      integer                 :: irploop,i,j,k,l,m,n
      integer                 :: n1,n2
      integer                 :: desrep,lencouplbl,neigenv
      integer                 :: icount,ndim,ialloc
      integer                 :: IDPT  !dummy pointer
      real*8, allocatable     :: ciajb(:),daibj(:)
      real*8, allocatable     :: cc(:)   ! work space for diag_lanc
!
!-----------------------interface section---------------------
!
      interface

        SUBROUTINE BUILD_X_MAIN (a,b,l)
          real*8, dimension(:)  :: a
          real*8, dimension(*)  :: b
          integer               :: l
        END SUBROUTINE

        SUBROUTINE SOLVE_HPHP(ra1,i1,i2)
          real*8, dimension(:)  :: ra1 
          integer               :: i1,i2
        END SUBROUTINE

        SUBROUTINE BUILD_X_COUSAT_R(h,a,i1,i2,i3,i4,lo,r2,i5)
          real*8, dimension(:)  :: h
          real*8, dimension(*)  :: a
          integer               :: i1,i2,i3,i4
          logical               :: lo
          real*8                :: r2
          integer               :: i5
        END SUBROUTINE

      end interface
!
!-----------------------end interface section---------------------
!

!
!  _________________________ INTRO SECTION __________________________
! |
! |
! |
      WRITE(IW,*)
      WRITE(IW,*) '***************************************************'
      WRITE(IW,*) '*** '
      WRITE(IW,*) '***                  STARTING THE'
      WRITE(IW,*) '*** '
      WRITE(IW,*) '***         RELATIVISTIC EXCITATION ADC(2x)'
      WRITE(IW,*) '*** '
      WRITE(IW,*) '***                     MODULE'
      WRITE(IW,*) '*** '
      WRITE(IW,*) '***************************************************'
      WRITE(IW,*)
      IF(reladc_doxext) THEN
        WRITE(IW,*) 'Doing extended ADC-2'
      ELSE
        WRITE(IW,*) 'Doing strict ADC-2'
      ENDIF

      IF(reladc_adcthr.gt.0.0d0) THEN
        WRITE(IW,*) 'A threshold of',reladc_adcthr,' will be applied to'
        WRITE(IW,*) 'the PP ADC matrix.'
      ENDIF
!
! determine symmetries to be calculated
!
      n = 0
      do i = 1,nrep
        n = n + reladc_excreps(i)
      enddo
      if(n.eq.0) then
        write(iw,*) 'Program checks symmetries automatically'
        do i=1,nrep
          reladc_excreps(i) = i
        enddo
      else
        write(iw,*) 'User has made the following symmetry choices:'
        write(iw,*)
        icount=1
        do i=1,nrep
          if(reladc_excreps(i).ne.0) then
            write(iw,*) '   No: ',icount,
     &                  '  symmetry: ',reladc_excreps(i)
            icount = icount + 1
          endif
        enddo
      endif
      write(iw,*)

      write(IW,'(10X,A)') 'Requested EXC eigenvectors:'
      n=0
      do i=1,32
        n=n+reladc_exceigv(i)
      enddo
      if(n.eq.0) then
        write(IW,'(10X,A)') 'None.'
      else
        write(IW,'(10X,16I4)') (reladc_exceigv(i),i=1,16)
        write(IW,'(10X,16I4)') (reladc_exceigv(i),i=17,32)
      endif

      INTBUF = 5*1024*1024     !we work with a 5 Megaword buffer.
! |
! |
! |
! |________________END INTRO SECTION__________________________
!
!
!  __________________  BUILD HP-HP (MAIN) BLOCK  _____________
! |
! |
! |
!   all symmetries in the HP-HP block are created at a time
!   this facilitates easy handling but costs a bit more memory.
!   allocate hphp block according to array:
!   IVOVO  (IREP+1) = IVOVO  (IREP) + NVO (IREP)*NVO (IREP)
!   This chunk of memory stays allocated also during the calc.
!   of the other blocks.
!
!   on exit of BUILD_X_MAIN the hp/hp block ist stored in the
!   BKC convention <|<|,|>|>, corresponding to the generic
!   VOVO storage mode.
!
      n1 = ivovo(nrep+1) * rcw
      allocate(ciajb(n1),stat=ialloc)
      if(ialloc.ne.0) stop 'allocation problem CIAJB'
!
! do the work
!
      CALL BUILD_X_MAIN (CIAJB,EPS,N1)
!
      n2 = jvovo(nrep+1) * rcw
      allocate(daibj(n2),stat=ialloc)
      if(ialloc.ne.0) stop 'allocation problem AUX'
!
!   We do one SRT16 resort and store everything as C_ai,bj
!   and not as C_ia,jb as written in the manuscript.
!   SRT16: <aj||bi> -->  ai,bj with <||>,|><|  BKC
!
      CALL SRT16(NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,
     &            JJVO,JJVO,ciajb,daibj)
!
! diagonalize main block for Davidson start vectors
!
      CALL SOLVE_HPHP(daibj,n2,iobase)
! |
! |
! |_____________  E_N_D__BUILD HP-HP BLOCK  _________________________
!

!  ____________ SYMMETRY LOOP FOR COUPLING + SAT BLOCK___________
! |
! |____________ TOGETHER WITH DIAGONALIZING IN EACH SYMMETRY ____
! |
! |
      NAMEMATFIL='ADCXPM'

      DO IRPLOOP = 1,NREP

        open(unit=iobase,file=NAMEMATFIL, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='unknown')

        IF(reladc_excreps(IRPLOOP).EQ.0) CYCLE
        IF(reladc_excreps(IRPLOOP).GT.NREP) THEN
           WRITE(IW,*) 'Chosen Symmetry',reladc_excreps(IRPLOOP),
     &     ' does not exist and is skipped!'
           CYCLE
        ENDIF

        DESREP = reladc_excreps(IRPLOOP)
        IF(MVO(DESREP).EQ.0) THEN
          WRITE(iw,*) 'No final states in symmetry',desrep,repna(desrep)
          CYCLE
        ENDIF
        WRITE(IW,*)
        WRITE(IW,*)
        WRITE(IW,*)
        WRITE(IW,'(A)') '********************************************'
        WRITE(IW,'(A,I4,X,A4)') 'Starting calculations in symmetry:',
     &         desrep,repna(desrep)
        WRITE(IW,'(A)') '********************************************'
        WRITE(IW,*)
        WRITE(IW,*)
!       CALL PST('Computing COUPLING & SAT blocks+')
        IF(CARITH) THEN
          stop 'Complex arithm. for excitations not yet available.'
        ELSE
          CALL BUILD_X_COUSAT_R(DAIBJ,EPS,DESREP,LENCOUPLBL,
     &                          INTBUF,NBUFS,reladc_doxext,
     &                          reladc_adcthr,iobase)
        ENDIF

        CLOSE(iobase)
!
! ______________________________________
!| prepare diagonalization              |
!|______________________________________|
!
! configuration information file
! is created within WCONDAT_X.
! file handle is reusable after call!
!
        NAMECNFFIL='ADCXCF'
!       WRITE(IW,*)
!       WRITE(IW,*) 'Creating configuration information'
!       WRITE(IW,*) 'For symmetry ',desrep
!       WRITE(IW,*)
        CALL WCONDAT_X(iobase,namecnffil,irecl,desrep)
!       WRITE(IW,*) 'Record length from WCONDAT_X:',IRECL
!
!  write matrix control file and call corresponding diagonalizer
!
        IONIZLEVEL = 3
        LADC = LENCOUPLBL + MVO(DESREP)

!
!  call full diagonalizer for new matrix format if requested.
!
! ______________________________________
!| Call full diagonalizer, if activated |
!|______________________________________|

        IF(reladc_dofull) THEN
          IF(CARITH) THEN
            CALL FULLDIAC2(IOBASE, NAMEMATFIL, INTBUF, NBUFS,
     &                     LADC, DESREP)
          ELSE
            CALL FULLDIAR2(IOBASE, NAMEMATFIL, INTBUF, NBUFS,
     &                     LADC, DESREP)
          ENDIF
        ENDIF
!
! ______________________________________
!| Call Lanczos diag.      if activated |
!|______________________________________|
!  if called from ADCEXC_MAIN we need to allocate the CC
!  array in order to make it functional.
!  Later, this will be abandoned completely.
!
        IF(reladc_dolanc) THEN
!
!** setup matrix descriptor for Lanczos
!
          neigenv = reladc_exceigv(IRPLOOP)

          reladc_md_iobase   = iobase
          reladc_md_ionizl   = IONIZLEVEL
          reladc_md_ioldnew  = 2
          reladc_md_intbuf   = INTBUF
          reladc_md_desrep   = DESREP
          reladc_md_rcw      = RCW
          reladc_md_lnzitr   = reladc_exciter
          reladc_md_matdim   = LADC
          reladc_md_irecl    = IRECL
          reladc_md_nmain    = MVO(DESREP)
          reladc_md_nbufs    = NBUFS
          reladc_md_neigenv  = neigenv
          reladc_md_fileadc  = NAMEMATFIL
          reladc_md_filediag = '      '
          reladc_md_filecnf  = NAMECNFFIL
          reladc_md_nmspec   = NMSPECF

          allocate(CC(reladc_lancmem),stat=ialloc)
          if(ialloc.ne.0) stop 'Memory allocation problem.'

          KFREE = 1; LFREE = reladc_lancmem; KFREEBASE = 1
!
! do one "allocation" and "deallocation" in order to 
! initialize mempkg routines.
! ask for LANCMEM-10 in order to allow space for tags.
!
          CALL MEMOPT('REAL',IDPT,reladc_lancmem-10,CC,
     &                KFREE,LFREE,.true.)
          CALL MEMREL('IDPT',CC,KFREEBASE,KFREEBASE,KFREE,LFREE)
!
!  perform Lanczos
!
          CALL DIAG_LANC(CC,KFREE,LFREE,IW,reladc_doincore)
!
! free all memory that was used for diagonalization
!
          deallocate(CC)
        ENDIF
! ______________________________________
!| Call Davidson diag.     if activated |
!|______________________________________|
!
        IF(reladc_dodavi) THEN
!
!  write current matrix specification into descriptor block
!  only the values **** needed for Davidson **** are written.
!  In the descriptor file a few more variables are contained.
!
          reladc_md_iobase    = iobase
          reladc_md_intbuf    = intbuf
          reladc_md_desrep    = desrep
          reladc_md_rcw       = rcw
          reladc_md_matdim    = ladc
          reladc_md_irecl     = irecl
          reladc_md_nmain     = mvo(desrep)
          reladc_md_nbufs     = nbufs
          reladc_md_davroots  = reladc_davroots
          reladc_md_maxdavsp  = reladc_davmaxsp
          reladc_md_maxdavit  = reladc_davmaxit
          reladc_md_convthr   = reladc_davconv
          reladc_md_fileadc   = NAMEMATFIL
          reladc_md_filecnf   = NAMECNFFIL
          reladc_md_nmspec    = NMSPECF 

          CALL DIAG_DAVI(IW)
        ENDIF
!
!  delete adc matrix file
!
        open(unit=iobase,file=NAMEMATFIL, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='unknown')
        close(iobase,status='DELETE')

      ENDDO    !IRPLOOP
! |
! |
! |
! |______________END___SYMMETRY + DIAGONALIZATION___LOOP_____________

      RETURN
      END
! |
! |
! |
! |________________E_N_D___M A I N___P R O G R A M __________________

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_X_MAIN (CIAJB,EPS,DIMCIAJB)
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
      REAL*8,dimension(:)         :: CIAJB  !assumed shape array. Was declared within F90
      REAL*8,dimension(*)         :: EPS(*) !ATT* must be declared as assumed size array
!                                            (with *)
!                                    otherwise the transfer as an argument to the old F77
!                                    routines will cause a program crash! This is not the
!                                    optimal convention. The compiler has *no* information
!                                    about the shape of the array. Therefore assumed size.
      INTEGER                     :: DIMCIAJB
!
!---------------Common Blocks, symmetry arrays from CCSETI --------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
!
!---------------Local variables--------------------------------------
!
      real*8, allocatable        :: sa_vvoott(:)
      real*8, allocatable        :: sb_vvoott(:)
      real*8, allocatable        :: buf1(:),buf2(:)
      real*8, allocatable        :: buf3(:),buf4(:)
      real*8, allocatable        :: buf5(:)
      real*8, allocatable        :: bab(:)
      real*8, allocatable        :: dij(:)

      integer                    :: ialloc
!
!---------------Executable code--------------------------------------
!
      IF(CARITH) THEN
        CALL PST('Building complex HP/HP block+')
      ELSE
        CALL PST('Building real HP/HP block+')
      ENDIF

!******************* preparations
!
!  we often need the v_ac,kl array which is V_ac,kl (VVOO)
!  divided by the energy denominator
!  e_a + e_c - e_k - e_l and the original V_ac,kl
!  we therefore calculate it once and store it in this subroutine.
!  the variables icX count the number of integrals. For the actual
!  storage requirements one has to multiply with RCW!
!
      ic1 = ivvoott(nrep+1)
      allocate(sa_vvoott(ic1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'allocation problem 1 in BUILD_X_MAIN'
      allocate(sb_vvoott(ic1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'allocation problem 2 in BUILD_X_MAIN'

      call getvvoo(sa_vvoott)
      call xcopy(ic1,sa_vvoott,1,sb_vvoott,1)
      call denomvvoo(eps,sa_vvoott,sa_vvoott)
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
      allocate(bab(ic2*rcw)); bab = 0.0d0

      ic3 = max(lvvoot(nrep+1),kvvoot(nrep+1))
      allocate(buf1(ic3*rcw))  ! array has to be slightly larger due to
      allocate(buf2(ic3*rcw))  ! A,BIJ storage and/or ABI,J storage.
      buf1 = 0.0d0; buf2 = 0.0d0

      ic4 = NFOO
      allocate(dij(ic4*rcw)); dij = 0.0d0

      ic5 = max(j2vovo(nrep+1),jvovo(nrep+1))
      allocate(buf3(ic5*rcw))  ! all arrays have to hold VOVO integrals
      allocate(buf4(ic5*rcw))  ! in MVO*MVO instead of NVO*NVO manner!
      allocate(buf5(ic5*rcw))
      buf3 = 0.0d0; buf4 = 0.0d0; buf5 = 0.0d0

      ic6 = ivovo(nrep+1)      ! size of the generic VOVO array in generic BKC
!
!  -------------------------------
!  Calculate    TERM C(2)A       |
!  result in bab  <||> BKC       |
!  -------------------------------
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
!  -------------------------------
!  Calculate    TERM C(2)B       |
!  result in dij  <||> = D*_ji   |
!  -------------------------------
!
!  sort v_cd,ki   -->  v_cdk,i    (buf1)
!  sort V_cd,kj   -->  V_cdk,j    (buf2)
!  V_cdk,j        -->  V*_cdk,j   (buf2)
!  v_cdk,i  *  V*_cdk,j  -->  Dij   -->  'hermitize'
!  c.c. Dij   -->  D_ji with <||>  BKC
!
      buf1 = 0.0d0; buf2 = 0.0d0
      call srt1t3 (nrep,multb,lfa,nvvt,no,no,nvvot,kvvoot,kkvvot,
     &             sa_vvoott,buf1)
      if(carith) call conjuga(ic1,sb_vvoott,1)
      call srt1t3 (nrep,multb,lfa,nvvt,no,no,nvvot,kvvoot,kkvvot,
     &             sb_vvoott,buf2)
      if(carith) call conjuga(ic1,sb_vvoott,1)
      call cntrct ('T','N',no,no,nvvot,a1,buf1,buf2,a1,dij,nrep)
      do krep=1,nrep
        ioff=jjoo(krep,krep)*rcw + 1
        if(carith) then
          call cmatsym(dij(ioff),no(krep),no(krep))
        else
          call matsym(dij(ioff),no(krep),no(krep))
        endif
      enddo
      if(carith) call conjuga(ic4,dij,1)
!
!  ----------------------------
!  Calculate    TERM C(2)C    |
!  result in buf3.            |
!  ----------------------------
!
!  sort v_ac,ik    -->  v_ai,ck    (buf3)
!  sort V*_bc,jk   -->  V*_bj,ck   (buf4)
!
      call srt1tt4 (nrep,multb,lfa,lfa,nv,nv,no,no,
     &              mvo,j2vovo,jjvo,jjvo,sa_vvoott,buf3)
      if(carith) call conjuga(ic1,sb_vvoott,1)
      call srt1tt4 (nrep,multb,lfa,lfa,nv,nv,no,no,
     &              mvo,j2vovo,jjvo,jjvo,sb_vvoott,buf4)
      if(carith) call conjuga(ic1,sb_vvoott,1)
!
!  here again special contraction, with 'T' in XGEMM
!  we end up as v_ai,ck * V*_bj,ck  -->  ai,bj intermediate (buf5)
!  having BKC:  <| |>,|><|.  
!
      ioff1 = 1
      ioff3 = 1
      do irep = 1, nrep
         jrep = multb(irep+nrep,1+nrep,2)
         m = mvo(irep)
         n = mvo(irep)
         k = mvo(jrep)
         ioff2 = j2vovo(jrep) * rcw + 1
         call xgemm ('N','T',m,n,k,a1,buf3(ioff1),m,buf4(ioff2),n,
     &               a0,buf5(ioff3),m)
         ioff1 = ioff1 + m * k * rcw
         ioff3 = ioff3 + m * n * rcw
      enddo
!
!  We have to backsort in order to arrive at
!  Y(1)_aj,bi order with <|<||>|>  BKC.
! do backsort   ATTention ***** deviation of array order in SRT 16 inverse mode !!!!!
!
!  ai,bj (buf5) ---->  Y(1)_aj,bi  (buf3)  with generic VOVO BKC.
!
!  in SRT 16 the *second* buffer is the source buffer when it comes to inverse !  sorting.
!
      call srt16 (nrep,multb,lfa,ltr,nv,no,nv,no,mvo,jvovo,jjvo,jjvo,
     &            buf3,buf5)
!
!  We have to symmetrize according to A_aj,bi = Y(1)_aj,bi + Y(2)*_bi,aj
!  Hereby we again interpret the Y(1)_aj,bi array as a two-dimensional matrix
!  of NVO(i) * NVO(i) submatrices contiguously stored (VOVO storage scheme!)
!
      do krep=1,nrep
        ioff=ivovo(krep)*rcw + 1
        if(carith) then
          call cmatsym(buf3(ioff),nvo(krep),nvo(krep))
        else
          call matsym(buf3(ioff),nvo(krep),nvo(krep))
        endif
      enddo
!
! now we have B_ab, D_ij and A_aj,bi. but without the prefactors!
!
! the array holding the hp/hp block is ciajb (name in the paper) and has the size
! of the VOVO integrals. All contribs are added to this array initialized with the
! VOVO integrals.
!
!  ----------------------------------
!  Calculate    TERM C(1)           |
!  which is the negative VOVO ints  |
!  ----------------------------------
!
      call getvovo(ciajb)
      call xnegate(dimciajb,ciajb,ciajb)
!
!  add -0.5 *  A_aj,bi array to it (VOVO storage)
!  this is the C(2C) contribution
!
      call xaxpy  (ic6,-ap5,buf3,1,ciajb,1)
!
!  'weave' into this array the C(0), the C(2A) and C(2B) contributions
!  hereby C(0) is added, the C(2A) and C(2B) contributions are formed from
!  the Bab and Dij!
!
      if(carith) then
        call bweaverc(ciajb,eps,bab,dij,ic6,ic2,ic4)
      else
        call bweaverr(ciajb,eps,bab,dij,ic6,ic2,ic4)
      endif
!
!  the ciajb array now contains the full hp/hp block in generic VOVO
!  storage convention.
!
!  we release all memory allocated so far and retain only the ciajb block.
!
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
      SUBROUTINE BWEAVERR(bhphp,eps,bab,dij,ic6,ic2,ic4)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Crawl through the VOVO structure and distribute the
!     C(0), C(2A) and C(2B) contributions to the hphp block
!
!---------------Calling variables--------------------------------------
!
      integer, intent(in)         :: ic6, ic2, ic4
      real*8, intent(inout)       :: bhphp(ic6)
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
!  the loop structure follows the VO||VO convention with <|<|  |>|>
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
              bhphp(aibj) = bhphp(aibj) + esum
            ENDIF
!
!  treat the C(2A) contribution: C(2A) = +0.5*Delta_ij*B_ab
!
            IF(IRP.eq.JRP) THEN
              IF(I.eq.J) THEN
                if(ARP.NE.BRP) call quit('bweaverr_error C(2A)')
                aboff = IVV(BRP) + (A-1)*NV(BRP) + B
                bhphp(aibj) = bhphp(aibj) + 0.5d0*bab(aboff)
              ENDIF
            ENDIF
!
!  treat the C(2B) contribution: C(2B) = +0.5*Delta_ab*D_ij
!
            IF(ARP.eq.BRP) THEN
              IF(A.eq.B) THEN
                if(IRP.NE.JRP) call quit('bweaverr_error C(2B)')
                ijoff = IOO(JRP) + (I-1)*NO(JRP) + J
                bhphp(aibj) = bhphp(aibj) + 0.5d0*dij(ijoff)
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
      SUBROUTINE BWEAVERC(bhphp,eps,bab,dij,ic6,ic2,ic4)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Count through the VOVO integrals and distribute the
!     C(0), C(2A) and C(2B) contributions to the hphp block
!
!---------------Calling variables--------------------------------------
!
      integer, intent(in)         :: ic6, ic2, ic4
      complex*16, intent(inout)   :: bhphp(ic6)
      real*8, intent(in)          :: eps(*)
      complex*16, intent(in)      :: bab(ic2)
      complex*16, intent(in)      :: dij(ic4)
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
!  the loop structure follows the VO||VO convention with <|<|  |>|>
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
              bhphp(aibj) = bhphp(aibj) + dcmplx(esum,0.0d0)
            ENDIF
!
!  treat the C(2A) contribution: C(2A) = +0.5*Delta_ij*B_ab
!
            IF(IRP.eq.JRP) THEN
              IF(I.eq.J) THEN
                if(ARP.NE.BRP) call quit('bweaverc_error C(2A)')
                aboff = IVV(BRP) + (A-1)*NV(BRP) + B
                bhphp(aibj) = bhphp(aibj) + 0.5d0*bab(aboff)
              ENDIF
            ENDIF
!
!  treat the C(2B) contribution: C(2B) = +0.5*Delta_ab*D_ij
!
            IF(ARP.eq.BRP) THEN
              IF(A.eq.B) THEN
                if(IRP.NE.JRP) call quit('bweaverc_error C(2B)')
                ijoff = IOO(JRP) + (I-1)*NO(JRP) + J
                bhphp(aibj) = bhphp(aibj) + 0.5d0*dij(ijoff)
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
      SUBROUTINE SOLVE_HPHP(daibj,n,iounit)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Diagonalizes the hphp block in each symmetry separately.
!     For the singles amplitudes a local coniguration analysis is
!     provided for convenience.
!
!---------------Calling variables--------------------------------------
!
      real*8,dimension(:)         :: daibj
      integer                     :: n,iounit
!
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
      write(*,*) 'Generating Davidson start vectors.'

!  MUTAB handling check
!  we do it in the **COLUMN** index pair because we
!  construct the matrix columnwise!
!
!  According to the notes, bj is the decisive pair.
!
      do desrep=1,nrep   !desrep is a bosonic loop (beta_col)
        do jrep = 1,nrep
          brep = multb(desrep+nrep,jrep,1)
          trep = multb(jrep,brep,2)
          if(desrep.ne.trep) STOP 'error in MUTAB handling!'
        enddo
      enddo
!
!  create local single configuration descriptor array
!  remember: first index is virtual, second is occupied
!  this corresponds to the storage mode!
!
      icount = 0
      nx = jvoxx(nrep+1)
      allocate(cfield(nx))

      do desrep=1,nrep   !desrep is a bosonic loop
        do jrep = 1,nrep   !jrep,brep are fermionic
          brep = multb(desrep+nrep,jrep,1)
          do jfun = 1,NO(jrep)
          do bfun = 1,NV(brep)
            jabs=IO(jrep) + jfun
            babs=IO(nrep+1) + IV(brep) + bfun
            icount = icount + 1
            write(cfield(icount),777) jabs,repna(jrep),babs,repna(brep)
!           write(*,*) cfield(icount)
          enddo
          enddo
        enddo
      enddo
 777  Format(I3,' (',A4,') --->',1X,I3,' (',A4,')')
      if(icount.ne.jvoxx(nrep+1)) then
         write(*,*) icount,jvoxx(nrep+1)
         stop 'conf array error!'
      endif
!
!  diagonalize hphp submatrices
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
          write(*,*) 'No states in symmetry',desrep,repna(desrep)
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
        if(ialloc.ne.0) STOP 'Allocation problem 1 in SOLVE_HPHP'
        allocate(xmat(nsiz*nsiz),stat=ialloc)
        if(ialloc.ne.0) STOP 'Allocation problem 2 in SOLVE_HPHP'
        amat = 0.0d0; xmat = 0.0d0
!
!  we convert the linear MVO*MVO array to quadratic array
!  hereby we need xmat, compiler does not understand it otherwise
!
        resh_array(1)=nsiz; resh_array(2)=nsiz
        xmat(1:nsiz*nsiz) = daibj(ioff+1:ioff + nsiz*nsiz)
        amat = reshape(xmat,resh_array)
!
!  check asymmetry of obtained HPHP matrix
!
        x2 = dabs(qmasymr(amat,nsiz))
        if(x2.gt.1.0E-10) then
          write(*,*) 'Warning: Asymmetry of HP/HP matrix',
     &         ' in symmetry',desrep,' amounts to',x2
        endif
!
!  amat array has ph/ph structure and is transferred to LAPACK
!
!       write(*,*) 'Diagonalizing HP/HP in symmetry:',desrep
!       write(*,*) 'Size of amat:',size(amat)
!
!  specific allocations and settings for entering LAPACK routine
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
!  on screen. This we comment out because the energies are inaccurate
!  and the eigenvectors are only used as start vector for Davidson.
!

!       write(*,*)
!       write(*,*) 'Eigenvalues of states in symm.',desrep,repna(desrep)
!       write(*,*)
!       do i=1,nsiz
!         write(*,'(I5,2F22.12,2X,A4,X,A2)') i,evl(i),evl(i)*autoev,
!    &             repna(desrep),'@X'
!         write(*,'(5X,A)') '-----------------------------------------'
!         do j=1,nsiz
!           if(dabs(amat(j,i)).gt.0.1d0) then
!             ioff=jvoxx(desrep)
!             write(*,*) j,':  ',cfield(ioff+j),amat(j,i)
!           endif
!         enddo
!         write(*,*)
!       enddo
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

      enddo
      deallocate(cfield)
!|
!|
!|
!|__________________________________________________________


      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_X_COUSAT_R (HPHP_BLOCK,EPS,DESREP,LENCOUPLBL,
     &                             INTBUF,NBUFS,DOXEXT,WTHR,IOCH)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!  Computes (real) 2H2P/HP and SAT block of the exc. ADC Matrix
!  In the corresponding symmetry DESREP.
!
!---------------Calling variables--------------------------------------
!
      REAL*8                           :: HPHP_BLOCK(:)
      REAL*8                           :: EPS(*)
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
!
!---------------Local variables--------------------------------------
!
      INTEGER              :: n1,irep
      INTEGER              :: jrep,brep,jfun,bfun
      INTEGER              :: ialloc,mjcol,mjrow,nmvod

      INTEGER              :: mxno,mxnv,ioff
      INTEGER, dimension (:,:,:), allocatable  :: oolt
      INTEGER, dimension (:,:,:), allocatable  :: vvlt
!
      REAL*8, dimension(:), allocatable :: voooz
      REAL*8, dimension(:), allocatable :: vovvz
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

        SUBROUTINE XCOLMAKE(i1,i2,ra1,
     &                      i3,i4,i5,i6,
     &                      i7,i8,i9,
     &                      ra2,ra3,ia1,ia2)
          INTEGER                         :: i1,i2
          REAL*8, dimension(:)            :: ra1
          INTEGER                         :: i3,i4,i5,i6
          INTEGER                         :: i7,i8,i9
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

        SUBROUTINE XCOLFLSH(i1,i2,i3,i4,ra1,ia1,ia2)
          INTEGER                         :: i1,i2,i3,i4
          REAL*8, dimension(:)            :: ra1
          INTEGER, dimension(:)           :: ia1,ia2
        END SUBROUTINE XCOLFLSH

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

      end interface
!
!---------------Executable code--------------------------------------
!
!________________________________________________________
!|     prepare ADC matrix write buffers
!|
!|
      allocate(wbufr(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer could not be allocated.'
      allocate(wbufi(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer could not be allocated.'
      allocate(wbufj(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer could not be allocated.'
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
!  make room for VOOO and VOVV integrals
!  and read them into memory
!
      n1 = ivooot(nrep+1)
      allocate(voooz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #1 in BUILD_X_COUSAT_R'
      n1 = ivovvt(nrep+1)
      allocate(vovvz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #2 in BUILD_X_COUSAT_R'
      call getvooo(voooz)
      call rdvovv(vovvz)
!     write(*,*) 'breakpoint 1'
!
! determine actual length of coupling block
! and allocate column buffer
!
      NMVOD = MVO(desrep)
      CALL XCOLTEST(desrep,lencouplbl)
      ALLOCATE(BUFCOL(lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #3 in BUILD_X_COUSAT_R'
      ALLOCATE(BUFCOLEX(NMVOD + lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #4 in BUILD_X_COUSAT_R'
!     WRITE(IW,*) 'Block length of main block:',NMVOD
!     WRITE(IW,*) 'Block length of coupling block:',lencouplbl
!     WRITE(IW,*) 'Matrix dimension:',lencouplbl + NMVOD
!
!  according to the SRT16 convention, 
!  the columns also have DESREP symmetry. 
!  We loop over the columns in the same manner and construct
!  the columns of the 2H2P/HP coupling block. 
!  after construction we write out the column together with the
!  components of the main block.
!
      mjcol = 1
      do jrep = 1,nrep   !jrep,brep are fermionic
        brep = multb(desrep+nrep,jrep,1)
        do jfun = 1,NO(jrep)
        do bfun = 1,NV(brep)

          CALL XCOLMAKE(desrep,lencouplbl,bufcol,
     &                  jrep,brep,jfun,bfun,
     &                  mxno,mxnv,nrep,
     &                  voooz,vovvz,oolt,vvlt)
!
! combine main block column with coupling block column
!
          istart = JVOVO(desrep) + (mjcol-1)*nmvod
          bufcolex(1:nmvod)=hphp_block(istart+1:istart+nmvod)
          bufcolex(nmvod+1:nmvod+lencouplbl)=
     &        bufcol(1:lencouplbl)
!
! write out this column, but only the i>=j part!
! tridiagonal storage of matrix due to symmetry!
!
          CALL XCOLWRIT(IOCH,wthr,bufcolex,wbufr,wbufi,wbufj,
     &                  IBFP,NBUFS,LENCOUPLBL+NMVOD,mjcol,INTBUF)
          mjcol = mjcol + 1
        enddo
        enddo
      enddo

      IF((mjcol-1).ne.mvo(desrep))
     &  stop 'error 3 in build_x_cousat_r!'
!________________________________________________________
!|
!|
!| from here the main and coupling block are written to
!| the external file. We now build the SAT block.
!| ATT**  we can release, and should release *all*
!| memory which we do not need anymore.
!| Due to the large size of the matrices we can not afford
!| a four-indexed lookup table as it was used in the two-particle
!| propagator!
!|
!|_______________________________________________________

      deallocate(bufcolex)
      deallocate(voooz)
      deallocate(vovvz)

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
      write(*,*) '---- Number of buffers:',NBUFS
      write(*,*) '--------------------------------------------'

      deallocate(bufcol)
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
      SUBROUTINE XCOLTEST (DESREP,LENCOL)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
! Loop through VVOO structure with DESREP as destination irrep
! count corresponding number of elements and return it
! ATT!! One gets more elements than there are contained in the
! original VVOOTT array because there the overall direct product is
! always of A symmetry. We therefore get the number IVVOOTT(nrep+1)
! for DESREP=1. This reflecs the fact that the VVOO stream can be longer
! in the ADC matrix than the ordinary VVOOTT stream.
!
      INTEGER                   :: DESREP
      INTEGER, intent(INOUT)      :: LENCOL
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Executable code--------------------------------------
!
!     write(*,*) 'In xcoltest...'
      desrepkr=multb(desrep+nrep,1+nrep,2)
      if(multb(desrepkr+nrep,desrep+nrep,1).ne.1) then
        stop 'Kramers partner wrongly determined!'
      endif
!     write(*,*) 'Kramers partner of:',desrep,' is',desrepkr

      CDKL = 0
      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through

      CDREP = MULTB(DESREP+NREP,KLREP+NREP,1)
!     write(*,'(A,4I4)') 'µ: cdrep,klrep,cdrep*klrep*,desrepkr',
!    &  cdrep,klrep,MULTB(CDREP+NREP,KLREP+NREP,2),desrepkr
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
!
!  do the CD side...
!  CDREP, KLREP and DESREP are bosonic!
!
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                   CDKL = CDKL + 1
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
!
      LENCOL = CDKL
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLMAKE (desrep,lencouplbl,bufcol,
     &                     jrep,brep,J,B,
     &                     mxno,mxnv,dimrep,
     &                     voooz,vovvz,oolt,vvlt)
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
      INTEGER                          :: desrep,lencouplbl
      REAL*8, dimension(:)             :: bufcol
      INTEGER                          :: jrep,brep,J,B
      INTEGER                          :: mxno,mxnv,dimrep
      REAL*8, dimension(:)             :: voooz,vovvz
      INTEGER, dimension(:,:,:)        :: oolt,vvlt
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      real*8   ::  S
!
!---------------Executable code--------------------------------------
!
      desrepkr=multb(desrep+nrep,1+nrep,2)
      if(multb(desrepkr+nrep,desrep+nrep,1).ne.1) then
        stop 'Kramers partner wrongly determined!'
      endif

      bufcol = 0.0d0
      CDKL = 0
      if(dimrep.ne.nrep) STOP 'Internal error in XCOLMAKE!'
!
!  do the |>|> side...
!
      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through

        CDREP = MULTB(DESREP+NREP,KLREP+NREP,1)
!       write(*,'(A,4I4)') '¶: cdrep,klrep,cdrep*klrep*,desrep',
!    &   cdrep,klrep,MULTB(CDREP+NREP,KLREP+NREP,2),desrepkr
!
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
!
!  do the CD side...
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
! we have to treat four delta cases. they are processed sequentially
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
        ZKLREP=MULTB(KREP,LREP,1)
        if(zklrep.ne.klrep) stop 'internal error #44'
        IF(KREP.EQ.LREP) THEN
          if(oolt(k,l,krep).eq.0) stop'oolt 0!'
          IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
        ELSE
          IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
        ENDIF
        IOFF=IOFF * NVO(ZKLREP) + IVOOOT(ZKLREP)
        IOFF=IOFF + IIVO(DREP,JREP) + (J-1)*NV(DREP) + D
        if(ioff.gt.ivooot(nrep+1).or.ioff.lt.1) 
     &     stop 'intern. error #66!'
        bufcol(cdkl) = bufcol(cdkl) - voooz(ioff)
      endif
!|
!|  Case B:  < cj||kl > delta_bd, we calc a VOOO-offset.
!|
      if(BREP.eq.DREP.and.B.eq.D) then
        ZKLREP=MULTB(KREP,LREP,1)
        if(zklrep.ne.klrep) stop 'internal error #45'
        IF(KREP.EQ.LREP) THEN
          if(oolt(k,l,krep).eq.0) stop'oolt 0!'
          IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
        ELSE
          IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
        ENDIF
        IOFF=IOFF * NVO(ZKLREP) + IVOOOT(ZKLREP)
        IOFF=IOFF + IIVO(CREP,JREP) + (J-1)*NV(CREP) + C
        if(ioff.gt.ivooot(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #67!'
        bufcol(cdkl) = bufcol(cdkl) + voooz(ioff)
      endif
!|
!|  Case C:  < bl||cd >* delta_jk, we calc a VOVV offset
!|
      if(JREP.eq.KREP.and.J.eq.K) then
        ZCDREP=MULTB(CREP,DREP,1)
        if(zcdrep.ne.cdrep) stop 'internal error #46'
        IF(CREP.EQ.DREP) THEN
          if(vvlt(c,d,crep).eq.0) stop'vvlt 0!'
          IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
        ELSE
          IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
        ENDIF
        IOFF=IOFF * NVO(ZCDREP) + IVOVVT(ZCDREP)
        IOFF=IOFF + IIVO(BREP,LREP) + (L-1)*NV(BREP) + B
        if(ioff.gt.ivovvt(nrep+1)) stop 'intern. error #68!'
        bufcol(cdkl) = bufcol(cdkl) - vovvz(ioff)
      endif
!|
!|  Case D:  < bk||cd >* delta_jl, we calc a VOVV offset
!|
      if(jrep.eq.lrep.and.J.eq.L) then
        ZCDREP=MULTB(CREP,DREP,1)
        if(zcdrep.ne.cdrep) stop 'internal error #47'
        IF(CREP.EQ.DREP) THEN
          if(vvlt(c,d,crep).eq.0) stop'vvlt 0!'
          IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
        ELSE
          IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
        ENDIF
        IOFF=IOFF * NVO(ZCDREP) + IVOVVT(ZCDREP)
        IOFF=IOFF + IIVO(BREP,KREP) + (K-1)*NV(BREP) + B
        if(ioff.gt.ivovvt(nrep+1)) stop 'intern. error #69!'
        bufcol(cdkl) = bufcol(cdkl) + vovvz(ioff)
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
!  writes out a column of the ADC matrix held in bufferin to the
!  buckets provided by the caller. If the buckets are full, they 
!  are emptied to the file in IOCH.
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
!         write(*,'(A,2I4,F20.12)') 'Z-->',i,mjcol,S
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
      REAL*8                               :: EA,EB,EI,EJ
      REAL*8                               :: S
      INTEGER                              :: JDUMMY
!
!---------------Start execution--------------------------------------
!     write(*,*) 'XSATMAKE wthr:',WTHR
!     write(*,*) 'In xsatmake...'
      desrepkr=multb(desrep+nrep,1+nrep,2)
      if(multb(desrepkr+nrep,desrep+nrep,1).ne.1) then
        stop 'Kramers partner wrongly determined!'
      endif
!     write(*,*) 'Kramers partner of:',desrep,' is',desrepkr

!
! allocate space and read in integrals
!
      n1 = ivvvvtt(nrep+1)
      allocate(vvvvz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'VVVV buffer could not be allocated.'

      n1 = ioooott(nrep+1)
      allocate(ooooz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'OOOO buffer could not be allocated.'

      n1 = ivovo(nrep+1)
      allocate(vovoz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOVO buffer could not be allocated.'

      CALL RDVVVV(vvvvz)
      CALL GETOOOO(ooooz)
      CALL GETVOVO(vovoz)
!
! CDKL is the column counter, has column symmetry DESREP
! ABIJ is the row counter has row symmetry DESREP* = KRP(DESREP)
! SAT elements are calculated only for the lower triangle which
! means that ABIJ >= CDKL !
! At the moment we have to examine this condition by counting through
! for we do not construct a lookup table (too big!)
!
      CDKL=0    !column counter
      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,1)
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
      DO IJREP = 1,NREP     ! KLREP (bosonic!) loops through
      ABREP = MULTB(DESREP+NREP,IJREP+NREP,1)
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
!        write(*,'(A,2I4,F20.12)') 'Z-->',ABIJ + MVO(DESREP),
!    &                             CDKL + MVO(DESREP),S
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
