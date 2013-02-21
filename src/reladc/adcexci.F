!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!      *****************************************************
!      ****        D H F - A D C (2x)  M O D U L E      ****
!      ****           POLARIZATION PROPAGATOR           ****
!      *****************************************************
!
!                       AUTHOR: M.Pernpointner
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADCEXC_MAIN (IW,EPS)
 
      IMPLICIT INTEGER(A_Z)

      integer iw      ! I/O file handle for screen output
      real*8  eps(*)  ! array of spinor energies
!
!---------------Common Blocks--------------------------------------
!
#include  "adcinpt.inc"
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
!      COMMON/RDMATND/NBUFS
!
!---------------Local variables memory requirements------------------
!
! file names for the Lanczos
!     CHARACTER*6 NAME2,NAME3
!     PARAMETER(NAME2='ADC2PF',NAME3='ADCCNF')
!     CHARACTER*7 INFOFN
!     PARAMETER(INFOFN='MATINFO')
!     CHARACTER*5 NMSPECF
!     PARAMETER(NMSPECF='DSPEC')
! variables for the lanczos
!     REAL*8 RDUMMY
!     INTEGER IONIZLEVEL
!     INTEGER IRECL
      integer                                ::  irploop,i,j,k,l,m,n
      integer                                ::  desrep,lensatbl
      integer                                ::  ICOLLEN(MXREP)
      integer                                ::  icount,ndim
      real*8, allocatable                    ::  ciajb(:)
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
      IF(DOFULL) THEN
        WRITE(IW,*) 'Full diagonalizer additionally requested.'
        WRITE(IW,*) 'Remark: Only executed up to a Matrix dimension'
        WRITE(IW,*) '     ****  of 5000    ****'
        WRITE(IW,*) 'and is for debugging purposes only.'
      ELSE
        WRITE(IW,*) 'Only Lanzcos diagonalizer will be called.'
      ENDIF
      IF(ADCTHR.gt.0.0d0) THEN
        WRITE(IW,*) 'A threshold of',ADCTHR,' will be applied to'
        WRITE(IW,*) 'the PP ADC matrix.'
      ENDIF
!
! determine symmetries to be calculated
!
      n = 0
      do i = 1,nrep
        n = n + excreps(i)
      enddo
      if(n.eq.0) then
        write(iw,*) 'User has not given explicit symmetries.'
        write(iw,*) 'Program calculates the following EXC symmetries:'
        write(iw,*)
        icount=1
        do i=1,nrep
          excreps(icount) = i
          write(iw,*) '   No: ',icount,'  symmetry: ',i
          icount = icount + 1
        enddo
        write(iw,*)
      else
        write(iw,*) 'User has made explicit symmetry choices.'
        write(iw,*) 'Program calculates the following EXC symmetries:'
        write(iw,*)
        icount=1
        do i=1,nrep
          if(excreps(i).ne.0) then
            write(iw,*) '   No: ',icount,'  symmetry: ',excreps(i)
            icount = icount + 1
          endif
        enddo
        write(iw,*) 'Program calculates',icount,' symmetries.'
      endif
      write(iw,*)
! |
! |
! |
! |________________END INTRO SECTION__________________________
!
!
!  __________________  BUILD HP-HP BLOCK  ___________________________
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
        n = ivovo(nrep+1) * rcw
        allocate(ciajb(n))
        CALL BUILD_X_HPHP (CIAJB,EPS)
! |
! |
! |_____________  E_N_D__BUILD HP-HP BLOCK  _________________________

!  _________________________ HP/2H2P SYMMETRY LOOP
! |
! |
      symm_loop_hp_2h2p:   DO IRPLOOP = 1,NREP

        IF(excreps(IRPLOOP).EQ.0) CYCLE
        IF(excreps(IRPLOOP).GT.NREP) THEN
           WRITE(IW,*) 'Chosen Symmetry',excreps(IRPLOOP),
     &     ' does not exist and is skipped!'
           CYCLE
        ENDIF

        IF(CARITH) THEN
          CALL PST('Building complex hp/2h2p coupling block+')
        ELSE
          CALL PST('Building real hp/2h2p coupling block+')
        ENDIF
!
!  determine hp/hp block in chosen symmetry
!
        DESREP = excreps(IRPLOOP)
        ICOLLEN(DESREP) = 0
        IF(NVO (DESREP).GT.0) THEN   ! *** clarify !!
          WRITE(IW,'(A,I4)')  
     &       ' *** Found HP final states in symm. ',DESREP
          IF(CARITH) THEN
            CALL BUILD_X_2H2P_HP_C(EPS,DESREP,LENSATBL)
          ELSE
            CALL BUILD_X_2H2P_HP_R(EPS,DESREP,LENSATBL)
          ENDIF
          ICOLLEN(DESREP) = LENSATBL
        ENDIF


      ENDDO symm_loop_hp_2h2p
! |
! |
! |_____________________ END HP/2H2P SYMMETRY LOOP
!
!

!     INTBUF = 5*1024*1024  !output buffer is now 5 MWORDS (40 MByte)

!
!  __________________  BUILD 2H2P-2H2P (SAT) BLOCK  _________________
! |
! |   This is done over individual symmetries. The packed 2h2p/hp
! |   matrices are already on disk.
! |

!     NDIM=(IOOOOTT(NREP+1)+1)*RCW
!     CALL MEMOPT('REAL',BFOOOO,NDIM,CC,KFREE,LFREE,ALLOC_ME)
!     NDIM=(IVOVO(NREP+1)+1)*RCW
!     CALL MEMOPT('REAL',BFVOVO,NDIM,CC,KFREE,LFREE,ALLOC_ME)
!     NDIM=IO(NREP+1)*2
!     CALL MEMOPT('INTE',BSPTRO,NDIM,CC,KFREE,LFREE,ALLOC_ME)
!     NDIM=IO(NREP+1)
!     NDIM=NDIM*NDIM*NDIM*NDIM
!     CALL MEMOPT('INTE',BLKUPO,NDIM,CC,KFREE,LFREE,ALLOC_ME)
!     NDIM=IV(NREP+1)*IO(NREP+1)
!     NDIM=NDIM*NDIM
!     CALL MEMOPT('INTE',BLKUPV,NDIM,CC,KFREE,LFREE,ALLOC_ME)
!
!  initialize transposition arrays
!
!     CALL MAKESPTRA(CC(BSPTRO),IO(NREP+1))
!
!  initialize lookup tables for OOOO and VOVO
!
!     CALL MAKELOOKUPO(CC(BLKUPO),IO(NREP+1))
!     CALL MAKELOOKUPV(CC(BLKUPV),IO(NREP+1),IV(NREP+1))

!     CALL MEMCHK('before BUILD_SAT',CC,KFRSAV)

c  loop over the symmetries
c  att: the additional, symmetry-specific memory is allocated
c  within this loop and immediately released afterwards.
c  att2: it was observed that the ADC2 matrices get very large and
c  therefore a diagonalization is required immediately after the
c  SAT block construction!
c      
      symm_loop_sat_block: DO IRPLOOP = 1,NREP
        IF(CARITH) THEN
          CALL PST('Building 1. order complex satellite block+')
        ELSE
          CALL PST('Building 1. order real satellite block+')
        ENDIF
        IF(excreps(IRPLOOP).EQ.0) CYCLE
        IF(excreps(IRPLOOP).GT.NREP) THEN
           WRITE(IW,*) 'Chosen Symmetry',excreps(IRPLOOP),
     &     ' does not exist and is skipped!'
           CYCLE
        ENDIF

        DESREP = excreps(IRPLOOP)

        IF(NVO(DESREP).EQ.0) CYCLE    ! ****  clarify !!!!
        WRITE(IW,*) ' *** Found final states in symm.',DESREP

        LENSATBL = ICOLLEN(DESREP)
        NDIM = LENSATBL
!       CALL MEMOPT('REAL',BSATEV,NDIM,CC,KFREE,LFREE,ALLOC_ME)
!       NDIM = LENSATBL*4
!       CALL MEMOPT('INTE',BFIRCD,NDIM,CC,KFREE,LFREE,ALLOC_ME)
!       NDIM = INTBUF
!       CALL MEMOPT('INTE',BUFIOI,NDIM,CC,KFREE,LFREE,ALLOC_ME)
!       CALL MEMOPT('INTE',BUFIOJ,NDIM,CC,KFREE,LFREE,ALLOC_ME)
!       NDIM = INTBUF*RCW
!       CALL MEMOPT('REAL',BUFBUF,NDIM,CC,KFREE,LFREE,ALLOC_ME)
!       NDIM = IO(NREP+1)*2
!       CALL MEMOPT('INTE',BIROOO,NDIM,CC,KFREE,LFREE,ALLOC_ME)

        IF(CARITH) THEN
          CALL BUILD_X_SAT_C(LENSATBL,DESREP,ADCTHR,INTBUF)
        ELSE
          CALL BUILD_X_SAT_R(LENSATBL,DESREP,ADCTHR,INTBUF)
        ENDIF
!
!  call hp/2h2p configuration creator, file handle is available
!  afterwards
!
!       LADC = NOOT(DESREP) + LENSATBL
        CALL WCONDAT_X(ITAPADC,NAME3,IRECL,DESREP,LADC,
     &                 IO(NREP+1),CC(BIROOO))
        WRITE(IW,*) 'Record length as determined by WCONDAT_X:',IRECL
!
!  write general matrix control file and perform Lanczos diagonalization
!  for each symmetry immediately.
!
        OPEN(ITAPADC,FILE=INFOFN,FORM='FORMATTED',STATUS='UNKNOWN')

        WRITE(ITAPADC,'(I1)')  IONIZLEVEL        !IONIZLEVEL
        WRITE(ITAPADC,'(I2)')  DESREP            !KREP
        WRITE(ITAPADC,'(I1)')  RCW               !RCW
        N=2
        WRITE(ITAPADC,'(I1)')  N                 !IOLDNEW
        WRITE(ITAPADC,'(I12)') LADC              !LADC
        WRITE(ITAPADC,'(A6)')  NAME2             !FILEADC
        WRITE(ITAPADC,'(A6)')  '      '          !FILEDIA
        WRITE(ITAPADC,'(A6)')  NAME3             !FILECNF
        WRITE(ITAPADC,'(A5)')  NMSPECF           !Name of SPEC
        WRITE(ITAPADC,'(I9)')  IRECL             !IRECL of config
        WRITE(ITAPADC,'(I9)')  DIPITER           !NITER
        WRITE(ITAPADC,'(I9)')  NOOT(DESREP)      !2h space size
        WRITE(ITAPADC,'(I9)')  NBUFS             !NBUFS
        WRITE(ITAPADC,'(I9)')  DIPEIGV(DESREP)   !NEIGENV
        WRITE(ITAPADC,'(F10.1)') DIPPRNT         !EVLLIM

        CLOSE(ITAPADC)
!
!  call full diagonalizer for new matrix format if requested.
!
        IF(DOFULL) THEN
          IF(CARITH) THEN
            CALL FULLDIAC2(ITAPADC, NAME2, INTBUF, NBUFS, LADC, DESREP)
          ELSE
            CALL FULLDIAR2(ITAPADC, NAME2, INTBUF, NBUFS, LADC, DESREP)
          ENDIF
        ENDIF
!
!  call Lanczos diagonalizer.
!
        OPEN(ITAPADC,FILE=INFOFN,FORM='FORMATTED',STATUS='OLD')
        CALL DIAGEXE(CC,KFRSAV2,LFREE,ITAPADC,IW,DOINCORE)
        CLOSE(ITAPADC,STATUS='DELETE')
      ENDDO symm_loop_sat_block
! |
! |
! |
! |________________E_N_D______S_E_C_T_I_O_N__________________________

      RETURN
      END
! |
! |
! |
! |________________E_N_D___M A I N___P R O G R A M __________________

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_X_HPHP (CIAJB,EPS)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Computes HP/HP block of the Excitation ADC Matrix
!     (polarization propagator). 
!
!---------------Calling variables--------------------------------------
!
      REAL*8,  intent(in)                           :: CIAJB(:)
      REAL*8,  intent(in)                           :: EPS(:)
!
!---------------Common Blocks, symmetry arrays from CCSETI --------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
!
!---------------Local variables--------------------------------------
!
      real*8, allocatable                           :: sa_vvoott(:)
      real*8, allocatable                           :: sb_vvoott(:)
      real*8, allocatable                           :: buf1(:),buf2(:)
      real*8, allocatable                           :: buf3(:),buf4(:)
      real*8, allocatable                           :: buf5(:)
      real*8, allocatable                           :: bab(:)
      real*8, allocatable                           :: dij(:)
!
!---------------Executable code--------------------------------------
!
      IF(CARITH) THEN
        CALL PST('Building complex HP/HP block+')
      ELSE
        CALL PST('Building real HP/HP block+')
      ENDIF
!
!  we often need the v_ac,kl array which is V_ac,kl (VVOO)
!  divided by the energy denominator
!  e_a + e_c - e_k - e_l and the original V_ac,kl
!  we therefore calculate it once and store it in this subroutine.
!  the variables icX count the number of integrals. For the actual
!  storage requirements one has to multiply with RCW!
!
      ic1 = ivvoott(nrep+1)
      allocate(sa_vvoott(ic1*rcw))
      allocate(sb_vvoott(ic1*rcw))
      call getvvoo(sa_vvoott)
      call xcopy(ic1,sa_vvoott,1,sb_vvoott,1)
      call denomvvoo(eps,sa_vvoott,sa_vvoott)
!
! we have v_ac,kl in sa_vvoott and   (storage a...)
! V_ac,kl         in sb_vvoott       (storage b...)
!

!******************************************
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

      ic5 = j2vovo(nrep+1)
      allocate(buf3(ic5*rcw))  ! all arrays have to hold VOVO integrals
      allocate(buf4(ic5*rcw))  ! in MVO*MVO instead of NVO*NVO manner!
      allocate(buf5(ic5*rcw))
      buf3 = 0.0d0; buf4 = 0.0d0; buf5 = 0.0d0

      ic6 = ivovo(nrep+1)      ! size of the generic VOVO array in generic BKC
!******************************************
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
!  **att** tha array bab is packed as nv1*nv1(beta1),...
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
         n = mvv(irep)
         k = mvo(jrep)
         ioff2 = j2vvvo(jrep) * rcw + 1
         call xgemm ('N','T',m,n,k,a1,buf3(ioff1),m,buf4(ioff2),n,
     &               a0,buf5(ioff3),m)
         ioff1 = ioff1 + m * k * rcw
         ioff3 = ioff3 + m * n * rcw
      enddo
!
!  We have to backsort in order to arrive at
!  Y(1)_aj,bi order with <|<||>|>  BKC.
! do backsort   ATTention ***** possible problems in SRT 16 inverse mode !!!!!
! do backsort   ATTention ***** possible problems in SRT 16 inverse mode !!!!!
! do backsort   ATTention ***** possible problems in SRT 16 inverse mode !!!!!
!
!  ai,bj (buf5) ---->  Y(1)_aj,bi  (buf3)  with generic VOVO BKC.
!
! remark: also in the inverse mode the first buffer remains the source and the
! second one the destination. However, the source structure of the integrals
! reflects the destination structure in a regular sort.
!
      call srt16 (nrep,multb,lfa,ltr,nv,no,nv,no,mvo,jvovo,jjvo,jjvo,
     &            buf5,buf3)
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
!  fetch the C(1) contribution to hp/hp block, these are the VOVO integrals.
!
      call getvovo(ciajb)
!
!  add -0.5 *  A_aj,bi array to it (VOVO storage), this is the C(2C)
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
!  the ciajb array now contains the full hp/hp block
!
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BWEAVERR(bhphp,eps,bab,dij,ic6,ic2,ic4)
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
      real*8, intent(inout)       :: bhphp(ic6)
      real*8, intent(in)          :: eps(:)
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
!               if(jrp.ne.krp) call quit('weaver_error')
!               IAOFF = IOO(KRP) + (J-1)*NO(KRP) + K
!               IAOFF = IOO(KRP) + (J-1)*NO(KRP) + K
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
      complex*16, intent(inout)       :: bhphp(ic6)
      real*8, intent(in)          :: eps(:)
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
      SUBROUTINE BUILD_3H1P_2H_R (MAXOCCSP,IXOOO,IROOO,ORIJK,
     &                            VOOO1,VOOO2,COL,GIJKL,EPS,DESREP,
     &                            LEN_RKLM)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes 3H1P/2H block of the Double ionization ADC Matrix
C     This is the real number implementation.
c
C---------------Calling variables--------------------------------------
C
      INTEGER MAXOCCSP
      INTEGER IXOOO(MAXOCCSP,MAXOCCSP,MAXOCCSP)
      INTEGER IROOO(MAXOCCSP,2)
      INTEGER ORIJK(MAXOCCSP*MAXOCCSP*MAXOCCSP)
      REAL*8 VOOO1(*),VOOO2(*),COL(*),GIJKL(*),EPS(*)
      INTEGER DESREP,LEN_RKLM
c
c  MAXOCCSP is the total number of occupied spinors over all ireps.
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 AUX
      CHARACTER*6 NAME3,NAME4
      PARAMETER(NAME3='ADC2PE',NAME4='ADC2PP')
      CHARACTER*9 FNAME3,FNAME4
      REAL*8 ENERGY_I,ENERGY_J,ENERGY_K,ENERGY_R
C
C---------------Executable code--------------------------------------
c
c  desrep is the destination symmetry of the two-hole final state
c  the packed matrix for the coupling block just consists of a
c  sequential packed file which will be reorganized in the routine
c  for the satellite block and transformed to its final form.
c
c_______  INITIALIZE ARRAYS FOR THE COUPL. BLOCK ____________________
c                                                                    |
c                                                                    |
c
c  all occupied spinors are counted contiguously through the
c  irreps. Assign the corresponding irep to this spinor number.
c  in the code the spinors are accessed through this unique number.
c  this unique number can be interpreted as *absolute* spinor number
c  in contrast to the *relative* spinor number per irep.
c
      IX=1
      Do KRP=1,NREP
        Do K=1,NO(KRP)
          IROOO(IX,1)=K
          IROOO(IX,2)=KRP
          IX = IX + 1
        Enddo
      Enddo
c
c  set up triple index array for totally antisymmetric storage of 
c  the Gamma_rklm,ij
c  The gamma array is totally antisymmetric in klm!
c  This means that I>J>K must be maintained throughout and
c  no I=J, J=K or I=K combination must occur 
c
      IX=1
      DO K=1,MAXOCCSP
        DO J=1,MAXOCCSP
          DO I=1,MAXOCCSP
            IXOOO(I,J,K) = 0
            ORIJK(IX)=0
            IX = IX + 1
          ENDDO
        ENDDO
      ENDDO
c
c  initialize column offset array.
c  In due course also calculate the energies and use the COL(*) 
c  buffer for storage. Later on the COL(*) buffer holds the PHHH 
c  column.
c  For the chosen storage mode I,J,K with I>J>K we additionally
c  know the corresponding IREPS and can then access the one-particle
c  spinor energies accordingly.
c
      IX = 1
      IEX = 1
      DO K=1,MAXOCCSP
        KRP = IROOO(K,2)
        DO J=K+1,MAXOCCSP
          JRP = IROOO(J,2)
          JKRP = MULTB(JRP,KRP,1)
          DO I=J+1,MAXOCCSP
            IRP = IROOO(I,2)
            IJKRP = MULTB(IRP,JKRP+NREP,1)
            IXOOO(I,J,K)=IX
            IX = IX + 1
            RRP = MULTB(DESREP+NREP,IJKRP,2)
c
c  we have determined G_r such that G_r* x G_ijk = G_dest
c  attention: DESREP is bosonic !
c
            ORIJK(IX) = ORIJK(IX-1) + NV(RRP)
c
c  energies. from here I,J,K and IRP,JRP,KRP,RRP are known.
c  We loop over R for the different particle indices and construct the
c  energy array.
c
            ENERGY_I = EPS(IO(IRP) + IROOO(I,1))
            ENERGY_J = EPS(IO(JRP) + IROOO(J,1))
            ENERGY_K = EPS(IO(KRP) + IROOO(K,1))
            AUX = ENERGY_I + ENERGY_J + ENERGY_K
            DO 17 IXR = 1,NV(RRP)
              ENERGY_R = EPS(IO(NREP+1) + IV(RRP) + IXR)
              COL(IEX) = ENERGY_R - AUX
              IEX = IEX + 1
 17         CONTINUE

          ENDDO
        ENDDO
      ENDDO
      LEN_RKLM=ORIJK(IX)   
      WRITE(IW,'(4X,A,I3,A,I8)') 'Length of coupling block in symmetry',
     &      DESREP,':',LEN_RKLM
c                                                                    |
c____________________________________________________________________|
c_____________________________________________________________________
c                                                                    |
c            OPEN CONTROL/MATRIX/ENERGY FILES FOR THE ADC-2          |
c                                                                    |
      IF(DESREP.GT.9) THEN
        WRITE(FNAME3,'(A6,A1,I2)') NAME3,'.',DESREP
        WRITE(FNAME4,'(A6,A1,I2)') NAME4,'.',DESREP
      ELSE
        WRITE(FNAME3,'(A6,A2,I1)') NAME3,'.0',DESREP
        WRITE(FNAME4,'(A6,A2,I1)') NAME4,'.0',DESREP
      ENDIF

c     WRITE(IW,*) '2P-ADC block data (packed) in       ',FNAME4
c     WRITE(IW,*) '2P-ADC matrix satellite energies in ',FNAME3

      OPEN(ITAPADC+3,FILE=FNAME3,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAPADC+4,FILE=FNAME4,FORM='UNFORMATTED',STATUS='UNKNOWN')
c
c       write satellite block diagonal energies immediately
c
      WRITE(ITAPADC+3) (COL(IX),IX=1,LEN_RKLM)
      CLOSE(ITAPADC+3)
c                                                                    |
c____________________________________________________________________|
c
c  read VOOO integrals and sort to VOO,O:  (I K>L, J)
c
      CALL GETVOOO(VOOO1)
      CALL SRT6 (NREP,MULTB,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           VOOO1,VOOO2)
c
C *****************************************************
C ************  MASTER SYMMETRY LOOP OVER THE 2H SYMS *
C *****************************************************
C
c
c  the externally given symmetry is DESREP. Of course, there are
c  multiple combinations of gamma_j and gamma_i to fulfill
c  gamma_i x gamma_j = DESREP. They all have to be included and
c  form the final array in the genuine order below.
c  this genuine order also defines the order of columns in GIJKL
c
      GCOL=0
      DO 90 JXREP = 1, NREP
        IXREP = MULTB(JXREP,DESREP+NREP,2)
        IF (IXREP.LT.JXREP) GOTO 90
        DO 85 JX = 1, NO(JXREP)
          IXMIN = 1
          IF (IXREP.EQ.JXREP) IXMIN = JX + 1
          DO 80 IX = IXMIN, NO(IXREP)

c
c  the following initializations and updates *HAVE* to take place
c  for *EACH* new column, determined by Gamma_i, Gamma_j, IX and JX:
c
c  - GCOL:  increase column counter by one.
c  - CLEAR intermediate column COL(*). Contributions are added on top
c          therefore it must be set to zero for each new column. Once
c          a col is finished it is written out immediately.
c
            GCOL = GCOL + 1
            DO IY=1,LEN_RKLM
              COL(IY) = 0.0D0 
            ENDDO
c
c  from here we know the ireps IREP,JREP and the corresponding spinor
c  indices i and j. Now we rumble through the sorted VOOOs and distribute the
c  integrals
c_______________________________________________________________
c|
c|      internal VOOO Loop start
c|  
c|  first treat the delta_i? cases
c|
      OFFOVOO = LOVOOT(JXREP) + (JX-1)*NVOOT(JXREP)
      ICNT = 1
      DO 55 LMREP=1,NREP
        DO 50 MREP=1,NREP
          LREP = MULTB(MREP,LMREP+NREP,2)
           IF (LREP.LT.MREP) GOTO 50
           DO M=1,NO(MREP)
             LMIN = 1
             IF (LREP.EQ.MREP) LMIN = M + 1
             DO L = LMIN, NO(LREP)
               RREP=MULTB(JXREP,LMREP+NREP,2)

c
c  from here all the rlm Ireps and rlm indices are known
c  we can assign a destination offset.
c
c  *** Attention:  IX, JX, L and M are relative spinor numbers
c  but the IXOOO array needs absolute ones ! 
c  the R spinor index stays relative...
c
c  the rklm values are created by looping through the VOO,O
c  integrals only once. Successively the loop variables
c  L and M are interpreted as L,M/K,M and K,L !
c  
               DO 45 R=1,NV(RREP)

                 OFFS = (OFFOVOO + (ICNT-1) + R)
                 AUX = VOOO2(OFFS)
 
                 XXK = IO(IXREP) + IX
                 XXL = IO(LREP)  + L
                 XXM = IO(MREP)  + M
                 XX1 = IXOOO(XXK,XXL,XXM)
                 IF(XX1.NE.0) THEN
                   OFFD = (ORIJK(XX1) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - AUX
                 ENDIF
                   
                 XXK = IO(LREP)  + L
                 XXL = IO(IXREP) + IX
                 XXM = IO(MREP)  + M
                 XX2 = IXOOO(XXK,XXL,XXM)
                 IF(XX2.NE.0) THEN
                   OFFD = (ORIJK(XX2) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + AUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(MREP)  + M
                 XXM = IO(IXREP) + IX
                 XX3 = IXOOO(XXK,XXL,XXM)
                 IF(XX3.NE.0) THEN
                   OFFD = (ORIJK(XX3) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - AUX
                 ENDIF

 45            CONTINUE
   
               ICNT = ICNT + NV(RREP)

             ENDDO
           ENDDO
 50     CONTINUE
 55   CONTINUE
c|
c|  second treat the delta_j? case
c|
      OFFOVOO = LOVOOT(IXREP) + (IX-1)*NVOOT(IXREP)
      ICNT = 1
      DO 75 LMREP=1,NREP
        DO 70 MREP=1,NREP
          LREP = MULTB(MREP,LMREP+NREP,2)
           IF (LREP.LT.MREP) GOTO 70
           DO M=1,NO(MREP)
             LMIN = 1
             IF (LREP.EQ.MREP) LMIN = M + 1
             DO L = LMIN, NO(LREP)
               RREP=MULTB(IXREP,LMREP+NREP,2)
c
c  from here all the rlm Ireps and rlm indices are known
c  we can assign a destination offset
c
               DO 65 R=1,NV(RREP)

                 OFFS = (OFFOVOO + (ICNT-1) + R)
                 AUX = VOOO2(OFFS)

                 XXK = IO(JXREP) + JX
                 XXL = IO(LREP)  + L
                 XXM = IO(MREP)  + M
                 XX1 = IXOOO(XXK,XXL,XXM)
                 IF(XX1.NE.0) THEN
                   OFFD = (ORIJK(XX1) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + AUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(JXREP) + JX
                 XXM = IO(MREP)  + M
                 XX2 = IXOOO(XXK,XXL,XXM)
                 IF(XX2.NE.0) THEN
                   OFFD = (ORIJK(XX2) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - AUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(MREP)  + M
                 XXM = IO(JXREP) + JX
                 XX3 = IXOOO(XXK,XXL,XXM)
                 IF(XX3.NE.0) THEN
                   OFFD = (ORIJK(XX3) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + AUX
                 ENDIF

 65            CONTINUE

               ICNT = ICNT + NV(RREP)

             ENDDO
           ENDDO
 70     CONTINUE
 75   CONTINUE
c|
c|      internal VOOO Loop end
c|  
c|______________________________________________________________
c
c  from here one column of the coupling block matrix is finished
c  and ranges from 1 to LEN_RKLM
c  Next the corresponding 2h column
c  will be placed in front of it and then written to disk.
c_______________________________________________________________
c   FORM COMPLETE COLUMN AND WRITE TO DISK IN PACKED FORM.     |
c  The packed matrix storage implies column counting starting  |
c  from 1 and the next column is indicated by a row index of   |
c  -1. EOF is indicated by a row index of -255.                |
c                                                              |
          JOFF = NOOT(DESREP)
          GOFFSET = IOOOOTT(DESREP) + (GCOL-1)*JOFF
          DO JXX = LEN_RKLM,1,-1
            COL(JXX+JOFF) = COL(JXX)
          ENDDO
          DO JXX = 1,JOFF
            COL(JXX) = GIJKL(GOFFSET + JXX)
          ENDDO
          JXXM = LEN_RKLM + JOFF
          DO JXX = GCOL, JXXM
            IF(COL(JXX).NE.0.0D0) THEN
              WRITE(ITAPADC+4) JXX,COL(JXX)
            ENDIF
          ENDDO
          JXX = -1
          WRITE(ITAPADC+4) JXX,0.0d0
c                                                              |
c______________________________________________________________|
c
c
  80      CONTINUE
  85    CONTINUE
  90  CONTINUE
c
c  close adc matrix file (packed and unpacked)
c  write end of file descriptor.
c
      JXX = -255
      WRITE(ITAPADC+4) JXX,0.0d0
      CLOSE(ITAPADC+4)

      RETURN
      END
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_3H1P_2H_C(MAXOCCSP,IXOOO,IROOO,ORIJK,
     &                           VOOO1,VOOO2,COL,ECOL,GIJKL,
     &                           EPS,DESREP,LEN_RKLM)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes 3H1P/2H block of the Double ionization ADC Matrix
C     This is the complex number implementation.
C
C---------------Calling variables--------------------------------------
C
      INTEGER MAXOCCSP
      INTEGER IXOOO(MAXOCCSP,MAXOCCSP,MAXOCCSP)
      INTEGER IROOO(MAXOCCSP,2)
      INTEGER ORIJK(MAXOCCSP*MAXOCCSP*MAXOCCSP)
      COMPLEX*16 VOOO1(*),VOOO2(*),COL(*),GIJKL(*)
      REAL*8 ECOL(*),EPS(*)
      INTEGER DESREP,LEN_RKLM
c
c  MAXOCCSP is the total number of occupied spinors over all ireps.
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 AUX
      COMPLEX*16 CAUX
      CHARACTER*6 NAME3,NAME4
      PARAMETER(NAME3='ADC2PE',NAME4='ADC2PP')
      CHARACTER*9 FNAME3,FNAME4
      REAL*8 ENERGY_I,ENERGY_J,ENERGY_K,ENERGY_R
C
C---------------Executable code--------------------------------------
c
c  desrep is the destination symmetry of the two-hole final state
c  the packed matrix for the coupling block just consists of a
c  sequential packed file which will be reorganized in the routine
c  for the satellite block and transformed to its final form.
c
c_______  INITIALIZE ARRAYS FOR THE COUPL. BLOCK ____________________
c                                                                    |
c                                                                    |
c
c  all occupied spinors are counted contiguously through the
c  irreps. Assign the corresponding irep to this spinor number.
c  in the code the spinors are accessed through this unique number.
c  this unique number can be interpreted as *absolute* spinor number
c  in contrast to the *relative* spinor number.
c
      IX=1
      Do KRP=1,NREP
        Do K=1,NO(KRP)
          IROOO(IX,1)=K
          IROOO(IX,2)=KRP
          IX = IX + 1
        Enddo
      Enddo
c
c  set up triple index array for totally antisymmetric storage of 
c  the Gamma_rklm,ij
c  The gamma array is totally antisymmetric in klm!
c  This means that I>J>K must be maintained throughout!
c  no I=J, J=K or I=K combination must occur 
c
      IX=1
      DO K=1,MAXOCCSP
        DO J=1,MAXOCCSP
          DO I=1,MAXOCCSP
            IXOOO(I,J,K) = 0
            ORIJK(IX)=0
            IX = IX + 1
          ENDDO
        ENDDO
      ENDDO
c
c  initialize column offset array.
c  In due course also calculate the energies and use the COL(*) 
c  buffer for storage. Later on the COL(*) buffer holds the PHHH 
c  column.
c  For the chosen storage mode I,J,K with I>J>K we additionally
c  know the corresponding IREPS and can then access the one-particle
c  spinor energies accordingly.
c
      IX = 1
      IEX = 1
      DO K=1,MAXOCCSP
        KRP = IROOO(K,2)
        DO J=K+1,MAXOCCSP
          JRP = IROOO(J,2)
          JKRP = MULTB(JRP,KRP,1)
          DO I=J+1,MAXOCCSP
            IRP = IROOO(I,2)
            IJKRP = MULTB(IRP,JKRP+NREP,1)
            IXOOO(I,J,K)=IX
            IX = IX + 1
            RRP = MULTB(DESREP+NREP,IJKRP,2)
c
c  we have determined G_r such that G_r* x G_ijk = G_dest
c  attention: DESREP is bosonic !
c
            ORIJK(IX) = ORIJK(IX-1) + NV(RRP)
c
c  energies. from here I,J,K and IRP,JRP,KRP,RRP are known.
c  We loop over R for the different particle indices and construct the
c  energy array.
c
            ENERGY_I = EPS(IO(IRP) + IROOO(I,1))
            ENERGY_J = EPS(IO(JRP) + IROOO(J,1))
            ENERGY_K = EPS(IO(KRP) + IROOO(K,1))
            AUX = ENERGY_I + ENERGY_J + ENERGY_K
            DO 17 IXR = 1,NV(RRP)
              ENERGY_R = EPS(IO(NREP+1) + IV(RRP) + IXR)
              ECOL(IEX) = ENERGY_R - AUX
              IEX = IEX + 1
 17         CONTINUE

          ENDDO
        ENDDO
      ENDDO
      LEN_RKLM=ORIJK(IX)   
      WRITE(IW,'(4X,A,I3,A,I8)') 'Length of coupling block in symmetry',
     &      DESREP,':',LEN_RKLM
c                                                                    |
c____________________________________________________________________|
c_____________________________________________________________________
c                                                                    |
c            OPEN CONTROL/MATRIX/ENERGY FILES FOR THE ADC-2          |
c                                                                    |
      IF(DESREP.GT.9) THEN
        WRITE(FNAME3,'(A6,A1,I2)') NAME3,'.',DESREP
        WRITE(FNAME4,'(A6,A1,I2)') NAME4,'.',DESREP
      ELSE
        WRITE(FNAME3,'(A6,A2,I1)') NAME3,'.0',DESREP
        WRITE(FNAME4,'(A6,A2,I1)') NAME4,'.0',DESREP
      ENDIF

      OPEN(ITAPADC+3,FILE=FNAME3,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAPADC+4,FILE=FNAME4,FORM='UNFORMATTED',STATUS='UNKNOWN')
c
c write satellite block diagonal energies immediately
c
      WRITE(ITAPADC+3) (ECOL(IX),IX=1,LEN_RKLM)
      CLOSE(ITAPADC+3)
c                                                                    |
c____________________________________________________________________|
c
c  read VOOO integrals and sort to VOO,O:  (I K>L, J)
c
      CALL GETVOOO(VOOO1)
      CALL SRT6 (NREP,MULTB,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           VOOO1,VOOO2)
C
C *****************************************************
C ************  MASTER SYMMETRY LOOP OVER THE 2H SYMS *
C *****************************************************
C
c
c  the externally given symmetry is DESREP. Of course, there are
c  multiple combinations of gamma_j and gamma_i to fulfill
c  gamma_i x gamma_j = DESREP. They all have to be included and
c  form the final array in the genuine order below.
c  this genuine order also defines the order of columns in GIJKL
c
      GCOL=0
      DO 90 JXREP = 1, NREP
        IXREP = MULTB(JXREP,DESREP+NREP,2)
        IF (IXREP.LT.JXREP) GOTO 90
        DO 85 JX = 1, NO(JXREP)
          IXMIN = 1
          IF (IXREP.EQ.JXREP) IXMIN = JX + 1
          DO 80 IX = IXMIN, NO(IXREP)

c
c  the following initializations and updates *HAVE* to take place
c  for *EACH* new column, determined by Gamma_i, Gamma_j, IX and JX:
c
c  - GCOL:  increase column counter by one.
c  - CLEAR intermediate column COL(*). Contributions are added on top
c          therefore it must be set to zero for each new column. Once
c          a col is finished it is written out immediately.
c
            GCOL = GCOL + 1
            DO IY=1,LEN_RKLM
              COL(IY) = (0.0D0,0.0D0) 
            ENDDO
c
c  from here we know the ireps IREP,JREP and the corresponding spinor
c  indices i and j. Now we rumble through the sorted VOOOs and distribute the
c  integrals
c_______________________________________________________________
c|
c|      internal VOOO Loop start
c|  
c|  first treat the delta_i? cases
c|
      OFFOVOO = LOVOOT(JXREP) + (JX-1)*NVOOT(JXREP)
      ICNT = 1
      DO 55 LMREP=1,NREP
        DO 50 MREP=1,NREP
          LREP = MULTB(MREP,LMREP+NREP,2)
           IF (LREP.LT.MREP) GOTO 50
           DO M=1,NO(MREP)
             LMIN = 1
             IF (LREP.EQ.MREP) LMIN = M + 1
             DO L = LMIN, NO(LREP)
               RREP=MULTB(JXREP,LMREP+NREP,2)

c       write(iw,'(A,4I3)') 'In LM LOOP:',LREP,MREP,L,M
c
c  from here all the rlm Ireps and rlm indices are known
c  we can assign a destination offset.
c
c  *** Attention:  IX, JX, L and M are relative spinor numbers
c  but the IXOOO array needs absolute ones ! 
c  the R spinor index stays relative...
c
c  the rklm values are created by looping through the VOO,O
c  integrals only once. Successively the loop variables
c  L and M are interpreted as L,M/K,M and K,L !
c  
               DO 45 R=1,NV(RREP)

                 OFFS = (OFFOVOO + (ICNT-1) + R)
                 CAUX = VOOO2(OFFS)
 
                 XXK = IO(IXREP) + IX
                 XXL = IO(LREP)  + L
                 XXM = IO(MREP)  + M
                 XX1 = IXOOO(XXK,XXL,XXM)
                 IF(XX1.NE.0) THEN
                   OFFD = (ORIJK(XX1) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - CAUX
                 ENDIF
                   
                 XXK = IO(LREP)  + L
                 XXL = IO(IXREP) + IX
                 XXM = IO(MREP)  + M
                 XX2 = IXOOO(XXK,XXL,XXM)
                 IF(XX2.NE.0) THEN
                   OFFD = (ORIJK(XX2) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + CAUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(MREP)  + M
                 XXM = IO(IXREP) + IX
                 XX3 = IXOOO(XXK,XXL,XXM)
                 IF(XX3.NE.0) THEN
                   OFFD = (ORIJK(XX3) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - CAUX
                 ENDIF

 45            CONTINUE
   
               ICNT = ICNT + NV(RREP)

             ENDDO
           ENDDO
 50     CONTINUE
 55   CONTINUE
c|
c|  second treat the delta_j? case
c|
      OFFOVOO = LOVOOT(IXREP) + (IX-1)*NVOOT(IXREP)
      ICNT = 1
      DO 75 LMREP=1,NREP
        DO 70 MREP=1,NREP
          LREP = MULTB(MREP,LMREP+NREP,2)
           IF (LREP.LT.MREP) GOTO 70
           DO M=1,NO(MREP)
             LMIN = 1
             IF (LREP.EQ.MREP) LMIN = M + 1
             DO L = LMIN, NO(LREP)
               RREP=MULTB(IXREP,LMREP+NREP,2)
c
c  from here all the rlm Ireps and rlm indices are known
c  we can assign a destination offset
c
               DO 65 R=1,NV(RREP)

                 OFFS = (OFFOVOO + (ICNT-1) + R)
                 CAUX = VOOO2(OFFS)

                 XXK = IO(JXREP) + JX
                 XXL = IO(LREP)  + L
                 XXM = IO(MREP)  + M
                 XX1 = IXOOO(XXK,XXL,XXM)
                 IF(XX1.NE.0) THEN
                   OFFD = (ORIJK(XX1) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + CAUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(JXREP) + JX
                 XXM = IO(MREP)  + M
                 XX2 = IXOOO(XXK,XXL,XXM)
                 IF(XX2.NE.0) THEN
                   OFFD = (ORIJK(XX2) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - CAUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(MREP)  + M
                 XXM = IO(JXREP) + JX
                 XX3 = IXOOO(XXK,XXL,XXM)
                 IF(XX3.NE.0) THEN
                   OFFD = (ORIJK(XX3) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + CAUX
                 ENDIF

 65            CONTINUE

               ICNT = ICNT + NV(RREP)

             ENDDO
           ENDDO
 70     CONTINUE
 75   CONTINUE
c|
c|      internal VOOO Loop end
c|  
c|______________________________________________________________
c
c  from here one column of the coupling block matrix is finished
c  and ranges from 1 to LEN_RKLM. Since the c.c. was created according
c  to the formulation of the equations we have
c  to complex conjugate the corresponding vector in the coupling block.
c  the remaining 2h/2h part must not be touched !

          DO JXX = 1,LEN_RKLM
            COL(JXX) = DCONJG(COL(JXX))
          ENDDO

c  Next the corresponding 2h column
c  will be placed in front of it and then written to disk.
c_______________________________________________________________
c   FORM COMPLETE COLUMN AND WRITE TO DISK IN PACKED FORM.     |
c  The packed matrix storage implies column counting starting  |
c  from 1 and the next column is indicated by a row index of   |
c  -1. EOF is indicated by a row index of -255.                |
c                                                              |
          JOFF = NOOT(DESREP)
          GOFFSET = IOOOOTT(DESREP) + (GCOL-1)*JOFF
          DO JXX = LEN_RKLM,1,-1
            COL(JXX+JOFF) = COL(JXX)
          ENDDO
          DO JXX = 1,JOFF
            COL(JXX) = GIJKL(GOFFSET + JXX)
          ENDDO
          JXXM = LEN_RKLM + JOFF
          DO JXX = GCOL, JXXM
            IF(COL(JXX).NE.0.0D0) THEN
              WRITE(ITAPADC+4) JXX,COL(JXX)
            ENDIF
          ENDDO
          JXX = -1
          WRITE(ITAPADC+4) JXX,(0.0d0,0.0d0)
c                                                              |
c______________________________________________________________|
c
c
  80      CONTINUE
  85    CONTINUE
  90  CONTINUE
c
c  close adc matrix file (packed and unpacked)
c  write end of file descriptor.
c
      JXX = -255
      WRITE(ITAPADC+4) JXX,(0.0d0,0.0d0)
      CLOSE(ITAPADC+4)

      RETURN
      END
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_SAT_R(OOOO,VOVO,IROOO,MAXOCCSP,MAXVIRSP,
     &                       SATEV,IRCD,IRCDLEN,LTO,LTV,DESREP,
     &                       ADCTHR,IOI,IOJ,BUF,INTBUF)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes 3H1P/3H1P block of the Double ionization ADC Matrix
C     real number implementation. Attention: row/column matrix storage
C     happens via packed integers ! One can safely assume that the
C     matrix dimension of the ADC matrix never will exceed 2^31-1 !!
C
C---------------Calling variables--------------------------------------
C
      REAL*8 OOOO(0:*),VOVO(0:*)
      INTEGER IROOO(MAXOCCSP,2),MAXOCCSP,MAXVIRSP
      INTEGER IRCD(4,IRCDLEN),IRCDLEN,DESREP
      INTEGER LTO(MAXOCCSP,MAXOCCSP,MAXOCCSP,MAXOCCSP)
      INTEGER LTV(MAXVIRSP,MAXOCCSP,MAXVIRSP,MAXOCCSP)
      REAL*8 SATEV(IRCDLEN)
      REAL*8 ADCTHR
c for the packed matrix output file
      INTEGER INTBUF
      INTEGER IOI(INTBUF),IOJ(INTBUF)
      REAL*8 BUF(INTBUF)
c
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/files.inc"

      COMMON/RDMATND/NBUFS
C
C---------------Local variables--------------------------------------
C
      REAL*8 G
      INTEGER IB,NBUFS,JDUMMY,NTOTELE
      CHARACTER*6 NAME2,NAME3,NAME4
      PARAMETER(NAME2='ADC2PE',NAME3='ADC2PP',NAME4='ADC2PF')
      CHARACTER*9 FNAME2,FNAME3,FNAME4
      INTEGER R,I,J,K,S,L,M,N
      CHARACTER*44 FIELD

      REAL*8 SPARS,XNZEL,XTOTEL
c
c  meaning of the file names:
c  ADC2PC:  control file with the dimensions and file infos
c  ADC2PE:  contains the 3h1p/3h1p diagonal energies
c  ADC2PP:  packed matrix of the 3h1p/2h block. contains also diags
c  ADC2PF:  full packed ADC matrix in the conventional nonzero element
c           storage. This file is created in this routine. Important:
c           In contrast to the one-particle case the diagonal elements
c           are already included in this matrix !
c  ADCDAT:  configuration information for the 2h and 3h1p vector
C
C---------------Executable code--------------------------------------
c
c
c  say hello
c

      CALL PST('Building real SAT block+')
      WRITE(IW,*) 'Threshold for ADC matrix:',ADCTHR

c  create master index table now with *absolute* spinor indices!
c  This index table assigns to a row or column number within the
c  satellite block a quadruple of absolute spinor indices in order
c  to enable the treatment of all the delta conditions.
c  The counting starts with 1 and means the first entry of the
c  3h1p/3h1p block not the start of the total ADC matrix!
c
c  take care with the conjuga routine! This routine *really* assumes a 
c  *complex argument*. It must not be applied to real arrays with the
c  implicit assumption that it will not affect anything.
c

      IXTOT = 1
      DO K=1,MAXOCCSP
        KRP = IROOO(K,2)
        DO J=K+1,MAXOCCSP
          JRP = IROOO(J,2)
          JKRP = MULTB(JRP,KRP,1)
          DO I=J+1,MAXOCCSP
            IRP = IROOO(I,2)
            IJKRP = MULTB(IRP,JKRP+NREP,1)
            RRP = MULTB(DESREP+NREP,IJKRP,2)
            DO R = 1,NV(RRP)
              IRCD(1,IXTOT) = IV(RRP) + R
              IRCD(2,IXTOT) = I
              IRCD(3,IXTOT) = J
              IRCD(4,IXTOT) = K
              IXTOT = IXTOT + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF( (IXTOT-1) .NE. IRCDLEN) stop 'internal error in BUILD_SAT!'
c
c  read in OOOO and VOVO integrals and cc the OOOO integrals
c
      CALL GETOOOO(OOOO(1))
      CALL GETVOVO(VOVO(1))
c
c  insert the neutral element for addition
c
      OOOO(0)=0.0d0
      VOVO(0)=0.0d0
c
c  create symmetry-specific filenames for 
c  sat-diag eigenvalues, for packed 3h1p/2h matrix 
c  ATT: final packed ADC matrix now only described by generic name!
c  These **huge** files will be deleted immediately after the
c  diagonalization.
c
      IF(DESREP.GT.9) THEN
        WRITE(FNAME2,'(A6,A1,I2)') NAME2,'.',DESREP
        WRITE(FNAME3,'(A6,A1,I2)') NAME3,'.',DESREP
      ELSE
        WRITE(FNAME2,'(A6,A2,I1)') NAME2,'.0',DESREP
        WRITE(FNAME3,'(A6,A2,I1)') NAME3,'.0',DESREP
      ENDIF

      WRITE(IW,*) '2P-ADC matrix satellite energies in ',FNAME2
      WRITE(IW,*) '2P-ADC 3h1p/2h (packed) data in     ',FNAME3
      WRITE(IW,*) '2P-ADC final packed ADC matrix in   ',NAME4

      ITAP2 = ITAPADC+2
      ITAP3 = ITAPADC+3
      ITAP4 = ITAPADC+4

      OPEN(ITAP2,FILE=FNAME2,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAP3,FILE=FNAME3,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAP4,FILE=NAME4, FORM='UNFORMATTED',STATUS='UNKNOWN')
c
c  read packed 3h1p/2h packed matrix created previously and transform
c  it to final format where only the nonzero numbers together with
c  row and column indices are stored.
c
      IB = 0
      NBUFS = 0
      JDUMMY = 0
      NTOTELE = 0
c
c  the column switching is implicitly done. Only the row index was
c  stored. The Column and row indices are absolute, no offset is
c  to be added in this case.
c
c
c  *ATT* *ATT* *ATT*  DSYEV yields eigenvalues sorted from the lowest
c  to the highest. Lanczos is not showing this behavior.
c  It converges fastest to the extremum eigenvalue which would be in our
c  case the highest possible two-hole final state and not the desired
c  lowest one!! Writing out the negative of the ADC matrix leads to the
c  desired behavior for Lanczos!
c
      SATCOL=1
 177  CONTINUE
      READ(ITAP3) SATROW,G
      IF(SATROW.EQ.-1) THEN
        SATCOL = SATCOL+1
        GOTO 177
      ENDIF
      IF(SATROW.EQ.-255) GOTO 180
      IF(G.EQ.0.0D0) STOP 'Read zero element in BUILD_SAT_R'
      IF(SATROW.LT.SATCOL) 
     &     STOP 'ROW/COL inconsistency in BUILD_SAT_R'
c
c  write matrix element only if greater than ADCTHR
c
      IF(ABS(G).GT.ADCTHR) THEN
        NTOTELE = NTOTELE + 1
        IB=IB+1
c
        BUF(IB)=G
        IOI(IB)=SATROW
        IOJ(IB)=SATCOL
        IF(IB.EQ.INTBUF) THEN
          WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &                 (IOI(IXX),IXX=1,INTBUF),
     &                 (IOJ(IXX),IXX=1,INTBUF),
     &                 INTBUF,JDUMMY
          IB = 0
          NBUFS = NBUFS + 1
        ENDIF
      ENDIF
      GOTO 177   !next element
 180  CONTINUE
      CLOSE(ITAP3,STATUS='DELETE')

c
c  from here the output buffer may remain partially filled
c  because more elements are about to come down below.
c

c
c  read in diagonal satellite eigenvalues (SATEV)
c
      READ(ITAP2) (SATEV(IZ),IZ=1,IRCDLEN)
      CLOSE(ITAP2,STATUS='DELETE')
c
c  loop through the columns and check each entry
c  The index conventions are as in the accomp. material
c
      SATOFF = NOOT(DESREP)
      DO 401 SATCOL=1,IRCDLEN

        S   = IRCD(1,SATCOL)
        L   = IRCD(2,SATCOL)
        M   = IRCD(3,SATCOL)
        N   = IRCD(4,SATCOL)

        DO 400 SATROW=SATCOL,IRCDLEN

          R   = IRCD(1,SATROW)
          I   = IRCD(2,SATROW)
          J   = IRCD(3,SATROW)
          K   = IRCD(4,SATROW)
c
c  initialize column entry
c
          G=0.0d0
c
c  start the delta blocks
c
c **** D_rs block
c
          IF(R.EQ.S) THEN

            IF(I.EQ.L) THEN
              G = G + OOOO(LTO(j,k,m,n))
            ENDIF
            IF(I.EQ.M) THEN
              G = G - OOOO(LTO(j,k,l,n))
            ENDIF
            IF(I.EQ.N) THEN
              G = G + OOOO(LTO(j,k,l,m))
            ENDIF

            IF(J.EQ.L) THEN
              G = G - OOOO(LTO(i,k,m,n))
            ENDIF
            IF(J.EQ.M) THEN
              G = G + OOOO(LTO(i,k,l,n))
            ENDIF
            IF(J.EQ.N) THEN
              G = G - OOOO(LTO(i,k,l,m))
            ENDIF

            IF(K.EQ.L) THEN
              G = G + OOOO(LTO(i,j,m,n))
            ENDIF
            IF(K.EQ.M) THEN
              G = G - OOOO(LTO(i,j,l,n))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + OOOO(LTO(i,j,l,m))
            ENDIF

          ENDIF
c
c **** D_occ,occ block, alpha
c
          IF(J.EQ.M) THEN
            IF(K.EQ.N) THEN
              G = G - VOVO(LTV(r,l,s,i))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, beta
c
          IF(J.EQ.L) THEN
            IF(K.EQ.M) THEN
              G = G - VOVO(LTV(r,n,s,i))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + VOVO(LTV(r,m,s,i))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, gamma
c
          IF(I.EQ.M) THEN
            IF(J.EQ.N) THEN
              G = G - VOVO(LTV(r,l,s,k))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + VOVO(LTV(r,l,s,j))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, delta
c
          IF(I.EQ.L) THEN
            IF(J.EQ.M) THEN
              G = G - VOVO(LTV(r,n,s,k))
            ENDIF
            IF(J.EQ.N) THEN
              G = G + VOVO(LTV(r,m,s,k))
            ENDIF
            IF(K.EQ.M) THEN
              G = G + VOVO(LTV(r,n,s,j))
            ENDIF
            IF(K.EQ.N) THEN
              G = G - VOVO(LTV(r,m,s,j))
            ENDIF
          ENDIF
c
c check if diagonal entry has to be added
c
          IF(SATROW.EQ.SATCOL) THEN
            G = G + SATEV(SATROW)
          ENDIF
c
c fill this element in the output buffer and flush, if necessary.
c the buffer counter is set outside the row/column loop
c reminder: satrow and satcol count from 1 to length of sat block,
c but the matrix starts with 2h2h block ==> add offset to these
c variables.
c
          IF(ABS(G).GT.ADCTHR) THEN
            NTOTELE = NTOTELE + 1
            IB=IB+1
c
            BUF(IB)=G
            IOI(IB)=SATROW + SATOFF
            IOJ(IB)=SATCOL + SATOFF
          
            IF(IB.EQ.INTBUF) THEN
              WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &                     (IOI(IXX),IXX=1,INTBUF),
     &                     (IOJ(IXX),IXX=1,INTBUF),
     &                     INTBUF,JDUMMY
              IB = 0
              NBUFS = NBUFS + 1
            ENDIF
          ENDIF


 400    CONTINUE
 401  CONTINUE
c
c  write out *remaining* elements of the file buffer
c  and count all nonzero elements in the ADC matrix.
c
      IF(IB.GT.0) THEN
        NBUFS = NBUFS + 1
        WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &               (IOI(IXX),IXX=1,INTBUF),
     &               (IOJ(IXX),IXX=1,INTBUF),
     &               IB,JDUMMY
      ENDIF
      CLOSE(ITAP4)
      WRITE(IW,*)
      WRITE(IW,*) 'Number of elements above threshold: ',NTOTELE
      WRITE(IW,*) 'Number of buffers in ADC matrix file   : ',NBUFS
      WRITE(IW,*) 'Dimension of ADC matrix in this symmetry:',
     6            IRCDLEN+SATOFF
      XNZEL=DBLE(NTOTELE)
      XTOTEL=DBLE(IRCDLEN+SATOFF)
      XTOTEL=XTOTEL*XTOTEL
      SPARS=(XNZEL/XTOTEL)*100.0D0
      SPARS=100.0D0 - SPARS
      WRITE(IW,'(A,F5.2,A)') ' Sparsity of ADC matrix:',SPARS,' %'
      WRITE(IW,*)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_SAT_C(OOOO,VOVO,IROOO,MAXOCCSP,MAXVIRSP,
     &                       SATEV,IRCD,IRCDLEN,LTO,LTV,DESREP,
     &                       ADCTHR,IOI,IOJ,BUF,INTBUF)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes 3H1P/3H1P block of the Double ionization ADC Matrix
C     real number implementation. Attention: row/column matrix storage
C     happens via packed integers ! One can safely assume that the
C     matrix dimension of the ADC matrix never will exceed 2^31-1 !!
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 OOOO(0:*),VOVO(0:*)
      INTEGER IROOO(MAXOCCSP,2),MAXOCCSP,MAXVIRSP
      INTEGER IRCD(4,IRCDLEN),IRCDLEN,DESREP
      INTEGER LTO(MAXOCCSP,MAXOCCSP,MAXOCCSP,MAXOCCSP)
      INTEGER LTV(MAXVIRSP,MAXOCCSP,MAXVIRSP,MAXOCCSP)
      REAL*8 SATEV(IRCDLEN)
      REAL*8 ADCTHR
c for the packed matrix output file
      INTEGER INTBUF
      INTEGER IOI(INTBUF),IOJ(INTBUF)
      COMPLEX*16 BUF(INTBUF)
c
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/files.inc"

      COMMON/RDMATND/NBUFS
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 G
      INTEGER IB,NBUFS,JDUMMY,NTOTELE
      CHARACTER*6 NAME2,NAME3,NAME4
      PARAMETER(NAME2='ADC2PE',NAME3='ADC2PP',NAME4='ADC2PF')
      CHARACTER*9 FNAME2,FNAME3,FNAME4
      INTEGER R,I,J,K,S,L,M,N
      CHARACTER*44 FIELD

      REAL*8 SPARS,XNZEL,XTOTEL
c
c  meaning of the file names:
c  ADC2PE:  contains the 3h1p/3h1p diagonal energies
c  ADC2PP:  packed matrix of the 3h1p/2h block. contains also diags
c  ADC2PF:  full packed ADC matrix in the conventional nonzero element
c           storage. This file is created in this routine. Important:
c           In contrast to the one-particle case the diagonal elements
c           are already included in this matrix !
c  ADCDAT:  configuration information for the 2h and 3h1p vector
C
C---------------Executable code--------------------------------------
c
c
c  say hello
c

      CALL PST('Building complex SAT block+')
      WRITE(IW,*) 'Threshold for ADC matrix:',ADCTHR

c  create master index table now with *absolute* spinor indices!
c  This index table assigns to a row or column number within the
c  satellite block a quadruple of absolute spinor indices in order
c  to enable the treatment of all the delta conditions.
c  The counting starts with 1 and means the first entry of the
c  3h1p/3h1p block not the start of the total ADC matrix!
c
c  take care with the conjuga routine! This routine *really* assumes a 
c  *complex argument*. It must not be applied to real arrays with the
c  implicit assumption that it will not affect anything.
c

      IXTOT = 1
      DO K=1,MAXOCCSP
        KRP = IROOO(K,2)
        DO J=K+1,MAXOCCSP
          JRP = IROOO(J,2)
          JKRP = MULTB(JRP,KRP,1)
          DO I=J+1,MAXOCCSP
            IRP = IROOO(I,2)
            IJKRP = MULTB(IRP,JKRP+NREP,1)
            RRP = MULTB(DESREP+NREP,IJKRP,2)
            DO R = 1,NV(RRP)
              IRCD(1,IXTOT) = IV(RRP) + R
              IRCD(2,IXTOT) = I
              IRCD(3,IXTOT) = J
              IRCD(4,IXTOT) = K
              IXTOT = IXTOT + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF( (IXTOT-1) .NE. IRCDLEN) stop 'internal error in BUILD_SAT!'
c
c  read in OOOO and VOVO integrals and cc the OOOO integrals
c  according to the equations on paper the VOVOs must not be
c  changed.
c
      CALL GETOOOO(OOOO(1))
      CALL GETVOVO(VOVO(1))
      DO K=1,IOOOOTT(NREP+1)
        OOOO(K)=DCONJG(OOOO(K))
      ENDDO
c
c  insert the neutral element for addition
c
      OOOO(0)=(0.0d0,0.0d0)
      VOVO(0)=(0.0d0,0.0d0)
c
c  create symmetry-specific filenames for
c  sat-diag eigenvalues, for packed 3h1p/2h matrix
c  ATT: final packed ADC matrix now only described by generic name!
c  These **huge** files will be deleted immediately after the
c  diagonalization.
c
      IF(DESREP.GT.9) THEN
        WRITE(FNAME2,'(A6,A1,I2)') NAME2,'.',DESREP
        WRITE(FNAME3,'(A6,A1,I2)') NAME3,'.',DESREP
      ELSE
        WRITE(FNAME2,'(A6,A2,I1)') NAME2,'.0',DESREP
        WRITE(FNAME3,'(A6,A2,I1)') NAME3,'.0',DESREP
      ENDIF

      WRITE(IW,*) '2P-ADC matrix satellite energies in ',FNAME2
      WRITE(IW,*) '2P-ADC 3h1p/2h (packed) data in     ',FNAME3
      WRITE(IW,*) '2P-ADC final packed ADC matrix in   ',NAME4

      ITAP2 = ITAPADC+2
      ITAP3 = ITAPADC+3
      ITAP4 = ITAPADC+4

      OPEN(ITAP2,FILE=FNAME2,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAP3,FILE=FNAME3,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAP4,FILE=NAME4, FORM='UNFORMATTED',STATUS='UNKNOWN')
c
c  read packed 3h1p/2h packed matrix created previously and transform
c  it to final format where only the nonzero numbers together with
c  row and column indices are stored.
c
      IB = 0
      NBUFS = 0
      JDUMMY = 0
      NTOTELE = 0
c
c  the column switching is implicitly done. Only the row index was
c  stored. The Column and row indices are absolute, no offset is
c  to be added in this case.
c
c
c  *ATT* *ATT* *ATT*  DSYEV yields eigenvalues sorted from the lowest
c  to the highest. Lanczos is not showing this behavior. It converges to
c  the highest eigenvalue first. This is not the behavior of Lanczos!
c  It converges fastest to the extremum eigenvalue which would be in our
c  case the highest possible two-hole final state and not the desired
c  lowest one!! Writing out the negative of the ADC matrix leads to the
c  desired behavior for Lanczos!
c
      SATCOL=1
 177  CONTINUE
      READ(ITAP3) SATROW,G
      IF(SATROW.EQ.-1) THEN
        SATCOL = SATCOL+1
        GOTO 177
      ENDIF
      IF(SATROW.EQ.-255) GOTO 180
      IF(G.EQ.(0.0D0,0.0D0)) STOP 'Read zero element in BUILD_SAT_C'
      IF(SATROW.LT.SATCOL) 
     &     STOP 'ROW/COL inconsistency in BUILD_SAT_C!'
      IF(ABS(G).GT.ADCTHR) THEN
        NTOTELE = NTOTELE + 1
        IB=IB+1
        BUF(IB)=G
        IOI(IB)=SATROW
        IOJ(IB)=SATCOL
        IF(IB.EQ.INTBUF) THEN
          WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &                 (IOI(IXX),IXX=1,INTBUF),
     &                 (IOJ(IXX),IXX=1,INTBUF),
     &                 INTBUF,JDUMMY
          IB = 0
          NBUFS = NBUFS + 1
        ENDIF
      ENDIF
      GOTO 177   !next element
 180  CONTINUE
      CLOSE(ITAP3,STATUS='DELETE')

c
c  from here the output buffer may remain partially filled
c  because more elements are about to come down below.
c

c
c  read in diagonal satellite eigenvalues (SATEV)
c
      READ(ITAP2) (SATEV(IZ),IZ=1,IRCDLEN)
      CLOSE(ITAP2,STATUS='DELETE')
c
c  loop through the columns and check each entry
c  The index conventions are as in the accomp. material
c
      SATOFF = NOOT(DESREP)
      DO 401 SATCOL=1,IRCDLEN

        S   = IRCD(1,SATCOL)
        L   = IRCD(2,SATCOL)
        M   = IRCD(3,SATCOL)
        N   = IRCD(4,SATCOL)

        DO 400 SATROW=SATCOL,IRCDLEN

          R   = IRCD(1,SATROW)
          I   = IRCD(2,SATROW)
          J   = IRCD(3,SATROW)
          K   = IRCD(4,SATROW)
c
c  initialize column entry
c
          G=(0.0d0,0.0d0)
c
c  start the delta blocks
c
c **** D_rs block
c
          IF(R.EQ.S) THEN

            IF(I.EQ.L) THEN
              G = G + OOOO(LTO(j,k,m,n))
            ENDIF
            IF(I.EQ.M) THEN
              G = G - OOOO(LTO(j,k,l,n))
            ENDIF
            IF(I.EQ.N) THEN
              G = G + OOOO(LTO(j,k,l,m))
            ENDIF

            IF(J.EQ.L) THEN
              G = G - OOOO(LTO(i,k,m,n))
            ENDIF
            IF(J.EQ.M) THEN
              G = G + OOOO(LTO(i,k,l,n))
            ENDIF
            IF(J.EQ.N) THEN
              G = G - OOOO(LTO(i,k,l,m))
            ENDIF

            IF(K.EQ.L) THEN
              G = G + OOOO(LTO(i,j,m,n))
            ENDIF
            IF(K.EQ.M) THEN
              G = G - OOOO(LTO(i,j,l,n))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + OOOO(LTO(i,j,l,m))
            ENDIF

          ENDIF
c
c **** D_occ,occ block, alpha
c
          IF(J.EQ.M) THEN
            IF(K.EQ.N) THEN
              G = G - VOVO(LTV(r,l,s,i))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, beta
c
          IF(J.EQ.L) THEN
            IF(K.EQ.M) THEN
              G = G - VOVO(LTV(r,n,s,i))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + VOVO(LTV(r,m,s,i))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, gamma
c
          IF(I.EQ.M) THEN
            IF(J.EQ.N) THEN
              G = G - VOVO(LTV(r,l,s,k))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + VOVO(LTV(r,l,s,j))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, delta
c
          IF(I.EQ.L) THEN
            IF(J.EQ.M) THEN
              G = G - VOVO(LTV(r,n,s,k))
            ENDIF
            IF(J.EQ.N) THEN
              G = G + VOVO(LTV(r,m,s,k))
            ENDIF
            IF(K.EQ.M) THEN
              G = G + VOVO(LTV(r,n,s,j))
            ENDIF
            IF(K.EQ.N) THEN
              G = G - VOVO(LTV(r,m,s,j))
            ENDIF
          ENDIF
c
c according to the formulation on paper the c.c. satellite
c block entry is created here. We therefore have to c.c. the 
c entry.
c
          G=DCONJG(G)
c
c check if diagonal entry has to be added
c
          IF(SATROW.EQ.SATCOL) THEN
            G = G + DCMPLX(SATEV(SATROW),0.0D0)
          ENDIF
c
c fill this element in the output buffer and flush, if necessary.
c the buffer counter is set outside the row/column loop
c reminder: satrow and satcol count from 1 to length of sat block,
c but the matrix starts with 2h2h block ==> add offset to these
c variables.
c
          IF(ABS(G).GT.ADCTHR) THEN
            NTOTELE = NTOTELE + 1
            IB=IB+1
c
c   !!!
c
            BUF(IB)=G
            IOI(IB)=SATROW + SATOFF
            IOJ(IB)=SATCOL + SATOFF
          
            IF(IB.EQ.INTBUF) THEN
              WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &                     (IOI(IXX),IXX=1,INTBUF),
     &                     (IOJ(IXX),IXX=1,INTBUF),
     &                     INTBUF,JDUMMY
              IB = 0
              NBUFS = NBUFS + 1
            ENDIF
          ENDIF


 400    CONTINUE
 401  CONTINUE
c
c  write out *remaining* elements of the file buffer
c  and count all nonzero elements in the ADC matrix.
c
      IF(IB.GT.0) THEN
        NBUFS = NBUFS + 1
        WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &               (IOI(IXX),IXX=1,INTBUF),
     &               (IOJ(IXX),IXX=1,INTBUF),
     &               IB,JDUMMY
      ENDIF
      CLOSE(ITAP4)
      WRITE(IW,*)
      WRITE(IW,*) 'Total of nonzero elements in ADC matrix: ',NTOTELE
      WRITE(IW,*) 'Number of buffers in ADC matrix file   : ',NBUFS
      WRITE(IW,*) 'Dimension of ADC matrix in this symmetry:',
     6            IRCDLEN+SATOFF
      XNZEL=DBLE(NTOTELE)
      XTOTEL=DBLE(IRCDLEN+SATOFF)
      XTOTEL=XTOTEL*XTOTEL
      SPARS=(XNZEL/XTOTEL)*100.0D0
      SPARS=100.0D0 - SPARS
      WRITE(IW,'(A,F5.2,A)') ' Sparsity of ADC matrix:',SPARS,' %'
      WRITE(IW,*)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKE_EIJKL (EPS,E)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Create energy array 0.5*(ei + ej + ek + el) in the same order as a
C     <OO||OO> integral stream is stored !
C     this energy array is always real, no complexation considerd.
C     ATT: one factor of 0.5 is absorbed due to summation r,s -->  r>s !
C
C    **** important remark ***
C
C  is is very error-prone to change all the indices from IJ to
C  KL and from AB to IJ because this is the denomvvoo source!
C  We therefore are lazy but safe and *interpret* just for this case
C  the AB variables as IJ and the IJ variable as KL. Therefore AB are
C  now *occupieds*  !!! We just adjust the energy offsets and signs (all
C  energies are *positive*), everything else stays identical
C  due to the twofold triangular storage.
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 E(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 FAC,FAC1,FAC2,FAC3
C
C---------------Executable code--------------------------------------
C
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = + EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = + EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IO(BRP)
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            AOFF = IO(ARP)
            DO B = 1, NO(BRP)
               BB = BB + 1
               FAC3 = FAC2 + EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NO(ARP)
                  AA = AOFF + A
                  FAC = FAC3 + EPS(AA)
                  ABIJ = ABIJ + 1

                  E(ABIJ) = 0.5D0 * FAC

               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      IF(ABIJ.NE.IOOOOTT(NREP+1)) CALL QUIT('Error in MAKE_EIJKL')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKELOOKUPO(LT,LENLT)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C    **** important remark ***
C
C  is is very error-prone to change all the indices from IJ to
C  KL and from AB to IJ because this is the denomvvoo source!
C  We therefore are lazy but safe and *interpret* just for this case
C  the AB variables as IJ and the IJ variable as KL. Therefore AB are
C  now *occupieds*  !!! 
C
C---------------Calling variables--------------------------------------
C
      INTEGER LT(LENLT,LENLT,LENLT,LENLT)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER AABSL,BABSL,IABSL,JABSL
C
C---------------Executable code--------------------------------------
C
c
c  initialize lookup table (extremely important due to null offsets)
c  loop structure follows storage mode of the OOOO integrals.
c
      do i=1,lenlt
      do j=1,lenlt
      do k=1,lenlt
      do l=1,lenlt
        LT(l,k,j,i) = 0
      enddo
      enddo
      enddo
      enddo

c     WRITE(IW,*) 'Preparing OOOO Lookup table...'

      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            DO 20 BRP = 1, NREP
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            DO B = 1, NO(BRP)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NO(ARP)
                  ABIJ = ABIJ + 1

                  AABSL = IO(ARP) + A
                  BABSL = IO(BRP) + B
                  IABSL = IO(IRP) + I
                  JABSL = IO(JRP) + J

                  LT(AABSL,BABSL,IABSL,JABSL)=ABIJ

               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      IF(ABIJ.NE.IOOOOTT(NREP+1)) CALL QUIT('Error in MAKELOOKUPO')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKELOOKUPV(LT,LENLTO,LENLTV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C    **** important remark ***
C
C  is is very error-prone to change all the indices from IJ to
C  KL and from AB to IJ because this is the denomvvoo source!
C  We therefore are lazy but safe and *interpret* just for this case
C  the AB variables as IJ and the IJ variable as KL. Therefore AB are
C  now *occupieds*  !!! 
C
C---------------Calling variables--------------------------------------
C
      INTEGER LT(LENLTV,LENLTO,LENLTV,LENLTO)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER AABSL,BABSL,IABSL,JABSL
C
C---------------Executable code--------------------------------------
C
c
c  initialize lookup table (extremely important due to null offsets)
c  loop structure follows storage mode of the OOOO integrals.
c
      do i=1,lenlto
      do j=1,lenltv
      do k=1,lenlto
      do l=1,lenltv
        LT(l,k,j,i) = 0
      enddo
      enddo
      enddo
      enddo

c     WRITE(IW,*) 'Preparing VOVO Lookup table...'

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
                  AIBJ = AIBJ + 1

                  AABSL = IV(ARP) + A
                  IABSL = IO(IRP) + I
                  BABSL = IV(BRP) + B
                  JABSL = IO(JRP) + J

                  LT(AABSL,IABSL,BABSL,JABSL)=AIBJ

               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      IF(AIBJ.NE.IVOVO(NREP+1)) CALL QUIT('Error in MAKELOOKUPV')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKESPTRA(IROOO,ILENO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Creates the auxiliary array for conversion of relative/absolute
C     spinor numbers for occupied spinors (only these are needed)
C
C---------------Calling variables--------------------------------------
C
      INTEGER IROOO(ILENO,2)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/files.inc"
C
C---------------Executable code--------------------------------------
c
      IX=1
      Do KRP=1,NREP
        Do K=1,NO(KRP)
          IROOO(IX,1)=K
          IROOO(IX,2)=KRP
          IX = IX + 1
        Enddo
      Enddo
      If((IX-1).ne.ileno) STOP 'Consistency violation in MAKESPTRA!'

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




