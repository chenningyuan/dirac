!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_TRAMO(desrep,ioch,eps)
 
      IMPLICIT INTEGER (A-Z)
!     IMPLICIT none
!
!---------------Description--------------------------------------------
!
!     Calculates transition moments from the transformed
!     x,y,z dipole moment integrals and the Davidson eigenvectors
!
!     Important: Also for real two-electron integrals (high symmetry and/or 
!     NONREL,SPINFREE Hamiltonians) the property integrals are in complex form.
!     These cases will be handled by the corresponding phase factors stored at the
!     beginning of the property integral file.
!
!---------------Calling variables--------------------------------------
!
      INTEGER                          :: desrep,ioch
      REAL*8,dimension(:)              :: eps
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/complex.inc"
#include "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!

! filename of property integrals

      Character*8,dimension(3)               :: trafile

! constants

      Complex*16                             :: A0,A1,AI

! index arrays for the nontotally symmetric (NTS) property arrays

      INTEGER,allocatable,dimension(:)       :: JPOOXX,JPOVXX
      INTEGER,allocatable,dimension(:)       :: JPVOXX,JPVVXX

! offset arrays for NTS two-electron contributions
                                            
      INTEGER,allocatable,dimension(:)       :: ABIJTT
      INTEGER,allocatable,dimension(:,:)     :: JAVVOO
      INTEGER,allocatable,dimension(:,:)     :: JBVVOO
      INTEGER,allocatable,dimension(:)       :: MVVOOT
      INTEGER,allocatable,dimension(:)       :: BVVOOT

! phase information for the property x,y,z arrays in case of real arithmetic

      Complex*16,allocatable,dimension(:)    :: pphase
      Complex*16                             :: pphasex
      Complex*16                             :: pphasey
      Complex*16                             :: pphasez

! intermediate and final F_ai, F_aibj matrices (as vectors)

      REAL*8,allocatable,dimension(:)        :: faipx,faipy,faipz
      Complex*16,allocatable,dimension(:)    :: faitx,faity,faitz
      Complex*16,allocatable,dimension(:)    :: fabijtx,fabijty,fabijtz

! dipole moment integrals in MO basis for all three components

      REAL*8,allocatable,dimension(:)        :: doox,dooy,dooz
      REAL*8,allocatable,dimension(:)        :: dovx,dovy,dovz
      REAL*8,allocatable,dimension(:)        :: dvox,dvoy,dvoz
      REAL*8,allocatable,dimension(:)        :: dvvx,dvvy,dvvz
                                            
! auxiliary arrays for contractions

      REAL*8,allocatable,dimension(:)        :: aac,aki,aak,abj
                                            
! integral storage arrays

      REAL*8,allocatable,dimension(:)        :: vvoo1
      REAL*8,allocatable,dimension(:)        :: vvoo2
      REAL*8,allocatable,dimension(:)        :: vvoo3
      REAL*8,allocatable,dimension(:)        :: vvoo4
      REAL*8,allocatable,dimension(:)        :: vvoop
                                            
      REAL*8,allocatable,dimension(:)        :: vooo1,vooo2
      REAL*8,allocatable,dimension(:)        :: vovo1,vovo2,vovo3
      REAL*8,allocatable,dimension(:)        :: vovv1,vovv2
                                            
      REAL*8,allocatable,dimension(:)        :: oooo

      REAL*8,allocatable,dimension(:)        :: vabijx
      REAL*8,allocatable,dimension(:)        :: vabijy
      REAL*8,allocatable,dimension(:)        :: vabijz

      Integer                                :: n1,n2,n3,n4
      Integer                                :: ialloc,iat

! some parameters

      Logical                                :: LFA=.false.
      Logical                                :: LTR=.true.

      integer*8                              :: n8
      character*9                            :: preconfn='VVVVXVVOO'
      character*30                           :: iorder
      integer                                :: incsat
      logical                                :: okflag
!
!---------------interface region --------------------------------------
!
      interface

        SUBROUTINE propfetch(i1,ca1,
     &                          i2,i3,i4,i5,
     &                          ra1,ra2,ra3,
     &                          ra4,ra5,ra6,
     &                          ra7,ra8,ra9,
     &                          ra10,ra11,ra12)
          integer                      :: i1
          complex*16, dimension(:)     :: ca1
          integer                      :: i2,i3,i4,i5
          real*8, dimension(:)         :: ra1,ra2,ra3
          real*8, dimension(:)         :: ra4,ra5,ra6
          real*8, dimension(:)         :: ra7,ra8,ra9
          real*8, dimension(:)         :: ra10,ra11,ra12
        END SUBROUTINE

        INTEGER FUNCTION XCOLDET(i1)
          integer                      :: i1
        END FUNCTION

        SUBROUTINE denomvots(ra1,ra2)
          real*8, dimension(:)         :: ra1,ra2
        END SUBROUTINE

        SUBROUTINE denomvovo(ra1,ra2)
          real*8, dimension(:)         :: ra1,ra2
        END SUBROUTINE

        SUBROUTINE rcaxpy(i1,l1,c1,ra1,ca1)
          integer                      :: i1
          logical                      :: l1
          complex*16                   :: c1
          real*8, dimension(:)         :: ra1
          complex*16, dimension(:)     :: ca1
        END SUBROUTINE

        SUBROUTINE rca_negate(i1,ra1)
          integer                      :: i1
          real*8, dimension(:)         :: ra1
        END SUBROUTINE

        SUBROUTINE TM_GEN(i1,i2,i3,ch1,
     &                    l1,
     &                    ca1,ca2,ca3,
     &                    ca4,ca5,ca6)
          integer                       :: i1,i2,i3
          character*30                  :: ch1
          logical                       :: l1
          complex*16, dimension(:)      :: ca1,ca2,ca3
          complex*16, dimension(:)      :: ca4,ca5,ca6
        END SUBROUTINE

      end interface
!
!---------------Executable code--------------------------------------
!
      call pst('Calculating transition moments+')
      write(*,*)
      write(*,*) '   Note: Even in SPINFREE/NONREL cases the TMs'
      write(*,*) '   can be complex due to the quaternionic structure'
      write(*,*) '   of the wave function.'
      write(*,*)
      write(*,*) '   **********************'
      write(*,*) '   Symmetry: ',desrep
      write(*,*) '   **********************'
      write(*,*)
!
! initialization of constant variables
!
      trafile = (/'XPROPINT','YPROPINT','ZPROPINT'/)
      do i=1,30
        iorder(i:i)=' '
      enddo
      A0=(0.0d0,0.0d0)
      A1=(1.0d0,0.0d0)
      AI=(0.0d0,1.0d0)
!
! Generate index arrays for OO,OV,VO,VV property arrays
!
      n1 = nrep + 1
      allocate(JPOOXX(n1))
      allocate(JPVOXX(n1))
      allocate(JPOVXX(n1))
      allocate(JPVVXX(n1))

      JPOOXX=0
      JPVOXX=0
      JPOVXX=0
      JPVVXX=0

      DO IREP = 1, NREP
        JPOOXX(IREP+1) = JPOOXX(IREP) + MOO(IREP)
        JPOVXX(IREP+1) = JPOVXX(IREP) + MOV(IREP)
        JPVOXX(IREP+1) = JPVOXX(IREP) + MVO(IREP)
        JPVVXX(IREP+1) = JPVVXX(IREP) + MVV(IREP)
      ENDDO
!
! Generate index arrays for NTS two-particle contribs
! and check for correct coupling block lengths
!
      allocate(ABIJTT(NREP))

      allocate(JAVVOO(NREP,NREP))
      allocate(JBVVOO(NREP,NREP))
      allocate(MVVOOT(NREP))
      allocate(BVVOOT(NREP))

      ABIJTT  = 0
      JAVVOO  = 0
      JBVVOO  = 0
      MVVOOT  = 0
      BVVOOT  = 0

      DO JREP = 1, NREP
        DO IREP = 1, NREP
           IJREP = MULTB(IREP,JREP,2)
           IJBREP = MULTB(IREP+NREP,JREP+NREP,2)  !ABIJTT refers to bosonic ireps!

           ABIJTT(IJBREP) = ABIJTT(IJBREP) + NVVT(IREP)*NOOT(JREP)
           JAVVOO(IREP,JREP) = MVVOOT(IJREP)
           JBVVOO(IREP,JREP) = BVVOOT(IJREP)
           MVVOOT(IJREP) = MVVOOT(IJREP) + NVVOT(IREP)*NO(JREP)
           BVVOOT(IJREP) = BVVOOT(IJREP) + NV(IREP)*NVOOT(JREP)

        ENDDO
      ENDDO

      DO irep = 1, NREP
        if(ABIJTT(irep).ne.xcoldet(irep)) then
          write(*,*) 'Internal error in coupling block lengths!'
          write(*,*) 'IREP:',irep,' ABIJTT(irep):',ABIJTT(irep),
     &               ' XCOLDET(irep):',xcoldet(irep)
          Call quit('**********')
        endif
      ENDDO

!*****************************************
! complex array for the phases.
!*****************************************

      allocate(pphase(3))
      pphase = A0
!
!***************************************************
! ****  Transition element vectors. Attention
! ****  Even in a real calculation they will be kept
! ****  complex!
! ****  the final transition moment arrays are always complex.
! ****  the intermediate transition moment arrays follow CARITH.
!***************************************************
!
      n1 = mvo(desrep)

      allocate(faitx(n1))
      allocate(faity(n1))
      allocate(faitz(n1))
      faitx = A0
      faity = A0
      faitz = A0
      
      allocate(faipx(rcw*n1))
      allocate(faipy(rcw*n1))
      allocate(faipz(rcw*n1))
      faipx = 0.0d0
      faipy = 0.0d0
      faipz = 0.0d0
!
!***************************************************
! ****  Arrays holding the property integrals. 
! ****  They are either
! ****  complex if CARITH=T or real, if CARITH=F.
! ****  However, in the latter case a phase information
! ****  indicating purely imaginary property integrals
! ****  has to be additionally stored and processed.
!***************************************************
!
      idaoo = JPOOXX(NREP+1)   !  dimensions
      idaov = JPOVXX(NREP+1)
      idavo = JPVOXX(NREP+1)
      idavv = JPVVXX(NREP+1)
      n1 = idaoo*rcw   ! lengths
      n2 = idaov*rcw
      n3 = idavo*rcw
      n4 = idavv*rcw
      allocate(doox(n1),dooy(n1),dooz(n1))
      allocate(dovx(n2),dovy(n2),dovz(n2))
      allocate(dvox(n3),dvoy(n3),dvoz(n3))
      allocate(dvvx(n4),dvvy(n4),dvvz(n4))
      doox = 0.0d0
      dooy = 0.0d0
      dooz = 0.0d0
      dovx = 0.0d0
      dovy = 0.0d0
      dovz = 0.0d0
      dvox = 0.0d0
      dvoy = 0.0d0
      dvoz = 0.0d0
      dvvx = 0.0d0
      dvvy = 0.0d0
      dvvz = 0.0d0

!*********************************************************************
!  VVOO-type integral arrays, they always remain allocated since they
!  hold integrals repeatedly needed in nearly every step.
!*********************************************************************

      n1 = ivvoott(nrep+1)
      allocate(vvoo1(n1*rcw))
      allocate(vvoo2(n1*rcw))
      vvoo1 = 0.0d0
      vvoo2 = 0.0d0

      n1 = max0(lvvoot(nrep+1),jvovo(nrep+1),
     &          ivvoott(nrep+1),ivvoot(nrep+1),
     &          JVOOV(nrep+1))

      allocate(vvoo3(n1*rcw))
      allocate(vvoo4(n1*rcw))
      vvoo3 = 0.0d0
      vvoo4 = 0.0d0
!
!  auxiliary arrays
!
      n1 = MVV(1)
      allocate(aac(n1*rcw))
      n1 = MOO(1)
      allocate(aki(n1*rcw))
      n1 = MVO(1)
      allocate(aak(n1*rcw))
      n1 = JPVOXX(NREP+1) 
      allocate(abj(n1*rcw))
      aac = 0.0d0
      aki = 0.0d0
      aak = 0.0d0
      abj = 0.0d0
!______________________________________________________
!|
!|
!|
!|
!|
!|       START CALCULATIONS
!|
!|
!|
!|
!|_____________________________________________________
!
!  Initialize VVOO integral arrays
!  v_abij ---> vvoo1  /  v*_cbjk ---> vvoo2
!
      n1 = ivvoott(nrep+1)
      call getvvoo(vvoo1)
      call xcopy(n1,vvoo1,1,vvoo2,1)
      if(carith) call conjuga(n1,vvoo2,1)
      call denomvvoo(eps,vvoo1,vvoo1)
      call denomvvoo(eps,vvoo2,vvoo2)
!
!  get the property integrals and phases from corresponding files.
!
      CALL propfetch(ioch,pphase,    !fetch properties, we transfer dimensions!
     &               idaoo,idaov,idavo,idavv,
     &               doox,dooy,dooz,
     &               dovx,dovy,dovz,
     &               dvox,dvoy,dvoz,
     &               dvvx,dvvy,dvvz)

      pphasex = pphase(1)
      pphasey = pphase(2)
      pphasez = pphase(3)
!
!-----------------------------------------------
!------ Determine relevant irep such that
!------ b * b-1 = A where b=desrep of the state
!------ and check consistency.
!-----------------------------------------------
!
      proprep = multb(desrep+nrep,1+nrep,2)
      if(MVO(desrep).ne.MVO(proprep))then
        CALL QUIT('Property dimension inconsistency!')
      endif
!
!------------------------
!------  F_ai(0)
!------------------------
!
      n1 = MVO(proprep)       !number of integrals
      off = JPVOXX(proprep) + 1       !offset from where to read the integrals

      call xcopy(n1,dvox(off),1,faipx,1)
      call xcopy(n1,dvoy(off),1,faipy,1)
      call xcopy(n1,dvoz(off),1,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th order                     '
      incsat = 0
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
!
!------------------------
!------  F_ai(1)
!------------------------
!
! sort v_ab,ij --> v_ai,bj
! the minus sign in the contraction is accounted for in the 
! XGEMV routines!
!
      vvoo3 = 0.0d0
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &                 MVO,JVOVO,JJVO,JJVO,vvoo1,vvoo3)
      off = JPVOXX(desrep) + 1
      n1 = MVO(desrep)
      if(carith) then
        call conjuga(n1,dvox(off),1)
        call conjuga(n1,dvoy(off),1)
        call conjuga(n1,dvoz(off),1)
      endif

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOVO(desrep) + 1
      off2 = JPVOXX(desrep) + 1
      M=MVO(desrep)
      K=MVO(desrep)
      CALL XGEMV ('N',M,K,-A1,vvoo3(off1),M,
     &            dvox(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,-A1,vvoo3(off1),M,
     &            dvoy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,-A1,vvoo3(off1),M,
     &            dvoz(off2),1,A0,faipz,1)
      if(carith) then
        call conjuga(n1,dvox(off),1)
        call conjuga(n1,dvoy(off),1)
        call conjuga(n1,dvoz(off),1)
      endif
!
! add result on top of the final arrays
!
      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)
!
!  calculate corresponding 0th + 1st order TMs
!
      iorder='0th + 1st order               '
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

!
!  enter calculation for second-order contributions
!
!------------------------
!------  F_ai(2)A
!------------------------
!
!  sort  v_ab,jk  -->  v_a,bjk    (vvoo3)
!  sort v*_cb,jk  -->  v*_c,bjk   (vvoo4)
!         v_a,bjk  *  v*_c,bjk   -->  aux1
!
      call srt1t2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      call srt1t2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo2,vvoo4)
      call cntrct ('N','T',nv,nv,nvoot,a1,vvoo3,vvoo4,A0,aac,nrep)
!
!  contract A_ac and D_ci(x,y,z) over c in special loop because
!  D_ci is not totally symmetric and the order of Y_c in D_ci does not
!  follow the one in A_ac.
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         crep = multb(desrep+nrep,irep,2)
         arep = crep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NV(crep)
         ioff1 = JJVV(arep,crep) + 1
         ioff2 = JPVOXX(desrep) + JJVO(crep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,aac(ioff1),m,dvox(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aac(ioff1),m,dvoy(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aac(ioff1),m,dvoz(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo
      n1 = MVO(desrep)


      CALL XSCAL(n1,-A1/2.0d0,faipx,1)
      CALL XSCAL(n1,-A1/2.0d0,faipy,1)
      CALL XSCAL(n1,-A1/2.0d0,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2)A finished'

      iorder='0th + 1st + 2ndA'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
!
!------------------------
!------  F_ai(2)B
!------------------------
!
!  sort  v*_bc,jk (vvoo2)  -->  v*_bcj,k  (vvoo3)  (1T3)
!  sort  v_bc,ji  (vvoo1)  -->  v_bcj,i   (vvoo4)  (1T3)
!
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             vvoo2,vvoo3)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             vvoo1,vvoo4)
      CALL CNTRCT('T','N',NO,NO,NVVOT,A1,vvoo3,vvoo4,A0,aki,NREP)
!
!  contract D_ak(x,y,z) and A_ki over k in special loop because
!  D_ak is not totally symmetric!
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         arep = multb(desrep+nrep,irep,2)
         krep = irep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NO(krep)
         ioff1 = JPVOXX(desrep) + JJVO(arep,krep) + 1
         ioff2 = JJOO(krep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,dvox(ioff1),m,aki(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvoy(ioff1),m,aki(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvoz(ioff1),m,aki(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      CALL XSCAL(n1,-A1/2.0d0,faipx,1)
      CALL XSCAL(n1,-A1/2.0d0,faipy,1)
      CALL XSCAL(n1,-A1/2.0d0,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2)B finished'

      iorder='0th + 1st + 2ndAB'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
!
!------------------------
!------  F_ai(2)C
!------------------------
!
!
!  sort  v_ab,ij  (vvoo1)  -->  v_ai,bj    (vvoo3)  (1TT4)
!  v*_cb,kj = v*_bc,jk !
!  sort  v*_bc,jk (vvoo2)  -->  v*_bj,ck   (vvoo4)  (1TT4)

!  additionally, allocate proper array
!
      vvoo3 = 0.0d0
      vvoo4 = 0.0d0
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &                 MVO,JVOVO,JJVO,JJVO,vvoo1,vvoo3)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &                 MVO,JVOVO,JJVO,JJVO,vvoo2,vvoo4)
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      vovo1 = 0.0d0
      if(ialloc.ne.0)  call quit('vovo1 allocation problem (tramo)')
!
!  do contraction of v_ai,bj  *  v*_bj,ck  =  A_ai,ck
!
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vvoo3(OFF1),M,vvoo4(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
!  do contraction of A_ai,ck  *  D_ck (determine proper symmetry!)
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      M = MVO(desrep)
      K = MVO(desrep)
      OFF1 = jvovo(desrep) + 1
      OFF2 = JPVOXX(desrep) + 1
      CALL XGEMV ('N',M,K,A1,vovo1(OFF1),M,dvox(OFF2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vovo1(OFF1),M,dvoy(OFF2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vovo1(OFF1),M,dvoz(OFF2),1,A0,faipz,1)
      deallocate(vovo1)
!
!  scale result by +0.5d0 and add to final values
!
      n1 = MVO(desrep)
      CALL XSCAL(n1,A1/2.0d0,faipx,1)
      CALL XSCAL(n1,A1/2.0d0,faipy,1)
      CALL XSCAL(n1,A1/2.0d0,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)
      
      write(*,*) '------  F_ai(2)C finished'

      iorder='0th + 1st + 2ndABC'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
!
!
!------------------------
!------  F_ai(2,1)
!------------------------
!
!  only temporarily allocate space for the VOVV integrals.
!
      iat = 0
      n1 = ivovvt(nrep+1)
      allocate(vovv1(n1*rcw),stat=ialloc); iat = iat + ialloc
      n1 = lvovvt(nrep+1)
      allocate(vovv2(n1*rcw),stat=ialloc); iat = iat + ialloc
      if(ialloc.ne.0) STOP 'VOVV allocation problem in TRAMO!'
!
!  sort V_aj,bc  --->  Va,jbc   (1S2)    vovv1 --> vovv2
!  sort v_bc,kj  --->  v_jbc,k  (SRT22)  vvoo1 --> vvoo3
!  then contract Va,jbc  *  v_jbc,k  ---> A_ak (tot symm)
!  divide A_ak by (eps_a  -  eps_k), and multiply by -1.0
!  and contract A_ak with D_ki(x,y,z) (special loop!)
!
      CALL RDVOVV(vovv1); vovv2 = 0.0d0
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NV,NO,NVVT,LVOVVT,LLOVVT,
     &             vovv1,vovv2)
      CALL SRT22 (NREP,MULTB,LFA,NVVT,NO,NO,NOVVT,LOOVVT,LLOVVT,
     &            vvoo1,vvoo3)
      CALL CNTRCT ('N','N',NV,NO,NOVVT,A1,vovv2,vvoo3,A0,aak,NREP)
      CALL DENOMVOTS (EPS,aak)
      CALL XSCAL(MVO(1),-A1,aak,1)
      deallocate(vovv1,vovv2)
!
!  contract A_ak/(e_a - e_k) and D_ki(x,y,z) over k in special loop because
!  D_ki is not totally symmetric and the order of Y_k in D_ki does not
!  follow the one in A_ak.
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         krep = multb(desrep+nrep,irep,2)
         arep = krep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NO(krep)
         ioff1 = JJVO(arep,krep) + 1
         ioff2 = JPOOXX(desrep) + JJOO(krep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,doox(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooy(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooz(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,1'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

      write(*,*) '------  F_ai(2,1) finished'
      
!
!------------------------
!------  F_ai(2,2)
!------------------------
!
      iat = 0
      n1 = ivooot(nrep+1)
      allocate(vooo1(n1*rcw),stat=ialloc); iat = iat + ialloc
      n1 = lovoot(nrep+1)
      allocate(vooo2(n1*rcw),stat=ialloc); iat = iat + ialloc
      if(ialloc.ne.0) STOP 'VOOO allocation problem in TRAMO!'
!
!  sort v_ab,jk  --->  v_a,bjk   (still in vvoo1 ---> vvoo3)
!  sort V*_bl,jk --->  V*_bjk,l  vooo1 ---> vooo2
!  contract  v_a,bjk  *  V*_bjk,l  ---> A_al (aak)
!  divide A_ak by (eps_a  -  eps_k), and multiply by -1.0
!  and contract A_ak with D_ki(x,y,z)
!
      call SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      CALL GETVOOO(vooo1);vooo2 = 0.0d0
      IF(CARITH) CALL CONJUGA(ivooot(nrep+1),vooo1,1)
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           vooo1,vooo2)
      CALL CNTRCT ('N','N',NV,NO,NVOOT,A1,vvoo3,vooo2,A0,aak,NREP)
      CALL DENOMVOTS (EPS,aak)
      CALL XSCAL(MVO(1),-A1,aak,1)
      deallocate(vooo1,vooo2)
!
!  special loop as above for contracting A_al and D_li (identical to
!  A_ak and D_ki)
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      do irep = 1, nrep
         krep = multb(desrep+nrep,irep,2)
         arep = krep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NO(krep)
         ioff1 = JJVO(arep,krep) + 1
         ioff2 = JPOOXX(desrep) + JJOO(krep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,doox(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooy(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooz(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,2'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      
      write(*,*) '------  F_ai(2,2) finished'
!
!-------------------------------------
!------  F_ai(2,3), a DVV contribution
!-------------------------------------
      n1 = ivooot(nrep+1)
      allocate(vooo1(n1*rcw),stat=ialloc)
      n1 = lovoot(nrep+1)
      allocate(vooo2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOOO allocation problem in TRAMO!'
!
!  sort v_bc,jk  --->  v_b,cjk   (still in vvoo1 ---> vvoo3)
!  sort V*_ci,jk --->  V*_cjk,i  vooo1 ---> vooo2
!  contract  v_b,cjk  *  V*_cjk,i  ---> A_bi (aak)
!  divide A_bi by (eps_b  -  eps_i)
!  and contract D_ab(x,y,z) with A_bi
!
      call SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      CALL GETVOOO(vooo1)
      IF(CARITH) CALL CONJUGA(ivooot(nrep+1),vooo1,1)
      vooo2 = 0.0d0
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           vooo1,vooo2)
      CALL CNTRCT ('N','N',NV,NO,NVOOT,A1,vvoo3,vooo2,A0,aak,NREP)
      deallocate(vooo1,vooo2)
      CALL DENOMVOTS (EPS,aak) !A_bi now in aak.
!
! do special contraction of D_ab and A_bi over b
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         arep = multb(desrep+nrep,irep,2)
         brep = irep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NV(brep)
         ioff1 = JPVVXX(desrep) + JJVV(arep,brep) + 1
         ioff2 = JJVO(brep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,dvvx(ioff1),m,aak(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvy(ioff1),m,aak(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvz(ioff1),m,aak(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,3'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

      write(*,*) '------  F_ai(2,3) finished'
!
!-------------------------------------
!------  F_ai(2,4), a DVV contribution
!-------------------------------------
!
!  only temporarily allocate space for the VOVV integrals.
!
      n1 = ivovvt(nrep+1)
      allocate(vovv1(n1*rcw),stat=ialloc)
      n1 = lvovvt(nrep+1)
      allocate(vovv2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOVV allocation problem in TRAMO!'
!
!  sort v_bc,ij  --->  v_jbc,i  (SRT22)  vvoo1 --> vvoo3
!  sort V_dj,bc  --->  Vd,jbc   (1S2)    vovv1 --> vovv2
!  then contract Vd,jbc  *  v_jbc,i  ---> A_di (tot symm)
!  divide A_di by (eps_d  -  eps_i)
!  and contract D_ad(x,y,z) with Adi
!
      CALL SRT22 (NREP,MULTB,LFA,NVVT,NO,NO,NOVVT,LOOVVT,LLOVVT,
     &            vvoo1,vvoo3)
      CALL RDVOVV(vovv1); vovv2 = 0.0d0
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NV,NO,NVVT,LVOVVT,LLOVVT,
     &             vovv1,vovv2)
      CALL CNTRCT ('N','N',NV,NO,NOVVT,A1,vovv2,vvoo3,A0,aak,NREP)
      deallocate(vovv1, vovv2)
      CALL DENOMVOTS (EPS,aak)
!
! do special contraction of D_ad and A_di over d! (loop from above)
! is the same as of D_ab and A_bi over b
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         arep = multb(desrep+nrep,irep,2)
         brep = irep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NV(brep)
         ioff1 = JPVVXX(desrep) + JJVV(arep,brep) + 1
         ioff2 = JJVO(brep,irep) + 1
         ioff3 = JJVO(arep,irep) + 1
         call xgemm ('N','N',m,n,k,a1,dvvx(ioff1),m,aak(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvy(ioff1),m,aak(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvz(ioff1),m,aak(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,4'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

      write(*,*) '------  F_ai(2,4) finished'
!
!-------------------------------------
!------  F_ai(2,5), a DOV contribution
!-------------------------------------
!
!  read Vaj,bi ---> Vai,bj (SRT16) (vovo1 ---> vovo2)
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
        if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      allocate(vovo2(n1*rcw),stat=ialloc)
        if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      allocate(vovo3(n1*rcw),stat=ialloc)
        if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      CALL getvovo(vovo1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo1,vovo2)
!
!  v_b>c,k>j  --->  -v_b>c,j>k  (RCA_NEGATE) vvoo1 --> -vvoo1
!  -v_b>c,j>k  --->  -v_bj,ck   (SRT1TT4) -vvoo1 --> vovo1
!  negate vvoo1 again in order to restore original values.

      n1 = ivvoott(nrep+1) 
      CALL RCA_NEGATE(n1,vvoo1)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vovo1)
      CALL RCA_NEGATE(n1,vvoo1)

!  contract Vai,bj (vovo2) and -v_bj,ck (vovo1) to Aai,ck (vovo3)
!  Aai,ck  attains  <||>,<||>  BKC  structure

      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = JVOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vovo2(OFF1),M,vovo1(OFF2),K,
     &               A0,vovo3(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO

! ************** ATT !  no backsort, instead denomvovo!
! divide A_ai,ck array (JVOVO structure) by the corresponding
! energy denominators e_a - e_i + e_c - e_k

      CALL DENOMVOVO(eps,vovo3)

!  contract A_ai,kc (vovo3)  and Dkc over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOVO(desrep) + 1
      off2 = JPVOXX(desrep) + 1
      M=MVO(desrep)
      K=MVO(desrep)

      if(carith) call conjuga(M,dvox(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dvox(off2),1,A0,faipx,1)
      if(carith) call conjuga(M,dvox(off2),1)

      if(carith) call conjuga(M,dvoy(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dvoy(off2),1,A0,faipy,1)
      if(carith) call conjuga(M,dvoy(off2),1)

      if(carith) call conjuga(M,dvoz(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dvoz(off2),1,A0,faipz,1)
      if(carith) call conjuga(M,dvoz(off2),1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      deallocate(vovo1,vovo2,vovo3)

      iorder='0th + 1st + 2ndABC + 2,5'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

      write(*,*) '------  F_ai(2,5) finished'
!
!-------------------------------------
!------  F_ai(2,6): a DOV contribution
!-------------------------------------
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      allocate(vovo2(n1*rcw),stat=ialloc)
      allocate(vovo3(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
!
!  sort the -v_ab,ji (vvoo1) ---> -v_ai,bj (vovo3)
!  and the V*_bk,cj ---> V*_bj,ck (SRT 16)
!  then contract -v_ai,bj and V*_bj,ck to  A_ai,ck
!       VO,VO  x  VO,VO  --->  VO,VO <||>,<||>
!
      n1 = ivvoott(nrep+1) 
      CALL RCA_NEGATE(n1,vvoo1)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vovo3)
      CALL RCA_NEGATE(n1,vvoo1)

      CALL getvovo(vovo1)
      n1 = ivovo(nrep+1) 
      if(carith) call conjuga(n1,vovo1,1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo1,vovo2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vovo3(OFF1),M,vovo2(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
!  A_ai,ck now in vovo1
!  divide by energy denominators
!
      CALL DENOMVOVO(eps,vovo1)

!  contract B_ai,ck (vovo1) and Dkc = D*ck over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOVO(desrep) + 1
      off2 = JPVOXX(desrep)+1
      M=MVO(desrep)
      K=MVO(desrep)

      if(carith) call conjuga(M,dvox(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dvox(off2),1,A0,faipx,1)
      if(carith) call conjuga(M,dvox(off2),1)

      if(carith) call conjuga(M,dvoy(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dvoy(off2),1,A0,faipy,1)
      if(carith) call conjuga(M,dvoy(off2),1)

      if(carith) call conjuga(M,dvoz(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dvoz(off2),1,A0,faipz,1)
      if(carith) call conjuga(M,dvoz(off2),1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      deallocate(vovo1,vovo2,vovo3)

      iorder='0th + 1st + 2ndABC + 2,6'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

      write(*,*) '------  F_ai(2,6) finished'
!
!-------------------------------------
!------  F_ai(2,7): a DOV contribution
!-------------------------------------
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      allocate(vovo2(n1*rcw),stat=ialloc)
      allocate(vovo3(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
!
!  sort the v_cb,ji (vvoo1) ---> v_cj,bi (1TT4), vovo3
!  and the V_aj,ck ---> V_ak,cj (SRT 16)
!  contract V_ak,cj  x  v_cj,bi  to A_ak,bi
!
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vovo3)
      CALL getvovo(vovo1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo1,vovo2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vovo2(OFF1),M,vovo3(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO

! A_ak,bi is now in vovo1
! divide it by the energy denominators --> B_ak,bi

      Call DENOMVOVO(eps,vovo1)

! next we tentatively do a SRT16 inverse sort in order to
! arrive at the B_ai,bk array. (watch the array order in SRT16
! inverse!) The B_ai,bk array is now in vovo2

      CALL SRT16 (NREP,MULTB,LFA,LTR,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo2,vovo1)

!  contract B_ai,bk (vovo2)  and Dbk* = Dkb over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOVO(desrep) + 1
      off2 = JPVOXX(desrep)+1
      M=MVO(desrep)
      K=MVO(desrep)
!**************************************** UTH
      if(carith) call conjuga(M,dvox(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dvox(off2),1,A0,faipx,1)
      if(carith) call conjuga(M,dvox(off2),1)

      if(carith) call conjuga(M,dvoy(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dvoy(off2),1,A0,faipy,1)
      if(carith) call conjuga(M,dvoy(off2),1)

      if(carith) call conjuga(M,dvoz(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dvoz(off2),1,A0,faipz,1)
      if(carith) call conjuga(M,dvoz(off2),1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,7'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

      write(*,*) '------  F_ai(2,7) finished'
      call quit('Temporary stop')
!
!-------------------------------------
!------  F_ai(2,8): a DOV contribution
!-------------------------------------
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      allocate(vovo2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
!
!  sort the v_ac,kj (vvoo1) ---> v_ak,cj (1TT4)
!  and the V*_ci,bj ---> V_cj,bk (SRT 16)
!  contract v_ak,cj  and   V*_cj,bi to A_ak,bi
!
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vvoo3)
      CALL getvovo(vovo1)
      if(carith) CALL conjuga(ivovo(nrep+1),vovo1,1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo1,vovo2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vvoo3(OFF1),M,vovo2(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
! A_ak,bi is now in vovo1, and we can proceed exactly in the same way
! as for F_ai(2,7)!
!
! Backsort A_ak,bi ---> B_ab,ki (1TT4 inv)
! Divide by energy denominator and negate in order to
! obtain b_ab,ik (note lowercase letter!)
!
      vvoo3 = 0.0d0
      CALL SRT1TT4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vovo1,vvoo3)
      CALL DENOMVVOO (eps,vvoo3,vvoo3)
      CALL RCA_NEGATE(ivvoott(nrep+1),vvoo3)
      deallocate(vovo1,vovo2)
!
! resort b_a>b,i>k to b_ab,i>k to C_ai,kb
!
      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,
     &              vvoo3,vvoo4)
      CALL SRT20_NTS(NREP,MULTB,NVV,NV,NV,NO,NO,MVO,
     &               JVOOV,JJVO,JJOV,vvoo4,vvoo3)

!  contract C_ai,kb (vvoo3)  and Dkb over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOOV(desrep) + 1
      off2 = JPOVXX(desrep)+1
      M=MVO(desrep)
      K=MOV(desrep)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,8'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

      write(*,*) '------  F_ai(2,8) finished'
!
!-------------------------------------
!------  F_ai(2,9): a DOV contribution
!-------------------------------------
!
! contraction of V_ab,cd and v_cd,ij is done in precontraction section
! of main part. We read in the VVVV x VVOO contracted integrals.
!
!
! Read the A_a>b,i>j to vvoop
!
      n1 = ivvoott(nrep+1)
      allocate(vvoop(n1*rcw))
      open(ioch,file=preconfn,form='unformatted',
     &     access='sequential',status='unknown')
      read(ioch) (vvoop(K),K=1,n1*RCW)
      close(ioch)
      CALL DENOMVVOO (eps,vvoop,vvoop)
!
! resort b_a>b,i>j to b_ab,i>j to C_ai,jb,
! contract with D_jb and store.
!
      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,
     &              vvoop,vvoo4)
      CALL SRT20_NTS(NREP,MULTB,NVV,NV,NV,NO,NO,MVO,
     &               JVOOV,JJVO,JJOV,vvoo4,vvoo3)
      deallocate(vvoop)

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOOV(desrep) + 1
      off2 = JPOVXX(desrep)+1
      M=MVO(desrep)
      K=MOV(desrep)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,9'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

      write(*,*) '------  F_ai(2,9) finished'
!
!-------------------------------------
!-----  F_ai(2,10): a DOV contribution
!-------------------------------------
!
!  contract v_ab,kl and V_kl,ij to A_ab,ij
!  from here we go on as in F_ai(2,9)
!
      n1=ioooott(nrep+1)
      allocate(oooo(n1*rcw))
      call getoooo(oooo)
      CALL CNTRCT('N','N',NVVT,NOOT,NOOT,A1,vvoo1,oooo,A0,vvoo3,NREP)
      deallocate(oooo)
!
      CALL DENOMVVOO (eps,vvoo3,vvoo3)
      CALL SRT1LS1 (NREP,MULTB,NVVT,NVV,NV,NV,NOOT,IVVOOT,IIVV,
     &              vvoo3,vvoo4)
      CALL SRT20M(NREP,MULTB,NVV,NV,NV,NO,NO,MVO,
     &            JVOOV,JJVO,JJOV,vvoo4,vvoo3)

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOOV(desrep) + 1
      off2 = JPOVXX(desrep)+1
      M=MVO(desrep)
      K=MOV(desrep)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      write(*,*) '------  F_ai(2,10) finished'
!
! ----------------------------------------------------------
! --- from here all second order h/p contributions are -----
! ---      done and we print the resulting TMs       -------
! ----------------------------------------------------------
!
      iorder='0th-2nd order complete ph'
 7744 CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
!
!  ----------------------------------------------------------
!  -------     2h2p second order contributions to F  --------
!  ----------------------------------------------------------
!
!  allocate required space, faibj arrays have to be complex!
!
      write(*,*) '------  calculating the 2h2p 2nd order contributions'
      n1=ABIJTT(desrep)  !ABIJTT is non-cumulative!
      write(*,*) 'Allocating F_abij buffers of length',n1
      allocate(fabijtx(n1))
      allocate(fabijty(n1))
      allocate(fabijtz(n1))
      fabijtx = A0
      fabijty = A0
      fabijtz = A0
      allocate(vvoop(n1*rcw))
      vvoop = 0.0d0


      n1=mvvoot(desrep)  !mvvoot is non-cumulative!
      write(*,*) 'Allocating A_abi,k buffers of length',n1
      allocate(vabijx(n1*rcw))
      allocate(vabijy(n1*rcw))
      allocate(vabijz(n1*rcw))
      vabijx = 0.0d0
      vabijy = 0.0d0
      vabijz = 0.0d0
!
!---------------------------------------------------
!-----  F_ab,ij(2,A+B and antisymmetrization) ------
!---------------------------------------------------
!
!  v_ab,ki = -v_ab,ik (vvoo1)
!  sort -v_ab,ik to -v_abi,k  (vvoo3)
!  contract v_abi,k and D(x,y,z)kj  ---> Babi,j **NTS**
!  backsort and antisymmetrize Babi,j ---> Bab,ij.
!
      n1=ivvoott(nrep+1)
      CALL RCA_NEGATE(n1,vvoo1)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             vvoo1,vvoo3)
      CALL RCA_NEGATE(n1,vvoo1)
!
! apply proper symmetry loop due to nontotally symmetric D_oo contributions
! 
      do krep = 1, nrep
        jrep = multb(krep,desrep+nrep,1)
        if(MULTB(krep,jrep,2).ne.desrep) then
          stop 'munch:-('
        endif
        m = NVVOT(krep)
        n = NO(jrep)
        k = NO(krep)
        ioff1 = KVVOOT(krep) + 1
        ioff2 = JPOOXX(desrep) + JJOO(krep,jrep) + 1
        ioff3 = JAVVOO(krep,jrep) + 1
        write(*,*) 'm,n,k:',m,n,k
!x
        call xgemm ('N','N',m,n,k,a1,vvoo3(ioff1),m,doox(ioff2),k,
     &              a0,vabijx(ioff3),m)
!y
        call xgemm ('N','N',m,n,k,a1,vvoo3(ioff1),m,dooy(ioff2),k,
     &              a0,vabijy(ioff3),m)
!z
        call xgemm ('N','N',m,n,k,a1,vvoo3(ioff1),m,dooz(ioff2),k,
     &              a0,vabijz(ioff3),m)
      enddo

      n1 = ABIJTT(desrep)
      write(*,*) 'ICEIL 1:',n1
      vvoop = 0.0d0
      CALL SRT1T3_NTSI (NREP,desrep,MULTB,NVVT,NO,NO,NVVOT,
     &                  JAVVOO,KKVVOT,vabijx,vvoop,n1)
! debug test resorted x array:
      do i=1,n1
        write(*,*) i,vvoop(i)
      enddo
! debug end
      call rcaxpy(n1,carith,pphasex,vvoop,fabijtx)
      vvoop = 0.0d0
      CALL SRT1T3_NTSI (NREP,desrep,MULTB,NVVT,NO,NO,NVVOT,
     &                  JAVVOO,KKVVOT,vabijy,vvoop,n1)
! debug test resorted x array:
      do i=1,n1
        write(*,*) i,vvoop(i)
      enddo
! debug end
      call rcaxpy(n1,carith,pphasey,vvoop,fabijty)
      vvoop = 0.0d0
      CALL SRT1T3_NTSI (NREP,desrep,MULTB,NVVT,NO,NO,NVVOT,
     &                  JAVVOO,KKVVOT,vabijz,vvoop,n1)
! debug test resorted x array:
      do i=1,n1
        write(*,*) i,vvoop(i)
      enddo
! debug end
      call rcaxpy(n1,carith,pphasez,vvoop,fabijtz)

      deallocate(vvoop)
      deallocate(vabijx)
      deallocate(vabijy)
      deallocate(vabijz)

      write(*,*) 'F_ab,ij(2,A+B) contributions finished.'
      write(*,*) 'skipping the F_ab,ij(2,C+D) contributions.'
      goto 7766
!
!---------------------------------------------------
!-----  F_ab,ij(2,C+D and antisymmetrization)  -----
!---------------------------------------------------
!
      n1=BVVOOT(desrep)
      write(*,*) 'Allocating three B_a,bij buffers of length',n1
      allocate(vabijx(n1*rcw)) !still called vabij? for simplicity
      allocate(vabijy(n1*rcw))
      allocate(vabijz(n1*rcw))
      vabijx = 0.0d0
      vabijy = 0.0d0
      vabijz = 0.0d0

      n1=ABIJTT(desrep)
      write(*,*) 'Allocating one more F_aibj buffer of length',n1
      allocate(vvoop(n1*rcw))
      vvoop = 0.0d0
!
!  sort v_cb,ij to v_c,bij (1T2)  -->  vvoo3 (dim = lvvoot(nrep+1)! )
!  contract Dac with v_c,bij to B_a,bij, keeping proper symmetry
!  backsort B_a,bij to B_ab,ij  including antisymmetrization.
!
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      do crep = 1, nrep
         arep = multb(desrep+nrep,crep,2)
         if(MULTB(arep,crep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NVOOT(crep)
         k = NV(crep)
         ioff1 = JPVVXX(desrep) + JJVV(arep,crep) + 1
         ioff2 = LVVOOT(crep) + 1
         ioff3 = JBVVOO(arep,crep) + 1
!x
         call xgemm ('N','N',m,n,k,a1,dvvx(ioff1),m,vvoo3(ioff2),k,
     &               a0,vabijx(ioff3),m)
!y
         call xgemm ('N','N',m,n,k,a1,dvvy(ioff1),m,vvoo3(ioff2),k,
     &               a0,vabijy(ioff3),m)
!z
         call xgemm ('N','N',m,n,k,a1,dvvz(ioff1),m,vvoo3(ioff2),k,
     &               a0,vabijz(ioff3),m)
      enddo

      n1 = ABIJTT(desrep)
      write(*,*) 'ICEIL 2:',n1

! debug: original call leading to boundary violations
!     CALL SRT1T2_NTSI (NREP,desrep,MULTB,NVVT,NV,NV,NOOT,
!    &                  JBVVOO,vabijx,vvoop,n1)
! debug end

      CALL SRT1T2_NTSI (NREP,desrep,MULTB,NVVT,NV,NV,NOOT,
     &                  JBVVOO,LLVOOT,vabijx,vvoop,n1)
      call rcaxpy(n1,carith,pphasex,vvoop,fabijtx)
      CALL SRT1T2_NTSI (NREP,desrep,MULTB,NVVT,NV,NV,NOOT,
     &                  JBVVOO,LLVOOT,vabijy,vvoop,n1)
      call rcaxpy(n1,carith,pphasey,vvoop,fabijty)
      CALL SRT1T2_NTSI (NREP,desrep,MULTB,NVVT,NV,NV,NOOT,
     &                  JBVVOO,LLVOOT,vabijz,vvoop,n1)
      call rcaxpy(n1,carith,pphasez,vvoop,fabijtz)


      deallocate(vvoop)
      deallocate(vabijx)
      deallocate(vabijy)
      deallocate(vabijz)
!
!  calculate total transition moment including all contribs.
!
 7766 iorder='0th-2nd order ph + 2p2h A+B'
      incsat = ABIJTT(desrep)

      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LTR,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      deallocate(fabijtx,fabijty,fabijtz)
!______________________________________________________
!|
!|
!|
!|
!|
!|  END CALCULATIONS, Deallocate all remaining arrays
!|
!|
!|
!|
!|_____________________________________________________
!
!
      deallocate(abj,aak,aki,aac)

      deallocate(vvoo3,vvoo4)
      deallocate(vvoo1,vvoo2)

      deallocate(doox,dooy,dooz)
      deallocate(dovx,dovy,dovz)
      deallocate(dvox,dvoy,dvoz)
      deallocate(dvvx,dvvy,dvvz)

      deallocate(faipx,faipy,faipz)
      deallocate(faitx,faity,faitz)

      deallocate(pphase)

      deallocate(ABIJTT)
      deallocate(JAVVOO)
      deallocate(JBVVOO)
      deallocate(MVVOOT)
      deallocate(BVVOOT)

      deallocate(JPOOXX,JPVOXX,JPOVXX,JPVVXX)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE TM_GEN(ioch,desrep,isatctl,iorder,
     &                  do2h2p,
     &                  faitx,faity,faitz,
     &                  fabijtx,fabijty,fabijtz)
      implicit none
!
!  this routine contracts the Davidson eigenvector with the 
!  F_ai/F_aibj matrix elements to form the three transition moments.
!  Remind that the F_ai and F_aibj elements together with the TMs
!  are kept complex!
!
!  due to high numerical effort in the calling routine a flag (do2h2p) 
!  determines if the 2h2p contributions should be included. In this
!  case the arrays fabijt[xyz] must be allocated! Otherwise they are not
!  referenced.
!
!  The routine can be called for each order of perturbation theory in the
!  F_ai matrix elements. It just needs access to the Davidson eigenvectors.
!
!  calling parameters
!
      integer                      :: ioch,desrep,isatctl
      character*30                 :: iorder
      logical                      :: do2h2p
      complex*16, dimension(:)     :: faitx,faity,faitz
      complex*16, dimension(:)     :: fabijtx,fabijty,fabijtz

!
!  includes
!
#include  "../relccsd/complex.inc"
!
!  locals
!
      character*8, parameter                  :: basename='DAVILONG'
      character*11                            :: davevcname
      logical                                 :: isthere

      real*8, allocatable, dimension(:)       :: rcoef
      complex*16, allocatable, dimension(:)   :: ccoef
      complex*16,allocatable,dimension(:)     :: tmx,tmy,tmz

      integer                                 :: ioresult
      integer                                 :: nroots,nmain,ladc
      real*8,allocatable,dimension(:)         :: evl,xnr
      real*8                                  :: ps
      real*8, parameter                       :: D0=0.0d0
      complex*16, parameter                   :: A0=(D0,D0)
      integer                                 :: i,k,ist
      complex*16                              :: ax,ay,az
      real*8, parameter                       :: autoev = 27.2113957D0
      integer                                 :: n2h2p
      character*7,dimension(2)                :: numtype
!
!  executable
!
      if(carith.and.RCW.ne.2) STOP 'Carith error'
      if(.not.carith.and.RCW.ne.1) STOP 'Carith error'
      numtype(1) = 'real'
      numtype(2) = 'complex'
      IF(desrep.GT.9) THEN
        WRITE(davevcname,'(A8,A1,I2)') basename,'.',desrep
      ELSE
        WRITE(davevcname,'(A8,A2,I1)') basename,'.0',desrep
      ENDIF
      INQUIRE(file=davevcname,exist=isthere)
      if(.not.isthere) then
        write(*,*) '***************************************'
        write(*,*) '*****  Warning:                   *****'
        write(*,*) '***************************************'
        write(*,*) 'No eigenvectors found in symmetry',desrep,'!'
        write(*,*) 'No TMs accessible!'
        write(*,*)
        return
      else
        write(*,*)
        write(*,*) 'Found ',numtype(rcw),' eigenvecs in symm.',desrep
      endif
      OPEN(ioch,file=davevcname,access='sequential',
     &     status='unknown',form='unformatted')
      REWIND(ioch)

!_________________________  perform the contractions
!|

      READ(ioch,iostat=ioresult) nroots
      READ(ioch,iostat=ioresult) nmain
      READ(ioch,iostat=ioresult) ladc

      write(*,*) 'Found',nroots,' states in the eigenvector file.'
      write(*,*) 'Length of ADC matrix (from file):',ladc
      write(*,*) 'Length of main space (from file):',nmain
      n2h2p = ladc-nmain
      write(*,*) 'Length of sat space:',n2h2p

      if(do2h2p) then
        write(*,*) 'Including 2h2p contributions to TM.'
        if (isatctl.ne.n2h2p) 
     &    call quit('Internal error for coupling block length')
      else
        write(*,*) 'Omitting 2h2p contributions to TM'
      endif
!
! allocate eigenvalue, transition moment and coefficient arrays
!
      allocate(evl(nroots))
      evl = D0

      allocate(tmx(nroots))
      allocate(tmy(nroots))
      allocate(tmz(nroots))
      tmx = A0
      tmy = A0
      tmz = A0
      if(carith) then
        allocate(ccoef(ladc))
      else
        allocate(rcoef(ladc))
      endif

      allocate(xnr(nroots))
      xnr = D0

      Do k=1,nroots
        read(ioch) ist
        read(ioch) evl(k)
        read(ioch) ps
        write(*,'(A,X,I3,X,A,F10.4,A,X,F10.4,A,F10.4)') 
     &    'State',ist,':exc. energy',evl(k),'(au)',evl(k)*autoev,
     &    ' (eV) and PS',ps
        if(carith) then
          read(ioch)(ccoef(i),i=1,ladc)
          do i=1,ladc
            xnr(k)=xnr(k) + 
     &         (real(ccoef(i))**2.0d0 + aimag(ccoef(i))**2.0d0)
          enddo
        else
          read(ioch)(rcoef(i),i=1,ladc)
          do i=1,ladc
            xnr(k)=xnr(k) + rcoef(i)**2.0d0
          enddo
        endif
        xnr(k) = dsqrt(xnr(k))
        
!
! the do2h2p flag controls the level of contraction.
!
        ax=A0
        ay=A0
        az=A0
        if(carith) then
          do i=1,nmain
            ax = ax + ccoef(i) * faitx(i)
            ay = ay + ccoef(i) * faity(i)
            az = az + ccoef(i) * faitz(i)
          enddo
          if(do2h2p) then
            do i=1,n2h2p
              ax = ax + ccoef(i+nmain) * fabijtx(i)
              ay = ay + ccoef(i+nmain) * fabijty(i)
              az = az + ccoef(i+nmain) * fabijtz(i)
            enddo
          endif
        else
          do i=1,nmain
            ax = ax + dcmplx(rcoef(i),D0)*faitx(i)
            ay = ay + dcmplx(rcoef(i),D0)*faity(i)
            az = az + dcmplx(rcoef(i),D0)*faitz(i)
!           write(*,'(F10.6,5X,2F10.6,5X,2F10.6)') rcoef(i),faitx(i),ax
          enddo
          if(do2h2p) then
            do i=1,n2h2p
              ax = ax + dcmplx(rcoef(i+nmain),D0)*fabijtx(i)
              ay = ay + dcmplx(rcoef(i+nmain),D0)*fabijty(i)
              az = az + dcmplx(rcoef(i+nmain),D0)*fabijtz(i)
            enddo
          endif
        endif
        tmx(k) = ax
        tmy(k) = ay
        tmz(k) = az

      Enddo   ! nroots

      CLOSE(ioch)
!|________________________ end reading section

!_________________________  print output
!|
!|
      write(*,*) '   Transition moments for states in symmetry',desrep
      write(*,*)
      write(*,*)
      write(*,*) '   *******************************'
      write(*,'(A,A)') '   Order of computation: ',iorder
      write(*,*) '   *******************************'
      write(*,*)
      write(*,*)
      write(*,*) '----------------------------------------------------'
      write(*,'(A,4X,A,9X,A,20X,A,20X,A)')
     &       'State','Energy(eV)','TM(x)','TM(y)','TM(z)'
      write(*,*) '----------------------------------------------------'
      write(*,*)
      do k=1,nroots
        write(*,'(I4,F14.7,F8.4,2F15.9,A,2F15.9,A,2F15.9)') k,
     &     evl(k)*autoev,xnr(k),tmx(k),' | ',tmy(k),' | ',tmz(k)
      enddo
      write(*,*) '----------------------------------------------------'
!|
!|________________________ end print output

      if(carith) then
        deallocate(ccoef)
      else
        deallocate(rcoef)
      endif
      deallocate(tmx)
      deallocate(tmy)
      deallocate(tmz)

      deallocate(evl)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
