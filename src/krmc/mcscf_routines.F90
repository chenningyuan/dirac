module mcscf_routines

public :: raddfq90,raddfq90_mlf, m2dnz390, raddfq90_test

private :: get_fq, geth2gu

contains

  subroutine raddfq90_test(fq,h2gu,h2ug,pv,nxy,ix,iy,nz,lupri,norbt,nasht, &
                                     NFSYM,NASH,ipqtoq)
  !***********************************************************************
  !
  !     Add up the infinite number of terms that constitutes the FQ fock
  !     matrix.
  !
  !     For proper documentation see Ph.D. thesis of Joern Thyssen where
  !     all the formulae below should be documented :-)
  !
  !     Formulas (B.40) - (B.43)
  !
  !     Input:
  !        IREPIJ: symmetry of integrals (currently unused)
  !        H2XY, H2YX: (p,v) and (v,p)^T integrals for compound index XY
  !        NXY   : the XY compound index
  !        IX,IY : x and y
  !        PV    : 2-particle density matrix in (NZ,3) format, see appendix B3
  !
  !     Output:
  !        FQ    : the aux. Fock matrix
  !
  !     Written by J. Thyssen - Nov 16 2000
  !
  !***********************************************************************
  implicit none

  ! External variables
  integer, intent(in) :: nxy, ix, iy,ipqtoq(4,0:7)
  integer, intent(in) :: nz,lupri, norbt, nasht
  real*8,  intent(in) :: h2gu(norbt,nasht,nz,nz), h2ug(nasht,norbt,nz,nz)
  real*8,  intent(in) :: pv(2*nasht,2*nasht,2*nasht,2*nasht,2)
  real*8,  intent(inout) :: fq(norbt,nasht,nz)

  ! Internal variables
  integer :: izz, npvcol, IV, IP,IU, IVTR,IUTR, IXTR, IYTR
  INTEGER :: I, J, K, L, M, IJ, II, JJ, KK, LL, ITR, JTR, KTR, LTR
  INTEGER :: IER, NCLASS, IFSYM,NFSYM, IZ, KR(-300:300), NASH(2)
  REAL*8  :: P 

  real*8, parameter   :: d1 = 1.0d0, dm1 = -1.0d0

  !------------------------------------------------------------------------------!
  !                               Print section                                  !
  !------------------------------------------------------------------------------!
!  if(iprint .ge. 20)then
    call header('Output from raddfq90_mlf',-1)
    write(lupri,'(a,i3)') 'nxy       = ',nxy
    write(lupri,'(a,i3)') 'ix        = ',ix
    write(lupri,'(a,i3)') 'iy        = ',iy
    call getpv(pv)
    call get_h2gu(h2gu, h2ug,norbt,nasht,nz)







      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO



  !------------------------------------------------------------------------------!
  !                          *** Real part of FQ ***                             !
  !------------------------------------------------------------------------------!
  ! (B.40): each "term" is a line of the formula

  !--------------------!
  ! Term 1 (real part) !
  !--------------------!

  ! \sum_{x >= y} \sum_v -[pv|xy]_{2,1} Im[P_{qv,xy}]

  IXTR = KR(-IX)
  IYTR = KR(-IY)
  DO IP = 1, NORBT
    DO IU= 1, NASHT
      DO IV = 1, NASHT
        IVTR = KR(-IV)
        IUTR = KR(-IU)
      !  WRITE(6,*) 'IU,IV,IX,IY,P ',IU,IV,IX,IY, PV(IU,IV,IX,IY,1)
      !  WRITE(6,*) 'IU,IV,IYTR,IXTR,P ',IU,IV,IYTR,IXTR,PV(IU,IV,IYTR,IXTR,1)
      !  WRITE(6,*) 'IVTR,IUTR,IX,IY,P ',IVTR,IUTR,IX,IY,PV(IVTR,IUTR,IX,IY,1)
      !  WRITE(6,*) 'IVTR,IUTR,IYTR,IXTR,P ',IVTR,IUTR,IYTR,IXTR,PV(IVTR,IUTR,IYTR,IXTR,1)
      !  WRITE(6,*) 'IX,IY,IU,IV,P ',IX,IY,IU,IV,PV(IX,IY,IU,IV,1)
      !  WRITE(6,*) 'IYTR,IXTR,IU,IV,P ',IYTR,IXTR,IU,IV,PV(IYTR,IXTR,IU,IV,1)
      !  WRITE(6,*) 'IX,IY,IVTR,IUTR,P ',IX,IY,IVTR,IUTR,PV(IX,IY,IVTR,IUTR,1)
      !  WRITE(6,*) 'IYTR,IXTR,IVTR,IUTR,P ',IYTR,IXTR,IVTR,IUTR,PV(IYTR,IXTR,IVTR,IUTR,1)

! correct code
! 1<->3 2<->4
!        P = PV(IU,IV,IX,IY,2) + PV(IU,IV,IYTR,IXTR,2) + PV(IVTR,IUTR,IX,IY,2) &
!       +   PV(IVTR,IUTR,IYTR,IXTR,2) + PV(IX,IY,IU,IV,2)+ PV(IYTR,IXTR,IU,IV,2)&
!       + PV(IX,IY,IVTR,IUTR,2)+ PV(IYTR,IXTR,IVTR,IUTR,2)


! 1<->2 3<->4
!        P = PV(IU,IV,IX,IY,2) + PV(IU,IV,IYTR,IXTR,2) + PV(IVTR,IUTR,IX,IY,2) &
!       +   PV(IVTR,IUTR,IYTR,IXTR,2) &
!       + PV(IV,IU,IY,IX,2) + PV(IV,IU,IXTR,IYTR,2) + PV(IUTR,IVTR,IY,IX,2) &
!       +   PV(IUTR,IVTR,IXTR,IYTR,2) 

! 1<->4 2<->3
        P = PV(IU,IV,IX,IY,2) + PV(IU,IV,IYTR,IXTR,2) + PV(IVTR,IUTR,IX,IY,2) &
       +   PV(IVTR,IUTR,IYTR,IXTR,2) &
       + PV(IY,IX,IV,IU,2) + PV(IXTR,IYTR,IV,IU,2) + PV(IY,IX,IUTR,IVTR,2) &
       +   PV(IXTR,IYTR,IUTR,IVTR,2) 


      !  P = PV(IU,IV,IX,IY,2) + PV(IU,IV,IYTR,IXTR,2) + PV(IVTR,IUTR,IX,IY,2) &
      ! +   PV(IVTR,IUTR,IYTR,IXTR,2) + PV(IX,IY,IU,IV,2)+ PV(IYTR,IXTR,IU,IV,2)&
      ! +  PV(IX,IY,IVTR,IUTR,2) +  PV(IYTR,IXTR,IVTR,IUTR,2)

      
        FQ(IP,IU,1) = FQ(IP,IU,1) - H2GU(IP,IV,1,2)*P
      !  IF (P.GT.1E-12) THEN
 !         WRITE(6,*) 'IU,IV,IX,IY,PTOTAL ', IU,IV,IX,IY,P
       ! ENDIF
  !      WRITE(6,*) 'IU,IV,IX,IY,PV ',IU,IV,IX,IY,P
      ENDDO
    ENDDO
  ENDDO

  if(ix /= iy)then
    ! \sum_{x > y} \sum_v -[vp|xy]_{2,1} Im[P_{vq,xy}]
    DO IP = 1, NORBT 
      DO IU = 1, NASHT 
        DO IV = 1, NASHT 
          IVTR = KR(-IV)
          IUTR = KR(-IU)

      !  WRITE(6,*) 'IV,IU,IX,IY,P ',IV,IU,IX,IY, PV(IV,IU,IX,IY,1)
      !  WRITE(6,*) 'IV,IU,IYTR,IXTR,P ',IV,IU,IYTR,IXTR,PV(IV,IU,IYTR,IXTR,1)
      !  WRITE(6,*) 'IUTR,IVTR,IX,IY,P ',IUTR,IVTR,IX,IY,PV(IUTR,IVTR,IX,IY,1)
      !  WRITE(6,*) 'IUTR,IVTR,IYTR,IXTR,P ',IUTR,IVTR,IYTR,IXTR,PV(IUTR,IVTR,IYTR,IXTR,1)
      !  WRITE(6,*) 'IX,IY,IV,IU,P ',IX,IY,IV,IU,PV(IX,IY,IV,IU,1)
      !  WRITE(6,*) 'IYTR,IXTR,IV,IU,P ',IYTR,IXTR,IV,IU,PV(IYTR,IXTR,IV,IU,1)
      !  WRITE(6,*) 'IX,IY,IUTR,IVTR,P ',IX,IY,IUTR,IVTR,PV(IX,IY,IUTR,IVTR,1)
      !  WRITE(6,*) 'IYTR,IXTR,IUTR,IVTR,P ',IYTR,IXTR,IUTR,IVTR,PV(IYTR,IXTR,IUTR,IVTR,1)
      !    P = PV(IV,IU,IX,IY,1) + PV(IV,IU,IYTR,IXTR,1) + PV(IUTR,IVTR,IX,IY,1) &
      !     + PV(IUTR,IVTR,IYTR,IXTR,1) + PV(IX,IY,IV,IU,1)+ PV(IYTR,IXTR,IV,IU,1)&
      !     + PV(IX,IY,IUTR,IVTR,1) + PV(IYTR,IXTR,IUTR,IVTR,1)

! correct code
! 1<->3 2<->4
!          P = PV(IV,IU,IX,IY,2) + PV(IV,IU,IYTR,IXTR,2) + PV(IUTR,IVTR,IX,IY,2) &
!           + PV(IUTR,IVTR,IYTR,IXTR,2) &
!           + PV(IX,IY,IV,IU,2)+ PV(IYTR,IXTR,IUTR,IVTR,2)&
!           +  PV(IYTR,IXTR,IV,IU,2) +  PV(IX,IY,IUTR,IVTR,2)
!




! 1<->2 3<->4
!          P = PV(IV,IU,IX,IY,2) + PV(IV,IU,IYTR,IXTR,2) + PV(IUTR,IVTR,IX,IY,2)&
!           + PV(IUTR,IVTR,IYTR,IXTR,2) &
!           + PV(IU,IV,IY,IX,2) + PV(IU,IV,IXTR,IYTR,2) + PV(IVTR,IUTR,IY,IX,2)&
!           + PV(IVTR,IUTR,IXTR,IYTR,2) 

! 1<->4 2<->3
          P = PV(IV,IU,IX,IY,2) + PV(IV,IU,IYTR,IXTR,2) + PV(IUTR,IVTR,IX,IY,2)&
           + PV(IUTR,IVTR,IYTR,IXTR,2) &
           + PV(IY,IX,IU,IV,2) + PV(IXTR,IYTR,IU,IV,2) + PV(IY,IX,IVTR,IUTR,2)&
           + PV(IXTR,IYTR,IVTR,IUTR,2) 




          FQ(IP,IU,1) = FQ(IP,IU,1) - H2UG(IV,IP,1,2)*P
      !  IF (P.GT.1E-12) THEN
  !        WRITE(6,*) 'IV,IU,IX,IY,PTOTAL ', IV,IU,IX,IY,P
       ! ENDIF
        ENDDO
      ENDDO
    ENDDO
  WRITE(6,*) 'OUTPUT FQ_TEST ON EXIT raddfq90_mlf'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)

  ENDIF





  end subroutine raddfq90_test

  subroutine raddfq90_mlf(fq,h2full,pvmolf,nz,lupri,norbt,nasht,maxash,idxu2g)
  !***********************************************************************
  !
  !     Add up the infinite number of terms that constitutes the FQ fock
  !     matrix.
  !
  !     For proper documentation see Ph.D. thesis of Joern Thyssen where
  !     all the formulae below should be documented :-)
  !
  !     Formulas (B.40) - (B.43)
  !
  !     Input:
  !        IREPIJ: symmetry of integrals (currently unused)
  !        H2XY, H2YX: (p,v) and (v,p)^T integrals for compound index XY
  !        NXY   : the XY compound index
  !        IX,IY : x and y
  !        PV    : 2-particle density matrix in (NZ,3) format, see appendix B3
  !
  !     Output:
  !        FQ    : the aux. Fock matrix
  !
  !     Written by J. Thyssen - Nov 16 2000
  !
  !***********************************************************************
  implicit none

  ! External variables
  integer, intent(in) :: nz, norbt, nasht,maxash,IDXU2G(MAXASH),lupri
  real*8,  intent(in) :: H2FULL(2*NORBT,2*NORBT,2*NASHT,2*NASHT,2)
  real*8,  intent(in) :: PVMOLF(2*NASHT,2*NASHT,2*NASHT,2*NASHT,2)
  real*8,  intent(inout) :: fq(2*norbt,2*nasht,2)

  ! Internal variables
  real*8 :: H2GAC(2*NORBT,2*NASHT,2*NASHT,2*NASHT,2)
  real*8 :: PVMOLF_SYM(2*NASHT,2*NASHT,2*NASHT,2*NASHT,2)
  integer :: IG, IV, IVG, IX, IY, IU, IZ, P, Q
  real*8  :: trace(2)
  real*8, parameter   :: d1 = 1.0d0, dm1 = -1.0d0, d0 = 0.0d0

  !------------------------------------------------------------------------------!
  !                               Print section                                  !
  !------------------------------------------------------------------------------!


  !    II = 0
  !    JJ = 0
  !    DO IFSYM = 1, NFSYM
  !       DO I = 1, NASH(IFSYM)
  !          II = II + 1
  !          JJ = JJ + 1
  !          KR(II) = JJ
  !          KR(-II) = JJ + NASH(IFSYM)
  !       ENDDO
  !       JJ = JJ + NASH(IFSYM)
  !    ENDDO



    DO IG = 1, 2*NORBT
      DO IV = 1, NASHT
        IVG = IDXU2G(IV)
        DO IX = 1, 2*NASHT
          DO IY = 1, 2*NASHT
            DO IZ = 1, MIN(2,NZ)
              H2GAC(IG,IV,IX,IY,IZ) = H2FULL(IG,IVG,IX,IY,IZ) 
              H2GAC(IG,IV+NASHT,IX,IY,IZ) = H2FULL(IG,IVG+NORBT,IX,IY,IZ) 
              !WRITE(LUPRI,*) 'H2GAC(IG,IV,IX,IY,IZ) ', &
              !IG,IV,IX,IY,IZ,H2GAC(IG,IV,IX,IY,IZ)
              !WRITE(LUPRI,*) 'H2GAC(IG,IV+NASHT,IX,IY,IZ) '&
              !,IG,IV+NASHT,IX,IY,IZ,H2GAC(IG,IV+NASHT,IX,IY,IZ)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO IV = 1, 2*NASHT
      DO IU = 1, 2*NASHT
        DO IX = 1, 2*NASHT
          DO IY = 1, 2*NASHT 
            PVMOLF_SYM(IV,IU,IX,IY,1) = PVMOLF(IX,IY,IV,IU,1) + PVMOLF(IV,IU,IX,IY,1)
            PVMOLF_SYM(IV,IU,IX,IY,2) = PVMOLF(IX,IY,IV,IU,2) + PVMOLF(IV,IU,IX,IY,2)
          ENDDO
        ENDDO
      ENDDO
    ENDDO


     CALL DZERO(FQ, 8*norbt*nasht)
 !   CALL GETPV(PVMOLF)
    DO IU = 1, 2*NASHT
      DO IG = 1, 2*NORBT
        DO IV = 1, 2*NASHT
          DO IX = 1, 2*NASHT
            DO IY = 1, 2*NASHT
              FQ(IG,IU,1) = FQ(IG,IU,1) + PVMOLF_SYM(IU,IV,IX,IY,1)*H2GAC(IG,IV,IX,IY,1) &
                           -PVMOLF_SYM(IU,IV,IX,IY,2)*H2GAC(IG,IV,IX,IY,2) 
              FQ(IG,IU,2) = FQ(IG,IU,2) + PVMOLF_SYM(IU,IV,IX,IY,2)*H2GAC(IG,IV,IX,IY,1) &
                           +PVMOLF_SYM(IU,IV,IX,IY,1)*H2GAC(IG,IV,IX,IY,2)
            !  write(lupri,*) 'IU,IG,IV,IX,IY', IU,IG,IV,IX,IY
            !  write(lupri,*) 'PVMOLF(IU,IV,IX,IY,1)*H2GAC(IG,IV,IX,IY,1)'&
            !                  ,PVMOLF_SYM(IU,IV,IX,IY,1)*H2GAC(IG,IV,IX,IY,1)
            !  write(lupri,*) '-PVMOLF(IU,IV,IX,IY,1)*H2GAC(IG,IV,IX,IY,1)'&
            !                  ,-PVMOLF_SYM(IU,IV,IX,IY,2)*H2GAC(IG,IV,IX,IY,2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO 
 9001 FORMAT(4I5,5X,4F20.10)

!        WRITE(LUPRI,*) 'RADDFQ90_MLF: PRINT H2GAC'
!         DO IU = 1,2*NASHT
!            DO IV = 1,2*NASHT
!               DO IX = 1,NASHT
!                  DO IY = 1,2*NORBT
!                    P = IDXU2G(IY)
!                    Q = IDXU2G(IX)
!                         WRITE(LUPRI,9001) &
!                         P,Q,IV,IU,(H2FULL(IY,Q,IV,IU,IZ),IZ=1,MIN(NZ,2))
!                  END DO
!               END DO
!            END DO
!         END DO
!



    trace = d0

    DO IV = 1, 2*NASHT
      IVG = IDXU2G(IV)
      TRACE(1) = TRACE(1) + FQ(IVG,IV,1)
      TRACE(2) = TRACE(2) + FQ(IVG,IV,2)
    ENDDO

    WRITE(LUPRI,*) 'OUTPUT: RADDFQ90_MLF'
    WRITE(LUPRI,*) 'The trace of ReFQ is: ', TRACE(1)
    WRITE(LUPRI,*) 'The trace of ImFQ is: ', TRACE(2)

    WRITE(LUPRI,*) 'BELOW IS FQ'
    DO IG = 1, 2*NORBT
      DO IV = 1, 2*NASHT
        WRITE(LUPRI,*) 'FQ, IG, IV, ',IG,IV,(FQ(IG,IV,IZ),IZ=1,2)
      ENDDO
    ENDDO




  end subroutine raddfq90_mlf

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
  subroutine raddfq90(fq,h2xy,h2yx,pv,nxy,ix,iy,irepij,iprint, &
                      nz,ipqtoq,lupri,norbt,nasht,nnashx,key)
  !***********************************************************************
  !
  !     Add up the infinite number of terms that constitutes the FQ fock
  !     matrix.
  !
  !     For proper documentation see Ph.D. thesis of Joern Thyssen where
  !     all the formulae below should be documented :-)
  !
  !     Formulas (B.40) - (B.43)
  !
  !     Input:
  !        IREPIJ: symmetry of integrals (currently unused)
  !        H2XY, H2YX: (p,v) and (v,p)^T integrals for compound index XY
  !        NXY   : the XY compound index
  !        IX,IY : x and y
  !        PV    : 2-particle density matrix in (NZ,3) format, see appendix B3
  !
  !     Output:
  !        FQ    : the aux. Fock matrix
  !
  !     Written by J. Thyssen - Nov 16 2000
  !
  !***********************************************************************
  implicit none

  ! External variables
  integer, intent(in), optional :: key
  integer, intent(in) :: nxy, ix, iy, irepij, iprint
  integer, intent(in) :: nz, ipqtoq(4,0:7), lupri, norbt, nasht, nnashx
  real*8,  intent(in) :: h2xy(norbt,nasht,nz,3), h2yx(norbt,nasht,nz,3)
  real*8,  intent(in) :: pv(nasht,nasht,nnashx,nz,3)
  real*8,  intent(inout) :: fq(norbt,nasht,nz)

  ! Internal variables
  integer :: i, j, k, izz, npvcol, ier

  real*8, parameter   :: d1 = 1.0d0, dm1 = -1.0d0

  !------------------------------------------------------------------------------!
  !                               Print section                                  !
  !------------------------------------------------------------------------------!
  if(iprint .ge. 20)then
    call header('Output from raddfq90',-1)
    write(lupri,'(a,i3)') 'nxy       = ',nxy
    write(lupri,'(a,i3)') 'ix        = ',ix
    write(lupri,'(a,i3)') 'iy        = ',iy
    write(lupri,'(a,i3)') 'irepij    = ',irepij
    write(lupri,'(/a)') ' fq on entry'
    call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
    npvcol = nasht * nnashx
    do i = 1, 3
      write(lupri,'(a,i2,a)') ' pv density matrix',i,' of 3'
      call prqmat(pv(1,1,nxy,1,i),nasht,nasht,nasht,npvcol,nz,ipqtoq(1,0),lupri)
    end do
    do i = 1, 3
      write(lupri,'(a,i2,a)') ' h2xy matrix',i,' of 3'
      call prqmat(h2xy(1,1,1,i),norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
    end do
    if (ix /= iy) then
      do i = 1, 3
        write(lupri,'(a,i2,a)') ' h2yx matrix',i,' of 3'
        call prqmat(h2yx(1,1,1,i),norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
      end do
    end if
  end if

  !------------------------------------------------------------------------------!
  !                          *** Real part of FQ ***                             !
  !------------------------------------------------------------------------------!
  ! (B.40): each "term" is a line of the formula

  !--------------------!
  ! Term 1 (real part) !
  !--------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;1,1}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,1,1,1,1)

  if(ix /= iy)then
    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;1,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,d1,1,1,1,1,1)
  endif
  write(lupri,*) 'term1 FQ ix, iy', ix, iy
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)

  !--------------------!
  ! Term 2 (real part) !
  !--------------------!

  if(nz >= 2)then

    ! - \sum_{x >= y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;2,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,2,1,2,1,1)

    if(ix /= iy)then
      ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;2,1}
      call get_fq('n','n',h2yx,pv,fq,nxy,dm1,2,1,2,1,1)
    endif

  endif
  write(lupri,*) 'term2 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  !--------------------------!
  ! Term 3 and 4 (real part) !
  !--------------------------!

  if(nz == 4 .and. ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;3,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,1,3,1,1)
  write(lupri,*) 'term3 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;4,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,4,1,4,1,1)
  write(lupri,*) 'term4 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  endif

  !--------------------!
  ! Term 5 (real part) !
  !--------------------!

  if(ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,2} P_{qv,xy;1,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,2,1,2,1)
  write(lupri,*) 'term5 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  endif

  !--------------------!
  ! Term 6 (real part) !
  !--------------------!

  if(nz >= 2 .and. ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,2} P_{qv,xy;2,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,2,2,2,2,1)
  write(lupri,*) 'term6 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  endif

  if(nz == 4)then

  !--------------------!
  ! Term 7 (real part) !
  !--------------------!

    ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{3,2} P_{qv,xy;3,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,2,3,2,1)

    if(ix /= iy)then

      ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,3} P_{vq,xy;3,3}
      call get_fq('n','n',h2yx,pv,fq,nxy,d1,3,3,3,3,1)

    endif
  write(lupri,*) 'term7 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  !--------------------!
  ! Term 8 (real part) !
  !--------------------!

    ! - \sum_{x >= y} \sum_v \Twoint{pv|xy}_{4,2} P_{qv,xy;4,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,4,2,4,2,1)

    if(ix /= iy)then

      ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,3} P_{vq,xy;4,3}
      call get_fq('n','n',h2yx,pv,fq,nxy,dm1,4,3,4,3,1)

    endif
  write(lupri,*) 'term8 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  endif

  !--------------------!
  ! Term 9 (real part) !
  !--------------------!

  if(ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;1,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,3,1,3,1)

  endif
  write(lupri,*) 'term9 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  !---------------------!
  ! Term 10 (real part) !
  !---------------------!

  if(nz >= 2 .and. ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;2,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,2,3,2,3,1)

  endif
  write(lupri,*) 'term10 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  if(nz == 4 .and. ix /= iy)then

  !---------------------!
  ! Term 11 (real part) !
  !---------------------!

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,1} P_{vq,xy;3,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,d1,3,1,3,1,1)
  write(lupri,*) 'term11 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  !---------------------!
  ! Term 12 (real part) !
  !---------------------!

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,1} P_{vq,xy;4,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,4,1,4,1,1)
  write(lupri,*) 'term12 FQ'
  call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  endif

  ! Skip remaining 36 terms if NZ = 1.
  if(nz > 1)then ! goto 100

  !------------------------------------------------------------------------------!
  !                         *** i-imag part of FQ ***                            !
  !------------------------------------------------------------------------------!
  ! (B.41)

  !----------------------!
  ! Term 1 (i-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;2,1}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,1,2,1,2)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;2,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,1,1,2,1,2)

  endif

  !----------------------!
  ! Term 2 (i-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;1,1}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,2,1,1,1,2)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;1,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,2,1,1,1,2)

  endif

  !----------------------------!
  ! Term 3 and 4 (i-imag part) !
  !----------------------------!

  if(nz == 4 .and. ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;4,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,1,4,1,2)

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;3,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,4,1,3,1,2)

  endif

  if(ix /= iy)then

  !----------------------!
  ! Term 5 (i-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,2} P_{qv,xy;2,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,2,2,2,2)

  !----------------------!
  ! Term 6 (i-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,2} P_{qv,xy;1,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,2,2,1,2,2)

  endif

  if(nz == 4)then

  !----------------------!
  ! Term 7 (i-imag part) !
  !----------------------!

    ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{3,2} P_{qv,xy;4,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,2,4,2,2)

    if(ix /= iy)then

      ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,3} P_{vq,xy;4,3}
      call get_fq('n','n',h2yx,pv,fq,nxy,dm1,3,3,4,3,2)

    endif

  !----------------------!
  ! Term 8 (i-imag part) !
  !----------------------!

    ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{4,2} P_{qv,xy;3,2}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,4,2,3,2,2)

    if(ix /= iy)then

      ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,3} P_{vq,xy;3,3}
      call get_fq('n','n',h2yx,pv,fq,nxy,dm1,4,3,3,3,2)

    endif

  endif

  if(ix /= iy)then

  !----------------------!
  ! Term 9 (i-imag part) !
  !----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;2,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,1,3,2,3,2)

  !-----------------------!
  ! Term 10 (i-imag part) !
  !-----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;1,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,2,3,1,3,2)

  endif

  if(nz == 4 .and. ix /= iy)then

  !-----------------------!
  ! Term 11 (i-imag part) !
  !-----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,1} P_{vq,xy;4,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,3,1,4,1,2)

  !-----------------------!
  ! Term 12 (i-imag part) !
  !-----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,1} P_{vq,xy;3,1}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,4,1,3,1,2)

  endif

  ! Skip remaining 24 terms if NZ = 2.
  if(nz > 2)then ! goto 101

  !------------------------------------------------------------------------------!
  !                         *** j-imag part of FQ ***                            !
  !------------------------------------------------------------------------------!
  ! (B.42)

  !----------------------!
  ! Term 1 (j-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;3,3}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,1,3,3,3)

  if(ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;3,2}
    call get_fq('n','n',h2yx,pv,fq,nxy,d1,1,1,3,2,3)

  endif

  !----------------------!
  ! Term 2 (j-imag part) !
  !----------------------!

  ! - \sum_{x >= y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;4,3}
  call get_fq('n','t',h2xy,pv,fq,nxy,dm1,2,1,4,3,3)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;4,2}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,2,1,4,2,3)

  endif

  !----------------------!
  ! Term 3 (j-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{3,2} P_{vq,xy;1,1}
  call get_fq('n','n',h2xy,pv,fq,nxy,d1,3,2,1,1,3)

  if(ix /= iy)then

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,3} P_{qv,xy;1,1}
    call get_fq('n','t',h2yx,pv,fq,nxy,d1,3,3,1,1,3)

  endif

  !----------------------!
  ! Term 4 (j-imag part) !
  !----------------------!

  ! - \sum_{x >= y} \sum_v \Twoint{pv|xy}_{4,2} P_{vq,xy;2,1}
  call get_fq('n','n',h2xy,pv,fq,nxy,dm1,4,2,2,1,3)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,3} P_{qv,xy;2,1}
    call get_fq('n','t',h2yx,pv,fq,nxy,dm1,4,3,2,1,3)

  endif

  if(ix /= iy)then

  !----------------------!
  ! Term 5 (j-imag part) !
  !----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;3,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,1,3,3,1,3)

  !----------------------!
  ! Term 6 (j-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;4,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,2,3,4,1,3)

  !----------------------!
  ! Term 7 (j-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,2} P_{vq,xy;3,1}
    call get_fq('n','n',h2xy,pv,fq,nxy,d1,1,2,3,1,3)

  !----------------------!
  ! Term 8 (j-imag part) !
  !----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,2} P_{vq,xy;4,1}
    call get_fq('n','n',h2xy,pv,fq,nxy,dm1,2,2,4,1,3)

  !----------------------!
  ! Term 9 (j-imag part) !
  !----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,1} P_{qv,xy;1,2}
    call get_fq('n','t',h2yx,pv,fq,nxy,dm1,3,1,1,2,3)

  !-----------------------!
  ! Term 10 (j-imag part) !
  !-----------------------!

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,1} P_{qv,xy;2,2}
    call get_fq('n','t',h2yx,pv,fq,nxy,d1,4,1,2,2,3)

  !-----------------------!
  ! Term 11 (j-imag part) !
  !-----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;1,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,1,1,3,3)

  !-----------------------!
  ! Term 12 (j-imag part) !
  !-----------------------!

    ! - \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;2,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,dm1,4,1,2,3,3)

  endif

  !------------------------------------------------------------------------------!
  !                         *** k-imag part of FQ ***                            !
  !------------------------------------------------------------------------------!
  ! (B.43)

  !----------------------!
  ! Term 1 (k-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;4,3}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,1,4,3,4)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;4,2}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,1,1,4,2,4)

  endif

  !----------------------!
  ! Term 2 (k-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;3,3}
  call get_fq('n','t',h2xy,pv,fq,nxy,d1,2,1,3,3,4)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;3,2}
    call get_fq('n','n',h2yx,pv,fq,nxy,dm1,2,1,3,2,4)

  endif

  !----------------------!
  ! Term 3 (k-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{3,2} P_{vq,xy;2,1}
  call get_fq('n','n',h2xy,pv,fq,nxy,d1,3,2,2,1,4)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,3} P_{qv,xy;2,1}
    call get_fq('n','t',h2yx,pv,fq,nxy,dm1,3,3,2,1,4)

  endif

  !----------------------!
  ! Term 4 (k-imag part) !
  !----------------------!

  ! \sum_{x >= y} \sum_v \Twoint{pv|xy}_{4,2} P_{vq,xy;1,1}
  call get_fq('n','n',h2xy,pv,fq,nxy,d1,4,2,1,1,4)

  if(ix /= iy)then

    ! - \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,3} P_{qv,xy;1,1}
    call get_fq('n','t',h2yx,pv,fq,nxy,dm1,4,3,1,1,4)

  endif

  if(ix /= iy)then

  !----------------------!
  ! Term 5 (k-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;4,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,1,3,4,1,4)

  !----------------------!
  ! Term 6 (k-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;3,1}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,2,3,3,1,4)

  !----------------------!
  ! Term 7 (k-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{1,2} P_{vq,xy;4,1}
    call get_fq('n','n',h2xy,pv,fq,nxy,d1,1,2,4,1,4)

  !----------------------!
  ! Term 8 (k-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{2,2} P_{vq,xy;3,1}
    call get_fq('n','n',h2xy,pv,fq,nxy,d1,2,2,3,1,4)

  !----------------------!
  ! Term 9 (k-imag part) !
  !----------------------!

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{3,1} P_{qv,xy;2,2}
    call get_fq('n','t',h2yx,pv,fq,nxy,d1,3,1,2,2,4)

  !-----------------------!
  ! Term 10 (k-imag part) !
  !-----------------------!

    ! \sum_{x > y} \sum_v \Twoint{vp|xy}_{4,1} P_{qv,xy;1,2}
    call get_fq('n','t',h2yx,pv,fq,nxy,d1,4,1,1,2,4)

  !-----------------------!
  ! Term 11 (k-imag part) !
  !-----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;2,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,3,1,2,3,4)

  !-----------------------!
  ! Term 12 (k-imag part) !
  !-----------------------!

    ! \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;1,3}
    call get_fq('n','t',h2xy,pv,fq,nxy,d1,4,1,1,3,4)

  endif

  endif ! 101
  endif ! 100

  !----------------!
  ! Output section !
  !----------------!

 ! mvh uncomment later if (iprint .ge. 20) then
    write(lupri,'(a)') 'fq on exit'
    call prqmat(fq,norbt,nasht,norbt,nasht,nz,ipqtoq(1,0),lupri)
  !end if

  !------------------------------------------------------------------------------!
  !           Write to the file ouput information for unit testing               !
  !------------------------------------------------------------------------------!
  if( present(key) )then

    ! open file
    open(60,file='raddfq90.out',status='unknown',form='formatted', &
                                access='sequential',iostat=ier)
    if(ier /= 0)then
      write(6,*) 'raddfq90: ERROR while creating raddfq90.out file'
      stop
    endif

    rewind 60

    ! write output data to the file
    do i=1,nz
    do j=1,nasht
    do k=1,norbt
    write(60,*) fq(k,j,i)
    enddo
    enddo
    enddo

    ! close file
    close(60,status='keep')

  endif

  end subroutine raddfq90
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_fq(ch1,ch2,h2,pv,fq,nxy,alpha,izh,ich,izp,icp,izf)
  implicit none

  ! External variables
  character, intent( in) :: ch1*1, ch2*1
  integer,   intent( in) :: nxy, izh, ich, izp, icp, izf
  real*8,    intent( in) :: alpha
  real*8,    intent( in) :: h2(:,:,:,:), pv(:,:,:,:,:)
  real*8,    intent(out) :: fq(:,:,:)

  ! Internal variables
  integer :: norbt, nasht

  norbt = size(h2,1)
  nasht = size(h2,2)

  call dgemm(ch1,ch2,norbt,nasht,nasht,alpha, &
       h2(1,1,izh,ich),norbt,                 &
       pv(1,1,nxy,izp,icp),nasht,             &
       1.0d0,fq(1,1,izf),norbt)

  end subroutine get_fq
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_h2gu(h2gu, h2ug,norbt,nasht,nz)
  implicit none

  ! External variables
  real*8, intent(in) :: h2gu(norbt,nasht,nz,nz), h2ug(nasht,norbt,nz,nz)
  integer, intent(in) :: nz
  integer :: IP, IU, IZ1, IZ2


  ! Internal variables
  integer :: norbt, nasht

  DO IP = 1, NORBT 
    DO IU = 1, NASHT 
      DO IZ1 = 1, NZ 
        DO IZ2 = 1, NZ     
          WRITE(6,*) 'IP,IU,H2GU,IZ1,IZ2 ', IP, IU, IZ1, IZ2, H2GU(IP, IU, IZ1, IZ2)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

    DO IP = 1, NORBT 
    DO IU = 1, NASHT 
      DO IZ1 = 1, NZ 
        DO IZ2 = 1, NZ     
          WRITE(6,*) 'IP,IU,H2UG,IZ1,IZ2 ', IP, IU, IZ1, IZ2, H2GU(IU, IP, IZ1, IZ2)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  end subroutine get_h2gu

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE M2DNZ390(TSYMM,KSYMM,PSYMM,VMUUUU,VDUUUU,MMOLF,IPRINT,   &
                      IPQTOQ,NASH,NASHT,NNASHX,NZ,NZ_IN_CI,NM5, &
                      NBSYM,NFSYM,LUPRI,SPINFR,LEVYLE)
!***********************************************************************
!
!     Transform a UUUU matrix in Molfdir format to a UUUU matrix
!     in Dirac (NZ,3) format.
!
!     Input:
!        TSYMM: transition density matrix symmetrization:
!           do <0|...|B> + <B|...0> symmetrization
!        KSYMM: do Kramers symmetrization (i.e., calculate P++)
!        PSYMM: particle symmetrization: pqrs + rspq
!
!     Written by J. Thyssen - Jul 21 2000
!     Last revision :
!
!***********************************************************************
!     USE UNIT_TEST_GENERATOR
      IMPLICIT NONE

      LOGICAL, INTENT( IN) :: PSYMM, KSYMM, TSYMM, SPINFR, LEVYLE
      INTEGER, INTENT( IN) :: NASHT, NNASHX, NZ, NZ_IN_CI, NM5, NBSYM, NFSYM, IPRINT, LUPRI
      INTEGER, INTENT( IN) :: NASH(2), IPQTOQ(4,0:7)

      REAL*8,  INTENT( IN) :: VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,NM5)
      REAL*8,  INTENT(OUT) :: VDUUUU(NASHT,NASHT,NNASHX,NZ,3)
      REAL*8,  INTENT(OUT) :: MMOLF(2*NASHT,2*NASHT,2*NASHT,2*NASHT,NM5)

      INTEGER :: I, J, K, L, M, IJ, II, JJ, KK, LL, ITR, JTR, KTR, LTR
      INTEGER :: IER, NCLASS, IFSYM, IZ, KR(-300:300)
      REAL*8  :: THR_PRINT

  !------------------------------------------------------------------------------!
  !           Write to the file input information for unit testing               !
  !------------------------------------------------------------------------------!

#ifdef CREATE_UNIT_TEST_M2DNZ390
      call set_suffix(psymm)
      call set_suffix(ksymm)
      call set_suffix(tsymm)
      call set_suffix(icount)
      call set_driver_routine('m2dnz390')
      call add_module('mcscf_routines')

      call add_inp_variable(psymm ,'psymm ')
      call add_inp_variable(ksymm ,'ksymm ')
      call add_inp_variable(tsymm ,'tsymm ')
      call add_inp_variable(spinfr,'spinfr')
      call add_inp_variable(levyle,'levyle')
      call add_inp_variable(nasht ,'nasht ')
      call add_inp_variable(nnashx,'nnashx')
      call add_inp_variable(nz    ,'nz    ')
      call add_inp_variable(nz_in_ci,'nz_in_ci')
      call add_inp_variable(nm5   ,'nm5   ')
      call add_inp_variable(nbsym ,'nbsym ')
      call add_inp_variable(nfsym ,'nfsym ')
      call add_inp_variable(iprint,'iprint')
      call add_inp_variable(lupri ,'lupri ')
      call add_inp_variable(nash  ,'nash  ')
      call add_inp_variable(ipqtoq,'ipqtoq')
      call add_inp_variable(vmuuuu,'vmuuuu')
#endif

      write(6,*) 'in M2DNZ390',psymm,ksymm,tsymm,NM5

  !------------------------------------------------------------------------------!
      NCLASS = NZ * NZ * NBSYM / NFSYM
      CALL DZERO(VDUUUU,NASHT * NASHT * NNASHX * NZ * 3)

      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO

 9001 FORMAT(4I5,5X,2F20.10)

!      IF (IPRINT .GE. 30) THEN
         THR_PRINT = 1.0D-12
         WRITE(LUPRI,'(/A,3(1X,L1)//A)') &
              ' (M2DNZ3) TSYMM, KSYMM, PSYMM = ',TSYMM,KSYMM,PSYMM, &
              ' (M2DNZ3) Input VUUUU matrix in Molfdir format'
         DO L = 1,2*NASHT
            DO K = 1,2*NASHT
               DO J = 1,2*NASHT
                  DO I = 1,2*NASHT
                   !  IF ( (ABS(VMUUUU(I,J,K,L,1)) .GT. THR_PRINT   &
                    !  .OR. ABS(VMUUUU(I,J,K,L,2)) .GT. THR_PRINT)) &
                          WRITE(LUPRI,9001) &
                          I,J,K,L,(VMUUUU(I,J,K,L,IZ),IZ=1,NZ_IN_CI)
                  END DO
               END DO
            END DO
         END DO
!      END IF

      CALL DZERO(MMOLF,NASHT*32)
       WRITE(LUPRI,'(/A,3(1X,L1)//A)') &
              ' (M2DNZ3) TSYMM, KSYMM, PSYMM = ',TSYMM,KSYMM,PSYMM, &
              ' (M2DNZ3) Input MMOLF matrix '
      DO L = 1,2*NASHT
         DO K = 1,2*NASHT
            DO J = 1,2*NASHT
               DO I = 1,2*NASHT
                  DO IZ = 1, NZ_IN_CI
                    MMOLF(I,J,K,L,IZ) = VMUUUU(I,J,K,L,IZ) + VMUUUU(K,L,I,J,IZ)
                  ENDDO
                  WRITE(LUPRI,9001) &
                          I,J,K,L,(MMOLF(I,J,K,L,IZ),IZ=1,NZ_IN_CI)
               END DO
            END DO
         END DO
      END DO



      IJ = 0
      DO II = 1, NASHT
         DO JJ = 1, II
            IJ = IJ + 1
            DO LL = 1, NASHT
               DO KK = 1, NASHT
                  I = KR(II)
                  J = KR(JJ)
                  K = KR(KK)
                  L = KR(LL)
                  ITR = KR(-II)
                  JTR = KR(-JJ)
                  KTR = KR(-KK)
                  LTR = KR(-LL)

!                 Calculate integral classes.

                  IF (TSYMM) THEN
                     IF (PSYMM) THEN
                        IF (KSYMM) THEN

!                          (1) Symm. transition density matrix
!                          (2) Do particle symmetrization
!                          (3) Do Kramer's symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) + VMUUUU(L,K,J,I,1) &
                                + VMUUUU(I,J,K,L,1) + VMUUUU(J,I,L,K,1) &
                                + VMUUUU(LTR,KTR,I,J,1) &
                                + VMUUUU(KTR,LTR,J,I,1) &
                                + VMUUUU(I,J,LTR,KTR,1) &
                                + VMUUUU(J,I,KTR,LTR,1) &
                                + VMUUUU(K,L,JTR,ITR,1) &
                                + VMUUUU(L,K,ITR,JTR,1) &
                                + VMUUUU(JTR,ITR,K,L,1) &
                                + VMUUUU(ITR,JTR,L,K,1) &
                                + VMUUUU(LTR,KTR,JTR,ITR,1) &
                                + VMUUUU(KTR,LTR,ITR,JTR,1) &
                                + VMUUUU(JTR,ITR,LTR,KTR,1) &
                                + VMUUUU(ITR,JTR,KTR,LTR,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(I,JTR,K,LTR,1) &
                                + VMUUUU(LTR,K,JTR,I,1) &
                                + VMUUUU(JTR,I,LTR,K,1) &
                                - VMUUUU(L,KTR,I,JTR,1) &
                                - VMUUUU(I,JTR,L,KTR,1) &
                                - VMUUUU(KTR,L,JTR,I,1) &
                                - VMUUUU(JTR,I,KTR,L,1) &
                                - VMUUUU(K,LTR,J,ITR,1) &
                                - VMUUUU(J,ITR,K,LTR,1) &
                                - VMUUUU(LTR,K,ITR,J,1) &
                                - VMUUUU(ITR,J,LTR,K,1) &
                                + VMUUUU(L,KTR,J,ITR,1) &
                                + VMUUUU(J,ITR,L,KTR,1) &
                                + VMUUUU(KTR,L,ITR,J,1) &
                                + VMUUUU(ITR,J,KTR,L,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(I,JTR,KTR,L,1) &
                                + VMUUUU(L,KTR,JTR,I,1) &
                                + VMUUUU(JTR,I,L,KTR,1) &
                                - VMUUUU(LTR,K,I,JTR,1) &
                                - VMUUUU(I,JTR,LTR,K,1) &
                                - VMUUUU(K,LTR,JTR,I,1) &
                                - VMUUUU(JTR,I,K,LTR,1) &
                                - VMUUUU(KTR,L,J,ITR,1) &
                                - VMUUUU(J,ITR,KTR,L,1) &
                                - VMUUUU(L,KTR,ITR,J,1) &
                                - VMUUUU(ITR,J,L,KTR,1) &
                                + VMUUUU(LTR,K,J,ITR,1) &
                                + VMUUUU(J,ITR,LTR,K,1) &
                                + VMUUUU(K,LTR,ITR,J,1) &
                                + VMUUUU(ITR,J,K,LTR,1)
                        ELSE

!                          (1) Symm. transition density matrix
!                          (2) Do particle symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) + VMUUUU(L,K,J,I,1) &
                                + VMUUUU(I,J,K,L,1) + VMUUUU(J,I,L,K,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(LTR,K,JTR,I,1) &
                                + VMUUUU(I,JTR,K,LTR,1) &
                                + VMUUUU(JTR,I,LTR,K,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(L,KTR,JTR,I,1) &
                                + VMUUUU(I,JTR,KTR,L,1) &
                                + VMUUUU(JTR,I,L,KTR,1)
                        END IF
                     ELSE
                        IF (KSYMM) THEN

!                          (1) Symm. transition density matrix
!                          (2) Kramers symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) &
                                + VMUUUU(L,K,J,I,1) &
                                + VMUUUU(LTR,KTR,I,J,1) &
                                + VMUUUU(KTR,LTR,J,I,1) &
                                + VMUUUU(K,L,JTR,ITR,1) &
                                + VMUUUU(L,K,ITR,JTR,1) &
                                + VMUUUU(LTR,KTR,JTR,ITR,1) &
                                + VMUUUU(KTR,LTR,ITR,JTR,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(LTR,K,JTR,I,1) &
                                - VMUUUU(L,KTR,I,JTR,1) &
                                - VMUUUU(KTR,L,JTR,I,1) &
                                - VMUUUU(K,LTR,J,ITR,1) &
                                - VMUUUU(LTR,K,ITR,J,1) &
                                + VMUUUU(L,KTR,J,ITR,1) &
                                + VMUUUU(KTR,L,ITR,J,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(L,KTR,JTR,I,1) &
                                - VMUUUU(LTR,K,I,JTR,1) &
                                - VMUUUU(K,LTR,JTR,I,1) &
                                - VMUUUU(KTR,L,J,ITR,1) &
                                - VMUUUU(L,KTR,ITR,J,1) &
                                + VMUUUU(LTR,K,J,ITR,1) &
                                + VMUUUU(K,LTR,ITR,J,1)
                        ELSE

!                          (1) Symm. transition density matrix

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) + VMUUUU(L,K,J,I,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(LTR,K,JTR,I,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(L,KTR,JTR,I,1)
                        END IF
                     END IF
                  ELSE
                     IF (PSYMM) THEN
                        IF (KSYMM) THEN

!                          (1) Particle symmetrization
!                          (2) Kramers symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) &
                                + VMUUUU(I,J,K,L,1) &
                                + VMUUUU(LTR,KTR,I,J,1) &
                                + VMUUUU(I,J,LTR,KTR,1) &
                                + VMUUUU(K,L,JTR,ITR,1) &
                                + VMUUUU(JTR,ITR,K,L,1) &
                                + VMUUUU(LTR,KTR,JTR,ITR,1) &
                                + VMUUUU(JTR,ITR,LTR,KTR,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(I,JTR,K,LTR,1) &
                                - VMUUUU(L,KTR,I,JTR,1) &
                                - VMUUUU(I,JTR,L,KTR,1) &
                                - VMUUUU(K,LTR,J,ITR,1) &
                                - VMUUUU(J,ITR,K,LTR,1) &
                                + VMUUUU(L,KTR,J,ITR,1) &
                                + VMUUUU(J,ITR,L,KTR,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(I,JTR,KTR,L,1) &
                                - VMUUUU(LTR,K,I,JTR,1) &
                                - VMUUUU(I,JTR,LTR,K,1) &
                                - VMUUUU(KTR,L,J,ITR,1) &
                                - VMUUUU(J,ITR,KTR,L,1) &
                                + VMUUUU(LTR,K,J,ITR,1) &
                                + VMUUUU(J,ITR,LTR,K,1)
                        ELSE

!                          (1) Particle symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) &
                                + VMUUUU(I,J,K,L,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                + VMUUUU(I,JTR,K,LTR,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                + VMUUUU(I,JTR,KTR,L,1)
                        END IF
                     ELSE
                        IF (KSYMM) THEN

!                          (1) Kramers symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = &
                                VMUUUU(K,L,I,J,1) &
                                + VMUUUU(LTR,KTR,I,J,1) &
                                + VMUUUU(K,L,JTR,ITR,1) &
                                + VMUUUU(LTR,KTR,JTR,ITR,1)
                           VDUUUU(KK,LL,IJ,1,2) = &
                                VMUUUU(K,LTR,I,JTR,1) &
                                - VMUUUU(L,KTR,I,JTR,1) &
                                - VMUUUU(K,LTR,J,ITR,1) &
                                + VMUUUU(L,KTR,J,ITR,1)
                           VDUUUU(KK,LL,IJ,1,3) = &
                                VMUUUU(KTR,L,I,JTR,1) &
                                - VMUUUU(LTR,K,I,JTR,1) &
                                - VMUUUU(KTR,L,J,ITR,1) &
                                + VMUUUU(LTR,K,J,ITR,1)
                        ELSE

!                          No symmetrization

                           VDUUUU(KK,LL,IJ,1,1) = VMUUUU(K,L,I,J,1)
                           VDUUUU(KK,LL,IJ,1,2) = VMUUUU(K,LTR,I,JTR,1)
                           VDUUUU(KK,LL,IJ,1,3) = VMUUUU(KTR,L,I,JTR,1)
                        END IF
                     END IF
                  END IF

!                 Complex point groups
!                 --------------------

                  IF (NZ .GE. 2) THEN
                     IF (TSYMM) THEN
                        IF (PSYMM) THEN
                           IF (KSYMM) THEN

!                             (1) Symm. transition density matrix
!                             (2) Do particle symmetrization
!                             (3) Do Kramer's symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   + VMUUUU(I,J,K,L,2) &
                                   - VMUUUU(L,K,J,I,2) &
                                   - VMUUUU(J,I,L,K,2) &
                                   + VMUUUU(LTR,KTR,I,J,2) &
                                   + VMUUUU(I,J,LTR,KTR,2) &
                                   - VMUUUU(KTR,LTR,J,I,2) &
                                   - VMUUUU(J,I,KTR,LTR,2) &
                                   + VMUUUU(K,L,JTR,ITR,2) &
                                   + VMUUUU(JTR,ITR,K,L,2) &
                                   - VMUUUU(L,K,ITR,JTR,2) &
                                   - VMUUUU(ITR,JTR,L,K,2) &
                                   + VMUUUU(LTR,KTR,JTR,ITR,2) &
                                   + VMUUUU(JTR,ITR,LTR,KTR,2) &
                                   - VMUUUU(KTR,LTR,ITR,JTR,2) &
                                   - VMUUUU(ITR,JTR,KTR,LTR,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   + VMUUUU(I,JTR,K,LTR,2) &
                                   - VMUUUU(LTR,K,JTR,I,2) &
                                   - VMUUUU(JTR,I,LTR,K,2) &
                                   - VMUUUU(L,KTR,I,JTR,2) &
                                   - VMUUUU(I,JTR,L,KTR,2) &
                                   + VMUUUU(KTR,L,JTR,I,2) &
                                   + VMUUUU(JTR,I,KTR,L,2) &
                                   - VMUUUU(K,LTR,J,ITR,2) &
                                   - VMUUUU(J,ITR,K,LTR,2) &
                                   + VMUUUU(LTR,K,ITR,J,2) &
                                   + VMUUUU(ITR,J,LTR,K,2) &
                                   + VMUUUU(L,KTR,J,ITR,2) &
                                   + VMUUUU(J,ITR,L,KTR,2) &
                                   - VMUUUU(KTR,L,ITR,J,2) &
                                   - VMUUUU(ITR,J,KTR,L,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   + VMUUUU(I,JTR,KTR,L,2) &
                                   - VMUUUU(L,KTR,JTR,I,2) &
                                   - VMUUUU(JTR,I,L,KTR,2) &
                                   - VMUUUU(LTR,K,I,JTR,2) &
                                   - VMUUUU(I,JTR,LTR,K,2) &
                                   + VMUUUU(K,LTR,JTR,I,2) &
                                   + VMUUUU(JTR,I,K,LTR,2) &
                                   - VMUUUU(KTR,L,J,ITR,2) &
                                   - VMUUUU(J,ITR,KTR,L,2) &
                                   + VMUUUU(L,KTR,ITR,J,2) &
                                   + VMUUUU(ITR,J,L,KTR,2) &
                                   + VMUUUU(LTR,K,J,ITR,2) &
                                   + VMUUUU(J,ITR,LTR,K,2) &
                                   - VMUUUU(K,LTR,ITR,J,2) &
                                   - VMUUUU(ITR,J,K,LTR,2)
                           ELSE

!                            (1) Symm. transition density matrix
!                            (2) Do particle symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   - VMUUUU(L,K,J,I,2) &
                                   + VMUUUU(I,J,K,L,2) &
                                   - VMUUUU(J,I,L,K,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   - VMUUUU(LTR,K,JTR,I,2) &
                                   + VMUUUU(I,JTR,K,LTR,2) &
                                   - VMUUUU(JTR,I,LTR,K,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   - VMUUUU(L,KTR,JTR,I,2) &
                                   + VMUUUU(I,JTR,KTR,L,2) &
                                   - VMUUUU(JTR,I,L,KTR,2)
                           END IF
                        ELSE
                           IF (KSYMM) THEN

!                            (1) Symm. transition density matrix
!                            (2) Kramers symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   - VMUUUU(L,K,J,I,2) &
                                   + VMUUUU(LTR,KTR,I,J,2) &
                                   - VMUUUU(KTR,LTR,J,I,2) &
                                   + VMUUUU(K,L,JTR,ITR,2) &
                                   - VMUUUU(L,K,ITR,JTR,2) &
                                   + VMUUUU(LTR,KTR,JTR,ITR,2) &
                                   - VMUUUU(KTR,LTR,ITR,JTR,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   - VMUUUU(LTR,K,JTR,I,2) &
                                   - VMUUUU(L,KTR,I,JTR,2) &
                                   + VMUUUU(KTR,L,JTR,I,2) &
                                   - VMUUUU(K,LTR,J,ITR,2) &
                                   + VMUUUU(LTR,K,ITR,J,2) &
                                   + VMUUUU(L,KTR,J,ITR,2) &
                                   - VMUUUU(KTR,L,ITR,J,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   - VMUUUU(L,KTR,JTR,I,2) &
                                   - VMUUUU(LTR,K,I,JTR,2) &
                                   + VMUUUU(K,LTR,JTR,I,2) &
                                   - VMUUUU(KTR,L,J,ITR,2) &
                                   + VMUUUU(L,KTR,ITR,J,2) &
                                   + VMUUUU(LTR,K,J,ITR,2) &
                                   - VMUUUU(K,LTR,ITR,J,2)
                           ELSE

!                            (1) Symm. transition density matrix

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   - VMUUUU(L,K,J,I,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   - VMUUUU(LTR,K,JTR,I,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   - VMUUUU(L,KTR,JTR,I,2)
                           END IF
                        END IF
                     ELSE
                        IF (PSYMM) THEN
                           IF (KSYMM) THEN

!                            (1) Particle symmetrization
!                            (2) Kramers symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   + VMUUUU(I,J,K,L,2) &
                                   + VMUUUU(LTR,KTR,I,J,2) &
                                   + VMUUUU(I,J,LTR,KTR,2) &
                                   + VMUUUU(K,L,JTR,ITR,2) &
                                   + VMUUUU(JTR,ITR,K,L,2) &
                                   + VMUUUU(LTR,KTR,JTR,ITR,2) &
                                   + VMUUUU(JTR,ITR,LTR,KTR,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   + VMUUUU(I,JTR,K,LTR,2) &
                                   - VMUUUU(L,KTR,I,JTR,2) &
                                   - VMUUUU(I,JTR,L,KTR,2) &
                                   - VMUUUU(K,LTR,J,ITR,2) &
                                   - VMUUUU(J,ITR,K,LTR,2) &
                                   + VMUUUU(L,KTR,J,ITR,2) &
                                   + VMUUUU(J,ITR,L,KTR,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   + VMUUUU(I,JTR,KTR,L,2) &
                                   - VMUUUU(LTR,K,I,JTR,2) &
                                   - VMUUUU(I,JTR,LTR,K,2) &
                                   - VMUUUU(KTR,L,J,ITR,2) &
                                   - VMUUUU(J,ITR,KTR,L,2) &
                                   + VMUUUU(LTR,K,J,ITR,2) &
                                   + VMUUUU(J,ITR,LTR,K,2)
                           ELSE

!                            (1) Particle symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   + VMUUUU(I,J,K,L,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   + VMUUUU(I,JTR,K,LTR,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   + VMUUUU(I,JTR,KTR,L,2)
                           END IF
                        ELSE
                           IF (KSYMM) THEN

!                            (1) Kramers symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2) &
                                   + VMUUUU(LTR,KTR,I,J,2) &
                                   + VMUUUU(K,L,JTR,ITR,2) &
                                   + VMUUUU(LTR,KTR,JTR,ITR,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2) &
                                   - VMUUUU(L,KTR,I,JTR,2) &
                                   - VMUUUU(K,LTR,J,ITR,2) &
                                   + VMUUUU(L,KTR,J,ITR,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2) &
                                   - VMUUUU(LTR,K,I,JTR,2) &
                                   - VMUUUU(KTR,L,J,ITR,2) &
                                   + VMUUUU(LTR,K,J,ITR,2)
                           ELSE

!                             No symmetrization

                              VDUUUU(KK,LL,IJ,2,1) = &
                                   VMUUUU(K,L,I,J,2)
                              VDUUUU(KK,LL,IJ,2,2) = &
                                   VMUUUU(K,LTR,I,JTR,2)
                              VDUUUU(KK,LL,IJ,2,3) = &
                                   VMUUUU(KTR,L,I,JTR,2)
                           END IF
                        END IF
                     END IF
                  END IF

          if(nz == 4)then  ! goto nz4

          if(tsymm)then
            if(psymm)then
              if(ksymm)then

!               (1) Symm. transition density matrix
!               (2) Particle symmetrization
!               (3) Kramer's symmetrization
!                   Order of the implementation: (3), (1), (2)

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) - vmuuuu(ltr,ktr,j  ,itr,1) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,1) - vmuuuu(k  ,l  ,j  ,itr,1) &
                                     + vmuuuu(l  ,k  ,jtr,i  ,1) - vmuuuu(ktr,ltr,itr,j  ,1) &
                                     + vmuuuu(ktr,ltr,jtr,i  ,1) - vmuuuu(l  ,k  ,itr,j  ,1) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,1) - vmuuuu(j  ,itr,ltr,ktr,1) &
                                     + vmuuuu(i  ,jtr,ltr,ktr,1) - vmuuuu(j  ,itr,k  ,l  ,1) &
                                     + vmuuuu(jtr,i  ,l  ,k  ,1) - vmuuuu(itr,j  ,ktr,ltr,1) &
                                     + vmuuuu(jtr,i  ,ktr,ltr,1) - vmuuuu(itr,j  ,l  ,k  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(ltr,ktr,j  ,itr,2) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,2) - vmuuuu(k  ,l  ,j  ,itr,2) &
                                     - vmuuuu(l  ,k  ,jtr,i  ,2) + vmuuuu(ktr,ltr,itr,j  ,2) &
                                     - vmuuuu(ktr,ltr,jtr,i  ,2) + vmuuuu(l  ,k  ,itr,j  ,2) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,2) - vmuuuu(j  ,itr,ltr,ktr,2) &
                                     + vmuuuu(i  ,jtr,ltr,ktr,2) - vmuuuu(j  ,itr,k  ,l  ,2) &
                                     - vmuuuu(jtr,i  ,l  ,k  ,2) + vmuuuu(itr,j  ,ktr,ltr,2) &
                                     - vmuuuu(jtr,i  ,ktr,ltr,2) + vmuuuu(itr,j  ,l  ,k  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) - vmuuuu(l  ,ktr,jtr,itr,1) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,1) + vmuuuu(k  ,ltr,jtr,itr,1) &
                                     + vmuuuu(ltr,k  ,j  ,i  ,1) - vmuuuu(ktr,l  ,itr,jtr,1) &
                                     - vmuuuu(ktr,l  ,j  ,i  ,1) + vmuuuu(ltr,k  ,itr,jtr,1) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,1) - vmuuuu(jtr,itr,l  ,ktr,1) &
                                     - vmuuuu(i  ,j  ,l  ,ktr,1) + vmuuuu(jtr,itr,k  ,ltr,1) &
                                     + vmuuuu(j  ,i  ,ltr,k  ,1) - vmuuuu(itr,jtr,ktr,l  ,1) &
                                     - vmuuuu(j  ,i  ,ktr,l  ,1) + vmuuuu(itr,jtr,ltr,k  ,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(l  ,ktr,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,2) + vmuuuu(k  ,ltr,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,j  ,i  ,2) + vmuuuu(ktr,l  ,itr,jtr,2) &
                                     + vmuuuu(ktr,l  ,j  ,i  ,2) - vmuuuu(ltr,k  ,itr,jtr,2) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,2) - vmuuuu(jtr,itr,l  ,ktr,2) &
                                     - vmuuuu(i  ,j  ,l  ,ktr,2) + vmuuuu(jtr,itr,k  ,ltr,2) &
                                     - vmuuuu(j  ,i  ,ltr,k  ,2) + vmuuuu(itr,jtr,ktr,l  ,2) &
                                     + vmuuuu(j  ,i  ,ktr,l  ,2) - vmuuuu(itr,jtr,ltr,k  ,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) - vmuuuu(ltr,k  ,jtr,itr,1) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,1) + vmuuuu(ktr,l  ,jtr,itr,1) &
                                     + vmuuuu(l  ,ktr,j  ,i  ,1) - vmuuuu(k  ,ltr,itr,jtr,1) &
                                     - vmuuuu(k  ,ltr,j  ,i  ,1) + vmuuuu(l  ,ktr,itr,jtr,1) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,1) - vmuuuu(jtr,itr,ltr,k  ,1) &
                                     - vmuuuu(i  ,j  ,ltr,k  ,1) + vmuuuu(jtr,itr,ktr,l  ,1) &
                                     + vmuuuu(j  ,i  ,l  ,ktr,1) - vmuuuu(itr,jtr,k  ,ltr,1) &
                                     - vmuuuu(j  ,i  ,k  ,ltr,1) + vmuuuu(itr,jtr,l  ,ktr,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(ltr,k  ,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,2) + vmuuuu(ktr,l  ,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,j  ,i  ,2) + vmuuuu(k  ,ltr,itr,jtr,2) &
                                     + vmuuuu(k  ,ltr,j  ,i  ,2) - vmuuuu(l  ,ktr,itr,jtr,2) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,2) - vmuuuu(jtr,itr,ltr,k  ,2) &
                                     - vmuuuu(i  ,j  ,ltr,k  ,2) + vmuuuu(jtr,itr,ktr,l  ,2) &
                                     - vmuuuu(j  ,i  ,l  ,ktr,2) + vmuuuu(itr,jtr,k  ,ltr,2) &
                                     + vmuuuu(j  ,i  ,k  ,ltr,2) - vmuuuu(itr,jtr,l  ,ktr,2)
              else

!               (1) Symm. transition density matrix
!               (2) Particle symmetrization
!                   Order of the implementation: (1), (2)

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) + vmuuuu(l  ,k  ,jtr,i  ,1) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,1) + vmuuuu(jtr,i  ,l  ,k  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(l  ,k  ,jtr,i  ,2) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,2) - vmuuuu(jtr,i  ,l  ,k  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) + vmuuuu(ltr,k  ,j  ,i  ,1) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,1) + vmuuuu(j  ,i  ,ltr,k  ,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(ltr,k  ,j  ,i  ,2) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,2) - vmuuuu(j  ,i  ,ltr,k  ,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) + vmuuuu(l  ,ktr,j  ,i  ,1) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,1) + vmuuuu(j  ,i  ,l  ,ktr,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(l  ,ktr,j  ,i  ,2) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,2) - vmuuuu(j  ,i  ,l  ,ktr,2)
              endif
            else
              if(ksymm)then

!               (1) Symm. transition density matrix
!               (3) Kramers symmetrization
!                   Order of the implementation: (3), (1)

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) - vmuuuu(ltr,ktr,j  ,itr,1) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,1) - vmuuuu(k  ,l  ,j  ,itr,1) &
                                     + vmuuuu(l  ,k  ,jtr,i  ,1) - vmuuuu(ktr,ltr,itr,j  ,1) &
                                     + vmuuuu(ktr,ltr,jtr,i  ,1) - vmuuuu(l  ,k  ,itr,j  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(ltr,ktr,j  ,itr,2) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,2) - vmuuuu(k  ,l  ,j  ,itr,2) &
                                     - vmuuuu(l  ,k  ,jtr,i  ,2) + vmuuuu(ktr,ltr,itr,j  ,2) &
                                     - vmuuuu(ktr,ltr,jtr,i  ,2) + vmuuuu(l  ,k  ,itr,j  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) - vmuuuu(l  ,ktr,jtr,itr,1) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,1) + vmuuuu(k  ,ltr,jtr,itr,1) &
                                     + vmuuuu(ltr,k  ,j  ,i  ,1) - vmuuuu(ktr,l  ,itr,jtr,1) &
                                     - vmuuuu(ktr,l  ,j  ,i  ,1) + vmuuuu(ltr,k  ,itr,jtr,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(l  ,ktr,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,2) + vmuuuu(k  ,ltr,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,j  ,i  ,2) + vmuuuu(ktr,l  ,itr,jtr,2) &
                                     + vmuuuu(ktr,l  ,j  ,i  ,2) - vmuuuu(ltr,k  ,itr,jtr,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) - vmuuuu(ltr,k  ,jtr,itr,1) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,1) + vmuuuu(ktr,l  ,jtr,itr,1) &
                                     + vmuuuu(l  ,ktr,j  ,i  ,1) - vmuuuu(k  ,ltr,itr,jtr,1) &
                                     - vmuuuu(k  ,ltr,j  ,i  ,1) + vmuuuu(l  ,ktr,itr,jtr,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(ltr,k  ,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,2) + vmuuuu(ktr,l  ,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,j  ,i  ,2) + vmuuuu(k  ,ltr,itr,jtr,2) &
                                     + vmuuuu(k  ,ltr,j  ,i  ,2) - vmuuuu(l  ,ktr,itr,jtr,2)
              else

!               (1) Symm. transition density matrix

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) + vmuuuu(l  ,k  ,jtr,i  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(l  ,k  ,jtr,i  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) + vmuuuu(ltr,k  ,j  ,i  ,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(ltr,k  ,j  ,i  ,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) + vmuuuu(l  ,ktr,j  ,i  ,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(l  ,ktr,j  ,i  ,2)
              endif
            endif
          else
            if(psymm)then
              if(ksymm)then

!               (2) Particle symmetrization
!               (3) Kramers symmetrization
!                   Order of the implementation: (3), (2)

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) - vmuuuu(ltr,ktr,j  ,itr,1) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,1) - vmuuuu(k  ,l  ,j  ,itr,1) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,1) - vmuuuu(j  ,itr,ltr,ktr,1) &
                                     + vmuuuu(i  ,jtr,ltr,ktr,1) - vmuuuu(j  ,itr,k  ,l  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(ltr,ktr,j  ,itr,2) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,2) - vmuuuu(k  ,l  ,j  ,itr,2) &
                                     + vmuuuu(i  ,jtr,k  ,l  ,2) - vmuuuu(j  ,itr,ltr,ktr,2) &
                                     + vmuuuu(i  ,jtr,ltr,ktr,2) - vmuuuu(j  ,itr,k  ,l  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) - vmuuuu(l  ,ktr,jtr,itr,1) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,1) + vmuuuu(k  ,ltr,jtr,itr,1) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,1) - vmuuuu(jtr,itr,l  ,ktr,1) &
                                     - vmuuuu(i  ,j  ,l  ,ktr,1) + vmuuuu(jtr,itr,k  ,ltr,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(l  ,ktr,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,2) + vmuuuu(k  ,ltr,jtr,itr,2) &
                                     + vmuuuu(i  ,j  ,k  ,ltr,2) - vmuuuu(jtr,itr,l  ,ktr,2) &
                                     - vmuuuu(i  ,j  ,l  ,ktr,2) + vmuuuu(jtr,itr,k  ,ltr,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) - vmuuuu(ltr,k  ,jtr,itr,1) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,1) + vmuuuu(ktr,l  ,jtr,itr,1) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,1) - vmuuuu(jtr,itr,ltr,k  ,1) &
                                     - vmuuuu(i  ,j  ,ltr,k  ,1) + vmuuuu(jtr,itr,ktr,l  ,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(ltr,k  ,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,2) + vmuuuu(ktr,l  ,jtr,itr,2) &
                                     + vmuuuu(i  ,j  ,ktr,l  ,2) - vmuuuu(jtr,itr,ltr,k  ,2) &
                                     - vmuuuu(i  ,j  ,ltr,k  ,2) + vmuuuu(jtr,itr,ktr,l  ,2)
              else

!               (2) Particle symmetrization

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) + vmuuuu(i  ,jtr,k  ,l  ,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) + vmuuuu(i  ,jtr,k  ,l  ,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) + vmuuuu(i  ,j  ,k  ,ltr,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) + vmuuuu(i  ,j  ,k  ,ltr,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) + vmuuuu(i  ,j  ,ktr,l  ,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) + vmuuuu(i  ,j  ,ktr,l  ,2)
              endif
            else
              if(ksymm)then

!               (3) Kramers symmetrization

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1) - vmuuuu(ltr,ktr,j  ,itr,1) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,1) - vmuuuu(k  ,l  ,j  ,itr,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2) - vmuuuu(ltr,ktr,j  ,itr,2) &
                                     + vmuuuu(ltr,ktr,i  ,jtr,2) - vmuuuu(k  ,l  ,j  ,itr,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1) - vmuuuu(l  ,ktr,jtr,itr,1) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,1) + vmuuuu(k  ,ltr,jtr,itr,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2) - vmuuuu(l  ,ktr,jtr,itr,2) &
                                     - vmuuuu(l  ,ktr,i  ,j  ,2) + vmuuuu(k  ,ltr,jtr,itr,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1) - vmuuuu(ltr,k  ,jtr,itr,1) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,1) + vmuuuu(ktr,l  ,jtr,itr,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2) - vmuuuu(ltr,k  ,jtr,itr,2) &
                                     - vmuuuu(ltr,k  ,i  ,j  ,2) + vmuuuu(ktr,l  ,jtr,itr,2)
              else

!               No symmetrization

                vduuuu(kk,ll,ij,3,1) = vmuuuu(k  ,l  ,i  ,jtr,1)
                vduuuu(kk,ll,ij,4,1) = vmuuuu(k  ,l  ,i  ,jtr,2)
                vduuuu(kk,ll,ij,3,2) = vmuuuu(k  ,ltr,i  ,j  ,1)
                vduuuu(kk,ll,ij,4,2) = vmuuuu(k  ,ltr,i  ,j  ,2)
                vduuuu(kk,ll,ij,3,3) = vmuuuu(ktr,l  ,i  ,j  ,1)
                vduuuu(kk,ll,ij,4,3) = vmuuuu(ktr,l  ,i  ,j  ,2)
              endif
            endif
          endif

          endif  ! label nz4

! END OF LOOPS FOR UUUU matrix transformation
        enddo
      enddo
    enddo
  enddo

!     remove integrals that should be zero in spinfree calculations but
!     are not due to numerical noise
      if(spinfr.or.levyle)then
        do i = 2,3
          call dzero(VDUUUU(1,1,1,1,i),NASHT*NASHT*NNASHX*NZ)
        end do
      end if

 !     IF (IPRINT .GE. 30) THEN
         CALL HEADER('VUUUU matrix in Dirac (NZ,3) format',-1)
         CALL PRDNZ3(VDUUUU,NASHT,NNASHX,NZ,IPQTOQ(1,0),LUPRI)
 !     END IF

  !------------------------------------------------------------------------------!
  !           Write to the file ouput information for unit testing               !
  !------------------------------------------------------------------------------!
#ifdef CREATE_UNIT_TEST_M2DNZ390
    call add_out_variable(vduuuu,'vduuuu')
    call generate_test()
#endif

  end subroutine m2dnz390



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck dnz32m */
      SUBROUTINE GGUUDNZ32M(VDGGUU,VMGGUU,NORBT,NORB,NASHT,NASH,NNASHX &
        ,NZ,NBSYM,NFSYM,IPQTOQ,IPRINT,idxu2g,maxash)
!***********************************************************************
!
!     Transform a UUUU matrix in Dirac (NZ,3) format to a UUUU matrix
!     in Molfdir format.
!
!     Written by J. Thyssen - Nov 20 2000
!     Last revision :
!
!***********************************************************************

      INTEGER, INTENT(IN) :: NORBT,NASHT,NNASHX,NZ,NBSYM,NFSYM,NORB(2),NASH(2)
      INTEGER, INTENT(IN) :: IPQTOQ(4,0:7),maxash,IDXU2G(MAXASH)
      REAL*8, INTENT(IN) :: VDGGUU(NORBT,NORBT,NNASHX,NZ,*)
      PARAMETER (THROUT = 1.0D-15)
      REAL*8, INTENT(OUT) :: VMGGUU(2*NORBT,2*NORBT,2*NASHT,2*NASHT,*)

      INTEGER :: P,J,Q,PP,QQ,JJ,IJ,PTR,JTR,QTR,KR(-300:300),IFSYM
      INTEGER :: KRG(-NORBT:NORBT)
      DIMENSION DINT(8,2)

!
      NCLASS = NZ * NZ * NBSYM / NFSYM
      CALL DZERO(VMGGUU, (2*NASHT)**2*(2*NORBT)**2 * MIN(NZ,2))
!
      iprint_local = iprint
#ifdef DMRG_DEBUG
!     debug
      iprint_local = 30
#endif
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO

      PP = 0
      QQ = 0 
      DO IFSYM = 1, NFSYM 
        DO I = 1, NORB(IFSYM)
          PP = PP + 1
          QQ = QQ + 1
          KRG(PP) = QQ
          KRG(-PP) = QQ + NORB(IFSYM)
        ENDDO
        QQ = QQ + NORB(IFSYM)
      ENDDO


!
! MvH      IF (iprint_local .GE. 30) THEN
!         WRITE(LUPRI,'(/1X,A)') '(DNZ32M) Input VUUUU matrix in Dirac (NZ,3) format'
!         CALL PRDNZ3(VDGGUU,NORBT,NNASHX,NZ,IPQTOQ(1,0),LUPRI)
!      END IF
!
      DO II = 1, NASHT
         DO JJ = 1, II
            IJ = (II * (II - 1))/2 + JJ
            DO QQ = 1, NORBT
               DO PP = 1, NORBT
                  I = KR(II)
                  J = KR(JJ)
                  P = KRG(PP)
                  Q = KRG(QQ)
                  ITR = KR(-II)
                  JTR = KR(-JJ)
                  PTR = KRG(-PP)
                  QTR = KRG(-QQ)

                 ! write(6,*) 'I,J,P,Q,ITR,JTR,PTR,QTR ',I,J,P,Q,ITR,JTR,PTR,QTR

!
!                 Canonical integrals (I>=J)
!                 (pl|ij), (pL|iJ), (Pl|iJ)
!
                  VMGGUU(P,Q,I,J,1)     = VDGGUU(PP,QQ,IJ,1,1)
                  VMGGUU(P,QTR,I,JTR,1) = VDGGUU(PP,QQ,IJ,1,2)
                  VMGGUU(PTR,Q,I,JTR,1) = VDGGUU(PP,QQ,IJ,1,3)
!
!                 The rest:
!
!                 (Pl|Ij) = (pL|iJ)^*

                  VMGGUU(PTR,Q,ITR,J,1) =  VDGGUU(PP,QQ,IJ,1,2)
!                 (PL|ij) = (LP|IJ)^* = (pl|ji)
                  VMGGUU(PTR,QTR,I,J,1) = VDGGUU(QQ,PP,IJ,1,1)
!                 (KL|IJ) = (lk|ji) = (kl|ij)^*
                  VMGGUU(PTR,QTR,ITR,JTR,1) = VDGGUU(PP,QQ,IJ,1,1)
!                 (kl|IJ) = (kl|ji) = (lk|ij)
                  VMGGUU(P,Q,ITR,JTR,1) = VDGGUU(QQ,PP,IJ,1,1)
!                 (kL|Ij) = - (kL|Ji) = - (Lk|iJ)^*
                  VMGGUU(P,QTR,ITR,J,1) = - VDGGUU(QQ,PP,IJ,1,3)
!
                  IF (II .NE. JJ) THEN
!
                     VMGGUU(Q,P,J,I,1)     = VDGGUU(PP,QQ,IJ,1,1)
                     VMGGUU(QTR,P,JTR,I,1) = VDGGUU(PP,QQ,IJ,1,2)
                     VMGGUU(Q,PTR,JTR,I,1) = VDGGUU(PP,QQ,IJ,1,3)
!
                     VMGGUU(Q,PTR,J,ITR,1) = - VDGGUU(QQ,PP,IJ,1,2)
                     VMGGUU(QTR,PTR,J,I,1) =   VDGGUU(QQ,PP,IJ,1,1)
                     VMGGUU(QTR,PTR,JTR,ITR,1) = VDGGUU(PP,QQ,IJ,1,1)
                     VMGGUU(Q,P,JTR,ITR,1) =   VDGGUU(QQ,PP,IJ,1,1)
                     VMGGUU(QTR,P,J,ITR,1) = - VDGGUU(QQ,PP,IJ,1,3)
                  END IF
!
!
                  IF (NZ .GE. 2) THEN
!
!                    Canonical integrals (I>=J)
!                    (kl|ij), (kL|iJ), (Kl|iJ)
!
                     VMGGUU(P,Q,I,J,2)     = VDGGUU(PP,QQ,IJ,2,1)
                     VMGGUU(P,QTR,I,JTR,2) = VDGGUU(PP,QQ,IJ,2,2)
                     VMGGUU(PTR,Q,I,JTR,2) = VDGGUU(PP,QQ,IJ,2,3)
!
!                    The rest:
!
!                    (Kl|Ij) = - (Kl|Ji) = - (lK|iJ)^*
                     VMGGUU(PTR,Q,ITR,J,2) = VDGGUU(QQ,PP,IJ,2,2)
!                    (KL|ij) = (lk|ij)
                     VMGGUU(PTR,QTR,I,J,2) = VDGGUU(QQ,PP,IJ,2,1)
!                    (KL|IJ) = (lk|ji) = (kl|ij)^*
                     VMGGUU(PTR,QTR,ITR,JTR,2) = - VDGGUU(PP,QQ,IJ,2,1)
!                    (kl|IJ) = (kl|ji) = (lk|ij)^*
                     VMGGUU(P,Q,ITR,JTR,2) = - VDGGUU(QQ,PP,IJ,2,1)
!                    (kL|Ij) = - (kL|Ji) = - (Lk|iJ)^*
                     VMGGUU(P,QTR,ITR,J,2) = VDGGUU(QQ,PP,IJ,2,3)
!
                     IF (II .NE. JJ) THEN
!
!                       (lk|ji) = (kl|ij)^*
                        VMGGUU(Q,P,J,I,2)     = - VDGGUU(PP,QQ,IJ,2,1)
!                       (Lk|Ji) = (kL|iJ)^*
                        VMGGUU(QTR,P,JTR,I,2) = - VDGGUU(PP,QQ,IJ,2,2)
!                       (lK|Ji) = (Kl|iJ)^*
                        VMGGUU(Q,PTR,JTR,I,2) = - VDGGUU(PP,QQ,IJ,2,3)
!
!                       (lK|jI) = (kL|iJ)
                        VMGGUU(Q,PTR,J,ITR,2) =   VDGGUU(PP,QQ,IJ,2,2)
!                       (LK|ji) = (kl|ji) = (lk|ij)^*
                        VMGGUU(QTR,PTR,J,I,2) = - VDGGUU(QQ,PP,IJ,2,1)
!                       (LK|JI) = (kl|ij)
                        VMGGUU(QTR,PTR,JTR,ITR,2) = VDGGUU(PP,QQ,IJ,2,1)
!                       (lk|JI) = (lk|ij)
                        VMGGUU(Q,P,JTR,ITR,2) = VDGGUU(QQ,PP,IJ,2,1)
!                       (Lk|jI) = (Kl|iJ)
                        VMGGUU(QTR,P,J,ITR,2) = VDGGUU(PP,QQ,IJ,2,3)
                     END IF
!
                  END IF
!
                  IF (NZ .EQ. 4) THEN
!
!                    (Kl|ij) (canonical)
!
                     VMGGUU(PTR,Q,I,J,1) = VDGGUU(PP,QQ,IJ,3,3)
                     VMGGUU(PTR,Q,I,J,2) = VDGGUU(PP,QQ,IJ,4,3)
                     IF (II .NE. JJ) THEN
!                       (lK|ji) = (Kl|ij)^*
                        VMGGUU(Q,PTR,J,I,1) =   VDGGUU(PP,QQ,IJ,3,3)
                        VMGGUU(Q,PTR,J,I,2) = - VDGGUU(PP,QQ,IJ,4,3)
                     END IF
!
!                    (kL|ij) (canonical)
!
                     VMGGUU(P,QTR,I,J,1) = VDGGUU(PP,QQ,IJ,3,2)
                     VMGGUU(P,QTR,I,J,2) = VDGGUU(PP,QQ,IJ,4,2)
                     IF (II .NE. JJ) THEN
!                       (Lk|ji) = (kL|ij)^*
                        VMGGUU(QTR,P,J,I,1) =   VDGGUU(PP,QQ,IJ,3,2)
                        VMGGUU(QTR,P,J,I,2) = - VDGGUU(PP,QQ,IJ,4,2)
                     END IF
!
!                    (kl|Ij) = - (kl|Ji) = - (lk|iJ)^*
!
                     VMGGUU(P,Q,ITR,J,1) = - VDGGUU(QQ,PP,IJ,3,1)
                     VMGGUU(P,Q,ITR,J,2) =   VDGGUU(QQ,PP,IJ,4,1)
                     IF (II .NE. JJ) THEN
!                       (lk|jI) = - (lk|iJ)
                        VMGGUU(Q,P,J,ITR,1) = - VDGGUU(QQ,PP,IJ,3,1)
                        VMGGUU(Q,P,J,ITR,2) = - VDGGUU(QQ,PP,IJ,4,1)
                     END IF
!
!                    (kl|iJ) (canonical)
!
                     VMGGUU(P,Q,I,JTR,1) = VDGGUU(PP,QQ,IJ,3,1)
                     VMGGUU(P,Q,I,JTR,2) = VDGGUU(PP,QQ,IJ,4,1)
                     IF (II .NE. JJ) THEN
!                       (lk|Ji) = (kl|iJ)^*
                        VMGGUU(Q,P,JTR,I,1) =   VDGGUU(PP,QQ,IJ,3,1)
                        VMGGUU(Q,P,JTR,I,2) = - VDGGUU(PP,QQ,IJ,4,1)
                     END IF
!
!
!                    (3 bars)
!
!                    (kL|IJ) = -(Kl|ij)^*
!
                     VMGGUU(P,QTR,ITR,JTR,1) = - VDGGUU(PP,QQ,IJ,3,3)
                     VMGGUU(P,QTR,ITR,JTR,2) =   VDGGUU(PP,QQ,IJ,4,3)
                     IF (II .NE. JJ) THEN
!                       (Lk|JI) = (kL|IJ)^* = -(Kl|ij)
                        VMGGUU(QTR,P,JTR,ITR,1) = - VDGGUU(PP,QQ,IJ,3,3)
                        VMGGUU(QTR,P,JTR,ITR,2) = - VDGGUU(PP,QQ,IJ,4,3)
                     END IF
!
!                    (Kl|IJ) = -(kL|ij)^*
!
                     VMGGUU(PTR,Q,ITR,JTR,1) = - VDGGUU(PP,QQ,IJ,3,2)
                     VMGGUU(PTR,Q,ITR,JTR,2) =   VDGGUU(PP,QQ,IJ,4,2)
                     IF (II .NE. JJ) THEN
!                       (lK|JI) = (Kl|IJ)^* = -(kL|ij)
                        VMGGUU(Q,PTR,JTR,ITR,1) = - VDGGUU(PP,QQ,IJ,3,2)
                        VMGGUU(Q,PTR,JTR,ITR,2) = - VDGGUU(PP,QQ,IJ,4,2)
                     END IF
!
!                    (KL|Ij) = -(kl|iJ)^*
!
                        VMGGUU(PTR,QTR,ITR,J,1) = - VDGGUU(PP,QQ,IJ,3,1)
                        VMGGUU(PTR,QTR,ITR,J,2) =   VDGGUU(PP,QQ,IJ,4,1)
                     IF (II .NE. JJ) THEN
!                       (KL|jI) = -(kl|Ji)^* = -(lk|iJ)
                        VMGGUU(PTR,QTR,J,ITR,1) = - VDGGUU(QQ,PP,IJ,3,1)
                        VMGGUU(PTR,QTR,J,ITR,2) = - VDGGUU(QQ,PP,IJ,4,1)
                     END IF
! MvH inserted my code here
                        VMGGUU(PTR,QTR,I,JTR,1) = VDGGUU(QQ,PP,IJ,3,1)
                        VMGGUU(PTR,QTR,I,JTR,2) = VDGGUU(QQ,PP,IJ,4,1)
                       IF (II .NE. JJ) THEN
!                       (KL|Ji) = -(kl|jI)^* = -(jI|kl)^*
                        VMGGUU(PTR,QTR,JTR,I,1) =  VDGGUU(QQ,PP,IJ,3,1)
                        VMGGUU(PTR,QTR,JTR,I,2) = -VDGGUU(QQ,PP,IJ,4,1)
                       END IF




! MVH original code below

!                     IF (KK .GE. LL) THEN
!                        KL = (KK*(KK-1))/2 + LL
!!
!!                    (KL|iJ) = -(kl|Ij)^* = -(Ij|kl)^*
!!
!                        VMUUUU(KTR,LTR,I,JTR,1) = - VDUUUU(II,JJ,KL,3,3)
!                        VMUUUU(KTR,LTR,I,JTR,2) =   VDUUUU(II,JJ,KL,4,3)
!                       IF (II .NE. JJ) THEN
!!                       (KL|Ji) = -(kl|jI)^* = -(jI|kl)^*
!                        VMUUUU(KTR,LTR,JTR,I,1) = - VDUUUU(JJ,II,KL,3,2)
!                        VMUUUU(KTR,LTR,JTR,I,2) =   VDUUUU(JJ,II,KL,4,2)
!                       END IF
!                     ELSE ! KK .lt. LL
!                        LK = (LL*(LL-1))/2 + KK
!!
!!                    (KL|iJ) = -(kl|Ij)^* = -(Ij|kl)^* = -(jI|lk)
!!
!                        VMUUUU(KTR,LTR,I,JTR,1) = - VDUUUU(JJ,II,LK,3,2)
!                        VMUUUU(KTR,LTR,I,JTR,2) = - VDUUUU(JJ,II,LK,4,2)
!                       IF (II .NE. JJ) THEN
!!                       (KL|Ji) = -(kl|jI)^* = -(jI|kl)^* = -(Ij|lk)
!                        VMUUUU(KTR,LTR,JTR,I,1) = - VDUUUU(II,JJ,LK,3,3)
!                        VMUUUU(KTR,LTR,JTR,I,2) = - VDUUUU(II,JJ,LK,4,3)
!                       END IF
!                     END IF
!
                  END IF
               END DO
            END DO
         END DO
      END DO
!


!
!!9001 FORMAT(1P,1X,4I3,4X,2(D14.7,1X))
! 9002 FORMAT(6I3,5X,4F20.10)
!! MvH      IF (iprint_local .GE. 30) THEN
!         WRITE(LUPRI,'(A,1P,D10.2)')  ' (GGUUDNZ32M) Input VGGUU matrix in Molfdir format,'&
!         //' print threshold =',THROUT
!!#ifdef old_code
!! hjaaj mar 2018: new code below groups the 2-electron integrals,
!! making it easier to spot some problems (particle symmetry, hermiticity)
!! but not Kramers conjugation problems - that requires 2**4=16 times lot more
!! lines for each quadruplet, but with "DO I = 1,NASHT" instead.
!         DO IZ = 1, NZ
!            DO I = 1,NNASHX
!              DO Q = 1,NORBT
!                DO P = 1,NORBT
!                  WRITE(LUPRI,9001) &
!                         P,Q,I,IZ,(VDGGUU(P,Q,I,IZ,I3),I3=1,3)
!                END DO
!              END DO
!            END DO
!        END DO
!         
!



!!
!!9001 FORMAT(1P,1X,4I3,4X,2(D14.7,1X))
! 9001 FORMAT(4I5,5X,4F20.10)
!! MvH      IF (iprint_local .GE. 30) THEN
!         WRITE(LUPRI,'(A,1P,D10.2)')  ' (GGUUDNZ32M) Output VGGUU matrix in Molfdir format,'&
!         //' print threshold =',THROUT
!!#ifdef old_code
!! hjaaj mar 2018: new code below groups the 2-electron integrals,
!! making it easier to spot some problems (particle symmetry, hermiticity)
!! but not Kramers conjugation problems - that requires 2**4=16 times lot more
!! lines for each quadruplet, but with "DO I = 1,NASHT" instead.
!         DO J = 1,2*NASHT
!            DO I = 1,2*NASHT
!               DO K = 1,NASHT
!                  DO L = 1,2*NORBT
!                    P = IDXU2G(L)
!                    Q = IDXU2G(K)
!                         WRITE(LUPRI,*) &
!                         L,Q,I,J,(VMGGUU(L,Q,I,J,IZ),IZ=1,MIN(NZ,2))
!                  END DO
!               END DO
!            END DO
!         END DO
!      ELSE IF (iprint_local .GE. 10) THEN
!#endif
!         DO I = 1, 2*NASHT
!            DO J = 1, I
!               DO L = 1, I
!                  DO K = 1, L
!                  WRITE(LUPRI,'(A)')
!     & 'New quadruplet re(ij|kl), re(kl|ij), im(ij|kl), im(kl|ij)'
!                  WRITE(LUPRI,9001)
!     & I,J,K,L,(VMUUUU(I,J,K,L,IZ),VMUUUU(K,L,I,J,IZ),IZ=1,NZ_in_CI)
!                  WRITE(LUPRI,9001)
!     & J,I,L,K,(VMUUUU(J,I,L,K,IZ),VMUUUU(L,K,J,I,IZ),IZ=1,NZ_in_CI)
!                  WRITE(LUPRI,9001)
!     & I,J,L,K,(VMUUUU(I,J,L,K,IZ),VMUUUU(L,K,I,J,IZ),IZ=1,NZ_in_CI)
!                  WRITE(LUPRI,9001)
!     & J,I,K,L,(VMUUUU(J,I,K,L,IZ),VMUUUU(K,L,J,I,IZ),IZ=1,NZ_in_CI)
!                  END DO
!               END DO
!            END DO
!         END DO
!      END IF
!
      RETURN
      END SUBROUTINE GGUUDNZ32M




end module mcscf_routines
