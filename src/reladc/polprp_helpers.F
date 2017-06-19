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
      SUBROUTINE EARR_VVOO (EPS,T2)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     This routine creates an energy denominator array in the
!     VVOO bra-ket storage <|<|,|>|> in order to be sorted with SRT16
!     Normally the VVOO integrals are tridiagonally stored. In this case
!     the full ranges are covered because SRT 16 down not handle tridiagonal
!     storage cases!
!     The energies are all real and for the complex case we set the imaginary
!     part to zero!
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
!     call PST('Creating ABIJ energy denominators+')
!     write(*,*) 'jvovo:',jvovo(nrep+1)
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      if(multb(irp,jrp,1).ne.ijrp) stop 'earr error'
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = - EPS(JJ)
         DO I = 1, NO(IRP)
            II = IOFF + I
            FAC2 = - EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,IJRP+NREP,2)
      if(multb(arp,brp,1).ne.ijrp) stop 'earr error'
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 + EPS(BB)
               DO A = 1, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 + EPS(AA)
                  ABIJ = ABIJ + 1
                  IF(CARITH) THEN
                    T2(2*ABIJ - 1) = 1.0/FAC  !real part
                    T2(2*ABIJ) = 0.0  ! imaginary part
                  ELSE
                    T2(ABIJ) = 1.0/FAC
                  ENDIF
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      IF(ABIJ.ne.JVOVO(nrep+1)) THEN
        WRITE(*,*) 'ABIJ/JVOVO:',ABIJ,JVOVO(nrep+1)
        STOP 'COUNTING ERROR IN EARR_VVOO!'
      ENDIF
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMVOVO(EPS,T)
 
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
! Divides the <a|i><b|j> integral stream by the corresponding
! denominators e_a - e_i + e_b - e_j.
! The stream has JVOVO organization.
!
!---------------Calling variables--------------------------------------
!
      REAL*8,dimension(:)              :: EPS
      REAL*8,dimension(:)              :: T
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/symm.inc"
#include "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8                           :: EA,EI,EB,EJ
!
!---------------Executable code--------------------------------------
!

!     write(*,*) '+++++++++++ entering DENOMVOVO +++++++++++'
      AIBJ = 0

      DO BJREP = 1, NREP   !bosonic irep
        DO 10 JREP = 1, NREP   !fermionic irep
          JOFF = IO(JREP)
          BREP = MULTB(JREP,BJREP+NREP,2)
          BOFF = IO(NREP+1) + IV(BREP)
          DO J = 1, NO(JREP)
            JJ = JOFF + J
            EJ = eps(JJ)
            DO B = 1, NV(BREP)
               BB = BOFF + B
               EB = eps(BB)
               DO 20 IREP = 1, NREP
                 IOFF = IO(IREP)
                 AREP = MULTB(IREP,BJREP+NREP,2)
                 AOFF = IO(NREP+1) + IV(AREP)
                 DO I = 1, NO(IREP)
                   II = IOFF + I
                   EI = eps(II)
                   DO A = 1, NV(AREP)
                     AA = AOFF + A
                     EA = eps(AA)

                     AIBJ = AIBJ + 1
                     IF (CARITH) THEN
                        stop 'complex not yet available'
                     ELSE
                        T(AIBJ) = T(AIBJ)/(EA+EB-EI-EJ)
                     ENDIF

                   ENDDO !A
                 ENDDO !I
 20           ENDDO  !IREP
            ENDDO  !B
          ENDDO  !J
 10     ENDDO !JREP
      ENDDO  !BJREP

      IF((AIBJ).ne.JVOVO(nrep+1)) Then
        write(*,*) 'Consistency error in DENOMVOVO.'
        write(*,*) 'Counter:',AIBJ
        write(*,*) 'Required value:',JVOVO(nrep+1)
        stop
      ENDIF

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMVOTS(eps,aak)
 
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Divides the totally symmetric A_ak array by eps_a - eps_k
!
!---------------Calling variables--------------------------------------
!
      REAL*8,dimension(:)              :: EPS
      REAL*8,dimension(:)              :: AAK
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8                   :: epsa,epsk,ediff
C
C---------------Executable code--------------------------------------
C
      AK = 0
      DO KREP = 1, NREP
        AREP = KREP    ! array is totally symmetric
        KOFF = IO(KREP)
        AOFF = IO(NREP+1) + IV(AREP)
        DO K = 1,NO(KREP)
          DO A = 1,NV(AREP)
            epsa = EPS(AOFF + A)
            epsk = EPS(KOFF + K)
            ediff = epsa - epsk
            AK = AK+1
            AAK(AK) = AAK(AK)/ediff
!           write(*,*) 'DENOMVOTS: ak,aak(ak):',ak,aak(ak)
          ENDDO
        ENDDO
      ENDDO
      IF(AK.ne.MVO(1)) STOP 'Inconsistency in DENOMVOTS!'

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE rcaxpy(n,carith,phase,as,ad)
      Implicit none
!
!---------------- description ----------------
!
! this routine xaxpys from as to ad depending on phase
! ad is always complex, as is either real or complex (following carith)
!
!---------------- formal parameters ----------------
!
      integer                      :: n
      logical                      :: carith
      complex*16                   :: phase
      real*8,dimension(:)          :: as
      complex*16,dimension(:)      :: ad
!
!---------------- local variables ----------------
!
      integer                      :: rcw
      integer                      :: i
!
!---------------- executable code ----------------
!
      rcw = 1
      if(carith) rcw = 2

!     write(*,*) ' *** in rcaxpy ***',carith,rcw

      if(.not.carith) then       ! ** as = real case + phase **
!       write(*,*) '                  * * * rcaxpy phase:',phase
        do i=1,n
          ad(i) = ad(i) + phase*dcmplx(as(i),0.0d0)
        enddo
      else                       ! ** as = complex case, no phase **
        do i=1,n
          ad(i) = ad(i) + dcmplx(as(2*i-1),as(2*i))
        enddo
      endif

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RCA_NEGATE(n,arr)
!
!  this routine changes sign of a real or a complex array
!  depending on the CARITH flag. In case of complex:
!  z -->  -z  not  z --> z*  !
!
      implicit none
      integer                      :: n
      real*8, dimension(:)         :: arr
!
! it is assumed that n counts the elements and not the array length!
!

#include  "../relccsd/complex.inc"

      integer                      :: i

      DO I=1,N
        IF(CARITH) THEN
          arr(2*i-1) = -1.0D0*arr(2*i-1)
          arr(2*i)   = -1.0D0*arr(2*i)
        ELSE
          arr(i) = -1.0D0*arr(i)
        ENDIF
      ENDDO
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE print_input_info(iw,nrep)
      use polprp_cfg     !these data are visible on master node only!

      implicit none

      integer               :: iw         ! I/O file handle for proper screen output
      integer               :: nrep       ! number of total symmetries in the calculation
!
      integer               :: i,n,icount

      write(iw,*) 'POLPRP control variables:'
      write(iw,*) '-------------------------'
      write(iw,*)
      write(iw,*) 'we do extended calculation: ',polprp_doextended
      write(iw,*) 'we calculate transition moments: ',polprp_dotrmo
      write(iw,'(A,F10.5)') 
     &   'threshold for elements in ADC matrix: ',polprp_writethr
      write(iw,*) 'we perform a Davidson diagonalization:',polprp_dodiag
      write(iw,*) 'The way final state symmetries are chosen is:'
!
! determine symmetries to be calculated
!
      n = 0
      do i = 1,nrep
        n = n + polprp_statesym(i)
      enddo
      if(n.eq.0) then
        write(iw,*) 'Program checks symmetries automatically'
        do i=1,nrep
          polprp_statesym(i) = i
        enddo
      else
        write(iw,*) 'User has selected the following final'
        write(iw,*) 'state symmetries:'
        icount=1
        do i=1,nrep
          if(polprp_statesym(i).ne.0) then
            write(iw,*) '   No: ',icount,
     &                  '  symmetry: ',polprp_statesym(i)
            icount = icount + 1
            if(polprp_statesym(i).gt.nrep) then
              write(iw,*) 'This state symmetry is higher than'
              write(iw,*) 'the number of available symmetries.'
              write(iw,*) 'Setting it to the maximum value.'
              polprp_statesym(i) = nrep
            endif
          endif
        enddo
      endif
      write(iw,*)
!
! report diagonalization parameters
!
      write(iw,*) 'Davidson control variables:'
      write(iw,*) '---------------------------'
      write(iw,*)

      if(polprp_dodiag) then
        WRITE(IW,*) 'Number of eigenstates to be calculated (roots)',
     &     polprp_davroots
        WRITE(IW,*) 'Size of Krylov space: ',polprp_davmaxsp
        WRITE(IW,*) 'Number of Davidson macro iterations: ',
     &     polprp_davmaxit
        WRITE(IW,'(A,E20.4)') ' Requested energy conv.: ',polprp_davconv
        IF(polprp_davooc.eqv..false.) THEN
           WRITE(IW,*) 'Trying full incore storage of Davidson vectors.'
           WRITE(IW,*) 'Will default to out-of-core if resources are'
           WRITE(IW,*) 'not available.'
        ENDIF
      endif

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     SUBROUTINE PRECONTRACT(iounit,iw,eps)
!
!
!  this routine is called by master only. If we have a parallel run
!  then the nodes provide integrals.
!  The routine provides <ab||cd> * <cd||ji>/(ec + ed - ej - ei)
!  contracted over cd.
!

!if defined (VAR_MPI)
!     use interface_to_mpi
!endif

!     implicit none

!     integer                  :: iounit         ! I/O file handle for writing to file
!     integer                  :: iw             ! file handle for writing to screen
!     real*8, dimension(:)     :: eps
!
! local variables
!
!     real*8, allocatable, dimension(:)        :: vvvvz
!     real*8, allocatable, dimension(:)        :: vvoo1
!     real*8, allocatable, dimension(:)        :: vvoop
!     integer                                  :: n1,ixx
!     integer                                  :: ialloc
!     integer                                  :: signal
!     integer*8                                :: n8,i8
!     COMPLEX*16                               :: A0,A1
!     character*9                              :: preconfn='VVVVXVVOO'
!     logical                                  :: isthere


!include  "../relccsd/symm.inc"
!include  "../relccsd/complex.inc"
!include  "../relccsd/ccpar.inc"
!
!   function interfaces
!
!     interface

!       SUBROUTINE denomvvoo_pp(ra1,ra2)
!         real*8, dimension(:)         :: ra1,ra2
!       END SUBROUTINE

!       SUBROUTINE mpi_vvvv_complete(ra1,i1)
!         real*8, dimension(:)            :: ra1
!         integer*8                       :: i1
!       END SUBROUTINE mpi_vvvv_complete

!     end interface
!
! execution
!
!     if(carith) then
!       write(iw,*) ' **) Performing (complex) precontraction'
!       if(rcw.ne.2) call quit ('intenal error in PRECONTRACT')
!     else
!       write(iw,*) ' **) Performing (real) precontraction'
!       if(rcw.ne.1) call quit ('intenal error in PRECONTRACT')
!     endif
!     A0=(0.0D00,0.0D00)
!     A1=(1.0D00,0.0D00)


!if defined (VAR_MPI)
!
!  send the slaves a start signal for int supply only
!
!     signal = 2  !start code.
!     write(iw,*) 'Master indicates PRECONTRACTION (signal=2)'
!     call interface_mpi_bcast(signal,1,MASTER,
!    &                             global_communicator)
!endif

! allocate buffers and read in integrals

!     n1 = ivvoott(nrep+1)
!     allocate(vvoo1(n1*rcw))
!     allocate(vvoop(n1*rcw))

!     n8 = ivvvvtt(nrep+1)
!     allocate(vvvvz(n8*rcw),stat=ialloc)
!     if(ialloc.ne.0) then
!       STOP 'VVVV buffer could not be allocated.'
!     else
!       write(*,*) 'Master allocated',n8,' words for VVVV buffer.'
!     endif

!if defined (VAR_MPI)
!     do i8=1,n8*rcw
!       vvvvz(i8) = -4096.0d0
!     enddo
!     call rdvvvv(vvvvz)  ! master reads its own chunks
!     write(*,*) 'masters chunk reading ... complete.'
!     call mpi_vvvv_complete(vvvvz,n8*rcw)  !.... and completes the stream.
!     write(*,*) 'VVVV check in PRECONTRACT...'
!     do i8=1,n8
!       if(vvvvz(i8).eq. -4096.0d0)
!    &    call quit('Internal error: gap in the VVVV stream!')
!     enddo
!     write(*,*) '.... success. No leaks.'
!else
!     call rdvvvv(vvvvz)  ! master reads complete stream in serial mode
!endif
!     call getvvoo(vvoo1)
!     call denomvvoo_pp(eps,vvoo1)
!
!  do the contraction over all symmetries
!
!     CALL CNTRCT('N','N',NVVT,NOOT,NVVT,A1,vvvvz,vvoo1,A0,vvoop,NREP)
!
!  check if file already present. If yes, delete it.
!
!     Inquire(file=preconfn,exist=isthere)
!     If(isthere) then
!        write(*,*) 'Found a precontraction file. Deleting...'
!        open(unit=iounit, file=preconfn, status='old')
!        close(iounit, status='delete')
!     Endif
!
!  open corresponding file and write out contracted integrals.
!
!     open(unit=iounit,file=preconfn, access='SEQUENTIAL',
!    &     form='UNFORMATTED',status='NEW')
!     write(iounit) (vvoop(ixx), ixx=1,n1*rcw)
!     close(iounit)
!
!  free buffers immediately 
!
!     deallocate(vvoop)
!     deallocate(vvoo1)
!     deallocate(vvvvz)

!     write(*,*) 'PRECONTRACTION finished. Buffers freed'
!     write(*,*)
!     write(*,*)

!     return
!     end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRECONTRACT_NEW(iounit,iw,eps)
!
!
!  this routine is called by master only. If we have a parallel run
!  then the nodes provide integrals.
!  The routine calculates <ab||cd> * <cd||ji>/(ec + ed - ej - ei)
!  contracted over cd.
!

#if defined (VAR_MPI)
      use interface_to_mpi
#endif

      implicit none

      integer                  :: iounit         ! I/O file handle for writing to file
      integer                  :: iw             ! file handle for writing to screen
      real*8, dimension(:)     :: eps
!
! local variables
!
      real*8, allocatable, dimension(:)        :: vvvvz
      real*8, allocatable, dimension(:)        :: vvoo1
      real*8, allocatable, dimension(:)        :: vvoop
      integer                                  :: n1,ixx
      integer                                  :: ialloc
      integer                                  :: signal
      integer*8                                :: n8,i8
      COMPLEX*16                               :: A0,A1
      character*9                              :: preconfn='VVVVXVVOO'
      logical                                  :: isthere


#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/ccpar.inc"
!
!   function interfaces
!
      interface

        SUBROUTINE denomvvoo_pp(ra1,ra2)
          real*8, dimension(:)         :: ra1,ra2
        END SUBROUTINE

        SUBROUTINE mpi_master_vvvv_complete(ra1,i1)
          real*8, dimension(:)            :: ra1
          integer*8                       :: i1
        END SUBROUTINE mpi_master_vvvv_complete

      end interface
!
! execution
!
      if(carith.and.rcw.eq.2) then
        write(iw,*) ' **) Performing (complex) precontraction'
      else if(.not.carith.and.rcw.eq.1) then
        write(iw,*) ' **) Performing (real) precontraction'
      else
        write(iw,*) ' Control error in precontract_new!'
        call quit('error')
      endif
      A0=(0.0D00,0.0D00)
      A1=(1.0D00,0.0D00)

! allocate buffers and read in integrals

      n1 = ivvoott(nrep+1)
      allocate(vvoo1(n1*rcw))
      allocate(vvoop(n1*rcw))

      n8 = ivvvvtt(nrep+1)
      allocate(vvvvz(n8*rcw),stat=ialloc)
      if(ialloc.ne.0) then
        STOP 'VVVV buffer could not be allocated.'
      else
        write(*,*) 'Master allocated',n8,' words for VVVV buffer.'
      endif
!____________________________________________________________________
!|
!|  fetching integrals in serial/parallel mode
!|
#if defined (VAR_MPI)
      do i8=1,n8*rcw
        vvvvz(i8) = -4096.0d0
      enddo
      call rdvvvv(vvvvz)  ! master reads its own chunks
      call mpi_master_vvvv_complete(vvvvz,n8*rcw)  !fetch integrals from slaves
      do i8=1,n8*rcw
        if(vvvvz(i8).eq. -4096.0d0)
     &    call quit('Internal error: gap in the VVVV stream!')
      enddo
      write(*,*) 'Master: <VV||VV> stream ready.'
#else
      call rdvvvv(vvvvz)  ! master reads complete stream in serial mode
#endif
!     write(*,*) '################# int list ##################'
!     do i8=1,n8*rcw
!        write(*,*) ' ### vvvv integral',i8,vvvvz(i8)
!     enddo
!     write(*,*) '################# int list ##################'
!     write(*,*) 'Master VVVV chksum:',dot_product(vvvvz,vvvvz)
!|
!|  
!|___________________________________________________________________

      call getvvoo(vvoo1)
      call denomvvoo_pp(eps,vvoo1)
!
!  do the contraction over all symmetries
!
      CALL CNTRCT('N','N',NVVT,NOOT,NVVT,A1,vvvvz,vvoo1,A0,vvoop,NREP)
!     write(*,*) 'Master: after CNTRCT.Chksum:',dot_product(vvoop,vvoop)
!
!  check if file already present. If yes, delete it.
!
      Inquire(file=preconfn,exist=isthere)
      If(isthere) then
         write(*,*) 'Deleting previous precontraction data.'
         open(unit=iounit, file=preconfn, status='old')
         close(iounit, status='delete')
      Endif
!
!  open corresponding file and write out contracted integrals.
!
      open(unit=iounit,file=preconfn, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='NEW')
      write(iounit) (vvoop(ixx), ixx=1,n1*rcw)
      close(iounit)
      write(*,*) 'File ',preconfn,' generated.'
!
!  free buffers immediately 
!
      deallocate(vvoop)
      deallocate(vvoo1)
      deallocate(vvvvz)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine compare_evc_r(io1,desrep,nroots,ladc,
     &                         davibase,fullbase)
      implicit none
!
! this routine evaluates eigenvectors for the available Davidson roots
! with respect to the full diagonalization results. It determines phases and
! absolute deviations
!
      integer        :: io1,desrep,nroots,ladc
      character*8    :: davibase
      character*7    :: fullbase 
!
! local variables
! 
      character*11                        :: daviname
      character*10                        :: fullname
      logical                             :: fullisthere
      integer                             :: io2
      integer                             :: i,j,k
      integer                             :: k1,k2
      integer                             :: nmain
      real*8                              :: lambda1,lambda2,dummy
      real*8,parameter                    :: autoev = 27.2113957D0
      real*8,allocatable,dimension (:)    :: e1,e2
      real*8                              :: devmod,devabs,devphs
!
! executable code
! 
      Call PST('Comparing D/F eigenvectors+')

      IF(desrep.GT.9) THEN
        WRITE(daviname,'(A8,A1,I2)') davibase,'.',desrep
        WRITE(fullname,'(A7,A1,I2)') fullbase,'.',desrep
      ELSE
        WRITE(daviname,'(A8,A2,I1)') davibase,'.0',desrep
        WRITE(fullname,'(A7,A2,I1)') fullbase,'.0',desrep
      ENDIF
!
! do we have an eigenvector file from fulldia ?
!
      fullisthere=.false.
      INQUIRE(FILE=fullname,EXIST=fullisthere)
      if(.not.fullisthere) then
        write(*,*) 'The eigenvector file from a full diagonalization'
        write(*,*) 'is not available. No comparison is made. Returning.'
        return
      endif
      write(*,*) 'Comparing eigenvectors for',nroots,' roots'
      write(*,*)
!
! comparing eigenvectors
!
      io2 = io1 + 1
      open(io1,file=daviname,access='sequential',
     &     status='unknown',form='unformatted')
      open(io2,file=fullname,access='sequential',
     &     status='unknown',form='unformatted')

      allocate(e1(ladc))
      allocate(e2(ladc))
      e1 = 0.0d0
      e2 = 0.0d0
!
! three integer reads from the Davidson file
!
      read(io1) k
      read(io1) nmain
      read(io1) k
      write(*,*) 'Main space has dimension',nmain
      write(*,*)

      write(*,'(18X,2A7,4X,3A8)') 'root1','root2',' devmod',
     &                          'devabs','devphs'
      do k=1,nroots
        read(io1) k1
        read(io2) k2
        if(.not.(k1.eq.k2 .and. k1.eq.k)) then
          write(*,*) 'Inconsistency in root numbering!'
          write(*,*) 'This is serious, continuation questionable!'
          write(*,*) 'k,k1,k2:',k,k1,k2
          exit  !leave this do loop
        endif
        read(io1) lambda1
        read(io2) lambda2
        read(io1) dummy
        read(io2) dummy
        read(io1) (e1(i),i=1,ladc)
        read(io2) (e2(i),i=1,ladc)
!
! compare vectors
!
        devmod = 0.0d0
        devabs = 0.0d0
        devphs = 0.0d0
        do i=1,ladc
          devmod = devmod + (abs(e1(i)) - abs(e2(i)))
          devabs = devabs + (e1(i) - e2(i))
          devphs = devphs + (e1(i) + e2(i))
        enddo
        write(*,'(A,I3,2X,2F12.6,2X,3F14.8)') 'Comparing root',k,
     &         lambda1,lambda2,devmod,devabs,devphs
        write(*,*)
        do i=1,min(10,nmain)
          write(*,*) '       ',e1(i),e2(i)
        enddo
        write(*,*) '..................'
        write(*,*)
      enddo
!
! release ressources
!
      deallocate(e1)
      deallocate(e2)

      close(io1)
      close(io2)

      return
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
      SUBROUTINE MATXMAT(n,a,b,c)
 
      IMPLICIT none
!
!---------------Description--------------------------------------------
!
! multiplies the matrix (vector) a with b and stores result in c
! This is done *elementwise*
!                           c = a * b
! If the real arrays carry complex numbers, we combine two numbers to
! one complex number
!
!---------------Calling variables--------------------------------------
!
      Integer                          :: n
      REAL*8,dimension(:)              :: a,b,c
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/complex.inc"
!
!---------------Local Variables--------------------------------------
!
      Integer         :: i
      Complex*16      :: s1,s2,s3
!
!---------------Executable --------------------------------------
!
!     call pst('in matxmat+')
      write(*,*) n,carith
      If(carith) then
        Do i=1,n
          s1=dcmplx(a(2*i-1),a(2*i))
          s2=dcmplx(b(2*i-1),b(2*i))
          s3 = s1 * s2
          c(2*i-1) = dble(s3)
          c(2*i)   = aimag(s3)
        Enddo
      Else
        Do i=1,n
          c(i) = a(i) * b(i)
        Enddo
      Endif

      Return
      End
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1T3_NTS (NREP,MULTB,DOINV,NPAIR,KFIE,LFIE,NTRIPL,
     &                       OFF,OFF1,BUF1,BUF2,DESREP,ICEIL)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     BUF(IJ,K>L:KLREP) <--> BUF(IJ K,L:LREP)
!     where the product of all four ireps can be non-totally symmetric
!     If (DOINV) inverse sort.
!
!                 <IJ| x |KL> = desrep
!
!---------------Calling variables--------------------------------------
!
      INTEGER                  NREP
      INTEGER                  MULTB(64,64,2)
      LOGICAL                  DOINV
      INTEGER                  NPAIR(NREP),KFIE(NREP)
      INTEGER                  LFIE(NREP),NTRIPL(NREP)
      INTEGER                  OFF(NREP+1),OFF1(32,32)
      REAL*8                   BUF1(*),BUF2(*)
      INTEGER                  DESREP,ICEIL
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
!
!---------------Local variables--------------------------------------
!
       integer summer,sumkl
!
!---------------Executable code--------------------------------------
!
      IF (.NOT.DOINV) THEN
         CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      ENDIF
      summer = 0
      IJKL = 1
      DO KLREP = 1, NREP

       IJREP = MULTB (DESREP+NREP,KLREP+NREP,2)
       IF(MULTB(ijrep+nrep,klrep+nrep,2).ne.desrep) 
     &        stop 'Symm error 1 in SRT_1T3_NTS'

       N = NPAIR(IJREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10

        IJKREP = MULTB (IJREP+NREP,KREP,2)
        IJLREP = MULTB (IJREP+NREP,LREP,2)
        if(multb(ijkrep,lrep,1).ne.desrep) stop 'Symm error 2'
        if(multb(ijlrep,krep,1).ne.desrep) stop 'Symm error 3'

        LOFF = (OFF(LREP) + OFF1(IJREP,KREP)) * RCW
        KOFF = (OFF(KREP) + OFF1(IJREP,LREP)) * RCW

        DO L = 1, LFIE(LREP)
         IJL  = (L-1) * N * RCW + 1
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
            IJK  = (K-1) * N * RCW + 1
            IJKL1 = LOFF + (L-1) * NTRIPL(IJKREP) * RCW + IJK 
            IJLK1 = KOFF + (K-1) * NTRIPL(IJLREP) * RCW + IJL 
            if(IJKL1+N.gt.iceil) write(*,*) '***alert A'
            if(IJLK1+N.gt.iceil) write(*,*) '***alert B'
            IF (DOINV) THEN
               CALL XCOPY (N,BUF1(IJKL1),1,BUF2(IJKL),1)
               CALL XAXPY (N,-A1,BUF1(IJLK1),1,BUF2(IJKL),1)
            ELSE
               CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJKL1),1)
               CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJLK1),1)
               CALL XSCAL (N,-A1,BUF2(IJLK1),1)
            ENDIF
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1T2_NTS (NREP,MULTB,DOINV,NPAIR1,IFIE,JFIE,
     &                       NPAIR2,OFF,OFF2,BUF1,BUF2,DESREP)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     BUF(I>J,KL:KLREP) <--> BUF(I,J KL:JKLREP)
!     If (DOINV) inverse sort.
!
!---------------Last modified------------------------------------------
!
!    by MP for nontotally symmetric cases.
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP,DESREP,MULTB(64,64,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF2(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
!
!---------------Executable code--------------------------------------
!
      IF (.NOT.DOINV) THEN
         CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      ENDIF
      IJKL = 0
      DO KLREP = 1, NREP

! new for nontotally symmetric cases:
       IJREP = MULTB (DESREP+NREP,KLREP+NREP,2)
       IF(MULTB(ijrep+nrep,klrep+nrep,2).ne.desrep) stop 'Error!'

       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IF (IREP.LT.JREP) GOTO 10

! new for nontotally symmetric cases:

        JKLREP = MULTB(JREP,KLREP+NREP,2)
        IKLREP = MULTB(IREP,KLREP+NREP,2)

        JKLOFF = (OFF(JKLREP)+OFF2(JREP,KLREP)*IFIE(IREP)) * RCW
        IKLOFF = (OFF(IKLREP)+OFF2(IREP,KLREP)*JFIE(JREP)) * RCW
        M2 = IFIE(IREP) * JFIE(JREP)
        DO J = 1, JFIE(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, IFIE(IREP)
            IJKL1 = JKLOFF + ((J-1) * IFIE(IREP) + I - 1) * RCW + 1
            JIKL1 = IKLOFF + ((I-1) * JFIE(JREP) + J - 1) * RCW + 1
            IF (DOINV) THEN
               CALL XCOPY(N,BUF1(IJKL1),M2,BUF2(IJKL+IJ),M)
               CALL XAXPY(N,-A1,BUF1(JIKL1),M2,BUF2(IJKL+IJ),M)
            ELSE
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(IJKL1),M2)
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(JIKL1),M2)
               CALL XSCAL(N,-A1,BUF2(JIKL1),M2)
            ENDIF
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Function get_file_unit()
!
! returns a unit number ifor file I/O that is not in use
!
      Implicit none

#define LUMIN 70
#define LUMAX 99

      Integer                   :: get_file_unit
      integer                   :: lu,iostat
      logical                   :: istaken
!
      istaken = .true.
      do lu = LUMIN, LUMAX
         inquire (unit=lu, opened=istaken, iostat=iostat)
         if (iostat.ne.0) cycle
         if (.not.istaken) exit
      end do
!
      if(istaken) then
        write(*,*) ' *** Attention. Reached',LUMAX
        write(*,*) ' *** Attention. No free file handle. Terminating.'
        call quit('I/O problem')
      else
        write(*,*) 'Obtained file handle',lu
      endif
      get_file_unit = lu
      return
      end Function
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

