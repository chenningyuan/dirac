!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end
!
#ifdef HAVE_MKL_BLAS
!radovan: workaround for strange crash when using MKL zdotc
#define ZDOTC z_dot_c
#define zdotc z_dot_c
#endif
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!      ****  Diagonalization Routines for the RELADC package  ****
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FULLDIAR(IOUNIT, FILEN4, FILEN5, INTBUF, LADC, KREP)

      implicit none

      INTEGER                       :: IOUNIT,INTBUF,LADC,KREP
      CHARACTER(len=6)              :: FILEN4, FILEN5
!
!---------------Description--------------------------------------------
!
!  For not too large symmetric matrices the full diagonalization is
!  performed for debugging purposes. The required arrays are allocated
!  dynamically in order to achieve maximal independence.
!
!  The only specific entities are two available file handles,
!  the length of the write chunks and the length of the ADC matrix.  
!
!---------------Local variables--------------------------------------
!
      REAL*8                          :: S
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      REAL*8, allocatable             :: A(:,:),WORK(:),EVL(:),BUF(:)
      INTEGER, allocatable            :: IOI(:),IOJ(:)
      INTEGER                         :: NBUFS,LWORK,TOTELE
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO
      INTEGER                         :: LUN_DIA,LUN_OFF,LUN_OUT
      CHARACTER(len=9)                :: longfilename
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!

      CALL PST('Entering FULL diagonalizer (debug)+')

      IF(LADC.gt.5000) THEN
        write(*,*) 'Matrix dimension too large for direct'
        write(*,*) 'diagonalization! Skipping this step...'
        return
      ENDIF

      LUN_DIA=IOUNIT
      LUN_OFF=IOUNIT+1
      LUN_OUT=IOUNIT+2




      IF(KREP.GT.9) THEN
        WRITE(longfilename,'(A7,I2)') 'FULLEV.',KREP
      ELSE
        WRITE(longfilename,'(A6,A2,I1)') 'FULLEV','.0',KREP
      ENDIF

      allocate(a(ladc,ladc))
      allocate(evl(ladc))
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      lwork=10*ladc
      allocate(work(lwork))
      
      A=0.0d0
!_________________________ reading section
!|
!|
!|
!|
      OPEN(LUN_DIA,FILE=FILEN4,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(LUN_OFF,FILE=FILEN5,FORM='UNFORMATTED',STATUS='UNKNOWN')

      REWIND(LUN_DIA); REWIND(LUN_OFF)
      READ(LUN_DIA,ERR=88,END=99) (A(IXX,IXX),IXX=1,LADC),NBUFS
      TOTELE = 0
      DO IREC = 1,NBUFS
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                        (IOI(IXX),IXX=1,INTBUF),
     &                        (IOJ(IXX),IXX=1,INTBUF),
     &                         NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          S    = BUF(K)

          IF(ICOL.GE.IROW) CALL QUIT('ICOL.GE.IROW should not happen')
 
          A(IROW,ICOL)=S
          TOTELE=TOTELE + 1
        ENDDO
      ENDDO
      WRITE(*,555) TOTELE
      CLOSE(LUN_DIA)
      CLOSE(LUN_OFF)
!|
!|
!|
!|________________________ end reading section

!
!  call symmetric diagonalizer
!  calculate eigenvectors being placed in A
!
      JOBZ='V'
      UPLO='L'

!     write(*,*) 'Calling DSYEV with lwork=',lwork
      CALL DSYEV( JOBZ, UPLO, LADC, A, LADC, EVL, WORK, LWORK, INFO )
!     write(*,*) 'AFTER  DSYEV with optimal lwork=',work(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in DSYEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAR (lapack dsyev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential', 
     &     form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)')
     &         '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
        do j=1,ladc
          write(lun_out,'(i6,3x,f17.12)') j,a(j,i)
        enddo
      enddo
      write(lun_out,'(A)') '     end of output'
      close(lun_out)

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

 88   STOP 'Error reading the diagonal (FULLDIAR)'
 99   STOP 'Unexpected end of diagonal file (FULLDIAR)'
 199  STOP 'Error reading off-diagonal elements (FULLDIAR)'

 555  FORMAT (1X,'Total number of elements in OD-file:',T40,I8)

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FULLDIAC(IOUNIT, FILEN4, FILEN5, INTBUF, LADC, KREP)
      implicit none

      INTEGER                       :: IOUNIT,INTBUF,LADC,KREP
      CHARACTER(len=6)              :: FILEN4, FILEN5
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
      COMPLEX*16                      :: S
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      COMPLEX*16,   allocatable       :: A(:,:),CWORK(:),BUF(:)
      REAL*8,       allocatable       :: EVL(:),RWORK(:)
      INTEGER,      allocatable       :: IOI(:),IOJ(:)
      INTEGER                         :: NBUFS,LWORK,TOTELE
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO
      INTEGER                         :: LUN_DIA,LUN_OFF,LUN_OUT
      CHARACTER(len=9)                :: longfilename
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!

      CALL PST('Entering FULL diagonalizer (debug)+')

      IF(LADC.gt.5000) THEN
        write(*,*) 'Matrix dimension too large for direct'
        write(*,*) 'diagonalization! Skipping this step...'
        return
      ENDIF

      LUN_DIA=IOUNIT
      LUN_OFF=IOUNIT+1
      LUN_OUT=IOUNIT+2



      IF(KREP.GT.9) THEN
        WRITE(longfilename,'(A7,I2)') 'FULLEV.',KREP
      ELSE
        WRITE(longfilename,'(A6,A2,I1)') 'FULLEV','.0',KREP
      ENDIF

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
!|
!|
!|
      OPEN(LUN_DIA,FILE=FILEN4,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(LUN_OFF,FILE=FILEN5,FORM='UNFORMATTED',STATUS='UNKNOWN')

      REWIND(LUN_DIA); REWIND(LUN_OFF)
      READ(LUN_DIA,ERR=88,END=99) (A(IXX,IXX),IXX=1,LADC),NBUFS
      TOTELE = 0
      DO IREC = 1,NBUFS
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                        (IOI(IXX),IXX=1,INTBUF),
     &                        (IOJ(IXX),IXX=1,INTBUF),
     &                         NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          S    = BUF(K)

          IF(ICOL.GE.IROW) CALL QUIT('ICOL.GE.IROW should not happen')
 
          A(IROW,ICOL)=S
          TOTELE=TOTELE + 1
        ENDDO
      ENDDO
      WRITE(*,555) TOTELE
      CLOSE(LUN_DIA)
      CLOSE(LUN_OFF)
!|
!|
!|
!|________________________ end reading section

!
!  call symmetric diagonalizer
!  calculate eigenvectors being placed in A
!
      JOBZ='V'
      UPLO='L'

!     write(*,*) 'Calling ZHEEV with lwork=',lwork
      CALL ZHEEV( JOBZ, UPLO, LADC, A, LADC, EVL, CWORK,
     &            LWORK, RWORK, INFO )
!     write(*,*) 'AFTER  ZHEEV with optimal lwork=',cwork(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in ZHEEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAC (lapack zheev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential',
     &     form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)')
     &         '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
        do j=1,ladc
          write(lun_out,'(i6,3x,2f17.12)') j,a(j,i)
        enddo
      enddo
      write(lun_out,'(A)') '     end of output'
      close(lun_out)
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

 88   STOP 'Error reading the diagonal (FULLDIAC)'
 99   STOP 'Unexpected end of diagonal file (FULLDIAC)'
 199  STOP 'Error reading off-diagonal elements (FULLDIAC)'

 555  FORMAT (1X,'Total number of elements in OD-file:',T40,I8)

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FULLDIAR2(IOUNIT, FILEN5, INTBUF, NBUFS, LADC, KREP)

      use memory_allocator
      use adc_fano_exchange
      use adc_mat

      implicit none

      INTEGER                       :: IOUNIT,INTBUF,NBUFS,LADC,KREP
      CHARACTER(len=6)              :: FILEN5
!
!---------------Description--------------------------------------------
!
!  For not too large symmetric matrices the full diagonalization is
!  performed for debugging purposes. The required arrays are allocated
!  dynamically in order to achieve maximal independence.
!
!  The only specific entities are two available file handles,
!  the length of the write chunks and the length of the ADC matrix.  
!
!---------------Local variables--------------------------------------
!
      REAL*8                          :: S
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      REAL*8, allocatable             :: A(:,:),WORK(:),EVL(:),BUF(:)
      INTEGER, allocatable            :: IOI(:),IOJ(:)
      INTEGER                         :: LWORK,TOTELE
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO
      INTEGER                         :: LUN_DIA,LUN_OFF,LUN_OUT
      CHARACTER(len=10)               :: longfilename
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!

      CALL PST('Entering (new format) REAL FULL diagonalizer+')
      write(*,*) 'Number of buffers to read:',nbufs

      IF(LADC.gt.5000) THEN
        write(*,*) 'Matrix dimension too large for direct'
        write(*,*) 'diagonalization! Exiting from this...'
        return
      ENDIF

      LUN_OFF=IOUNIT
      LUN_OUT=IOUNIT+1



      IF(KREP.GT.9) THEN
        WRITE(longfilename,'(A8,I2)') 'FULLEVD.',KREP
      ELSE
        WRITE(longfilename,'(A7,A2,I1)') 'FULLEVD','.0',KREP
      ENDIF

      allocate(a(ladc,ladc))
      allocate(evl(ladc))
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      lwork=10*ladc
      allocate(work(lwork))
      
      A=0.0d0
!_________________________ reading section
!|
!|
!|
!|
      OPEN(LUN_OFF,FILE=FILEN5,FORM='UNFORMATTED',STATUS='UNKNOWN')

      REWIND(LUN_OFF)
      TOTELE = 0
      DO IREC = 1,NBUFS
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                        (IOI(IXX),IXX=1,INTBUF),
     &                        (IOJ(IXX),IXX=1,INTBUF),
     &                         NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          S    = BUF(K)
!         write(*,'(A,2I4,F20.12)') 'Q<--',irow,icol,S

          IF(ICOL.GT.IROW) CALL QUIT('ICOL.GT.IROW should not happen')
 
          A(IROW,ICOL)=S
          TOTELE=TOTELE + 1
        ENDDO
      ENDDO
      WRITE(*,555) TOTELE
      CLOSE(LUN_OFF)
!|
!|
!|
!|________________________ end reading section

!
!  call symmetric diagonalizer
!  calculate eigenvectors being placed in A
!
      JOBZ='V'
      UPLO='L'

!     write(*,*) 'Calling DSYEV with lwork=',lwork
      CALL DSYEV( JOBZ, UPLO, LADC, A, LADC, EVL, WORK, LWORK, INFO )
!     write(*,*) 'AFTER  DSYEV with optimal lwork=',work(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in DSYEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAR (lapack dsyev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 20):'
      DO I=1,min(20,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential',
     &     form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)')
     &         '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
        do j=1,ladc
          write(lun_out,'(i6,3x,f17.12)') j,a(j,i)
        enddo
      enddo
      write(lun_out,'(A)') '     end of output'
      close(lun_out)
!
!  in case of a Fano run, write eigenvectors to array
!
      IF(reladc_md_isfano) THEN
        CALL alloc(fin_evecs,ladc,ladc,id='array of Fano final states')
        CALL alloc(fin_energies,ladc,id='2h1p energies of final states')
        WRITE(*,*) 'Saving final state vectors and energies to arrays'
        
        fin_energies(1:ladc)     = evl(1:ladc)
        fin_evecs(1:ladc,1:ladc) = a(1:ladc,1:ladc)

      END IF
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

 199  STOP 'Error reading off-diagonal elements (FULLDIAR)'

 555  FORMAT (1X,'Total number of elements in OD-file:',T40,I8)

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FULLDIAC2(IOUNIT, FILEN5, INTBUF, NBUFS, LADC, KREP)

      use memory_allocator
      use adc_fano_exchange
      use adc_mat

      implicit none

      INTEGER                       :: IOUNIT,INTBUF,NBUFS,LADC,KREP
      CHARACTER(len=6)              :: FILEN5
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
      COMPLEX*16                      :: S
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      COMPLEX*16,   allocatable       :: A(:,:),CWORK(:),BUF(:)
      REAL*8,       allocatable       :: EVL(:),RWORK(:)
      INTEGER,      allocatable       :: IOI(:),IOJ(:)
      INTEGER                         :: LWORK,TOTELE
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO
      INTEGER                         :: LUN_DIA,LUN_OFF,LUN_OUT
      CHARACTER(len=10)                :: longfilename
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!

      CALL PST('Entering (new format) FULL complex diagonalizer+')
      write(*,*) 'Number of buffers to read:',nbufs

      IF(LADC.gt.5000) THEN
        write(*,*) 'Matrix dimension too large for direct'
        write(*,*) 'diagonalization! Exiting from this ...'
        return
      ENDIF

      LUN_OFF=IOUNIT
      LUN_OUT=IOUNIT+1



      IF(KREP.GT.9) THEN
        WRITE(longfilename,'(A8,I2)') 'FULLEVD.',KREP
      ELSE
        WRITE(longfilename,'(A7,A2,I1)') 'FULLEVD','.0',KREP
      ENDIF

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
!|
!|
!|
      OPEN(LUN_OFF,FILE=FILEN5,FORM='UNFORMATTED',STATUS='UNKNOWN')

      REWIND(LUN_OFF)
      TOTELE = 0
      DO IREC = 1,NBUFS
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                        (IOI(IXX),IXX=1,INTBUF),
     &                        (IOJ(IXX),IXX=1,INTBUF),
     &                         NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          S    = BUF(K)

          IF(ICOL.GT.IROW) CALL QUIT('ICOL.GT.IROW should not happen')
 
          A(IROW,ICOL)=S
          TOTELE=TOTELE + 1
        ENDDO
      ENDDO
      WRITE(*,555) TOTELE
      CLOSE(LUN_OFF)
!|
!|
!|
!|________________________ end reading section

!
!  call symmetric diagonalizer
!  calculate eigenvectors being placed in A
!
      JOBZ='V'
      UPLO='L'

!     write(*,*) 'Calling ZHEEV with lwork=',lwork
      CALL ZHEEV( JOBZ, UPLO, LADC, A, LADC, EVL, CWORK,
     &            LWORK, RWORK, INFO )
!     write(*,*) 'AFTER  ZHEEV with optimal lwork=',cwork(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in ZHEEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAC (lapack zheev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential',
     &     form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)')
     &         '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
        do j=1,ladc
          write(lun_out,'(i6,3x,2f17.12)') j,a(j,i)
        enddo
      enddo
      write(lun_out,'(A)') '     end of output'
      close(lun_out)
!
!  in case of a Fano run, write eigenvectors to array
!
      IF(reladc_md_isfano) THEN
        CALL alloc(cfin_evecs,ladc,ladc,id='array of Fano final states')
        CALL alloc(fin_energies,ladc,id='2h1p energies of final states')
        WRITE(*,*) 'Saving final state vectors and energies to arrays'
        
        fin_energies(1:ladc)     = evl(1:ladc)
        cfin_evecs(1:ladc,1:ladc) = a(1:ladc,1:ladc)

      END IF
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

 199  STOP 'Error reading off-diagonal elements (FULLDIAC)'

 555  FORMAT (1X,'Total number of elements in OD-file:',T40,I8)

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DIAG_LANC(CC,KFREE,LFREE,IW,DOINCORE)

      use adc_mat      !obtain matrix parameters via descriptor file
      use adc_cfg

      IMPLICIT INTEGER (A-Z)

      REAL*8  CC(*)
      INTEGER KFREE,LFREE
      INTEGER IW
      LOGICAL DOINCORE
C
C---------------Description--------------------------------------------
C
C  This routine executes diagonalizer on any matrix (SIP/DIP) coming
C  
C  Infohandle is the base file handle.
C  the subsequent handles INFOHANDLE+1, +2 ... are available and
C  can be used. The allocation is as follows:
c
C  INFOHANDLE+1:   packed ADC matrix
C  INFOHANDLE+2:   old format: ADC matrix diagonal elements.
C  INFOHANDLE+3:   a free handle for the XLANCZOS routines
C  INFOHANDLE+4:   a free handle

C
C---------------Local variables--------------------------------------
C
c *** for infofile
      INTEGER INFOHANDLE
      INTEGER IONIZLEVEL
      INTEGER INTBUF
      INTEGER KREP
      INTEGER RCW
      INTEGER IOLDNEW
      INTEGER LADC
      CHARACTER*6 FILEADC
      CHARACTER*6 FILEDIA
      CHARACTER*6 FILECNF 
      CHARACTER*5 NMSPEC 
      INTEGER IRECL
      INTEGER NITER
      INTEGER NOCC
      INTEGER NEIGENV
      REAL*8 EVLLIM
      INTEGER IVALID   ! number of valid elements for Lanczos incore diagonalization.
c
c *** others
      LOGICAL ALLOC_ME,LFLAG,CARITH
      COMPLEX*16 CDUMMY
      REAL*8 RDUMMY
      INTEGER NINCWORDS
C
C---------------Local variables--------------------------------------
C
      COMMON/MATCHUNKS/NBUFS
c
C
C---------------Executable code--------------------------------------
C
!     WRITE(IW,*) 'Memory pointers on DIAG_LANC entry:',KFREE,LFREE
      KFRSAV=KFREE
      CALL MEMCHK('Beginning of diagonalization',CC,KFRSAV)
c
c default to a predefined buffer length. Will be changed
c in the future because we need a flexible buffer size,
c
      INTBUF = 5*1024*1024  !output buffer is now 5 MWORDS (40 MByte)
c
c  transfer required parameters from matrix descriptor module
c  to actual modules.
c
      IONIZLEVEL   =  reladc_md_ionizl
      INTBUF       =  reladc_md_intbuf
      INFOHANDLE   =  reladc_md_iobase
      KREP         =  reladc_md_desrep
      RCW          =  reladc_md_rcw
      IOLDNEW      =  reladc_md_ioldnew
      LADC         =  reladc_md_matdim
      FILEADC      =  reladc_md_fileadc
      FILEDIA      =  reladc_md_filediag
      FILECNF      =  reladc_md_filecnf
      NMSPEC       =  reladc_md_nmspec
      IRECL        =  reladc_md_irecl
      NITER        =  reladc_md_lnzitr
      NOCC         =  reladc_md_nmain
      NBUFS        =  reladc_md_nbufs    
      NEIGENV      =  reladc_md_neigenv
      IF (IONIZLEVEL.eq.1) THEN
        EVLLIM     =  reladc_sipprnt
      ELSE IF (IONIZLEVEL.eq.2) THEN
        EVLLIM     =  reladc_dipprnt
      ELSE IF (IONIZLEVEL.eq.3) THEN
        EVLLIM     =  reladc_excprnt
      END IF
c
c  open the ADC Matrix and configuration files and do diagonalization
c
      OPEN(INFOHANDLE+1,
     %     FILE=FILEADC,FORM='UNFORMATTED',STATUS='UNKNOWN')
      IF(IOLDNEW.EQ.1) OPEN(INFOHANDLE+2,
     %     FILE=FILEDIA,FORM='UNFORMATTED',STATUS='UNKNOWN')

c
c this is the absolute lower bound to which we will return as soon as the
c diagonalization has been finished. We must not go further down because
c some buffers are needed in the calling code.
c
      KFREESAVE=KFREE
      ALLOC_ME=.TRUE.
      NP = NOCC
      CARITH=.FALSE.
      IF(RCW.EQ.2) CARITH=.TRUE.

      CALL PST('    Initializing iterative diagonalizer+')
!
!  if repacking was successful we work from here with smaller ADC matrix
!  (in the double ionization case).
!
      IF(DOINCORE) THEN
        IF(IOLDNEW.eq.1) THEN
          WRITE(IW,*) 'Incore calculation for the old-format'
          WRITE(IW,*) 'ADC matrix not possible. Switching back to'
          WRITE(IW,*) 'Out of core diagonalization'
          DOINCORE = .FALSE.
        ELSE
          WRITE(IW,*) 'Trying (partial) incore diagonalization.'
        ENDIF
      ENDIF

      WRITE(IW,*) 'Memory pointers:',KFREE,LFREE
      WRITE(IW,*) 'Words available:',LFREE-KFREE

      WRITE(IW,*)
      WRITE(IW,*) 'ADC matrix is (1=SIP,2=DIP,3=EXC):',IONIZLEVEL
      WRITE(IW,*) 'Write buffer length is:           ',INTBUF
      WRITE(IW,*) 'Final states have symmetry:       ',KREP
      WRITE(IW,*) 'ADC matrix is (1=real,2=complex): ',RCW
      WRITE(IW,*) 'ADC matrix format (1=old/2=new):  ',IOLDNEW
      WRITE(IW,*) 'ADC matrix dimension (# of rows): ',LADC
      WRITE(IW,*) 'Packed matrix stored in:          ',FILEADC
      WRITE(IW,*) 'Diagonal elem. (old) stored in:   ',FILEDIA
      WRITE(IW,*) 'Configuration data stored in:     ',FILECNF
      WRITE(IW,*) 'Eigenvalues written to:           ',NMSPEC,'.#IREP'
      WRITE(IW,*) 'Record length of Config file:     ',IRECL
      WRITE(IW,*) 'Number of requ. Lanc iterations:  ',NITER
      WRITE(IW,*) 'Half band width (1h/2h/hp states):',NOCC
      WRITE(IW,*) 'Number of chunks in ADC file:     ',NBUFS
      WRITE(IW,*) 'Number of requested eigenvectors: ',NEIGENV
      WRITE(IW,*) 'Printing limit for eigenvalues:   ',EVLLIM
C
C  allocate required work space for band Lanczos iterations
C  real/complex case are accounted for by the RCW variable
C  wrk_diag is always real since we have a hermitian matrix
C  disk buffers are already allocated.
C
C  the OLD/NEW matrix formats are accounted for by the IOLDNEW
C  variable. In this case the file handle for the diagfile carries
C  the number -255
C

      ITAP1=INFOHANDLE+1    ! handle for the packed ADC matrix
      ITAP2=-255
      IF(IOLDNEW.EQ.1) ITAP2=INFOHANDLE+2  ! old matrix format !!


      ITAP3=INFOHANDLE+3    ! free handle
      ITAP5=INFOHANDLE+5    ! free handle
      ITAP6=INFOHANDLE+6    ! free handle

      N=2*NP+1
      CALL MEMOPT('INTE',WRK_IVX,N,CC,KFREE,LFREE,ALLOC_ME)
      write(*,*) 'allocation  1 done.' 
      N=(2*NP+1)*LADC*RCW
      CALL MEMOPT('REAL',WRK_P,N,CC,KFREE,LFREE,ALLOC_ME)
      write(*,*) 'allocation  2 done.' 
      N=(NP+1)*NITER*RCW
      CALL MEMOPT('REAL',WRK_T,N,CC,KFREE,LFREE,ALLOC_ME)
      write(*,*) 'allocation  3 done.' 
      N=LADC*RCW
      CALL MEMOPT('REAL',WRK_R,N,CC,KFREE,LFREE,ALLOC_ME)
      write(*,*) 'allocation  4 done.' 
      CALL MEMOPT('REAL',WRK_DIAG,N,CC,KFREE,LFREE,ALLOC_ME)
      write(*,*) 'allocation  5 done.' 
      N=INTBUF*RCW
      CALL MEMOPT('REAL',WRK_BUF,N,CC,KFREE,LFREE,ALLOC_ME)
      write(*,*) 'allocation  6 done.' 
      N=INTBUF
      CALL MEMOPT('INTE',WRK_OI,N,CC,KFREE,LFREE,ALLOC_ME)
      write(*,*) 'allocation  7 done.' 
      CALL MEMOPT('INTE',WRK_OJ,N,CC,KFREE,LFREE,ALLOC_ME)
      write(*,*) 'allocation  8 done.' 
      CALL MEMCHK('BEFORE  LANCZOS',CC,KFRSAV)
      write(*,*) 'allocation  9 done.' 
c
c  remember memory pointers before we do the huge incore parts.
c
      IF(DOINCORE) THEN
        KINCSAV = KFREE
        WRITE(IW,*) 'Memory pointers before incore diag:',
     &               KFREE,LFREE
        WRITE(IW,*) 'Saving free memory:',KINCSAV
        WRITE(IW,*) 'Available buffer for incore diag:'
        WRITE(IW,*) (LFREE-10),' words'
        if(RCW.eq.1) then
           N=(LFREE - 10)/3
           WRITE(IW,*) 'We can store',N,' real numbers'
        else
           N=(LFREE - 10)/4
           WRITE(IW,*) 'We can store',N,' complex numbers'
        endif
        IF(N.lt.INTBUF) THEN
          WRITE(IW,*) 'Buffer for inc. diag too small!'
          WRITE(IW,*) 'Switching back to out-of-core diag.'
          DOINCORE = .FALSE.
        ELSE
          CALL MEMOPT('REAL',WRK_IDR,N*RCW,CC,KFREE,LFREE,ALLOC_ME)
          CALL MEMOPT('INTE',WRK_IDI,N,CC,KFREE,LFREE,ALLOC_ME)
          CALL MEMOPT('INTE',WRK_IDJ,N,CC,KFREE,LFREE,ALLOC_ME)
          CALL MEMCHK('ICD',CC,KFRSAV)
          NINCWORDS = N
          IF(CARITH) THEN
            CALL INITINCORE_C(ITAP1,ITAP3,NBUFS,NINCWORDS,INTBUF,
     &                        CC(WRK_BUF),CC(WRK_OI),CC(WRK_OJ),
     &                        CC(WRK_IDR),CC(WRK_IDI),CC(WRK_IDJ),
     &                        IVALID)
          ELSE
            CALL INITINCORE_R(ITAP1,ITAP3,NBUFS,NINCWORDS,INTBUF,
     &                        CC(WRK_BUF),CC(WRK_OI),CC(WRK_OJ),
     &                        CC(WRK_IDR),CC(WRK_IDI),CC(WRK_IDJ),
     &                        IVALID)
          ENDIF
        ENDIF
      ENDIF
C
C  enter real resp. complex Band Lanczos iterations
C  for symmetric resp. Hermitian matrices
C  In the first run: no eigenvectors because no eigenvalues are
C  available
C
      LFLAG = .FALSE.

      IF(DOINCORE.EQV..true.) THEN
        IF(CARITH.EQV..TRUE.) THEN
        CALL ZLANCZOS_IC(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
     &                   CC(WRK_P), CC(WRK_T), CC(WRK_R),
     &                   CC(WRK_DIAG), CDUMMY, CDUMMY, CC(WRK_BUF),
     &                   CC(WRK_OI), CC(WRK_OJ), CDUMMY,LFLAG,
     &                   INTBUF,ITAP1, ITAP2, ITAP3, NINCWORDS,
     &                   cc(wrk_idr), cc(wrk_idi), cc(wrk_idj),ivalid)
        CALL MEMCHK('AFTER COMPLEX INCORE-LANCZOS',CC,KFRSAV)
        ELSE
        CALL RLANCZOS_IC(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
     &          CC(WRK_P), CC(WRK_T), CC(WRK_R),
     &          CC(WRK_DIAG), RDUMMY, RDUMMY, CC(WRK_BUF),
     &          CC(WRK_OI), CC(WRK_OJ), RDUMMY,LFLAG,
     &          INTBUF,ITAP1, ITAP2, ITAP3, NINCWORDS,
     &          cc(wrk_idr), cc(wrk_idi), cc(wrk_idj),ivalid)
        CALL MEMCHK('AFTER REAL INCORE-LANCZOS',CC,KFRSAV)
        ENDIF
        CALL MEMREL('After incore Lanczos',CC,1,KINCSAV,KFREE,LFREE)
        WRITE(IW,*)  'Memory after incore Lanczos released.'
      ELSE 
        IF(CARITH.EQV..TRUE.) THEN
        CALL ZLANCZOS(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
     &          CC(WRK_P), CC(WRK_T), CC(WRK_R),
     &          CC(WRK_DIAG), CDUMMY, CDUMMY, CC(WRK_BUF),
     &          CC(WRK_OI), CC(WRK_OJ), CDUMMY,LFLAG,
     &          INTBUF,ITAP1, ITAP2, ITAP3)
        CALL MEMCHK('AFTER COMPLEX LANCZOS',CC,KFRSAV)
        ELSE
        CALL RLANCZOS(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
     &          CC(WRK_P), CC(WRK_T), CC(WRK_R),
     &          CC(WRK_DIAG), RDUMMY, RDUMMY, CC(WRK_BUF),
     &          CC(WRK_OI), CC(WRK_OJ), RDUMMY,LFLAG,
     &          INTBUF,ITAP1, ITAP2, ITAP3)
        CALL MEMCHK('AFTER REAL LANCZOS',CC,KFRSAV)
        ENDIF
      ENDIF
C
C  from here the Lanczos matrix T is available. we acquire work
C  space for the LAPACK Band-diagonalizer.
C  After the BANDDIA runs the corresponding eigenVECTORS of the
C  Lanczos matrix are available in ZEV  resp. in CC(WRK_ZEV)
C  ITAP5 will be provided as a free file handle. 
C  It is only used locally and reusable after finishing RBANDDIA/ZBANDDIA
C  This file handle is used for storing the T eigenvector matrix
C
      N=NITER
      CALL MEMOPT('REAL',WRK_EVL,N,CC,KFREE,LFREE,ALLOC_ME)
      N=NITER*NITER*RCW
      CALL MEMOPT('REAL',WRK_ZEV,N,CC,KFREE,LFREE,ALLOC_ME)
      N=4*NITER*RCW
      CALL MEMOPT('REAL',WRK_WORK1,N,CC,KFREE,LFREE,ALLOC_ME)
      N=6*NITER
      CALL MEMOPT('REAL',WRK_RWORK1,N,CC,KFREE,LFREE,ALLOC_ME)

      IF(CARITH) THEN
        CALL ZBANDDIA(ITAP5,NITER,NP,NP+1,NEIGENV,CC(WRK_T),CC(WRK_EVL),
     &           CC(WRK_ZEV),CC(WRK_WORK1),CC(WRK_RWORK1))
        CALL MEMCHK('AFTER ZBANDDIA',CC,KFRSAV)
      ELSE
        CALL RBANDDIA(ITAP5,NITER,NP,NP+1,NEIGENV,CC(WRK_T),CC(WRK_EVL),
     &           CC(WRK_ZEV),CC(WRK_RWORK1))
        CALL MEMCHK('AFTER RBANDDIA',CC,KFRSAV)
      ENDIF
c
c remove spurious/multiple eigenvalues and
c analyze eigenvectors of band matrix
c open corresponding configuration file
c
      OPEN(ITAP6,FILE=FILECNF,
     &     ACCESS='DIRECT',RECL=IRECL,STATUS='UNKNOWN')

      N=2*NP*NITER*RCW
      CALL MEMOPT('REAL',WRK_Z,N,CC,KFREE,LFREE,ALLOC_ME)
      N=NITER
      CALL MEMOPT('INTE',WRK_IND,N,CC,KFREE,LFREE,ALLOC_ME)
      CALL MEMOPT('INTE',WRK_MLT,N,CC,KFREE,LFREE,ALLOC_ME)
      CALL MEMOPT('INTE',WRK_LSP,N,CC,KFREE,LFREE,ALLOC_ME)

      IF(CARITH) THEN
         CALL ZBLANA(IW,KREP,niter,np,nocc,ionizlevel,
     &        cc(wrk_evl),cc(wrk_zev),cc(wrk_z),
     &        cc(wrk_t),cc(wrk_ind),cc(wrk_mlt),cc(wrk_lsp),
     &        ITAP5,ITAP6,EVLLIM,NMSPEC)
      ELSE
         CALL RBLANA(IW,KREP,niter,np,nocc,ionizlevel,
     &        cc(wrk_evl),cc(wrk_zev),cc(wrk_z),
     &        cc(wrk_t),cc(wrk_ind),cc(wrk_mlt),cc(wrk_lsp),
     &        ITAP5,ITAP6,EVLLIM,NMSPEC)
      ENDIF
c
c  close configuration file
c
      CLOSE(ITAP6)

      CALL MEMCHK('AFTER BLANA routines',CC,KFRSAV)
c
c check if user wants to calculate eigenvectors of full ADC matrix
c NGT is the actual number of nonspurious/unique eigenvalues.
c *** important: the corresponding intensities and corrections for
c the multiple eigenvalues are accounted for in the short eigenvector
c array z. This array has to be transferred to the eigenvector
c calculation.
C The corresponding indices sit in the array IND2 which has to be
C transferred to the eigenvector routine.
c
c
      IF(NEIGENV.GT.0) THEN

        N=LADC*NEIGENV*RCW
        CALL MEMOPT('REAL',WRK_AEV,N,CC,KFREE,LFREE,ALLOC_ME)
        N=LADC*RCW
        CALL MEMOPT('REAL',WRK_OEV,N,CC,KFREE,LFREE,ALLOC_ME)
c
c  next redo Lanczos now with the FULL eigenvectors of T and create
c  the LONG eigenvectors for the TOTAL ADC matrix
c  this is executed if asked for long eigenvectors.
c  file 1 and 2 (the ADC matrix/diag files) are still open.
c
        LFLAG = .TRUE.

        IF(DOINCORE.EQV..true.) THEN
          IF(CARITH.EQV..TRUE.) THEN
            CALL ZLANCZOS_IC(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
     &                       CC(WRK_P), CC(WRK_T), CC(WRK_R),
     &                       CC(WRK_DIAG), CC(WRK_ZEV), CC(WRK_AEV),
     &                       CC(WRK_BUF),CC(WRK_OI), CC(WRK_OJ),
     &                       CC(WRK_EVL),LFLAG, INTBUF,
     &                       ITAP1, ITAP2, ITAP3, NINCWORDS,
     &                       cc(wrk_idr), cc(wrk_idi), cc(wrk_idj),
     &                       ivalid)
          ELSE
            CALL RLANCZOS_IC(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
     &                       CC(WRK_P), CC(WRK_T), CC(WRK_R),
     &                       CC(WRK_DIAG), CC(WRK_ZEV), CC(WRK_AEV),
     &                       CC(WRK_BUF),CC(WRK_OI), CC(WRK_OJ),
     &                       CC(WRK_EVL),LFLAG, INTBUF,
     &                       ITAP1, ITAP2, ITAP3, NINCWORDS,
     &                       cc(wrk_idr), cc(wrk_idi), cc(wrk_idj),
     &                       ivalid)
          ENDIF
          CALL MEMREL('After incore Lanczos',CC,1,KINCSAV,KFREE,LFREE)
          WRITE(IW,*)  'Memory after incore Lanczos released.'
        ELSE
          IF(CARITH.EQV..TRUE.) THEN
            CALL ZLANCZOS(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
     &                    CC(WRK_P), CC(WRK_T), CC(WRK_R),
     &                    CC(WRK_DIAG), CC(WRK_ZEV), CC(WRK_AEV),
     &                    CC(WRK_BUF),CC(WRK_OI), CC(WRK_OJ),
     &                    CC(WRK_EVL),LFLAG, INTBUF,
     &                    ITAP1, ITAP2, ITAP3)
          ELSE
            CALL RLANCZOS(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
     &                    CC(WRK_P), CC(WRK_T), CC(WRK_R),
     &                    CC(WRK_DIAG), CC(WRK_ZEV), CC(WRK_AEV),
     &                    CC(WRK_BUF),CC(WRK_OI), CC(WRK_OJ),
     &                    CC(WRK_EVL),LFLAG, INTBUF,
     &                    ITAP1, ITAP2, ITAP3)
          ENDIF
        ENDIF
      ENDIF
c
c  close and delete the (diag) and the packed ADC matrix file
c
      CLOSE(ITAP1,STATUS='DELETE')
      IF(IOLDNEW.EQ.1) CLOSE(ITAP2,STATUS='DELETE')

      CALL MEMREL('ADC-memrel',CC,KFRSAV,KFREESAVE,KFREE,LFREE)
c
c  from here the file handles itap1 - itap3 are available again and we
c  can perform the detailed analysis of the long ADC eigenvectors together
c  with their configurations
c
c  All the memory is released and we do the analysis with a F90 routine.
c  irecl still contains the record length of the configuration file.
c
      itap1=infohandle+1
      itap2=infohandle+2
      itap3=infohandle+3
c
c  call long vector analyzer
c  handles 1-3 are exclusively used there, all files are closed
c  on exit.
c
      IF(NEIGENV.GT.0) THEN
        IF(CARITH) THEN
          CALL LVECANALY_C(itap1,itap2,itap3,krep,ladc,
     &                     nocc,irecl,filecnf,ionizlevel)
        ELSE
          CALL LVECANALY_R(itap1,itap2,itap3,krep,ladc,
     &                     nocc,irecl,filecnf,ionizlevel)
        ENDIF
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RLANCZOS(IRREP, LADC, NP, NITER, NEIGENV, IVX,
     &                    P, T, R, DIAG, Z, AEV, BUF, OI, OJ, CEVL, 
     &                    CFLAG, INTBUF, IO_OFF, IO_DIA, IO_LAN)
C
      IMPLICIT INTEGER (A-Z)
 
      INTEGER IRREP,LADC,NP,NITER,NEIGENV,IVX(2*NP+1)
      INTEGER INTBUF,IO_OFF,IO_DIA,IO_LAN
      REAL*8 P(LADC,2*NP+1), T(NP+1,NITER), R(LADC), DIAG(LADC)
      REAL*8 Z(NITER,NITER), AEV(LADC,NEIGENV), BUF(INTBUF)
      REAL*8 CEVL(NITER)
      INTEGER OI(INTBUF),OJ(INTBUF)
      LOGICAL CFLAG
      REAL*8 CPUTIM1,CPUTIM2,timer1,timer2,timertot
C
C---------------Description--------------------------------------------
C
C  Perform the real Lanczos iteration, generate matrix T_ij
C  the arrays BUF,OI,OJ are buffers for the adc matrix reading process
C  the file handles are IO_DIA and IO_OFF
C  If only eigenvalues are to be calculated, resp. in the first call of
C  this routine the Z-array is empty and will not be referenced !
C  In an eigenvector run Z contains the eigenvectors of the band matrix
C  diagonalization and is multiplied with the Lanczos vectors in order
C  to obtain the approximative eigenvectors of the original problem.
C  the variable .flag. determines if the additional eigenvector run
C  should be performed which is equivalent to the presence of the Z
C  vectors.
C  We write out the Lanczos matrix to file because it will be needed
C  for possible eigenvector calculations.
C  If an eigenvalue calculation is asked for we have to read in the
C  matrix of purged eigenvectors from ZMAT.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables and function definitions ------------
C
      REAL*8 Y,DDOT,Z1
      CHARACTER*10 LEVECFN
      real*8 rdummy
C
C---------------Executable code--------------------------------------
c
c  create symmetry-specific name for the long eigenvector file.
c  The eigenvector file of the Lanczos matrix needed for the
c  creation of the long eigenvectors is generically stored.
c
      IF(IRREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',IRREP
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',IRREP
      ENDIF

      IF(CFLAG) THEN
        CALL PST('Entering real Lanczos for Eigenvectors+')
      ELSE
        CALL PST('Entering real Lanczos for Eigenvalues+')
      ENDIF

      WRITE(IW,'(A30,I7)') ' Symmetry:',IRREP
      WRITE(IW,'(A30,I7)') ' Matrix dimension:',LADC
      WRITE(IW,'(A30,I7)') ' Band width:',NP
      WRITE(IW,'(A30,I7)') ' Number of iterations:',NITER
      IF(CFLAG) THEN
        WRITE(IW,'(A30,I7)') ' Requested eigenvectors:',NEIGENV
c       WRITE(IW,'(A30,A)') ' Long eigenvectors in:',LEVECFN
      ENDIF
      WRITE(IW,*)
c
c  ... clear large eigenvector array if requested
c
      IF(CFLAG) THEN
        DO I=1,NEIGENV
          DO J=1,LADC
            AEV(J,I)=0.0D0
          ENDDO
        ENDDO
      ENDIF
c
c  ... read in purged eigenvector matrix TMATEVC for calculation of
c      the long eigenvectors. These are folded into the Lanczos
c      iterations. The third available file handle is used for that.
c      After the read process the handle can be reused and file is
c      deleted because this is the last position in the code where
c      it is needed.
c
      IF(CFLAG) THEN
        OPEN(IO_LAN,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IO_LAN)
        READ(IO_LAN) N1
        READ(IO_LAN) N2
        WRITE(IW,*) 'Found',N1,' Lanczos eigenvectors of length',
     &              N2,' on TMATEVC'
        IF(NEIGENV.gt.N1) THEN
          WRITE(*,*) 'More EVCs requested than available.'
          WRITE(*,*) 'Adjusting to the value given above.'
          NEIGENV = N1
        ENDIF
        IF(N2.NE.NITER)
     &    CALL QUIT('Incompatible T matrix eigenvector lengths!')
        DO I=1,NEIGENV
          READ(IO_LAN) CEVL(I)
          READ(IO_LAN) (Z(IX,I),IX=1,N2)
        ENDDO
        CLOSE(IO_LAN,STATUS='DELETE')
      ENDIF

c
c  ... clear Lanczos matrix
c
      DO I=1,NITER
        DO J=1,NP+1
          T(J,I)=0.0D0
        ENDDO
      ENDDO
c
c  ... clear r/p-vector array and initialize cartesian start vectors
c      start vectors are orthonormal and sit in the fields
c      np+1,np+2,...,2*np
c
      DO I=1,2*NP+1
        DO J=1,LADC
          P(J,I)=0.0D0
        ENDDO
      ENDDO
      DO J=1,LADC
        R(J)=0.0D0
      ENDDO
      DO I=1,NP
        P(I,I+NP)=1.0D0
      ENDDO

c
c  ... initialize vector mapping array
c
      DO I=1,2*NP+1
        IVX(I)=I
      ENDDO
c
c  ... and time count
c
      timertot=0.0d0
      CALL CPUUSED(CPUTIM1)


c----------------------------- lanczos loop start ----------------------

      DO 99 IL = 1,NITER
c
c
      IF(mod(il,100).eq.0) THEN
        WRITE(IW,*) '---  Band Lanczos iteration ',IL,' of',niter
        CALL CPUUSED(CPUTIM2)
        timer1 = cputim2 - cputim1
        write(*,*) 'Time per 100 iterations:',timer1
        timertot = timertot + timer1
        cputim1 = cputim2
      ENDIF
c
c  form r = Hq_i - sum t_ji*q_j
c  the j index in the sum runs from i-p up to i-1
c  and the loop is only executed if previous vectors exist.
c  the actual vector to be multiplied by A always sits in
c  P(*,np+1)
c
c
c  since for the double ionization the diagonal elements are
c  not stored separately, we define an alternative matrix-vector
c  multiplication routine which directly multiplies the matrix.
c  If IO_DIA has value of -255 (which is surely no valid file handle)
c  we know that the alternative routine should be accessed.
c
c  att: we need the buffer number in RMATBLK_ND! In order to avoid a lot
c  of parameter rewriting we communicate this via a common block.
c
        IF(IO_DIA.eq.-255) THEN
          CALL RMATBLK_ND(LADC, P(1,IVX(NP+1)), R, BUF,
     &                    OI, OJ, INTBUF, IO_OFF)
        ELSE
          CALL RMATBLK(LADC, P(1,IVX(NP+1)), R, DIAG, BUF,
     &                 OI, OJ, INTBUF, IO_DIA, IO_OFF)
        ENDIF

c  subtract all the previous components

        ICC = 1
        DO IP = MAX0(IL-NP,1), IL-1
          ITROW = ICC+1
          ITCOL = IL - ICC
          DO J=1,LADC
            R(J) = R(J)-T(ITROW,ITCOL)*P(J,IVX(NP-ICC+1))
          ENDDO
          ICC = ICC + 1
        ENDDO
c
c  form the t_ji = q_j^+ * r and r=r-t_ji*q_j for all
c  j=i...i+p-1
c
        DO IP = 1,NP
          Y = DDOT(LADC,P(1,IVX(NP+IP)),1,R,1)
          T(IP,IL) = Y
          DO J=1,LADC
            R(J) = R(J) - Y*P(J,IVX(NP+IP))
          ENDDO
        ENDDO
c
c  now r contains the whole RHS of the band-lanczos equation
c  and will be normalized
c  form t_i+p,i = |r|, q_i+p,i = r/|r|
c
        Y = DDOT(LADC,R,1,R,1)
        Y = DSQRT(Y)
        T(NP+1,IL) = Y

        Y = 1.0d0/Y
        DO J=1,LADC
          P(J,IVX(2*NP+1)) = R(J) * Y
        ENDDO
c
c  if applicable form the approximate eigenvector
c  with the first Lanczos vector Q_1
c  Important: the ith ready Lanczos vector is available at
c  P(*,IVX(NP+1)) exactly at this position of the
c  iteration. In case of an eigenvector run its contributions
c  to ALL eigenvectors have to be taken into account. After
c  cycling the contributions of the i+1. lanczos vector will
C  be folded into the final eigenvectors
c
        IF(CFLAG) THEN
          DO IOM = 1,NEIGENV
            Z1 = Z(IL,IOM)
            DO ICOL = 1,LADC
              AEV(ICOL,IOM) = AEV(ICOL,IOM) + Z1*P(ICOL,IVX(NP+1))
            ENDDO
          ENDDO
        ENDIF
c
c  now cycle the mapping array IVX
c  and start new main iteration
c
        CALL IVXROT(IVX,2*NP+1)

 99   CONTINUE
c
c  write out long ADC eigenvectors when eigenvector run finished
c
      IF(CFLAG) THEN
        OPEN(IO_LAN, FILE=LEVECFN, FORM='UNFORMATTED', 
     &       ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
        WRITE(IO_LAN) LADC,NEIGENV
        DO J=1,NEIGENV
          WRITE(IO_LAN) CEVL(J)
          WRITE(IO_LAN) (AEV(I,J),I=1,LADC)
        ENDDO
        CLOSE(IO_LAN)
        WRITE(IW,*) ' ADC LONG EIGENVECTORS WRITTEN.'
      ENDIF

      write(*,*) 'Total time spent in Matrix * Vector:',timertot
      write(*,*) 'Average time for Matrix * Vector:',
     &            timertot/dble(niter)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZLANCZOS(IRREP, LADC, NP, NITER, NEIGENV, IVX,
     &                    P, T, R, DIAG, Z, AEV, BUF, OI, OJ, CEVL, 
     &                    CFLAG, INTBUF, IO_OFF, IO_DIA, IO_LAN)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER IRREP,LADC,NP,NITER,NEIGENV,IVX(2*NP+1)
      INTEGER INTBUF,IO_OFF,IO_DIA,IO_LAN
      COMPLEX*16 P(LADC,2*NP+1), T(NP+1,NITER), R(LADC), DIAG(LADC)
      COMPLEX*16 Z(NITER,NITER), AEV(LADC,NEIGENV), BUF(INTBUF)
      REAL*8 CEVL(NITER)
      INTEGER OI(INTBUF),OJ(INTBUF)
      LOGICAL CFLAG
C
C---------------Description--------------------------------------------
C
C  Perform the complex Lanczos iteration, generate matrix T_ij
C  the arrays BUF,OI,OJ are buffers for the adc matrix reading process
C  the file handles are IO_DIA and IO_OFF
C  If only eigenvalues are to be calculated, resp. in the first call of
C  this routine the Z-array is empty and will not be referenced !
C  In an eigenvector run Z contains the eigenvectors of the band matrix
C  diagonalization and is multiplied with the Lanczos vectors in order
C  to obtain the approximative eigenvectors of the original problem.
C  the variable .flag. determines if the additional eigenvector run
C  should be performed which is equivalent to the presence of the Z
C  vectors.
C  We write out the Lanczos matrix to file because it will be needed
C  for possible eigenvector calculations.
C  If an eigenvalue calculation is demanded we have to read in the
C  matrix of purged eigenvectors from ZMAT.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables and function definitions ------------
C
      COMPLEX*16 Y,ZDOTC,Z1
      CHARACTER*10 LEVECFN
      real*8 rdummy
C
C---------------Executable code--------------------------------------
C
      IF(IRREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',IRREP
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',IRREP
      ENDIF

      IF(CFLAG) THEN
        CALL PST('Entering complex Lanczos for Eigenvectors+')
      ELSE
        CALL PST('Entering complex Lanczos for Eigenvalues+')
      ENDIF

      WRITE(IW,'(A30,I7)') ' Symmetry:',IRREP
      WRITE(IW,'(A30,I7)') ' Matrix dimension:',LADC
      WRITE(IW,'(A30,I7)') ' Band width:',NP
      WRITE(IW,'(A30,I7)') ' Number of iterations:',NITER
      IF(CFLAG) THEN
        WRITE(IW,'(A30,I7)') ' Requested eigenvectors:',NEIGENV
c       WRITE(IW,'(A30,A)') ' Long eigenvectors in:',LEVECFN
      ENDIF
      WRITE(IW,*)
c
c  ... clear approx eigenvector array if requested
c
      IF(CFLAG) THEN
        DO I=1,NEIGENV
          DO J=1,LADC
            AEV(J,I)=(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDIF
c
c  ... read in purged eigenvector matrix TEVC for calculation of
c      the long eigenvectors. These are folded into the Lanczos
c      iterations. The third available file handle is used for that.
c      After the read process the handle can be reused and file is
c      deleted because this is the last position in the code where
c      it is needed.
c
      IF(CFLAG) THEN
        OPEN(IO_LAN,FILe='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IO_LAN)
        READ(IO_LAN) N1
        READ(IO_LAN) N2
        WRITE(IW,*) 'Found',N1,' Lanczos eigenvectors of length',
     &              N2,' on TMATEVC'
        IF(NEIGENV.gt.N1) THEN
          WRITE(*,*) 'More EVCs requested than available.'
          WRITE(*,*) 'Adjusting to the value given above.'
          NEIGENV = N1
        ENDIF
        IF(N2.NE.NITER)
     &    CALL QUIT('Incompatible T matrix eigenvector lengths!')
        DO I=1,NEIGENV
          READ(IO_LAN) CEVL(I)
          READ(IO_LAN) (Z(IX,I),IX=1,N2)
        ENDDO
        CLOSE(IO_LAN,STATUS='DELETE')
      ENDIF
c
c  ... clear Lanczos matrix
c
      DO I=1,NITER
        DO J=1,NP+1
          T(J,I)=(0.0D0,0.0D0)
        ENDDO
      ENDDO
c
c  ... clear r/p-vector array and initialize start vectors
c      start vectors are orthonormal and sit in the fields
c      np+1,np+2,...,2*np
c
      DO I=1,2*NP+1
        DO J=1,LADC
          P(J,I)=(0.0D0,0.0D0)
        ENDDO
      ENDDO
      DO J=1,LADC
        R(J)=(0.0D0,0.0D0)
      ENDDO
      DO I=1,NP
        P(I,I+NP)=(1.0D0,0.0D0)
      ENDDO

c
c  ... initialize vector mapping array
c
      DO I=1,2*NP+1
        IVX(I)=I
      ENDDO

c----------------------------- lanczos loop start ----------------------

      DO 99 IL = 1,NITER
c
c
      IF(mod(il,100).eq.0) THEN
        WRITE(IW,*) '---  Band Lanczos iteration ',IL,' of',niter
      ENDIF
c
c  form r = Hq_i - sum t_ji*q_j
c  the j index in the sum runs from i-p up to i-1
c  and the loop is only executed if previous vectors exist.
c  the actual vector to be multiplied by A always sits in
c  P(*,np+1)
c
        IF(IO_DIA.eq.-255) THEN
          CALL CMATBLK_ND(LADC, P(1,IVX(NP+1)), R, BUF,
     &                 OI, OJ, INTBUF, IO_OFF)
        ELSE
          CALL CMATBLK(LADC, P(1,IVX(NP+1)), R, DIAG, BUF,
     &                 OI, OJ, INTBUF, IO_DIA, IO_OFF)
        ENDIF


c  subtract all the previous components

        ICC = 1
        DO IP = MAX0(IL-NP,1), IL-1
          ITROW = ICC+1
          ITCOL = IL - ICC
          DO J=1,LADC
            R(J) = R(J)-DCONJG(T(ITROW,ITCOL))*P(J,IVX(NP-ICC+1))
          ENDDO
          ICC = ICC + 1
        ENDDO
c
c  form the t_ji = q_j^+ * r and r=r-t_ji*q_j for all
c  j=i...i+p-1
c
        DO IP = 1,NP
          Y = ZDOTC(LADC,P(1,IVX(NP+IP)),1,R,1)
          T(IP,IL) = Y
          DO J=1,LADC
            R(J) = R(J) - Y*P(J,IVX(NP+IP))
          ENDDO
        ENDDO
c
c  now r contains the whole RHS of the band-lanczos equation
c  and will be normalized
c  form t_i+p,i = |r|, q_i+p,i = r/|r|
c
        Y = ZDOTC(LADC,R,1,R,1)
c
        Y = SQRT(Y)
        T(NP+1,IL) = Y

        Y = (1.0d0,0.0d0)/Y
        DO J=1,LADC
          P(J,IVX(2*NP+1)) = R(J) * Y
        ENDDO
c
c  if applicable form the approximate eigenvector
c  with the first Lanczos vector Q_1
c  Important: the ith ready Lanczos vector is available at
c  P(*,IVX(NP+1)) exactly at this position of the
c  iteration. In case of an eigenvector run its contributions
c  to ALL eigenvectors have to be taken into account. After
c  cycling the contributions of the i+1. lanczos vector will
C  be folded into the final eigenvectors
c
        IF(CFLAG) THEN
          DO IOM = 1,NEIGENV
            Z1 = Z(IL,IOM)
            DO ICOL = 1,LADC
              AEV(ICOL,IOM) = AEV(ICOL,IOM) + Z1*P(ICOL,IVX(NP+1))
            ENDDO
          ENDDO
        ENDIF
c
c  now cycle the mapping array IVX
c  and start new main iteration
c
        CALL IVXROT(IVX,2*NP+1)

 99   CONTINUE
c
c  write out long ADC eigenvectors when eigenvector run finished
c
      IF(CFLAG) THEN
        OPEN(IO_LAN, FILE=LEVECFN, FORM='UNFORMATTED',
     &       ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
        WRITE(IO_LAN) LADC,NEIGENV
        DO J=1,NEIGENV
          WRITE(IO_LAN) CEVL(J)
          WRITE(IO_LAN) (AEV(I,J),I=1,LADC)
        ENDDO
        CLOSE(IO_LAN)
        WRITE(IW,*) ' ADC LONG EIGENVECTORS WRITTEN.'
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RLANCZOS_IC(IRREP, LADC, NP, NITER, NEIGENV, IVX,
     &                       P, T, R, DIAG, Z, AEV, BUF, OI, OJ, CEVL, 
     &                       CFLAG, INTBUF, IO_OFF, IO_DIA, IO_LAN,
     &                       NINCW,BUFC,OIC,OJC,IVALID)
C
      IMPLICIT INTEGER (A-Z)
 
      INTEGER IRREP,LADC,NP,NITER,NEIGENV,IVX(2*NP+1)
      INTEGER INTBUF,IO_OFF,IO_DIA,IO_LAN
      REAL*8 P(LADC,2*NP+1), T(NP+1,NITER), R(LADC), DIAG(LADC)
      REAL*8 Z(NITER,NITER), AEV(LADC,NEIGENV), BUF(INTBUF)
      REAL*8 CEVL(NITER)
      INTEGER OI(INTBUF),OJ(INTBUF)
      LOGICAL CFLAG
      INTEGER NINCW
      REAL*8 BUFC(NINCW)
      INTEGER OIC(NINCW),OJC(NINCW),IVALID
      REAL*8 CPUTIM1,CPUTIM2,timer1,timer2,timertot
C
C---------------Description--------------------------------------------
C
C  Perform the real Lanczos iteration, generate matrix T_ij
C  the arrays BUF,OI,OJ are buffers for the adc matrix reading process
C  the file handles are IO_DIA and IO_OFF
C  If only eigenvalues are to be calculated, resp. in the first call of
C  this routine the Z-array is empty and will not be referenced !
C  In an eigenvector run Z contains the eigenvectors of the band matrix
C  diagonalization and is multiplied with the Lanczos vectors in order
C  to obtain the approximative eigenvectors of the original problem.
C  the variable .flag. determines if the additional eigenvector run
C  should be performed which is equivalent to the presence of the Z
C  vectors.
C  We write out the Lanczos matrix to file because it will be needed
C  for possible eigenvector calculations.
C  If an eigenvalue calculation is asked for we have to read in the
C  matrix of purged eigenvectors from ZMAT.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables and function definitions ------------
C
      REAL*8 Y,DDOT,Z1
      CHARACTER*10 LEVECFN
      real*8 rdummy
C
C---------------Executable code--------------------------------------
c
c  create symmetry-specific name for the long eigenvector file.
c  The eigenvector file of the Lanczos matrix needed for the
c  creation of the long eigenvectors is generically stored.
c
      IF(IRREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',IRREP
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',IRREP
      ENDIF

      IF(CFLAG) THEN
        CALL PST('Entering real IC-Lanczos for Eigenvectors+')
      ELSE
        CALL PST('Entering real IC-Lanczos for Eigenvalues+')
      ENDIF

      WRITE(IW,'(A40,I12)') ' Symmetry:',IRREP
      WRITE(IW,'(A40,I12)') ' Matrix dimension:',LADC
      WRITE(IW,'(A40,I12)') ' Band width:',NP
      WRITE(IW,'(A40,I12)') ' Number of iterations:',NITER
      WRITE(IW,'(A40,I12)') ' Number of words for incore:',NINCW
      WRITE(IW,'(A40,I12)') ' Number of valid elements incore:',IVALID
      IF(CFLAG) THEN
        WRITE(IW,'(A30,I7)') ' Requested eigenvectors:',NEIGENV
c       WRITE(IW,'(A30,A)') ' Long eigenvectors in:',LEVECFN
      ENDIF
      WRITE(IW,*)

      IF(IO_DIA.ne.-255)
     &  CALL QUIT('No incore version for old-format matrices!')
c
c  ... clear large eigenvector array if requested
c
      IF(CFLAG) THEN
        DO I=1,NEIGENV
          DO J=1,LADC
            AEV(J,I)=0.0D0
          ENDDO
        ENDDO
      ENDIF
c
c  ... read in purged eigenvector matrix TMATEVC for calculation of
c      the long eigenvectors. These are folded into the Lanczos
c      iterations. The third available file handle is used for that.
c      After the read process the handle can be reused and file is
c      deleted because this is the last position in the code where
c      it is needed.
c
      IF(CFLAG) THEN
        OPEN(IO_LAN,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IO_LAN)
        READ(IO_LAN) N1
        READ(IO_LAN) N2
        WRITE(IW,*) 'Found',N1,' Lanczos eigenvectors of length',
     &              N2,' on TMATEVC'
        IF(NEIGENV.gt.N1) THEN
          WRITE(*,*) 'More EVCs requested than available.'
          WRITE(*,*) 'Adjusting to the value given above.'
          NEIGENV = N1
        ENDIF
        IF(N2.NE.NITER)
     &    CALL QUIT('Incompatible T matrix eigenvector lengths!')
        DO I=1,NEIGENV
          READ(IO_LAN) CEVL(I)
          READ(IO_LAN) (Z(IX,I),IX=1,N2)
        ENDDO
        CLOSE(IO_LAN,STATUS='DELETE')
      ENDIF

c
c  ... clear Lanczos matrix
c
      DO I=1,NITER
        DO J=1,NP+1
          T(J,I)=0.0D0
        ENDDO
      ENDDO
c
c  ... clear r/p-vector array and initialize cartesian start vectors
c      start vectors are orthonormal and sit in the fields
c      np+1,np+2,...,2*np
c
      DO I=1,2*NP+1
        DO J=1,LADC
          P(J,I)=0.0D0
        ENDDO
      ENDDO
      DO J=1,LADC
        R(J)=0.0D0
      ENDDO
      DO I=1,NP
        P(I,I+NP)=1.0D0
      ENDDO

c
c  ... initialize vector mapping array
c
      DO I=1,2*NP+1
        IVX(I)=I
      ENDDO
c
c  ... and time count
c
      timertot=0.0d0
      CALL CPUUSED(CPUTIM1)

c----------------------------- lanczos loop start ----------------------

      DO 99 IL = 1,NITER
c
c
      IF(mod(il,100).eq.0) THEN
        WRITE(IW,*) '---  Band Lanczos iteration ',IL,' of',niter
        CALL CPUUSED(CPUTIM2)
        timer1 = cputim2 - cputim1
        write(*,*) 'Time per 100 iterations:',timer1
        timertot = timertot + timer1
        cputim1 = cputim2
      ENDIF
c
c  form r = Hq_i - sum t_ji*q_j
c  the j index in the sum runs from i-p up to i-1
c  and the loop is only executed if previous vectors exist.
c  the actual vector to be multiplied by A always sits in
c  P(*,np+1)
c
c
c  since for the double ionization the diagonal elements are
c  not stored separately, we define an alternative matrix-vector
c  multiplication routine which directly multiplies the matrix.
c  If IO_DIA has value of -255 (which is surely no valid file handle)
c  we know that the alternative routine should be accessed.
c
c  att: we need the buffer number in RMATBLK_ND! In order to avoid a lot
c  of parameter rewriting we communicate this via a common block.
c  We transfer the huge memory chunk to the incore Mat/Vec routine.
c
        CALL RMATBLK_IC(LADC, P(1,IVX(NP+1)), R, BUF,
     &                  OI, OJ, INTBUF, IO_OFF,
     &                  NINCW,BUFC,OIC,OJC,IVALID)

c  subtract all the previous components

        ICC = 1
        DO IP = MAX0(IL-NP,1), IL-1
          ITROW = ICC+1
          ITCOL = IL - ICC
          DO J=1,LADC
            R(J) = R(J)-T(ITROW,ITCOL)*P(J,IVX(NP-ICC+1))
          ENDDO
          ICC = ICC + 1
        ENDDO
c
c  form the t_ji = q_j^+ * r and r=r-t_ji*q_j for all
c  j=i...i+p-1
c
        DO IP = 1,NP
          Y = DDOT(LADC,P(1,IVX(NP+IP)),1,R,1)
          T(IP,IL) = Y
          DO J=1,LADC
            R(J) = R(J) - Y*P(J,IVX(NP+IP))
          ENDDO
        ENDDO
c
c  now r contains the whole RHS of the band-lanczos equation
c  and will be normalized
c  form t_i+p,i = |r|, q_i+p,i = r/|r|
c
        Y = DDOT(LADC,R,1,R,1)
        Y = DSQRT(Y)
        T(NP+1,IL) = Y

        Y = 1.0d0/Y
        DO J=1,LADC
          P(J,IVX(2*NP+1)) = R(J) * Y
        ENDDO
c
c  if applicable form the approximate eigenvector
c  with the first Lanczos vector Q_1
c  Important: the ith ready Lanczos vector is available at
c  P(*,IVX(NP+1)) exactly at this position of the
c  iteration. In case of an eigenvector run its contributions
c  to ALL eigenvectors have to be taken into account. After
c  cycling the contributions of the i+1. lanczos vector will
C  be folded into the final eigenvectors
c
        IF(CFLAG) THEN
          DO IOM = 1,NEIGENV
            Z1 = Z(IL,IOM)
            DO ICOL = 1,LADC
              AEV(ICOL,IOM) = AEV(ICOL,IOM) + Z1*P(ICOL,IVX(NP+1))
            ENDDO
          ENDDO
        ENDIF
c
c  now cycle the mapping array IVX
c  and start new main iteration
c
        CALL IVXROT(IVX,2*NP+1)

 99   CONTINUE
c
c  write out long ADC eigenvectors when eigenvector run finished
c
      IF(CFLAG) THEN
        OPEN(IO_LAN, FILE=LEVECFN, FORM='UNFORMATTED', 
     &       ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
        WRITE(IO_LAN) LADC,NEIGENV
        DO J=1,NEIGENV
          WRITE(IO_LAN) CEVL(J)
          WRITE(IO_LAN) (AEV(I,J),I=1,LADC)
        ENDDO
        CLOSE(IO_LAN)
        WRITE(IW,*) ' ADC LONG EIGENVECTORS WRITTEN.'
      ENDIF

      write(*,*) 'Total time spent in Matrix * Vector:',timertot
      write(*,*) 'Average time for Matrix * Vector:',
     &            timertot/dble(niter)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZLANCZOS_IC(IRREP, LADC, NP, NITER, NEIGENV, IVX,
     &                       P, T, R, DIAG, Z, AEV, BUF, OI, OJ, CEVL, 
     &                       CFLAG, INTBUF, IO_OFF, IO_DIA, IO_LAN,
     &                       NINCW,BUFC,OIC,OJC,IVALID)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER IRREP,LADC,NP,NITER,NEIGENV,IVX(2*NP+1)
      INTEGER INTBUF,IO_OFF,IO_DIA,IO_LAN
      COMPLEX*16 P(LADC,2*NP+1), T(NP+1,NITER), R(LADC), DIAG(LADC)
      COMPLEX*16 Z(NITER,NITER), AEV(LADC,NEIGENV), BUF(INTBUF)
      REAL*8 CEVL(NITER)
      INTEGER OI(INTBUF),OJ(INTBUF)
      LOGICAL CFLAG
      INTEGER NINCW
      COMPLEX*16 BUFC(NINCW)
      INTEGER OIC(NINCW),OJC(NINCW),IVALID
      REAL*8 CPUTIM1,CPUTIM2,timer1,timer2,timertot

C
C---------------Description--------------------------------------------
C
C  Perform the complex Lanczos iteration, generate matrix T_ij
C  the arrays BUF,OI,OJ are buffers for the adc matrix reading process
C  the file handles are IO_DIA and IO_OFF
C  If only eigenvalues are to be calculated, resp. in the first call of
C  this routine the Z-array is empty and will not be referenced !
C  In an eigenvector run Z contains the eigenvectors of the band matrix
C  diagonalization and is multiplied with the Lanczos vectors in order
C  to obtain the approximative eigenvectors of the original problem.
C  the variable .flag. determines if the additional eigenvector run
C  should be performed which is equivalent to the presence of the Z
C  vectors.
C  We write out the Lanczos matrix to file because it will be needed
C  for possible eigenvector calculations.
C  If an eigenvalue calculation is demanded we have to read in the
C  matrix of purged eigenvectors from ZMAT.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables and function definitions ------------
C
      COMPLEX*16 Y,ZDOTC,Z1
      CHARACTER*10 LEVECFN
      real*8 rdummy
C
C---------------Executable code--------------------------------------
C
      IF(IRREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',IRREP
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',IRREP
      ENDIF

      IF(CFLAG) THEN
        CALL PST('Entering complex IC Lanczos for Eigenvectors+')
      ELSE
        CALL PST('Entering complex IC Lanczos for Eigenvalues+')
      ENDIF

      WRITE(IW,'(A30,I7)') ' Symmetry:',IRREP
      WRITE(IW,'(A30,I7)') ' Matrix dimension:',LADC
      WRITE(IW,'(A30,I7)') ' Band width:',NP
      WRITE(IW,'(A30,I7)') ' Number of iterations:',NITER
      WRITE(IW,'(A40,I12)') ' Number of words for incore:',NINCW

      IF(CFLAG) THEN
        WRITE(IW,'(A30,I7)') ' Requested eigenvectors:',NEIGENV
      ENDIF

      IF(IO_DIA.ne.-255)
     &  CALL QUIT('No incore version for old-format matrices!')
c
c  ... clear approx eigenvector array if requested
c
      IF(CFLAG) THEN
        DO I=1,NEIGENV
          DO J=1,LADC
            AEV(J,I)=(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDIF
c
c  ... read in purged eigenvector matrix TMATEVC for calculation of
c      the long eigenvectors. These are folded into the Lanczos
c      iterations. The third available file handle is used for that.
c      After the read process the handle can be reused and file is
c      deleted because this is the last position in the code where
c      it is needed.
c
      IF(CFLAG) THEN
        OPEN(IO_LAN,FILe='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IO_LAN)
        READ(IO_LAN) N1
        READ(IO_LAN) N2
        WRITE(IW,*) 'Found',N1,' Lanczos eigenvectors of length',
     &              N2,' on TMATEVC'
        IF(NEIGENV.gt.N1) THEN
          WRITE(*,*) 'More EVCs requested than available.'
          WRITE(*,*) 'Adjusting to the value given above.'
          NEIGENV = N1
        ENDIF
        IF(N2.NE.NITER)
     &    CALL QUIT('Incompatible T matrix eigenvector lengths!')
        DO I=1,NEIGENV
          READ(IO_LAN) CEVL(I)
          READ(IO_LAN) (Z(IX,I),IX=1,N2)
        ENDDO
        CLOSE(IO_LAN,STATUS='DELETE')
      ENDIF
c
c  ... clear Lanczos matrix
c
      DO I=1,NITER
        DO J=1,NP+1
          T(J,I)=(0.0D0,0.0D0)
        ENDDO
      ENDDO
c
c  ... clear r/p-vector array and initialize start vectors
c      start vectors are orthonormal and sit in the fields
c      np+1,np+2,...,2*np
c
      DO I=1,2*NP+1
        DO J=1,LADC
          P(J,I)=(0.0D0,0.0D0)
        ENDDO
      ENDDO
      DO J=1,LADC
        R(J)=(0.0D0,0.0D0)
      ENDDO
      DO I=1,NP
        P(I,I+NP)=(1.0D0,0.0D0)
      ENDDO

c
c  ... initialize vector mapping array
c
      DO I=1,2*NP+1
        IVX(I)=I
      ENDDO
c
c  ... and time count
c
      timertot=0.0d0
      CALL CPUUSED(CPUTIM1)

c----------------------------- lanczos loop start ----------------------

      DO 99 IL = 1,NITER
c
c
      IF(mod(il,100).eq.0) THEN
        WRITE(IW,*) '---  Band Lanczos iteration ',IL,' of',niter
        CALL CPUUSED(CPUTIM2)
        timer1 = cputim2 - cputim1
        write(*,*) 'Time per 100 iterations:',timer1
        timertot = timertot + timer1
        cputim1 = cputim2
      ENDIF
c
c  form r = Hq_i - sum t_ji*q_j
c  the j index in the sum runs from i-p up to i-1
c  and the loop is only executed if previous vectors exist.
c  the actual vector to be multiplied by A always sits in
c  P(*,np+1)
c
      CALL CMATBLK_IC(LADC, P(1,IVX(NP+1)), R, BUF,
     &                OI, OJ, INTBUF, IO_OFF,
     &                NINCW,BUFC,OIC,OJC,IVALID)

c  subtract all the previous components

        ICC = 1
        DO IP = MAX0(IL-NP,1), IL-1
          ITROW = ICC+1
          ITCOL = IL - ICC
          DO J=1,LADC
            R(J) = R(J)-DCONJG(T(ITROW,ITCOL))*P(J,IVX(NP-ICC+1))
          ENDDO
          ICC = ICC + 1
        ENDDO
c
c  form the t_ji = q_j^+ * r and r=r-t_ji*q_j for all
c  j=i...i+p-1
c
        DO IP = 1,NP
          Y = ZDOTC(LADC,P(1,IVX(NP+IP)),1,R,1)
          T(IP,IL) = Y
          DO J=1,LADC
            R(J) = R(J) - Y*P(J,IVX(NP+IP))
          ENDDO
        ENDDO
c
c  now r contains the whole RHS of the band-lanczos equation
c  and will be normalized
c  form t_i+p,i = |r|, q_i+p,i = r/|r|
c
        Y = ZDOTC(LADC,R,1,R,1)
c
        Y = SQRT(Y)
        T(NP+1,IL) = Y

        Y = (1.0d0,0.0d0)/Y
        DO J=1,LADC
          P(J,IVX(2*NP+1)) = R(J) * Y
        ENDDO
c
c  if applicable form the approximate eigenvector
c  with the first Lanczos vector Q_1
c  Important: the ith ready Lanczos vector is available at
c  P(*,IVX(NP+1)) exactly at this position of the
c  iteration. In case of an eigenvector run its contributions
c  to ALL eigenvectors have to be taken into account. After
c  cycling the contributions of the i+1. lanczos vector will
C  be folded into the final eigenvectors
c
        IF(CFLAG) THEN
          DO IOM = 1,NEIGENV
            Z1 = Z(IL,IOM)
            DO ICOL = 1,LADC
              AEV(ICOL,IOM) = AEV(ICOL,IOM) + Z1*P(ICOL,IVX(NP+1))
            ENDDO
          ENDDO
        ENDIF
c
c  now cycle the mapping array IVX
c  and start new main iteration
c
        CALL IVXROT(IVX,2*NP+1)

 99   CONTINUE
c
c  write out long ADC eigenvectors when eigenvector run finished
c
      IF(CFLAG) THEN
        OPEN(IO_LAN, FILE=LEVECFN, FORM='UNFORMATTED',
     &       ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
        WRITE(IO_LAN) LADC,NEIGENV
        DO J=1,NEIGENV
          WRITE(IO_LAN) CEVL(J)
          WRITE(IO_LAN) (AEV(I,J),I=1,LADC)
        ENDDO
        CLOSE(IO_LAN)
        WRITE(IW,*) ' ADC LONG EIGENVECTORS WRITTEN.'
      ENDIF

      write(*,*) 'Total time spent in Matrix * Vector:',timertot
      write(*,*) 'Average time for Matrix * Vector:',
     &            timertot/dble(niter)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RBLANA(IW,IRREP,N,NP,NHOLES,IONIZLEVEL,
     &                  D,ZLONG,Z,T,IND,MLT,LSP,
     &                  IOUNIT,IOUNITC,ELIM,NMSPEC)

      IMPLICIT INTEGER (A-Z)

      INTEGER IW,IRREP,N,NP,NHOLES,IONIZLEVEL
      REAL*8 D(N)
      REAL*8 ZLONG(N,N),Z(2*NP,N),T(NP+1,N)
      INTEGER IND(N),MLT(N),LSP(N)
      INTEGER IOUNIT,IOUNITC
      REAL*8 ELIM
      CHARACTER*5 NMSPEC
C
C---------------Description--------------------------------------------
C
C  Analyses the eigenvalues/vectors of the Band matrix T and
C  treats multiple/spurious eigenvalues. All this for the
C  current main program symmetry stored in IRREP
c  (H.-D. Meyer's and Thomas Sommerfeld's formulation)
c
c  IOUNIT is a free handle for creating the XSPEC files
c  IOUNITC is an occupied handle referring to the config file.
c
c  eigenvalues come in D(*), complete eigenvectors in ZLONG(N,N) !
c  the first and last NP components will be needed for the analysis 
c  only. in order to use the unmodified code we copy those components
c  to the short vector array z(2*np,n)
c
c  the number N corresponds to the number of L-iterations == NITER
c
c  *** attention *** in the analysis some components of z are added
c  together or multiplied by some scaling factor stemming from the
c  multiplicity of some eigenvales. These modifications should also be applied
c  to the LONG eigenvector as well in order to get the right results.
c  If approx. eigenvectors are calculated the LONG eigenvector array
c  will be used because we need EVERY component of Z not only the
c  first and last NP ones. This causes trouble because according to
c  Dieter's analysis only the hole parts (the first NP ones) are weighted
c  and for the other components no strict modification scheme
c  is given. This means that we should determine the corrected long
c  Z vector by INVERSE INTERATION. There no artefacts will occur because
c  multiple and spurious eigenvalues are artefacts from LANCZOS and not
c  from the physical problem.
c
c
c  the band-matrix T is needed for the error analysis
c  according to Parlett: "The symmetric eigenvalue problem",
c  Chapter 13, sect 13-10-2 where a formula for the accuracy
c  of a Rayleigh-Ritz pair is given for the case of band-Lanczos.
C  From that one sees that the unused lower right triangle of T
C  together with the last NP components of the calculated eigenvectors
c  are needed in order to estimate the quality and convergence of 
c  the corresponding eigenvector.
c
c  the variable cntev actually counts the printed eigenvalues.
c  This is of use because each eigenvalue sequence then gets a contiguous
c  numbering. Numbering according to threshold will have gaps.
C
      REAL*8 AA,AK,A1K,SUM,SUP
      REAL*8 AUEV,THRS,OVL
      REAL*8 X1,Y1

      REAL*8 XM,DDOT
      DATA AUEV/27.2113957D0/

      CHARACTER*8 SPECFILE
      REAL*8 THTHR,WTHR
      CHARACTER*44 FIELD
C
C---------------Executable code--------------------------------------
C


      CALL PST('Processing the spectrum (real)+')
      WTHR=0.01D0
      WRITE(IW,*) 'Upper limit for printing eigenvalues on screen:',ELIM
      WRITE(IW,*) 'Threshold for considering eigenvectors:',WTHR
      WRITE(IW,*) 'Number of iterations:',N
      WRITE(IW,*) 'Band width:',NP

c
c  copy the first and last NP components of the long eigenvectors
c  to Z in order to apply dieters code easily.
c
      do i=1,n
        IND(i)=0
        MLT(i)=0
        LSP(i)=0
      enddo

c  prepare short array z
         
      do i=1,n
      do j=1,np
        z(j,i)=zlong(j,i)
        z(j+np,i)=zlong(n-np+j,i)
      enddo
      enddo

C     ----------------------------------------------------------
C     DETERMINE THE INDEX-VECTOR 'IND' IN ORDER TO REMOVE
C     MULTIPLE AND SPURIOUS EIGENVALUES.
C     ----------------------------------------------------------

      NP1 = NP + 1
      NP2 = 2 * NP

      NSP = 0
      NMU = 0
      ISP = 0
      JND = 0
      J   = 0
C
  50  J   = J + 1
      IF( J .GT. n )                               GOTO  40

      AA=DABS(DDOT(NP,Z(1,J),1,Z(1,J),1))

      NSP = NSP + 1
      ISP = ISP + 1
C.....TEST OF SPURIOUS EIGENVALUE.
      IF( AA .LT. 1.D-14 )                     GOTO  50
      NSP = NSP - 1
      ISP = ISP - 1
      JND = JND + 1
      IND(JND) = J
      MLT(JND) = 0
      LSP(JND) = ISP
      ISP      = 0
C
      IF( J .EQ. n )                               GOTO  40
      THRS = 5.D-7*( 1.D0+DABS(d(J)) )
C.....TEST OF MULTIPLE EIGENVALUE.
      IF( DABS( d(J) - d(J+1) ) .GT. THRS )           GOTO  50
      SUM = AA
      SUP = AA
      KK  = 0
      IA  = 0
C
      DO  41  K = 1, n - J
      OVL = 1.D-5
      IF( DABS( d(J+IA) - d(J+K) ) .GT. THRS )        GOTO  42

      AK=DABS(DDOT(NP,Z(1,J+K),1,Z(1,J+K),1))

      IF( AK .LE. 1.D-14 )  THEN
         KK     = K
         GOTO  41
      ELSE IF( AK/SUP + SUP/AK  .LT.  2.D+6 )    THEN
           OVL = 1.D-7
           IF( DABS(d(J+IA)-d(J+K)) .GT. 1.D-2*THRS )   GOTO  42
      END IF

      A1K=dabs(ddot(NP,Z(1,J+K),1,Z(1,J+IA),1))

      A1K = A1K*A1K/(SUP*AK)
      IF( DABS( A1K-1.D0 ) .LE. OVL )            THEN
         SUM = SUM + AK
         IF( AK .GT. SUP )  THEN
           SUP = AK
           IA  =  K
         END IF
         KK = K
      ELSE
        GOTO  42
      END IF
  41  CONTINUE
C
  42  d(J) = d(J+IA)
      DO  43  I = 1, NP
         Z(   I,J) = Z(   I,J+IA)*DSQRT(SUM/SUP)
         Z(NP+I,J) = Z(NP+I,J+IA)
  43  CONTINUE
C
      NMU = NMU + KK
      J   =  J  + KK
      MLT(JND) =  KK
      GOTO  50
C
  40  CONTINUE
C
C.....REMOVE SPURIOUS AND MULTIPLE EIGENVALUES.
C  ATT!
C  The complete eigenvectors also have to be purged!
C  In the Z-array only the upper and lower NP components
C  are available!
C
      DO  48  J = 1, JND
         I = IND(J)
         d(J) = d(I)
         DO  49  K = 1, NP2
            Z(K,J) = Z(K,I)
 49      continue
         DO  67  K = 1, N
            ZLONG(K,J) = ZLONG(K,I)
 67      continue
 48   CONTINUE

      write (IW,'(/,A)') '  Information on purging:'
      PRINT 905, JND, N, int(100.0d0*dble(jnd)/dble(n))
 905  FORMAT(4x,'There are',I5,' relevant residues out of',I5, 
     $     '   (',I3,'%)')
      write(6,'(4x,A,I4,6x,A,I4)') 'N_spurios = ', NSP,
     $     'N_multiple = ', NMU
c     write(6,*) 'Threshold for vector printing:',WTHR
C
C.....TRANSFORM THE EIGENVECTORS(FIRST NP ROWS) TO OBTAIN THE RESIDUES.
C.....STORE THE RESULT IN THE FIRST NP ROWS OF THE EIGENVECTOR MATRIX Z.
C.....MULTIPLY :  EIGENVECTORS(LAST ROWS) * 'ERROR-TRIANGLE'
C.....STORE THIS RESULT IN THE LAST ROWS OF THE EIGENVECTOR MATRIX Z.
C
      DO  30  J  = 1, JND
       DO  35  I  = 1, NP
         XM = 0.0D0
         DO  36  K = I, NP
           KK = K + NP
           IC = K + n - NP
           JC = I - K + NP1
           XM = XM + Z(KK,J)*T(JC,IC)
  36     CONTINUE
  35   Z(I+NP,J) = XM
  30  CONTINUE     
c
c    * * * * * *     S P E C T R U M     * * * * * * * * * * * *
c  spectrum is purged and available: all eigenvalues with corresponding
c  intensities obtained as squares of the hole/hole components
c  (pole strengths). IMPORTANT: the band width of the Lanczos matrix
c  can be larger than the number of actual hole states due to stability
c  reasons. If we sum the intensities we just look at the hole states !
c
c  print eigenvalues and accuracy in a more convenient format
c  threshold applies only to the screen output. pole strength = square of
c  eigenvector component. Write all eigenvalues and pole strength
c  to file named in SPECFILE
c
      if(irrep.gt.9) then
        write(SPECFILE,'(a5,a1,i2)') NMSPEC,'.',irrep
      else
        write(SPECFILE,'(a5,a2,i1)') NMSPEC,'.0',irrep
      endif
      OPEN(IOUNIT,FILE=SPECFILE,FORM='FORMATTED',STATUS='UNKNOWN')

      cntev=0
      DO J = 1, JND 
c
c  pole strength calculation
c
        x1 = 0.0d0
        do i = 1, nholes
           x1 = x1 + z(i,j)*z(i,j)
        enddo
c
c  error analysis
c
        y1 = 0.0d0
        do i = 1, np
          y1 = y1 + z(np+i,j)*z(np+i,j)
        enddo
        y1 = dsqrt(y1)
c
c  everything goes to file, no threshold for pole strength and
c  no upper limit for the eigenvalues. That means: all eigenvalues
c  together with their pole strengths and error estimates are written
c  to the XSPEC.XX files.
c
        IF(IONIZLEVEL.EQ.1) THEN
c
c  take care of single ionizations
c
           isipcount=0
           DO iths=1,nholes
             ththr=z(iths,j)*z(iths,j)
             if(ththr.gt.0.001) isipcount = isipcount + 1
           ENDDO 
           write(iounit,'(2x,F14.8,4x,2F12.6,4x,A1,I4,A,I3)')
     &        d(j)*AUEV,x1,y1,'@',isipcount,'  File: ',irrep
           DO iths=1,nholes
             ththr=z(iths,j)*z(iths,j)
             if(ththr.gt.0.001) then
               read(iounitc,rec=iths) field
               write(iounit,'(10X,I5,A3,F12.5,1X,"(",F7.5,")",A8,A22)')
     &         iths,':  ',z(iths,j),z(iths,j)*z(iths,j),'   conf:',field
             endif
           ENDDO

        ELSE
c
c  take care of double ionizations and excitations
c
           idipcount=0
           DO iths=1,nholes
             ththr=z(iths,j)*z(iths,j)
             if(ththr.gt.0.001) idipcount = idipcount + 1
           ENDDO 
           write(iounit,'(2x,F14.8,4x,2F12.6,4x,A1,I4,A,I3)') 
     &        d(j)*AUEV,x1,y1,'@',idipcount,'  File: ',irrep
           DO iths=1,nholes
             ththr=z(iths,j)*z(iths,j)
             if(ththr.gt.0.001) then
               read(iounitc,rec=iths) field
               write(iounit,'(10X,I5,A3,F12.5,1X,"(",F7.5,")",A8,A44)') 
     &         iths,':  ',z(iths,j),z(iths,j)*z(iths,j),'  conf: ',field
             endif
           ENDDO          
        ENDIF
c
c  print on screen only those eigenvalues that are below the given limit.
c  use cntev for contiguous numbering !
c
c  ATT !!! For the two-hole final states an output of the squared
c  coefficient is shrouding information! We therefore write out the
c  coefficients themselves!
c
        if(d(j)*AUEV.le.ELIM) THEN
          cntev = cntev + 1
          write(iw,388)
 388      format(//2x,60('-')/'  no.    eigenvalue (eV)   mean-res.',
     &    3x,'error estimate'/2x,60('-')/)
          write(iw,389) cntev,'@', d(j)*AUEV, x1, y1, ind(j)
 389      format(I5,A1,2x,F14.8,4x,F8.4,4x,1PE11.2,'  (',I5,')')
          IF(IONIZLEVEL.EQ.2) THEN
            write(iw,*) '       -----Contributing 2h configurations:'
            DO iths=1,nholes
              ththr=z(iths,j)*z(iths,j)
              if(ththr.gt.0.0001) then
                read(iounitc,rec=iths) field
                write(iw,'(12X,I5,A3,F12.5,A22)') 
     &             iths,':  ',z(iths,j),field
              endif
            ENDDO          
          ELSE
            write(iw,'(100I8)') (ix,ix=1,nholes)
            write(iw,'(100F8.4)') 
     &           (z(ix,j)*z(ix,j),ix=1,nholes)
          ENDIF
        endif
      ENDDO

      CLOSE(IOUNIT)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZBLANA(IW,IRREP,N,NP,NHOLES,IONIZLEVEL,
     &                  D,ZLONG,Z,T,IND,MLT,LSP,
     &                  IOUNIT,IOUNITC,ELIM,NMSPEC)

      IMPLICIT INTEGER (A-Z)

      INTEGER IW,IRREP,N,NP,NHOLES,IONIZLEVEL
      REAL*8 D(N)
      COMPLEX*16 ZLONG(N,N),Z(2*NP,N),T(NP+1,N)
      INTEGER IND(N),MLT(N),LSP(N)
      INTEGER IOUNIT,IOUNITC
      REAL*8 ELIM
      CHARACTER*5 NMSPEC
C
C---------------Description--------------------------------------------
C
C  Analyses the eigenvalues/vectors of the Band matrix T and
C  treats multiple/spurious eigenvalues. All this for the
C  current main program symmetry stored in IRREP
c  (H.-D. Meyer's and Thomas Sommerfeld's formulation)
c
c  IOUNIT is a free handle for creating the XSPEC files
c  IOUNITC is an occupied handle referring to the config file.
c
c  eigenvalues come in D(*), complete eigenvectors in ZLONG(N,N) !
c  the first and last NP components will be needed for the analysis 
c  only. in order to use the unmodified code we copy those components
c  to the short vector array z(2*np,n)
c
c  the number N corresponds to the number of L-iterations == NITER
c
c  *** attention *** in the analysis some components of z are added
c  together or multiplied by some scaling factor stemming from the
c  multiplicity of some eigenvales. These modifications should also be applied
c  to the LONG eigenvector as well in order to get the right results.
c  If approx. eigenvectors are calculated the LONG eigenvector array
c  will be used because we need EVERY component of Z not only the
c  first and last NP ones. This causes trouble because according to
c  Dieter's analysis only the hole parts (the first NP ones) are weighted
c  and for the other components no strict modification scheme
c  is given. This means that we should determine the corrected long
c  Z vector by INVERSE INTERATION. There no artefacts will occur because
c  multiple and spurious eigenvalues are artefacts from LANCZOS and not
c  from the physical problem.
c
c
c  the band-matrix T is needed for the error analysis
c  according to Parlett: "The symmetric eigenvalue problem",
c  Chapter 13, sect 13-10-2 where a formula for the accuracy
c  of a Rayleigh-Ritz pair is given for the case of band-Lanczos.
C  From that one sees that the unused lower right triangle of T
C  together with the last NP components of the calculated eigenvectors
c  are needed in order to estimate the quality and convergence of 
c  the corresponding eigenvector.
c
c  the variable cntev actually counts the printed eigenvalues.
c  This is of use because each eigenvalue sequence then gets a contiguous
c  numbering. Numbering according to threshold will have gaps.
C
C---------------Local variables--------------------------------------
C
      REAL*8 AA,AK,A1K,SUM,SUP
      REAL*8 AUEV,THRS,OVL
      REAL*8 X1,Y1

      COMPLEX*16 XM,ZDOTC
      DATA AUEV/27.2113957D0/

      CHARACTER*8 SPECFILE
c variables for hole-configuration analysis
      REAL*8 THTHR,WTHR
      CHARACTER*44 FIELD
C
C---------------Executable code--------------------------------------
C


      CALL PST('Processing the spectrum (complex)+')
      WTHR=0.01D0
      WRITE(IW,*) 'Upper limit for printing eigenvalues on screen',ELIM
      WRITE(IW,*) 'Threshold pole strength for eigenvectors:',WTHR
      WRITE(IW,*) 'Number of iterations:',N
      WRITE(IW,*) 'Band width:',NP
c
c  copy the first and last NP components of the long eigenvectors
c  to Z in order to apply dieters code easily.
c
      do i=1,n
        IND(i)=0
        MLT(i)=0
        LSP(i)=0
      enddo

c  prepare short array z
         
      do i=1,n
      do j=1,np
        z(j,i)=zlong(j,i)
        z(j+np,i)=zlong(n-np+j,i)
      enddo
      enddo

C     ----------------------------------------------------------
C     DETERMINE THE INDEX-VECTOR 'IND' IN ORDER TO REMOVE
C     MULTIPLE AND SPURIOUS EIGENVALUES.
C     ----------------------------------------------------------

      NP1 = NP + 1
      NP2 = 2 * NP

      NSP = 0
      NMU = 0
      ISP = 0
      JND = 0
      J   = 0
C
  50  J   = J + 1
      IF( J .GT. n )                               GOTO  40

      AA=ABS(ZDOTC(NP,Z(1,J),1,Z(1,J),1))

      NSP = NSP + 1
      ISP = ISP + 1
C.....TEST OF SPURIOUS EIGENVALUE.
      IF( AA .LT. 1.D-14 )                     GOTO  50
      NSP = NSP - 1
      ISP = ISP - 1
      JND = JND + 1
      IND(JND) = J
      MLT(JND) = 0
      LSP(JND) = ISP
      ISP      = 0
C
      IF( J .EQ. n )                               GOTO  40
      THRS = 5.D-7*( 1.D0+DABS(d(J)) )
C.....TEST OF MULTIPLE EIGENVALUE.
      IF( DABS( d(J) - d(J+1) ) .GT. THRS )           GOTO  50
      SUM = AA
      SUP = AA
      KK  = 0
      IA  = 0
C
      DO  41  K = 1, n - J
      OVL = 1.D-5
      IF( DABS( d(J+IA) - d(J+K) ) .GT. THRS )        GOTO  42

      AK=ABS(ZDOTC(NP,Z(1,J+K),1,Z(1,J+K),1))

      IF( AK .LE. 1.D-14 )  THEN
         KK     = K
         GOTO  41
      ELSE IF( AK/SUP + SUP/AK  .LT.  2.D+6 )    THEN
           OVL = 1.D-7
           IF( DABS(d(J+IA)-d(J+K)) .GT. 1.D-2*THRS )   GOTO  42
      END IF

      A1K=abs(zdotc(NP,Z(1,J+K),1,Z(1,J+IA),1))

      A1K = A1K*A1K/(SUP*AK)
      IF( DABS( A1K-1.D0 ) .LE. OVL )            THEN
         SUM = SUM + AK
         IF( AK .GT. SUP )  THEN
           SUP = AK
           IA  =  K
         END IF
         KK = K
      ELSE
        GOTO  42
      END IF
  41  CONTINUE
C
  42  d(J) = d(J+IA)
      DO  43  I = 1, NP
         Z(   I,J) = Z(   I,J+IA)*DCMPLX(DSQRT(SUM/SUP),0.0d0)
         Z(NP+I,J) = Z(NP+I,J+IA)
  43  CONTINUE
C
      NMU = NMU + KK
      J   =  J  + KK
      MLT(JND) =  KK
      GOTO  50
C
  40  CONTINUE
C
C.....REMOVE SPURIOUS AND MULTIPLE EIGENVALUES.
C  ATT!
C  The complete eigenvectors also have to be purged!
C  In the Z-array only the upper and lower NP components
C  are available!
C
      DO  48  J = 1, JND
         I = IND(J)
         d(J) = d(I)
         DO  49  K = 1, NP2
            Z(K,J) = Z(K,I)
 49      CONTINUE
         DO  67  K = 1, N
            ZLONG(K,J) = ZLONG(K,I)
 67      CONTINUE
 48   CONTINUE

      write (IW,'(/,A)') '  Information on purging:'
      PRINT 905, JND, N, int(100.0d0*dble(jnd)/dble(n))
 905  FORMAT(4x,'There are',I5,' relevant residues out of',I5, 
     $     '   (',I3,'%)')
      write(IW,'(4x,A,I4,6x,A,I4)') 'N_spurios = ', NSP,
     $     'N_multiple = ', NMU
c     write(IW,*) 'Threshold for vector printing:',WTHR
C
C
C.....TRANSFORM THE EIGENVECTORS(FIRST NP ROWS) TO OBTAIN THE RESIDUES.
C.....STORE THE RESULT IN THE FIRST NP ROWS OF THE EIGENVECTOR MATRIX Z.
C.....MULTIPLY :  EIGENVECTORS(LAST ROWS) * 'ERROR-TRIANGLE'
C.....STORE THIS RESULT IN THE LAST ROWS OF THE EIGENVECTOR MATRIX Z.
C
      DO  30  J  = 1, JND
       DO  35  I  = 1, NP
         XM = (0.0D0,0.0D0)
         DO  36  K = I, NP
           KK = K + NP
           IC = K + n - NP
           JC = I - K + NP1
           XM = XM + Z(KK,J)*T(JC,IC)
  36     CONTINUE
  35   Z(I+NP,J) = XM
  30  CONTINUE     

c
c    * * * * * *     S P E C T R U M     * * * * * * * * * * * *
c  spectrum is purged and available: all eigenvalues with corresponding
c  intensities obtained as squares of the hole/hole components
c  (pole strengths). IMPORTANT: the band width of the Lanczos matrix
c  can be larger than the number of actual hole states due to stability
c  reasons. If we sum the intensities we just look at the hole states !
c
c
c  print eigenvalues and accuracy in a more convenient format
c  threshold applies only to the screen output. pole strength = square of
c  eigenvector component. Write all eigenvalues and pole strength
c  to file named in SPECFILE
c
      if(irrep.gt.9) then
        write(SPECFILE,'(a5,a1,i2)') NMSPEC,'.',irrep
      else
        write(SPECFILE,'(a5,a2,i1)') NMSPEC,'.0',irrep
      endif
      open(iounit,file=SPECFILE,form='FORMATTED',status='UNKNOWN')
      
      cntev=0
      DO J = 1, JND 
c
c  pole strength calculation
c
        x1 = 0.0d0
        do i = 1, nholes
           x1 = x1 + DBLE( z(i,j)*dconjg(z(i,j)) )
        enddo
c
c  error analysis
c
        y1 = 0.0d0
        do i = 1, np
          y1 = y1 + DBLE( z(np+i,j)*dconjg(z(np+i,j)) )
        enddo
        y1 = dsqrt(y1)
c
c  everything goes to file, no threshold for pole strength and
c  no upper limit for the eigenvalues. That means: all eigenvalues
c  together with their pole strengths and error estimates are written
c  to the XSPEC.XX files.
c
        IF(IONIZLEVEL.EQ.1) THEN
c
c  take care of single ionizations
c
           isipcount=0
           DO iths=1,nholes
             ththr=z(iths,j)*dconjg(z(iths,j))
             if(ththr.gt.0.001) isipcount = isipcount + 1
           ENDDO 
           write(iounit,'(2x,F14.8,4x,2F12.6,4x,A1,I4,A,I3)')
     &        d(j)*AUEV,x1,y1,'@',isipcount,'  File: ',irrep
           DO iths=1,nholes
             ththr=dble(z(iths,j)*dconjg(z(iths,j)))
             if(ththr.gt.0.001) then
               read(iounitc,rec=iths) field
               write(iounit,'(10X,I5,A3,2F12.5,1X,"(",F7.5,")",A8,A22)')
     &         iths,':  ',z(iths,j),dble(z(iths,j)*dconjg(z(iths,j))),
     &         '   conf:',field
             endif
           ENDDO

        ELSE
c
c  take care of double ionizations
c
           idipcount=0
           DO iths=1,nholes
             ththr=z(iths,j)*dconjg(z(iths,j))
             if(ththr.gt.0.001) idipcount = idipcount + 1
           ENDDO 
           write(iounit,'(2x,F14.8,4x,2F12.6,4x,A1,I4,A,I3)')
     &        d(j)*AUEV,x1,y1,'@',idipcount,'  File: ',irrep
           DO iths=1,nholes
             ththr=dble(z(iths,j)*dconjg(z(iths,j)))
             if(ththr.gt.0.001) then
               read(iounitc,rec=iths) field
               write(iounit,'(10X,I5,A3,2F12.5,1X,"(",F7.5,")",A8,A44)')
     &            iths,':  ',z(iths,j),ththr,'  conf: ',field
             endif
           ENDDO
        ENDIF
c
c  print on screen only those eigenvalues that are below the given limit.
c  use cntev for contiguous numbering ! Again: we now write out the
c  coefficients themselves for the two-particle case!
c
        if(d(j)*AUEV.le.ELIM) then
          cntev = cntev + 1
          write(iw,388)
 388      format(//2x,60('-')/'  no.    eigenvalue (eV)   mean-res.',
     &    3x,'error estimate'/2x,60('-')/)
          write(iw,389) cntev,'@', d(j)*AUEV, x1, y1, ind(j)
 389      format(I5,A1,2x,F14.8,4x,F8.4,4x,1PE11.2,'  (',I5,')')

          IF(IONIZLEVEL.EQ.2) THEN
            write(iw,*) '       -----Contributing 2h configurations:'
            DO iths=1,nholes
              ththr=z(iths,j)*dconjg(z(iths,j))
              if(ththr.gt.0.0001) then
                read(iounitc,rec=iths) field
                write(iw,'(12X,I5,A3,2F12.5,A22)') 
     &             iths,':  ',z(iths,j),field
              endif
            ENDDO    
          ELSE
            write(iw,'(100I8)') (ix,ix=1,nholes)
            write(iw,'(100F8.4)') 
     &           (dble(z(ix,j)*dconjg(z(ix,j))),ix=1,nholes)
          ENDIF
        endif

      ENDDO

      CLOSE(IOUNIT)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RBANDDIA(IOUNIT,N,NP,LDAB,NEIGENV,AB,EVL,Z,RWORK)
C
      use adc_mat
      use memory_allocator
      use adc_fano_diag

      IMPLICIT INTEGER (A-Z)

      INTEGER IOUNIT,N,NP,LDAB
      INTEGER NEIGENV
      REAL*8 AB(LDAB,*),Z(N,*)
      REAL*8 RWORK(6*N),EVL(N)
      
C
C---------------Description--------------------------------------------
C
C  Perform the diagonalization of the real symmetric band matrix
C  by calling the appropriate LAPACK routine. This routine
C  serves as an interface in order to avoid too many command line
C  arguments.
C  The Lanczos matrix T from the caller comes in the AB array.
C     ****  ATTENTION ****  due to the usual programming conventions
C  of the LAPACK routines the original matrix T is DESTROYED !!!
C
C  important: we used to save the eigenvectors of the Lanc matrix
C  *after* purging. But when one analyzes the rebuilding of the long ADC
C  vectors inconsistencies may arise. We therefore store *all* evecs of
C  L and deal with multitudes separately.

C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER INFO,LDZ
      CHARACTER*1 JOBZ,UPLO
      REAL*8 Y
      REAL*8 AUEV
      INTEGER                            :: counter

      INTEGER, ALLOCATABLE, DIMENSION(:) :: positions
      INTEGER                            :: nr_vecs
      INTEGER                            :: j, k
      REAL                               :: upper, lower

C
C---------------Executable code--------------------------------------
C
C we always compute the eigenvectors of the band matrix, because
C they are needed in the further analysis
C These eigenvectors have the FULL length of the Lanczos matrix. 
C but not the length of the original ADC matrix!
C In the original formulation only
C the first and last NP components were available and the whole
C eigenvector was calculated by inverse iteration. We do not pursue this
C direction due to serious problems with the inverse interation
C algorithm
C
C N is the dimension of the band matrix (not the number of sidebands!)

      JOBZ='V'
      UPLO='L'
      LDZ=N

      CALL DSBEV(JOBZ, UPLO, N, NP, AB, LDAB, EVL, Z, LDZ,
     &           RWORK, INFO)

      WRITE(iw,*) 'info = ', info

      IF(INFO.NE.0) THEN
        CALL QUIT('Error in DSBEV')
      ENDIF

      WRITE(IW,*)
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*) 'First 100 RAW Lanczos EVAL/EVECS + intensity.'
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*)
      WRITE(IW,*) ' No.         (a.u.)      (eV)'


      AUEV=27.2113957D0

      DO I=1,MIN(100,N)
        Y = 0.0D0
        DO J=1,NP
          Y = Y + Z(J,I)*Z(J,I)
        ENDDO
        WRITE(IW,111) I,EVL(I),EVL(I)*AUEV,DABS(Y)
      ENDDO

 111  FORMAT(I5,2x,2F12.6,3x,'mean-res',2x,F10.4)
c
c eigenvectors of the Lanczos matrix are ready. Save them *all* to disk
c Since each final state symmetry is treated immediately after
c calculation, TMATVEC does not get an individual symmetry label.
c

      IF (reladc_md_isfano) THEN
        CALL alloc(positions,N,id="array of position numbers fano Lanc")
        CALL fano_find_in_vec(N,EVL,Z,positions,nr_vecs,iw)

        NEIGENV = nr_vecs
        OPEN(IOUNIT,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        WRITE(IOUNIT) nr_vecs   ! number of computed eigenvectors
        WRITE(IOUNIT) N   ! length of eigenvectors (= matrix dimension)
        DO counter=1, nr_vecs
          j = positions(counter)
          WRITE(IOUNIT) EVL(J) ! eigenvalue
          WRITE(IOUNIT) (Z(IX,J),IX=1,N)  !corresponding eigenvector
        ENDDO
        CLOSE(IOUNIT)
      
        WRITE(*,*) 'Fano set of Lanczos eigenvectors saved to TMATEVC'

        CALL dealloc(positions)

      ELSE
        CALL alloc(positions,N,id='chosen eigenvectors for analysis')
        nr_vecs = 0

        upper = reladc_md_eeigv_upper
        lower = reladc_md_eeigv_lower

        DO j = 1, N
          IF ((EVL(j).GE.lower).AND.(EVL(j).LE.upper)) THEN
            nr_vecs = nr_vecs + 1
            positions(nr_vecs) = j
          END IF
        END DO

        NEIGENV = nr_vecs

        OPEN(IOUNIT,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        WRITE(IOUNIT) nr_vecs   ! number of computed eigenvectors (all!)
        WRITE(IOUNIT) N   ! length of eigenvectors (= matrix dimension)
        DO k = 1, nr_vecs
          j = positions(k)
          WRITE(IOUNIT) EVL(J) ! eigenvalue
          WRITE(IOUNIT) (Z(IX,J),IX=1,N)  !corresponding eigenvector
        END DO
        CLOSE(IOUNIT)
    
        WRITE(iw,*) NEIGENV, 'Lanczos eigenvectors saved to TMATEVC'
        WRITE(iw,*) 'Energy range in a.u.: ', lower, '--', upper

        CALL dealloc(positions)
      END IF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZBANDDIA(IOUNIT,N,NP,LDAB,NEIGENV,AB,EVL,Z,WORK,RWORK)
C
      use adc_mat
      use memory_allocator
      use adc_fano_diag
!
      IMPLICIT INTEGER (A-Z)

      INTEGER N,NP,LDAB
      INTEGER NEIGENV
      COMPLEX*16 AB(LDAB,*),Z(N,*),WORK(4*N)
      DOUBLE PRECISION RWORK(6*N),EVL(N)
      
C
C---------------Description--------------------------------------------
C
C  Perform the diagonalization of the Hermitian band matrix
C  by calling the appropriate LAPACK routine. This routine
C  is interfacing in order to avoid too many command line
C  arguments.
C  The Lanczos matrix T from the caller comes in the AB array.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER INFO,LDZ
      CHARACTER*1 JOBZ,UPLO
      COMPLEX*16 Y
      REAL*8 AUEV

      INTEGER                            :: nr_vecs            
      INTEGER                            :: counter
      INTEGER, ALLOCATABLE, DIMENSION(:) :: positions
      INTEGER                            :: j, k
      REAL                               :: upper, lower

C
C---------------Executable code--------------------------------------
C
C we always compute the eigenvectors of the band matrix, because
C they are needed in the further analysis
C These eigenvectors have FULL length. In the original formulation only
C the first and last NP components were available and the whole
C eigenvector was calculated by inverse iteration.
C

      JOBZ='V'
      UPLO='L'
      LDZ=N

      CALL ZHBEV(JOBZ, UPLO, N, NP, AB, LDAB, EVL, Z, LDZ,
     &           WORK, RWORK, INFO)

      IF(INFO.NE.0) THEN
        CALL QUIT('Error in ZHBEV')
      ENDIF

      WRITE(IW,*)
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*) 'First 100 RAW Lanczos EVAL/EVECS + intensity.'
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*)
      WRITE(IW,*) ' No.         (a.u.)      (eV)'


      AUEV=27.2113957D0

      DO I=1,MIN(100,N)
        Y = (0.0D0,0.0D0)
        DO J=1,NP
          Y = Y + Z(J,I)*DCONJG(Z(J,I))
        ENDDO
        WRITE(IW,111) I,EVL(I),EVL(I)*AUEV,ABS(Y)
      ENDDO

 111  FORMAT(I5,2x,2F12.6,3x,'mean-res',2x,F10.4)
c
c eigenvectors of the Lanczos matrix are ready. Save the chosen to disk
c Since each final state symmetry is treated immediately after
c calculation, TMATVEC does not get an individual symmetry label.
c
                                                                  
      IF (reladc_md_isfano) THEN                                  
        CALL alloc(positions,N,id="array of position numbers fano Lanc")
        CALL cfano_find_in_vec(N,EVL,Z,positions,nr_vecs,iw)    
                                                                  
        NEIGENV = nr_vecs                                      
        OPEN(IOUNIT,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        WRITE(IOUNIT) nr_vecs   ! number of computed eigenvectors
        WRITE(IOUNIT) N   ! length of eigenvectors (= matrix dimension)
        DO counter=1, nr_vecs                                  
          j = positions(counter)                                  
          WRITE(IOUNIT) EVL(J) ! eigenvalue                       
          WRITE(IOUNIT) (Z(IX,J),IX=1,N)  !corresponding eigenvector
        ENDDO                                                     
        CLOSE(IOUNIT)                                             
                                                                  
        WRITE(*,*) 'Fano set of Lanczos eigenvectors saved to TMATEVC'
                                                                  
        CALL dealloc(positions)                                   
                                                                  
      ELSE                                                        
        CALL alloc(positions,N,id='chosen eigenvectors for analysis')
        nr_vecs = 0

        upper = reladc_md_eeigv_upper
        lower = reladc_md_eeigv_lower

        DO j = 1, N
          IF ((EVL(j).GE.lower).AND.(EVL(j).LE.upper)) THEN
            nr_vecs = nr_vecs + 1
            positions(nr_vecs) = j
          END IF
        END DO

        NEIGENV = nr_vecs

        OPEN(IOUNIT,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        WRITE(IOUNIT) nr_vecs   ! number of computed eigenvectors (all!)
        WRITE(IOUNIT) N   ! length of eigenvectors (= matrix dimension)
        DO k = 1, nr_vecs
          j = positions(k)
          WRITE(IOUNIT) EVL(J) ! eigenvalue
          WRITE(IOUNIT) (Z(IX,J),IX=1,N)  !corresponding eigenvector
        END DO
        CLOSE(IOUNIT)
    
        WRITE(iw,*) NEIGENV, 'Lanczos eigenvectors saved to TMATEVC'
        WRITE(iw,*) 'Energy range in a.u.: ', lower, '--', upper

        CALL dealloc(positions)
      END IF  


      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE IVXROT(IVX,N)
      IMPLICIT INTEGER (A-Z)
C
C  This routine cycles the mapping function according to a new
C  macro iteration IL in the calling program
C
      INTEGER IVX(N)
      INTEGER ITEMP

      ITEMP = IVX(1)
      DO I=1,N-1
        IVX(I)=IVX(I+1)
      ENDDO
      IVX(N)=ITEMP

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RMATBLK(LADC,PO,PN,DIAG,BUF,IOI,IOJ,INTBUF,
     &                   IO_DIA,IO_OFF)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LADC,INTBUF,IO_DIA,IO_OFF
      REAL*8 PO(LADC), PN(LADC), DIAG(LADC)
      REAL*8 BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)
C
C---------------Description--------------------------------------------
C
C  Perform the real symmetric multiplication of the ADC matrix
C  with the Krylov vectors. The ADC matrix is real symmetric
C  This has to be taken into account.
C
C---------------Local variables--------------------------------------
C
      INTEGER NBUFS,IREC,NACT,JDUMMY,IROW,ICOL
      REAL*8 A
C
C---------------Executable code--------------------------------------

C
C  the diagonal elements
C
      REWIND(IO_DIA)
      READ(IO_DIA,ERR=88,END=99) (DIAG(IXX),IXX=1,LADC),NBUFS

      DO J = 1, LADC
        PN(J) = PO(J) * DIAG(J)
      ENDDO

C  for the off-diagonal elements we note that only the lower triangular
C  part of the ADC-matrix is stored. the missing contribution is from
C  the symmetric upper triangular part !

c  .. read off-diagonal elements and multiply on the fly

      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)
            PN(IROW) = PN(IROW) + PO(ICOL) * A
            PN(ICOL) = PN(ICOL) + PO(IROW) * A
        ENDDO
      ENDDO

      RETURN

 88   CALL QUIT('Diagonal read error in RMATBLK')
 99   CALL QUIT('Unexpected end in RMATBLK')
 199  CALL QUIT('Off-diagonal read error in RMATBLK')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RMATBLK_ND(LADC,PO,PN,BUF,IOI,IOJ,INTBUF,IO_OFF)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LADC,INTBUF,IO_OFF
      REAL*8 PO(LADC), PN(LADC)
      REAL*8 BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)

C
C---------------Description--------------------------------------------
C
C  Perform the real symmetric multiplication of the ADC matrix
C  with the Krylov vectors. Here *NO* separate diagonal entries are
C  processed ! They are contained in the matrix already and
C  a single multiplication suffices.
c
c  att! the required number of buffers is communicated through this
c  common block.
C
C---------------Common Blocks--------------------------------------
C
      INTEGER NBUFS
      COMMON/MATCHUNKS/NBUFS
C
C---------------Local variables--------------------------------------
C
      INTEGER IREC,NACT,JDUMMY,IROW,ICOL
      REAL*8 A
C
C---------------Executable code--------------------------------------

C  for the matrix elements we note that only the lower triangular
C  part *AND* the diagonal of the ADC-matrix is stored. 
C  the missing contribution is from the symmetric upper triangular part !

c  the if takes care of not double counting the diagonal.
c  since the diagonal is now *contained* in the matrix
c
c  new vector has to be zeroed!
      DO J = 1, LADC
        PN(J) =  0.0d0
      ENDDO

c  .. read off-diagonal elements and multiply on the fly

      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
        ENDDO
      ENDDO

      RETURN

 889  CALL QUIT('Matrix read error in RMATBLK_ND')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RMATBLK_IC(LADC,PO,PN,BUF,IOI,IOJ,INTBUF,IO_OFF,
     &                      nincw,bufc,ioic,iojc,ivalid)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LADC,INTBUF,IO_OFF
      REAL*8  PO(LADC), PN(LADC)
      REAL*8  BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)
      integer nincw
      real*8  bufc(nincw)
      integer ioic(nincw),iojc(nincw),ivalid

C
C---------------Description--------------------------------------------
C
C  Perform the real symmetric multiplication of the ADC matrix
C  with the Krylov vectors. Here *NO* separate diagonal entries are
C  processed ! They are contained in the matrix already and
C  a single multiplication suffices.
c
c  att! we have the icore version here.  Matmul splits in an incore
c  part and an external part for those buffers that can not be held
c  in memory.
C
C---------------Common Blocks--------------------------------------
C
      INTEGER NBUFS
      COMMON/MATCHUNKS/NBUFS
C
C---------------Local variables--------------------------------------
C
      INTEGER IREC,NACT,JDUMMY,IROW,ICOL
      REAL*8 A
C
C---------------Executable code--------------------------------------

C  for the matrix elements we note that only the lower triangular
C  part *AND* the diagonal of the ADC-matrix is stored. 
C  the missing contribution is from the symmetric upper triangular part !

c  the if takes care of not double counting the diagonal.
c  since the diagonal is now *contained* in the matrix

c  new vector has to be zeroed!

      DO J = 1, LADC
        PN(J) =  0.0d0
      ENDDO
c
c  do incore part. In this part of the program there *IS* an
c  incore contribution. So we do not need to check this condition
c  any further.
c
      DO K = 1, IVALID
        IROW = ioic(K)
        ICOL = iojc(K)
        A    = bufc(K)
        PN(IROW) = PN(IROW) + PO(ICOL) * A
        IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
      ENDDO
c
c  next do out of core part (if there is any)
c  remark: if matrix is completely incore the following loop
c  will never be executed. (nbufs = 0).
c  Under the IO_OFF file handle we find the *truncated* ADC
c  matrix only containing the ooc buffers. (there is no lseek
c  for this type of sequential files)
c
      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
        ENDDO
      ENDDO

      RETURN

 889  CALL QUIT('Matrix read error in RMATBLK_IC')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMATBLK(LADC,PO,PN,DIAG,BUF,IOI,IOJ,INTBUF,
     &                   IO_DIA,IO_OFF)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LADC,INTBUF,IO_DIA,IO_OFF
      COMPLEX*16 PO(LADC), PN(LADC), DIAG(LADC)
      COMPLEX*16 BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)

C
C---------------Description--------------------------------------------
C
C  Perform the complex multiplication of the ADC matrix
C  with the Krylov vectors. The ADC matrix is Hermitian !
C  This has to be taken into account. The diagonal elements
C  are necessarily real.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER NBUFS,IREC,NACT,JDUMMY,IROW,ICOL
      COMPLEX*16 A
C
C---------------Executable code--------------------------------------

C
C  the diagonal elements
C
      REWIND(IO_DIA)
      READ(IO_DIA,ERR=88,END=99) (DIAG(IXX),IXX=1,LADC),NBUFS

      DO J = 1, LADC
        PN(J) = PO(J) * DIAG(J)
      ENDDO

C  for the off-diagonal elements we note that only the lower triangular
C  part of the ADC-matrix is stored. the missing contribution is from
C  the complex conjugate upper triangular part !

c  .. read off-diagonal elements

      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)
            PN(IROW) = PN(IROW) + PO(ICOL) * A
            PN(ICOL) = PN(ICOL) + PO(IROW) * DCONJG(A)
        ENDDO
      ENDDO

      RETURN

 88   CALL QUIT('Diagonal read error in CMATBLK')
 99   CALL QUIT('Unexpected end in CMATBLK')
 199  CALL QUIT('Off-diagonal read error in CMATBLK')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMATBLK_ND(LADC,PO,PN,BUF,IOI,IOJ,INTBUF,IO_OFF)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LADC,INTBUF,IO_OFF
      COMPLEX*16 PO(LADC), PN(LADC)
      COMPLEX*16 BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)

C
C---------------Description--------------------------------------------
C
C  Perform the complex hermitian multiplication of the ADC matrix
C  with the Krylov vectors. Here *NO* separate diagonal entries are
C  processed ! They are contained in the matrix and
C  a single multiplication suffices.
c
c  att! the required number of buffers is communicated through this
c  common block.
C
C---------------Common Blocks--------------------------------------
C
      INTEGER NBUFS
      COMMON/MATCHUNKS/NBUFS
C
C---------------Local variables--------------------------------------
C
      INTEGER IREC,NACT,JDUMMY,IROW,ICOL
      COMPLEX*16 A
C
C---------------Executable code--------------------------------------

C  for the matrix elements we note that only the lower triangular
C  part *AND* the diagonal of the ADC-matrix is stored. 
C  the missing contribution is from the symmetric upper triangular part !

c  the if takes care of not double counting the diagonal.
c
c  new vector has to be zeroed!
      DO J = 1, LADC
        PN(J) =  (0.0d0,0.0d0)
      ENDDO


c  .. read off-diagonal elements and multiply on the fly

      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * DCONJG(A)
        ENDDO
      ENDDO

      RETURN

 889  CALL QUIT('Matrix read error in CMATBLK_ND')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMATBLK_IC(LADC,PO,PN,BUF,IOI,IOJ,INTBUF,IO_OFF,
     &                      nincw,bufc,ioic,iojc,ivalid)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER     LADC,INTBUF,IO_OFF
      complex*16  PO(LADC), PN(LADC)
      complex*16  BUF(INTBUF)
      INTEGER     IOI(INTBUF),IOJ(INTBUF)
      integer     nincw
      complex*16  bufc(nincw)
      integer     ioic(nincw),iojc(nincw),ivalid

C
C---------------Description--------------------------------------------
C
C  Perform the complex hermitian multiplication of the ADC matrix
C  with the Krylov vectors. Here *NO* separate diagonal entries are
C  processed ! They are contained in the matrix already.
c
c  att! we have the icore version here.  Matmul splits in an incore
c  part and an external part for those buffers that can not be held
c  in memory.
C
C---------------Common Blocks--------------------------------------
C
      INTEGER NBUFS
      COMMON/MATCHUNKS/NBUFS
C
C---------------Local variables--------------------------------------
C
      INTEGER IREC,NACT,JDUMMY,IROW,ICOL
      complex*16 A
C
C---------------Executable code--------------------------------------

C  for the matrix elements we note that only the lower triangular
C  part *AND* the diagonal of the ADC-matrix is stored. 
C  the missing contribution is from the symmetric upper triangular part !

c  the if takes care of not double counting the diagonal.
c  since the diagonal is now *contained* in the matrix

c  new vector has to be zeroed!

      DO J = 1, LADC
        PN(J) = (0.0d0,0.0d0)
      ENDDO
c
c  do incore part. In this part of the program there *IS* an
c  incore contribution. So we do not need to check this condition
c  any further.
c
      DO K = 1, IVALID
        IROW = ioic(K)
        ICOL = iojc(K)
        A    = bufc(K)
        PN(IROW) = PN(IROW) + PO(ICOL) * A
        IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * DCONJG(A)
      ENDDO
c
c  next do out of core part (if there is any)
c  remark: if matrix is completely incore the following loop
c  will never be executed. (nbufs = 0).
c  Under the IO_OFF file handle we find the *truncated* ADC
c  matrix only containing the ooc buffers. (there is no lseek
c  for this type of sequential files)
c
      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * DCONJG(A)
        ENDDO
      ENDDO

      RETURN

 889  CALL QUIT('Matrix read error in CMATBLK_IC')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DIAG_DAVI(iw)

      use adc_cfg   !make input variables accessible here
!                   !abandons the old adc_input common block.
      use adc_mat   !contains matrix descriptor for the matrix
                    !in the actual run. Abandons INFOFILE.

      IMPLICIT NONE

      INTEGER                      :: iw
!
!---------------Description--------------------------------------------
!
! entry point for Davidson diagonalisator
!
!
! *** variables for Matrix descriptor file
      INTEGER         :: iobase
      INTEGER         :: intbuf,desrep,rcw
      INTEGER         :: ladc,irecl,nmain,nbufs
      CHARACTER*6     :: fileadc,filecnf
      CHARACTER*5     :: nmspec 
! *** variables for Davidson control
      INTEGER         :: nroots,maxdavsp,maxdavit
      REAL*8          :: convthr
! *** end definition variables for Matrix descriptor file
!
!---------------Local variables----------------------------------------
!
      integer         :: i,j,k,ialloc
      integer         :: nsiz
      integer         :: iroot,locdim
      integer         :: lwork,info
      integer         :: idavit,iobase1
      integer, allocatable, dimension(:)  :: ievec
      logical         :: notallconv,proxflag,orthoviolflag

      real*8, allocatable, dimension(:,:)   :: bold,bnew,aij,evint
      real*8, allocatable, dimension(:)     :: lambda,reserr,work
      real*8, allocatable, dimension(:)     :: qm,adiag,xi,dmp1
      real*8                                :: qmn,a1,a2,ortho
      real*8,parameter                      :: autoev = 27.2113957D0
      character*11                          :: davstartname
      character*11                          :: davevcname
      character*8                           :: dsroot ='DAVSTART'
      character*8                           :: dsroot2='DAVILONG'
      character*8                           :: dsroot3='DAVISHRT'
      character*60                          :: listout
      character*56                          :: field
      logical, allocatable, dimension(:)    :: notconverged

!
!---------------Executable code----------------------------------------
!

      Call PST('Entering (real) Davidson diagonalizer+')
!
!  fill current descriptor variables with transferred values
!
      iobase   =  reladc_md_iobase      !file handle for matrix handling
      intbuf   =  reladc_md_intbuf      !buffer length in matrix
      desrep   =  reladc_md_desrep      !current symmetry
      rcw      =  reladc_md_rcw         !real/complex matrix
      ladc     =  reladc_md_matdim      !matrix dimension
      fileadc  =  reladc_md_fileadc     !file name of ADC matrix
      filecnf  =  reladc_md_filecnf     !file name of configuration file
      nmspec   =  reladc_md_nmspec      !file name of spectral output file
      irecl    =  reladc_md_irecl       !record length of configuration data
      nmain    =  reladc_md_nmain       !size of ADC matrix main block 
      nbufs    =  reladc_md_nbufs       !number of buffers of len "intbuf"
      nroots   =  reladc_md_davroots    !number of requested Davidson roots
      maxdavsp =  reladc_md_maxdavsp    !maximum dimension of Davidson space
      maxdavit =  reladc_md_maxdavit    !maximum number of Dav. iterations
      convthr  =  reladc_md_convthr     !convergence of Davidson eigenvectors
!
!  write out matrix  and Davidson parameters
!
      WRITE(iw,*) ' ------ Matrix parameters -----------'
      WRITE(iw,*)
      WRITE(iw,*) 'Write buffer length is:           ',intbuf
      WRITE(iw,*) 'Final states have symmetry:       ',desrep
      WRITE(iw,*) 'ADC matrix is (1=real,2=complex): ',rcw
      WRITE(iw,*) 'ADC matrix dimension (# of rows): ',ladc
      WRITE(iw,*) 'Packed matrix stored in:          ',fileadc
      WRITE(iw,*) 'Configuration data stored in:     ',filecnf
      WRITE(iw,*) 'Eigenvalues written to:           ',nmspec,'.XX'
      WRITE(iw,*) 'Record length of Config file:     ',irecl
      WRITE(iw,*) 'Size of main space:               ',nmain
      WRITE(iw,*) 'Number of chunks in ADC file:     ',nbufs
      WRITE(iw,*)
      WRITE(iw,*) ' ------ Davidson parameters -----------'
      WRITE(iw,*)
      WRITE(iw,*) 'Number of Davidson roots:         ',nroots
      WRITE(iw,*) 'Maximum Davidson subspace:        ',maxdavsp
      WRITE(iw,*) 'Number of Macro iterations:       ',maxdavit
      WRITE(iw,*) 'Eigenvalue convergence up to:     ',convthr
!
!   allocate work space arrays
!
      lwork = 3*maxdavsp
      allocate(aij(maxdavsp,maxdavsp))   !Dav. overlap matrix (DOM)
      allocate(lambda(maxdavsp))         !Eigenvalues of DOM
      allocate(reserr(maxdavsp))         !Residual errors of eigenvalues
      allocate(work(lwork))              !workspace for DSYEV
      allocate(qm(ladc))                 !intermediate vector
      allocate(adiag(ladc))              !vector of ADC diagonal elements
      allocate(xi(ladc))                 !intermediate vector
      allocate(ievec(ladc))
      allocate(dmp1(ladc))
      allocate(notconverged(maxdavsp))

      allocate(bold(ladc,maxdavsp),stat=ialloc)
      if(ialloc.ne.0) then
        call quit('Allocation problem in DO_DAVI!')
      endif

      allocate(bnew(ladc,maxdavsp),stat=ialloc)
      if(ialloc.ne.0) then
        call quit('Allocation problem in DO_DAVI!')
      endif

      allocate(evint(ladc,nroots),stat=ialloc)
      if(ialloc.ne.0) then
        call quit('Allocation problem in DO_DAVI!')
      endif
!
! clear relevant arrays
!
      aij    = 0.0d0
      reserr = 0.0d0
      qm     = 0.0d0
      adiag  = 0.0d0
      xi     = 0.0d0
      dmp1   = 0.0d0
      bold   = 0.0d0
      bnew   = 0.0d0
      evint  = 0.0d0
!__________________________________________________________
!| open file of main block eigenvectors and read them in 
!| select the first 1,...,nroots vectors from the file
!| which are the lowest ones in the main space.
!|
      IF(desrep.GT.9) THEN
        WRITE(davstartname,'(A8,A1,I2)') dsroot,'.',desrep
      ELSE
        WRITE(davstartname,'(A8,A2,I1)') dsroot,'.0',desrep
      ENDIF

      OPEN(iobase,file=davstartname,access='sequential',
     &     status='unknown',form='unformatted')
      READ(iobase) nsiz
      WRITE(*,*) 'Found',nsiz,' Start vectors on file ',davstartname

      If(2*nroots.gt.maxdavsp) then
         Call Quit('Increase Davidson space!')
      Endif

      If(nroots.gt.nsiz) then
        Write(*,*) 'You request more states than there',
     &             ' are in the Main space.'
        Write(*,*) 'Adapting requested roots',nroots,' to',nsiz
        nroots = nsiz
      Endif

      Do i=1,nroots
        Read(iobase) (bold(j,i),j=1,nsiz)
      Enddo
      Close(iobase,status='delete')
!|
!|  end reading initial vectors and diagonal.
!|___________________________________________________________

      WRITE(*,*) '  '
!
!  from now on the ADC matrix file handle is used and the
!  file is open. Diagonal elements are fetched once.
!
      orthoviolflag = .false.

      Open(iobase,FILE=fileadc,FORM='UNFORMATTED',STATUS='UNKNOWN')
      Call XMATINI_R(iobase,intbuf,nbufs,ladc,adiag) 
!__________________________________________________________
!|
!|
!|             DAVIDSON  MACRO ITERATION LOOP
!|
      Do idavit = 1,maxdavit

      Write(*,*) '  ****   MACRO iteration',idavit
      Write(*,*)
!
!      STEP  A.
!
      locdim = nroots    ! we start with the smallest dav space
      notconverged = .true.  !all vectors are at first not converged
!
! check convergence of the start vectors
! 
      ortho = 0.0d0
      Do i=1,nroots
      Do j=1,nroots
        a1 = dot_product(bold(:,i),bold(:,j))
        ortho = ortho + abs(a1)
      Enddo
      Enddo
      ortho = ortho - dble(nroots)
      WRITE(*,*) 'Orthogonality of the Davidson start vectors:',ortho
      WRITE(*,*)
      WRITE(*,'(3X,A,4X,A,7X,A,7X,A,5X,A)') 'root','dim',
     &           'exc.energy','error','status'
      WRITE(*,'(56A1)') ('-',i=1,56)

      Call XMVMULBL_R(iobase,intbuf,nbufs,ladc,locdim,
     &                bold,bnew)
    
!_______________________________________________________________
!|
!|         Davidson micro iteration loop
!|
      DO      ! MICRO iteration loop. Terminated by an exit statement
!
!  form  b_i*A*b_j and diagonalize. Matrix is replaced by eigenvectors.
!  Compiler knows length of matrix columns
!  AIJ is symmetric because ADC matrix is symmetric.
!
        Do j=1,locdim
        Do i=j,locdim
          aij(i,j) = dot_product(bold(:,i),bnew(:,j))
          aij(j,i) = aij(i,j)   !symmetrize
        Enddo
        Enddo
        CALL DSYEV('V','L',locdim,aij,maxdavsp,lambda,work,lwork,info)
        If(INFO.NE.0) Then
          Call Quit('Error in DSYEV (diag_davi routine)')
        Endif
!
!  STEP B + C.  Form residual vector q_m for each eigenvalue k
!
!***************************** loop over the roots
        Do k = 1,nroots   ! loop over sought eigenvalues lambda_k
!*****************************

          proxflag = .false.

          qm = 0.0d0
          Do i = 1,locdim
            qm(:) = qm(:) + aij(i,k)*bnew(:,i)
          Enddo

          evint(:,k) = qm(:)   !if converged, then we construct the normalized eigenvector.
          a1 = dot_product(evint(:,k),evint(:,k))
          a1 = sqrt(a1)
          evint(:,k) = evint(:,k)/a1

          Do i = 1,locdim
            qm(:) = qm(:) - aij(i,k)*lambda(k)*bold(:,i)
          Enddo
          qmn = sqrt(dot_product(qm,qm))
          reserr(k) = qmn
         
          if(qmn.le.convthr) then
            notconverged(k) = .false.
!           write(*,*) '+++++ Root',k,' converged:',qmn,convthr
!           Write(*,*) '+++++',lambda(k),' Res.norm:',qmn
          else
!           Write(*,*)  k,lambda(k),' Res.norm:',qmn
            notconverged(k) = .true.
          endif
          write(listout,'(A3,I3,A3,I5,6X,F10.7,6X,F10.7)') '   ',
     &           k,'   ',locdim,lambda(k),qmn
          if(notconverged(k).eqv..false.) then
            listout(47:48) =  ' *'
          else
            listout(47:48) =  ' o'
          endif
!
!  STEP D. form new preconditioner Xi
!  here one has to take care of the division!!!
!  It may happen that lambda(k) is approx. eq. to adiag(j)!
!
          Do j=1,ladc
            a1=lambda(k) - adiag(j)
            if(dabs(a1).lt.1.0e-03) then
              proxflag = .true.
              listout(49:51) =  ' PX'
              xi(j) = 1.0d0
            else
              xi(j)=qm(j)/a1
            endif
          Enddo
!         write(*,*) '   @Norm of preconditioner xi:',
!    &           sqrt(dot_product(xi(:),xi(:)))
!
!  STEP E. form  vector d_(M+1) for the kth eigenvalue
!  orthogonalize to *all* previous b_i directions from xi
!  Importance of new vector is decided on its norm !
!
          dmp1 = xi
          Do i = 1,locdim + k - 1
            a1=dot_product(xi,bold(:,i))
            dmp1(:) = dmp1(:) - a1*bold(:,i)
          Enddo
          a1 = sqrt(dot_product(dmp1,dmp1))
!         write(*,*) '   @Norm of preconditioner d(M+1):',a1
!
!  STEP F. form  vector b_(M+1): b_(M+1)=d_(M+1)/||d_(M+1)||
!  this vector has to be included in the array of old vectors!
!  important. This is useful only if the vector b_(M+1) adds
!  new information (should be linearly independent)
!
          dmp1(:)=dmp1(:)/a1
          bold(:,locdim+k)=dmp1(:)
          xi = bold(:,locdim+k)
          a2 = sqrt(dot_product(xi,xi))
!
!  intermediate ortho check of new b vector.
!
          ortho=0.0d0
          xi = bold(:,locdim+k)
          Do i=1,locdim+k-1
            a1 = dot_product(xi,bold(:,i))
            ortho = ortho + a1
          Enddo
          If(abs(ortho).lt.1.0E-06) then
            listout(52:54) = ' OK'
          else
            listout(52:54) = ' VI'
            orthoviolflag = .true.
          endif
!         WRITE(*,*) 'Ortho. of the newly added b-vector',ortho
!
!  STEP G.  Form A*b_(M+1)
!  this vector has to be included in the array of new vectors!
!
          Call XMVMUL_R(iobase,intbuf,nbufs,ladc,
     &                  bold(:,locdim+k),bnew(:,locdim+k))

          Write(*,*)  listout
!
! check orthogonality violation of newly added vectors.
! If it can not be maintained we have to exit in order to
! avoid intruder states.
!
          If(orthoviolflag.eqv..true.) EXIT   !leave this do loop


!*****************************  end loop over the roots
        Enddo
!*****************************
        Write(*,*) '    '        !spacer

!
! check vector convergence in microiterations
!
        notallconv = .false.
        do k=1,nroots
          notallconv=notallconv.or.notconverged(k)
        enddo

        if(notallconv.eqv..false.) then
          write(*,*) 'All vectors converged. Leaving',
     &               ' microiteration loop.'
          exit
        endif

        If(orthoviolflag.eqv..true.) exit
      
        locdim = locdim + nroots
!       write(*,*) 'New local dimension:',locdim
!
! check for macroiteration, transfer new start vector set.
! Attention!! It was observed that orthogonality slowly gets
! lost when macroiterations are entered. We therefore
! reorthogonalize the new start vector set.
!
        if( (locdim + nroots).gt.maxdavsp) then
          write(*,*) 'Davidson space exhausted.'
          write(*,*) 'Generating new (orthogonalized) start vectors.'
          write(*,*) 'For the next macroiteration.'
          do k=1,nroots
            bold(:,k)=evint(:,k)
          enddo
          call reortho(bold,ladc,nroots)
          bnew = 0.0d0
          exit
        endif

      ENDDO
!|
!|         End Davidson micro iteration loop
!|
!|______________________________________________________________
!
!  we have to check convergence again to distinguish between next
!  macroiteration or end of davidson
!
      if(notallconv.eqv..false.) exit

      If(orthoviolflag.eqv..true.) exit

      ENDDO   !idavit loop

      write(*,*) 'Leaving Davidson macro iteration loop.'

      If(orthoviolflag.eqv..true.) then
        write(*,*) 'Orthogonality of new vectors can not'
        write(*,*) 'be maintained any longer. Exiting...'
      Endif

      if(notallconv.eqv..false.) then
        write(*,*) 'All vectors converged to desired threshold'
      else
        write(*,*) '**WARNING: Not all vectors converged to desired'
        write(*,*) '**         Threshold. Check residual errors!'
      endif
!|
!|         End Davidson MACRO iteration loop
!|
!|______________________________________________________________
!
!   close ADC matrix file and delete it
!
      close(iobase,status='delete')
!
!   open configuration file
!
      iobase1 = iobase + 1
      OPEN(iobase1,FILE=FILECNF,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
!
!   write out results as LONG output (symmetry specific)
!
      IF(desrep.GT.9) THEN
        WRITE(davevcname,'(A8,A1,I2)') dsroot2,'.',desrep
      ELSE
        WRITE(davevcname,'(A8,A2,I1)') dsroot2,'.0',desrep
      ENDIF
      OPEN(iobase,file=davevcname,access='sequential',
     &     status='unknown',form='formatted')
!
!__________________ loop over the roots and write out ** LONG ** results
!|
!|
!|
      Do K = 1,nroots
!
! determine pole strength
!
        a1 = 0.0d0
        do j=1,nmain
          a1 = a1 + evint(j,k)*evint(j,k)
        enddo

        Write(iobase,'(A14,X,I4)') 'Excited state:',K
        Write(iobase,'(A19)') '-----------------'
        Write(iobase,'(A12,F10.6,A4,F10.6,A8,F10.6,A5,F10.6,A5,I4)') 
     &   'Exc. energy:',
     &   lambda(K),' au ',lambda(K)*autoev,
     &   ' eV  PS:',a1,' Err:',reserr(k),' Sym:',desrep
        Write(iobase,'(A28,F15.6)') 'Norm^2 of double exc. contr:',
     &         1.0d0 - a1
        Write(iobase,'(A)') ' '
!
! write out all vector components
!
        do j=1,ladc
          Read(iobase1,rec=j) field
          Write(iobase,'(A56,F15.6,3X,F15.6)') field,
     &       evint(j,k),evint(j,k)*evint(j,k)
        enddo
        Write(iobase,'(A17)') '-----------------'
        Write(iobase,'(A)') ' '

      Enddo
      CLOSE(iobase)
!|
!|
!|_______________________________________________
!
!   write out results as SHORT output (symmetry specific)
!
      IF(desrep.GT.9) THEN
        WRITE(davevcname,'(A8,A1,I2)') dsroot3,'.',desrep
      ELSE
        WRITE(davevcname,'(A8,A2,I1)') dsroot3,'.0',desrep
      ENDIF
      OPEN(iobase,file=davevcname,access='sequential',
     &     status='unknown',form='formatted')
!
!__________________ loop over the roots and write out ** SHORT ** results
!|
!|
!|
      Do K = 1,nroots
!
! determine pole strength
!
        a1 = 0.0d0
        do j=1,nmain
          a1 = a1 + evint(j,k)*evint(j,k)
        enddo

        Write(iobase,'(A14,X,I4)') 'Excited state:',K
        Write(iobase,'(A19)') '-----------------'
        Write(iobase,'(A12,F10.6,A4,F10.6,A8,F10.6,A5,F10.6,A5,I4)') 
     &   'Exc. energy:',
     &   lambda(K),' au ',lambda(K)*autoev,
     &   ' eV  PS:',a1,' Err:',reserr(k),' Sym:',desrep
        Write(iobase,'(A28,F15.6)') 'Norm^2 of double exc. contr:',
     &         1.0d0 - a1
        Write(iobase,'(A)') ' '
!
! determine 20 most relevant configurations
! reuse xi-array for modulus calculation and sort
! the corresponding eigenvector with respect to coeffs.
!
        do j=1,ladc;ievec(j)=j;enddo
        xi=-abs(evint(:,k))
        Call DISORT(ladc,xi,ievec)  !sorts from smallest to largest
        Do j=1,min(20,ladc)
          Read(iobase1,rec=ievec(j)) field
          a1 = evint(ievec(j),k)
          Write(iobase,'(A56,F15.6,3X,F15.6)') field,a1,a1*a1
        Enddo
        Write(iobase,'(A17)') '-----------------'
        Write(iobase,'(A)') ' '

      Enddo
      CLOSE(iobase)
!|
!|
!|_______________________________________________

      CLOSE(iobase1)
!
!  deallocate all intermediate arrays
!
      deallocate(aij)
      deallocate(lambda)
      deallocate(work)
      deallocate(qm)
      deallocate(adiag)
      deallocate(xi)
      deallocate(dmp1)
      deallocate(notconverged)
      deallocate(bold)
      deallocate(bnew)
      deallocate(evint)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XMATINI_R(iobase,intbuf,nbufs,ladc,adiag)
C
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      INTEGER                        :: iobase,intbuf,nbufs,ladc
      REAL*8                         :: adiag(ladc)
!
!------------------------ local variables -------------------------
!
      REAL*8,  allocatable, dimension(:)      :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      REAL*8                                  :: A

!
!---------------Description--------------------------------------------
!
!  Fetch the diagonal entries of the sparse matrix A
!  This routine will simplify considerably in the future.
!---------------Executable code--------------------------------------
!
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))

c  .. read elements

      REWIND(iobase)
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

      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 889  CALL QUIT('Matrix read error in XMATINI_R')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XMVMULBL_R(iobase,intbuf,nbufs,ladc,nvecs,PO,PN)
C
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      INTEGER                        :: iobase,intbuf,nbufs,ladc,nvecs
      REAL*8                         :: PO(LADC,nvecs), PN(LADC,nvecs)
!
!------------------------ local variables -------------------------
!
      REAL*8,  allocatable, dimension(:)      :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      INTEGER                                 :: I,J,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      REAL*8                                  :: A

!
!---------------Description--------------------------------------------
!
!  Perform the real symmetric multiplication of the external ADC matrix
!  with the vector PO and store the result in the vector PN
!
!  Matrix format is assumed to include the diagonal (new format)
!
!  Intermediate buffers are allocated internally !
!
!---------------Executable code--------------------------------------
!
!  for the matrix elements we note that only the lower triangular
!  part *AND* the diagonal of the ADC-matrix is stored. 
!  the if takes care of not double counting the diagonal.
!  since the diagonal is now *contained* in the matrix

      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
!
!  new vector has to be zeroed!
!
      DO J = 1, nvecs
      DO I = 1, ladc
        PN(I,J) =  0.0d0
      ENDDO
      ENDDO

c  .. read off-diagonal elements and multiply on the fly

      REWIND(iobase)
      DO I = 1,NBUFS
        READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          DO J = 1,NVECS
         
            PN(IROW,J) = PN(IROW,J) + PO(ICOL,J) * A
            IF(IROW.NE.ICOL) PN(ICOL,J) = PN(ICOL,J) + PO(IROW,J) * A

          ENDDO
        ENDDO
      ENDDO

      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 889  CALL QUIT('Matrix read error in XMVMULBL_R')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XMVMUL_R(iobase,intbuf,nbufs,ladc,PO,PN)
C
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      INTEGER                        :: iobase,intbuf,nbufs,ladc
      REAL*8                         :: PO(LADC), PN(LADC)
!
!------------------------ local variables -------------------------
!
      REAL*8,  allocatable, dimension(:)      :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      REAL*8                                  :: A

!
!---------------Description--------------------------------------------
!
!  Perform the real symmetric multiplication of the external ADC matrix
!  with the vector PO and store the result in the vector PN
!
!  Matrix format is assumed to include the diagonal (new format)
!
!  Intermediate buffers are allocated internally !
!
!---------------Executable code--------------------------------------
!
!  for the matrix elements we note that only the lower triangular
!  part *AND* the diagonal of the ADC-matrix is stored. 
!  the if takes care of not double counting the diagonal.
!  since the diagonal is now *contained* in the matrix

      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
!
!  new vector has to be zeroed!
!
      DO I = 1, LADC
        PN(I) =  0.0d0
      ENDDO

c  .. read off-diagonal elements and multiply on the fly

      REWIND(iobase)
      DO I = 1,NBUFS
        READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
        ENDDO
      ENDDO

      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 889  CALL QUIT('Matrix read error in XMVMUL_R')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE REORTHO(v,ladc,nroots)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      INTEGER                        :: ladc,nroots
      REAL*8                         :: v(ladc,nroots)
!
!------------------------ local variables -------------------------
!
      REAL*8,  allocatable, dimension(:)     :: x
      REAL*8                                 :: a1
      INTEGER                                :: i,k

!
!---------------Description--------------------------------------------
!
!     write(*,*) 'In reortho ***********************'
      allocate(x(ladc))
!
! reorthogonalize nvec vectors stored columnwise in the array a
!
      Do k=1,nroots
        x = v(:,k)
        Do i=1,k-1
          a1 = dot_product(x(:),v(:,i))
!         write(*,*) i,k,a1
          x = x - a1*v(:,i)
        Enddo
        a1 = sqrt(dot_product(x,x))
!       write(*,*) 'norm k:',k,a1
        v(:,k)  =  x/a1
      Enddo
!     write(*,*) 'Exiting  reortho ***********************'
!
!---------------Executable part----------------------------------------
!
      deallocate(x)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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
      SUBROUTINE WCONDAT_S(LUN,FILENAME,IRECL,DESREP,LENADCM)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Writes explicit irep names and spinor numbers for the various
C     final state configurations belonging to symmetry DESREP.
C     Avoiding too complicated things we fix IRECL in this routine
c     no matter which value is transferred to this routine. Therefore
c     IRECL is an OUTPUT parameter.
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*6 FILENAME
      INTEGER LUN,IRECL,DESREP,LENADCM
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER lcount
      CHARACTER*44 FIELD44
      INTEGER bl_conf
      CHARACTER*4 bl_irna,bl_abs
      
C
C---------------Executable code--------------------------------------
C
C  Code for ** single ionization **:  1h/2h1p
C
c fix record length for the direct access file and set some
c blank variables. Then open file and initialize configuration
c counter.
c *ATT* Is has turned out advantageous that storing the *absolute*
c spinor numbers is superior to storing the relative ones for the
c subsequent postprocessing!
c
c this will be followed by the single ionization case for the moment
c (1h and 2h1p) configurations
c
      IRECL = 44
      BL_CONF = 0
      BL_IRNA = '    '
      BL_ABS  = 'abs.'
      OPEN(LUN,FILE=FILENAME,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
      LCOUNT = 1

C_________________________ SIPS ______________________________
C|
C|
C|
        WRITE(IW,*) 'Writing SIP configs to ',FILENAME
        WRITE(IW,*) 'Final state symmetry: ',DESREP

c ** start with 1h configurations

        DO k=1,no(DESREP)
          koff=io(DESREP)+k
          WRITE(FIELD44,381) 
     &            k,repna(DESREP),
     &            koff,repna(DESREP),
     &            bl_conf,bl_irna,
     &            bl_conf,bl_irna 
          WRITE(LUN,REC=lcount) FIELD44
          lcount = lcount + 1
        ENDDO

c ** continue with 2h1p configurations
c ** write absolute spinor numbers

        DO 60 KLREP=1,NREP
          AREP=MULTB(DESREP,KLREP+NREP,2)
          DO 45 LREP=1,NREP
           KREP=MULTB(LREP,KLREP+NREP,2)
           IF(KREP.LT.LREP) GOTO 45
           DO L=1,NO(LREP)
             LOFF=IO(LREP)+L
             KMIN=1
             IF(KREP.EQ.LREP) KMIN = L + 1
             DO K=KMIN,NO(KREP)
               KOFF=IO(KREP)+K
               DO A=1,NV(AREP)
                 AOFF=IO(NREP+1)+IV(AREP)+A
c                WRITE(FIELD44,381) 
c    &             K,REPNA(KREP),
c    &             L,REPNA(LREP),
c    &             A,REPNA(AREP),
c    &             BL_CONF,BL_IRNA
                 WRITE(FIELD44,381) 
     &             KOFF,REPNA(KREP),
     &             LOFF,REPNA(LREP),
     &             AOFF,REPNA(AREP),
     &             BL_CONF,BL_IRNA
                 WRITE(LUN,REC=lcount) FIELD44
                 LCOUNT=LCOUNT+1
               ENDDO
             ENDDO
           ENDDO
 45       CONTINUE
 60     CONTINUE

C ** file is open only within this subroutine

        CLOSE(LUN)

C ** correct configuration counter

        LCOUNT = LCOUNT - 1
        IF(LCOUNT.NE.LENADCM) THEN
          WRITE(IW,*) 'LCOUNT,LENADCM:',LCOUNT,LENADCM
          CALL QUIT('Error in WCONDAT_S for single ionizations')
        ENDIF
C|
C|
C|____________________________________________________________

 381  FORMAT(4(I3,'(',A4,')  '))

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WCONDAT_D(LUN,FILENAME,IRECL,DESREP,LENADCM,
     &                     MAXOCCSP,IROOO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Writes explicit irep names and spinor numbers for the various
C     final state configurations belonging to  symmetry DESREP.
C     Avoiding too complicated things we fix IRECL in this routine
c     no matter which value is transferred to this routine. Therefore
c     IRECL is an OUTPUT parameter.
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*6 FILENAME
      INTEGER LUN,IRECL,DESREP,LENADCM
      INTEGER MAXOCCSP,IROOO(MAXOCCSP,2)
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER lcount
      CHARACTER*44 FIELD44
      INTEGER bl_conf
      CHARACTER*4 bl_irna,bl_abs
      
C
C---------------Executable code--------------------------------------
C
C  Do the hole space configuration analysis for the
C  DIP case exclusively in this routine.
c  Initialize arrays:
c
      IX=1
      Do KRP=1,NREP
        Do K=1,NO(KRP)
          IROOO(IX,1)=K
          IROOO(IX,2)=KRP
          IX = IX + 1
        Enddo
      Enddo
C
c fix record length for the direct access file and set some
c blank variables. Then open file and initialize configuration
c counter.
 
      IRECL = 44
      BL_CONF = 0
      BL_IRNA = '    '
      BL_ABS  = 'abs.'
      OPEN(LUN,FILE=FILENAME,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
      LCOUNT = 1

C____________ START WITH DIPS ______________________________
C|
C|
C|
      WRITE(IW,*) 'Writing DIP configs to ',FILENAME
      WRITE(IW,*) 'Final state symmetry: ',DESREP

      DO 14 JRP = 1, NREP
        IRP = MULTB(JRP,DESREP+NREP,2)
        IF (IRP.LT.JRP) GOTO 14
        DO J = 1, NO(JRP)
          IMIN = 1
          IF (IRP.EQ.JRP) IMIN = J + 1
          DO I = IMIN, NO(IRP)
c
c  -- additionally determine absolute spinor numbers
c  -- in the active space for mulpop
c
            ISPINOR=IO(IRP)+I
            JSPINOR=IO(JRP)+J
            WRITE(FIELD44,381)
     &         I,REPNA(IRP),
     &         J,REPNA(JRP),
     &         ISPINOR,BL_ABS,
     &         JSPINOR,BL_ABS 
            WRITE(LUN,REC=LCOUNT) FIELD44
            LCOUNT = LCOUNT + 1
c
          ENDDO
        ENDDO
  14  CONTINUE
      IF( (LCOUNT-1).NE.NOOT(DESREP)) THEN
        WRITE(IW,*) 'LCOUNT,LENADCM:',LCOUNT,NOOT(DESREP)
        CALL QUIT('Error in WCONDAT_D!')
      ENDIF
c
c  determine satellite configurations
c  Since we work in a Kramers-paired formalism (closed shell)
c  then number of virtuals in RRP is equal to the number in its
c  Kramers irrep (RRP*). For the output RRP has to be the Kramers
c  partner since internally the cc of this orbital is used !
c
      DO K=1,MAXOCCSP
        KSP = IROOO(K,1)
        KRP = IROOO(K,2)
        DO J=K+1,MAXOCCSP
          JSP = IROOO(J,1)
          JRP = IROOO(J,2)
          JKRP = MULTB(JRP,KRP,1)
          DO I=J+1,MAXOCCSP
            ISP = IROOO(I,1)
            IRP = IROOO(I,2)
            IJKRP = MULTB(IRP,JKRP+NREP,1)
            RRP = MULTB(DESREP+NREP,IJKRP,2)
c determine Kramers-partner of RRP (is cc of RRP)
            RRPKR = MULTB(RRP,1+NREP,2)
            DO RSP = 1,NV(RRP)
            WRITE(FIELD44,381) ISP,REPNA(IRP),
     &                         JSP,REPNA(JRP),
     &                         KSP,REPNA(KRP),
     &                         RSP,REPNA(RRPKR)
            WRITE(LUN,REC=LCOUNT) FIELD44
            LCOUNT = LCOUNT + 1

            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CLOSE(LUN)
C|
C|
C|____________________________________________________________

 381  FORMAT(4(I3,'(',A4,')  '))

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WCONDAT_X(LUN,FILENAME,IRECL,DESREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Writes explicit configuration information for the excited
C     states in symmetry DESREP.
C
C---------------Calling variables--------------------------------------
C
      INTEGER LUN,IRECL,DESREP
      CHARACTER*6 FILENAME
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
      CHARACTER*56 FIELD
      INTEGER bl_conf
      CHARACTER*4 bl_irna
      INTEGER RECCOUNT
      
C
C---------------Executable code--------------------------------------
C
C open config file
C
      IRECL = 56
      OPEN(LUN,FILE=FILENAME,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
C
C First do the H--->P part. Column counter is the decisive
C IRREP.
C
      BL_CONF = 0
      BL_IRNA = '    '
      RECCOUNT = 1

C____________ START WITH H ---> P __________________________
C|
C|
C|

      do jrep = 1,nrep   !jrep,brep are fermionic, desrep is bosonic
        brep = multb(desrep+nrep,jrep,2)
        do jfun = 1,NO(jrep)
        do bfun = 1,NV(brep)
          jabs=IO(jrep) + jfun
          babs=IO(nrep+1) + IV(brep) + bfun

          WRITE(FIELD,333)
     &         jabs,REPNA(jrep),
     &         babs,REPNA(brep),
     &         BL_CONF,BL_IRNA,
     &         BL_CONF,BL_IRNA 

          WRITE(LUN,REC=RECCOUNT) FIELD
          RECCOUNT = RECCOUNT + 1

        enddo
        enddo
      enddo
      if( (reccount-1).ne.MVO(desrep))
     &   stop 'Inconsistent rec counter in WCONDAT_X!'
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
     &         kabs,REPNA(krep),
     &         ccabs,REPNA(crep),
     &         labs,REPNA(lrep),
     &         ddabs,REPNA(drep)

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

 333  FORMAT(2(I3,' (',A4,')  ---> ',I3,' (',A4,') '))
 
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
