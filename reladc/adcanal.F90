!
!
!  ************************************************************
!  ************************************************************
!  *******   ADCANAL                                 **********
!  *******   In this file all the analysis routines  **********
!  *******   for RELADC are combined.                **********
!  *******                                           **********
!  ************************************************************
!  ************************************************************
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine RADCVECANALYZE_F90(io1,io2,io3,desrep,ladc,nocc, & 
                                    irecl,filecnf)

      implicit none
!
! transferred variables
!
      integer, intent(in) :: io1,io2,io3,desrep,ladc,nocc,irecl
      character*6, intent(in) :: filecnf
!
! local variables
!
      CHARACTER*10                :: LEVECFN,ANALOUT
      CHARACTER*44                :: field
      CHARACTER*11                :: blankf='           '
      logical                     :: isthere
      integer                     :: neigenv,j,ix,ladcx,ivcount,nv,idum
      integer, allocatable        :: indarr(:), indarr2(:)
      real*8, allocatable         :: evec(:),eval(:)
      real*8                      :: enorm,polestr,polethr,diffthr,seval
      real*8                      :: autoev = 27.2113957d0
!
! ------------------------ execution starts ----------------------
!
      call pst('Analyzer for long ADC real eigenvectors+')
      write(*,*) 'Calling arguments:'
      write(*,*) 'Symmetry of final state:',DESREP
      write(*,*) 'Length of ADC matrix:   ',LADC
      write(*,*) 'Length of main space:   ',NOCC
      write(*,*) 'Record length of config file:',IRECL
      write(*,*) 'Generic config file  name: ',FILECNF

! constructing eigenvector and analysis output file name

      IF(DESREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',DESREP
        WRITE(ANALOUT,'(A8,I2)') 'EVCANAL.',DESREP
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',DESREP
        WRITE(ANALOUT,'(A9,I1)') 'EVCANAL.0',DESREP
      ENDIF
      write(*,*) 'Fetching eigenvectors from: ',LEVECFN
      write(*,*) 'Writing analysis results to: ',ANALOUT

! check the existence of eigenvector and configuration file

      INQUIRE(file=levecfn,exist=isthere)
      if(.not.isthere) then
        write(*,*) '******* ATTENTION ********'
        write(*,*) 'Eigenvector file with ADC vectors not present'
        write(*,*) 'for symmetry',desrep
        write(*,*) 'Skipping analysis'
        return
      endif
      INQUIRE(file=filecnf,exist=isthere)
      if(.not.isthere) then
        write(*,*) '******* ATTENTION ********'
        write(*,*) 'Configuration file not present'
        write(*,*) 'for symmetry',desrep
        write(*,*) 'Skipping analysis'
        return
      endif

! opening eigenvector, configuration and output file

      OPEN(io1, FILE=LEVECFN, FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL', STATUS='old')
      OPEN(io2,FILE=FILECNF,ACCESS='DIRECT',RECL=IRECL, &
          STATUS='UNKNOWN')
      OPEN(io3, FILE=ANALOUT, FORM='FORMATTED', &
          ACCESS='SEQUENTIAL', STATUS='UNKNOWN')

! start analysis

      read(io1) ladcx,neigenv
      write(*,*) 'found',neigenv,' eigenvectors on ',levecfn
      write(*,*) 'of length',ladcx
      if(ladcx.ne.ladc) then
        call quit('internal error in adcanalyze_f90!')
      endif


      allocate(evec(ladcx))
      allocate(eval(neigenv))
      allocate(indarr(neigenv))

      polethr = 0.001d0; indarr=0; indarr(1)=1

! make a first browse for the polestrengths

      ivcount=1
      DO J=1,NEIGENV
        read(io1) eval(J)
        read(io1) (evec(ix),ix=1,LADCX)
        enorm=dot_product(evec,evec)
        evec=evec/sqrt(enorm)
        polestr=dot_product(evec(1:nocc),evec(1:nocc))
        write(*,'(A,2F16.8)') 'found eigenvector for', &
                              eval(j)*autoev,polestr
        IF(polestr.GE.polethr) then
          indarr(ivcount)=J; ivcount = ivcount + 1
!         write(*,'(A,I5,2F16.8)') '  ++++ including', &
!               J,eval(j)*autoev,polestr
        ENDIF
      ENDDO
!     write(*,'(A,I5)') '  I grabbed ',ivcount
      nv = ivcount - 1

! make a second browse through the eigenvalues for duplicates
! and modify index array accordingly.

      diffthr = 1.0D-06
      IF(nv.gt.1) THEN
        allocate(indarr2(nv))
        ivcount=1
        DO J=1,nv-1
!         write(*,*) 'poooo:',eval(indarr(j+1))*autoev, &
!                             eval(indarr(j))*autoev
          IF(dabs(eval(indarr(j+1)) - eval(indarr(j))).GE.diffthr) THEN
            indarr2(ivcount) = indarr(J)
            ivcount = ivcount + 1           
          ENDIF
        ENDDO
        indarr2(ivcount) = indarr(J)
        nv = ivcount
        write(*,*) 'Finally after purging:',nv,' vectors:'
        DO J=1,nv
          write(*,*) j,indarr2(j),eval(indarr2(j))*autoev
        ENDDO
      ENDIF

! now loop over all eigenvectors again and write out the relevant ones
! whose indices are stored in indarr2.

      REWIND(io1)
      read(io1) idum,idum  ! read off the first two integers
      ivcount=1
      DO J=1,NEIGENV

        read(io1) seval
        read(io1) (evec(ix),ix=1,LADCX)

        IF(J.ne.indarr2(ivcount)) CYCLE   ! skip printout if not right evec.

        write(*,*) 'writing eigenvector for',seval,seval*autoev
        write(io3,'(A,I4,A,F16.10,A)') 'Eigenvector nr. ', &
              ivcount,'   IP = ',seval*autoev,' eV'
        enorm=dot_product(evec,evec)
        write(*,*) 'norm of eigenvector:',enorm
        write(*,*) 'renormalizing eigenvector.'
        evec=evec/sqrt(enorm)
        polestr=dot_product(evec(1:nocc),evec(1:nocc))
        write(io3,'(A,F16.10)') 'Pole strength: ',polestr
        write(io3,'(A)') ' '
        write(io3,'(A)') ' '
        write(io3,'(A8,2X,A33,A20)') &
          'No.',' 1h  and  2h1p configurations  ',' Coeff(real)'
        write(io3,'(A)') ' '

        DO ix=1,ladcx
          read(io2,REC=ix) field
          IF(ix.gt.nocc) then
!           write(*,'(I8,A33,F20.12)') ix,field,evec(ix)
            write(io3,'(I8,2X,A33,F20.12)') ix,field,evec(ix)
          ELSE
!           write(*,'(I8,3A11,F20.12)') ix,field,blankf,blankf, &
!                                        evec(ix)
            write(io3,'(I8,2X,3A11,F20.12)') ix,field,blankf,blankf, &
                                             evec(ix)
          ENDIF
        ENDDO
        ivcount = ivcount + 1
        write(io3,'(A)') ' '
        write(io3,'(A)') '--------------------------------------------'
        write(io3,'(A)') ' '

      ENDDO

! release memory

      deallocate(indarr2)
      deallocate(indarr)
      deallocate(eval)
      deallocate(evec)

      close(io1)
      close(io2)
      close(io3)
      
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine ZADCVECANALYZE_F90(io1,io2,io3,desrep,ladc,nocc, & 
                                    irecl,filecnf)

      implicit none
!
! transferred variables
!
      integer, intent(in) :: io1,io2,io3,desrep,ladc,nocc,irecl
      character*6, intent(in) :: filecnf
!
! local variables
!
      CHARACTER*10                :: LEVECFN,ANALOUT
      CHARACTER*44                :: field
      CHARACTER*11                :: blankf='           '
      logical                     :: isthere
      integer                     :: neigenv,j,ix,ladcx,ivcount,nv,idum
      integer, allocatable        :: indarr(:), indarr2(:)
      real*8, allocatable         :: eval(:)
      complex*16, allocatable     :: evec(:)
      real*8                      :: enorm,polestr,polethr,diffthr,seval
      real*8                      :: autoev = 27.2113957d0
!
!  function declaration
!
      real*8                      :: dot_productz
!
! ------------------------ execution starts ----------------------
!
      call pst('Analyzer for long ADC complex eigenvectors+')
      write(*,*) 'Calling arguments:'
      write(*,*) 'Symmetry of final state:',DESREP
      write(*,*) 'Length of ADC matrix:   ',LADC
      write(*,*) 'Length of main space:   ',NOCC
      write(*,*) 'Record length of config file:',IRECL
      write(*,*) 'Generic config file  name: ',FILECNF

! constructing eigenvector and analysis output file name

      IF(DESREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',DESREP
        WRITE(ANALOUT,'(A8,I2)') 'EVCANAL.',DESREP
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',DESREP
        WRITE(ANALOUT,'(A9,I1)') 'EVCANAL.0',DESREP
      ENDIF
      write(*,*) 'Fetching eigenvectors from: ',LEVECFN
      write(*,*) 'Writing analysis results to: ',ANALOUT

! check the existence of eigenvector and configuration file

      INQUIRE(file=levecfn,exist=isthere)
      if(.not.isthere) then
        write(*,*) '******* ATTENTION ********'
        write(*,*) 'Eigenvector file with ADC vectors not present'
        write(*,*) 'for symmetry',desrep
        write(*,*) 'Skipping analysis'
        return
      endif
      INQUIRE(file=filecnf,exist=isthere)
      if(.not.isthere) then
        write(*,*) '******* ATTENTION ********'
        write(*,*) 'Configuration file not present'
        write(*,*) 'for symmetry',desrep
        write(*,*) 'Skipping analysis'
        return
      endif

! opening eigenvector, configuration and output file

      OPEN(io1, FILE=LEVECFN, FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL', STATUS='old')
      OPEN(io2,FILE=FILECNF,ACCESS='DIRECT',RECL=IRECL, &
          STATUS='UNKNOWN')
      OPEN(io3, FILE=ANALOUT, FORM='FORMATTED', &
          ACCESS='SEQUENTIAL', STATUS='UNKNOWN')

! start analysis

      read(io1) ladcx,neigenv
      write(*,*) 'found',neigenv,' eigenvectors on ',levecfn
      write(*,*) 'of length',ladcx
      if(ladcx.ne.ladc) then
        call quit('internal error in adcanalyze_f90!')
      endif


      allocate(evec(ladcx))
      allocate(eval(neigenv))
      allocate(indarr(neigenv))

      polethr = 0.001d0; indarr=0; indarr(1)=1

! make a first browse for the polestrengths

      ivcount=1
      DO J=1,NEIGENV
        read(io1) eval(J)
        read(io1) (evec(ix),ix=1,LADCX)
        enorm=dot_productz(evec,ladcx)
        evec=evec/cmplx(sqrt(enorm),0.0d0)
        write(*,*) 'Normchk: ',enorm
        polestr=dot_productz(evec,nocc)
        write(*,'(A,2F16.8)') 'found eigenvector for', &
                              eval(j)*autoev,polestr
        IF(polestr.GE.polethr) then
          indarr(ivcount)=J; ivcount = ivcount + 1
!         write(*,'(A,I5,2F16.8)') '  ++++ including', &
!               J,eval(j)*autoev,polestr
        ENDIF
      ENDDO
!     write(*,'(A,I5)') '  I grabbed ',ivcount
      nv = ivcount - 1

! make a second browse through the eigenvalues for duplicates
! and modify index array accordingly.

      diffthr = 1.0D-06
      IF(nv.gt.1) THEN
        allocate(indarr2(nv))
        ivcount=1
        DO J=1,nv-1
!         write(*,*) 'poooo:',eval(indarr(j+1))*autoev, &
!                             eval(indarr(j))*autoev
          IF(dabs(eval(indarr(j+1)) - eval(indarr(j))).GE.diffthr) THEN
            indarr2(ivcount) = indarr(J)
            ivcount = ivcount + 1           
          ENDIF
        ENDDO
        indarr2(ivcount) = indarr(J)
        nv = ivcount
        write(*,*) 'Finally after purging:',nv,' vectors:'
        DO J=1,nv
          write(*,*) j,indarr2(j),eval(indarr2(j))*autoev
        ENDDO
      ENDIF

! now loop over all eigenvectors again and write out the relevant ones
! whose indices are stored in indarr2.

      REWIND(io1)
      read(io1) idum,idum  ! read off the first two integers
      ivcount=1
      DO J=1,NEIGENV

        read(io1) seval
        read(io1) (evec(ix),ix=1,LADCX)

        IF(J.ne.indarr2(ivcount)) CYCLE   ! skip printout if not right evec.

        write(*,*) 'writing eigenvector for',seval,seval*autoev
        write(io3,'(A,I4,A,F16.10,A)') 'Eigenvector nr. ', &
              ivcount,'   IP = ',seval*autoev,' eV'
        enorm=dot_productz(evec,ladcx)
        write(*,*) 'norm of eigenvector:',enorm
        write(*,*) 'renormalizing eigenvector.'
        evec=evec/cmplx(sqrt(enorm),0.0d0)
        polestr=dot_productz(evec,nocc)
        write(io3,'(A,F16.10)') 'Pole strength: ',polestr
        write(io3,'(A)') ' '
        write(io3,'(A)') ' '
        write(io3,'(A8,2X,A33,2A20)') &
          'No.',' 1h  and  2h1p configurations  ',' Coeff(real)', &
          ' Coeff(imag)'
        write(io3,'(A)') ' '

        DO ix=1,ladcx
          read(io2,REC=ix) field
          IF(ix.gt.nocc) then
!           write(*,'(I8,A33,F20.12)') ix,field,evec(ix)
            write(io3,'(I8,2X,A33,2F20.12)') ix,field,evec(ix)
          ELSE
!           write(*,'(I8,3A11,F20.12)') ix,field,blankf,blankf, &
!                                        evec(ix)
            write(io3,'(I8,2X,3A11,2F20.12)') ix,field,blankf,blankf, &
                                             evec(ix)
          ENDIF
        ENDDO
        ivcount = ivcount + 1
        write(io3,'(A)') ' '
        write(io3,'(A)') '--------------------------------------------'
        write(io3,'(A)') ' '

      ENDDO

! release memory

      deallocate(indarr2)
      deallocate(indarr)
      deallocate(eval)
      deallocate(evec)

      close(io1)
      close(io2)
      close(io3)
      
      end
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
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF), &
                              (IOI(IXX),IXX=1,INTBUF), &
                              (IOJ(IXX),IXX=1,INTBUF), &
                               NACT,JDUMMY
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

      WRITE(*,*) 'Eigenvalues in eV (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential', &
           form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)') &
               '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
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
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF), &
                              (IOI(IXX),IXX=1,INTBUF), &
                              (IOJ(IXX),IXX=1,INTBUF), &
                               NACT,JDUMMY
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
      CALL ZHEEV( JOBZ, UPLO, LADC, A, LADC, EVL, CWORK, &
                  LWORK, RWORK, INFO )
!     write(*,*) 'AFTER  ZHEEV with optimal lwork=',cwork(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in ZHEEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAC (lapack zheev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues in eV (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential', &
           form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)') &
               '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
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
      write(*,*) 'Numer of buffers to read:',nbufs

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
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF), &
                              (IOI(IXX),IXX=1,INTBUF), &
                              (IOJ(IXX),IXX=1,INTBUF), &
                               NACT,JDUMMY
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

!     write(*,*) 'Calling DSYEV with lwork=',lwork
      CALL DSYEV( JOBZ, UPLO, LADC, A, LADC, EVL, WORK, LWORK, INFO )
!     write(*,*) 'AFTER  DSYEV with optimal lwork=',work(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in DSYEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAR (lapack dsyev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues in eV (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential', &
           form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)') &
               '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
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
      SUBROUTINE FULLDIAC2(IOUNIT, FILEN5, INTBUF, NBUFS, LADC, KREP)
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
      write(*,*) 'Numer of buffers to read:',nbufs

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
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF), &
                              (IOI(IXX),IXX=1,INTBUF), &
                              (IOJ(IXX),IXX=1,INTBUF), &
                               NACT,JDUMMY
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
      CALL ZHEEV( JOBZ, UPLO, LADC, A, LADC, EVL, CWORK, &
                  LWORK, RWORK, INFO )
!     write(*,*) 'AFTER  ZHEEV with optimal lwork=',cwork(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in ZHEEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAC (lapack zheev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues in eV (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential', &
           form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)') &
               '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
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
      subroutine TWOHPOP_R(lupri,io_lan,ibeta,nrep,multb,no,noot)
!
!  Purpose: Prepare the propagator contributions to the reduced density matrix.
!
!  the Delta_mn matrix is then written to a file 'DADCPOP'
!  in the format
!              integer dimension
!              real(:,:) dmn   :  real part, vector 1
!              real(:,:) dmn   :  imaginary part, vector 1
!                       ....
!              real(:,:) dmn   :  real part, vector N
!              real(:,:) dmn   :  imaginary part, vector N
!
      implicit none
      integer, intent(in)  :: lupri    ! file handle for printing to stdout
      integer, intent(in)  :: io_lan   ! file handle start
      integer, intent(in)  :: ibeta    ! symmetry of the two-hole final state 
      integer, intent(in)  :: nrep     ! number of ireps in the overall run
      integer, intent(in)  :: multb(64,64,2)     ! fixed dimension multiplication table
      integer, intent(in)  :: no(nrep)    ! number of occupied orbitals in each irep
      integer, intent(in)  :: noot(nrep)  ! number of OO-combinations in each irep
!
!  local variables
!
      real*8,  allocatable  :: dmn(:,:), xmna(:,:), xmn2(:,:),xmnc(:,:)
      real*8,  allocatable  :: xmnd(:,:),thvecs(:,:),evec(:)
      real*8,  allocatable  :: cevl(:), polestr(:)
      real*8,  allocatable  :: dmnext(:,:,:)
      real*8                :: polethr = 0.001d0
      real*8                :: autoev = 27.2113957D0
      real*8                :: A0=0.0d0, A1=1.0d0
      integer, allocatable  :: no_offs(:)
      integer               :: i,j,ix,ivcount,iwwcount
      integer               :: n1,n2
      integer               :: nocc,ipoint,no_sum
      integer               :: irp,jrp,imin,irow,jcol
      logical               :: isthere
      character*7           :: dmn_filename='DADCPOP'

!  say hello and check basic information

      call PST('Performing 2H population analysis (real vectors)+')
      nocc=noot(ibeta)
      write(lupri,*) 'Final state symmetry:',ibeta
      write(lupri,*) 'Number of 2h states in this symm:',nocc

!__________________________________________________________________
!|
!|  read the file of Lanczos eigenvectors and extract the 2h/2h block
!|        ----->   in this symmetry !!  <-----
!|  from it. 
!|  Fill the corresponding arrays with the read vectors and determine
!|  the polestrengths. The raw Lanc eigenvector goes to *evec* and
!|  the 2h/2h block goes to *thvecs*

      Inquire(file='TMATEVC',exist=isthere)
      IF(.not.isthere) THEN
        write(lupri,*) '**************'
        write(lupri,*) '**** Attention'
        write(lupri,*) '**************'
        write(lupri,*) 'TMATEVC file not present.'
        write(lupri,*) 'No population analysis possible'
        return
      ENDIF

      OPEN(IO_LAN,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
      REWIND(IO_LAN)
      READ(IO_LAN) N1
      READ(IO_LAN) N2
      WRITE(lupri,*) 'Found',N1,' Lanczos eigenvectors of length', &
                   N2,' on TMATEVC'
!     WRITE(lupri,*) nocc,' of them are 2H states!'

      allocate(thvecs(nocc,nocc), &
               evec(n2),          &
               cevl(nocc),        &
               polestr(nocc))

      thvecs = 0.0d0; evec = 0.0d0; cevl = 0.0d0; polestr = 0.0d0

      iwwcount = 0
      DO I=1,nocc
        READ(IO_LAN) cevl(i)
        READ(IO_LAN) (evec(ix),ix=1,n2)
        polestr(i)=dot_product(evec(1:nocc),evec(1:nocc))
        if(polestr(i).gt.polethr) then
              write(lupri,'(A,I5,A,2F16.6,A)')           &
                  'Taking vector ',i,                    &
                  ' with', cevl(i)*autoev,polestr(i),    &
                  ' for population analysis'
              iwwcount = iwwcount + 1
              thvecs(:,i) = evec(1:nocc)  ! copy corresponding components
        endif 
      ENDDO
      write(lupri,*) iwwcount, &
        ' vectors have enough pole strength for 2H analysis'
      CLOSE(IO_LAN,STATUS='DELETE')
!|
!|
!|
!|
!|_______ 2h/2h block in memory together with the polestrengths ______


!  prepare occupied orbital offset table

      allocate(no_offs(nrep+1))
      no_offs = 0; no_sum = 0
      DO i=1,nrep
        no_offs(i+1) = no_offs(i) + no(i)
        no_sum = no_sum + no(i)
      ENDDO
!     write(lupri,*) 'offset table:',no_offs
!     write(lupri,*) 'total number of occupied orbitals:',no_sum

!  allocate space for Delta_mn matrix and clear matrix

      allocate(dmn(no_sum,no_sum))
      allocate(dmnext(no_sum,no_sum,2))
      allocate(xmna(no_sum,no_sum))
      allocate(xmnc(no_sum,no_sum))
      allocate(xmnd(no_sum,no_sum))
      allocate(xmn2(no_sum,no_sum))

!  open DADCPOP file for sequential access and write 
!  number of matrices and corresponding matrix dimension
!
!  ATT We check if file exists already and append the next symmetry
!  to it. This liberates us from writing symmetry-specific files altogether.
 
      inquire(file=dmn_filename,exist=isthere)

      select case(isthere)

      case (.FALSE.)      ! Delta_mn file was not created
      
        open(unit=IO_LAN,          & 
             file=dmn_filename,    &
             status='new',         &
             access='sequential',  &
             form='unformatted')

      case (.TRUE.)       ! Delta_mn file is there.

        open(unit=IO_LAN,          & 
             file=dmn_filename,    &
             status='old',         &
             access='sequential',  &
             form='unformatted',   &
             position='append')

      end select

!  write relevant data in this symmetry to the DNM file
!  symmetry, # of analyzed vectors, matrix dimension

      write(IO_LAN) ibeta, iwwcount, no_sum

!__________________________________________________________________
!|
!|  Loop over all vectors in the 2h/2h block possessing high enough 
!|  pole strength

      DO ivcount=1,nocc 

        IF(polestr(ivcount).lt.polethr) cycle ! too small, no POPANA

        xmna = 0.0d0;  xmn2 = 0.0d0; xmnc = 0.0d0; xmnd = 0.0d0

        write(lupri,*) '   ** doing analysis for vector:',ivcount
        write(lupri,*) '   ** eigenvalue: ',cevl(ivcount)*autoev

!
!  found a vector, construct x_mn matrix
!
        ipoint = 0
        DO 14 JRP = 1, NREP
          IRP = MULTB(JRP,ibeta+NREP,2)
          IF (IRP.LT.JRP) GOTO 14
          DO J = 1, NO(JRP)
            IMIN = 1
            IF (IRP.EQ.JRP) IMIN = J + 1
            DO I = IMIN, NO(IRP)
              ipoint = ipoint + 1   ! get next eigenvector component
              irow   = no_offs(irp) + I
              jcol   = no_offs(jrp) + J
              xmna(irow,jcol) = thvecs(ipoint,ivcount)
!             write(*,'(A,F10.6,2I5)') '   read in components:', &
!       thvecs(ipoint,ivcount),irow,jcol
!             write(*,'(A,4I4,A,2I4,F10.6)') & 
!                   ' irp,jrp,i,j',irp,jrp,i,j, &
!                   ' row/col/entry:',irow,jcol,xmna(irow,jcol)
            ENDDO
          ENDDO 
  14    CONTINUE

        IF( (ipoint).NE.NOOT(ibeta)) THEN
          WRITE(lupri,*) 'ipoint,noot',ipoint,noot(ibeta)
          CALL QUIT('Error in TWOHPOP for 2h configurations!')
        ENDIF

!       write(lupri,*) 'Contributions from xmna:'
!       do i=1,no_sum
!         write(lupri,'(40F10.5)') (xmna(i,j),j=1,no_sum)
!       enddo

!  contract x_ml and x*_nl. XGEMM can be used because all irrelevant entries are
!  filled with zeros and will not contribute to the matrix product.
!  Second call: contract x_lm and x_ln. Transpose the first matrix.

! first part: contract x_ml . x_ln* over l<m,n
! result goes into xmnc
! 
        xmn2 = transpose(xmna)  ! form x_ln from x_nl. 
        CALL XGEMM ('N','N',no_sum,no_sum,no_sum,A1,xmna,no_sum, &
                    xmn2,no_sum,A0,xmnc,no_sum)

!       write(lupri,*) 'Contributions from xmnc (real):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (REAL(xmnc(i,j)),j=1,no_sum)
!       enddo
!
! second part: form x_lm from x_ml (in xmna) and then
! form x_ln* result goes into xmnd
!
        xmn2 = transpose(xmna) 
        CALL XGEMM ('N','N',no_sum,no_sum,no_sum,A1,xmn2,no_sum, &
                    xmna,no_sum,A0,xmnd,no_sum)
        xmnd = transpose(xmnd)
        xmn2 = xmnc + xmnd

!       write(lupri,*) 'Contributions from part 1+2 (real only):'
!       do i=1,no_sum
!         write(lupri,'(40F10.5)') (xmn2(i,j),j=1,no_sum)
!       enddo
!       stop

!  now form Delta_mn matrix. Part 1 is just one on the diagonal,
!  part 2 is just minus the pole strength on the dagonal.
        
        dmn = 0.0d0
        DO i=1,no_sum
          dmn(i,i) = 1.0d0 - polestr(ivcount)
!         write(*,*) 'real polestrength:',polestr(ivcount)
        ENDDO
        dmn = dmn + xmn2

        write(lupri,*) 'Delta_mn:'
        do i=1,no_sum
          write(lupri,'(40F10.5)') (dmn(i,j),j=1,no_sum)
        enddo

!  write delta_mn matrix to file in the format requested by dirac conversion
!  routine.

        dmnext=0.0d0
        dmnext(:,:,1)=dmn(:,:)
        dmn=0.0d0
        dmnext(:,:,2)=dmn(:,:)
        write(io_lan) dmnext

      ENDDO   !ivcount

!  perform all deallocations in reverse order of the allocations
!  and close file

!     write(lupri,*) 'entering deallocation section'

      deallocate(xmn2)
      deallocate(xmnd)
      deallocate(xmnc)
      deallocate(xmna)
      deallocate(dmnext)
      deallocate(dmn)
      deallocate(no_offs)
      deallocate(polestr)
      deallocate(cevl)
      deallocate(evec)
      deallocate(thvecs)

!     write(lupri,*) 'deallocation finished'

      close(io_lan)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine TWOHPOP_C(lupri,io_lan,ibeta,nrep,multb,no,noot)
!
!  Purpose: Prepare the propagator contributions to the reduced density matrix.
!
!  the Delta_mn matrix is then written to a file 'DADCPOP'
!  in the format
!              integer dimension
!              real(:,:) dmn   :  real part, vector 1
!              real(:,:) dmn   :  imaginary part, vector 1
!                       ....
!              real(:,:) dmn   :  real part, vector N
!              real(:,:) dmn   :  imaginary part, vector N
!
      implicit none
      integer, intent(in)  :: lupri    ! file handle for printing to stdout
      integer, intent(in)  :: io_lan   ! file handle start
      integer, intent(in)  :: ibeta    ! symmetry of the two-hole final state 
      integer, intent(in)  :: nrep     ! number of ireps in the overall run
      integer, intent(in)  :: multb(64,64,2)     ! fixed dimension multiplication table
      integer, intent(in)  :: no(nrep)    ! number of occupied orbitals in each irep
      integer, intent(in)  :: noot(nrep)  ! number of OO-combinations in each irep
!
!  local variables
!
      complex*16, allocatable   :: dmn(:,:),        &
                                   xmna(:,:),       &
                                   xmn2(:,:),       &
                                   xmnc(:,:),       &
                                   xmnd(:,:)
      complex*16, allocatable   :: thvecs(:,:), evec(:)
      real*8, allocatable       :: cevl(:), polestr(:), dmnext(:,:,:)
      real*8                    :: dot_productz      ! is a function
      real*8                    :: polethr = 0.001d0
      real*8                    :: autoev = 27.2113957D0
      complex*16                :: A0=(0.0d0,0.0d0), A1=(1.0d0,0.0d0)
      integer, allocatable      :: no_offs(:)
      integer                   :: i,j,ix,ivcount,iwwcount
      integer                   :: n1,n2
      integer                   :: nocc,ipoint,no_sum
      integer                   :: irp,jrp,imin,irow,jcol
      logical                   :: isthere
      character*7               :: dmn_filename='DADCPOP'

!  say hello and check basic information

      call PST('Performing 2H population analysis (complex vectors)+')
      nocc=noot(ibeta)
      write(lupri,*) 'Final state symmetry:',ibeta
      write(lupri,*) 'Number of 2h states in this symm:',nocc

!__________________________________________________________________
!|
!|  read the file of Lanczos eigenvectors and extract the 2h/2h block
!|        ----->   in this symmetry !!  <-----
!|  from it. 
!|  Fill the corresponding arrays with the read vectors and determine
!|  the polestrengths. The raw Lanc eigenvector goes to *evec* and
!|  the 2h/2h block goes to *thvecs*

      Inquire(file='TMATEVC',exist=isthere)
      IF(.not.isthere) THEN
        write(lupri,*) '**************'
        write(lupri,*) '**** Attention'
        write(lupri,*) '**************'
        write(lupri,*) 'TMATEVC file not present.'
        write(lupri,*) 'No population analysis possible'
        return
      ENDIF

      OPEN(IO_LAN,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='OLD')
      REWIND(IO_LAN)
      READ(IO_LAN) N1
      READ(IO_LAN) N2
      WRITE(lupri,*) 'Found',N1,' Lanczos eigenvectors of length', &
                   N2,' on TMATEVC'
!     WRITE(lupri,*) nocc,' of them are 2H states!'

      allocate(thvecs(nocc,nocc), &
               evec(n2),          &
               cevl(nocc),        &
               polestr(nocc))

      thvecs = A0; evec = A0; cevl = 0.0d0; polestr = 0.0d0

      iwwcount = 0
      DO I=1,nocc
        READ(IO_LAN) cevl(i)
        READ(IO_LAN) (evec(ix),ix=1,n2)
!       write(*,*) 'cevl,evec:',cevl,'\\',evec
!       pause
        polestr(i)=dot_productz(evec,nocc)
        if(polestr(i).gt.polethr) then
              write(lupri,'(A,I5,A,2F16.6,A)')           &
                  'Taking vector ',i,                    &
                  ' with', cevl(i)*autoev,polestr(i),    &
                  ' for population analysis'
              iwwcount = iwwcount + 1
              thvecs(:,i) = evec(1:nocc)   ! copy corresponding components
!             write(*,*) 'corr ev:',evec(1:nocc)
        endif 
      ENDDO
      write(lupri,*) iwwcount, &
        ' vectors have enough pole strength for 2H analysis'
      CLOSE(IO_LAN,STATUS='DELETE')
!|
!|
!|
!|
!|_______ 2h/2h block in memory together with the polestrengths ______


!  prepare occupied orbital offset table

      allocate(no_offs(nrep+1))
      no_offs = 0; no_sum = 0
      DO i=1,nrep
        no_offs(i+1) = no_offs(i) + no(i)
        no_sum = no_sum + no(i)
      ENDDO
!     write(lupri,*) 'offset table:',no_offs
!     write(lupri,*) 'total number of occupied orbitals:',no_sum

!  allocate space for Delta_mn matrix and clear matrix

      allocate(dmn(no_sum,no_sum))
      allocate(dmnext(no_sum,no_sum,2))
      allocate(xmna(no_sum,no_sum))
      allocate(xmnc(no_sum,no_sum))
      allocate(xmnd(no_sum,no_sum))
      allocate(xmn2(no_sum,no_sum))

!  open DADCPOP file for sequential access and write 
!  number of matrices and corresponding matrix dimension
!
!  ATT We check if file exists already and append the next symmetry
!  to it. This liberates us from writing symmetry-specific files altogether.
 
      inquire(file=dmn_filename,exist=isthere)

      select case(isthere)

      case (.FALSE.)      ! Delta_mn file was not created
      
        open(unit=IO_LAN,          & 
             file=dmn_filename,    &
             status='new',         &
             access='sequential',  &
             form='unformatted')

      case (.TRUE.)       ! Delta_mn file is there.

        open(unit=IO_LAN,          & 
             file=dmn_filename,    &
             status='old',         &
             access='sequential',  &
             form='unformatted',   &
             position='append')

      end select

!  write relevant data in this symmetry to the DNM file
!  symmetry, # of analyzed vectors, matrix dimension

      write(IO_LAN) ibeta, iwwcount, no_sum

!__________________________________________________________________
!|
!|  Loop over all vectors in the 2h/2h block possessing high enough 
!|  pole strength

      DO ivcount=1,nocc 

        IF(polestr(ivcount).lt.polethr) cycle ! too small, no POPANA

        xmna = A0; xmn2 = A0; xmnc = A0; xmnd = A0

        write(lupri,*) '   ** doing analysis for vector:',ivcount
        write(lupri,*) '   ** eigenvalue: ',cevl(ivcount)*autoev

!
!  found a vector, construct x_mn matrix
!
        ipoint = 0
        DO 14 JRP = 1, NREP
          IRP = MULTB(JRP,ibeta+NREP,2)
          IF (IRP.LT.JRP) GOTO 14
          DO J = 1, NO(JRP)
            IMIN = 1
            IF (IRP.EQ.JRP) IMIN = J + 1
            DO I = IMIN, NO(IRP)
              ipoint = ipoint + 1   ! get next eigenvector component
              irow   = no_offs(irp) + I
              jcol   = no_offs(jrp) + J
              xmna(irow,jcol) = thvecs(ipoint,ivcount)
!             write(*,*) '   read in components:',thvecs(ipoint,ivcount)
!             write(*,*) 'check: irow,icol:',irow,jcol ! irow must always be greater than jcol
!             write(*,'(A,4I4,A,2I4,F10.6)') & 
!                   ' irp,jrp,i,j',irp,jrp,i,j, &
!                   ' row/col/entry:',irow,jcol,xmna(irow,jcol)
            ENDDO
          ENDDO 
  14    CONTINUE

        IF( (ipoint).NE.NOOT(ibeta)) THEN
          WRITE(lupri,*) 'ipoint,noot',ipoint,noot(ibeta)
          CALL QUIT('Error in TWOHPOP for 2h configurations!')
        ENDIF

!       write(lupri,*) 'Contributions from xmna (real):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (REAL(xmna(i,j)),j=1,no_sum)
!       enddo
!       write(lupri,*) 'Contributions from xmna (imag):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (AIMAG(xmna(i,j)),j=1,no_sum)
!       enddo

!  contract x_ml and x*_nl. XGEMM can be used because all irrelevant entries are
!  filled with zeros and will not contribute to the matrix product.
!  Second call: contract x_lm and x_ln. Transpose the first matrix.

! first part: contract x_ml . x_ln* over l<m,n
! result goes into xmnc

        xmn2 = transpose(xmna)  ! form x_ln from x_nl. 
        xmn2 = dconjg(xmn2)     ! form x_ln* from x_nl. 
        CALL XGEMM ('N','N',no_sum,no_sum,no_sum,A1,xmna,no_sum, &
                    xmn2,no_sum,A0,xmnc,no_sum)

!       write(lupri,*) 'Contributions from xmnc (real):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (REAL(xmnc(i,j)),j=1,no_sum)
!       enddo
!
! second part: form x_lm from x_ml (in xmna) and then
! form x_ln* result goes into xmnd
!
        xmn2 = dconjg(transpose(xmna))  ! form x_ln* from x_ml. 
        CALL XGEMM ('N','N',no_sum,no_sum,no_sum,A1,xmn2,no_sum, &
                    xmna,no_sum,A0,xmnd,no_sum)
        xmnd = transpose(xmnd)
        xmn2 = xmnc + xmnd
!
!
!       write(lupri,*) 'Contributions from part1+2 (real):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (REAL(xmn2(i,j)),j=1,no_sum)
!       enddo
!       write(lupri,*) 'Contributions from part1+2 (imag):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (AIMAG(xmn2(i,j)),j=1,no_sum)
!       enddo
!       stop


!  now form Delta_mn matrix. Part 1 is just one on the diagonal,
!  part 2 is just minus the pole strength on the dagonal.
        
        dmn = 0.0d0
        DO i=1,no_sum
          dmn(i,i) = A1 - DCMPLX(polestr(ivcount),0.0d0)
        ENDDO
        dmn = dmn + xmn2

        write(lupri,*) 'Delta_mn (real part):'
        do i=1,no_sum
          write(lupri,'(40F10.5)') (REAL(dmn(i,j)),j=1,no_sum)
        enddo
        write(lupri,*) 'Delta_mn (imag part):'
        do i=1,no_sum
          write(lupri,'(40F10.5)') (AIMAG(dmn(i,j)),j=1,no_sum)
        enddo
        write(lupri,*)
        write(lupri,'(A)') '------------------------------------------'
        write(lupri,*)

!  write delta_mn matrix to file in the format requested by dirac conversion
!  routine. dmnext can be a real array since real part is first layer and 
!  imaginary part is second layer written out.

        dmnext=0.0d0
        dmnext(:,:,1)=real(dmn(:,:))
        dmnext(:,:,2)=aimag(dmn(:,:))
        write(io_lan) dmnext

      ENDDO   !ivcount

!  perform all deallocations in reverse order of the allocations
!  and close file

!     write(lupri,*) 'entering deallocation section'

      deallocate(xmn2)
      deallocate(xmnd)
      deallocate(xmnc)
      deallocate(xmna)
      deallocate(dmnext)
      deallocate(dmn)
      deallocate(no_offs)
      deallocate(polestr)
      deallocate(cevl)
      deallocate(evec)
      deallocate(thvecs)

!     write(lupri,*) 'deallocation finished'

      close(io_lan)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine DNMMATCHK(lupri,io_lan)
!
      implicit none
      integer :: lupri    ! file handle for printing to stdout
      integer :: io_lan   ! file handle start
      logical :: errflg   ! error flag
!
!  local variables
!
      real*8,  allocatable  :: dmnext(:,:,:)
      integer               :: i,j,ivcount,ibeta,iwwcount,no_sum
      character*7           :: dmn_filename='DADCPOP'

      call PST('Performing DMN matrix check+')

!  open DNM matrix and read in everything

      open(unit=io_lan,             &
           file=dmn_filename,       &
           status='old',            &
           access='sequential',     &
           form='unformatted')

      errflg=.false.
 12   read(io_lan,end=800,err=900) ibeta,iwwcount,no_sum
!     write(lupri,*) 'Record:',ibeta,iwwcount,no_sum
      allocate(dmnext(no_sum,no_sum,2))
      dmnext = 0.0d0
      do ivcount=1,iwwcount
        read(io_lan) dmnext
      enddo

!     write(lupri,*) 'Delta_mn (real):'
!     do i=1,no_sum
!       write(lupri,'(40F10.5)') (dmnext(i,j,1),j=1,no_sum)
!     enddo
!     write(lupri,*) 'Delta_mn (complex):'
!     do i=1,no_sum
!       write(lupri,'(40F10.5)') (dmnext(i,j,2),j=1,no_sum)
!     enddo

      write(lupri,*) 'Symmetry',ibeta,' passed.'
      deallocate(dmnext)
      goto 12

 800  return

 900  errflg=.false.
      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION DOT_PRODUCTZ(a,n)
      real*8 dot_productz
!
!  This little function computes the squared norm of a *complex*
!  vector according to
!        Dot_productz = conjg(a) * a
!  and therefore yields a real value.
!
      integer,    intent(in)    ::  n
      complex*16, intent(in)    ::  a(n)
!
!  local variables
!
      complex*16 :: ca1
      real*8     :: rsum,thrimag=1.0D-12
      integer    :: i


      rsum=0.0d0
      DO i=1,n
        ca1 = conjg(a(i))*a(i)
        if(aimag(ca1).gt.thrimag) then
          stop 'unreasonably high imaginary part in DOT_PRODUCTZ !'
        endif
        rsum = rsum + dble(ca1)
      ENDDO

      dot_productz = rsum
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

