C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C
C      ****  Diagonalization Routines for the RELADC package  ****
C
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FULLDIAR(LUN_DIA,LUN_OFF,INTBUF,LADC,A,EVL,WORK,LWORK,
     &                   IOI,IOJ,BUF)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LUN_DIA,LUN_OFF,INTBUF,LADC,LWORK
      REAL*8 A(LADC,LADC),WORK(LWORK),EVL(*),BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)
C
C---------------Description--------------------------------------------
C
C  For not too large Hermitian matrices the full diagonalization is
C  started.
C
C---------------Local variables--------------------------------------
C
      REAL*8 S
      CHARACTER*1 JOBZ,UPLO

      REAL*8 AUTOEV
      PARAMETER(AUTOEV = 27.2113957D0)
C
C---------------Executable code--------------------------------------
C
C
C  construct ADC matrix by reading nonzero elements from file
C
      DO I=1,LADC
      DO J=1,LADC
        A(J,I)=0.0D0
      ENDDO
      ENDDO
      REWIND(LUN_DIA)
      READ(LUN_DIA,ERR=88,END=99) (A(IXX,IXX),IXX=1,LADC),NBUFS
      REWIND(LUN_OFF)
      TOTELE1 = 0
      DO IREC = 1,NBUFS
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                        (IOI(IXX),IXX=1,INTBUF),
     &                        (IOJ(IXX),IXX=1,INTBUF),
     &                         NACT,JDUMMY
C       WRITE(*,*) ' ## REC:',irec,' ## NACT:',nact,' ## JDUMMY:',jdummy
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          S    = BUF(K)

          IF(ICOL.GE.IROW) CALL QUIT('ICOL.GE.IROW should not happen')
 
          A(IROW,ICOL)=S
          TOTELE1=TOTELE1 + 1
        ENDDO
      ENDDO
      WRITE(*,555) TOTELE1
C
C  call symmetric diagonalizer
C  calculate eigenvectors being placed in A
C
      JOBZ='V'
      UPLO='L'

      CALL DSYEV( JOBZ, UPLO, LADC, A, LADC, EVL, WORK, LWORK, INFO )

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in DSYEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAR (lapack dsyev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues in eV (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO


      RETURN

 88   STOP 'Error reading the diagonal (FULLDIAR)'
 99   STOP 'Unexpected end of diagonal file (FULLDIAR)'
 199  STOP 'Error reading off-diagonal elements (FULLDIAR)'

 555  FORMAT (1X,'Total number of elements in OD-file:',T40,I8)

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FULLDIAC(LUN_DIA,LUN_OFF,INTBUF,LADC,A,EVL,CWORK,RWORK,
     &                    LWORK,IOI,IOJ,BUF)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LUN_DIA,LUN_OFF,INTBUF,LADC,LWORK
      COMPLEX*16 A(LADC,LADC),CWORK(*)
      REAL*8 EVL(*),RWORK(LWORK)
      INTEGER IOI(INTBUF),IOJ(INTBUF)
      COMPLEX*16 BUF(INTBUF)
C
C---------------Description--------------------------------------------
C
C  For not too large Hermitian matrices the full diagonalization is
C  started.
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "../relccsd/files.inc"
      INCLUDE "../relccsd/complex.inc"
      INCLUDE "../relccsd/param.inc"
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 S
      CHARACTER*1 JOBZ,UPLO

      REAL*8 AUTOEV
      PARAMETER(AUTOEV = 27.2113957D0)
C
C---------------Executable code--------------------------------------
C
C
C  construct complex matrix by reading nonzero elements from file
C
      DO I=1,LADC
      DO J=1,LADC
        A(J,I)=A0
      ENDDO
      ENDDO
      REWIND(LUN_DIA)
      READ(LUN_DIA,ERR=88,END=99) (A(IXX,IXX),IXX=1,LADC),NBUFS
c     write(iw,*) (A(IXX,IXX),IXX=1,LADC)
      REWIND(LUN_OFF)
      TOTELE1 = 0
      DO IREC = 1,NBUFS
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                        (IOI(IXX),IXX=1,INTBUF),
     &                        (IOJ(IXX),IXX=1,INTBUF),
     &                         NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          S    = BUF(K)

          IF(ICOL.GE.IROW) STOP 'ICOL.GE.IROW should not happen'
 
          A(IROW,ICOL)=S
          TOTELE1=TOTELE1 + 1
        ENDDO
      ENDDO
      WRITE(*,555) TOTELE1
C
C  call hermitian diagonalizer, no eigenvectors to be calculated
c  ** ATT ** LWORK is here the size of CWORK and NOT of RWORK!
c  therefore we set this variable to the same value as in the calling
c  program!
C
      JOBZ='V'
      UPLO='L'
      LWORK = (LADC+1)*LADC

      CALL ZHEEV( JOBZ, UPLO, LADC, A, LADC, EVL, CWORK, LWORK, RWORK,
     &                  INFO )
      IF(INFO.NE.0)
     &    CALL QUIT ('*** Error in FULLDIAC (lapack zheev) ***')

      WRITE(*,*) 'Eigenvalues in au/eV (lowest 20):'
      DO I=1,min(20,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(IW,*) 'Corresponding eigenvectors (first 100 entries)'
      do iev=1,min(20,ladc)
        WRITE(IW,*) iev,EVL(iev)
        WRITE(IW,*) '----------------------'
        do i=1,min(100,ladc)
          write(iw,*) i,a(i,iev)
        enddo
        write(iw,*)
        write(iw,*)
      enddo
 
      RETURN

 88   CALL QUIT('Error reading the diagonal (FULLDIAC)')
 99   CALL QUIT('Unexpected end of diagonal file (FULLDIAC)')
 199  CALL QUIT('Error reading off-diagonal elements (FULLDIAC)')
 555  FORMAT (1X,'Total number of elements in Offdiag-file:',T45,I8)

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RVECANALYZE(IOUNIT,IOUNIT2,DAFILE,KREP,LADC,
     &                       NEIGENV,EVL,EVC,THR)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Analyzes eigenvectors of the ADC matrix
C     Writes pole strengths and corresponding 1p-2h configuration
C     to stdout if eigenvector entry is larger than THR.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*9 DAFILE
      INTEGER IOUNIT,IOUNIT2,KREP,LADC,NEIGENV
      REAL*8  EVL(LADC),THR
      REAL*8 EVC(LADC,NEIGENV)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "../relccsd/symm.inc"
      INCLUDE "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 AUEV, WH, WP
      DATA AUEV /27.2113957D0/
      CHARACTER*7 EVECFN
      REAL*8 B
      character*78 field
C
C---------------Executable code--------------------------------------
C
      WRITE(IW,800) KREP,THR
 800  FORMAT (//" ==================================================="/
     &        " ===    Eigenvector analysis for symmetry ",I2/
     &        " ===    Configuration component threshold ",F10.4/
     &        " ==================================================="//)
c    &        " ===    From minimum eigenvalue of ",
c    &           F10.4,' (Hartree)',3X,F10.4,' (eV)'/
c    &        " ===    to maximum eigenvalue of   ",
c    &           F10.4,' (Hartree)',3X,F10.4,' (eV)'/

c
c  open eigenvector file
c
      write(evecfn,'(a5,i1)') 'EVEC.',krep
      open(iounit, file=EVECFN, form='UNFORMATTED',
     &     access='SEQUENTIAL', status='UNKNOWN')
      read(iounit) iladc,neigenv
      if(iladc.ne.ladc) then
        call quit('Incompatible eigenvector files!')
      endif
      write(iw,*) 'There are',neigenv,' eigenvectors on file.'
      if(neigenv.gt.ladc) then
        call quit('More eigenvectors than matrix dimension ?')
      endif
c
c  reading eigenvectors for eigenvalues
c
      do j=1,neigenv
        read(iounit) evl(j) 
        read(iounit) (evc(ix,j),ix=1,ladc)
      enddo
      close(iounit)
c
c  open direct access file with configuration data
c
      WRITE(IW,*) 'ADC configuration data read from ',DAFILE
      IRECL = 78
      OPEN(iounit2,FILE=DAFILE,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
c
c  starting analysis
c
      DO 300 I = 1,NEIGENV

        WRITE(iw,*)
        WRITE(iw,*) '----------------------------------------'
        WRITE(iw,*) 'Eigenvalue ',I,':',EVL(I),' (a.u.)',
     &               EVL(I)*AUEV,' (eV)'
C
C  determine weight of hole and particle space
C
        WH=0.0D0
        WP=0.0D0
        DO J=1,LADC
            B=EVC(J,I)
            IF(J.LE.NO(KREP)) THEN
              WH=WH + B*B
            ELSE
              WP=WP + B*B
            ENDIF
        ENDDO
        WRITE(iw,860) WH*100.0D0,WP*100.0D0
        WRITE(iw,*) '----------------------------------------'

        DO 150 J=1,LADC
          B=EVC(J,I)
          B = B*B
          IF(DABS(B).GE.THR) THEN
            READ(iounit2,REC=J) field
            WRITE(iw,'(A5,I6,A2,A78)') 'Conf:',J,' |',field
          ENDIF
 150    CONTINUE

 300  CONTINUE

      CLOSE(iounit2) 

      RETURN
C
C  formats...
C:
 850  FORMAT (1X,'number of eigenpairs on output file : ',I4,/,
     &        1X,'eigenpairs analysis requested for   : ',I4,' roots',/,
     &        1X,'analysis will be performed    for   : ',I4,' roots')
 860  FORMAT (/" Weight of 1h-space:",T25,F6.1,' %'/
     &         " Weight of 1p-2h-space:",T25,F6.1,' %')
 880  FORMAT (/" 1h space"//"   Configuration ",I3,2X,"Spinor ",I3,
     &        " Irep ",A4," Coeff ",F10.4/)
 890  FORMAT (/" 1p-2h space"//"   Configuration ",I4,"   (",
     &        I2,A4,3X,I2,A4,3X,I2,A4,")"," Coeff ",F10.4)
 1100 FORMAT (1X,/,1X,49('-'),/,1X,'# ',I4,' )',1X,
     &        F12.6,' (eV)',1X,F12.6,' (Hartree)',/,1X,49('-'),/)

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZVECANALYZE(IOUNIT,IOUNIT2,DAFILE,KREP,LADC,
     &                       NEIGENV,EVL,EVC,THR)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Analyzes eigenvectors of the ADC matrix
C     Writes pole strengths and corresponding 1p-2h configuration
C     to stdout if eigenvector entry is larger than THR.
C
C---------------Last modified------------------------------------------
C
C     Author : MP
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*9 DAFILE
      INTEGER IOUNIT,IOUNIT2,KREP,LADC,NEIGENV
      REAL*8  EVL(LADC),THR
      COMPLEX*16 EVC(LADC,NEIGENV)
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "../relccsd/symm.inc"
      INCLUDE "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 AUEV, WH, WP
      DATA AUEV /27.2113957D0/
      CHARACTER*7 EVECFN
      COMPLEX*16 B
      character*78 field
C
C---------------Executable code--------------------------------------
C
      WRITE(IW,800) KREP,THR
 800  FORMAT (//" ==================================================="/
     &        " ===    Eigenvector analysis for symmetry ",I2/
     &        " ===    Configuration component threshold ",F10.4/
     &        " ==================================================="//)
c    &        " ===    From minimum eigenvalue of ",
c    &           F10.4,' (Hartree)',3X,F10.4,' (eV)'/
c    &        " ===    to maximum eigenvalue of   ",
c    &           F10.4,' (Hartree)',3X,F10.4,' (eV)'/

c
c  open eigenvector file
c
      write(evecfn,'(a5,i1)') 'EVEC.',krep
      open(iounit, file=EVECFN, form='UNFORMATTED',
     &     access='SEQUENTIAL', status='UNKNOWN')
      read(iounit) iladc,neigenv
      if(iladc.ne.ladc) then
        call quit('Incompatible eigenvector files!')
      endif
      write(iw,*) 'There are',neigenv,' eigenvectors on file.'
      if(neigenv.gt.ladc) then
        call quit('More eigenvectors than matrix dimension ?')
      endif
c
c  reading eigenvectors for eigenvalues
c
      do j=1,neigenv
        read(iounit) evl(j) 
        read(iounit) (evc(ix,j),ix=1,ladc)
      enddo
      close(iounit)
c
c  open direct access file with configuration data
c
      WRITE(IW,*) 'ADC configuration data read from ',DAFILE
      IRECL = 78
      OPEN(iounit2,FILE=DAFILE,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
c
c  starting analysis
c
      DO 300 I = 1,NEIGENV

        WRITE(iw,*)
        WRITE(iw,*) '----------------------------------------'
        WRITE(iw,*) 'Eigenvalue ',I,':',EVL(I),' (a.u.)',
     &               EVL(I)*AUEV,' (eV)'
C
C  determine weight of hole and particle space
C
        WH=0.0D0
        WP=0.0D0
        DO J=1,LADC
            B=EVC(J,I)
            IF(J.LE.NO(KREP)) THEN
              WH=WH + dble(B*dconjg(B))
            ELSE
              WP=WP + dble(B*dconjg(B))
            ENDIF
        ENDDO
        WRITE(iw,860) WH*100.0D0,WP*100.0D0
        WRITE(iw,*) '----------------------------------------'

        DO 150 J=1,LADC
          B=EVC(J,I)
          B = B*dconjg(B)
          IF(ABS(B).GE.THR) THEN
            READ(iounit2,REC=J) field
            WRITE(iw,'(A5,I6,A2,A78)') 'Conf:',J,' |',field
          ENDIF
 150    CONTINUE

 300  CONTINUE

      CLOSE(iounit2) 

      RETURN
C
C  formats...
C:
 850  FORMAT (1X,'number of eigenpairs on output file : ',I4,/,
     &        1X,'eigenpairs analysis requested for   : ',I4,' roots',/,
     &        1X,'analysis will be performed    for   : ',I4,' roots')
 860  FORMAT (/" Weight of 1h-space:",T25,F6.1,' %'/
     &         " Weight of 1p-2h-space:",T25,F6.1,' %')
 880  FORMAT (/" 1h space"//"   Configuration ",I3,2X,"Spinor ",I3,
     &        " Irep ",A4," Coeff ",F10.4/)
 890  FORMAT (/" 1p-2h space"//"   Configuration ",I4,"   (",
     &        I2,A4,3X,I2,A4,3X,I2,A4,")"," Coeff ",F10.4)
 1100 FORMAT (1X,/,1X,49('-'),/,1X,'# ',I4,' )',1X,
     &        F12.6,' (eV)',1X,F12.6,' (Hartree)',/,1X,49('-'),/)

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RLANCZOS(IRREP, LADC, NP, NITER, NEIGENV, IVX,
     &                    P, T, R, DIAG, Z, AEV, BUF, OI, OJ, CEVL, 
     &                    CFLAG, INTBUF, IO_DIA, IO_OFF, IO_LAN)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER IRREP,LADC,NP,NITER,NEIGENV,IVX(2*NP+1)
      INTEGER INTBUF,IO_DIA,IO_OFF,IO_LAN
      REAL*8 P(LADC,2*NP+1), T(NP+1,NITER), R(LADC), DIAG(LADC)
      REAL*8 Z(NITER,NITER), AEV(LADC,NEIGENV), BUF(INTBUF)
      REAL*8 CEVL(NITER)
      INTEGER OI(INTBUF),OJ(INTBUF)
      LOGICAL CFLAG
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
      INCLUDE "../relccsd/files.inc"
      INCLUDE "../relccsd/complex.inc"
C
C---------------Local variables and function definitions ------------
C
      REAL*8 Y,DDOT,Z1
      CHARACTER*7 EVECFN
C
C---------------Executable code--------------------------------------
C
      WRITE(EVECFN,'(A5,I1)') 'EVEC.',IRREP

      CALL PST('Executing real symmetric Band-Lanczos+')
      WRITE(IW,*)
      WRITE(IW,*) '*** Symmetry ',IRREP
      WRITE(IW,*)
      WRITE(IW,*) ' Matrix dimension:',LADC
      WRITE(IW,*) ' Band width:      ',NP
      WRITE(IW,*) ' Number of iterations:',NITER
      WRITE(IW,*) ' Unit numbers dia,off:',IO_DIA,IO_OFF
      WRITE(IW,*) ' Eigenvector run:     ',CFLAG
      WRITE(IW,*) ' Eigenvector file:    ',EVECFN
      WRITE(IW,*)
c
c  ... clear approx eigenvector array if requested
c
      IF(CFLAG) THEN
        DO I=1,NEIGENV
          DO J=1,LADC
            AEV(J,I)=0.0D0
          ENDDO
        ENDDO
      ENDIF
c
c  ... read in purged eigenvector matrix ZVEC for calculation of
c      the long eigenvectors.
c
      if(cflag) then
        open(io_lan,filE='ZVEC',form='UNFORMATTED',status='UNKNOWN')
        rewind(io_lan)
        read(io_lan) n1
        read(io_lan) n2
        write(iw,*) 'There are',n1,' Lanczos eigenvectors available.'
        write(iw,*) 'Reading ZVEC ...'
        if(n1.ne.neigenv.or.n2.ne.niter)
     &    call quit('incompatible eigenvector parameters!')
        do i=1,n1
          read(io_lan) 
          read(io_lan) (z(ix,i),ix=1,n2)
        enddo
        close(io_lan,status='DELETE')
      endif

c
c  ... clear Lanczos matrix
c
      DO I=1,NITER
        DO J=1,NP+1
          T(J,I)=0.0D0
        ENDDO
      ENDDO
c
c  ... clear r/p-vector array and initialize start vectors
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

c----------------------------- lanczos loop start ----------------------

      DO 99 IL = 1,NITER
c
c
      IF(mod(il,50).eq.0) THEN
        WRITE(IW,*) '---  Band Lanczos iteration ',IL,' of',niter
      ENDIF
c
c  form r = Hq_i - sum t_ji*q_j
c  the j index in the sum runs from i-p up to i-1
c  and the loop is only executed if previous vectors exist.
c  the actual vector to be multiplied by A always sits in
c  P(*,np+1)
c
        CALL RMATBLK(LADC, P(1,IVX(NP+1)), R, DIAG, BUF,
     &               OI, OJ, INTBUF, IO_DIA, IO_OFF)

c  subtract all the previous components

        ICC = 1
        DO IP = MAX0(IL-NP,1), IL-1
c     write(iw,*) ' ++++ executing ORTHO LOOP +++++++',IL,ip
          ITROW = ICC+1
          ITCOL = IL - ICC
c     write(iw,*) ' ++++ itrow,itcol,np-icc+1:',itrow,itcol,np-icc+1
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
c  now write out Lanczos matrix when called for the FIRST time !
c
      IF(.NOT.CFLAG) THEN
c       write(iw,*) ' writing out Lanczos matrix...'
        open(io_lan,file='TMAT',form='UNFORMATTED',access='SEQUENTIAL',
     &              status='UNKNOWN')
        write(io_lan) np+1,niter
        do j=1,niter
          write(io_lan) (t(i,j),i=1,np+1)
        enddo
        close(io_lan)
        write(iw,*) ' Lanczos matrix written.'
      ENDIF
c
c  write out long ADC eigenvectors when eigenvector run finished
c
      IF(CFLAG) THEN
        open(io_lan, file=EVECFN, form='UNFORMATTED', 
     &       access='SEQUENTIAL', status='UNKNOWN')
        write(io_lan) LADC,NEIGENV
        do j=1,neigenv
          write(io_lan) cevl(j)
          write(io_lan) (aev(i,j),i=1,ladc)
        enddo
        close(io_lan)
        write(iw,*) ' ADC eigenvectors written.'

      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZLANCZOS(IRREP, LADC, NP, NITER, NEIGENV, IVX,
     &                    P, T, R, DIAG, Z, AEV, BUF, OI, OJ, CEVL, 
     &                    CFLAG, INTBUF, IO_DIA, IO_OFF, IO_LAN)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER IRREP,LADC,NP,NITER,NEIGENV,IVX(2*NP+1)
      INTEGER INTBUF,IO_DIA,IO_OFF,IO_LAN
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
      INCLUDE "../relccsd/files.inc"
      INCLUDE "../relccsd/complex.inc"
C
C---------------Local variables and function definitions ------------
C
      COMPLEX*16 Y,ZDOTC,Z1
      CHARACTER*7 EVECFN
C
C---------------Executable code--------------------------------------
C
      WRITE(EVECFN,'(A5,I1)') 'EVEC.',IRREP

      WRITE(IW,*)
      WRITE(IW,*) '*************************************'
      WRITE(IW,*) '*** executing Hermitian Band-Lanczos '
      WRITE(IW,*) '*** in symmetry ',IRREP
      WRITE(IW,*) '*************************************'
      WRITE(IW,*)
      WRITE(IW,*) ' Matrix dimension:',LADC
      WRITE(IW,*) ' Band width:      ',NP
      WRITE(IW,*) ' Number of iterations:',NITER
      WRITE(IW,*) ' Unit numbers dia,off:',IO_DIA,IO_OFF
      WRITE(IW,*) ' Eigenvector run:     ',CFLAG
      WRITE(IW,*) ' Eigenvector file:    ',EVECFN
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
c  ... read in purged eigenvector matrix ZVEC for calculation of
c      the long eigenvectors.
c
      if(cflag) then
        open(io_lan,filE='ZVEC',form='UNFORMATTED',status='UNKNOWN')
        rewind(io_lan)
        read(io_lan) n1
        read(io_lan) n2
        write(iw,*) 'There are',n1,' Lanczos eigenvectors available.'
        write(iw,*) 'Reading ZVEC ...'
        if(n1.ne.neigenv.or.n2.ne.niter)
     &    call quit('incompatible eigenvector parameters!')
        do i=1,n1
          read(io_lan) 
          read(io_lan) (z(ix,i),ix=1,n2)
        enddo
        close(io_lan,status='DELETE')
      endif

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
      IF(mod(il,50).eq.0) THEN
        WRITE(IW,*) '---  Band Lanczos iteration ',IL,' of',niter
      ENDIF
c
c  form r = Hq_i - sum t_ji*q_j
c  the j index in the sum runs from i-p up to i-1
c  and the loop is only executed if previous vectors exist.
c  the actual vector to be multiplied by A always sits in
c  P(*,np+1)
c
        CALL CMATBLK(LADC, P(1,IVX(NP+1)), R, DIAG, BUF,
     &               OI, OJ, INTBUF, IO_DIA, IO_OFF)

c  subtract all the previous components

        ICC = 1
        DO IP = MAX0(IL-NP,1), IL-1
c     write(iw,*) ' ++++ executing ORTHO LOOP +++++++',IL,ip
          ITROW = ICC+1
          ITCOL = IL - ICC
c     write(iw,*) ' ++++ itrow,itcol,np-icc+1:',itrow,itcol,np-icc+1
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
c       Y = ZSQRT(Y)
c  att: the PGI compiler does not know zsqrt
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
c  now write out Lanczos matrix when called for the FIRST time !
c
      IF(.NOT.CFLAG) THEN
        open(io_lan,file='TMAT',form='UNFORMATTED',access='SEQUENTIAL',
     &              status='UNKNOWN')
        write(io_lan) np+1,niter
        do j=1,niter
          write(io_lan) (t(i,j),i=1,np+1)
        enddo
        close(io_lan)
        write(iw,*) ' Lanczos matrix written.'
c  debug
c       write(iw,*) ' Lanczos matrix:'
c       do j=1,niter
c         write(*,*) (t(i,j),i=1,np+1)
c       enddo
c  debug end
      ENDIF
c
c  write out long ADC eigenvectors when eigenvector run finished
c
      IF(CFLAG) THEN
        open(io_lan, file=EVECFN, form='UNFORMATTED', 
     &       access='SEQUENTIAL', status='UNKNOWN')
        write(io_lan) LADC,NEIGENV
        do j=1,neigenv
          write(io_lan) cevl(j)
          write(io_lan) (aev(i,j),i=1,ladc)
        enddo
        close(io_lan)
        write(iw,*) ' ADC eigenvectors written.'
c
c  debug purpose: write some eigenvectors to output file
c  (gets too large !)
c
c       do i=1,min(20,neigenv)
c         write(iw,*)
c         write(iw,*) i,cevl(i)
c         write(iw,*)
c         do j=1,ladc
c           write(iw,*) j,aev(j,i)
c         enddo
c         write(iw,*) '------------------------------------------'
c       enddo

      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RINVERSIT(NDIM,NP,NP1,MM,REW,T,A,A1,B,BX,
     &                     IPVT,NINVI,SUCCESS)
c
c  this routine is based on numerical recipes
c

c     ndim: number of iterations in the lanczos. corresponds to niter
c           the true length of a T-eigenvector is niter
c     mmm:  this is just the blown-up dimension of T
c           corresponding to twice the bandwidth + 1.
c           remember: in the Numerical Recipes routine
c           the upper and lower part of the band-matrix are stored.
c    ninvi: maximum number of allowes inverse iterations. reason:
c           if an eigenvalue is converged the corresponding eigenvector
c           is immediately obtainable by 1 or two iterations. but it was
c           observed that in the case of multiple Lanczos vectors with a
c           corresp. eigenvector close to zero norm inverse iteration
c           does not lead to convergence due to numerical problems.
c           this case must be caught for stable running conditions.

      IMPLICIT INTEGER (A-Z)
c
      INTEGER NDIM, NP, NP1, MM, IPVT(NDIM), NINVI
      REAL*8 REW
      REAL*8 T(NP1,NDIM), A(NDIM,MM), A1(NDIM,MM), B(NDIM), BX(NDIM)
      logical success
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "../relccsd/files.inc"
C
C--------------- Local variables ----------------------------------
C     
      real*8 EW
      real*8 dum,y
      integer i, j,numit

      REAL*8 AUEV
      DATA AUEV/27.2113957D0/
C
C--------------- Executable code -----------------------------------
C     

c     ---------------------------------------------------
c     construct A from T according to the needs of NUMREC
c     ---------------------------------------------------

      EW=REW

      DO J = 1, MM
         DO I = 1, NDIM
            A(I,J) = 0.0D0
         ENDDO
      ENDDO
c

      success = .true.
c
      do i = 1, ndim
         a(i,np1) = t(1,i)
      enddo
      do j = 1, np
         do i = 1, ndim - j
            a(i,np1+j) = t(j+1,i)
            a(i+j,np1-j) = t(j+1,i)
         enddo
      enddo
c     ------------------------
c     form A - lambda*I
c     ------------------------
      do i = 1, ndim
         a(i,np1) = a(i,np1) - ew
      enddo
c     -----------------------------------------------
c     LU Decomposition of a - ew
c     use second vector BX for convergence control of
c     the inverse iteration
c     -----------------------------------------------
      call rbandec(a,ndim,np,np,ndim,mm,a1,mm,ipvt,dum)
c
      call rnormalize(ndim,b)

      numit = 0
 380  continue

      call dcopy(ndim,b,1,bx,1)
      call rbanbks(a,ndim,np,np,ndim,mm,a1,mm,ipvt,b)
      call rnormalize(ndim,b)
      call rnormdiff(ndim,b,bx,y)
      if(y.gt.1.0D-10) then
        numit = numit + 1
        if(numit.ge.ninvi) then
          write(iw,*) ' *** ATT ***'
          write(iw,*) ' The inverse iteration for eigenvalue',REW*AUEV,
     &                ' (eV)'
          write(iw,*) ' did not converge. The corresponding eigen-'
          write(iw,*) ' vector is not to be trusted!'
          write(iw,*) ' Final norm:',y
          write(iw,*) ' *** Remedy ***'
          write(iw,*) ' Decrease number of Lanczos iterations',
     &                ' or increase band width'
          
          success = .false.

          goto 382
        endif
        goto 380
      endif

 382  continue

      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RNORMALIZE (NDIM,B)
      implicit none
c     
      integer ndim
      real*8 b(ndim)
c     
      integer i
      real*8  y,norm
c     
      norm = 0.0d0
      do 10 i = 1, ndim
         norm = norm + b(i)*b(i)
 10   continue
c     
      y = 1.0d0 / dsqrt(norm)
      do 20 i = 1, ndim
         b(i) = b(i) * y
 20   continue
c     
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RNORMDIFF (NDIM,V,W,DIFF)
c
c  this routine does not produce the ordinary vector difference 
c  norm but just compares the absolute values of the vector
c  components. this is due to possible phase changes in the inverse
c  iteration.
c
      implicit none
c     
      integer ndim
      real*8 v(ndim),w(ndim)
      real*8 diff
c     
      integer i
      real*8 sum,d1
c     
      sum = 0.0d0
      do 10 i = 1, ndim
         d1 = abs(v(i)) - abs(w(i))
         sum = sum + d1*d1
 10   continue
      diff = dsqrt(sum)
c     
      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RBANDEC(A,N,M1,M2,NP,MP,A1,MP1,INDX,D)
c     ------------------------------------------------
c     LU decomposition of a band martix a
c     numerical recipes --> chapter 2, banddiagonal systems
c     a is destroyed and contatins the upper triangle
c     a1 contains the lower triangle
c     indx records the row permutations
c     -----------------------------------------------------
      implicit none
c
      integer m1, m2, mp, mp1, n, np, indx(n)
      real*8 tiny, d
      real*8 a(np,mp), a1(np,mp1)
c
      parameter (tiny = 1.d-20)
c
      integer i, j, k, l, mm
      real*8 dum
c
      mm = m1 + m2 + 1
      if (mm.gt.mp .or. m1.gt.mp1 .or. n.gt.np) then
         write(6,*) 'bad args in zbandec'
         stop
      endif
c
      l = m1
      do 13 i = 1, m1
         do 11 j = m1 + 2 - i, mm
            a(i,j-l) = a(i,j)
 11      continue
         l = l - 1
         do 12 j = mm - l, mm
            a(i,j) = 0.0d0
 12      continue
 13   continue
c
      d = 1.0d0
      l = m1
c
c     for each row ...
c
      do 18 k = 1, n
         dum = a(k,1)
         i = k
         if (l .lt. n) then
            l = l + 1
         endif
c
c     Find the pivot
c     
         do 14 j = k+1, l
            if (dabs(a(j,1)) .gt. dabs(dum)) then
               dum = a(j,1)
               i = j
            endif
 14      continue
c     
         indx(k) = i
         if (dabs(dum) .eq. 0.0d0) then
            a(k,1) = tiny
         endif
c     
         if (i .ne. k) then
c     Interchange rows
            d = -d
            do 15 j = 1, mm
               dum = a(k,j)
               a(k,j) = a(i,j)
               a(i,j) = dum
 15         continue
         endif
c
c     Do the elimination
c     
         do 17 i = k + 1, l
            dum = a(i,1) / a(k,1)
            a1(k,i-k) = dum
            do 16 j = 2, mm
               a(i,j-1) = a(i,j) - dum*a(k,j)
 16         continue
            a(i,mm) = 0.0d0
 17      continue
c     
 18   continue
c     
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RBANBKS(A,N,M1,M2,NP,MP,A1,MP1,INDX,B)
c     ------------------------------------------------------
c     solves Ax=b for x
c     A is a bandmatrix and has been decomposed in upper(a)
c     and lower(a1) triangle with bandec.
c     x is returned in b
c     ------------------------------------------------------
      implicit none
c
      integer m1, m2, mp, mp1, n, np, indx(n)
      real*8 a(np,mp), a1(np,mp1), b(n)
c
      integer i, k, l, mm
      real*8 dum
c
      mm = m1 + m2 + 1
      if (mm.gt.mp .or. m1.gt.mp1 .or. n.gt.np) then
         write(6,*) 'bad args in banbks'
         stop
      endif
c
c     Forward substitution unscrambling permuted rows
c
      l = m1
      do 12 k = 1, n
         i = indx(k)
         if (i .ne. k) then
            dum = b(k)
            b(k) = b(i)
            b(i) = dum
         endif
         if (l .lt. n) then 
            l = l + 1
         endif
         do 11 i = k + 1, l
            b(i) = b(i) - a1(k,i-k)*b(k)
 11      continue
 12   continue
c
c     Backsubstitution
c
      l = 1
      do 14 i = n, 1, -1
         dum = b(i)
         do 13 k = 2, l
            dum = dum - a(i,k)*b(k+i-1)
 13      continue
         b(i) = dum/a(i,1)
         if (l .lt. mm) then
            l = l + 1
         endif
 14   continue
c
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZINVERSIT(NDIM,NP,NP1,MM,REW,T,A,A1,B,BX,
     &                     IPVT,NINVI,SUCCESS)
c
c  this routine is based on numerical recipes and is
c  adapted for the Hermitian matrices
c

c     ndim: number of iterations in the lanczos. corresponds to niter
c           the true length of a T-eigenvector is niter
c     mmm:  this is just the blown-up dimension of T
c           corresponding to twice the bandwidth + 1.
c           remember: in the Numerical Recipes routine
c           the upper and lower part of the band-matrix are stored.
c    ninvi: maximum number of allowes inverse iterations. reason:
c           if an eigenvalue is converged the corresponding eigenvector
c           is immediately obtainable by 1 or two iterations. but it was
c           observed that in the case of multiple Lanczos vectors with a
c           corresp. eigenvector close to zero norm inverse iteration
c           does not lead to convergence due to numerical problems.
c           this case must be caught for stable running conditions.

      IMPLICIT INTEGER (A-Z)
c
      INTEGER NDIM, NP, NP1, MM, IPVT(NDIM), NINVI
      REAL*8 REW
      COMPLEX*16 T(NP1,NDIM), A(NDIM,MM), A1(NDIM,MM), B(NDIM), BX(NDIM)
      logical success
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "../relccsd/files.inc"
C
C--------------- Local variables ----------------------------------
C     
      complex*16 EW
      real*8 dum,y
      integer i, j,numit

      REAL*8 AUEV
      DATA AUEV/27.2113957D0/
C
C--------------- Executable code -----------------------------------
C     

c     ---------------------------------------------------
c     construct A from T according to the needs of NUMREC
c     REMEMBER that T is HERMITIAN and not
c     COMPLEX SYMMETRIC !!!!
c     Therefore we complex conjugate the missing triagonal
c     part in the following lines.
c     ---------------------------------------------------

      EW=DCMPLX(REW)
c     write(iw,*) ' +++++ zinversit parameter:'
c     write(iw,*) ' ++ NDIM, NP, NP1, MM, NINVI',
c    &            NDIM, NP, NP1, MM, NINVI
c     write(iw,*) ' +++++ EW:',EW


      DO J = 1, MM
         DO I = 1, NDIM
            A(I,J) = (0.0D0, 0.0D0)
         ENDDO
      ENDDO
c

      success = .true.
c
      do i = 1, ndim
         a(i,np1) = t(1,i)
      enddo
      do j = 1, np
         do i = 1, ndim - j
            a(i,np1+j) = DCONJG(t(j+1,i))
            a(i+j,np1-j) = t(j+1,i)
         enddo
      enddo
c     ------------------------
c     form A - lambda*I
c     ------------------------
      do i = 1, ndim
         a(i,np1) = a(i,np1) - ew
      enddo
c     -----------------------------------------------
c     LU Decomposition of a - ew
c     use second vector BX for convergence control of
c     the inverse iteration
c     -----------------------------------------------
      call zbandec(a,ndim,np,np,ndim,mm,a1,mm,ipvt,dum)
c
      call znormalize(ndim,b)

      numit = 0
 380  continue

      call zcopy(ndim,b,1,bx,1)
      call zbanbks(a,ndim,np,np,ndim,mm,a1,mm,ipvt,b)
      call znormalize(ndim,b)
      call znormdiff(ndim,b,bx,y)
      if(y.gt.1.0D-10) then
        numit = numit + 1
        if(numit.ge.ninvi) then
          write(iw,*) ' *** ATT ***'
          write(iw,*) ' The inverse iteration for eigenvalue',REW*AUEV,
     &                ' (eV)'
          write(iw,*) ' did not converge. The corresponding eigen-'
          write(iw,*) ' vector is not to be trusted!'
          write(iw,*) ' Final norm:',y
          write(iw,*) ' *** Remedy ***'
          write(iw,*) ' Decrease number of Lanczos iterations',
     &                ' or increase band width'
          
          success = .false.

          goto 382
        endif
        goto 380
      endif

 382  continue

      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZNORMALIZE (NDIM,B)
      implicit none
c     
      integer ndim
      complex*16 b(ndim)
c     
      integer i
      complex*16 y,norm
c     
      norm = (0.0d0, 0.0d0)
      do 10 i = 1, ndim
         norm = norm + b(i)*dconjg(b(i))
 10   continue
c     
      y = (1.0d0, 0.0d0) / sqrt(norm)
      do 20 i = 1, ndim
         b(i) = b(i) * y
 20   continue
c     
      return
      end
c
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZNORMDIFF (NDIM,V,W,DIFF)
c
c  this routine does not produce the ordinary vector difference 
c  norm but just compares the absolute values of the vector
c  components. this is due to possible phase changes in the inverse
c  iteration.
c
      implicit none
c     
      integer ndim
      complex*16 v(ndim),w(ndim)
      real*8 diff
c     
      integer i
      real*8 sum,d1
c     
      sum = 0.0d0
      do 10 i = 1, ndim
         d1 = abs(v(i)) - abs(w(i))
         sum = sum + d1*d1
 10   continue
      diff = dsqrt(sum)
c     
      return
      end
c     
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZBANDEC(A,N,M1,M2,NP,MP,A1,MP1,INDX,D)
c     ------------------------------------------------
c     LU decomposition of a band martix a
c     numerical recipes --> chapter 2, banddiagonal systems
c     a is destroyed and contatins the upper triangle
c     a1 contains the lower triangle
c     indx records the row permutations
c     -----------------------------------------------------
      implicit none
c
      integer m1, m2, mp, mp1, n, np, indx(n)
      real*8 tiny, d
      complex*16 a(np,mp), a1(np,mp1)
c
      parameter (tiny = 1.d-20)
c
      integer i, j, k, l, mm
      complex*16 dum
c
      mm = m1 + m2 + 1
      if (mm.gt.mp .or. m1.gt.mp1 .or. n.gt.np) then
         write(6,*) 'bad args in zbandec'
         stop
      endif
c
      l = m1
      do 13 i = 1, m1
         do 11 j = m1 + 2 - i, mm
            a(i,j-l) = a(i,j)
 11      continue
         l = l - 1
         do 12 j = mm - l, mm
            a(i,j) = (0.0d0, 0.0d0)
 12      continue
 13   continue
c
      d = 1.0d0
      l = m1
c
c     for each row ...
c
      do 18 k = 1, n
         dum = a(k,1)
         i = k
         if (l .lt. n) then
            l = l + 1
         endif
c
c     Find the pivot
c     
         do 14 j = k+1, l
            if (abs(a(j,1)) .gt. abs(dum)) then
               dum = a(j,1)
               i = j
            endif
 14      continue
c     
         indx(k) = i
         if (abs(dum) .eq. 0.0d0) then
            a(k,1) = dcmplx(tiny, 0.0d0)
         endif
c     
         if (i .ne. k) then
c     Interchange rows
            d = -d
            do 15 j = 1, mm
               dum = a(k,j)
               a(k,j) = a(i,j)
               a(i,j) = dum
 15         continue
         endif
c
c     Do the elimination
c     
         do 17 i = k + 1, l
            dum = a(i,1) / a(k,1)
            a1(k,i-k) = dum
            do 16 j = 2, mm
               a(i,j-1) = a(i,j) - dum*a(k,j)
 16         continue
            a(i,mm) = (0.0d0, 0.0d0)
 17      continue
c     
 18   continue
c     
      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZBANBKS(A,N,M1,M2,NP,MP,A1,MP1,INDX,B)
c     ------------------------------------------------------
c     solves Ax=b for x
c     A is a bandmatrix and has been decomposed in upper(a)
c     and lower(a1) triangle with bandec.
c     x is returned in b
c     ------------------------------------------------------
      implicit none
c
      integer m1, m2, mp, mp1, n, np, indx(n)
      complex*16 a(np,mp), a1(np,mp1), b(n)
c
      integer i, k, l, mm
      complex*16 dum
c
      mm = m1 + m2 + 1
      if (mm.gt.mp .or. m1.gt.mp1 .or. n.gt.np) then
         write(6,*) 'bad args in banbks'
         stop
      endif
c
c     Forward substitution unscrambling permuted rows
c
      l = m1
      do 12 k = 1, n
         i = indx(k)
         if (i .ne. k) then
            dum = b(k)
            b(k) = b(i)
            b(i) = dum
         endif
         if (l .lt. n) then
            l = l + 1
         endif
         do 11 i = k + 1, l
            b(i) = b(i) - a1(k,i-k)*b(k)
 11      continue
 12   continue
c
c     Backsubstitution
c
      l = 1
      do 14 i = n, 1, -1
         dum = b(i)
         do 13 k = 2, l
            dum = dum - a(i,k)*b(k+i-1)
 13      continue
         b(i) = dum/a(i,1)
         if (l .lt. mm) then
            l = l + 1
         endif
 14   continue
c
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RBLANA(IRREP,N,NP,NEIGENV,NHOLES,
     &                  D,ZLONG,Z,T,IND,MLT,LSP,
     &                  A,A1,B,BX,IPIV,IOUNIT,WTHR,ELIM)

      IMPLICIT INTEGER (A-Z)

      INTEGER IRREP,N,NP,NEIGENV,NHOLES
      REAL*8 D(N)
      REAL*8 ZLONG(N,N),Z(2*NP,N),T(NP+1,N)
      INTEGER IND(N),MLT(N),LSP(N)
      REAL*8 A(N,2*NP+1),A1(N,2*NP+1),B(N),BX(N)
      INTEGER IPIV(N),IOUNIT
      REAL*8 WTHR,ELIM
C
C---------------Description--------------------------------------------
C
C  Analyses the eigenvalues/vectors of the Band matrix T and
C  treats multiple/spurious eigenvalues. All this for the
C  current main program symmetry stored in IRREP
c  (H.-D. Meyer's and Thomas Sommerfeld's formulation)
c
c  the corresponding eigenvalues are written out in the file
c  EVAL.# for later analysis.
c
c  eigenvalues come in D(*), complete eigenvectors in ZLONG(N,N) !
c  the first and last NP components will be needed for the analysis 
c  only. in order to use the unmodified code we copy those components
c  to the short vector array z(2*np,n)
c
c  the number N corresponds to the number of L-iterations == NITER
c  the number NEIGENV is the number of desired eigenvectors
c
c  The array a, a1, b and IPIV  are needed for the numerical recipes routines. 
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
C---------------Common Blocks--------------------------------------
C
      INCLUDE "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C

      REAL*8 AA,AK,A1K,SUM,SUP
      REAL*8 AUEV,THRS,OVL
      REAL*8 X1,Y1

      REAL*8 XM,DDOT
      DATA AUEV/27.2113957D0/

      CHARACTER*7 EVALFN
      LOGICAL INVSUC
C
C---------------Executable code--------------------------------------
C


      CALL PST('Processing the spectrum+')
      WRITE(IW,*) 'Upper limit for printing eigenvalues:',ELIM
      WRITE(IW,*) 'Threshold for considering eigenvectors:',WTHR

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
      DO  48  J = 1, JND
         I = IND(J)
         d(J) = d(I)
         DO  49  K = 1, NP2
            Z(K,J) = Z(K,I)
 49      continue
 48   CONTINUE

      write (6,'(/,A)') '  Information on purging:'
      PRINT 905, JND, N, int(100.0d0*dble(jnd)/dble(n))
 905  FORMAT(4x,'There are',I5,' relevant residues out of',I5, 
     $     '   (',I3,'%)')
      write(6,'(4x,A,I4,6x,A,I4)') 'N_spurios = ', NSP,
     $     'N_multiple = ', NMU
      write(6,*) 'Threshold for vector printing:',WTHR
C
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
c
c  print eigenvalues and accuracy in a more convenient format
c  threshold applies only to the screen output. pole strength = square of
c  eigenvector component. Write all eigenvalues and pole strength
c  to file SPEC.X
c
      if(irrep.gt.9) then
        write(evalfn,'(a5,i2)') 'SPEC.',irrep
      else
        write(evalfn,'(a6,i1)') 'SPEC.0',irrep
      endif
      open(iounit,file=evalfn,form='FORMATTED',status='UNKNOWN')
      
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
c  to the SPEC.XX files.
c
        write(iounit,'(I5,2x,F10.4,4x,2F12.6)') j,d(j)*AUEV,x1,y1
c
c  print on screen only those eigenvalues that are below the given limit.
c  use cntev for contiguous numbering !
c
        if(d(j)*AUEV.le.ELIM) THEN
          cntev = cntev + 1
          write(iw,388)
 388      format(//2x,60('-')/'  no.    eigenvalue (eV)   mean-res.',
     &    3x,'error estimate'/2x,60('-')/)
          write(iw,389) cntev,'@', d(j)*AUEV, x1, y1, ind(j)
 389      format(I5,A1,2x,F14.8,4x,F8.4,4x,1PE11.2,'  (',I5,')')
          write(iw,'(100I8)') (ix,ix=1,nholes)
          write(iw,'(100F8.4)') 
     &         (z(ix,j)*z(ix,j),ix=1,nholes)
        endif

      ENDDO

      close(iounit)

      WRITE(IW,*)
      WRITE(IW,*)
         
C
C.....PRINT MEAN RESIDUES (X), EIGENVALUES (W), ERROR-EST.(Y) AND MULT.
C----------------------------------------------------------------

c     ----------------------------------------------------------------
c     ----  from here the spectrum is purged and the number of valid
c     ----  eigenvectors now sits in JND. ==> in case of an eigenvector
c     ----  run we have all information for the file.
c     ----  in the following another control output is generated and
c     ----  the inverse iterator called if required.
c     ----------------------------------------------------------------

      CONTINUE
      IF(NEIGENV.GT.0) THEN
        IF(NEIGENV.GT.N) NEIGENV = N
        WRITE(IW,*)
        WRITE(IW,*) '          --------------------------------'
        WRITE(IW,*) '          ---  doing inverse iteration for'
        WRITE(IW,*) '          --- ',NEIGENV,' requested eigenvectors'
        WRITE(IW,*) '          --------------------------------'
        WRITE(IW,*)
        WRITE(IW,*) '          Reading Lanzcos matrix'
c
c  read in Lanczos matrix (it was destroyed during diagonalization)
c  the Lanczos matrix is not used after a specific symmetry is treated
c  and can therefore be deleted afterwards.
c
        open(iounit,file='TMAT',form='UNFORMATTED',status='UNKNOWN')
        rewind(iounit)

        read(iounit) k1,k2
        write(iw,*) ' control parameter:',k1,k2
        if(k1.ne.(np+1).or.k2.ne.N)
     &    call quit('Differing parameter in Lanczos matrix.')
        do j=1,n
          read(iounit) (t(i,j),i=1,np+1)
        enddo
        close(iounit,status='DELETE')
c
c  open file for the Lanczos eigenvectors
c
        OPEN(IOUNIT,FILE='ZVEC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IOUNIT)
        WRITE(IOUNIT) NEIGENV
        WRITE(IOUNIT) N
        NINVI = 10
c
c  start inverse iteration
c
        WRITE(IW,400)
 400    FORMAT(//2x,60('-')/'  no.    mean-res.     eigenvalue (eV)',3x,
     &     'error estimate'/2x,60('-')/)
        MM  = 2*NP + 1
        DO J = 1, NEIGENV 
          x1 = 0
          do i = 1, nholes
             x1 = x1 + z(i,j)*z(i,j)
          enddo
c         x1 = dsqrt(x1)         
  
          y1 = 0.0d0
          do i = 1, np
            y1 = y1 + z(np+i,j)*z(np+i,j)
          enddo
          y1 = dsqrt(y1)
          write(6,420) j, x1, d(j)*AUEV, y1, ind(j)
 420        format(I5,2x,F8.4,4x,F14.8,4x,1PE11.2,'  (',I5,')')
c
c  initialize B
c
          DO I = 1, N
             B(I) = 0.0D0
          ENDDO
          DO I = 1, NP
             B(I)          = Z(I,J)
             B(N - NP + I) = Z(NP + I,J)
          ENDDO
          CALL RINVERSIT(N,NP,NP1,MM,D(J),T,A,A1,B,BX,IPIV,NINVI,INVSUC)
          IF(.NOT.INVSUC) THEN
            write(iw,*) ' *** Could not create eigenvector for eval',
     &                  J,'  ',D(J)*AUEV,' (eV)'
            do ix=1,n
              b(ix)=0.0d0
            enddo
          ENDIF
            
          WRITE(IOUNIT) D(J)
          WRITE(IOUNIT) (B(IX),IX=1,N)
        ENDDO

        CLOSE(IOUNIT)

        WRITE(IW,*)
        WRITE(IW,*) '    ------------------------------------------'
        WRITE(IW,*) '    ---  inverse iteration finished '
        WRITE(IW,*) '    ---  T-matrix eigenvectors written to file'
        WRITE(IW,*) '    ------------------------------------------'

      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZBLANA(IRREP,N,NP,NEIGENV,NHOLES,
     &                  D,ZLONG,Z,T,IND,MLT,LSP,
     &                  A,A1,B,BX,IPIV,IOUNIT,WTHR,ELIM)

      IMPLICIT INTEGER (A-Z)

      INTEGER IRREP,N,NP,NEIGENV,NHOLES
      REAL*8 D(N)
      COMPLEX*16 ZLONG(N,N),Z(2*NP,N),T(NP+1,N)
      INTEGER IND(N),MLT(N),LSP(N)
      COMPLEX*16 A(N,2*NP+1),A1(N,2*NP+1),B(N),BX(N)
      INTEGER IPIV(N),IOUNIT
      REAL*8 WTHR,ELIM
C
C---------------Description--------------------------------------------
C
C  Analyses the eigenvalues/vectors of the Band matrix T and
C  treats multiple/spurious eigenvalues. All this for the
C  current main program symmetry stored in IRREP
c  (H.-D. Meyer's and Thomas Sommerfeld's formulation)
c
c  the corresponding eigenvalues are written out in the file
c  EVAL.# for later analysis.
c
c  eigenvalues come in D(*), complete eigenvectors in ZLONG(N,N) !
c  the first and last NP components will be needed for the analysis 
c  only. in order to use the unmodified code we copy those components
c  to the short vector array z(2*np,n)
c
c  the number N corresponds to the number of L-iterations == NITER
c  the number NEIGENV is the number of desired eigenvectors
c
c  The array a, a1, b and IPIV  are needed for the numerical recipes routines. 
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
C---------------Common Blocks--------------------------------------
C
      INCLUDE "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C

      REAL*8 AA,AK,A1K,SUM,SUP
      REAL*8 AUEV,THRS,OVL
      REAL*8 X1,X2,Y1,Y2

      COMPLEX*16 XM,ZDOTC
      DATA AUEV/27.2113957D0/

      CHARACTER*7 EVALFN
      LOGICAL INVSUC
C
C---------------Executable code--------------------------------------
C


      CALL PST('Processing the spectrum+')
      WRITE(IW,*) 'Upper limit for printing eigenvalues:',ELIM
      WRITE(IW,*) 'Threshold for considering eigenvectors:',WTHR
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
      DO  48  J = 1, JND
         I = IND(J)
         d(J) = d(I)
         DO  49  K = 1, NP2
            Z(K,J) = Z(K,I)
 49      continue
 48   CONTINUE

      write (6,'(/,A)') '  Information on purging:'
      PRINT 905, JND, N, int(100.0d0*dble(jnd)/dble(n))
 905  FORMAT(4x,'There are',I5,' relevant residues out of',I5, 
     $     '   (',I3,'%)')
      write(6,'(4x,A,I4,6x,A,I4)') 'N_spurios = ', NSP,
     $     'N_multiple = ', NMU
      write(6,*) 'Threshold for vector printing:',WTHR
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
c  spectrum is purged and available (all eigenvalues with corresponding
c  intensities obtained as squares of the hole/hole components
c  (pole strengths). IMPORTANT: the band width of the Lanczos matrix
c  can be larger than the number of actual hole states due to stability
c  reasons. If we sum the intensities we just look at the hole states !
c
c
c  print eigenvalues and accuracy in a more convenient format
c  threshold applies only to the screen output. pole strength = square of
c  eigenvector component. Write all eigenvalues and pole strength
c  to file SPEC.X
c
      if(irrep.gt.9) then
        write(evalfn,'(a5,i2)') 'SPEC.',irrep
      else
        write(evalfn,'(a6,i1)') 'SPEC.0',irrep
      endif
      open(iounit,file=evalfn,form='FORMATTED',status='UNKNOWN')
      
      cntev=0
      DO J = 1, JND 
c
c  pole strength
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
c  to the SPEC.XX files.
c
        write(iounit,'(I5,2x,F10.4,4x,2F12.6)') j,d(j)*AUEV,x1,y1
c
c  screen print with threshold and hole components.
c  use cntev for contiguous numbering !
c
        if(d(j)*AUEV.le.ELIM) then
          cntev = cntev + 1
          write(iw,388)
 388      format(//2x,60('-')/'  no.    eigenvalue (eV)   mean-res.',
     &    3x,'error estimate'/2x,60('-')/)
          write(iw,389) cntev,'@', d(j)*AUEV, x1, y1, ind(j)
 389      format(I5,A1,2x,F14.8,4x,F8.4,4x,1PE11.2,'  (',I5,')')
          write(iw,'(100I8)') (ix,ix=1,nholes)
          write(iw,'(100F8.4)') 
     &         (dble(z(ix,j)*dconjg(z(ix,j))),ix=1,nholes)
        endif

      ENDDO

      close(iounit)

      WRITE(IW,*)
      WRITE(IW,*)
         
C
C.....PRINT MEAN RESIDUES (X), EIGENVALUES (W), ERROR-EST.(Y) AND MULT.
C----------------------------------------------------------------
C *********** ATT *********** we skip this clumsy output !!!
      GOTO 712
C


      ICNV = 0
      JCNV = 0
      WRITE(IW,900)
 900  FORMAT(//2X,'E I G E N V A L U E S'/2X,125('-')//2
     $     ('   NO.  MEAN-RES.*100  EIGEN-VL (eV)',3x
     $     ,'ERROR     M   S',9x) // 2X,125('-')/ )
      DO  330  J = 1, (JND+1)/2
         Y1 = 0.D0
         Y2 = 0.D0
         K  = JND + 1 - J
         DO  332  I = 1, NP
            Y1 = Y1 + DBLE(DCONJG(Z(NP+I,J))*Z(NP+I,J))
            Y2 = Y2 + DBLE(DCONJG(Z(NP+I,K))*Z(NP+I,K))
 332     continue
         Y1 = dsqrt(Y1)
         Y2 = dsqrt(Y2)
         IF(Y1.LT.1.D-9 .AND. J.NE.JND/2+1) ICNV = ICNV + 1
         IF(Y2.LT.1.D-9)  ICNV = ICNV + 1
         IF(Y1.LT.1.D-6 .AND. J.NE.JND/2+1) JCNV = JCNV + 1
         IF(Y2.LT.1.D-6)  JCNV = JCNV + 1
         X1 = 0.D0
         X2 = 0.D0
         DO  331  I = 1, NP
            X1 = X1 + DBLE(DCONJG(Z(I,J))*Z(I,J))
            X2 = X2 + DBLE(DCONJG(Z(I,K))*Z(I,K))
 331     continue
c
c
         X1 = DSQRT(X1)
         X2 = DSQRT(X2)

         WRITE(IW,901) J,X1,d(J)*AUEV,Y1,MLT(J),LSP(J), 
     &         K,X2,d(K)*AUEV,Y2,MLT(K),LSP(K)
 901     FORMAT(1X, 2( I4, 2PF12.6, 0PF16.9, 1PE11.2, 2I4, 11X ) ) 
 330  CONTINUE
      WRITE(IW,990)
 990  FORMAT(/2('   NO.  MEAN-RES.*100  EIGEN-VL (eV)',3x,
     &     'ERROR     M   S         '),//,2X,125('-')/ )
      WRITE(IW,902) ICNV, 100.D0*ICNV/DBLE(JND), 100.D0*ICNV/DBLE(n)
      WRITE(IW,903) JCNV, 100.D0*JCNV/DBLE(JND), 100.D0*JCNV/DBLE(n)
 902  FORMAT(/,'    NUMBER OF CONVERGED (1.D-9) EIGENVALUES =',
     &         I4, F11.2, ' %  (JND) ', F9.2,' %  (JEND) ' )
 903  FORMAT(/,'    NUMBER OF CONVERGED (1.D-6) EIGENVALUES =',
     &         I4, F11.2, ' %  (JND) ', F9.2,' %  (JEND) ' )     

c     ----------------------------------------------------------------
c     ----  from here the spectrum is purged and the number of valid
c     ----  eigenvectors now sits in JND. ==> in case of an eigenvector
c     ----  run we have all information for the file.
c     ----  in the following another control output is generated and
c     ----  the inverse iterator called if required.
c     ----------------------------------------------------------------

 712  CONTINUE
      IF(NEIGENV.GT.0) THEN
        IF(NEIGENV.GT.N) NEIGENV = N
        WRITE(IW,*)
        WRITE(IW,*) '          --------------------------------'
        WRITE(IW,*) '          ---  doing inverse iteration for'
        WRITE(IW,*) '          --- ',NEIGENV,' requested eigenvectors'
        WRITE(IW,*) '          --------------------------------'
        WRITE(IW,*)
        WRITE(IW,*) '          Reading Lanzcos matrix'
c
c  read in Lanczos matrix (it was destroyed during diagonalization)
c  the Lanczos matrix is not used after a specific symmetry is treated
c  and can therefore be deleted afterwards.
c
        open(iounit,file='TMAT',form='UNFORMATTED',status='UNKNOWN')
        rewind(iounit)

        read(iounit) k1,k2
        write(iw,*) ' control parameter:',k1,k2
        if(k1.ne.(np+1).or.k2.ne.N)
     &    call quit('Differing parameter in Lanczos matrix.')
        do j=1,n
          read(iounit) (t(i,j),i=1,np+1)
        enddo
        close(iounit,status='DELETE')
c
c  open file for the Lanczos eigenvectors
c
        OPEN(IOUNIT,FILE='ZVEC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IOUNIT)
        WRITE(IOUNIT) NEIGENV
        WRITE(IOUNIT) N
        NINVI = 10
c
c  start inverse iteration
c
        WRITE(IW,400)
 400    FORMAT(//2x,60('-')/'  no.    mean-res.     eigenvalue (eV)',3x,
     &     'error estimate'/2x,60('-')/)
        MM  = 2*NP + 1
        DO J = 1, NEIGENV 
          x1 = 0
          do i = 1, nholes
             x1 = x1 + DBLE(z(i,j)*dconjg(z(i,j)) )
          enddo
c         x1 = dsqrt(x1)         
  
          y1 = 0.0d0
          do i = 1, np
            y1 = y1 + DBLE(z(np+i,j)*dconjg(z(np+i,j)) )
          enddo
          y1 = dsqrt(y1)
          write(6,420) j, x1, d(j)*AUEV, y1, ind(j)
 420        format(I5,2x,F8.4,4x,F14.8,4x,1PE11.2,'  (',I5,')')
c
c  initialize B
c
          DO I = 1, N
             B(I) = (0.0D0,0.0D0)
          ENDDO
          DO I = 1, NP
             B(I)          = Z(I,J)
             B(N - NP + I) = Z(NP + I,J)
          ENDDO
          CALL ZINVERSIT(N,NP,NP1,MM,D(J),T,A,A1,B,BX,IPIV,NINVI,INVSUC)
          IF(.NOT.INVSUC) THEN
            write(iw,*) ' *** Could not create eigenvector for eval',
     &                  J,'  ',D(J)*AUEV,' (eV)'
            do ix=1,n
              b(ix)=(0.0d0,0.0d0)
            enddo
          ENDIF
            
          WRITE(IOUNIT) D(J)
          WRITE(IOUNIT) (B(IX),IX=1,N)
        ENDDO

        CLOSE(IOUNIT)

        WRITE(IW,*)
        WRITE(IW,*) '    ------------------------------------------'
        WRITE(IW,*) '    ---  inverse iteration finished '
        WRITE(IW,*) '    ---  T-matrix eigenvectors written to file'
        WRITE(IW,*) '    ------------------------------------------'

      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RBANDDIA(N,NP,LDAB,AB,EVL,Z,RWORK)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER N,NP,LDAB
      REAL*8 AB(LDAB,*),Z(N,*)
      REAL*8 RWORK(6*N),EVL(N)
      
C
C---------------Description--------------------------------------------
C
C  Perform the diagonalization of the real symmetric band matrix
C  by calling the appropriate LAPACK routine. This routine
C  is interfacing in order to avoid too many command line
C  arguments.
C  The Lanczos matrix T from the caller comes in the AB array.
C     ****  ATTENTION ****  due to the usual conventions in the
C    numerical recipes routines the original matrix T is DESTROYED !!!
C
C---------------Common Blocks--------------------------------------
C
      INCLUDE "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER INFO,LDZ
      CHARACTER*1 JOBZ,UPLO
      REAL*8 Y
      REAL*8 AUEV

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

      CALL PST('Diagonalizing real tridiagonal matrix+')


      WRITE(IW,*) 'Parameters for real band diagonalizer'
      WRITE(IW,*) 'N,NP,LDAB:',N,NP,LDAB

      CALL DSBEV(JOBZ, UPLO, N, NP, AB, LDAB, EVL, Z, LDZ,
     &           RWORK, INFO)
      WRITE(IW,*) 'Code from DSBEV:',INFO
      IF(INFO.NE.0) THEN
        CALL QUIT('Error in DSBEV')
      ENDIF

      WRITE(IW,*)
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*) 'First 20 RAW Lanczos EVAL/EVECS + intensity.'
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*)
      WRITE(IW,*) ' No.         (a.u.)      (eV)'


      AUEV=27.2113957D0

      DO I=1,MIN(20,N)
        Y = 0.0D0
        DO J=1,NP
          Y = Y + Z(J,I)*Z(J,I)
        ENDDO
        WRITE(IW,111) I,EVL(I),EVL(I)*AUEV,DABS(Y)
      ENDDO

 111  FORMAT(I5,2x,2F12.6,3x,'mean-res',2x,F10.4)


      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZBANDDIA(N,NP,LDAB,AB,EVL,Z,WORK,RWORK)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER N,NP,LDAB
      COMPLEX*16 AB(LDAB,*),Z(N,*),WORK(4*N)
      REAL*8 RWORK(6*N),EVL(N)
      
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
      INCLUDE "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER INFO,LDZ
      CHARACTER*1 JOBZ,UPLO
      COMPLEX*16 Y
      REAL*8 AUEV

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

      CALL PST('Diagonalizing hermitian tridiagonal matrix+')

      WRITE(IW,*) 'Parameters for complex band diagonalizer'
      WRITE(IW,*) 'N,NP,LDAB:',N,NP,LDAB

      CALL ZHBEV(JOBZ, UPLO, N, NP, AB, LDAB, EVL, Z, LDZ,
     &           WORK, RWORK, INFO)
      WRITE(IW,*) 'Code from ZHBEV:',INFO
      IF(INFO.NE.0) THEN
        CALL QUIT('Error in ZHBEV')
      ENDIF

      WRITE(IW,*)
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*) 'First 20 RAW Lanczos EVAL/EVECS + intensity.'
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*)
      WRITE(IW,*) ' No.         (a.u.)      (eV)'


      AUEV=27.2113957D0

      DO I=1,MIN(20,N)
        Y = (0.0D0,0.0D0)
        DO J=1,NP
          Y = Y + Z(J,I)*DCONJG(Z(J,I))
        ENDDO
        WRITE(IW,111) I,EVL(I),EVL(I)*AUEV,ABS(Y)
      ENDDO

c     WRITE(IW,'(/,/)')
c     WRITE(IW,*) 'Corresponding eigenvectors (first 30 entries)'
c     do iev=1,min(1,N)
c       WRITE(IW,*) iev,EVL(iev)
c       WRITE(IW,*) '----------------------'
c       do i=1,min(30,N)
c         write(iw,*) z(i,iev)
c       enddo
c       write(iw,*)
c       write(iw,*) 'norm :  ',zdotc(n,z(1,iev),1,z(1,iev),1)
c       write(iw,*)
c     enddo

 111  FORMAT(I5,2x,2F12.6,3x,'mean-res',2x,F10.4)

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
      INCLUDE "../relccsd/files.inc"
      INCLUDE "../relccsd/complex.inc"
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
C---------------Common Blocks--------------------------------------
C
      INCLUDE "../relccsd/files.inc"
      INCLUDE "../relccsd/complex.inc"
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
