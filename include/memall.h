      IERR   = 0
      NBYTES = (LWRK+2) * 8
      XBYTES = NBYTES
      XBYTES = XBYTES / (1024*1024)
#if defined (SYS_AIX) || defined (SYS_IRIX) || defined (SYS_SUN) || defined (SYS_HPUX)
#if defined (SYS_AIX)
      WADDR = MALLOC(%VAL(NBYTES))
#endif
#if defined (SYS_IRIX) || defined (SYS_SUN) || defined (SYS_HPUX)
      WADDR = MALLOC(NBYTES)
#endif
      IF (WADDR .EQ. 0) THEN
        WRITE(LUPRI,1000)
        WRITE(LUERR,1000)
        IERR = 0
      END IF
#endif
#if defined (SYS_CONVEX)
C
C     Convex dynamic memory allocation
C     ================================
C
      CALL NALLOC(NBYTES,IPTR,IERR)
      IF (IERR .LT. 0) THEN
        WRITE(LUERR,*) ' >>> Convex NALLOC failure; IERR = ', IERR
        WRITE(LUPRI,*) ' >>> Convex NALLOC failure; IERR = ', IERR
        IF (IERR .EQ. -3) THEN
           WRITE(LUPRI,1000)
           IF (LUPRI .NE. 0) WRITE(0,1000)
        END IF
      ENDIF
      IWOFFB = IPTR - LOC(WORK(1))
      IWOFF  = IWOFFB / 8
      IF (IWOFF*8 .NE. IWOFFB) THEN
        WRITE (LUERR,'(A,3(/A,I12))')
     &      'Convex alignment problem',
     &      '   pointer from NALLOC    :',IPTR,
     &      '   off-set to WORK, IWOFFB=',IWOFFB,
     &      '   mod(IWOFFB,8) =         ',IWOFFB - IWOFF*8
        IERR = -1
      ENDIF
#endif
#if defined (SYS_CRAY)
C     
C     Cray dynamic memory allocation
C     
C     HPALLOC IERR codes:
C     -1 : length allocated smaller than zero
C     -2 : no memory left for the system
C     -7 : memory corrupted.
C     HPALLOC IABORT .ne. 0: abort on error
        IABORT = 0
        CALL HPALLOC(WADDR,LWRK+2,IERR,IABORT)
        IF (IERR .LT. 0) THEN
          WRITE (LUERR,'(/A,I3)')
     &        ' >>> Cray HPALLOC failure; IERR   =', IERR
          IF (LUPRI .NE. 0) WRITE (LUPRI,'(/A,I3)')
     &        ' >>> Cray HPALLOC failure; IERR   =', IERR
          IF (IERR .EQ. -2) THEN
             WRITE(LUPRI,1000)
             WRITE(LUERR,1000)
          END IF
        ENDIF
#endif

