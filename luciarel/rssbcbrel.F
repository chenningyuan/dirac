      SUBROUTINE RSSBCBREL(ONEINTS,INTSPACE,
     &     IASM,IATP,IOCTPA,IAOC,IAGRP,NIA,
     &     IBSM,IBTP,IOCTPB,IBOC,IBGRP,NIB,
     &     JASM,JATP,JOCTPA,JAOC,JAGRP,NJA,
     &     JBSM,JBTP,JOCTPB,JBOC,JBGRP,NJB,
     &     NGAS,SB,CB,IDOH2,
     &     ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &     NOBPTS,IOBPTS,NOBPTS2,IOBPTS2,
     &     MXPNGAS,MAXI,MAXK,
     &     SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &     NSMOB,NSMST,NSMSX,NSMDX,
     &     MXPOBS,IDC,
     &     CJRES,SIRES,I3,XI3S,I4,XI4S,
     &     MXSXBL,MXSXST,IPRNT,ISPINFREE,
     &     IPTSPC,JPTSPC,IJOP,SCLFAC,
     &     IRI,IRICO,IELCI,
     &     MXTSOB,SOMESCR,XIJILS)
*
*     Contributions to sigma block (iasm iatp, ibsm ibtp ) from
*     C block (jasm jatp , jbsm, jbtp)
*
*
*     Cleaned and shaved by Jeppe + Timo, March19, 1998
*     =====
*     Input
*     =====
*
*     IASM,IATP : Symmetry and type of alpha strings in sigma
*     IBSM,IBTP : Symmetry and type of beta  strings in sigma
*     JASM,JATP : Symmetry and type of alpha strings in C
*     JBSM,JBTP : Symmetry and type of beta  strings in C
*     NGAS      : Number of active spaces in calculation
*     IAOC,IBOC : Number of electrons in each AS for sigma supergroups
*     JAOC,JBOC : Number of electrons in each AS for C     supergroups
*     NAEL : Number of alpha electrons
*     NBEL : Number of  beta electrons
*     IJAGRP    : IA and JA belongs to this group of strings
*     IJBGRP    : IB and JB belongs to this group of strings
*     CB : Input c block
*     IDOH2 : = 0 => no two electron operator
*     IDOH2 : = 1 =>    two electron operator
*     ADASX : sym of a+, a => sym of a+a
*     ADSXA : sym of a+, a+a => sym of a
*     SXSTST : Sym of sx,!st> => sym of sx !st>
*     STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*     is nonvanishing by symmetry
*     DXSTST : Sym of dx,!st> => sym of dx !st>
*     STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*     is nonvanishing by symmetry
*     NTSOB  : Number of orbitals per type and symmetry
*     IBTSOB : base for orbitals of given type and symmetry
*     IBORB  : Orbitals of given type and symmetry
*     MAXI   : Largest Number of ' spectator strings 'treated simultaneously
*     MAXK   : Largest number of inner resolution strings treated at simult.
*
CTF
CTF   We need separate calls for real/imaginary parts according to
CTF   IRI: r/i part of sigma block
CTF   IRICO: r/i part of C block
CTF   (the appropriate integral list is then determined automatically)
CTF
CTF   IUB is given directly before routines are called for update of S
CTF   same for IPM1
CTF
*
*     ======
*     Output
*     ======
*     SB : fresh sigma block
*
*     =======
*     Scratch
*     =======
*     SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*     largest number of orbital pairs of given symmetries and
*     types.
*     I1, XI1S   : at least MXSTSO : Largest number of strings of given
*     type and symmetry
*     I1, XI1S   : at least MXSTSO : Largest number of strings of given
*     type and symmetry
*     C2 : Must hold largest STT block of sigma or C
*
*     XINT : Scratch space for integrals.
*
*     Jeppe Olsen , Winter of 1991
*
*     modified for relativistic double group calculations
*     Timo Fleig + Jeppe Olsen, March 1998
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
*
CTF DELTA_MK2 should be an integer!!
      integer DELTA_MK2
      INTEGER JKL_AC(3)
*
*     Input
      dimension CB(*)
*     Input (from WORK(KPICI1) in ... sblock)
      dimension ONEINTS(INTSPACE)
*
*     . Output
      DIMENSION SB(*)
*     . Scratch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*)
      DIMENSION C2(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION IBLOCK(8)
*
      call qenter('RSSBCBREL')
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
*
      XNORM = INPROD(CB,CB,NJA*NJB)
C     IF(XNORM.EQ.0.0D0) GOTO 1001
*
      if (NTEST.ge.100) call lmemchk
*
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' ====================='
         WRITE(6,*) ' RSSBCBREL :  C block '
         WRITE(6,*) ' ====================='
         CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
         WRITE(6,*) ' ============================='
         WRITE(6,*) ' RSSBCBREL : Initial  S block '
         WRITE(6,*) ' ============================='
         CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
         WRITE(6,*) ' Overall scalefactor ',SCLFAC
      END IF
*
      if (NTEST.ge.2000) then
        write(6,*)
        write(6,*) 'Checking integrals in rssbcbrel '
        do INT=1,INTSPACE,1
          write(6,'(I7,1F18.14)') INT,ONEINTS(INT)
        end do
      end if
*
      IF(NTEST.GE.500) THEN
         WRITE(6,*) ' IAOC and IBOC '
         CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
         CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
         WRITE(6,*) ' JAOC and JBOC  : '
         CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
         CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
         WRITE(6,*) ' IASM IATP JASM JATP ', IASM,IATP,JASM,JATP
         WRITE(6,*) ' IBSM IBTP JBSM JBTP ', IBSM,IBTP,JBSM,JBTP
         WRITE(6,*) ' NJA NJB NIA NIB ', NJA,NJB,NIA,NIB
      END IF
*     . No use of alpha-beta symmetry at the moment
      IUSEAB = 0
CTF
* Should the corresponding Hamiltonian matrix block be
* calculated exactly or approximately
*       NSEL2E = NNSEL2E
*  This is how it could be done somewhen (by passing NNSEL2E down)
*. Operator specified by input (for now, all excitations included)
      NSEL2E = 0
CTF
*     . Number of alpha- and beta- electrons in strings
      NIAEL = IELSUM(IAOC,NGAS)
      NIBEL = IELSUM(IBOC,NGAS)
      NJAEL = IELSUM(JAOC,NGAS)
      NJBEL = IELSUM(JBOC,NGAS)
*
      NIEL = NIAEL + NIBEL
      NJEL = NJAEL + NJBEL
*
      IF(NIEL.NE.NJEL) THEN
         WRITE(6,*)
     &        ' Total number of electrons in IaIb and JaJb differs'
         WRITE(6,*) ' NIAEL NIBEL NJAEL NJBEL', NIAEL,NIBEL,NJAEL,NJBEL
         Call Abend2
     &   ('RSSBCBREL: Tot. number of els. in IaIb and JaJb differs')
      END IF
*     . Difference in unbarred (A) and barred (B) occupation
*     ( number in sigma - number in C)
      DELTA_MK2 = 2*(NIAEL - NJAEL)
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' DELTA_MK2 = ', DELTA_MK2
      END IF
*
      IF (DELTA_MK2.EQ.0) THEN
         if (NTEST.ge.100) then
            write(6,*) '.........................'
            write(6,*) 'Entering DMK2=0 part ...'
            write(6,*) '.........................'
         end if
*
*     =================================
*     Sigma (delta Mk = 0) contribution
*     =================================
*
*     Sigma bb(IA,IB) = sum(i.gt.k,j.gt.l)<IB!Eb(ij)Eb(kl)!JB>
*     * ((ij!kl)-(il!kj)) C(IA,JB)
*     + sum(ij) <IB!Eb(ij)!JB> H(ij) C(IA,JB)
*     .One electron part beta
         JACTIVE = 0
         IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.NIBEL.GE.1) THEN
            if (NTEST.ge.100) then
               write(6,*) '........................'
               write(6,*) 'I have barred electrons'
               write(6,*) '........................'
            end if
            JACTIVE = 1
            CALL TRPMT3(SB,NIB,NIA,C2)
            CALL COPVEC(C2,SB,NIA*NIB)
            CALL TRPMT3(CB,NJB,NJA,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
*
            IF(NTEST.GE.500) THEN
               WRITE(6,*) ' SB before dmk0one'
               call wrtmat(sb,nia,nib,nia,nib)
            END IF
            IF(NTEST.GE.101)
     &           WRITE(6,*) ' I am going to call dmk0one'
*
*     Call new relativistic routines
C?          WRITE(6,*) 'Do memory check ...before dmk0one'
            if (NTEST.ge.100) call lmemchk
            IUB = 2
            IPM1 = 0
            call dmk0one(IBSM,IBTP,JBSM,JBTP,
     &           IBGRP,NIA,
     &           NGAS,IBOC,JBOC,
     &           SB,CB,
     &           ADSXA,STSTSX,
     &           NOBPTS2,IOBPTS2,
     &           MAXI,MAXK,
     &           SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &           NSMOB,NSMST,MOCAA,
     &           SCLFAC,
     &           ONEINTS,INTSPACE,
     &           IRI,IRICO,IPM1,IUB,NTEST)
*
            IF(NTEST.GE.500) THEN
               WRITE(6,*) ' SB after dmk0one'
               call wrtmat(sb,nia,nib,nia,nib)
            END IF
*
CTF  First call for two electron routines ...
CTF  First call for two electron routines ...
CTF  First call for two electron routines ...
CTF  First call for two electron routines ...
*
            IF(IDOH2.NE.0.AND.NIBEL.GE.2.and.IELCI.ge.2) THEN
               if (NTEST.ge.100) then
                  write(6,*) '...........................'
                  write(6,*) 'I even have > 1 barred el.'
                  write(6,*) '...........................'
               end if
*     . Two electron part beta-beta
               IF(NTEST.GE.101) then
                  WRITE(6,*) ' I am going to call dmk0two'
                  write(6,*) ' (Work is actually done by RSBB2A)'
               end if
*
C?             WRITE(6,*) 'Do memory check ...before RSBB2A'
C?             call memchk
*
*     . barred excitations ( no they are not barred)
               IBU = 1
               IF(IRI.NE.IRICO) THEN
                  IRI_INT = 2
               ELSE
                  IRI_INT = 1
               END IF
               IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag*imag => real gives minus
                  SCLFAC2 = - SCLFAC
               ELSE
                  SCLFAC2 = SCLFAC
               END IF
               CALL RSBB2A(IBSM,IBTP,JBSM,JBTP,IBGRP,NIA,NIB,
     &              NGAS,IBOC,JBOC,
     &              SB,CB,
     &              ADSXA,DXSTST,STSTDX,SXDXSX,
     &              NOBPTS2,IOBPTS2,ITSOB,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &              NSMOB,NSMST,NSMSX,NSMDX,
     &              CJRES,SIRES,MXSXST,MXSXBL,MOCAA,SCLFAC2,NTEST,
     &              NSEL2E,ISEL2E,IRI_INT,IBU,SOMESCR)
               IF(NTEST.GE.500) THEN
                  WRITE(6,*) ' SB after dmk0two'
                  call wrtmat(sb,nia,nib,nia,nib)
               END IF
*
            END IF
*           ^End if active 2 electron part
            CALL TRPMT3(SB,NIA,NIB,C2)
            CALL COPVEC(C2,SB,NIA*NIB)
            CALL TRPMT3(CB,NJA,NJB,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
         END IF
*        ^ End if IASM = JASM and IATP = JATP and atleast one electron
*
            IF(NIAEL.GE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
               if (NTEST.ge.100) then
                  write(6,*) '.......................'
                  write(6,*) 'I have unbarred elecs.'
                  write(6,*) '.......................'
               end if
*
*     alpha single excitation
*
               IF(NTEST.GE.101)
     &            WRITE(6,*) 'I am going to call dmk0one (last time)'
*
               IUB = 1
               IPM1 = 0
C?             WRITE(6,*) 'Do memory check ...before dmk0one'
C?             call memchk
*
               call dmk0one(IASM,IATP,JASM,JATP,
     &              IAGRP,NIB,
     &              NGAS,IAOC,JAOC,
     &              SB,CB,
     &              ADSXA,STSTSX,
     &              NOBPTS,IOBPTS,
     &              MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &              NSMOB,NSMST,MOCAA,
     &              SCLFAC,
     &              ONEINTS,INTSPACE,
     &              IRI,IRICO,IPM1,IUB,NTEST)
*
               IF(NTEST.GE.500) THEN
                  WRITE(6,*) ' SB transposed  after dmk0one'
                  call wrtmat(SB,nib,nia,nib,nia)
               END IF
            end if
*           ^ End NIAEL.GE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM
*
*     alpha double excitation
*
*
CTF  Second call for two electron routines ...
CTF  Second call for two electron routines ...
CTF  Second call for two electron routines ...
CTF  Second call for two electron routines ...
*
            if (IDOH2.ne.0.and.NIAEL.ge.2.and.IELCI.ge.2.AND.
     &          IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) then
               if (NTEST.ge.100) then
                WRITE(6,*) ' I am going to call RSBB2A (last time )'
                write(6,*) '..................................'
                write(6,*) 'meaning I have >1 unbarred elecs.'
                write(6,*) '..................................'
               end if
*
C?             WRITE(6,*) 'Do memory check ...before RSBB2A'
C?             call memchk
*
*     . unbarred excitations
               IBU = 2
               IF(IRI.NE.IRICO) THEN
                  IRI_INT = 2
               ELSE
                  IRI_INT = 1
               END IF
               IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
                  SCLFAC2 = - SCLFAC
               ELSE
                  SCLFAC2 = SCLFAC
               END IF
               CALL RSBB2A(IASM,IATP,JASM,JATP,IAGRP,NIB,NIA,
     &              NGAS,IAOC,JAOC,
     &              SB,CB,
     &              ADSXA,DXSTST,STSTDX,SXDXSX,
     &              NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &              NSMOB,NSMST,NSMSX,NSMDX,
     &              CJRES,SIRES,MXSXST,MXSXBL,MOCAA,SCLFAC2,NTEST,
     &              NSEL2E,ISEL2E,IRI_INT,IBU,SOMESCR)
            end if
*           ^ End if 2-electron aa is active
*
*
            IF(NTEST.GE.500) THEN
               WRITE(6,*) ' SB after dmk0two'
               call wrtmat(sb,nia,nib,nia,nib)
            END IF
*
CE?      END IF
*        ^End of IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM
*
*     Alpha-beta excitation
*
CTF  Third call for two electron routines ...
CTF  Third call for two electron routines ...
CTF  Third call for two electron routines ...
CTF  Third call for two electron routines ...
*
         if (NIAEL.ge.1.and.NIBEL.ge.1.and.IELCI.ge.2) then
            if (NTEST.ge.100) then
               write(6,*) '.................................'
               write(6,*) 'I have barred AND unbarred elecs.'
               write(6,*) '.................................'
            end if
*
C?          WRITE(6,*) 'Do memory check ...before RSBB2B'
C?          call memchk
*
            IF(IRI.NE.IRICO) THEN
               IRI_INT = 2
            ELSE
               IRI_INT = 1
            END IF
            IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
               SCLFAC2 = - SCLFAC
            ELSE
               SCLFAC2 = SCLFAC
            END IF
            CALL RSBB2B(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &           JASM,JATP,JBSM,JBTP,NJA,NJB,
     &           IAGRP,IBGRP,NGAS,
     &           IAOC,IBOC,JAOC,JBOC,
     &           SB,CB,
     &           ADSXA,STSTSX,
     &           NOBPTS,IOBPTS,NOBPTS2,IOBPTS2,
     &           ITSOB,MAXK,
     &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           XINT,NSMOB,NSMST,NSMSX,NSMDX,
     &           IUSEAB,ICJKAIB,CJRES,SIRES,C2,SCLFAC,NTEST,
     &           NSEL2E,ISEL2E,IRI_INT,SOMESCR)
         end if
      END IF
*     ^End of delta_mk2 = 0 part
*
      IF(NTEST.GE.500) THEN
         WRITE(6,*) ' SB after DMK2=0'
         call wrtmat(sb,nia,nib,nia,nib)
      END IF
*
      IF(DELTA_MK2 .EQ. 2.and.ISPINFREE.eq.0) THEN
         if (NTEST.ge.100) then
            write(6,*) '........................'
            write(6,*) 'Entering DMK2=2 part ...'
            write(6,*) '........................'
         end if
*
*     ===================================
*     Sigma (delta Mk = +1) contribution
*     ===================================
*
         IF(NTEST.GE.101)
     &        WRITE(6,*) ' I am going to call DMK1ONE'
*
         CALL TRPMT3(SB,NIB,NIA,C2)
         CALL COPVEC(C2,SB,NIA*NIB)
         CALL TRPMT3(CB,NJB,NJA,C2)
         CALL COPVEC(C2,CB,NJA*NJB)
*
         IPM1 = 1
         IUB = 0
         if (NTEST.ge.100) then
            write(6,*) '...................'
            write(6,*) 'dmk1one call with :'
            write(6,*) 'IPM1 =  ',IPM1
            write(6,*) 'IUB  =  ',IUB
            write(6,*) '...................'
         end if
*
C?       WRITE(6,*) 'Do memory check ...before dmk1one'
C?       call memchk
*
         call dmk1one(IASM,IBSM,IATP,IBTP,JASM,JBSM,
     &        JATP,JBTP,IAGRP,IBGRP,JAGRP,JBGRP,NIA,NIB,
     &        NJA,NJB,NGAS,IAOC,IBOC,JAOC,JBOC,
     &        SB,CB,
     &        ADSXA,STSTSX,NOBPTS,NOBPTS2,
     &        IOBPTS,IOBPTS2,MAXI,MAXK,
     &        SSCR,CSCR,I1,XI1S,I2,XI2S,
     &        XINT,NSMOB,NSMST,
     &        ONEINTS,INTSPACE,SCLFAC,
     &        IRI,IRICO,IPM1,IUB,NTEST)
*
C?       WRITE(6,*) 'Do memory check ...after dmk1one'
C?       call memchk
*
         if (NTEST.ge.500) then
            if (IPM1.eq.-1) write(6,*) 'MARKER-1'
            if (IPM1.eq.1) write(6,*) 'MARKER+1'
            write(6,*) 'Updated sigma after dmk1one'
            call wrtmat(SB,NIA,NIB,NIA,NIB)
         end if
*
*
CTF  Fourth call for two electron routines ...
CTF  Fourth call for two electron routines ...
CTF  Fourth call for two electron routines ...
CTF  Fourth call for two electron routines ...
*
         if (NIEL.ge.2.and.IELCI.ge.2) then
            if (NTEST.ge.100) then
               write(6,*) '.................'
               write(6,*) '>1 sigma electron'
               write(6,*) '.................'
            end if
*     . a+ ia a+ kb a lb a jb
            I_AC = 2
            JKL_AC(1) = 2
            JKL_AC(2) = 1
            JKL_AC(3) = 1
*
C?          WRITE(6,*) 'Do memory check ...before DMK1TWO'
C?          call memchk
*
            IF(IRI.NE.IRICO) THEN
               IRI_INT = 2
            ELSE
               IRI_INT = 1
            END IF
            IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
               SCLFAC2 = - SCLFAC
            ELSE
               SCLFAC2 = SCLFAC
            END IF
            IBU = 2
            JBU = 1
            KBU = 2
            LBU = 2
*
CTF
*  These routines are not tested yet!
*  In an initial relativistic CI run, all integrals
*  required in this part are = 0. This is because there
*  is no spin mixing of barred and unbarred spinors
*  (unbarred spinors are plain alpha, barred spinors
*   plain beta.) and therefore, spin integration makes
*  the integrals vanish.
CTF
*
c           CALL DMK1TWO(NIA,IASM,IATP,IAGRP,IAOC,
c    &           NIB,IBSM,IBTP,IBGRP,IBOC,
c    &           NJA,JASM,JATP,JAGRP,JAOC,
c    &           NJB,JBSM,JBTP,JBGRP,JBOC,
c    &           I_AC,JKL_AC,NGAS,SB,CB,
c    &           NOBPTS,NOBPTS2,MAXK,
c    &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
c    &           SCR,NSMOB,NSMST,IUSEAB,CJRES,SIRES,
c    &           SCLFAC2,NTEST,IUSE_PH,IPHGAS,IRI_INT,
c    &           IBU,JBU,KBU,LBU,
c    &           XIJILS)
*
C     DMK1TWO(NIA,IASM,IATP,IAGRP,IAOC,
C     &                    NIB,IBSM,IBTP,IBGRP,IBOC,
C     &                    NJA,JASM,JATP,JAGRP,JAOC,
C     &                    NJB,JBSM,JBTP,JBGRP,JBOC,
C     &          I_AC,JKL_AC,NGAS,SB,CB,MXPNGAS,
C     &          NALPTS,NBEPTS,MAXK,
C     &          SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
C     &          XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
C     &          CJRES,SIRES,SCLFAC,NTEST,
C     &          IUSE_PH,IPHGAS,
C     &          XIJILS)
*
         end if
*     . Restore order
         CALL TRPMT3(SB,NIA,NIB,C2)
         CALL COPVEC(C2,SB,NIA*NIB)
         CALL TRPMT3(CB,NJA,NJB,C2)
         CALL COPVEC(C2,CB,NJA*NJB)
*
         if (NTEST.ge.500) then
            write(6,*) 'Updated sigma after dmk1two'
            call wrtmat(SB,NIB,NIA,NIB,NIA)
         end if
*
CTF  Fifth call for two electron routines ...
CTF  Fifth call for two electron routines ...
CTF  Fifth call for two electron routines ...
CTF  Fifth call for two electron routines ...
*
         if (NIEL.ge.2.and.IELCI.ge.2) then
*     . a ib a+  ka a+ la a ja term
            I_AC = 1
            JKL_AC(1) = 2
            JKL_AC(2) = 2
            JKL_AC(3) = 1
            IBU = 1
            JBU = 2
            KBU = 1
            LBU = 1
*     . Proceed with blocks in the form C(Jb,Ja) and S(Ib,Ia) ( as on entry  )
C?          WRITE(6,*) 'Do memory check ...before DMK1TWO'
C?          call memchk
*
            IF(IRI.NE.IRICO) THEN
               IRI_INT = 2
            ELSE
               IRI_INT = 1
            END IF
            IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
               SCLFAC2 = - SCLFAC
            ELSE
               SCLFAC2 = SCLFAC
            END IF
c           CALL DMK1TWO(NIB,IBSM,IBTP,IBGRP,IBOC,
c    &           NIA,IASM,IATP,IAGRP,IAOC,
c    &           NJB,JBSM,JBTP,JBGRP,JBOC,
c    &           NJA,JASM,JATP,JAGRP,JAOC,
c    &           I_AC,JKL_AC,NGAS,SB,CB,
c    &           NOBPTS2,NOBPTS,MAXK,
c    &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
c    &           SCR,NSMOB,NSMST,IUSEAB,CJRES,SIRES,
c    &           SCLFAC2,NTEST,IUSE_PH,IPHGAS,IRI_INT,
c    &           IBU,JBU,KBU,LBU,XIJILS)
         end if
      END IF
*     ^End of DELTA_MK = 1 part
*
      IF(NTEST.GE.500) THEN
         WRITE(6,*) ' SB after DELTA_MK = +1'
         call wrtmat(sb,nia,nib,nia,nib)
      END IF
*
      IF(DELTA_MK2 .EQ. -2.and.ISPINFREE.eq.0) THEN
*
*     ===================================
*     Sigma (delta Mk = -1) contribution
*     ===================================
*
         IF(NTEST.GE.101)
     &        WRITE(6,*) ' I am going to call DMK1ONE'
*
         CALL TRPMT3(SB,NIB,NIA,C2)
         CALL COPVEC(C2,SB,NIA*NIB)
         CALL TRPMT3(CB,NJB,NJA,C2)
         CALL COPVEC(C2,CB,NJA*NJB)
*
         IPM1 =-1
         IUB = 0
         if (NTEST.ge.100) then
            write(6,*) '...................'
            write(6,*) 'dmk1one call with :'
            write(6,*) 'IPM1 =  ',IPM1
            write(6,*) 'IUB  =  ',IUB
            write(6,*) '...................'
         end if
*
C?       WRITE(6,*) 'Do memory check ...before dmk1one'
C?       call memchk
*
         call dmk1one(IASM,IBSM,IATP,IBTP,JASM,JBSM,
     &        JATP,JBTP,IAGRP,IBGRP,JAGRP,JBGRP,NIA,NIB,
     &        NJA,NJB,NGAS,IAOC,IBOC,JAOC,JBOC,
     &        SB,CB,
     &        ADSXA,STSTSX,NOBPTS,NOBPTS2,
     &        IOBPTS,IOBPTS2,MAXI,MAXK,
     &        SSCR,CSCR,I1,XI1S,I2,XI2S,
     &        XINT,NSMOB,NSMST,
     &        ONEINTS,INTSPACE,SCLFAC,
     &        IRI,IRICO,IPM1,IUB,NTEST)
*
         if (NTEST.ge.500) then
            if (IPM1.eq.-1) write(6,*) 'MARKER-1'
            if (IPM1.eq.1) write(6,*) 'MARKER+1'
            write(6,*) 'Updated sigma after dmk1one'
            call wrtmat(SB,NIA,NIB,NIA,NIB)
         end if
*
CTF  Sixth call for two electron routines ...
CTF  Sixth call for two electron routines ...
CTF  Sixth call for two electron routines ...
CTF  Sixth call for two electron routines ...
*
         if (NIEL.ge.2.and.IELCI.ge.2) then
*     . a  ia a+ jb a+ kb a lb
            IBU = 1
            JBU = 2
            KBU = 2
            LBU = 2
            I_AC = 1
            JKL_AC(1) = 2
            JKL_AC(2) = 2
            JKL_AC(3) = 1
C?          WRITE(6,*) 'Do memory check ...before DMK1TWO'
C?          call memchk
*
            IF(IRI.NE.IRICO) THEN
               IRI_INT = 2
            ELSE
               IRI_INT = 1
            END IF
            IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
               SCLFAC2 = - SCLFAC
            ELSE
               SCLFAC2 = SCLFAC
            END IF
c           CALL DMK1TWO(NIA,IASM,IATP,IAGRP,IAOC,
c    &           NIB,IBSM,IBTP,IBGRP,IBOC,
c    &           NJA,JASM,JATP,JAGRP,JAOC,
c    &           NJB,JBSM,JBTP,JBGRP,JBOC,
c    &           I_AC,JKL_AC,NGAS,SB,CB,
c    &           NOBPTS,NOBPTS2,MAXK,
c    &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
c    &           SCR,NSMOB,NSMST,IUSEAB,CJRES,SIRES,
c    &           SCLFAC2,NTEST,IUSE_PH,IPHGAS,IRI_INT,
c    &           IBU,JBU,KBU,LBU,XIJILS)
C     DMK1TWO(NIA,IASM,IATP,IAGRP,IAOC,
C     &                    NIB,IBSM,IBTP,IBGRP,IBOC,
C     &                    NJA,JASM,JATP,JAGRP,JAOC,
C     &                    NJB,JBSM,JBTP,JBGRP,JBOC,
C     &          I_AC,JKL_AC,NGAS,SB,CB,MXPNGAS,
C     &          NALPTS,NBEPTS,MAXK,
C     &          SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
C     &          XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
C     &          CJRES,SIRES,SCLFAC,NTEST,
C     &          IUSE_PH,IPHGAS,XIJILS)
*
         end if
*     . Restore order
         CALL TRPMT3(SB,NIA,NIB,C2)
         CALL COPVEC(C2,SB,NIA*NIB)
         CALL TRPMT3(CB,NJA,NJB,C2)
         CALL COPVEC(C2,CB,NJA*NJB)
*
CTF  Seventh call for two electron routines ...
CTF  Seventh call for two electron routines ...
CTF  Seventh call for two electron routines ...
CTF  Seventh call for two electron routines ...
*
         if (NIEL.ge.2.and.IELCI.ge.2) then
*     . a+ ib a+  ja a ka a la term
            IBU  = 1
            JBU = 2
            KBU =  2
            LBU = 2
            I_AC = 2
            JKL_AC(1) = 2
            JKL_AC(2) = 1
            JKL_AC(3) = 1
*     . Proceed with blocks in the form C(Jb,Ja) and S(Ib,Ia) ( as on entry  )
C?          WRITE(6,*) 'Do memory check ...before DMK1TWO'
C?          call memchk
*
            IF(IRI.NE.IRICO) THEN
               IRI_INT = 2
            ELSE
               IRI_INT = 1
            END IF
            IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
               SCLFAC2 = - SCLFAC
            ELSE
               SCLFAC2 = SCLFAC
            END IF
c           CALL DMK1TWO(NIB,IBSM,IBTP,IBGRP,IBOC,
c    &           NIA,IASM,IATP,IAGRP,IAOC,
c    &           NJB,JBSM,JBTP,JBGRP,JBOC,
c    &           NJA,JASM,JATP,JAGRP,JAOC,
c    &           I_AC,JKL_AC,NGAS,SB,CB,
c    &           NOBPTS2,NOBPTS,MAXK,
c    &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
c    &           SCR,NSMOB,NSMST,IUSEAB,CJRES,SIRES,
c    &           SCLFAC2,NTEST,IUSE_PH,IPHGAS,IRI_INT,
c    &           IBU,JBU,KBU,LBU,XIJILS)
         end if
      END IF
*     ^End of DELTA_MK = -1 part
*
      IF(NTEST.GE.500) THEN
         WRITE(6,*) ' SB after DELTA_MK = -1'
         call wrtmat(sb,nia,nib,nia,nib)
      END IF
*
CTF  Eighth call for two electron routines ...
CTF  Eighth call for two electron routines ...
CTF  Eighth call for two electron routines ...
CTF  Eighth call for two electron routines ...
*
      if (DELTA_MK2.eq.4.and.NIEL.ge.2.and.IELCI.ge.2.
     &                    and.ISPINFREE.eq.0) then
         if (NTEST.ge.100) then
            write(6,*) '........................'
            write(6,*) 'Entering DMK2=4 part ...'
            write(6,*) '........................'
         end if
*
*     ===================================
*     Sigma (delta Mk = +2) contribution
*     ===================================
*
*     +2 contribution:
C?       WRITE(6,*) 'Do memory check ...before dmk2two'
C?       call memchk
*
         IF(IRI.NE.IRICO) THEN
            IRI_INT = 2
         ELSE
            IRI_INT = 1
         END IF
         IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
            SCLFAC2 = - SCLFAC
         ELSE
            SCLFAC2 = SCLFAC
         END IF
CTF
         IIDMK2 = 2
*  See what's going on here in detail when routines are
*  tested ...
CTF
c        call dmk2two(IBSM,IBTP,NIB,IOCTPB,
c    &        IASM,IATP,NIA,IOCTPA,
c    &        JBSM,JBTP,NJB,JOCTPB,
c    &        JASM,JATP,NJA,JOCTPA,
c    &        NGAS,IBOC,IAOC,JBOC,JAOC,
c    &        SB,CB,
c    &        ADSXA,
c    &        NOBPTS,NOBPTS2,MAXI,MAXK,
c    &        SSCR,CSCR,I1,XI1S,I2,XI2S,
c    &        XINT,
c    &        NSMOB,NSMST,
c    &        SCLFAC2,NTEST,IRI_INT,IIDMK2,
c    &        SOMESCR)
      end if
*
*
*
CTF  Ninth call for two electron routines ...
CTF  Ninth call for two electron routines ...
CTF  Ninth call for two electron routines ...
CTF  Ninth call for two electron routines ...
*
      if (DELTA_MK2.eq.-4.and.NIEL.ge.2.and.IELCI.ge.2.
     &                    and.ISPINFREE.eq.0) then
         if (NTEST.ge.100) then
            write(6,*) '.........................'
            write(6,*) 'Entering DMK2=-4 part ...'
            write(6,*) '.........................'
         end if
*
*     ===================================
*     Sigma (delta Mk = -2) contribution
*     ===================================
*
C?       WRITE(6,*) 'Do memory check ...before dmk2two'
C?       call memchk
*
         IF(IRI.NE.IRICO) THEN
            IRI_INT = 2
         ELSE
            IRI_INT = 1
         END IF
         IF(IRI.EQ.1.AND.IRICO.EQ.2) THEN
*     imag**imag => real gives minus
            SCLFAC2 = - SCLFAC
         ELSE
            SCLFAC2 = SCLFAC
         END IF
         IIDMK2 = -2
c        call dmk2two(IASM,IATP,NIA,IOCTPA,
c    &        IBSM,IBTP,NIB,IOCTPB,
c    &        JASM,JATP,NJA,JOCTPA,
c    &        JBSM,JBTP,NJB,JOCTPB,
c    &        NGAS,IAOC,IBOC,JAOC,JBOC,
c    &        SB,CB,
c    &        ADSXA,
c    &        NOBPTS2,NOBPTS,MAXI,MAXK,
c    &        SSCR,CSCR,I1,XI1S,I2,XI2S,
c    &        XINT,
c    &        NSMOB,NSMST,
c    &        SCLFAC2,NTEST,IRI_INT,IIDMK2,
c    &        SOMESCR)
      end if
*
*
*
      IF(NTEST.GE.50) THEN
         WRITE(6,*) ' =========================='
         WRITE(6,*) ' RSSBCBREL : Final S block '
         WRITE(6,*) ' =========================='
         CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
      END IF
*
 1001 CONTINUE
      call qexit('RSSBCBREL')
      RETURN
      END
