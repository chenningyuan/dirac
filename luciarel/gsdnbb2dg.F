      SUBROUTINE GSDNBB2DG(I12,
     &     RHO1UUR,RHO1UUI,
     &     RHO1UBR,RHO1UBI,
     &     RHO1BBR,RHO1BBI,
     &     RHO1BUR,RHO1BUI,
     &     IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &     NGAS,IAOC,IBOC,JAOC,JBOC,
     &     IAGRP,IBGRP,JAGRP,JBGRP,
     &     SB,CB,C2,
     &     ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &     MXPNGAS,NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &     MAXI,MAXK,
     &     SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &     NSMOB,NSMST,NSMSX,NSMDX,
     &     NIA,NIB,NJA,NJB,MXPOBS,IPRNT,NACOB,NACOB2,
     &     RHO1S,SCLFAC,IRIL,IRIR)
*
* Contributions to density matrix from sigma block (iasm iatp, ibsm ibtp ) and
* C block (jasm jatp , jbsm, jbtp)
*
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha strings in sigma
* IBSM,IBTP : Symmetry and type of beta  strings in sigma
* JASM,JATP : Symmetry and type of alpha strings in C
* JBSM,JBTP : Symmetry and type of beta  strings in C
* NGAS : Number of As'es
* IAOC : Occpation of each AS for alpha strings in L
* IBOC : Occpation of each AS for beta  strings in L
* JAOC : Occpation of each AS for alpha strings in R
* JBOC : Occpation of each AS for beta  strings in R
* IAGRP    : IA belongs to this group of strings
* IBGRP    : IB belongs to this group of strings
* IBGRP    : IB belongs to this group of strings
* IBGRP    : IB belongs to this group of strings
* CB : Input c block
* SB : Input S block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* MXPNGAS : Largest number of As'es allowed by program
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO1BBR/I : Real/Imaginary part of Barred-Barred     1-particle density
* RHO1UUR/I : ...................... Unbarred-Unbarred ..................
* RHO1UBR/I : ...................... Unbarred-Barred ....................
* RHO1BUR/I : ...................... Barred-Unbarred ....................
*
* =======
* Scratch
* =======
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* C2 : Must hold largest STT block of sigma or C
* RHO1S : Space for intermediate multiplication
*
*------------------------------------------------------------
* Jeppe Olsen , Winter of 1991
*------------------------------------------------------------
*
*------------------------------------------------------------
* modified for relativistic double group CI
*  Timo Fleig, April/May 1998
*              December 1998 / January 1999
*------------------------------------------------------------
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
      integer DELTA_MK2
*. Input
      DIMENSION CB(*),SB(*)
      dimension IOBPTS(MXPNGAS,*),IOBPTS2(MXPNGAS,*)
      dimension NOBPTS(MXPNGAS,*),NOBPTS2(MXPNGAS,*)
*     . Output
      dimension RHO1UUR(*),RHO1UUI(*),RHO1UBR(*),RHO1UBI(*)
      dimension RHO1BBR(*),RHO1BBI(*),RHO1BUR(*),RHO1BUI(*)
*
*. Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),I2(*),I3(*),I4(*)
      dimension XI1S(*),XI2S(*),XI3S(*),XI4S(*)
      DIMENSION C2(*),RHO1S(*)
*
      CALL QENTER('GSDNBB2DG')
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' GSDNBB2DG :  R block '
        WRITE(6,*) ' ====================='
        CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' GSDNBB2DG :  L block '
        WRITE(6,*) ' ====================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
*
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in L '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in L '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in R '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in R '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
*
        WRITE(6,'(A16,3I5)') ' MAXI,MAXK,NSMOB',MAXI,MAXK,NSMOB
*
        WRITE(6,*) 'SCLFAC =',SCLFAC
      END IF
      IACTIVE = 0
*
* Determine which part of this density block is being calculated
* real part or imaginary part
      if (IRIL.eq.IRIR) then
         if (NTEST.ge.200)
     &        write(6,*) 'Calculating real part of density matrix ...'
         IRIDEN = 1
         SCALEF = 1
      else if (IRIL.ne.IRIR) then
         if (NTEST.ge.200)
     &        write(6,*)
     $        'Calculating imaginary part of density matrix ...'
         IRIDEN = -1
         SCALEF = 1
      end if
*
*. Determine number of electrons in left and right strings and their type (UB)
      NLAEL = IELSUM(IAOC,NGAS)
      NLBEL = IELSUM(IBOC,NGAS)
      NRAEL = IELSUM(JAOC,NGAS)
      NRBEL = IELSUM(JBOC,NGAS)
*
      NLEL = NLAEL + NLBEL
      NREL = NRAEL + NRBEL
*
      if (NTEST.ge.200) then
         write(6,*) 'NLAEL =  ',NLAEL
         write(6,*) 'NLBEL =  ',NLBEL
         write(6,*) 'NRAEL =  ',NRAEL
         write(6,*) 'NRBEL =  ',NRBEL
         write(6,*) 'total left   : ',NLEL
         write(6,*) 'total right  : ',NREL
      end if
*
      if (NLEL.ne.NREL) then
         write(6,*) 'Total number of electrons ',
     &              'in left and right strings differs.'
         write(6,'(A25,4I4)') 'NLAEL,NLBEL,NRAEL,NRBEL  ',
     &                         NLAEL,NLBEL,NRAEL,NRBEL
         Call Abend2( 'gsdnbb2dg: problem with #els. in strings' )
      end if
*
      DELTA_MK2 = 2*(NLAEL-NRAEL)
      if (NTEST.ge.200)
     &     write(6,*) 'DELTA_MK2 =  ',DELTA_MK2
      IPM1 = DELTA_MK2 / 2
      if (NTEST.ge.200)
     &     write(6,*) 'IPM1 =  ',IPM1
*
      if (DELTA_MK2.eq.0) then
         if (NTEST.ge.100) then
            write(6,*) '............................'
            write(6,*) ' Entering DELTA_MK2 = 0 part'
            write(6,*) '............................'
         end if
         IF(NLBEL.GE.1.AND.IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* ===================================
*  barred-barred contribution to RHO1
* ===================================
*
            if (IRIDEN.eq.1) then
               if (NTEST.ge.1000) then
                  write(6,*) 'gsbbd1dg will be called'
                  write(6,*) 'barred-barred contribution'
                  write(6,*) 'IRIDEN is  ',IRIDEN
                  write(6,*) 'REDUNDANT! Is being skipped ...'
               end if
*
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1BBR'
                 write(6,*) '////////////////////'
               end if
               CALL GSBBD1DG(RHO1BBR,NACOB2,
     &              IBSM,IBTP,JBSM,JBTP,IBGRP,NIA,
     &              NGAS,IBOC,JBOC,
     &              SB,CB,
     &              NOBPTS2,IOBPTS2,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,
     &              NSMOB,RHO1S,
     &              SCALEF,SCLFAC,IPRNT)
            else if (IRIDEN.eq.-1) then
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1BBI'
                 write(6,*) '////////////////////'
               end if
               CALL GSBBD1DG(RHO1BBI,NACOB2,
     &              IBSM,IBTP,JBSM,JBTP,IBGRP,NIA,
     &              NGAS,IBOC,JBOC,
     &              SB,CB,
     &              NOBPTS2,IOBPTS2,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,
     &              NSMOB,RHO1S,
     &              SCALEF,SCLFAC,IPRNT)
            end if
         end if
*
         IF(NLAEL.GE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* =======================================
*  unbarred-unbarred contribution to RHO1
* =======================================
*
            CALL TRPMT3(CB,NJA,NJB,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
            CALL TRPMT3(SB,NIA,NIB,C2)
            CALL COPVEC(C2,SB,NIA*NIB)
            if (IRIDEN.eq.1) then
               if (NTEST.ge.1000) then
                  write(6,*) 'gsbbd1dg will be called'
                  write(6,*) 'unbarred-unbarred contribution'
                  write(6,*) 'IRIDEN is  ',IRIDEN
               end if
*
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1UUR'
                 write(6,*) '////////////////////'
               end if
               CALL GSBBD1DG(RHO1UUR,NACOB,
     &              IASM,IATP,JASM,JATP,IAGRP,NIB,
     &              NGAS,IAOC,JAOC,
     &              SB,CB,
     &              NOBPTS,IOBPTS,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,
     &              NSMOB,RHO1S,
     &              SCALEF,SCLFAC,IPRNT)
            else if (IRIDEN.eq.-1) then
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1UUI'
                 write(6,*) '////////////////////'
               end if
               CALL GSBBD1DG(RHO1UUI,NACOB,
     &              IASM,IATP,JASM,JATP,IAGRP,NIB,
     &              NGAS,IAOC,JAOC,
     &              SB,CB,
     &              NOBPTS,IOBPTS,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,
     &              NSMOB,RHO1S,
     &              SCALEF,SCLFAC,IPRNT)
            end if
         end if
*
C?      WRITE(6,*) ' GSBBD1DG was called '
C?      WRITE(6,*) ' Memory check '
C?      CALL LMEMCHK
*
* ================================
* beta-beta contribution to RHO2
* ================================
*
         IF(I12.EQ.2.AND.NBEL.GE.2) THEN
C?        WRITE(6,*) ' GSBBD2A will be called (beta)'
CTF       CALL GSBBD2A(RHO2,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
CTF  &         NGAS,IBOC,JBOC,SB,CB,
CTF  &         ADSXA,SXSTST,STSTSX,SXDXSX,
CTF  &         NOBPTS,IOBPTS,MAXI,MAXK,
CTF  &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
CTF  &         NSMOB,NSMST,NSMSX,SCLFAC)
C?        WRITE(6,*) ' GSBBD2A was called '
*
C              GSBBD2A(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
C    &         NGAS,ISEL,ICEL,SB,CB,
C    &         ADSXA,SXSTST,STSTSX,SXDXSX,
C    &         NOBPTS,IOBPTS,MAXI,MAXK,
C    &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
C    &         NSMOB,NSMST,NSMSX)
         END IF
         IF(I12.EQ.2.AND.NAEL.GE.2) THEN
*
*
* ===================================
*  alpha-alpha contribution to RHO2
* ===================================
*
C?        WRITE(6,*) ' GSBBD2A will be called (alpha)'
CTF       CALL GSBBD2A(RHO2,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
CTF  &         NGAS,IAOC,JAOC,SB,CB,
CTF  &         ADSXA,SXSTST,STSTSX,SXDXSX,
CTF  &         NOBPTS,IOBPTS,MAXI,MAXK,
CTF  &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
CTF  &         NSMOB,NSMST,NSMSX,SCLFAC)
C?      WRITE(6,*) ' GSBBD2A was called '
            CALL TRPMT3(CB,NJB,NJA,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
            CALL TRPMAT(SB,NIB,NIA,C2)
            CALL COPVEC(C2,SB,NIB*NIA)
         END IF
*
      else if (DELTA_MK2.eq.-2) then
         if (NTEST.ge.100) then
            write(6,*) '.............................'
            write(6,*) ' Entering DELTA_MK2 = -2 part'
            write(6,*) '.............................'
         end if

         if (NLAEL.ge.1.and.NLBEL.ge.1) then
*
* =====================================
*  barred-unbarred contribution to RHO1
* =====================================
*
            if (IRIDEN.eq.1) then
              if (NTEST.ge.1000) then
                 write(6,*) 'gsubd1dg will be called'
                 write(6,*) 'barred-unbarred contribution'
                 write(6,*) 'IRIDEN is  ',IRIDEN
                 write(6,*) 'IPM1 is  ',IPM1
                 write(6,*) 'REDUNDANT! Is being skipped ...'
              end if
*
              if (NTEST.ge.100) then
                write(6,*) '////////////////////'
                write(6,*) 'Writing to RHO1BUR'
                write(6,*) '////////////////////'
              end if
              call gsubd1dg(RHO1BUR,NACOB2,
     &             IASM,IBSM,IATP,IBTP,
     &             JASM,JBSM,JATP,JBTP,IAGRP,IBGRP,
     &             JAGRP,JBGRP,NIA,NIB,NJA,NJB,NGAS,
     &             IAOC,IBOC,JAOC,JBOC,
     &             SB,CB,
     &             NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &             MAXI,MAXK,SSCR,CSCR,
     &             I1,XI1S,I2,XI2S,NSMOB,
     &             RHO1S,IPM1,SCALEF,SCLFAC,IPRNT)
           else if (IRIDEN.eq.-1) then
              if (NTEST.ge.100) then
                write(6,*) '////////////////////'
                write(6,*) 'Writing to RHO1BUI'
                write(6,*) '////////////////////'
              end if
              call gsubd1dg(RHO1BUI,NACOB2,
     &             IASM,IBSM,IATP,IBTP,
     &             JASM,JBSM,JATP,JBTP,IAGRP,IBGRP,
     &             JAGRP,JBGRP,NIA,NIB,NJA,NJB,NGAS,
     &             IAOC,IBOC,JAOC,JBOC,
     &             SB,CB,
     &             NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &             MAXI,MAXK,SSCR,CSCR,
     &             I1,XI1S,I2,XI2S,NSMOB,
     &             RHO1S,IPM1,SCALEF,SCLFAC,IPRNT)
           end if
         end if
*
      else if (DELTA_MK2.eq.2) then
         if (NTEST.ge.100) then
            write(6,*) '............................'
            write(6,*) ' Entering DELTA_MK2 = 2 part'
            write(6,*) '............................'
         end if
*
* =====================================
*  unbarred-barred contribution to RHO1
* =====================================
*
         if (NTEST.ge.1000) then
            write(6,*) 'gsubd1dg will be called'
            write(6,*) 'unbarred-barred contribution'
            write(6,*) 'IRIDEN is  ',IRIDEN
            write(6,*) 'IPM1 is  ',IPM1
         end if
*
         if (IRIDEN.eq.1) then
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1UBR'
                 write(6,*) '////////////////////'
               end if
            call gsubd1dg(RHO1UBR,NACOB,
     &           IASM,IBSM,IATP,IBTP,
     &           JASM,JBSM,JATP,JBTP,IAGRP,IBGRP,
     &           JAGRP,JBGRP,NIA,NIB,NJA,NJB,NGAS,
     &           IAOC,IBOC,JAOC,JBOC,
     &           SB,CB,
     &           NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &           MAXI,MAXK,SSCR,CSCR,
     &           I1,XI1S,I2,XI2S,NSMOB,
     &           RHO1S,IPM1,SCALEF,SCLFAC,IPRNT)
         else if (IRIDEN.eq.-1) then
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1UBI'
                 write(6,*) '////////////////////'
               end if
            call gsubd1dg(RHO1UBI,NACOB,
     &           IASM,IBSM,IATP,IBTP,
     &           JASM,JBSM,JATP,JBTP,IAGRP,IBGRP,
     &           JAGRP,JBGRP,NIA,NIB,NJA,NJB,NGAS,
     &           IAOC,IBOC,JAOC,JBOC,
     &           SB,CB,
     &           NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &           MAXI,MAXK,SSCR,CSCR,
     &           I1,XI1S,I2,XI2S,NSMOB,
     &           RHO1S,IPM1,SCALEF,SCLFAC,IPRNT)
         end if
*
      else
         write(6,*) 'DELTA_MK2 is  ',DELTA_MK2
         write(6,*) 'These cases are not debugged yet'
         write(6,*) 'I am skipping the relevant density routines'
      end if
*
      idelabs = abs(DELTA_MK2)
      if (IDELABS.ge.4) then
         write(6,*) 'CAUTION !!!'
         write(6,*) 'CAUTION !!!'
         write(6,*) 'Skipping RHO2 routines ...'
         write(6,*) 'CAUTION !!!'
         write(6,*) 'CAUTION !!!'
         goto 666
      end if
*
      if (DELTA_MK2.eq.-4) then
         if (NTEST.ge.100) then
            write(6,*) '.............................'
            write(6,*) ' Entering DELTA_MK2 = -4 part'
            write(6,*) '.............................'
         end if
*
* ===================================
*  alpha-beta contribution to RHO2
* ===================================
*
         IF(I12.EQ.2.AND.NAEL.GE.1.AND.NBEL.GE.1) THEN
*. Routine uses transposed blocks
            CALL TRPMT3(CB,NJA,NJB,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
            CALL TRPMT3(SB,NIA,NIB,C2)
            CALL COPVEC(C2,SB,NIA*NIB)
C?      WRITE(6,*) ' GSBBD2B will be called '
CTF     CALL GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
CTF  &                    JASM,JATP,JBSM,JBTP,NJA,NJB,
CTF  &                    IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
CTF  &                    SB,CB,ADSXA,STSTSX,
CTF  &                    NOBPTS,IOBPTS,MAXK,
CTF  &                    I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
CTF  &                    NSMOB,NSMST,NSMSX,NSMDX,IUSEAB,
CTF  &                    SSCR,CSCR,NACOB,NTEST,SCLFAC)
C?      WRITE(6,*) ' GSBBD2B was called '
CTF  &
C     GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
C    &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
C    &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
C    &                  SB,CB,ADSXA,STSTSX,
C    &                  NOBPTS,IOBPTS,MAXK,
C    &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
C    &                  NSMOB,NSMST,NSMSX,NSMDX,IUSEAB,
C    &                  CJRES,SIRES,NORB,NTEST)
            CALL TRPMT3(CB,NJB,NJA,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
            CALL TRPMAT(SB,NIB,NIA,C2)
            CALL COPVEC(C2,SB,NIB*NIA)
         end if
      else if (DELTA_MK2.eq.4) then
         if (NTEST.ge.100) then
            write(6,*) '............................'
            write(6,*) ' Entering DELTA_MK2 = 4 part'
            write(6,*) '............................'
         end if
*
         write(6,*) 'DELTA_MK2 is  ',DELTA_MK2
         write(6,*) 'This part is still missing ...'
      END IF
*
 666  continue
      CALL QEXIT('GSDNBB2DG')
      RETURN
      END
