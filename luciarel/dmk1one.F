      subroutine dmk1one(ISCSMA,ISCSMB,ISCTPA,ISCTPB,
     &                   ICCSMA,ICCSMB,ICCTPA,ICCTPB,
     &                   IAGRP,IBGRP,JAGRP,JBGRP,
     &                   NIROWA,NIROWB,
     &                   NJROWA,NJROWB,
     &                   NGAS,
     &                   ISELA,ISELB,ICELA,ICELB,
     &                   SB,CB,
     &                   ADSXA,STSTSX,
     &                   NOBPTS_U,NOBPTS_B,
     &                   IOBPTS_U,IOBPTS_B,MAXL,MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,
     &                   H,
     &                   NSMOB,NSMST,
     &                   ONEINTS,INTSPACE,SCLFAC,
     &                   IRI,IRICO,IPM1,IUB,IPRNT)
*
*  Info for alpha strings has been added. We need both info for
*  alpha and beta at the same time!!
*  This has to be taken care of in the calling routine,
*  which is: rssbcbrel
*  Be careful: The CA strings are assumed to be N-1 strings already.
*  But the CB strings are assumed to be N+1 strings.
*
* =====
* Input
* =====
*
* ISCSMA,ISCTPA : Symmetry and type of sigma columns for alpha type
* ICCSMA,ICCTPA : Symmetry and type of C     columns for alpha type
* IAGRP : String group of columns, unbarred strings
* IBGRP : String group of columns, barred strings
* NGAS : Number of active sets
* ISELA : Occupation in each active set for sigma block for alpha type
* ICELA : Occupation in each active set for C     block for alpha type
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST : Number of symmetries of orbitals and strings
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* IRI : real or imaginary part of sigma vector?
*       real :1
*       imaginary :2
* IRICO : real or imaginary part of coefficients as input?
*       real :1
*       imaginary :2
* IPM1 : either +1 (delta Mk=+1) or -1 (delta Mk=-1)
* IUB : should be set to 0 here (only relevant in dmk0one)
*
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* H : Space for one electron integrals
* (The program will itself decide what type (real or imaginary) of
*   integrals will be needed and these fetched)
*
*  Jeppe Olsen, Winter of 1991
*
*  changed for relativistic double group calculation
*  Timo Fleig, Summer of 1997
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
*  symmetry info for double group CI
#include "symm.inc"
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS_U(MXPNGAS,*),NOBPTS_B(MXPNGAS,*)
      INTEGER IOBPTS_U(MXPNGAS,*),IOBPTS_B(MXPNGAS,*)
*.Specific Input
      INTEGER ISELA(NGAS),ICELA(NGAS),ISELB(NGAS),ICELB(NGAS)
      character*72 ENVIRO
*
* Input
      dimension ONEINTS(INTSPACE)
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*),H(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
      dimension SSCRH1(MXNSPIPDMGRP*MXNSPIPDMGRP)
*.Local arrays ( assume MPNGAS = 16 ) !!!
C      DIMENSION ITP(16),JTP(16)
* Type of single excitations that connects the two column strings
*
      call qenter('DMK1ONE')
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.500)THEN
        write(6,*)
        write(6,*) ' ========================= '
        write(6,*) ' Information from DMK1ONE '
        write(6,*) ' ========================= '
        write(6,*)
        write(6,*) ' IPM1 = ',IPM1
        write(6,*) ' ISELA : '
        CALL IWRTMA(ISELA,1,NGAS,1,NGAS)
        write(6,*) ' ICELA : '
        CALL IWRTMA(ICELA,1,NGAS,1,NGAS)
        write(6,*) ' ISELB : '
        call iwrtma(ISELB,1,NGAS,1,NGAS)
        write(6,*) ' ICELB : '
        call iwrtma(ICELB,1,NGAS,1,NGAS)
*
        write(6,*) ' NIROWA, NIROWB ', NIROWA, NIROWB
        write(6,*) ' NJROWA, NJROWB ', NJROWA, NJROWB
*
        if (NTEST.ge.2000) then
          WRITE(6,*)  ' Input C block '
          CALL WRTMAT(CB, NJROWA, NJROWB, NJROWA, NJROWB)
          WRITE(6,*)  ' Input S blocks '
          CALL WRTMAT(SB,NIROWA, NIROWB,NIROWA, NIROWB)
        end if
*
        WRITE(6,*) ' ISCSMA,ISCSMB ', ISCSMA,ISCSMB
        WRITE(6,*) ' ICCSMA,ICCSMB ', ICCSMA,ICCSMB
*
      END IF
*
      IFRST = 1
      JFRST = 1
*
*  1) Determine C(Ja,Ib,j) = sum(Jb) <Jb| j_+ |Ib> . C(Ja,Jb)
*     (if IPM1 = -1 then
*               C(Jb,Ia,j) = sum(Ja) <Ja| j+ |Ia> . C(Jb,Ja) )
*    a) find single excitations connecting spectator |Ib> = ISELB
*                                         and column <Jb| = ICELB
*       same for <Ia| i+ |Ja>
*               (<Ib| i_+ |Jb>)
*
* ===========================
* The algorithm for IPM1 = 1
* ===========================
* Sigma(Ia,Ib) = <Ia|i+|Ja><Ib|j_|Jb>C(Ja,Jb) H(i,j)
*
* Part 1 : C(Ja,Kb,j) = C(Ja,Ib,j) <Ib|j_|Jb>C(Ja,Jb)
*          (Proceed as normal 1-electron case)
*
* Part 2 : S(Ja,Ib,i) = sum(j) C(Ja,Ib,j) h(i,j)
*          (as normal 1-electron case )
*
* Part 3 : S(Ia,Ib) = S(Ia,Ib) + <Ia|i+|Ja>S(Ja,Ib,i)
*          Notice this is ROW scattering, to use our normal
*          routines we do the calculation as updates to transposed block
*           Loop over orbital i
*             transpose S(Ja,Ib,i) to S(Ib,Ja,i)
*             Update S(Ib,Ia) = S(Ib,Ia) + <Ia|i+|Ja>S(Ja,Ib,i)
*           End of loop over orbital I
*
*
* =============================
* The algorithm for IPM1 = -1
* =============================
*
* Sigma(Ia,Ib) = <Ia|j|Ja><Ib|i_+|Jb>C(Ja,Jb) h(j,i)
*
* Part 1 : C(         ) = <Ia|j|Ja>C(Ja,Jb)
*    Notice this is row scattering, we therefore proceed as
*    a) C(Jb,Ia,j) = <Ia|j|Ja> C(Jb,Ja)
*    b) C(Jb,Ia,j) => C(Ia,Jb,j)
*
* Part2 : S(Ia,Jb,i) = C(Ia,Jb,j) h(j,i)
*         (as normal 1-electron block)
*
* Part3 : S(Ia,Ib) = S(Ia,Ib) + <Ib|i_+|Jb>S(Ia,Jb,i)
*         (as normal 1-electron block )
*
      IF(IPM1.EQ.1) THEN
*. Transpose sigma block
        CALL trpmt3(SB,NIROWA,NIROWB,CSCR)
        CALL COPVEC(CSCR,SB,NIROWA*NIROWB)
      ELSE
*. Transpose C block
        CALL trpmt3(CB,NJROWA,NJROWB,CSCR)
        CALL COPVEC(CSCR,CB,NJROWA*NJROWB)
      END IF
*
      if (IPM1.eq.1) then
         call adtyp(NSXTPB,JTP,NGAS,ICELB,ISELB)
         call adtyp(NSXTPA,ITP,NGAS,ISELA,ICELA)
      else if (IPM1.eq.-1) then
         call adtyp(NSXTPA,JTP,NGAS,ICELA,ISELA)
         call adtyp(NSXTPB,ITP,NGAS,ISELB,ICELB)
      else
         write(6,*) ' IPM1 out of range '
         write(6,*) ' IPM1 = ',IPM1
         Call Abend2( ' Problem in dmk1one ..... ' )
      end if
*
      if (NTEST.ge.1000) then
         write(6,*) ' ITP,JTP ',ITP,JTP
      end if
*
      if (NSXTPB.eq.0 .or. NSXTPA.eq.0) goto 2000
*
* determine sym(J) for non-vanishing <Jb| j_+ |Ib>
*                                   (<Ja| j+ |Ia>)
* determine sym(I) for non-vanishing <Ia| i+ |Ja>
*                                   (<Ib| i_+ |Jb>)
      if (IPM1.eq.1) then
         JSM = IDBGMULT(ICCSMB,IADJSYM(ISCSMB))
         ISM = IDBGMULT(ISCSMA,IADJSYM(ICCSMA))
         NJORB = NOBPTS_B(JTP,JSM)
         NIORB = NOBPTS_U(ITP,ISM)
      else
         JSM = IDBGMULT(ICCSMA,IADJSYM(ISCSMA))
         ISM = IDBGMULT(ISCSMB,IADJSYM(ICCSMB))
         NJORB = NOBPTS_U(JTP,JSM)
         NIORB = NOBPTS_B(ITP,ISM)
      end if
      if (NTEST.ge.500) write(6,*) ' JSM ISM ', JSM,ISM
      if (NJORB.eq.0 .or. NIORB.eq.0) goto 2000
*
* start with batches of K strings according to sum(Kb) <Jb| j_+ |Kb><Kb|1|Ib>
*                                             (sum(Ka) <Ja| j+ |Ka><Ka|1|Ia>)
      KBOT = 1 - MAXK
      KTOP = 0
 700  continue
         KBOT = KBOT + MAXK
         KTOP = KTOP + MAXK
*
* determine mapping <Jb| j_+ |Kb>
*                  (<Ja| j+ |Ka>)
*
*  For internal affairs (unclear at the moment, to me at least)
*  (Bug fix from middle of march (2000?), set to 2 instead of 1)
         J12 = 2
         K12 = 2
*
         if (JFRST.eq.1) KFRST = 1
         if (IPM1.eq.1) then
            JBORB = IOBPTS_B(JTP,JSM)
            call adst_gas(JBORB,NJORB,ICCTPB,ICCSMB,JBGRP,KBOT,KTOP,
     &                    I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,
     &                    J12,K12,KKJACT,SCLFAC)
         else
            JBORB = IOBPTS_U(JTP,JSM)
            call adst_gas(JBORB,NJORB,ICCTPA,ICCSMA,JAGRP,KBOT,KTOP,
     &                    I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,
     &                    J12,K12,KKJACT,SCLFAC)
         end if
         if (NKBTC.le.0) then
           write(6,*) 'NKBTC after adst: ',NKBTC
           goto 2000
         end if
*
         if (NTEST.ge.1000) write(6,*) 'NKBTC = ',NKBTC
*
* finding out whether e.g. <Kb|Ib> is non-vanishing is done
* "behind the curtain".
*
         JFRST = 0
         KFRST = 0
* loop over batches of La, according to sum(La) <Ia| i+ |La><La|Ja>
*                                      (sum(Lb) <Ib| i_+ |Lb><Lb|Jb>)
         LBOT = 1 - MAXL
         LTOP = 0
 600     continue
            LBOT = LBOT + MAXL
            LTOP = LTOP + MAXL
*
*  determine mapping <Ia| i+ |La>
*                   (<Ib| i_+ |Lb>)
            I12 = 1
            L12 = 1
            if (IFRST.eq.1) LFRST = 1
            if (IPM1.eq.1) then
               IBORB = IOBPTS_U(ITP,ISM)
               call adst_gas(IBORB,NIORB,ISCTPA,ISCSMA,IAGRP,LBOT,LTOP,
     &                       I2,XI2S,MAXL,NLBTC,LEND,IFRST,LFRST,
     &                       I12,L12,LLIACT,SCLFAC)
            else
               IBORB = IOBPTS_B(ITP,ISM)
               call adst_gas(IBORB,NIORB,ISCTPB,ISCSMB,IBGRP,LBOT,LTOP,
     &                       I2,XI2S,MAXL,NLBTC,LEND,IFRST,LFRST,
     &                       I12,L12,LLIACT,SCLFAC)
            end if
            if (NLBTC.le.0) then
              write(6,*) 'NLBTC after adst: ',NLBTC
              goto 2000
            end if
*
            if (NTEST.ge.1000) write(6,*) 'NLBTC = ',NLBTC
*
*  obtain C(Ja,Kb,j) = sum(Jb) <Jb| j_+ |Kb><Kb|Ib> C(Ja,Jb)
*        (C(Jb,Ka,j) = sum(Ja) <Ja| j+ |Ka><Ka|Ia>)
*
            IF(LBOT.LE.0.OR.KBOT.LE.0) THEN
               WRITE(6,*) ' LBOT, KBOT = ',  LBOT, KBOT
               Call Abend2( ' Negative LBOT or KBOT in dmk1one ' )
            END IF
*
            if (IPM1.eq.1) then
               do JJORB = 1,NJORB,1
                  ICGOFF = 1 + (JJORB-1)*NKBTC*NLBTC
                  call matcg(CB,CSCR(ICGOFF),NJROWA,NLBTC,LBOT,NKBTC,
     &                       I1(1,JJORB),XI1S(1,JJORB))
               end do
            else
               do JJORB = 1,NJORB,1
                  ICGOFF = 1 + (JJORB-1)*NKBTC*NLBTC
                  call matcg(CB,SSCR(ICGOFF),NJROWB,NLBTC,LBOT,NKBTC,
     &                       I1(1,JJORB),XI1S(1,JJORB))
                  CALL COPVEC(SSCR(ICGOFF),CSCR(ICGOFF),NKBTC*NLBTC)
               end do
            end if
*
            if (NKBTC.ne.0.and.NLBTC.ne.0) then
               call geth1rel3(ONEINTS,H,ISM,ITP,JSM,JTP,
     &                        IRI,IRICO,IPM1,IUB,IRETURN,SCALEF,
     &                        INTSPACE,ENVIRO)
*
*  Transposition in BU case; structure is always adapted to "alpha-beta"
*  Transposition is required in both cases.
*  Reason is unclear at the moment!
c              if (IPM1.eq.1) then
                 call trpmt3(H,NIORB,NJORB,SSCRH1(1))
                 call copvec(SSCRH1(1),H,NIORB*NJORB)
c              end if
*
               if (NTEST.ge.500) then
               write(6,*) ' dmkp1one H block '
                  if ((IRI.eq.1.and.IRICO.eq.1).or.
     &                (IRI.eq.2.and.IRICO.eq.2)) then
                     write(6,*) ' Real integrals '
                     call wrtmat(H,NIORB,NJORB,NIORB,NJORB)
                  else if ((IRI.eq.1.and.IRICO.eq.2).or.
     &                     (IRI.eq.2.and.IRICO.eq.1)) then
                     write(6,*) ' Imaginary integrals '
                     call wrtmat(H,NIORB,NJORB,NIORB,NJORB)
                  else
                     write(6,*) 'Error in real/imaginary structure'
                     write(6,*) 'Stopping in dmk1one'
                     Call Abend()
                  end if
                  if (NTEST.ge.1000) then
                     write(6,*) 'Checking integral vector'
                     write(6,*) 'in dmk1one ...'
                     do I=1,NIORB*NJORB,1
                        write(6,*) 'H(',I,') = ',H(I)
                     end do
                  end if
               end if
            end if
*
            NLK = NLBTC*NKBTC
            FACTORC = 0.0D0
            FACTORAB = SCALEF
*
*  obtain S(Ja,Kb,i) = sum(j) C(Ja,Kb,j) h(j_,i)
*        (S(Jb,Ka,i) = sum(j) C(Jb,Ka,j) h(j,i_))
            call matml7(SSCR,CSCR,H,NLK,NIORB,NLK,NJORB,NJORB,
     &                  NIORB,FACTORC,FACTORAB,0)
      CALL MEMMAN(0,0,'CHECK ',0,'MATML7')
*
*  transpose matrix S(La,Kb,i) --> S(Kb,La,i)
*  determine S(Ia,Kb) = S(Ia,Kb) + <Ia| i+ |La> S(Ja,Kb,i)
*           (S(Ib,Ka) = S(Ib,Ka) + <Ib| i_+ |Lb> S(Jb,Ka,i))
*
*  write transpose matrix to CSCR array which not in use at this point
            if (IPM1.eq.1) then
               do IIORB = 1,NIORB,1
                  KOFF = (IIORB-1)*NKBTC*NLBTC + 1
                  call trpmt3(SSCR(KOFF),NLBTC,NKBTC,CSCR(KOFF))
                  call matcas(CSCR(KOFF),SB,NKBTC,NIROWB,KBOT,NLBTC,
     &                        I2(1,IIORB),XI2S(1,IIORB))
               end do
            else
               IF(NKBTC.LE.0.OR.NLBTC.LE.0) THEN
                 WRITE(6,*) ' NKBTC, NLBTC = ', NKBTC, NLBTC
                 Call Abend2( ' Negative batch sizes ' )
               END IF
               do IIORB = 1,NIORB,1
                  KOFF = (IIORB-1)*NKBTC*NLBTC + 1
                  call trpmt3(SSCR(KOFF),NLBTC,NKBTC,CSCR(KOFF))
                  call matcas(CSCR(KOFF),SB,NKBTC,NIROWA,KBOT,NLBTC,
     &                        I2(1,IIORB),XI2S(1,IIORB))
C                      MATCAS(CIN,COUT,NROWI,NROWO,IROWO1,NGCOL,ISCA,SCASGN
               end do
            end if
*
*  final result is on array SB (updated sigma block)
         if (LEND.eq.0) goto 600
      if (KEND.eq.0) goto 700
2000  continue
*. Transpose back to obtain usual sigma block
*
      IF(IPM1.EQ.1) THEN
*. Transpose sigma block
        call trpmt3(SB,NIROWB,NIROWA,CSCR)
        call COPVEC(CSCR,SB,NIROWA*NIROWB)
      ELSE IF( IPM1.EQ.-1) THEN
*. Transpose C block
        call trpmt3(CB,NJROWB,NJROWA,CSCR)
        call COPVEC(CSCR,CB,NJROWB*NJROWA)
      END IF
*
      call qexit('DMK1ONE')
      return
      end
