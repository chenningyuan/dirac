      SUBROUTINE GSUBD1DG(RHO1,NACOB,
     &     ISCSMA,ISCSMB,ISCTPA,ISCTPB,
     &     ICCSMA,ICCSMB,ICCTPA,ICCTPB,
     &     IAGRP,IBGRP,JAGRP,JBGRP,
     &     NIROWA,NIROWB,NJROWA,NJROWB,
     &     NGAS,ISELA,ISELB,ICELA,ICELB,
     &     SB,CB,
     &     NOBPTS_U,NOBPTS_B,IOBPTS_U,IOBPTS_B,
     &     MAXL,MAXK,
     &     SSCR,CSCR,I1,XI1S,I2,XI2S,
     &     NSMOB,RHO1S,
     &     IPM1,SCALEF,SCLFAC,IPRNT)
*
* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen
*
*  modified for relativistic double group calculation
*  summer/fall of '97
*  april 1998, Timo Fleig
*
* =====
* Input
* =====
* RHO1  : One body density matrix to be calculated/updated
* NACOB : Number of active unbarred orbitals
* NACOB2 : Number of active barred orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IAGRP,IBGRP : String group of S columns
* JAGRP,JBGRP : String group of C columns
* NIROWA,NIROWB : Number of rows in S block
* NJROWA,NJROWB : Number of rows in C block
* NGAS : Number of active spaces
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* SB   : Input S block
* CB   : Input C block
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* NSMOB : Number of symmetries of spinors
* MAXL   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
* IPM1 : distinguish barred-unbarred from unbarred-barred part
*         +1 : operator i+ j_ (unbarred-barred)
*         -1 : operator i_+ j (barred-unbarred)
*
* ======
* Output
* ======
* RHO1 : Calculated/updated density block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXLJ*MAXL*MAXK, where MAXLJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : MAXK*Max number of orbitals of given type and symmetry
* I2, XI2S   : MAXK*Max number of orbitals of given type and symmetry
*              type and symmetry
* RHO1S : Space for one electron density
*
* Jeppe Olsen, Winter of 1991
* Updated for GAS , August '95
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
*. General input
      INTEGER NOBPTS_U(MXPNGAS,*),NOBPTS_B(MXPNGAS,*)
      integer IOBPTS_U(MXPNGAS,*),IOBPTS_B(MXPNGAS,*)
*.Input
      INTEGER ISELA(NGAS),ICELA(NGAS)
      integer ISELB(NGAS),ICELB(NGAS)
*
#include "symm.inc"
*
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO1(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),RHO1S(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*)
      DIMENSION I2(MAXK,*),XI2S(MAXK,*)
*
*.Local arrays
      DIMENSION ITP(16),JTP(16)
*
      NTESTL = 0000
      NTEST = max(NTESTL,IPRNT)
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' ===================='
         WRITE(6,*) ' GSUBD1DG in action '
         WRITE(6,*) ' ===================='
         WRITE(6,*)
         WRITE(6,*) ' Occupation of active left unbarred strings '
         CALL IWRTMA(ISELA,1,NGAS,1,NGAS)
         WRITE(6,*) ' Occupation of active Right unbarred strings '
         CALL IWRTMA(ICELA,1,NGAS,1,NGAS)
         WRITE(6,*) ' Occupation of active left barred strings '
         CALL IWRTMA(ISELB,1,NGAS,1,NGAS)
         WRITE(6,*) ' Occupation of active Right barred strings '
         CALL IWRTMA(ICELB,1,NGAS,1,NGAS)
*
         WRITE(6,*) ' GSBBD1DG, sclfac ',SCLFAC
      END IF
*
      IFRST = 1
      JFRST = 1
*
*  determine C(Ja,Ib,j) = sum(Jb) <Jb| j_+ |Ib> C(Ja,Jb)
*  (if IPM1 = -1 determine C(Ia,Jb,j) = sum(Ja) <Ja| j+ |Ia> C(Ja,Jb) )
      if (IPM1.eq.1) then
         call adtyp(NSXTPB,JTYP,NGAS,ICELB,ISELB)
         call adtyp(NSXTPA,ITYP,NGAS,ISELA,ICELA)
      else if (IPM1.eq.-1) then
         call adtyp(NSXTPA,JTYP,NGAS,ICELA,ISELA)
         call adtyp(NSXTPB,ITYP,NGAS,ISELB,ICELB)
      else
         write(6,*) ' IPM1 out of range '
         write(6,*) ' IPM1 = ',IPM1
         goto 2000
      end if
      if (NSXTPA.eq.0 .or. NSXTPB.eq.0) goto 2000
*
* determine sym(J) for non-vanishing <Jb| j_+ |Ib>
*                                   (<Ja| j+ |Ia>)
* determine sym(I) for non-vanishing <Ia| i+ |Ja>
*                                   (<Ib| i_+ |Jb>)
      if (IPM1.eq.1) then
* transpose sigma block
         call trpmt3(SB,NIROWA,NIROWB,CSCR)
         call copvec(CSCR,SB,NIROWA*NIROWB)
         JSM = IDBGMULT(ICCSMB,IADJSYM(ISCSMB))
         ISM = IDBGMULT(ISCSMA,IADJSYM(ICCSMA))
         NJORB = NOBPTS_B(JTYP,JSM)
         NIORB = NOBPTS_U(ITYP,ISM)
      else
* transpose C block
         call trpmt3(CB,NJROWA,NJROWB,CSCR)
         call copvec(CSCR,CB,NJROWA*NJROWB)
         JSM = IDBGMULT(ICCSMA,IADJSYM(ISCSMA))
         ISM = IDBGMULT(ISCSMB,IADJSYM(ICCSMB))
         NJORB = NOBPTS_U(JTYP,JSM)
         NIORB = NOBPTS_B(ITYP,ISM)
      end if
      if (NTEST.ge.1000)
     &     write(6,*) ' JSM ISM ', JSM,ISM
*
* start with batches of K strings according to sum(Kb) <Jb| j_+ |Kb><Kb|1|Ib>
*                                             (sum(Ka) <Ja| j+ |Ka><Ka|1|Ia> )
      if (NJORB.eq.0 .or. NIORB.eq.0) goto 2000
      KBOT = 1 - MAXK
      KTOP = 0
 700  continue
      KBOT = KBOT + MAXK
      KTOP = KTOP + MAXK
      J12 = 2
      K12 = 2
      if (JFRST.eq.1) KFRST = 1
*
*  determine mapping <Jb| j_+ |Kb>
*                   ( <Ja| j+ |Ka> )
      if (IPM1.eq.1) then
         JBORB = IOBPTS_B(JTYP,JSM)
         if (NTEST.ge.1000) write(6,*) 'offset right orb ',JBORB
         call adst_gas(JBORB,NJORB,ICCTPB,ICCSMB,JBGRP,KBOT,KTOP,
     &        I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,
     &        J12,K12,KKJACT,SCLFAC)
      else if (IPM1.eq.-1) then
         JBORB = IOBPTS_U(JTYP,JSM)
         if (NTEST.ge.1000) write(6,*) 'offset right orb ',JBORB
         call adst_gas(JBORB,NJORB,ICCTPA,ICCSMA,JAGRP,KBOT,KTOP,
     &        I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,
     &        J12,K12,KKJACT,SCLFAC)
      end if
*
      JFRST = 0
      KFRST = 0
* loop over batches of La, according to sum(La) <Ia| i+ |La><La|Ja>
      LBOT = 1 - MAXL
      LTOP = 0
 600  continue
      LBOT = LBOT + MAXL
      LTOP = LTOP + MAXL
      I12 = 1
      L12 = 1
      if (IFRST.eq.1) LFRST = 1
*
*  determine mapping <Ia| i+ |La>
*                   (<Ib| i_+ |Lb>)
      if (IPM1.eq.1) then
         IBORB = IOBPTS_U(ITYP,ISM)
         if (NTEST.ge.1000) write(6,*) 'offset left orb ',IBORB
         call adst_gas(IBORB,NIORB,ISCTPA,ISCSMA,IAGRP,LBOT,LTOP,
     &        I2,XI2S,MAXL,NLBTC,LEND,IFRST,LFRST,
     &        I12,L12,LLIACT,SCLFAC)
      else
         IBORB = IOBPTS_B(ITYP,ISM)
         if (NTEST.ge.1000) write(6,*) 'offset left orb ',IBORB
         call adst_gas(IBORB,NIORB,ISCTPB,ISCSMB,IBGRP,LBOT,LTOP,
     &        I2,XI2S,MAXL,NLBTC,LEND,IFRST,LFRST,
     &        I12,L12,LLIACT,SCLFAC)
      end if
*
*  obtain CSCR(Ja,Kb,j) = sum(Jb) <Jb| j_+ |Kb><Kb|Ib> C(Ja,Jb)
*        (CSCR(Jb,Ka,j) = sum(Ja) <Ja| j+ |Ka><Ka|Ia> C(Jb,Ja) )
*
*  obtain SSCR(La,Ib,i) = sum(Ia) <Ia| i+ |La><La|Ja> C(Ib,Ia)
*        (SSCR(Lb,Ia,i) = sum(Ib) <Ib| i_+ |Lb><Lb|Jb> C(Ia,Ib) )
      if (NTEST.ge.1000) write(6,*) 'IPM1 =  ',IPM1
      if (IPM1.eq.1) then
         do JJORB = 1,NJORB,1
            ICGOFF = 1 + (JJORB-1)*NKBTC*NLBTC
            call matcg(CB,CSCR(ICGOFF),NJROWA,NLBTC,LBOT,
     &           NKBTC,I1(1,JJORB),XI1S(1,JJORB))
         end do
         do IIORB = 1,NIORB,1
            ISGOFF = 1 + (IIORB-1)*NKBTC*NLBTC
            call matcg(SB,SSCR(ISGOFF),NIROWB,NKBTC,KBOT,
     &           NLBTC,I2(1,IIORB),XI2S(1,IIORB))
         end do
      else
         do JJORB = 1,NJORB,1
            ICGOFF = 1 + (JJORB-1)*NKBTC*NLBTC
            call matcg(CB,CSCR(ICGOFF),NJROWB,NLBTC,LBOT,
     &           NKBTC,I1(1,JJORB),XI1S(1,JJORB))
         end do
         do IIORB = 1,NIORB,1
            ISGOFF = 1 + (IIORB-1)*NKBTC*NLBTC
            call matcg(SB,SSCR(ISGOFF),NIROWA,NKBTC,KBOT,
     &           NLBTC,I2(1,IIORB),XI2S(1,IIORB))
         end do
      end if
*
      NKI = NKBTC*NLBTC
      FACTORC = 0.0D0
      FACTORAB = SCALEF
*
* Check transposing scheme before matrix multiply:
*  if IPM1.eq.1   SB has been transposed
*  if IPM1.eq.-1  CB has been transposed
* transpose back where required
*
      if (IPM1.eq.1) then
         call matml7(RHO1S,SSCR,CSCR,NIORB,NJORB,NIORB,NKI,
     &        NKI,NJORB,FACTORC,FACTORAB,0)
      else if (IPM1.eq.-1) then
         call trpmt3(SSCR,NIROWB,NIROWA,SB)
         call copvec(SB,SSCR,NIROWA*NIROWB)
         call trpmt3(CSCR,NJROWB,NJROWA,CB)
         call copvec(CB,CSCR,NJROWA*NJROWB)
         call matml7(RHO1S,SSCR,CSCR,NIORB,NJORB,NIORB,NKI,
     &        NKI,NJORB,FACTORC,FACTORAB,0)
      end if
*
*  Scatter out to complete matrix
*  We have 2 sets of spinors, and we are calculating density matrix
*  blocks for each of the four resulting blocks separately. This means
*  that the offsets and numbers of spinors are absolute values in each
*  of the blocks in question. The number of active orbitals corresponds
*  to the number of spinors in the row index of the block we are
*  calculating. This is taken care of in the form of separate variables
*  passed from above.
*
      do 900 JJORB = 1,NJORB,1
         JORB = JBORB - 1 + JJORB
         do 800 IIORB = 1,NIORB,1
            IORB = IBORB - 1 + IIORB
            if (NTEST.ge.1000) then
              write(6,'(A13,1X,F10.8,2X,A16,I5,A3,F10.8)')
     &              'Adding RHO1S =',
     &              RHO1S((JJORB-1)*NIORB+IIORB),
     &              'to element RHO1(',
     &              (JORB-1)*NACOB+IORB,')= ',
     &              RHO1((JORB-1)*NACOB+IORB)
            end if
            RHO1((JORB-1)*NACOB+IORB) = RHO1((JORB-1)*NACOB+IORB)
     &           + RHO1S((JJORB-1)*NIORB+IIORB)
            if (NTEST.ge.1000) then
                write(6,'(3X,A3,F16.14)') ' = ',
     &             RHO1((JORB-1)*NACOB+IORB)
            end if
 800     continue
 900  continue
*
      if (KEND.eq.0) goto 700
      if (LEND.eq.0) goto 600
*
 2000 continue
      return
      end
