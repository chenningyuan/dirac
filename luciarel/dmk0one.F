CTF          call dmk0one(IBSM,IBTP,JBSM,JBTP,
CTF                       IJBGRP,NIA,
CTF                       NGAS,IBOC,JBOC,
CTF                       SB,CB,
CTF                       ADSXA,STSTSX,
CTF                       NOBPTS,IOBPTS,
CTF                       MAXI,MAXK,
CTF                       SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
CTF                       NSMOB,NSMST,MOCAA,
CTF                       SCLFAC,ONEINTS,INTSPACE,
CTF                       IRI,IRICO,IPM1,IUB,
CTF                       NTEST)
*
*
      subroutine dmk0one(ISCSM,      ! sym of sigma columns
     &                   ISCTP,      ! type of sigma columns
     &                   ICCSM,      ! sym of C columns
     &                   ICCTP,      ! type of C columns
     &                   IGRP,       ! string group of columns
     &                   NROW,       ! number of rows in S and C block
     &                   NGAS,       ! number of active sets
     &                   ISEL,       ! occupation in each active set for S block
     &                   ICEL,       ! occupation in each active set for C block
     &                   SB,         ! output sigma block
     &                   CB,         ! input C block
     &                   ADSXA,      ! Sym(a+a) = (Sym(a+),Sym(a))
     &                   STSTSX,     ! Sym of !st>,sx!st'> => sym of sx
*                                      so <st|sx|st'> is totally symmetric
     &                   NOBPTS,     ! # orbitals per type/symmetry
     &                   IOBPTS,     ! offsets for orbs of type/symmetry
     &                   MAXI,MAXK,  ! max # spectator/inner resol. strings
*                                      treated simultaneously
     &                   SSCR,CSCR,
     &                   I1,XI1S,
     &                   I2,XI2S,
     &                   H,          ! list of transformed integrals
     &                   NSMOB,      ! number of orbital syms
     &                   NSMST,      ! number of string syms
     &                   MOC,        ! use MOC method (inst. of N-1 res.)
     &                   SCLFAC,     ! scaling factor (Jeppe)
     &                   ONEINTS,    ! transformed onel-ints from WORK(KPICI1)
     &                   INTSPACE,   ! length of this vector
     &                   IRI,        ! real/imag. part of sigma vector
     &                   IRICO,      ! real/imag. part of coefficients
     &                   IPM1,       ! sign of Delta Mk (0 here!)
     &                   IUB,        ! UU part (1) or BB (2); 0 if IPM1.ne.0
     &                   IPRNT)     ! Print flag (IPRSIG in input)
*
* One electron excitations on column strings
*. If IH2TRM .ne. 0 then the diagonal and one-electron
*  excitations arising from the two body operator are also
* included
*
*  Real/imaginary logics:
*    S(r) = +h(r) . C(r) - h(i) . C(i)
*    iS(i) = +ih(i) . C(r) - h(r) . iC(i)
*     |      |       |     |       |
*    IRI    SCALEF  IRICO  SCALEF  IRICO
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
      INTEGER ADSXA(MXPOBS,2*MXPOBS),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*)
*  we don't have 2 sets here, as only the relevant one is passed
*  to the routine.
      INTEGER IOBPTS(MXPNGAS,*)
*.Specific Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      character*72 ENVIRO
*
      dimension ONEINTS(INTSPACE)
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),H(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
*.Local arrays
      DIMENSION ITP(MXPNGAS),JTP(MXPNGAS)
*
      call qenter('DMK0ONE')
* Type of single excitations that connects the two column strings
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.500)THEN
        WRITE(6,*)
        WRITE(6,*) ' ========================== '
        WRITE(6,*) ' Information from DMK0ONE '
        WRITE(6,*) ' ========================== '
        WRITE(6,*)
CTF     WRITE(6,*) ' DMK0ONE : MOC,IH2TRM ', MOC,IH2TRM
        WRITE(6,*) ' ISEL : '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' ICEL : '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF
*
      IPM1 = 0
*  (fortunately, +0 = -0)
*
*  The mapping of strings is only done the first time around, i.e.
*  in the first loop for the first Kbatch of strings.
*  After that, IFRST, JFRST (and KFRST, LFRST in dmk1one) are set
*  to zero. In adst_gas, the arc weights if the first index (I) is
*  .ne.0 and the occupation of strings if the second index (here J)
*  is .ne.0.
      IFRST = 1
      JFRST = 1
*
*  1) determine D(Ja,Ib,j+,i) = sum(Jb) <Ib|j_+ i_|Jb> . C(Ja,Jb)
*    a) find single excitations connecting spectator <Ib| = ISEL
*                                         and column |Jb> = ICEL
*
      CALL SXTYP_GAS(NSXTP,JTP,ITP,NGAS,MXPNGAS,ISEL,ICEL)
*
*    b) determine S(j_+ i_) = S*(|Ib>) x S*(|Jb>)
*       i.e. we need adjoint symmetries of the strings
*
      ISCASM = IADJSYM(ISCSM)
      ICCASM = IADJSYM(ICCSM)
      JISM = IDBGMULT(ICCSM,IADJSYM(ISCSM))
      if (JISM.eq.0) goto 1000
*
*  loop over types of single excitations
*
      do 900 JITP = 1,NSXTP,1
         JTYP = JTP(JITP)
         ITYP = ITP(JITP)
         if (NTEST.ge.1000) write(6,*) ' ITYP JTYP ', ITYP,JTYP
         do 800 JSM = 1,NSMOB,1
*
*  loop over symmetries of creator J, resulting in symmetry of creator I (ISM)
*  (sym of creator = adjsym of annihilator)
*
            ISM = IDBGMULT(IADJSYM(JSM),JISM)
            ISM = IADJSYM(ISM)
*  New intermediate strings will be accessed, so
            KFRST = 1
            if (ISM.eq.0) goto 800
            if (NTEST.ge.1000) write(6,*) ' JSM ISM ', JSM,ISM
            NJORB = NOBPTS(JTYP,JSM)
            NIORB = NOBPTS(ITYP,ISM)
            if (NJORB.eq.0.or.NIORB.eq.0) goto 800
*
* Using N-1 resolution method (there is also MOC)
*
* loop over partitionings of N-1 strings
* start with negative values up to zero
* MAXK is given as maximum number of strings per K batch
* inserting batches of N-1 strings according to
*   sum(Kb) <Ib| j_+ |Kb><Kb| i_ |Jb>
*
            KBOT = 1 - MAXK
            KTOP = 0
 700        continue
               KBOT = KBOT + MAXK
               KTOP = KTOP + MAXK
*  offset number for GAS ITYP and symmetry ISM, annihilator
               IBORB = IOBPTS(ITYP,ISM)
               I12 = 1
               K12 = 1
               if (IFRST.eq.1) KFRST = 1
*  determine mapping <Kb| i_ |Jb>, meaning <Jb| i_+ |Kb>, which is real
*
*  Input
*    IBORB : offset number for first spinor of ITYP and ISM
*    NIORB : number of spinors in block ITYP,ISM
*    ICCTP : Type of column vectors Jb
*    ICCSM : Symmetry of column vectors Jb
*    IGRP  : String group of columns
*    KBOT  : First string from K batch
*    KTOP  : Last string from K batch
*    MAXK  : max. no. of K strings per batch
*    IFRST : 1 at start, then 0
*    KFRST : 1 at start, then 0
*    I12   : for local scratch
*    K12   : for local scratch
*    SCLFAC: multiplied with sign in ads1_gas; scaling factor
*
*  Output
*    NKBTC : no. of K strings in present batch
*            in subroutine: NK, no. of K strings checked
*    I1    : I1(KSTR,IORB) depending on ket string K from batch and
*                          annihilator (transferred to creator) IORB
*            value is either 0 (no coupling of Jb and Kb)
*                            +-1 (coupling, strings are connected by IORB)
*    XI1S  : +- a sign might turn up for successful coupling
*    KEND  : if .eq.0 not all K strings have been accessed
*    KKIACT: if some nonvanishing excitations were found, KACT is set to 1,
*            else it is zero.
*
               call ADST_GAS(IBORB,NIORB,ICCTP,ICCSM,IGRP,KBOT,KTOP,
     &                       I1,XI1S,MAXK,NKBTC,KEND,IFRST,KFRST,
     &                       I12,K12,KKIACT,SCLFAC)
*
               IFRST = 0
               KFRST = 0
*. For testing
               IFRST = 1
               KFRST = 1
*
*. Determine mapping <Ib| j_+ |Kb>
*. Set up I2(KSTR,JORB)
               JBORB = IOBPTS(JTYP,JSM)
               J12 = 2
               if (JFRST.eq.1) KFRST = 1
               ONE = 1.0D0
               call ADST_GAS(JBORB,NJORB,ISCTP,ISCSM,IGRP,KBOT,KTOP,
     &                       I2,XI2S,MAXK,NKBTC,KEND,JFRST,KFRST,
     &                       J12,K12,KKJACT,ONE)
               JFRST = 0
               KFRST = 0
* for testing
               JFRST = 1
               KFRST = 1
*
*. loop over partitions of spectator strings
*
               NIPART = NROW/MAXI
               if (NIPART*MAXI.ne.NROW) NIPART = NIPART + 1
               do 701 IPART = 1,NIPART,1
                  IBOT = (IPART-1)*MAXI + 1
                  ITOP = MIN(IBOT+MAXI-1,NROW)
                  NIBTC = ITOP - IBOT + 1
*
*. obtain C(Ib,Kb,i_+) = sum(Jb) <Jb| i_+ |Kb> . C(Ib,Jb)
*  as an intermediate for D, as we still have a partitioning of
*  strings in K batches
*
                  do IIORB = 1,NIORB,1
                     ICGOFF = 1 + (IIORB-1)*NKBTC*NIBTC
                     call matcg(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                          NKBTC,I1(1,IIORB),XI1S(1,IIORB))
                  end do
*
*. obtain one-electron integrals (ISM,ITP,JSM,JTP)
*
                  if (NKBTC.ne.0) then
                     ZERO = 0.D0
                     call setvec(H,ZERO,NIORB*NJORB)
                     call geth1rel3(ONEINTS,H,ISM,ITYP,JSM,JTYP,
     &                       IRI,IRICO,IPM1,IUB,IRETURN,SCALEF,
     &                       INTSPACE,ENVIRO)
*
*  Some important comments:
*
                     if (NTEST.ge.1000) then
                        write(6,*) ' dmk0one H block '
                        if ((IRI.eq.1.and.IRICO.eq.1) .or.
     &                      (IRI.eq.2.and.IRICO.eq.2)) then
                           write(6,*) ' Real part'
                           call wrtmat(H,NIORB,NJORB,NIORB,NJORB)
                        else if ((IRI.eq.1.and.IRICO.eq.2) .or.
     &                           (IRI.eq.2.and.IRICO.eq.1)) then
                           write(6,*) 'Imaginary part'
                           call wrtmat(H,NIORB,NJORB,NIORB,NJORB)
                        else
                           write(6,*) 'Error in r/i structure'
                           Call Abend()
                        end if
                     end if
                  end if
*
*. multiply C(Ib,Kb,i_+) with h(j_,i_) to obtain Sscr(Ib,Kb,j_)
*
                  NIK = NIBTC*NKBTC
                  FACTORC = 0.0D0
                  FACTORAB = SCALEF
                  call matml7(SSCR,CSCR,H,NIK,NJORB,NIK,NIORB,NIORB,
     &                        NJORB,FACTORC,FACTORAB,0)
*. should be modified later, as we have real and imaginary parts of integrals
                  do JJORB = 1,NJORB,1
                     JBSOFF = 1 + (JJORB-1)*NKBTC*NIBTC
                     IF(IBOT.LE.0) THEN
                       WRITE(6,*) ' Start of rows (IBOT) le. 0', IBOT
                       Call Abend2( ' Problem in dmk0one .......' )
                     END IF
*
                     call matcas(SSCR(JBSOFF),SB,NIBTC,NROW,IBOT,NKBTC,
     &                           I2(1,JJORB),XI2S(1,JJORB))
                  end do
 701           continue
*              ^ end of loop over I batches
            if (KEND.eq.0) goto 700
*           ^ end of loop over K batches
 800     continue
*        ^ end of loop over index symmetries
 900  continue
*     ^ end of loop over types of single excitations
1000  continue
*     ^ JISM = 0 condition
2000  continue
      call qexit('DMK0ONE')
      return
      end
