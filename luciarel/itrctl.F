*  Control section for integral transformation
*
      subroutine itrctl(INTSPACE,ITRA_TYPE,INTITER,
     &                  ENVIRO,ISPINFREE,IFOCKEX,NORB_FULL,
     &                  NL2D,IUB2D,ISM2D,ISYM_T,SIGDEN_ROUTE,
     &                  F1,F2,IBOSYM,MCSCF,CIRUN)
*
*---------------------------------------------------
*  Written by Timo Fleig (including all subroutines)
*---------------------------------------------------
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "symm.inc"
#include "glbbas.inc"
#include "input.inc"
#include "clunit.inc"
#include "oper.inc"
#include "ctcc.inc"
*
      character*72 ENVIRO
      character*6 SIGDEN_ROUTE
      character*1 COORD
      logical MCSCF
      dimension INTEG(1000)
      dimension NORB_FULL(MXPIRR)
*
* ITRA_TYPE       Defines the type of transformation
*---------------------------------------------------
* ITRA_TYPE = 1 : Normal nonrelativistic transformation
* ITRA_TYPE = 2 : 2/4-component relativistic transformation for CI,
*                 all nonredundant integral lists
* ITRA_TYPE = 3 : 2/4-component relativistic transformation for CI,
*                 Excluding integrals with an odd number of bars
* ITRA_TYPE = 4 : 2/4-component relativistic transformation for MCSCF,
*                 or density matrix spinor optimization
*                 all nonredundant integral lists
* ITRA_TYPE = 5 : 2/4-component relativistic transformation for MCSCF,
*                 or density matrix spinor optimization
*                 Excluding integrals with an odd number of bars
*
*
      if (ENVIRO(1:6).eq.'MOLCAS') then
        FNONEAO = 'ONEINT'
        FNONEMO = 'TRAONE'
        INTIMP = 5
      else if (ENVIRO(1:6).eq.'DIRAC ') then
        DIRACONE = 'MRCONEE'
        DIRACTWO = 'MDCINT'
        INTIMP = 6
        if (MCSCF) INTIMP = 7
      end if
*
      call qenter('ITRCTRL')
*
      NTESTL = 00
      NTEST = max(IPRINT,NTESTL)
*
      WRITE(6,*) ' *************************************************'
      WRITE(6,*) ' *                                               *'
      WRITE(6,*) ' *  Transformer man                              *'
      WRITE(6,*) ' *  You run the show                             *'
      WRITE(6,*) ' *  Remote control                               *'
      WRITE(6,*) ' *  Direct the action with the push of a button  *'
      WRITE(6,*) ' *                                               *'
      WRITE(6,*) ' *                               Neil Young      *'
      WRITE(6,*) ' *                                               *'
      WRITE(6,*) ' *************************************************'
*
*  CI or MCSCF transformation?
*  Note: We will do the one-particle transformation always with
*        the full orbital space, it doesn't kill us. Afterwards,
*        a routine will select a matrix of integrals over the
*        active spaces for the GAS routines.
*
      NTRIAEL = ((SYMDIM(1)+SYMDIM(2))**2
     &          +(SYMDIM(1)+SYMDIM(2))) /2
*
      NTRSM_FULL = 0
      NDIM_FULL = 0
      NDIM_TOT = 0
      do ISYM=1,NIRR_PN,1
        NDIM_FULL = NDIM_FULL + NORB_FULL(ISYM)
        NDIM_TOT = NDIM_TOT + SYDI(ISYM)
      end do
      NTRSM_FULL = (NDIM_FULL**2 + NDIM_FULL) / 2
*
      call isetvc(SYMDIM_F,0,NIRR_DG)
      if (DOUGRP.lt.7.and.DOUGRP.ge.4) then
        I = 1
        do IRR=1,NIRR_PN,1
          if (IRR.eq.3.or.IRR.eq.5.or.IRR.eq.7) I = I + 1
          SYMDIM_F(I) = SYMDIM_F(I) + NORB_FULL(IRR)
        end do
      else if (DOUGRP.eq.7) then
        SYMDIM_F(1) = 2*NORB_FULL(1)
        SYMDIM_F(2) = 2*NORB_FULL(2)
      else if (DOUGRP.eq.8) then
        SYMDIM_F(1) = 2*NORB_FULL(1)
      end if 
*
      if (IFOCKEX.eq.1.and.ENVIRO(1:6).eq.'DIRAC ') then
        write(6,*) 'Problem in itrctl.'
        write(6,*) 'IFOCKEX is  ',IFOCKEX
        write(6,*) 'and environment is ',ENVIRO(1:6)
        write(6,*) 'This is not allowed. Omit FOCKEX variable.'
        Call Abend2( 'Quitting.' )
      end if
*
*
*  1) Static memory
*
*  Transformed integrals, 32 lists
      call findhigh(NGSOBPI,MXPIRR*MXPNGAS,MXGSDIM)
      write(6,*) 'Highest GAS dimension is  ',MXGSDIM
      call findhigh(SYMDIM(1),NIRR_DG,IMXDDIM)
      call findhigh(NORB_FULL(1),NIRR_PN,IMXDIM_F)
      if (IMXDDIM.gt.MXNSPIPDMGRP) then
        write(6,*) 'Spinor expansion too long for double dim. arrays.'
        write(6,*) 'MXNSPIPDMGRP is  ',MXNSPIPDMGRP
        write(6,*) 'IMXDDIM is ',IMXDDIM
        write(6,*) 'Adjust MXNSPIPDMGRP in mxpdim.inc or reconsider.'
        Call Abend2( 'Quitting.' )
      end if
      call findhigh(SYDI(1),NIRR_DG,IMXDIM)
      INTSPACE = (2*MXGSDIM)**2*NGASPI**2*32
      ITRALINE = (2*MXGSDIM)**2*NGASPI**2
*
      if (ITRA_TYPE.le.3) then
        call memman(KPICI1,INTSPACE,'ADDS  ',2,'IPI1  ')
*  Space for incorporated sf integrals and symmetry info
        if (ENVIRO(1:6).ne.'DIRAC ') then
          call memman(KPICAOT,ITRIASUM,'ADDS  ',2,'IPISFT')
*  Incorporated spin-orbit integrals
          call memman(KPICLX,NTRIAEL,'ADDS  ',2,'IPILX ')
          call memman(KPICLY,NTRIAEL,'ADDS  ',2,'IPILY ')
          call memman(KPICLZ,NTRIAEL,'ADDS  ',2,'IPILZ ')
          if (IFOCKEX.eq.1) then
            call memman(KPICLX_F,NTRSM_FULL,'ADDS  ',2,'IPILXF')
            call memman(KPICLY_F,NTRSM_FULL,'ADDS  ',2,'IPILYF')
            call memman(KPICLZ_F,NTRSM_FULL,'ADDS  ',2,'IPILZF')
          end if
        else
          call memman(KBOSYM,2*NDIM_TOT,'ADDS  ',1,'IBOSUM')
          call memman(KISREO,2*NDIM_TOT,'ADDS  ',1,'IISREO')
        end if
      else
        if (ENVIRO(1:6).eq.'DIRAC ') then
          call memman(KBOSYM,2*NDIM_TOT,'ADDS  ',1,'IBOSUM')
          call memman(KISREO,2*NDIM_TOT,'ADDS  ',1,'IISREO')
        end if
      end if
*
* When integral transformation section is generalized to traditional
* or sigden route, the allocation should be revised. At the moment:
      if (SIGDEN_ROUTE.eq.'SIGDEN') then
        print*,'NINT_CC',NINT_CC
        call memman(KINDEX,4*NINT_CC,'ADDS  ',1,'IINDEX')
      end if
*
************************************************
* Set mark for flushing local memory
************************************************
      call memman(IDUM,IDUM,'MARK  ',IDUM,'ITRCTL')
*
*  2) Local memory
*  Some index and GAS info
      call memman(KIGSOFF,2*NGASPI,'ADDL  ',1,'IGSOFF')
*  Scratch space, HELPMAT2
      call memman(KPICHE2,(2*IMXDDIM)**2,'ADDL  ',2,'IPIHE2')
      if (ENVIRO(1:6).eq.'DIRAC ') then
        NUM_ONEI = (2*NDIM_TOT)**2
        call memman(KDUMAR,NDIM_TOT*2,'ADDL  ',2,'DUMARR')
        call memman(KIDUMAR,NDIM_TOT*2,'ADDL  ',1,'IDUMAR')
        call memman(KDIONEE,2*NUM_ONEI,'ADDL  ',2,'IDIONE')
        NUM_TWEI = (2*NDIM_TOT)**4
        call memman(KDITWOE,2*NUM_TWEI,'ADDL  ',2,'IDITWO')
        call memman(KINVREO,2*NDIM_TOT+1,'ADDL  ',1,'INVREO')
        call memman(KKR,(2*NDIM_TOT+1),'ADDL  ',1,'IKR   ')
        call memman(KINDK,NUM_TWEI,'ADDL  ',2,'INDK  ')
        call memman(KINDL,NUM_TWEI,'ADDL  ',2,'INDL  ')
        call memman(KIOFF12,4*(2*NDIM_TOT)**2,'ADDL  ',1,'IOFF12')
        call memman(KNUMINT,4*(2*NDIM_TOT)**2,'ADDL  ',1,'NUMINT')
        call memman(KTPLOC,(2*NDIM_TOT)**2,'ADDL  ',1,'ITPLOC')
      else
        if (IFOCKEX.eq.1) then
          call memman(KPICHE2_F,(4*IMXDIM_F)**2,'ADDL  ',2,'IPIHEF')
          call memman(KPICHE3_F,(4*IMXDIM_F)**2,'ADDL  ',2,'IPIHE3')
        end if
      end if
*
*  for C2H reordering:
      IHELEN = 0
      do I=2,3,1
         IHELEN = IHELEN + (ITOTDM(I)**2 + ITOTDM(I)) / 2
      end do
*
      if (ITRA_TYPE.le.3) then
*
*******************************************************
*  Integral reading section
*******************************************************
*
        if ((IOPER.eq.1.or.IOPER.eq.3).and.IHAM12.le.2) then
          if (ENVIRO(1:6).eq.'MOLECU') then
            if (IFOCKEX.eq.0) then
              call rdint_molecu(WORK(KPICAOT),WORK(KPICHE2),
     &                          LU2,LU9,LUSGI1,ITRIASUM,IHELEN)
            else if (IFOCKEX.eq.1) then
              write(6,*) 'Case not implemented yet.'
              Call Abend2( 'Quitting.' )
            end if
          else if (ENVIRO(1:6).eq.'MOLCAS') then
            if (IFOCKEX.eq.0) then
              call rdint_molcas(WORK(KPICAOT),FNONEAO,LU2,NATOMS,
     &                          IFOCKEX)
            else if (IFOCKEX.eq.1) then
              call rd_traone(WORK(KPICAOT),LU2,FNONEMO,IPRINT)
            end if
          else if (ENVIRO(1:6).eq.'DIRAC ') then
            if (IOPER.eq.1) then
              write(6,*) 'Itrctl: Spinfree Dirac Hamiltonian!'
              call rdone_dirac(WORK(KDIONEE),LU2,DIRACONE,MORB,
     &                         NDIM_TOT,ISPINFREE,WORK(KBOSYM),
     &                         WORK(KIDUMAR),WORK(KDUMAR),NTEST)
              MORB_DIRAC = MORB
              if (IHAM12.ge.2) then
                IMOFFSET = (2*NDIM_TOT)**4
                call rdtwo_dirac(WORK(KDITWOE),
     &                           WORK(KDITWOE+IMOFFSET),
     &                           LU40,DIRACTWO,
     &                           NDIM_TOT,ISPINFREE,MORB,
     &                           WORK(KKR),WORK(KIOFF12),WORK(KNUMINT),
     &                           WORK(KTPLOC),
     &                           WORK(KINDK),WORK(KINDL),NTEST)
              end if
            end if
          else
            write(6,*) 'Problem in itrctl'
            write(6,'(A21,A6)') 'ENVIRO variable is  :',ENVIRO
            Call Abend2( 'Stopping execution, fatal error.' )
          end if
        end if
*
*  Read SO integrals to packed vector format (all symms on one triangle)
        if ((IOPER.eq.2.or.IOPER.eq.3).and.IHAM12.le.2) then
          if (IFOCKEX.eq.0) then
            if (ENVIRO(1:6).eq.'MOLECU') then
              call rdsoint_molecu(WORK(KPICLX),WORK(KPICLY),
     &                     WORK(KPICLZ),WORK(KPICHE2),
     &                     NTRIAEL,
     &                     LU69,LU70,LU71,LUSGI2,LUSGI3,LUSGI4)
            else if (ENVIRO(1:6).eq.'MOLCAS') then
              if (DOUGRP.ge.4.and.DOUGRP.le.6) then
                LENTRIA_SO = NTRIAEL
              else if (DOUGRP.ge.7.and.DOUGRP.le.8) then
                LENTRIA_SO = ITRIA(1) + ITRIA(2) + SYDI(1)*SYDI(2)
              end if
              call rdsoint_molcas(WORK(KPICLX),WORK(KPICLY),
     &                     WORK(KPICLZ),LENTRIA_SO,
     &                     LU69,LU70,LU71,LUSGI2,LUSGI3,LUSGI4)
            else if (ENVIRO(1:6).eq.'DIRAC ') then
              write(6,*) 'Itrctl: Full Dirac Hamiltonian!'
              call rdone_dirac(WORK(KDIONEE),LU2,DIRACONE,MORB,
     &                         NDIM_TOT,ISPINFREE,WORK(KBOSYM),
     &                         WORK(KIDUMAR),WORK(KDUMAR),NTEST)
              MORB_DIRAC = MORB
              if (IHAM12.ge.2) then
                IMOFFSET = (2*NDIM_TOT)**4
                call rdtwo_dirac(WORK(KDITWOE),
     &                           WORK(KDITWOE+IMOFFSET),
     &                           LU40,DIRACTWO,
     &                           NDIM_TOT,ISPINFREE,MORB,
     &                           WORK(KKR),WORK(KIOFF12),WORK(KNUMINT),
     &                           WORK(KTPLOC),
     &                           WORK(KINDK),WORK(KINDL),NTEST)
              end if
            end if
          else if (IFOCKEX.eq.1) then
            if (ENVIRO(1:6).eq.'MOLECU') then
              write(6,*) 'Case not implemented yet.'
              Call Abend2( 'Quitting.' )
c             stop 'Timo forced me to stop after reading di integrals.'
            else if (ENVIRO(1:6).eq.'MOLCAS') then
              call rdsoint_molcas(WORK(KPICLX_F),WORK(KPICLY_F),
     &                     WORK(KPICLZ_F),NTRSM_FULL,
     &                     LU69,LU70,LU71,LUSGI2,LUSGI3,LUSGI4)
*
            end if
          end if
        end if
        if (IOPER.lt.1.or.IOPER.gt.3) then
          write(6,*) 'Severe problem in itrctl!'
          write(6,*) 'Type of operator not specified correctly.'
          write(6,*) 'IOPER is  ',IOPER
          Call Abend2( 'Quitting.' )
        end if
      end if
*
*
*******************************************************
*  Distribution and transformation section
*******************************************************
*
*
      if (IFOCKEX.eq.1) then
*  External Fock matrix generation activated
*-------------------------------------------------------------
* Distribute to double group GAS lists (integrals transformed)
*-------------------------------------------------------------
        if (IOPER.eq.1.or.IOPER.eq.3) then
          if (IHAM12.le.2.and.ITRA_TYPE.le.3) then
            write(6,*)
            write(6,*) '*****************************************'
            write(6,*) '*                                       *'
            write(6,*) '* Distributing sf integrals to DG lists *'
            write(6,*) '*   (1-electron integrals)              *'
            write(6,*) '*                                       *'
            write(6,*) '*****************************************'
*
            call dist_1eints(WORK(KPICHE2_F),WORK(KPICAOT),
     &                       WORK(KPICI1),ITRALINE,IMXDIM_F,
     &                       IPRINT)
          end if
        end if
*
        if ((IOPER.eq.2.or.IOPER.eq.3).and.IHAM12.le.2) then
*
          write(6,*)
          write(6,*) '*****************************************'
          write(6,*) '*                                       *'
          write(6,*) '* Running tra_so for SO transformation  *'
          write(6,*) '*                                       *'
          write(6,*) '*****************************************'
*
*  Transform the spin-orbit integrals and distribute to dbg lists
          do 20 IOUTLST=1,NIRR_PN,1
             if (IOUTLST.le.2) then
                IDIMGRP = 1
             else
                IDIMGRP = 2
             end if
*
             if (DOUGRP.eq.7.or.DOUGRP.eq.8) then
               IDIMGRP = IOUTLST
             end if
*
             write(*,*) 'There are ',NIRR_PN,' integral lists.'
             write(*,*) 'Symmetry output list:'
             write(*,*) 'IOUTLST =  ',IOUTLST
             write(*,*) 'IDIMGRP is ',IDIMGRP
             write(*,*) 'SYMDIM(IDIMGRP) is ',SYMDIM(IDIMGRP)
             if (SYMDIM(IDIMGRP).eq.0) then
                write(6,*) 'Zero dimension group.'
                write(6,*) 'Skipping transformation section.'
                goto 777
             end if
*
*
*  2) Local memory (scratch)
*
           if ((DOUGRP.ge.4.and.DOUGRP.le.6).or.DOUGRP.eq.8) then
             IALLOCDIM = SYMDIM(IDIMGRP)
             IALLOCDIM_F = SYMDIM_F(IDIMGRP)
           else if (DOUGRP.eq.7) then
             call findhigh(SYMDIM,2,IALLDIM_CI)
             call findhigh(SYMDIM_F,2,IALLDIM_CI_F)
             IALLOCDIM = IALLDIM_CI
             IALLOCDIM_F = IALLDIM_CI_F
           end if
*  Untransformed sf AO integrals
           call memman(KPICAOD,IALLOCDIM_F**2,'ADDL  ',2,'IPISFD')
*  Transformed integrals, sf, 3xSO
           call memman(KPICTRA,4*IALLOCDIM**2,'ADDL  ',2,'IPITRA')
           call memman(KPICTRA_F,4*IALLOCDIM_F**2,'ADDL  ',2,'IPITRF')
*  Space for MOAO coefficients and adjoint
*  (This can be saved later on; inefficient)
           call memman(KPICCMO,2*IALLOCDIM**2,'ADDL  ',2,'IPICMO')
           call memman(KPICCMO_F,2*IALLOCDIM_F**2,'ADDL  ',2,'IPICMF')
           call memman(KPICCAD_F,2*IALLOCDIM_F**2,'ADDL  ',2,'IPICAD')
*  Some more scratch space ...
           call memman(KPICPRO,2*IALLOCDIM**2,'ADDL  ',2,'IPIPRO')
           call memman(KPICPRO_F,2*IALLOCDIM_F**2,'ADDL  ',2,'IPIPRF')
*  Unfolded SO integrals
           call memman(KPICSO,IALLOCDIM_F**2,'ADDL  ',2,'IPISO ')
*  Some index and GAS info
           call memman(KISMIND,IALLOCDIM,'ADDL  ',1,'ISMIND')
           call memman(KIGSIND,NIRR_PN*IALLOCDIM,'ADDL  ',1,'IGSIND')
*
           if (NTEST.ge.3) then
              write(6,*) 'Checking allocations in itrctl ...'
              write(6,*) 'Name       Length       Start address'
              write(6,*) 'Transformed integral lists: ',INTSPACE,KPICI1
              write(6,*) 'AOINT: ',IALLOCDIM_F**2,KPICAOD
              write(6,*) 'HELPMAT2: ',(IMXDDIM*2)**2,KPICHE2
              write(6,*) 'TRAINT: ',4*IALLOCDIM**2,KPICTRA
*  For correct dimensioning in subroutine:
           if (NIRR_PN.eq.4) then
             NORB_FULL3 = NORB_FULL(3)
             NORB_FULL2 = NORB_FULL(2)
           else if (NIRR_PN.eq.2) then
             NORB_FULL3 = 0
             NORB_FULL2 = NORB_FULL(2)
           else if (NIRR_PN.eq.1) then
             NORB_FULL3 = 0
             NORB_FULL2 = 0
           end if
              write(6,*) 'COEF: ',2*IALLOCDIM_F**2,KPICCMO_F
              write(6,*) 'ADJCO: ',2*IALLOCDIM_F**2,KPICCAD_F
              write(6,*) 'PROD: ',2*IALLOCDIM**2,KPICPRO
              write(6,*) 'Unfolded SO ints: ',IALLOCDIM_F**2,KPICSO
              write(6,*) 'ISMIND: ',IALLOCDIM,KISMIND
              write(6,*) 'IGSIND: ',NIRR_PN*IALLOCDIM,KIGSIND
              write(6,*) 'IGSOFF: ',2*NGASPI,KIGSOFF
           end if
           call memman(KITRAN,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
           if (NTEST.ge.2) then
              write(6,*) 'First free address after itrctl'
              write(6,*) KITRAN
           end if
*
*  Read coefficients either from file 19 (initial CI run)
*  or from coefficient file unit (LUMOUT) for optimized
*  coefficients.
           FNIPH = 'JOBIPH'
           ICHECK = 987654
           if (ENVIRO(1:6).eq.'MOLECU') then
             call rdcoef_molecu(WORK(KPICCAD_F),WORK(KPICCMO_F),
     &                   WORK(KPICHE2),WORK(KPICPRO_F),
     &                   ICHECK,SYMDIM_F(IDIMGRP),
     &                   NORB_FULL(IMXDIM_F),NORB_FULL2,NORB_FULL3,
     &                   NORB_FULL(1),
     &                   NIRR_PN,NIRR_DG,LU19,INTITER,IMXDIM_F,
     &                   IOUTLST,DOUGRP,ITRA_TYPE)
           else if (ENVIRO(1:6).eq.'MOLCAS') then
             call rdcoef_molcas(WORK(KPICCMO_F),WORK(KPICHE2_F),
     &                   WORK(KPICPRO_F),ICHECK,SYMDIM_F(IDIMGRP),
     &                   SYMDIM_F(1),NORB_FULL(IMXDIM_F),
     &                   NORB_FULL2,NORB_FULL3,
     &                   NORB_FULL(1),NIRR_PN,NIRR_DG,FNIPH,LU19,
     &                   INTITER,IMXDIM_F,IOUTLST,DOUGRP,
     &                   IFOCKEX,LU19,LUMOUT,ITRA_TYPE,
     &                   LUMOUTR,LUMOUTI)
           else
             write(6,'(A16,A6)') 'ENVIRO(1:6) is  ',ENVIRO(1:6)
             Call Abend2
     &       ( 'Unknown environment in oneltra. Quitting.' )
           end if
*
           do ICOORD=1,3,1
             if (ICOORD.eq.1) then
               KPICLIN_F = KPICLX_F
               COORD = 'X'
             else if (ICOORD.eq.2) then
               KPICLIN_F = KPICLY_F
               COORD = 'Y'
             else if (ICOORD.eq.3) then
               KPICLIN_F = KPICLZ_F
               COORD = 'Z'
             end if
*
             write(6,*)
             write(6,*) '-------------------------------------------'
             write(6,'(A6,A1,A19,I2)')
     &             'Doing ',COORD,' integrals for list',IOUTLST
             write(6,*) '-------------------------------------------'
*
             call tra_so(WORK(KPICI1),WORK(KPICHE2_F),
     &                   WORK(KPICAOD),WORK(KPICLIN_F),
     &                   WORK(KPICHE2),WORK(KPICTRA),
     &                   WORK(KPICTRA_F),
     &                   WORK(KPICCMO_F),
     &                   WORK(KPICCAD_F),
     &                   WORK(KPICPRO),WORK(KPICPRO_F),
     &                   WORK(KPICSO),
     &                   WORK(KISMIND),WORK(KIGSOFF),
     &                   WORK(KIGSIND),
     &                   ITRALINE,IMXDIM,IMXDIM_F,MXGSDIM,
     &                   IOUTLST,SYMDIM(IDIMGRP),SYMDIM_F(IDIMGRP),
     &                   NTRIAEL,INTITER,ITRA_TYPE,
     &                   INTSPACE,ENVIRO,COORD,
     &                   NORB_FULL,
     &                   NTRSM_FULL,IFOCKEX)
           end do
*
777        continue
           INTITER = INTITER + 1
20        continue
*         ^ end do over LISTS
        end if
*       ^ End if Hamiltonian type
*
        if (IHAM12.ge.2) then
          if (ITRA_TYPE.le.3) then
            call dist_2eints(NL2D,IUB2D,ISM2D,IPRINT)
          else if (ITRA_TYPE.ge.4) then
*
*  Transform two-electron integrals
*  (Has to be revised completely. When appropriate?)
*  When activated, a general solution for reading the
*  coefficients has to be found. At the moment, this
*  would cause trouble!
*
C           NLISTO = NIRR_PN
c           MAX_NLISTO = 20
*
*
* New allocation with help of memman needed here!
*
CTF
*  This routine will be needed, as soon as one goes beyond
*  an initial CI calculation, i.e. a natural spinor or an
*  MCSCF calculation
CTF
c           call tweltra(SOFT(IPOIS(5)),SOFT(IPOIS(6)),SOFT(IPOIS(7)),
c    &                   SOFT(IPOIS(8)),SOFT(IPOIS(9)),SOFT(IPOIS(10)),
c    &                   SOFT(IPOIS(11)),SOFT(IPOIS(12)),
c    &                   ITRA_TYPE,INTITER)
*
*  Names of variables:
*	5: I12O(MAX_NLISTO)
*	6: I34O(MAX_NLISTO)
*	7: I1234O(MAX_NLISTO)
*	8: ICCSMO(MAX_NLISTO)
*	9: IRIO(MXBLKO)
*      10: ISPITRA_RC(MAX_NLISTO)
*      11: ISPITRA(4,MAX_NLISTO)
*      12: LSPINOR(N_SPINR_EXP)
*
*
          end if
*         ^ End if transformation type
        end if
*       ^ End if Hamiltonian type
*
      else if (IFOCKEX.eq.0) then
        if (IHAM12.ge.2.and.ENVIRO(1:6).ne.'DIRAC ') then
          if (ITRA_TYPE.le.3) then
            write(6,*)
            write(6,*) '*****************************************'
            write(6,*) '*                                       *'
            write(6,*) '* Distributing sf integrals to DG lists *'
            write(6,*) '*   (only 2-electron integrals)         *'
            write(6,*) '*                                       *'
            write(6,*) '*****************************************'
*
            call dist_2eints(NL2D,IUB2D,ISM2D,IPRINT)
*
          else if (ITRA_TYPE.ge.4) then
*
*  Transform two-electron integrals
*  (Has to be revised completely. When appropriate?)
*  When activated, a general solution for reading the
*  coefficients has to be found. At the moment, this
*  would cause trouble!
*
C           NLISTO = NIRR_PN
c           MAX_NLISTO = 20
*
*
* New allocation with help of memman needed here!
*
             if (NIRR_PN.eq.4) then
               ITOTDM3 = ITOTDM(3)
               ITOTDM2 = ITOTDM(2)
             else if (NIRR_PN.eq.2) then
               ITOTDM3 = 0
               ITOTDM2 = ITOTDM(2)
             else if (NIRR_PN.eq.1) then
               ITOTDM3 = 0
               ITOTDM2 = 0
             end if
CTF
*  This routine will be needed, as soon as one goes beyond
*  an initial CI calculation, i.e. a natural spinor or an
*  MCSCF calculation
CTF
c           call tweltra(SOFT(IPOIS(5)),SOFT(IPOIS(6)),SOFT(IPOIS(7)),
c    &                   SOFT(IPOIS(8)),SOFT(IPOIS(9)),SOFT(IPOIS(10)),
c    &                   SOFT(IPOIS(11)),SOFT(IPOIS(12)),
c    &                   ITRA_TYPE,INTITER)
*
*  Names of variables:
*	5: I12O(MAX_NLISTO)
*	6: I34O(MAX_NLISTO)
*	7: I1234O(MAX_NLISTO)
*	8: ICCSMO(MAX_NLISTO)
*	9: IRIO(MXBLKO)
*      10: ISPITRA_RC(MAX_NLISTO)
*      11: ISPITRA(4,MAX_NLISTO)
*      12: LSPINOR(N_SPINR_EXP)
*
*
          end if
*         ^ End if transformation type
*
        end if
*       ^ End if Hamiltonian type
*
        if (IHAM12.le.2.and.ENVIRO(1:6).ne.'DIRAC ') then
*-------------------------------------------------------------
* Transform full SAO integrals; distribute to GAS lists
*-------------------------------------------------------------
*  Loop over integral output lists:
*  (which of course corresponds to the dg irreps of the spinors)
*  1: C1
*  2: C2, Cs
*  4: C2h
          do 10 IOUTLST=1,NIRR_PN,1
             if (IOUTLST.le.2) then
                IDIMGRP = 1
             else
                IDIMGRP = 2
             end if
*
*  Treat special cases Ci and C1
*
             if (DOUGRP.eq.7.or.DOUGRP.eq.8) then
                IDIMGRP = IOUTLST
             end if
             write(*,*) 'There are ',NIRR_PN,' integral lists.'
             write(*,*) 'Symmetry output list:'
             write(*,*) 'IOUTLST =  ',IOUTLST
             write(*,*) 'IDIMGRP is ',IDIMGRP
             write(*,*) 'SYMDIM(IDIMGRP) is ',SYMDIM(IDIMGRP)
             if (SYMDIM(IDIMGRP).eq.0) then
                write(6,*) 'Zero dimension group.'
                write(6,*) 'Skipping transformation section.'
                goto 666
             end if
*
*
*  2) Local memory (scratch)
*
*  Mark first free address
             call memman(IDUM,IDUM,'MARK  ',IDUM,'ITRCTL')
*
             IALLOCDIM = SYMDIM(IDIMGRP)
*  Untransformed sf AO integrals
             call memman(KPICAOD,IALLOCDIM**2,'ADDL  ',2,'IPISFD')
*  Transformed integrals, sf, 3xSO
             call memman(KPICTRA,4*IALLOCDIM**2,'ADDL  ',2,'IPITRA')
*  Space for MOAO coefficients and adjoint
*  (This can be saved later on; inefficient)
             call memman(KPICCMO,2*IALLOCDIM**2,'ADDL  ',2,'IPICMO')
             call memman(KPICCAD,2*IALLOCDIM**2,'ADDL  ',2,'IPICAD')
*  Some more scratch space ...
C            call memman(KPICDAO,IALLOCDIM**2,'ADDL  ',2,'IPIDAO')
             call memman(KPICPRO,2*IALLOCDIM**2,'ADDL  ',2,'IPIPRO')
*  Unfolded SO integrals
             call memman(KPICSOX,IALLOCDIM**2,'ADDL  ',2,'IPISOX')
             call memman(KPICSOY,IALLOCDIM**2,'ADDL  ',2,'IPISOY')
             call memman(KPICSOZ,IALLOCDIM**2,'ADDL  ',2,'IPISOZ')
*  Some index and GAS info
             call memman(KISMIND,IALLOCDIM,'ADDL  ',1,'ISMIND')
             call memman(KIGSIND,NIRR_PN*IALLOCDIM,'ADDL  ',1,'IGSIND')
*
             if (NTEST.ge.3) then
                write(6,*) 'Checking allocations in itrctl ...'
                write(6,*) 'Name       Length       Start address'
                write(6,*) 'Transformed integral lists: ',
     &                      INTSPACE,KPICI1
                write(6,*) 'Read AO integrals: ',ITRIASUM,KPICAOT
                write(6,*) 'Spin-orbit integrals x: ',NTRIAEL,KPICLX
                write(6,*) 'Spin-orbit integrals y: ',NTRIAEL,KPICLY
                write(6,*) 'Spin-orbit integrals z: ',NTRIAEL,KPICLZ
                write(6,*) 'AOINT: ',IALLOCDIM**2,KPICAOD
                write(6,*) 'HELPMAT2: ',(IMXDDIM*2)**2,KPICHE2
                write(6,*) 'TRAINT: ',4*IALLOCDIM**2,KPICTRA
                write(6,*) 'COEF: ',2*IALLOCDIM**2,KPICCMO
                write(6,*) 'ADJCO: ',2*IALLOCDIM**2,KPICCAD
                write(6,*) 'PROD: ',2*IALLOCDIM**2,KPICPRO
                write(6,*) 'Unfolded SO ints x: ',IALLOCDIM**2,KPICSOX
                write(6,*) 'Unfolded SO ints y: ',IALLOCDIM**2,KPICSOY
                write(6,*) 'Unfolded SO ints z: ',IALLOCDIM**2,KPICSOZ
                write(6,*) 'ISMIND: ',IALLOCDIM,KISMIND
                write(6,*) 'IGSIND: ',NIRR_PN*IALLOCDIM,KIGSIND
                write(6,*) 'IGSOFF: ',2*NGASPI,KIGSOFF
             end if
             call memman(KITRAN,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
             if (NTEST.ge.2) then
                write(6,*) 'First free address after itrctl'
                write(6,*) KITRAN
             end if
*
             FNIPH = 'JOBIPH'
             ICHECK = 987654
             if (ENVIRO(1:6).eq.'MOLECU') then
               call rdcoef_molecu(WORK(KPICCAD),WORK(KPICCMO),
     &                     WORK(KPICHE2),WORK(KPICPRO),
     &                     ICHECK,SYMDIM(IDIMGRP),
     &                     ITOTDM(IMXDIM),ITOTDM2,ITOTDM3,ITOTDM(1),
     &                     NIRR_PN,NIRR_DG,LU19,INTITER,IMXDIM,
     &                     IOUTLST,DOUGRP,ITRA_TYPE)
             else if (ENVIRO(1:6).eq.'MOLCAS') then
               call rdcoef_molcas(WORK(KPICCMO),WORK(KPICHE2),
     &                     WORK(KPICPRO),ICHECK,SYMDIM(IDIMGRP),
     &                     SYMDIM(1),ITOTDM(IMXDIM),ITOTDM2,ITOTDM3,
     &                     ITOTDM(1),NIRR_PN,NIRR_DG,FNIPH,LU19,
     &                     INTITER,IMXDIM,IOUTLST,DOUGRP,
     &                     IFOCKEX,LU19,LUMOUT,ITRA_TYPE,
     &                     LUMOUTR,LUMOUTI)
             else
               write(6,'(A16,A6)') 'ENVIRO(1:6) is  ',ENVIRO(1:6)
               Call Abend2
     &         ( 'Unknown environment in oneltra. Quitting.' )
             end if
*
             call oneltra(WORK(KPICI1),WORK(KPICAOT),
     &                    WORK(KPICLX),WORK(KPICLY),WORK(KPICLZ),
     &                    WORK(KPICAOD),WORK(KPICHE2),WORK(KPICTRA),
     &                    WORK(KPICCMO),WORK(KPICCAD),
     &                    WORK(KPICPRO),
     &                    WORK(KPICSOZ),WORK(KPICSOX),
     &                    WORK(KPICSOY),WORK(KISMIND),WORK(KIGSOFF),
     &                    WORK(KIGSIND),
     &                    ITRALINE,IMXDIM,MXGSDIM,
     &                    IOUTLST,SYMDIM(IDIMGRP),
     &                    NTRIAEL,INTITER,ITRA_TYPE,
     &                    INTSPACE,ENVIRO,IFOCKEX)
*
666          continue
             INTITER = INTITER + 1
10        continue
        end if
        if (ENVIRO(1:6).eq.'DIRAC ') then
          write(6,*)
          write(6,*) '********************************************'
          write(6,*) '* Distributing MOLTRA transformed          *'
          write(6,*) '* dirac one- (and two-) electron integrals *'
          write(6,*) '* using input integral reordering          *'
          write(6,*) '********************************************'
          IRIOFF = (2*NDIM_TOT)**2
          NDIM_TOT2 = NDIM_TOT*2
          if (SIGDEN_ROUTE.eq.'TRADIT') then
            call dist_ints_dirac(WORK(KPICHE2),WORK(KDIONEE),
     &                           WORK(KDIONEE+IRIOFF),
     &                           WORK(KPICI1),WORK(KBOSYM),
     &                           WORK(KISREO),WORK(KINVREO),
     &                           NL2D,IUB2D,ISM2D,
     &                           ITRALINE,NDIM_TOT,NDIM_TOT2,
     &                           IHAM12,ISPINFREE,IMXDIM,DOUGRP)
          else if (SIGDEN_ROUTE.eq.'SIGDEN') then
            call memman(KISTR_CA,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_CA')
            call memman(KISTR_CB,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_CB')
            call memman(KISTR_AA,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_AA')
            call memman(KISTR_AB,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STR_AB')
            call dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,
     &                            WORK(KLSOBEX),WORK(KT_CC),
     &                            WORK(KISTR_CA),WORK(KISTR_CB),
     &                            WORK(KISTR_AA),WORK(KISTR_AB),
     &                            WORK(KINDEX),
     &                            NDIM_TOT2,F1,F2,IBOSYM,
     &                            WORK(KLIBSOBEX),WORK(KLIHIND),
     &                            WORK(KLABEXTP),CIRUN,NTEST)
C     subroutine dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,
C                                 ISOBEX,T,
C    &                            ISTR_CA,ISTR_CB,
C                                 ISTR_AA,ISTR_AB,
C    &                            INDEX,
C                                 NASHT,F1,F2,IBOSYM,
C    &                            IBSPOBEX_TP,IHINDX,
C                                 LABEL,CIRUN,IPRNT)
c           stop 'AFTER DIST_INTS_SIGDEN.'
          end if
        end if
      else
        write(6,*) 'Problem in itrctl.'
        write(6,*) 'IFOCKEX is  ',IFOCKEX
        Call Abend2( 'Quitting.' )
      end if
*     ^ End if external Fock matrix active or not
*
*  Flush local memory to mark set at beginning of do loop
      call memman(IDUM,IDUM,'FLUSM ',IDUM,'ITRCTL')
*
      call qexit('ITRCTRL')
      return
      end
