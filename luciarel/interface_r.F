!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Distribute external Fock type one-electron integrals
*  to double group lists
*  Pick active spin-orbit integrals
*
      subroutine dist_1eints(HELPMAT,ONEINT,SPININTS,ITRALINE,
     &                       IMXDIM,IPRINT)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "symm.inc"
#include "onel.inc"
*
      dimension ISTART(MXPIRR),IDIM(MXPIRR),ILEN(MXPIRR)
      dimension ONEINT(ITRIASUM),HELPMAT((4*IMXDIM)**2),
     &          SPININTS(ITRALINE,32)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
*  Start addresses and lengths:
*
      ISTART(1) = 1
      do ISYM=1,NIRR_PN,1
        ILEN(ISYM) = (SYDI(ISYM)**2 + SYDI(ISYM)) / 2
        if (ISYM.ne.NIRR_PN)
     &      ISTART(ISYM+1) = ISTART(ISYM) + ILEN(ISYM)
      end do
*
      call findhigh(NGSOBPI,MXPIRR*MXPNGAS,MXGSDIM)
      do LIST=1,NIRR_PN,1
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'DIST1E')
        if (LIST.le.2) then
          IDIMGRP = 1
        else
          IDIMGRP = 2
        end if
*
        if (DOUGRP.eq.7.or.DOUGRP.eq.8) then
          IDIMGRP = LIST
        end if
*
        if (NTEST.ge.1) then
          write(6,*) 'Doing DG list    ',LIST
          write(6,*) 'Dimension group  ',IDIMGRP
          write(6,*) 'SYMDIM(IDIMGRP)  ',SYMDIM(IDIMGRP)
        end if
*
        ITRALEN = SYMDIM(IDIMGRP)**2
      call memmar(KPICTRA,ITRALEN,'ADDL  ',2,'IPITRA')
        call setvec(WORK(KPICTRA),0.D0,ITRALEN)
      call memmar(KIGSIND,NIRR_PN*SYMDIM(IDIMGRP),
     &              'ADDL  ',2,'IGSIND')
      call memmar(KISMIND,SYMDIM(IDIMGRP),'ADDL  ',2,'ISMIND')
      call memmar(KIGSOFF,2*NGASPI,'ADDL  ',2,'IGSOFF')
*
        if (NTEST.ge.3) then
          write(6,*) 'KPICTRA ... ',KPICTRA,ITRALEN
          write(6,*) 'KIGSIND ... ',KIGSIND,NIRR_PN*SYMDIM(IDIMGRP)
          write(6,*) 'KISMIND ... ',KISMIND,SYMDIM(IDIMGRP)
          write(6,*) 'KIGSOFF ... ',KIGSOFF,2*NGASPI
        end if
*
        call evenodd(IEOLST,LIST)
        if (IEOLST.eq.1) then
          ISTA = LIST
        else if (IEOLST.eq.2) then
          ISTA = LIST - 1
        end if
        IEND = ISTA + 1
        ICT = 1
        do ISYM=ISTA,IEND,1
          IMEDSYM = ISYM
          if (DOUGRP.eq.7.or.DOUGRP.eq.8) then
            if (ISYM.eq.1.and.IEOLST.eq.1) IMEDSYM = 1
            if (ISYM.eq.2.and.IEOLST.eq.1) IMEDSYM = 1
            if (ISYM.eq.1.and.IEOLST.eq.2) IMEDSYM = 2
            if (ISYM.eq.2.and.IEOLST.eq.2) IMEDSYM = 2
          end if
          call unfold(HELPMAT,ONEINT(ISTART(IMEDSYM)),
     &                SYDI(IMEDSYM),1,ILEN(IMEDSYM),1)
          ILENSQ = ILEN(IMEDSYM)**2
          IROWOFF = ICT
          ICOLOFF = ICT
          call vectomat(HELPMAT,WORK(KPICTRA),ILENSQ,
     &                  SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &                  SYDI(IMEDSYM),SYDI(IMEDSYM),
     &                  IROWOFF,ICOLOFF)
          if (NTEST.ge.3) then
            write(6,*)
            write(6,*) 'Checking WORK(KPICTRA) after vectomat'
            write(6,*) 'LIST =  ',LIST
            write(6,*) 'IMEDSYM =  ',IMEDSYM
            IACDM = SYMDIM(IDIMGRP)
            call wrtmat(WORK(KPICTRA),IACDM,IACDM,IACDM,IACDM)
          end if
          ICT = ICT + SYDI(IMEDSYM)
        end do
        IACSYM = 1
        if (LIST.gt.2) IACSYM = 3
        INTLIST = 2*LIST - 1
        call intgas3(SPININTS(1,INTLIST),WORK(KPICTRA),
     &               NGSOBPI(1,1),NGASPI,NIRR_PN,
     &               INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &               WORK(KISMIND),WORK(KIGSOFF),WORK(KIGSIND),
     &               IACSYM,LIST,MXGSDIM,ISPSF_MO,DOUGRP)
*
*  Copy block lengths of IUB and GAS blocks to unused lists.
*  The imaginary lists 2,4,6,8 are needed in the integral
*  fetcher, where Y integrals are added to these blocks !!
        do IGAS=1,MXPNGAS,1
          do JGAS=1,MXPNGAS,1
            do I=1,4,1
              IPUBPGS(I,IGAS,JGAS,INTLIST+1)
     &        = IPUBPGS(I,IGAS,JGAS,INTLIST)
            end do
          end do
        end do
*
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'DIST1E')
      end do
*
      if (NTEST.ge.10) then
        do LIST=1,2*NIRR_PN,1
          call evenodd(IEOLST,LIST)
          write(6,*)
          write(6,*) 'Check in dist_1eints'
          write(6,*) '--------------------'
          if (IEOLST.eq.1) then
            write(6,*) 'Writing real integrals of list  ',LIST
          else if (IEOLST.eq.2) then
            write(6,*) 'Writing imaginary integrals of list  ',LIST
          end if
          do INT=1,ITRALINE,1
            write(6,'(I6,1F18.14)') INT,SPININTS(INT,LIST)
          end do
        end do
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Distribute externally transformed (MOLTRA) dirac one-electron
*  integrals to double group lists
*   
*----------------------------------------------
*  Timo Fleig, Oct. 2000
*----------------------------------------------
*
      subroutine dist_1eints_dirac(HELPMAT,ONEINT,IRI,INPOFF,
     &                             SPININTS,ISREO,
     &                             ITRALINE,ITOTDIM,IDIM,
     &                             LIST,NTEST)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "symm.inc"
*
      dimension ONEINT(ITOTDIM,ITOTDIM),HELPMAT(IDIM,IDIM),
     &          SPININTS(ITRALINE,32)
      dimension ISREO(IDIM)
*
      NTESTL = 00
      NTEST = max(NTESTL,NTEST)
*
*  Zero out help matrix
      ZERO = 0.D0
      do I=1,IDIM,1
        call setvec(HELPMAT(1,I),ZERO,IDIM)
      end do
*
      do IR=1,IDIM,1
        do IL=1,IDIM,1
          HELPMAT(ISREO(IL),ISREO(IR)) 
     &                   = ONEINT(IL+INPOFF,IR+INPOFF)
        end do
      end do
*
      if (NTEST.ge.1) then
        write(6,*) 'Checking reordered integrals:'
        if (IRI.eq.1) write(6,*) 'Real list:'
        if (IRI.eq.2) write(6,*) 'Imag. list:'
        call wrtmat(HELPMAT,IDIM,IDIM,IDIM,IDIM)
      end if
*
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'DISTDI')
      call memmar(KIGSIND,NIRR_PN*SYMDIM(IDIMGRP),
     &              'ADDL  ',2,'IGSIND')
      call memmar(KISMIND,SYMDIM(IDIMGRP),'ADDL  ',2,'ISMIND')
      call memmar(KIGSOFF,2*NGASPI,'ADDL  ',2,'IGSOFF')
*
      IACSYM = 1
      if (LIST.gt.2) IACSYM = 3
      INTLIST = 2*LIST - 2 + IRI
      call findhigh(NGSOBPI,MXPIRR*MXPNGAS,MXGSDIM)
      call intgas3(SPININTS(1,INTLIST),HELPMAT,
     &             NGSOBPI(1,1),NGASPI,NIRR_PN,
     &             INTLIST,IDIM,ITOTDM(IACSYM),
     &             WORK(KISMIND),WORK(KIGSOFF),WORK(KIGSIND),
     &             IACSYM,LIST,MXGSDIM,ISPSF_MO,DOUGRP)
*
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'DISTDI') 
*
      write(6,'(A,I3,A,I6)') 
     &      ' List no. ',LIST,' of 1-el. ints processed. DIM ',IDIM
      if (NTEST.ge.5) then
        call evenodd(IEOLST,LIST)
        write(6,*)
        write(6,*) 'Check in dist_1eints_dirac'
        write(6,*) '--------------------------'
        write(6,*)
        write(6,*) 'This is total list ',INTLIST
        if (IEOLST.eq.1) then
          write(6,*) 'Writing real integrals of list  ',LIST
        else if (IEOLST.eq.2) then
          write(6,*) 'Writing imaginary integrals of list  ',LIST
        end if
        do INT=1,ITRALINE,1
          write(6,'(I6,1F18.14)') INT,SPININTS(INT,INTLIST)
        end do
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE DIST_2EINTS(NLIST,IBU,ISM,IPRTRA)
*
* Distribute nonrelativistic integrals into
* double group format
*
*---------------------------------------------
* Jeppe Olsen and Timo Fleig, 1997/98
*---------------------------------------------
*
* Distribute DIRAC integrals from symmetry reduced list
* on WORK(KINT2) to double group lists
*
*---------------------------------------------
* Timo Fleig, October 2000
*---------------------------------------------
*
* Input :
*
* NLIST : Number of two electron integral lists
*
* IBU :
*
* IBU(INDEX,ILIST) = 1 => Index INDEX  in list ILIST is   barred
*                  = 2 => Index INDEX  in list ILIST is unbarred
*
* ISM :
*
* ISM(1,ILIST) = 1 => Symmetry between index 1 and 2 in list ILIST
*              = 0 => No Symmetry between index 1 and 2 in list ILIST
*
* ISM(2,ILIST) = 1 => Symmetry between index 3 and 4 nd  in list ILIST
*              = 0 => No Symmetry between index 3 and 4  in list ILIST
*
* ISM(3,ILIST) = 1 => Symmetry between index (12) and (34)  in list ILIST
*              = 0 => No Symmetry between  index (12) and (34) in  ILIST
*
* Note : at the moment all symmetry flags are supposed to be zero !!
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "cintfo.inc"
#include "glbbas.inc"
*
      INTEGER IBU(4,*), ISM(3,*)
*
      dimension NREOSM(2,MXPNGAS)
*
      NTESTL = 0
      NTEST = max(IPRTRA,NTESTL)
*
      IPNTGRP_SAVE = PNTGRP
      IDUM = 1
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'DIST2E')
*. Block holding largest block of integrals
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
        DO IOBSM = 1, NSMOB
          MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
        END DO
      END DO
      INTBLK = MXTSOB_L**4
      call memmar(KLINT2,INTBLK,'ADDL  ',2,'KLINT2')
*. Loop over block of normal integrals, fetch them and dispose them
      DO LSYM = 1, NSMOB/2
       DO KSYM = 1, NSMOB/2
        DO JSYM = 1, NSMOB/2
*. Point group symmetry
          PNTGRP = 1
          CALL SYMCOM_REL(3,1,LSYM,KSYM,LKSYM)
          CALL SYMCOM_REL(2,1,LKSYM,JISYM,1)
          CALL SYMCOM_REL(2,1,JSYM,ISYM,JISYM)
*. ISYM, JSYM, KSYM and LSYM are now known
          IF(NTEST.GE.3) THEN
            WRITE(6,*) ' ISYM, JSYM, KSYM, LSYM',
     &      ISYM,JSYM,KSYM,LSYM
          END IF
          PNTGRP = IPNTGRP_SAVE
*. Total number of orbitals of these  symmetries
          NITOT = NGSSHT(ISYM)
          NJTOT = NGSSHT(JSYM)
          NKTOT = NGSSHT(KSYM)
          NLTOT = NGSSHT(LSYM)
          IF(NTEST.GE.3) THEN
            WRITE(6,*) ' NITOT is  ',NITOT
            WRITE(6,*) ' NJTOT is  ',NJTOT
            WRITE(6,*) ' NKTOT is  ',NKTOT
            WRITE(6,*) ' NLTOT is  ',NLTOT
          END IF
*. Corresponding double group symmetries
          ISYM_U = IMOSF_SP(ISYM,1)
          ISYM_B = IMOSF_SP(ISYM,2)
          JSYM_U = IMOSF_SP(JSYM,1)
          JSYM_B = IMOSF_SP(JSYM,2)
          KSYM_U = IMOSF_SP(KSYM,1)
          KSYM_B = IMOSF_SP(KSYM,2)
          LSYM_U = IMOSF_SP(LSYM,1)
          LSYM_B = IMOSF_SP(LSYM,2)
*
          IF(NTEST.GE.3) THEN
            WRITE(6,*) ' ISYM,ISYM_B,ISYM_U', ISYM,ISYM_B,ISYM_U
            WRITE(6,*) ' JSYM,JSYM_B,JSYM_U', JSYM,JSYM_B,JSYM_U
            WRITE(6,*) ' KSYM,KSYM_B,KSYM_U', KSYM,KSYM_B,KSYM_U
            WRITE(6,*) ' LSYM,LSYM_B,LSYM_U', LSYM,LSYM_B,LSYM_U
          END IF
*
          NGSOBCT = 0
*
*
          DO LGAS = 1, NGAS
           NL = NGSSH(LSYM,LGAS)
           DO LLGAS = 1, LGAS
             IF(LLGAS.EQ.1) THEN
               LOFF = 1
             ELSE
               LOFF = LOFF + NGSSH(LSYM,LLGAS)
             END IF
           END DO
*
*
           DO KGAS = 1, NGAS
            NK = NGSSH(KSYM,KGAS)
            DO KKGAS = 1, KGAS
              IF(KKGAS.EQ.1) THEN
                KOFF = 1
              ELSE
                KOFF = KOFF + NGSSH(KSYM,KKGAS)
              END IF
            END DO
*
*
            DO JGAS = 1, NGAS
             NJ = NGSSH(JSYM,JGAS)
             DO JJGAS = 1, JGAS
               IF(JJGAS.EQ.1) THEN
                JOFF = 1
              ELSE
                JOFF = JOFF + NGSSH(JSYM,JJGAS)
              END IF
             END DO
*
*
             DO IGAS = 1, NGAS
              NI = NGSSH(ISYM,IGAS)
              DO IIGAS = 1, IGAS
                IF(IIGAS.EQ.1) THEN
                  IOFF = 1
                ELSE
                  IOFF = IOFF + NGSSH(ISYM,IIGAS)
                END IF
              END DO
*
*
              IGSOBCT = 0
              JGSOBCT = 0
              KGSOBCT = 0
              LGSOBCT = 0
*
              do IIGAS = 1,IGAS-1,1
                 IGSOBCT = IGSOBCT + NGSSH(ISYM,IIGAS)
              end do
*
              do JJGAS = 1,JGAS-1,1
                 JGSOBCT = JGSOBCT + NGSSH(JSYM,JJGAS)
              end do
              JGSOBCT = JGSOBCT * NITOT
*
              do KKGAS = 1,KGAS-1,1
                 KGSOBCT = KGSOBCT + NGSSH(KSYM,KKGAS)
              end do
              KGSOBCT = KGSOBCT * NITOT * NJTOT
*
              do LLGAS = 1,LGAS-1,1
                 LGSOBCT = LGSOBCT + NGSSH(LSYM,LLGAS)
              end do
              LGSOBCT = LGSOBCT * NITOT * NJTOT * NKTOT
*
              NGSOBCT = IGSOBCT + JGSOBCT + KGSOBCT + LGSOBCT
*
*
              if (NI.eq.0.or.NJ.eq.0.or.NK.eq.0.or.NL.eq.0)
     &           goto 666
*
              if (NTEST.ge.3) then
                write(6,'(A23,4I4)') ' IGAS, JGAS, KGAS, LGAS ',
     &          IGAS, JGAS, KGAS, LGAS
                write(6,'(A12,4I5)') ' NI,NJ,NK,NL', NI,NJ,NK,NL
c               write(6,*) 'IGSOBCT,JGSOBCT,KGSOBCT,LGSOBCT ='
c    &                     ,IGSOBCT,JGSOBCT,KGSOBCT,LGSOBCT
c               write(6,*) 'NGSOBCT =  ',NGSOBCT
              end if
*
*. Fetch nonrelativistic or scalar relativistic integrals
*  They are spin-free, that's for sure!
              if (NTEST.ge.2) then
                 write(6,*) 'Input for GETINTN:'
                 write(6,*) 'ISYM =  ',ISYM
                 write(6,*) 'JSYM =  ',JSYM
                 write(6,*) 'KSYM =  ',KSYM
                 write(6,*) 'LSYM =  ',LSYM
              end if
*
              CALL GETINTN(WORK(KLINT2),
     &           IGAS,ISYM,JGAS,JSYM,KGAS,KSYM,LGAS,LSYM,
     &           0,0,0,0,IRI)
C                  GETINTN(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &                    IXCHNG,IKSM,JLSM,ICOUL,IRI)
*
*(IJ!KL) is in WORK(KLINT2) as full (IK!JL) matrix
*
              if (NTEST.ge.2) then
                 write(6,*) 'What has GETINTN fetched?'
                 call wrtmat(WORK(KLINT2),NI*NK*NJ*NL,1,NI*NK*NJ*NL,1)
              end if
CTF
*. Distribute to integral lists
              DO JLIST = 1, NLIST
*. Spin integration : index 1,2 and index 3,4 should be identical
*    (well, actaully rather the barred or unbarred type should be
*     identical ...)
*     This is implemented for later generalization
CTF
                IF(IBU(1,JLIST).EQ.IBU(2,JLIST).AND.
     &             IBU(3,JLIST).EQ.IBU(4,JLIST)) THEN
CTF
* This means that the third list is accessed all the time, as it is
* the only one where all indices have the same IBU value (2), meaning
* all four indices are unbarred, corresponding to NR integrals.
CTF
*. nonvanishing spin integration
*. Offset to given symmetryblock in JLIST
                IF(IBU(1,JLIST).EQ.1) THEN
                 ISYM_D = ISYM_B
                ELSE
                 ISYM_D = ISYM_U
                END IF
*
                IF(IBU(2,JLIST).EQ.1) THEN
                 JSYM_D = JSYM_B
                ELSE
                 JSYM_D = JSYM_U
                END IF
*
                IF(IBU(3,JLIST).EQ.1) THEN
                 KSYM_D = KSYM_B
                ELSE
                 KSYM_D = KSYM_U
                END IF
*
                IF(IBU(4,JLIST).EQ.1) THEN
                 LSYM_D = LSYM_B
                ELSE
                 LSYM_D = LSYM_U
                END IF
*
                if (NTEST.ge.3) then
                   write(6,*) 'ISYM_D =  ',ISYM_D
                   write(6,*) 'JSYM_D =  ',JSYM_D
                   write(6,*) 'KSYM_D =  ',KSYM_D
                   write(6,*) 'LSYM_D =  ',LSYM_D
                end if
*
                IOFF_DBG = IOFF_DBG_LIST(
     &                     JLIST,ISYM_D,JSYM_D,KSYM_D,LSYM_D)
                if (NTEST.ge.3) write(6,*)
     &          'IOFF_DBG is  ',IOFF_DBG
*. And now : Distribute to double group list
                DO L = LOFF, LOFF+NL-1
                 DO K = KOFF, KOFF+NK-1
                  DO J = JOFF, JOFF+NJ-1
                   DO I = IOFF, IOFF+NI-1
*
                    IREL = I - IOFF + 1
                    JREL = J - JOFF + 1
                    KREL = K - KOFF + 1
                    LREL = L - LOFF + 1
*
*. Address in input list (has form (IK!JL))
                    IJKL_IN =
     &               ((LREL-1)*NJ+JREL-1)*NI*NK+(KREL-1)*NI+IREL
*
*. Address in output list (has usual structure (IJ!KL))
                    IJKL_OUT = IOFF_DBG - 1
     &                + ((LREL-1)*NKTOT+KREL-1)*NITOT*NJTOT
     &                + (JREL-1)*NITOT+IREL
     &                + NGSOBCT
*. and now
                    WORK(KINT2D(1,JLIST)-1+IJKL_OUT) =
     &              WORK(KLINT2-1+IJKL_IN)
                    if (NTEST.ge.5) then
                      write(6,*) 'IJKL_IN = ',IJKL_IN
                      write(6,*) 'WORK(KLINT2-1+IJKL_IN) ='
                      write(6,*) '   ',WORK(KLINT2-1+IJKL_IN)
                      write(6,*) 'IJKL_OUT = ',IJKL_OUT
                      write(6,*) 'WORK(KINT2D(1,JLIST)-1+IJKL_OUT) ='
                      write(6,*) '   ',WORK(KINT2D(1,JLIST)-1+IJKL_OUT)
                    end if
                   END DO
                  END DO
                 END DO
                END DO
*               ^ End of loop over orbital of given TS
               END IF
*              ^ End of check for nonvanishing spin integration
              END DO
*             ^End of loop over output lists
 666          continue
             END DO
            END DO
           END DO
          END DO
*          ^ End of loop over gasspaces
        END DO
       END DO
      END DO
*     ^ End of loop over symmetries
*( Oh Yes, I do, I do I do)
*. Vanishing integral lists
      DO JLIST = 1, NLIST
       LENGTH = L2D(JLIST)
       ZERO = 0.0D0
* 1 : all imaginary terms
       IOFF = KINT2D(2,JLIST)
C?     WRITE(6,*) ' JLIST, KINT2D(2,JLIST)',JLIST,  KINT2D(2,JLIST)
       CALL SETVEC(WORK(IOFF),ZERO,LENGTH)
* 2 : zero by spin integration
       IF(.NOT.(IBU(1,JLIST).EQ.IBU(2,JLIST).AND.
     &          IBU(3,JLIST).EQ.IBU(4,JLIST))     ) THEN
C?     WRITE(6,*) 'JLIST,  KINT2D(1,JLIST)',JLIST,  KINT2D(1,JLIST)
         IOFF = KINT2D(1,JLIST)
         CALL SETVEC(WORK(IOFF),ZERO,LENGTH)
       END IF
      END DO
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*)
        WRITE(6,*) ' Integral lists over spinors(real part) '
        WRITE(6,*) ' ---------------------------------------'
        WRITE(6,*) ' Distributed 2-electron integrals       '
        WRITE(6,*)
        DO ILIST = 1, NLIST
          WRITE(6,*)
          LEN_T = L2D(ILIST)
          IOFF = KINT2D(1,ILIST)
          write(6,*) 'List and length :'
          WRITE(6,'(I9)') ILIST,LEN_T
          WRITE(6,*)
          CALL WRTMAT(WORK(IOFF),LEN_T,1,LEN_T,1)
        END DO
      END IF
*
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'DIST2E')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Driver routine for distribution of DIRAC integrals to
* double group GAS lists
*  Creates some generally required information
*  Loops over spinor lists and real/imaginary parts
*
*----------------------------------------------
*  Timo Fleig, Oct. 2000
*----------------------------------------------
*
      subroutine dist_ints_dirac(HELPMAT,ONEINTR,ONEINTI,
     &                           SPININTS,IBOSYM,ISREO,
     &                           INVREO,
     &                           NL2D,IUB2D,ISM2D,
     &                           ITRALINE,ITOTDIM,ITOTDIM2,IHAM12,
     &                           ISPINFREE,IMXDIM,IDBG)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "symm.inc"
#include "onel.inc"
#include "input.inc"
*
      dimension ISREO(ITOTDIM2),INVREO(ITOTDIM),IBOSYM(ITOTDIM2)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
* What type of double group?
*  Real				 	IDGTYP = 1
*  Complex 		matrix group    IDGTYP = 2
*  Quaternionic				IDGTYP = 3
CTF  Only complex and real (inefficient) implemented for now.
*
      if (IDBG.le.3) IDGTYP = 1
      if (IDBG.ge.4.and.IDBG.le.6) IDGTYP = 2
      if (IDBG.ge.7.and.IDBG.le.8) IDGTYP = 3
*
* Check for inversion symmetry of group
      INVSYM = 0
      if (IDBG.eq.7.or.IDBG.eq.4.or.IDBG.eq.1) INVSYM = 1
*
* Check for matching dimensions of boson symmetry functions
      if (ISPINFREE.eq.1) then
        do I=1,NIRR_PN,1
          IBCT = 0
          IDGRP = 1
          IDOFF = 0
          if (I.ge.3) then 
            IDGRP = 2
            IDOFF = SYMDIM(1)
          end if
          do J=1+IDOFF,SYMDIM(IDGRP)+IDOFF,1
            if ((IBOSYM(J)+1).eq.I) IBCT = IBCT + 1
          end do
          if (IBCT.ne.SYDI(I)) then
            write(6,*)
            write(6,*) 'Problem with IBSPI array.'
            write(6,*) 'Boson symmetry dimensions do not match'
            write(6,*) 'with GAS matrix input.'
            write(6,*) 'Symmetry ',I
            write(6,*) 'No. of input functions ',SYDI(I)
            write(6,*) 'No. of functions on IBSPI array ',IBCT
            call abend('Quitting.')
          end if
        end do
        write(6,*)
        write(6,*) 'Boson symm. dimension check passed.'
      end if
*
*  Create symmetry reordering array for this list
*  ... and the inverted array
      ICT = 1
      JCT = 1
      IOFFBS = 0
      do INSYM=1,NIRR_PN,1
        if (INSYM.le.2) then
          IDIMGRP = 1
        else
          IDIMGRP = 2
        end if
        if (IDGTYP.eq.3) IDIMGRP = INSYM
*
* ... continue with offsets and skipping of even spinor
*     symmetries (these have redundant, i.e. Kramers-paired
*     classification
* For the moment, we hardwire the irrep order of C2h to the
*  generators XY and XYZ (necessarily in this order !!).
*  I.e., the irreps are classified 
*  as Ag, Bg, Au, Bu. This should be made more flexible later on.
        ILEN = SYMDIM(IDIMGRP)
        if (INSYM.ge.3) IOFFBS = 2*SYMDIM(1)
        if (IDBG.eq.4) then
*  INSYM_SEA corresponds to the irrep ordering generated by
*  the choice of symmetry generators in the integral program.
          if (INSYM.eq.1) INSYM_SEA = 1
          if (INSYM.eq.2) INSYM_SEA = 2
          if (INSYM.eq.3) INSYM_SEA = 3
          if (INSYM.eq.4) INSYM_SEA = 4
        else
          INSYM_SEA = INSYM
        end if
        if (INSYM.eq.3) ICT = 1
        if (NTEST.ge.10) then
          write(6,*) 'Current run       : ',INSYM
          write(6,*) 'Seeking symmetry  : ',INSYM_SEA
          write(6,*) 'Doing IND from ',1+IOFFBS, 'to',IOFFBS+ILEN
        end if
        do IND=1+IOFFBS,IOFFBS+ILEN,1
          if ((IBOSYM(IND)+1).eq.INSYM_SEA) then
            ISREO(IND) = ICT
            IND_REL = IND
            if (INSYM.ge.3) IND_REL = IND - SYMDIM(1)
            INVREO(JCT) = IND_REL
            JCT = JCT + 1
            ICT = ICT + 1
          end if
        end do
      end do
*  ISREO(orbital number DIRAC) = boson symmetry ordered orbital number
*        for one-electron part
*  INVREO(boson symmetry ordered orbital number) = orbital number DIRAC
*        for two-electron part
*
*
      if (NTEST.ge.3) then
        write(6,*) 'Symmetry reordering array DIRAC:'
        do IGRP=1,NIRR_PN/2,1
          if (IGRP.eq.1) then
            write(6,*) 'First dimension group:'
            call iwrtma(ISREO(1),1,SYMDIM(IGRP),1,SYMDIM(IGRP))
          else if (IGRP.eq.2) then
            IFEROFF = 2*SYMDIM(1)
            write(6,*) 'Second dimension group:'
            call iwrtma(ISREO(1+IFEROFF),1,SYMDIM(IGRP),1,SYMDIM(IGRP))
            write(6,*) 'Ferm. index offset for second IDIMGRP: ',IFEROFF
          end if
        end do
        write(6,*) 'Inverse reordering array:'
        call iwrtma(INVREO,1,ITOTDIM,1,ITOTDIM)
      end if
*
* Loop over spinor dg lists
      do LIST=1,NIRR_PN,1
        if (LIST.le.2) then
          IDIMGRP = 1
          IFEROFF = 0
        else
          IDIMGRP = 2
          IFEROFF = 2*SYMDIM(1)
        end if
*
        if (IDGTYP.eq.3) then
          IDIMGRP = LIST
        end if
*
*  Determine offsets for input indices
        INPOFF = 0
        do ILS=1,LIST-1,1
          INPOFF = INPOFF + SYMDIM(ILS)
        end do
*
        if (NTEST.ge.1) then
          write(6,*) 'Doing DG list    ',LIST
          write(6,*) 'Dimension group  ',IDIMGRP
          write(6,*) 'SYMDIM(IDIMGRP)  ',SYMDIM(IDIMGRP)
          write(6,*) 'Input index offset ',INPOFF
          write(6,*) 'Offset for reordering info array ',IFEROFF
        end if
*
* Number of 2-el. integrals on symmetry reduced list
*  and number of symmetry blocks
c       if (LIST.eq.1) call numints_red_rel(IGLSUM,NSMBLK,NTEST)
*
* Loop over real/imaginary part
        IRILP = 1
        if (ISPINFREE.eq.0) IRILP = 2
        do IRI=1,IRILP,1
          if (IRI.eq.1) then
            if (IHAM12.le.2) then
              call dist_1eints_dirac(HELPMAT,ONEINTR,IRI,INPOFF,
     &                               SPININTS,ISREO(1+IFEROFF),
     &                               ITRALINE,ITOTDIM2,SYMDIM(IDIMGRP),
     &                               LIST,NTEST)
            end if
            if (IHAM12.ge.2.and.LIST.eq.1) then
              call dist_2eints_rel(NL2D,IUB2D,ISM2D,IRI,NTEST)
            end if
          else if (IRI.eq.2) then
            if (IHAM12.le.2) then
              call dist_1eints_dirac(HELPMAT,ONEINTI,IRI,INPOFF,
     &                               SPININTS,ISREO(1+IFEROFF),
     &                               ITRALINE,ITOTDIM2,SYMDIM(IDIMGRP),
     &                               LIST,NTEST)
            end if
            if (IHAM12.ge.2.and.LIST.eq.1) then
              call dist_2eints_rel(NL2D,IUB2D,ISM2D,IRI,NTEST)
            end if
          end if
*
        end do
*       ^ End loop over real/imaginary part
      end do
*     ^ End loop over double group lists
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine dist_2eints_rel(NLIST,IUB,ISM,IRI,IPRTRA)
*
* Distribute relativistic integrals into
* double group format
*
* Modified good old non-relativistic routine
* dist_2eints
*
*---------------------------------------------
* Timo Fleig, December 2000
*---------------------------------------------
*
* Input :
*
* NLIST : Number of two electron integral lists
*
* IUB :
* IUB(INDEX,ILIST) = 1 => Index INDEX  in list ILIST is unbarred
*                  = 2 => Index INDEX  in list ILIST is   barred
*
* ISM :
* ISM(1,ILIST) = 1 => Symmetry between index 1 and 2 in list ILIST
*              = 0 => No Symmetry between index 1 and 2 in list ILIST
*
* ISM(2,ILIST) = 1 => Symmetry between index 3 and 4 nd  in list ILIST
*              = 0 => No Symmetry between index 3 and 4  in list ILIST
*
* ISM(3,ILIST) = 1 => Symmetry between index (12) and (34) in list ILIST
*              = 0 => No Symmetry between  index (12) and (34) in  ILIST
*
* Note : at the moment all symmetry flags are supposed to be zero !!
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "cintfo.inc"
#include "glbbas.inc"
*
      character*1 YINDEX
      dimension YINDEX(4)
      INTEGER IUB(4,*),ISM(3,*)
*
      NTESTL = 2
      NTEST = max(IPRTRA,NTESTL)
*
      if (NTEST.ge.1) then
        write(6,*)
        write(6,*) ' ================================ '
        write(6,*) '  Welcome to DIST_2EINTS_DIRAC    '
        write(6,*) ' ================================ '
        write(6,'(A,I3,A)') '    processing ',NLIST,' LISTS'
      end if
*
      IDUM = 1
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'DIST2E')
*
* Some additional information:
*. Block holding largest block of integrals
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
        DO IOBSM = 1, NSMOB
          MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
        END DO
      END DO
      INTBLK = MXTSOB_L**4
      call memmar(KLINT2,INTBLK,'ADDL  ',2,'KLINT2')
*
*
*. Distribute to double group integral lists
*  Loop over non-redundant lists
      do JLIST=1,NLIST,1
        if (NTEST.ge.2) then
          write(6,*) '+++++++++++++++++++++'
          write(6,'(A,I3)') '  This is list  ',JLIST
          write(6,*) '+++++++++++++++++++++'
        end if
*
* Barred/unbarred type of index (1 is unbarred)
        I_UB = IUB(1,JLIST)
        J_UB = IUB(2,JLIST)
        K_UB = IUB(3,JLIST)
        L_UB = IUB(4,JLIST)
*
*. Loop over double group symmetries 
*  determine corresponding point group symmetry
*  and fetch relevant integral block.
        do LFSYM=NIRR_PN+1,NSMOB,1
         do KFSYM=NIRR_PN+1,NSMOB,1
          do JFSYM=NIRR_PN+1,NSMOB,1
*
*. Double group symmetry
            LKFSYM = IDBGMULT(KFSYM,LFSYM)
            JIFSYM = IDBGMULT(1,IADJSYM(LKFSYM))
            IFSYM = IDBGMULT(JIFSYM,IADJSYM(JFSYM))
*
*. IFSYM, JFSYM, KFSYM and LFSYM are now known
            if (NTEST.ge.3) then
              write(6,*)
              write(6,*)
              write(6,'(A)') 'Double group index symmetries:'
              write(6,'(A,4I3)') ' IFSYM, JFSYM, KFSYM, LFSYM',
     &        IFSYM,JFSYM,KFSYM,LFSYM
            end if
*
* Determine boson symmetries of indices
            IBSYM = ISPSF_MO(IFSYM,I_UB)
            JBSYM = ISPSF_MO(JFSYM,J_UB)
            KBSYM = ISPSF_MO(KFSYM,K_UB)
            LBSYM = ISPSF_MO(LFSYM,L_UB)
*
            if (NTEST.ge.3) then
              write(6,'(A)') 'Corresponding point group symmetries:'
              write(6,'(A,4I3)') ' IBSYM, JBSYM, KBSYM, LBSYM',
     &        IBSYM,JBSYM,KBSYM,LBSYM
            end if
*
* Number of orbitals in this symmetry
            NITOT = NGSSHT(IBSYM)
            NJTOT = NGSSHT(JBSYM)
            NKTOT = NGSSHT(KBSYM)
            NLTOT = NGSSHT(LBSYM)
            if (NTEST.ge.3) then
              write(6,'(A)') 'Index dimensions:'
              write(6,'(A,4I3)') ' NITOT, NJTOT, NKTOT, NLTOT',
     &        NITOT,NJTOT,NKTOT,NLTOT
            end if

            NGSOBCT = 0
*
            DO LGAS = 1, NGAS
             NL = NGSSH(LBSYM,LGAS)
             DO LLGAS = 1, LGAS
               IF(LLGAS.EQ.1) THEN
                 LOFF = 1
               ELSE
                 LOFF = LOFF + NGSSH(LBSYM,LLGAS)
               END IF
             END DO
*
*
             DO KGAS = 1, NGAS
              NK = NGSSH(KBSYM,KGAS)
              DO KKGAS = 1, KGAS
                IF(KKGAS.EQ.1) THEN
                  KOFF = 1
                ELSE
                  KOFF = KOFF + NGSSH(KBSYM,KKGAS)
                END IF
              END DO
*
*
              DO JGAS = 1, NGAS
               NJ = NGSSH(JBSYM,JGAS)
               DO JJGAS = 1, JGAS
                 IF(JJGAS.EQ.1) THEN
                  JOFF = 1
                ELSE
                  JOFF = JOFF + NGSSH(JBSYM,JJGAS)
                END IF
               END DO
*
*
               DO IGAS = 1, NGAS
                NI = NGSSH(IBSYM,IGAS)
                DO IIGAS = 1, IGAS
                  IF(IIGAS.EQ.1) THEN
                    IOFF = 1
                  ELSE
                    IOFF = IOFF + NGSSH(IBSYM,IIGAS)
                  END IF
                END DO
*
*
                IGSOBCT = 0
                JGSOBCT = 0
                KGSOBCT = 0
                LGSOBCT = 0
*
                do IIGAS = 1,IGAS-1,1
                   IGSOBCT = IGSOBCT + NGSSH(IBSYM,IIGAS)
                end do
*
                do JJGAS = 1,JGAS-1,1
                   JGSOBCT = JGSOBCT + NGSSH(JBSYM,JJGAS)
                end do
                JGSOBCT = JGSOBCT * NITOT
*
                do KKGAS = 1,KGAS-1,1
                   KGSOBCT = KGSOBCT + NGSSH(KBSYM,KKGAS)
                end do
                KGSOBCT = KGSOBCT * NITOT * NJTOT
*
                do LLGAS = 1,LGAS-1,1
                   LGSOBCT = LGSOBCT + NGSSH(LBSYM,LLGAS)
                end do
                LGSOBCT = LGSOBCT * NITOT * NJTOT * NKTOT
*
                NGSOBCT = IGSOBCT + JGSOBCT + KGSOBCT + LGSOBCT
*
*
                if (NI.eq.0.or.NJ.eq.0.or.NK.eq.0.or.NL.eq.0) goto 666
*
                if (NTEST.ge.3) then
                  write(6,*)
                  write(6,'(A,4I4)') ' IGAS, JGAS, KGAS, LGAS ',
     &            IGAS, JGAS, KGAS, LGAS
                  write(6,'(A,4I5)') ' NI,NJ,NK,NL', NI,NJ,NK,NL
                  if (NTEST.ge.5) then
                    write(6,'(A,4I4)') 
     &                    'IGSOBCT,JGSOBCT,KGSOBCT,LGSOBCT ='
     &                    ,IGSOBCT,JGSOBCT,KGSOBCT,LGSOBCT
                    write(6,'(A,I5)') 'NGSOBCT =  ',NGSOBCT
                  end if
                end if
*
*. Fetch relativistic integrals from full real/imaginary lists
*  all supported environments (MOLCAS, DIRAC, ...) implemented.
                if (NTEST.ge.2) then
                   write(6,*) 'Input for GETINTN:'
                   write(6,*) 'IBSYM =  ',IBSYM
                   write(6,*) 'JBSYM =  ',JBSYM
                   write(6,*) 'KBSYM =  ',KBSYM
                   write(6,*) 'LBSYM =  ',LBSYM
                end if
*
                call getintn(WORK(KLINT2),
     &             IGAS,IBSYM,JGAS,JBSYM,KGAS,KBSYM,LGAS,LBSYM,
     &             0,0,0,0,IRI)
C                    getintn(XINT,
C                           ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &                      IXCHNG,IKSM,JLSM,ICOUL,IRI)
*
*(IJ!KL) is in WORK(KLINT2) as full (IK!JL) matrix
*
                if (NTEST.ge.2) then
                  write(6,*) 'What has GETINTN fetched?'
                  call wrtmat(WORK(KLINT2),NI*NK*NJ*NL,1,NI*NK*NJ*NL,1)
                end if
*
* Offset for double group output block
                IOFF_DBG = IOFF_DBG_LIST(JLIST,IFSYM,JFSYM,KFSYM,LFSYM)
                if (NTEST.ge.3) 
     &              write(6,'(A,I6)') 'IOFF_DBG is  ',IOFF_DBG
*
*. And now : Distribute to double group list
                do L=LOFF,LOFF+NL-1
                  do K=KOFF,KOFF+NK-1
                    do J=JOFF,JOFF+NJ-1
                      do I=IOFF,IOFF+NI-1
*
                        IREL = I - IOFF + 1
                        JREL = J - JOFF + 1
                        KREL = K - KOFF + 1
                        LREL = L - LOFF + 1
*
*. Address in input list (has form (IK!JL))
                        IJKL_IN = ((LREL-1)*NJ+JREL-1)*NI*NK
     &                           + (KREL-1)*NI
     &                           + IREL
*
*. Address in output list (has usual structure (IJ!KL))
                        IJKL_OUT = IOFF_DBG - 1
     &                    + ((LREL-1)*NKTOT+KREL-1)*NITOT*NJTOT
     &                    + (JREL-1)*NITOT+IREL
     &                    + NGSOBCT
*
*. and now process integral
                        WORK(KINT2D(IRI,JLIST)-1+IJKL_OUT) =
     &                  WORK(KLINT2-1+IJKL_IN)
                        if (NTEST.ge.5) then
                          write(6,'(A,I6)') 'IJKL_IN = ',IJKL_IN
                          write(6,'(A,1F18.12)') 
     &                          'WORK(KLINT2-1+IJKL_IN) =',
     &                           WORK(KLINT2-1+IJKL_IN)
                          write(6,'(A,I6)') 'IJKL_OUT = ',IJKL_OUT
                          write(6,'(A,1F18.12)') 
     &                          'WORK(KINT2D(1,JLIST)-1+IJKL_OUT) =',
     &                           WORK(KINT2D(1,JLIST)-1+IJKL_OUT)
                        end if
                      end do
                    end do
                  end do
                end do
*               ^ End of loop over orbital of given TS
 666            continue
               end do
              end do
             end do
            end do
*           ^ End of loop over gasspaces
          end do
         end do
        end do
*       ^ End of loop over symmetries
      end do
*     ^End of loop over output lists
*
      if (NTEST.ge.1) then
        WRITE(6,*)
        WRITE(6,*) ' Integral lists over spinors            '
        if (IRI.eq.1) write(6,*) '  (real part) '
        if (IRI.eq.2) write(6,*) '  (imaginary part) '
        WRITE(6,*) ' ---------------------------------------'
        WRITE(6,*) ' Distributed 2-electron integrals       '
        WRITE(6,*) ' ---------------------------------------'
        WRITE(6,*)
        do ILIST =1,NLIST,1
          WRITE(6,*)
          LEN_T = L2D(ILIST)
          IOFF = KINT2D(IRI,ILIST)
          do IL=1,4,1
            if (IUB(IL,ILIST).eq.1) YINDEX(IL) = 'U'
            if (IUB(IL,ILIST).eq.2) YINDEX(IL) = 'B'
          end do
          write(6,*) 'List and length :'
          write(6,'(I2,3X,A,4A1,A,3X,I6)') 
     &          ILIST,
     &          '(',(YINDEX(M),M=1,4,1),')',
     &          LEN_T
          write(6,*) ' ----------- '
          call wrtmat(WORK(IOFF),LEN_T,1,LEN_T,1)
        end do
      end if
*
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'DIST2E')
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Driver routine for distribution of DIRAC integrals to
* double group GAS lists
*
*   Routine using excitation operator sequence for integral
*   storage definition (SIGDEN route)
*
*----------------------------------------------
*  Timo Fleig, Mar. 2001 (Ellen Tabea is here!)
*  Jeppe Olsen modifies it. (Apr. 2001)
*  In complete form, June 2001, Timo Fleig
*
*  Implemented real groups (DIRAC) treated like complex groups
*  with zero imaginary parts.
*  Timo Fleig, August 25, 2003
*
*  implemented correct symmetry handling of one-electron property 
*  operator elements wrt to time-reversal symmetry
*  Stefan Knecht - Nov 2008
*
*----------------------------------------------
*
*  Implementation of T_BUFF use for MPIs "shared" memory version
* 
*
      subroutine dist_ints_sigden(INTIMP,ISPINFREE,ISYM_T,ISOBEX,T,
     &                            T_BUFF,
     &                            ISTR_CA,ISTR_CB,ISTR_AA,ISTR_AB,
     &                            NASHT,F1,F2,MDIASV,IBOSYM,
     &                            LABEL,CIRUN,IPRNT)
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "ctcc.inc"
#include "symm.inc"
#include "cgas.inc"
#include "orbinp.inc"
#include "glbbas.inc"
#include "clunit.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) IDISP
#endif
#include "parluci.h"
#include "files_r.inc"
#include "integrals_off.inc"
*
#include "dgroup.h"
*
      character*6 CIRUN
      dimension NICA_STR(MXPNGAS),NICB_STR(MXPNGAS),
     &          NIAA_STR(MXPNGAS),NIAB_STR(MXPNGAS)
      dimension ISOBEX(NGAS,4,NSPOBEX_TP),
     &          LABEL(4,NSPOBEX_TP)
      dimension I_ADD(4),INDEX(4),IND_SAVE(4),
     &          MDIASV(4,*)           
*
* Group notation of the strings:
      dimension IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS),
     &          IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
      dimension T(NALLINT), T_BUFF(*)
*
* Space for the construction of SUBSTRINGS
* dimensioned outside as integers MX_ST_TSOSO_BLK_MX
      dimension ISTR_CA(*),ISTR_CB(*),ISTR_AA(*),ISTR_AB(*)
      INTEGER*8 ICOUNT_I8
      LOGICAL NOINT_FILE
*
      NTESTL = 0000
      NTEST = max(NTESTL,IPRNT)
!     NTEST = 0
*
* Constants
      IZERO = 0
      ZERO = 0.D0
!
!     control writing of file INDICES.x (x =0,... last slave tag)
      NOINT_FILE = .FALSE.
      if( lowsrt_ijkl .or. reord_ijkl ) NOINT_FILE = .TRUE.
C     NOINT_FILE = .TRUE.
      IF( NOINT_FILE )THEN
        write(LUWRT,'(/A)') ' NOTICE from DIST_INTS_SIGDEN: '
        write(LUWRT,'( A)') ' INDICES file is not written as'
        write(LUWRT,'(/A)') ' requested by ".IJKLRO".       '
        write(LUWRT,'( A)') ' The calculation of natural    '
        write(LUWRT,'( A)') ' orbital occupations numbers   '
        write(LUWRT,'( A)') ' is therefore not allowed.     '
      END IF
*
* Double group type
      if (NZ.eq.1) then
        IDGTYP = 1
        IF( NTEST .gt. 1001 )THEN
          write(LUWRT,'(/A)') ' NOTICE from DIST_INTS_SIGDEN: '
          write(LUWRT,'( A)') ' Running in real algebra.      '
          write(LUWRT,'(A/)') ' Imaginary loops eliminated.   '
        END IF
      else if (NZ.eq.2) then
        IDGTYP = 2
      else if (NZ.eq.4) then
        IDGTYP = 3
      end if
*
* Check for inversion symmetry of group
      INVSYM = 0
      if (DOUGRP.eq.7.or.DOUGRP.eq.4.or.DOUGRP.eq.1) INVSYM = 1
*
* Get reordering info depending on environment
      if (INTIMP.eq.6) then
        call reo_info_dirac(NDIM_TOT,NIRR_PN,SYMDIM,IBOSYM,
     &                      SYDI,DOUGRP,WORK(KISREO),WORK(KINVREO),
     &                      INVSYM,IDGTYP,ISPINFREE,NTEST)
*
        if (ISPINFREE.eq.1.or.IDGTYP.eq.1) then
          IRILP = 1
        else
          IRILP = 2
        end if
      else if (INTIMP.eq.7) then
        call reo_info_dirac(NDIM_TOT,NIRR_PN,SYMDIM,IBOSYM,
     &                      SYDI,DOUGRP,WORK(KISREO),WORK(KINVREO),
     &                      INVSYM,IDGTYP,ISPINFREE,NTEST)
        if (ISPINFREE.eq.1.or.IDGTYP.eq.1) then
          IRILP = 1
        else
          IRILP = 2
        end if
      end if
!
!     open ASCII file for checking integral sorting information
!     ntest = 1
      if (NTEST.ge.1) then
        if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
          LUINGR = LUSC94
          open(unit=LUINGR,file=NLUSC94_MPI(1:LUSC94_MPIL),
     &         status='UNKNOWN',form='FORMATTED')
          rewind(LUINGR)
        end if
      end if
!
!     open file for saving index sequence
      LUINSV = LUSC95
      open(unit=LUINSV,file=NLUSC95_MPI(1:LUSC95_MPIL),
     &     status='UNKNOWN',form='FORMATTED')
      rewind(LUINSV)
!
!     open file with pre-fetched integrals
      if( lowsrt_ijkl ) then
#if defined (VAR_MPI2)
        call mpi_file_open(mpi_comm_self,"IJKL_PRESORT",
     &                     mpi_mode_rdonly + mpi_mode_delete_on_close,
     &                     mpi_info_null,ijkl_pre,ierr)
        idisp = 0
        call mpi_file_set_view(ijkl_pre,idisp,mpi_real8,mpi_real8,
     &                         "native",mpi_info_null,ierr)
#endif
      end if
!
!     set counter to start value and initialize
      ICOUNT_I8 = 1
      NDIAINT   = 0
      N1DIAINT  = 0
      N2DIAINT  = 0
!
!     Loop over real/imaginary part
!
      do IRI=1,IRILP,1
        if (NTEST.ge.1) then
          if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
            if (IRI.eq.1.and.NTEST.ge.2) write(LUWRT,*) 'REAL LOOP'
            if (IRI.eq.1) write(LUINGR,'(A)') 'Real  integrals.'
            if (IRI.eq.2.and.NTEST.ge.2) write(LUWRT,*) 'IMAG LOOP'
            if (IRI.eq.2) write(LUINGR,'(A)') 'Imag. integrals.'
          end if
        end if
        ITPSV = 0
        if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
          if (NTEST.ge.1) 
     &      write(LUINGR,'(A)') 'INTEGRAL FILE'
        end if
*
* Loop over excitation types
        do IXTP=1,NSPOBEX_TP,1
*
* Check for correct pointers to type blocks
          if (IRI.eq.1.and.ICOUNT_I8.ne.IBTSOSO_TP(IXTP)) then
            write(LUWRT,*) 'Pointers to types do not match numbers'
            write(LUWRT,*) 'of integrals generated.'
            write(LUWRT,*) 'LABEL array: ',(LABEL(I,IXTP),I=1,4,1)
            write(LUWRT,*) 'int.no.,offset: ',ICOUNT_I8,IBTSOSO_TP(IXTP)
            write(LUWRT,*) 'Current type: ',IXTP
            call abend2('Quitting.')
          end if
          if (NTEST.ge.1) then
            write(LUWRT,*) '======================'
            write(LUWRT,*) 'Processing type ',IXTP
            write(LUWRT,*) '======================'
            write(LUWRT,*) 'ISOBEX partition:'
            call iwrtma(ISOBEX(1,1,IXTP),NGAS,4,NGAS,4)
          end if
*
* write statement about next integral type (to file)
          if (CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 ') then
            if (NTEST.ge.1) 
     &          call integral_type(ISOBEX(1,1,IXTP),NGAS,ITPSV,LUINGR)
            if (NTEST.ge.1) 
     &        write(LUINGR,'(A,I5)') 'Type ',IXTP
          end if
*
*. Change from Occupation to group notation of strings
          call occ_to_grp(ISOBEX(1,1,IXTP),IGRP_CA,1)
          call occ_to_grp(ISOBEX(1,2,IXTP),IGRP_CB,1)
          call occ_to_grp(ISOBEX(1,3,IXTP),IGRP_AA,1)
          call occ_to_grp(ISOBEX(1,4,IXTP),IGRP_AB,1)
          if (NTEST.ge.5.and.IRI.eq.1) then
            write(LUWRT,*) 'Group information:'
            write(LUWRT,'(A,4I6)') 'IGRP_CA : ',(IGRP_CA(I),I=1,NGAS,1)
            write(LUWRT,'(A,4I6)') 'IGRP_CB : ',(IGRP_CB(I),I=1,NGAS,1)
            write(LUWRT,'(A,4I6)') 'IGRP_AA : ',(IGRP_AA(I),I=1,NGAS,1)
            write(LUWRT,'(A,4I6)') 'IGRP_AB : ',(IGRP_AB(I),I=1,NGAS,1)
          end if
*
* Corresponding number of electrons
          NEL_CA = ielsum(ISOBEX(1,1,IXTP),NGAS)
          NEL_CB = ielsum(ISOBEX(1,2,IXTP),NGAS)
          NEL_AA = ielsum(ISOBEX(1,3,IXTP),NGAS)
          NEL_AB = ielsum(ISOBEX(1,4,IXTP),NGAS)
          if (NTEST.ge.2.and.IRI.eq.1) then
            write(LUWRT,*) 'NEL_CA = ',NEL_CA
            write(LUWRT,*) 'NEL_CB = ',NEL_CB
            write(LUWRT,*) 'NEL_AA = ',NEL_AA
            write(LUWRT,*) 'NEL_AB = ',NEL_AB
          end if
*
* Symmetry loops
* (Annihilator symm. is assumed in terms of corresponding creators!)
*  Symmetry restrictions according to number of electrons
          NEL_C = NEL_CA + NEL_CB
          NEL_A = NEL_AA + NEL_AB
          call evenodd(IEONELC,NEL_C)
          if (IEONELC.eq.1) then
            ISMCSTA = NIRR_PN + 1
            ISMCEND = NIRR_DG
          else if (IEONELC.eq.2) then
            ISMCSTA = 1
            ISMCEND = NIRR_PN
          end if
*
          call evenodd(IEONELCA,NEL_CA)
          if (IEONELCA.eq.1) then
            ISMCASTA = NIRR_PN + 1
            ISMCAEND = NIRR_DG
          else if (IEONELCA.eq.2) then
            ISMCASTA = 1
            ISMCAEND = NIRR_PN
          end if
*
          if (NTEST.ge.5.and.IRI.eq.1) then
            write(LUWRT,*) 'C  symmetry loop: ',ISMCSTA,ISMCEND
            write(LUWRT,*) 'CA symmetry loop: ',ISMCASTA,ISMCAEND
          end if
*
* Symmetry loops
          do ISM_C=ISMCSTA,ISMCEND,1
            ISM_A = IDBGMULT(ISM_C,IADJSYM(ISYM_T))
            do ISM_CA=ISMCASTA,ISMCAEND,1
              ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
              do ISM_AA=1,NIRR_DG,1
                ISM_AB = IDBGMULT(ISM_A,INVELM(ISM_AA))
                if (NTEST.ge.2000.and.IRI.eq.1) then
                  write(LUWRT,*)
                  write(LUWRT,*)
                  write(LUWRT,*) 'Printing symmetries:'
                  write(LUWRT,*) 'Creator symm.      :',ISM_C
                  write(LUWRT,*) 'Annihilator symm.  :',ISM_A
                  write(LUWRT,*) 'Creator alpha symm.:',ISM_CA
                  write(LUWRT,*) 'Creator beta symm. :',ISM_CB
                  write(LUWRT,*) 'Annihi. alpha symm.:',ISM_AA
                  write(LUWRT,*) 'Annihi. beta symm. :',ISM_AB
                end if
*
*. Generate strings of these symmetries
*. obtain strings
                IUB = 1
                call getstr2_totsm_spgp_rel(IUB,IGRP_CA,NGAS,ISM_CA,
     &                                      NEL_CA,NSTR_CA,ISTR_CA,
     &                                      NORBT,0,IDUM,IDUM)
                IUB = 2
                call getstr2_totsm_spgp_rel(IUB,IGRP_CB,NGAS,ISM_CB,
     &                                      NEL_CB,NSTR_CB,ISTR_CB,
     &                                      NORBT,0,IDUM,IDUM)
                IUB = 1
                call getstr2_totsm_spgp_rel(IUB,IGRP_AA,NGAS,ISM_AA,
     &                                      NEL_AA,NSTR_AA,ISTR_AA,
     &                                      NORBT,0,IDUM,IDUM)
                IUB = 2
                call getstr2_totsm_spgp_rel(IUB,IGRP_AB,NGAS,ISM_AB,
     &                                      NEL_AB,NSTR_AB,ISTR_AB,
     &                                      NORBT,0,IDUM,IDUM)
*
* Loops over strings for the given symmetries and types
* Loop over T elements as matrix T(I_CA,I_CB,I_AA,I_AB)
*  Determine possible connections
                if ((NSTR_CA.ge.1.or.NSTR_CB.ge.1).and.
     &              (NSTR_AA.ge.1.or.NSTR_AB.ge.1)) then
*
                  do I_AB=1,NSTR_AB,1
                    do I_AA=1,NSTR_AA,1
                      do I_CB=1,NSTR_CB,1
                        do I_CA=1,NSTR_CA,1
*
* Notes: 
*  1) Standard notation in second quantization for 2-el. ints is:
*     a+I a+K_ aL aJ_ (IJ|KL) = - a+I a+K_ aJ_ aL (IJ|KL)
*     The code proceeds according to the right hand side. The
*     signs are located in WORK(KSIGNNHX) for each excitation
*     class. So first indices 2 and 3 are interchanged, then
*     3 and 4 (implicit in the code).
*  2) Integral classes with 2 alpha or 2 beta creators (or both)
*     are treated differently: (due to exchange character)
*     a+I a+K aL aJ [(IJ|KL) - (IL|KJ)]
*  3) In the full spin-dependent formalism, extra classes of
*     integrals are required. 
* 
*. Indices of STRING I_CA are stored in IOCC_CA and first element is 
*  IOCC_CA(1+NEL_CA*(I_CA-1)) 
*
                          if (NEL_C.eq.1) then
                            IXCH = 0
*
                            if (NEL_CA.eq.1.and.NEL_CB.eq.0) then
                              INDEX(1) 
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                              INDEX(2) = 0
                            else if (NEL_CA.eq.0.and.NEL_CB.eq.1) then
                              INDEX(1) = 0
                              INDEX(2) 
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of creators.'
                              write(LUWRT,*) 
     &                              'NEL_CA,NEL_CB = ',NEL_CA,NEL_CB
                              call abend2('Quitting.')
                            end if
*
                            if (NEL_AA.eq.1.and.NEL_AB.eq.0) then
                              INDEX(3) 
     &                              = ISTR_AA(1+NEL_AA*(I_AA-1))
                              INDEX(4) = 0
                            else if (NEL_AA.eq.0.and.NEL_AB.eq.1) then
                              INDEX(3) = 0
                              INDEX(4) 
     &                              = ISTR_AB(1+NEL_AB*(I_AB-1))
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of annihil.'
                              write(LUWRT,*) 
     &                              'NEL_AA,NEL_AB = ',NEL_AA,NEL_AB
                              call abend2('Quitting.')
                            end if
*
                          else if (NEL_C.eq.2) then
*
                            if (NEL_CA.eq.2.and.NEL_CB.eq.0) then
                              IXCH = 0
                              if (NEL_AA.eq.2) then
*                               Exchange integrals
                                IXCH = 1
                                do IEL=1,NEL_CA,1
                                  INDEX(IEL) = 
     &                                ISTR_CA(IEL+NEL_CA*(I_CA-1))
                                end do
                              else if (NEL_AB.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(1)
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                                INDEX(2)
     &                              = ISTR_CA(2+NEL_CA*(I_CA-1))
                              else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
*                               Spin-orbit integral
                                IXCH = 0
                                do IEL=1,NEL_CA,1
                                  INDEX(IEL) =
     &                                ISTR_CA(IEL+NEL_CA*(I_CA-1))
                                end do
                              end if
                            else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
                              IXCH = 0
                              INDEX(1) 
     &                              = ISTR_CA(1+NEL_CA*(I_CA-1))
                              INDEX(3) 
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                            else if (NEL_CA.eq.0.and.NEL_CB.eq.2) then
                              IXCH = 0
                              if (NEL_AB.eq.2) then
*                               Exchange integrals
                                IXCH = 1
                                do IEL=1,NEL_CB,1
                                  INDEX(IEL) = 
     &                                ISTR_CB(IEL+NEL_CB*(I_CB-1))
                                end do
                              else if (NEL_AA.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(1)
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                                INDEX(2)
     &                              = ISTR_CB(2+NEL_CB*(I_CB-1))
                              else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
*                               Spin-orbit integral
                                INDEX(1)
     &                              = ISTR_CB(1+NEL_CB*(I_CB-1))
                                INDEX(3)
     &                              = ISTR_CB(2+NEL_CB*(I_CB-1))
                              end if
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of creators.'
                              write(LUWRT,*) 
     &                              'NEL_CA,NEL_CB = ',NEL_CA,NEL_CB
                              call abend2('Quitting.')
                            end if
*
                            if (NEL_AA.eq.2.and.NEL_AB.eq.0) then
                              if (NEL_CA.eq.2) then
                                IXCH = 1
                                do IEL=1,NEL_AA,1
                                 INDEX(IEL+2) = 
     &                                ISTR_AA(IEL+NEL_AA*(I_AA-1))
                                end do
                              else if (NEL_CB.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                IXCH = 1
                                INDEX(3)
     &                              = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(4)
     &                              = ISTR_AA(2+NEL_AA*(I_AA-1))
                              else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spin-orbit integral
                                do IEL=1,NEL_AA,1
                                 INDEX(IEL*2) =
     &                                ISTR_AA(IEL+NEL_AA*(I_AA-1))
                                end do
                              end if
                            else if (NEL_AA.eq.1.and.NEL_AB.eq.1) then
                              IXCH = 0
                              if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spinfree integral
                                INDEX(2) 
     &                                = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(4) 
     &                                = ISTR_AB(1+NEL_AB*(I_AB-1))
                              else
*                               Spin-orbit integral
                                INDEX(4)
     &                                = ISTR_AA(1+NEL_AA*(I_AA-1))
                                INDEX(2)
     &                                = ISTR_AB(1+NEL_AB*(I_AB-1))
                              end if
                            else if (NEL_AA.eq.0.and.NEL_AB.eq.2) then
                              if (NEL_CB.eq.2) then
                                IXCH = 1
                                do IEL=1,NEL_AB,1
                                 INDEX(IEL+2) = 
     &                                ISTR_AB(IEL+NEL_AB*(I_AB-1))
                                end do
                              else if (NEL_CA.eq.2) then
*                               Spin-orbit integrals (no exchange)
                                INDEX(3)
     &                              = ISTR_AB(1+NEL_AB*(I_AB-1))
                                INDEX(4)
     &                              = ISTR_AB(2+NEL_AB*(I_AB-1))
                              else if (NEL_CA.eq.1.and.NEL_CB.eq.1) then
*                               Spin-orbit integral
                                do IEL=1,NEL_AB,1
                                 INDEX(IEL*2) =
     &                                ISTR_AB(IEL+NEL_AB*(I_AB-1))
                                end do
                              end if
                            else
                              write(LUWRT,*) 
     &                              'Impossible number of annihil.'
                              write(LUWRT,*) 
     &                              'NEL_AA,NEL_AB = ',NEL_AA,NEL_AB
                              call abend2('Quitting.')
                            end if
*
                          else
                         write(LUWRT,*) 'Impossible number of creators.'
                         write(LUWRT,*) 'NEL_C = ',NEL_C
                         call abend2('Quitting.')
                          end if
*
*  Stored integrals (DIRAC) refer to symmetry-type ordering.
*  Internal ordering is type-symmetry, so rename:
                          IR_IND = 0
                          JR_IND = 0
                          KR_IND = 0
                          LR_IND = 0
                          if (INDEX(1).ne.IZERO) 
     &                              IR_IND = IREOTS(INDEX(1))
                          if (INDEX(2).ne.IZERO) 
     &                              JR_IND = IREOTS(INDEX(2))
                          if (INDEX(3).ne.IZERO) 
     &                              KR_IND = IREOTS(INDEX(3))
                          if (INDEX(4).ne.IZERO) 
     &                              LR_IND = IREOTS(INDEX(4))
                          INDEX(1) = IR_IND
                          INDEX(2) = JR_IND
                          INDEX(3) = KR_IND
                          INDEX(4) = LR_IND
*
* Index reordering and second integral in exchange case:
* (Reason: CACA string generated e.g. for two alpha excitations,
*          but this does not correspond to index ordering in integral!)
*  Apply only to 2-electron integrals:
*
                          YXINT = ZERO
                          YINT  = ZERO
                          if (IXCH.eq.1) then
* Index sequence reads:
*  (IJ|KL) (which is not the operator sequence for the desired int!)
* Reorder for Coulomb integral:
*  (IJ|KL) --> (IK|LJ)
* Now create exchange integral:
*  (IJ|KL) --> (IJ|LK)
* (Overall, indices L and J are interchanged and a "-" introduced,
*  due to I+ K+ L J = - I+ K+ J L
* Coulomb:
                            I_SV = INDEX(1)
                            L_SV = INDEX(2)
                            K_SV = INDEX(3)
                            J_SV = INDEX(4)
                            INDEX(1) = I_SV
                            INDEX(2) = K_SV
                            INDEX(3) = L_SV
                            INDEX(4) = J_SV
* Exchange:
                            I_ADD(1) = I_SV
                            I_ADD(2) = J_SV
                            I_ADD(3) = L_SV
                            I_ADD(4) = K_SV
                            do IRUN = 1,4,1
                              I_ADD(IRUN) = 
     &                              I_ADD(IRUN) * 
     &                              LABEL(IRUN,IXTP)
                            end do
                            if (IRI.eq.1.or.(IRI.eq.2.and.NZ.ne.1)) then
                              if (CIRUN.ne.'DENS1 '
     &                            .and.CIRUN.ne.'DENS2 ')
     &                          call getint_sigden(YXINT,I_ADD,IRI,
     &                                    NDIM_TOT,NASHT,F1,F2,NTEST)
                            end if
                          end if
*
                          if (CIRUN.ne.'DENS1 '.and.
     &                        CIRUN.ne.'DENS2 ') then
*
*  Modify if integrals are over relativistic spinors:
                            if ((INTIMP.eq.6.or.INTIMP.eq.7).and.
     &                           ISPINFREE.eq.0) then
                              if (INDEX(1).ne.0.and.
     &                            INDEX(2).ne.0) then
                                do IRUN = 1,4,1
                                    INDEX(IRUN) = 
     &                                  INDEX(IRUN) * 
     &                                  LABEL(IRUN,IXTP)
                                end do
*
* Correct internal ordering of two-bar exchange type integrals:
* (This might be implemented in some more elegant way later on.)
                                if (INDEX(1).gt.0.and.
     &                              INDEX(2).lt.0.and.
     &                              INDEX(3).lt.0.and.
     &                              INDEX(4).gt.0) then
                                  J_SV = INDEX(2)
                                  L_SV = INDEX(4)
                                  INDEX(2) = -L_SV
                                  INDEX(4) = -J_SV
                                end if
                              end if
                            end if
*
* Check for diagonal integrals
                            if (IRI.eq.1) then
                              ICNT = 0
                              do IR = 1,4,1
                                if (INDEX(IR).ne.0) then
                                  ICNT = ICNT + 1
                                  IND_SAVE(ICNT) = INDEX(IR)
                                end if
                              end do
                              if (ICNT.eq.2) then
                                if (IND_SAVE(1).eq.IND_SAVE(2)) then
                                  if (INDEX(1).ne.IZERO) IUBDIA = 1
                                  if (INDEX(2).ne.IZERO) IUBDIA = 2
                                  if (IUBDIA.eq.1) then
                                    NDIAINT = NDIAINT + 1
                                    N1DIAINT = N1DIAINT + 1
                                    IPOSDIA(NDIAINT) = ICOUNT_I8
                                    MDIASV(1,NDIAINT) = IND_SAVE(1)
                                    MDIASV(2,NDIAINT) = IZERO
                                    MDIASV(3,NDIAINT) = IND_SAVE(2)
                                    MDIASV(4,NDIAINT) = IZERO
C           print*,'D INT ',NDIAINT,' is ',(MDIASV(JR,NDIAINT),JR=1,4),
C    &             ' at ',ICOUNT_I8
                                  else if (IUBDIA.eq.2) then
                                    NDIAINT = NDIAINT + 1
                                    N1DIAINT = N1DIAINT + 1
                                    IPOSDIA(NDIAINT) = ICOUNT_I8
                                    MDIASV(1,NDIAINT) = IZERO
                                    MDIASV(2,NDIAINT) = IND_SAVE(1)
                                    MDIASV(3,NDIAINT) = IZERO
                                    MDIASV(4,NDIAINT) = IND_SAVE(2)
C           print*,'D INT ',NDIAINT,' is ',(MDIASV(JR,NDIAINT),JR=1,4),
C    &             ' at ',ICOUNT_I8
                                  end if
                                end if
                              else if (ICNT.eq.4) then
                                if ((INDEX(1).eq.INDEX(4)).and.
     &                              (INDEX(2).eq.INDEX(3)).or.
     &                              (INDEX(1).eq.INDEX(2)).and.
     &                              (INDEX(3).eq.INDEX(4)).or.
     &                              (INDEX(1).eq.-INDEX(3)).and.
     &                              (INDEX(2).eq.-INDEX(4))) then
                                  NDIAINT = NDIAINT + 1
                                  N2DIAINT = N2DIAINT + 1
                                  IPOSDIA(NDIAINT) = ICOUNT_I8
                                  do IR = 1,4,1
                                    MDIASV(IR,NDIAINT) = INDEX(IR)
                                  end do
C           print*,'D INT ',NDIAINT,' is ',(MDIASV(JR,NDIAINT),JR=1,4),
C    &             ' at ',ICOUNT_I8
                                end if
                              else
                                write(6,*) ' Diagonal check failed in '
                                write(6,*) ' dist_ints_sigden.        '
                                write(6,*) '  ICNT = ',ICNT
                                stop 'Quitting.'
                              end if
                            end if
*
*                           Write INDEX sequence to line in file
*                           for use in diagonal routine (saving memory):
                            IF(.NOT. NOINT_FILE )
     &                         write(LUINSV,'(4I4)') (INDEX(JJ),JJ=1,4)
*
                            if (IRI.eq.1.or.(IRI.eq.2.and.NZ.ne.1)) then
                              call getint_sigden(YINT,INDEX(1),
     &                                  IRI,NDIM_TOT,NASHT,F1,F2,NTEST)
                            end if
C
                            T(ICOUNT_I8) = YINT - YXINT
C FIXME !! Uncommented sign here!!
                            if (IXCH.eq.1) T(ICOUNT_I8) = - T(ICOUNT_I8)
C FIXME !! Uncommented sign here!!
                              if (NTEST.ge.1) then
                                 write(LUWRT,'(A,I4,A,4I3,A,1F18.14)') 
     &                              'INTEGRAL no. ',ICOUNT_I8,'(',
     &                              (INDEX(M),M=1,4),
     &                              ') Value: ',T(ICOUNT_I8)
                                 write(LUWRT,'(A,3F18.14)') 
     &                              'INTEGRALs: YINT,YXINT,YINT-YXINT',
     &                                          YINT,YXINT,YINT-YXINT
*                             write integral to file (for debugging):
                                if (IXCH.eq.0) then
                                  write(LUINGR,'(I4,4X,4I3,1F18.14)') 
     &                                  ICOUNT_I8,INDEX(1),
     &                                         INDEX(2),
     &                                         INDEX(3),
     &                                         INDEX(4),T(ICOUNT_I8)
                                else
                        write(LUINGR,'(I4,4X,4I3,1X,A,1X,4I3,1F18.14)')
     &                                  ICOUNT_I8,INDEX(1),
     &                                         INDEX(2),
     &                                         INDEX(3),
     &                                         INDEX(4),'-',
     &                              I_ADD(1),I_ADD(2),I_ADD(3),I_ADD(4),
     &                                         T(ICOUNT_I8)
                                end if
                              end if
                          end if
                          ICOUNT_I8 = ICOUNT_I8 + 1
*
                        end do
*                       ^ End loop over alpha creation strings
                      end do
*                     ^ End loop over beta creation strings
                    end do
*                   ^ End loop over alpha annihilation strings
                  end do
*                 ^ End loop over beta annihilation strings
                end if
*               ^ End if connection condition
              end do
*             ^ End loop over alpha annihilator symmetry
            end do
*           ^ End loop over alpha creator symmetry
          end do
*         ^ End loop over creation part symmetry
        end do
*       ^ End loop over excitation types
      end do
*     ^ End loop over real/imaginary
      if (NTEST.ge.1.and.(CIRUN.ne.'DENS1 '.and.CIRUN.ne.'DENS2 '))
     &   close(LUINGR)
      close(LUINSV)
      NALLINT = ICOUNT_I8 - 1
C
 300  CONTINUE
C
      if (NTEST.ge.1) then
        write(LUWRT,*) ' DIST_INTS_SIGDEN reporting   : '
        write(LUWRT,*) '  Total number of integrals   : ',NALLINT
        write(LUWRT,*) '  Number of diagonal integrals: ',NDIAINT
        write(LUWRT,*) '    Number of diagonal 1-el. integrals: ',
     &                      N1DIAINT
        write(LUWRT,*) '    Number of diagonal 2-el. integrals: ',
     &                      N2DIAINT
        write(LUWRT,*) '    sorted integrals: '
        CALL WRTMATMN8(T,1,NALLINT,1,NALLINT,LUWRT)
      end if
      NTEST = 0
!
!     close (and destroy) file with pre-fetched integrals
      if( lowsrt_ijkl ) then
#if defined (VAR_MPI2)
        call mpi_file_close(ijkl_pre,ierr)
#endif
      end if
*
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_ORB_DIM_ENV_REL(ECORE_ENV)
*
* Obtain number of orbitals and basis functions from the
* programming environment.
* results stored in NAOS_ENV, NMOS_ENV
*
* Obtain environments CORE energy, ECORE_ENV
*
* Jeppe Olsen, December 97
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "crun.inc"
#include "lucinp.inc"
*
      IF(ENVIRO(1:6).EQ.'MOLCAS') THEN
CC        CALL GETOBS2(ECORE_ENV,NAOS_ENV,NMOS_ENV)
C            Getobs2(ECORE,NAOS_ENV,NMOS_ENV)
      ELSE IF(ENVIRO(1:6).EQ.'DALTON' ) THEN
CC        CALL GETOBS_DALTON(ECORE_ENV,NAOS_ENV,NMOS_ENV)
      ELSE IF(ENVIRO(1:4).EQ.'NONE') THEN
*. No environment,
        WRITE(6,*) 'GET_ORB_DIM_ENV  in problems '
        WRITE(6,*) 'No ENVIRO parameter defined '
        IZERO = 0
        CALL ISETVC(NAOS_ENV,IZERO,NSMOB)
        CALL ISETVC(NMOS_ENV,IZERO,NSMOB)
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' From GET_ORB_FROM_ENV : '
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' NAOS_ENV'
        CALL IWRTMA(NAOS_ENV,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' NMOS_ENV'
        CALL IWRTMA(NMOS_ENV,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' ECORE_ENV=', ECORE_ENV
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  get desired block of one-electron integrals
*
*
*
*================================================
*================================================
*     CAUTION !!!!!
*
* At the moment, geth1rel is only designed for an
* initial CI calculation. As soon as the transformed
* integral matrices contain mixed off-diagonal and
* real/imaginary contributions, the summation has
* to be expanded to all lists !!!!
*================================================
*================================================
*
*
      subroutine geth1rel3(SPININTS,
     &                     H,ISM,ITYP,JSM,JTYP,IRISIG,IRICO,IPM1,
     &                     IJUB,IRETURN,SCALIN,INTSPACE,ENVIRO)
*
*-----------------------------------
*   Timo Fleig
*-----------------------------------
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
#include "onel.inc"
*
      character*72 ENVIRO
*
      dimension SPININTS(INTSPACE)
      dimension H(*)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      if (NTEST.ge.3) then
        write(6,*)
        write(6,*) 'Checking input in geth1rel3 '
        write(6,*) 'ISM       ',ISM
        write(6,*) 'JSM       ',JSM
        write(6,*) 'ITYP      ',JTYP
        write(6,*) 'IRISIG    ',IRISIG
        write(6,*) 'IRICO     ',IRICO
        write(6,*) 'IPM1      ',IPM1
        write(6,*) 'IJUB      ',IJUB
      end if
*
      ITRALINE = INTSPACE / 32
      ITRACHK = ITRALINE*32
      if (ITRACHK.ne.INTSPACE) then
         write(6,*) 'Problem with int space in geth1rel3.'
         Call Abend2( 'quitting' )
      end if
*
      IOFFTRA(1) = 1
      do I=2,32,1
         IOFFTRA(I) = IOFFTRA(I-1) + ITRALINE
      end do
*
      call findhigh(NGSOBPI,MXPIRR*MXPNGAS,MXGSDIM)
      if (NTEST.ge.5) then
         write(6,*) 'Highest GAS dimension is  ',MXGSDIM
      end if
*
      I_AM_OKAY = 1
      if (ISM.ne.JSM) then
         write(6,*) 'No integrals in this symmetry block'
         I_AM_OKAY = 0
         goto 100
      else
         if (NTEST.eq.5) then
            write(6,*) 'ISM equals JSM. I am fetching integrals'
         end if
*
*  determine real or imaginary integral parts
         SCALIN = 1.D0
         if (IRISIG.eq.1.and.IRICO.eq.1) then
            IRIINT = 1
         else if (IRISIG.eq.1.and.IRICO.eq.2) then
            IRIINT = 2
            SCALIN = -1.D0
         else if (IRISIG.eq.2.and.IRICO.eq.1) then
            IRIINT = 2
         else if (IRISIG.eq.2.and.IRICO.eq.2) then
            IRIINT = 1
         else
            write(6,*) 'IRI problem in geth1rel - returning'
            Call Abend('Quitting')
         end if
*
         if (ISUMNEW.eq.1.and.ENVIRO(1:6).ne.'DIRAC ') then
            if (IOPER.eq.1) then
               continue
            else if (IOPER.eq.2) then
               do LIST=1,NIRR_DG,1
                  call evenodd(IEOLIST,LIST)
                  do I=1,ITRALINE,1
                     if (IEOLIST.eq.1) then
                        IRUN = I + IOFFTRA(LIST+8) - 1
                        IRUN1 = I + IOFFTRA(LIST+16) - 1
                        SPININTS(IRUN) =
     &                     SPININTS(IRUN1)
                     else if (IEOLIST.eq.2) then
                        IRUN = I + IOFFTRA(LIST+8) - 1
                        IRUN2 = I + IOFFTRA(LIST+24) - 1
                        SPININTS(IRUN) =
     &                     SPININTS(IRUN)  +
     &                     SPININTS(IRUN2)
                     end if
                  end do
               end do
            else if (IOPER.eq.3) then
               do LIST=1,NIRR_DG,1
                  call evenodd(IEOLIST,LIST)
                  do I=1,ITRALINE,1
                     if (IEOLIST.eq.1) then
                        IRUN = I + IOFFTRA(LIST) - 1
                        IRUN2 = I + IOFFTRA(LIST+16) - 1
                        SPININTS(IRUN) =
     &                        SPININTS(IRUN) +
     &                        SPININTS(IRUN2)
                     else if (IEOLIST.eq.2) then
                        IRUN = I + IOFFTRA(LIST) - 1
                        IRUN1 = I + IOFFTRA(LIST+8) - 1
                        IRUN3 = I + IOFFTRA(LIST+24) - 1
                        SPININTS(IRUN) =
     &                        SPININTS(IRUN1) +
     &                        SPININTS(IRUN3)
                     end if
                  end do
               end do
            end if
         end if
*
*  Determine integral list
         INTLIST = 2*(ISM-NIRR_PN) + IRIINT - 2
         if (IOPER.eq.2) INTLIST = INTLIST + 8
*  List offset
         LSTOFF = IOFFTRA(INTLIST)
*  Offset for ITYP (inner index)
         ITPOFF = 0
         do I=1,ITYP-1,1
            ITPOFF = ITPOFF + 4*MXGSDIM**2
         end do
*  Offset for JTYP (outer index)
         JTPOFF = 0
         do J=1,JTYP-1,1
            JTPOFF = JTPOFF + 4*MXGSDIM**2*NGASPI
         end do
*  Final offset
         IFINOFF = LSTOFF + JTPOFF + ITPOFF
*
         if (IPM1.ne.0.and.IJUB.ne.0) then
            write(6,*) 'Combination of IJUB and IPM1 impossible'
            write(6,*) 'I am quitting'
            Call Abend('Quitting')
         end if
*
         if (NTEST.ge.5) then
           write(6,*)
           write(6,*) 'List and offsets in geth1rel3 '
           write(6,*) 'INTLIST    .... ',INTLIST
           write(6,*) 'LSTOFF     .... ',LSTOFF
           write(6,*) 'ITPOFF     .... ',ITPOFF
           write(6,*) 'JTPOFF     .... ',JTPOFF
           write(6,*) 'IFINOFF    .... ',IFINOFF
         end if
*
*  Determine relevant offsets and copy desired block to H
*
*  ISEARCH in the off-diagonal case depends on the symmetry of
*  the desired spinors. It has to be permuted, if the spinor
*  irrep is within the second half of the fermion irreps:
*
         if (IPM1.eq.1) then
            ISEARCH = 2
         else if (IPM1.eq.-1) then
            ISEARCH = 3
         else if (IPM1.eq.0) then
            if (IJUB.eq.1) then
               ISEARCH = 1
            else if (IJUB.eq.2) then
               ISEARCH = 4
            else
               write(6,*) 'IJUB out of range. Quitting.'
               Call Abend('Quitting')
            end if
         else
            write(6,*) 'IPM1 out of range. Quitting.'
            Call Abend('Quitting')
         end if
         ILEN = IPUBPGS(ISEARCH,ITYP,JTYP,INTLIST)
*
         ISTART = 1
         do I=1,ISEARCH-1,1
            ISTART = ISTART + IPUBPGS(I,ITYP,JTYP,INTLIST)
         end do
         IEND = ISTART + ILEN - 1
         K = 1
         ISTART = ISTART + IFINOFF - 1
         IEND = IEND + IFINOFF - 1
*
         if (NTEST.ge.10) then
           write(6,*) 'ISEARCH block  ',ISEARCH
           write(6,*) 'Block length   ',ILEN
           write(6,*) 'Start address  ',ISTART
           write(6,*) 'End address    ',IEND
         end if
*
* Now do the work ...
         do J=ISTART,IEND,1
            H(K) = SPININTS(J)
            K = K + 1
         end do
*
         if (NTEST.ge.10) then
            write(6,*) 'The integrals are:'
             do L=1,ILEN,1
                write(6,*) H(L)
             end do
            write(6,*) 'ISM is  ',ISM
            write(6,*) 'JSM is  ',JSM
            write(6,*) 'ITYP is  ',ITYP
            write(6,*) 'JTYP is  ',JTYP
            write(6,*) 'I have fetched integrals successfully'
            write(6,*) 'and am now returning'
         end if
         IRETURN = 1
      end if
100   continue
      if (I_AM_OKAY.eq.0) then
         write(6,*) 'I am returning without integrals'
         IRETURN = 0
      end if
*
* Run variable for summation of integrals. Only in the first call
* for geth1rel3, the summation should be carried out.
      ISUMNEW = ISUMNEW + 1
      if (NTEST.ge.5)
     &  write(6,*) 'ISUMNEW has current value  ',ISUMNEW
*
500   return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETINCN(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,INTLST,IJKLOF,NSMOB,I2INDX,
     &                  ICOUL)
*
* Obtain integrals
*
*     ICOUL = 0 :
*                  XINT(IK,JL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 1 :
*                  XINT(IJ,KL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
* Storing for ICOUL = 1 not working if IKSM or JLSM .ne. 0
*
*
* Version for integrals stored in INTLST
*
*  Debugged to some extent ...     Timo Fleig, August '98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "csm.inc"
*. Integral list
      Real * 8 Intlst(*)
      Dimension IJKLof(NSMSH,NSMSH,NSMSH)
*. Pair of orbital indeces ( symmetry ordered ) => address in symmetry packed
*. matrix
      Dimension I2INDX(*)
*.Output
      DIMENSION XINT(*)
*. Local scratch
      DIMENSION IJARR(MXPORB)
*
      NTEST = 000
*
      iOrb = NGSSH(ISM,ITP)
      jOrb = NGSSH(JSM,JTP)
      kOrb = NGSSH(KSM,KTP)
      lOrb = NGSSH(LSM,LTP)
*
      if (iOrb.eq.0.or.jOrb.eq.0.or.kOrb.eq.0
     &    .or.lOrb.eq.0) then
         if (NTEST.ge.100) write(6,*) ' No integrals fetched!'
         goto 666
      end if
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
     &   ' GETINCN : IORB,JORB,KORB,LORB',IORB,JORB,KORB,LORB
      END IF
*. Offsets relative to start of all orbitals, symmetry ordered
      IOFF = IBSH(ISM)
      DO IITP = 1, ITP -1
        IOFF = IOFF + NGSSH(ISM,IITP)
      END DO
*
      JOFF = IBSH(JSM)
      DO JJTP = 1, JTP -1
        JOFF = JOFF + NGSSH(JSM,JJTP)
      END DO
*
      KOFF = IBSH(KSM)
      DO KKTP = 1, KTP -1
        KOFF = KOFF + NGSSH(KSM,KKTP)
      END DO
*
      LOFF = IBSH(LSM)
      DO LLTP = 1, LTP -1
        LOFF = LOFF + NGSSH(LSM,LLTP)
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
     &   ' GETINCN : IOFF,jOFF,KOFF,LOFF',IOFF,jOFF,KOFF,LOFF
      END IF

*
*     Collect Coulomb terms
*
      ijblk = max(ism,jsm)*(max(ism,jsm)-1)/2 + min(ism,jsm)
      klblk = max(ksm,lsm)*(max(ksm,lsm)-1)/2 + min(ksm,lsm)
*
      IF(IJBLK.GT.KLBLK) THEN
       IJRELKL = 1
       IBLOFF=IJKLOF(MAX(ISM,JSM),MIN(ISM,JSM),MAX(KSM,LSM))
      ELSE IF (IJBLK.EQ.KLBLK) THEN
       IJRELKL = 0
       IBLOFF=IJKLOF(MAX(ISM,JSM),MIN(ISM,JSM),MAX(KSM,LSM))
      ELSE IF (IJBLK.LT.KLBLK) THEN
       IJRELKL = -1
       IBLOFF = IJKLOF(MAX(KSM,LSM),MIN(KSM,LSM),MAX(ISM,JSM))
      END IF
*
      itOrb=NGSSHT(iSm)
      jtOrb=NGSSHT(jSm)
      ktOrb=NGSSHT(kSm)
      ltOrb=NGSSHT(lSm)
*
      If(ISM.EQ.JSM) THEN
       IJPAIRS = ITORB*(ITORB+1)/2
      ELSE
       IJPAIRS = ITORB*JTORB
      END IF
*
      IF(KSM.EQ.LSM) THEN
        KLPAIRS = KTORB*(KTORB+1)/2
      ELSE
        KLPAIRS = KTORB*LTORB
      END IF
*
      iInt=0
      Do lJeppe=lOff,lOff+lOrb-1
        jMin=jOff
        If ( JLSM.ne.0 ) jMin=lJeppe
        Do jJeppe=jMin,jOff+jOrb-1
*
*
*. Set up array IJ*(IJ-1)/2
          IF(IJRELKL.EQ.0) THEN
            DO II = IOFF,IOFF+IORB-1
              IJ = I2INDX((JJEPPE-1)*NTOOB+II)
              IJARR(II) = IJ*(IJ-1)/2
            END DO
          END IF
*
          Do kJeppe=kOff,kOff+kOrb-1
            iMin = iOff
            kl = I2INDX(KJEPPE+(LJEPPE-1)*NTOOB)
            If(IKSM.ne.0) iMin = kJeppe
            IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                IINT = (LJEPPE-LOFF)*Jorb*Korb*Iorb
     &                + (KJEPPE-KOFF)*Jorb*Iorb
     &                + (JJEPPE-JOFF)*Iorb
            END IF
*
            IF(IJRELKL.EQ.1) THEN
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) > (ksm,lsm)
              IJKL0 = IBLOFF-1+(kl-1)*ijPairs
              IJ0 = (JJEPPE-1)*NTOOB
              Do iJeppe=iMin,iOff+iOrb-1
                  ijkl = ijkl0 + I2INDX(IJEPPE+IJ0)
                  iInt=iInt+1
                  Xint(iInt) = Intlst(ijkl)
              End Do
            END IF
*
*. block (ISM JSM !ISM JSM)
            IF(IJRELKL.EQ.0) THEN
              IJ0 = (JJEPPE-1)*NTOOB
              KLOFF = KL*(KL-1)/2
              IJKL0 = (KL-1)*IJPAIRS-KLOFF
              Do iJeppe=iMin,iOff+iOrb-1
                ij = I2INDX(IJEPPE+IJ0   )
                If ( ij.ge.kl ) Then
C                 ijkl=ij+(kl-1)*ijPairs-klOff
                  IJKL = IJKL0 + IJ
                Else
                  IJOFF = IJARR(IJEPPE)
                  ijkl=kl+(ij-1)*klPairs-ijOff
                End If
                iInt=iInt+1
                Xint(iInt) = Intlst(iblOff-1+ijkl)
              End Do
            END IF
*
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) < (ksm,lsm)
            IF(IJRELKL.EQ.-1) THEN
              ijkl0 = IBLOFF-1+KL - KLPAIRS
              IJ0 = (JJEPPE-1)*NTOOB
              Do iJeppe=iMin,iOff+iOrb-1
                IJKL = IJKL0 + I2INDX(IJEPPE + IJ0)*KLPAIRS
                iInt=iInt+1
                Xint(iInt) = Intlst(ijkl)
              End Do
            END IF
*
          End Do
        End Do
      End Do
*
*     Collect Exchange terms
*
      If ( IXCHNG.ne.0 ) Then
*
      IF(ISM.EQ.LSM) THEN
       ILPAIRS = ITORB*(ITORB+1)/2
      ELSE
       ILPAIRS = ITORB*LTORB
      END IF
*
      IF(KSM.EQ.JSM) THEN
        KJPAIRS = KTORB*(KTORB+1)/2
      ELSE
        KJPAIRS = KTORB*JTORB
      END IF
*
        ilblk = max(ism,lsm)*(max(ism,lsm)-1)/2 + min(ism,lsm)
        kjblk = max(ksm,jsm)*(max(ksm,jsm)-1)/2 + min(ksm,jsm)
        IF(ILBLK.GT.KJBLK) THEN
          ILRELKJ = 1
          IBLOFF = IJKLOF(MAX(ISM,LSM),MIN(ISM,LSM),MAX(KSM,JSM))
        ELSE IF(ILBLK.EQ.KJBLK) THEN
          ILRELKJ = 0
          IBLOFF = IJKLOF(MAX(ISM,LSM),MIN(ISM,LSM),MAX(KSM,JSM))
        ELSE IF(ILBLK.LT.KJBLK) THEN
          ILRELKJ = -1
          IBLOFF = IJKLOF(MAX(KSM,JSM),MIN(KSM,JSM),MAX(ISM,LSM))
        END IF
*
        iInt=0
        Do lJeppe=lOff,lOff+lOrb-1
          jMin=jOff
          If ( JLSM.ne.0 ) jMin=lJeppe
*
          IF(ILRELKJ.EQ.0) THEN
           DO II = IOFF,IOFF+IORB-1
             IL = I2INDX(II+(LJEPPE-1)*NTOOB)
             IJARR(II) = IL*(IL-1)/2
           END DO
          END IF
*
          Do jJeppe=jMin,jOff+jOrb-1
            Do kJeppe=kOff,kOff+kOrb-1
              KJ = I2INDX(KJEPPE+(JJEPPE-1)*NTOOB)
              KJOFF = KJ*(KJ-1)/2
              iMin = iOff
*
              IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                  IINT = (LJEPPE-LOFF)*Jorb*Korb*Iorb
     &                  + (KJEPPE-KOFF)*Jorb*Iorb
     &                  + (JJEPPE-JOFF)*Iorb
              END IF
*
              If(IKSM.ne.0) iMin = kJeppe
*
              IF(ILRELKJ.EQ.1) THEN
                ILKJ0 = IBLOFF-1+( kj-1)*ilpairs
                IL0 = (LJEPPE-1)*NTOOB
                Do iJeppe=iMin,iOff+iOrb-1
                  ILKJ = ILKJ0 + I2INDX(IJEPPE + IL0)
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(ilkj)
                End Do
              END IF
*
              IF(ILRELKJ.EQ.0) THEN
                IL0 = (LJEPPE-1)*NTOOB
                ILKJ0 = (kj-1)*ilPairs-kjOff
                Do iJeppe=iMin,iOff+iOrb-1
                  IL = I2INDX(IJEPPE + IL0 )
                  If ( il.ge.kj ) Then
C                     ilkj=il+(kj-1)*ilPairs-kjOff
                      ILKJ = IL + ILKJ0
                    Else
                      ILOFF = IJARR(IJEPPE)
                      ilkj=kj+(il-1)*kjPairs-ilOff
                    End If
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(iBLoff-1+ilkj)
                End Do
              END IF
*
              IF(ILRELKJ.EQ.-1) THEN
                ILKJ0 = IBLOFF-1+KJ-KJPAIRS
                IL0 = (LJEPPE-1)*NTOOB
                Do iJeppe=iMin,iOff+iOrb-1
                  ILKJ = ILKJ0 + I2INDX(IJEPPE+ IL0)*KJPAIRS
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(ilkj)
                End Do
              END IF
*
            End Do
          End Do
        End Do
      End If
*
      if (NTEST.ge.100) then
         call findhigh(NGSSH,MXPIRR*MXPNGAS,MAXDIM)
         ILENGTH = MAXDIM**4
         write(6,*) 'Checking XINT after getincn:'
         write(6,*) 'ISM  JSM  KSM  LSM'
         write(6,*) ISM,JSM,KSM,LSM
         do I=1,ILENGTH,1
            write(6,*) I,XINT(I)
         end do
      end if
 666  continue
*
      Return
      End
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
* Get 2-electron integrals from DIRAC environment structure
*
      subroutine getind(XOUT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  INVREO,KR,NUMINT,
     &                  INDK,INDL,IOFF12,
     &                  XIN)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
*--------------------------------------
*  Timo Fleig, December 2000
*--------------------------------------
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "input.inc"
*
      dimension XIN((2*NDIM_TOT)**4),XOUT(*)
      dimension KR(-NDIM_TOT:NDIM_TOT),INVREO(-NDIM_TOT:NDIM_TOT),
     &          INDK((2*NDIM_TOT)**4),INDL((2*NDIM_TOT)**4),
     &          IOFF12(MORB_DIRAC,MORB_DIRAC),
     &          NUMINT(MORB_DIRAC,MORB_DIRAC)
*
* 1) Determine order of output integrals (GAS on SYM block)
* 2) Determine index quadruples for subsequent integrals
* 3) Translate to input (DIRAC) index scheme
* 4) Fetch integral and add to list
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
* Type offset within one symmetry
      IGSOFF = 0
      JGSOFF = 0
      KGSOFF = 0
      LGSOFF = 0
      do IGAS=1,ITP-1,1
        IGSOFF = IGSOFF + NGSSH(ISM,IGAS)
      end do
      do JGAS=1,JTP-1,1
        JGSOFF = JGSOFF + NGSSH(JSM,JGAS)
      end do
      do KGAS=1,KTP-1,1
        KGSOFF = KGSOFF + NGSSH(KSM,KGAS)
      end do
      do LGAS=1,LTP-1,1
        LGSOFF = LGSOFF + NGSSH(LSM,LGAS)
      end do
*
* Loop boundaries (including total symmetry offset)
      ISTA = IOFPGSM(ISM) + IGSOFF
      JSTA = IOFPGSM(JSM) + JGSOFF
      KSTA = IOFPGSM(KSM) + KGSOFF
      LSTA = IOFPGSM(LSM) + LGSOFF
      IEND = IOFPGSM(ISM) + IGSOFF + NGSSH(ISM,ITP)
      JEND = IOFPGSM(JSM) + JGSOFF + NGSSH(JSM,JTP)
      KEND = IOFPGSM(KSM) + KGSOFF + NGSSH(KSM,KTP)
      LEND = IOFPGSM(LSM) + LGSOFF + NGSSH(LSM,LTP)
*
      if (NTEST.ge.5) then
        write(6,'(A,4I4)') 
     &        'Start of loop per index: ',ISTA+1,JSTA+1,KSTA+1,LSTA+1
        write(6,'(A,4I4)') 
     &        'End   of loop per index: ',IEND,JEND,KEND,LEND
      end if
*
* Loop over input indices in required order
      ICT = 1
      ZERO = 0.D0
      do L=LSTA+1,LEND,1
        do J=JSTA+1,JEND,1
          do K=KSTA+1,KEND,1
            do I=ISTA+1,IEND,1
*
              IIN = INVREO(I)
              JIN = INVREO(J)
              KIN = INVREO(K)
              LIN = INVREO(L)
*
              if (NTEST.ge.5) then
                write(6,*) 'Kramers indices:'
                write(6,'(A,4I6)') 
     &            'KR(NIN) ',KR(IIN),KR(JIN),KR(KIN),KR(LIN)
              end if
*
*  Find non-vanishing Kramers block and define sign (KRFAC)
              IRIFAC = 1
              KRFAC = 1   ! always here
              if (IDGTYP.eq.1.or.IDGTYP.eq.2) then
                IBLLEN = NUMINT(KR(IIN),KR(JIN))
                if (JIN.gt.IIN) then
*  use real/imaginary symmetry   
                  if (IRI.eq.2) IRIFAC = -1
*  integral is imaginary, i.e.
*   (rr|ri) or
*   (ri|ii)
*   or all kinds of permutations of the latter.
*   In any case: (IJ|KL) = -(JI|LK) here.
                  IIN_SAVE = IIN
                  JIN_SAVE = JIN
                  IIN = JIN_SAVE
                  JIN = IIN_SAVE
                  KIN_SAVE = KIN
                  LIN_SAVE = LIN
                  KIN = LIN_SAVE
                  LIN = KIN_SAVE
                  IBLLEN = NUMINT(KR(IIN),KR(JIN))
                end if
                if (IBLLEN.eq.0) then
*  pick (IJ_|KL_)
                  JIN = -JIN
                  LIN = -LIN
                  IBLLEN = NUMINT(KR(IIN),KR(JIN))
                end if
                if (IBLLEN.eq.0) then
*  pick (I_J|K_L)
                  IIN = -IIN
                  JIN = -JIN
                  KIN = -KIN
                  LIN = -LIN
                  IBLLEN = NUMINT(KR(IIN),KR(JIN))
                end if
                if (IBLLEN.eq.0) then
*  pick (I_J_|K_L_)
                  JIN = -JIN
                  LIN = -LIN
                  IBLLEN = NUMINT(KR(IIN),KR(JIN))
                end if
                if (IBLLEN.eq.0) then
                  write(6,*) 'Dist_2eints_dirac: '
                  write(6,*) 'Seeking indices ',IIN,JIN
                  write(6,*) 'Vanishing Kramers block.'
                  call abend2('Stop.')
                end if
*
                if (NTEST.ge.5) then
                  write(6,*) 
     &            'Kramers indices after symm. reordering:'
                  write(6,'(A,4I6)') 
     &            'KR(NIN) ',KR(IIN),KR(JIN),KR(KIN),KR(LIN)
                end if
*
*  Fetch integral
                IBLKOFF = IOFF12(KR(IIN),KR(JIN))
                do INRUN=1,IBLLEN,1
                  if (KR(INDK(INRUN+IBLKOFF)).eq.KR(KIN).
     &            and.KR(INDL(INRUN+IBLKOFF)).eq.KR(LIN))
     &            then
                    XOUT(ICT) = IRIFAC * KRFAC * XIN(INRUN+IBLKOFF)
                    if (NTEST.ge.10) then
                      write(6,'(A,3I4,1F18.12)') 
     &                    'ICT,IBLKOFF,INRUN,XOUT(ICT) '
     &                    ,ICT,IBLKOFF,INRUN,XOUT(ICT)
                    end if
                    ICT = ICT + 1
                    goto 10
                  end if
                end do
*
*  Integral has not been found on existing symmetry (IJ| block
*  of input integrals. It has been deleted in DIRAC due to small
*  value. Insert ZERO in its place.
                XOUT(ICT) = ZERO
                if (NTEST.ge.10) then
                  write(6,'(10X,A,I9,A,1F18.12)') 
     &                'Integral ',ICT,' set to ',ZERO
                end if
                ICT = ICT + 1
              else if (IDGTYP.eq.3) then
                write(6,*) 'Complex matrix groups'
                write(6,*) 'not implemented yet.'
                call abend2('Stop.')
              end if
   10         continue
*             ^ Generate next output index quadruple
            end do
          end do
        end do
      end do
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Get 2-electron integral with selected indices
      subroutine getind_index(YOUT,IIN,JIN,KIN,LIN,NUMINT,
     &                        KR,INDK,INDL,IOFF12,IRI,XIN,IPRNT)
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "symm.inc"
#include "crun.inc"
*
      dimension NUMINT(MORB_DIRAC,MORB_DIRAC,4),KR(-NDIM_TOT:NDIM_TOT),
     &          INDK((2*NDIM_TOT)**4),INDL((2*NDIM_TOT)**4),
     &          IOFF12(MORB_DIRAC,MORB_DIRAC,4)
      dimension XIN((2*NDIM_TOT)**4)
*
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
      if (NTEST.ge.10) then
        write(6,*) '========================'
        write(6,*) ' GETIND_INDEX called    '
        write(6,*) '========================'
        write(6,*) 'Double group      ',DOUGRP
        write(6,*) 'Double group type ',IDGTYP
        write(6,'(A,4I4)') 'Indices : ',IIN,JIN,KIN,LIN
      end if
*
*  Find non-vanishing Kramers block and define sign (KRFAC)
*  Initialize
      ZERO = 0.D0
      IRIFAC = 1
      KRFAC = 1
*
**************************************************
* Spinfree formalism of DIRAC environment        *
**************************************************
      if (ISPINFREE.eq.1) then
        IT = 1
        IBLLEN = NUMINT(KR(IIN),KR(JIN),IT)
        if (JIN.gt.IIN) then
*  use real/imaginary symmetry   
          if (IRI.eq.2) then
            IRIF = -1
          else
            IRIF = 1
          end if
*  integral is imaginary, i.e.
*   (rr|ri) or
*   (ri|ii)
*   or all kinds of permutations of the latter.
*   In any case: (IJ|KL) = -(JI|LK) here.
          call index_exch(IIN,JIN,IRIFAC,1,1)
          call index_exch(KIN,LIN,IRIFAC,1,IRIF)
          IBLLEN = NUMINT(KR(IIN),KR(JIN),IT)
        end if
        if (IBLLEN.eq.0) then
* The following symmetry relations can only be used in SPINFREE
* calculations. In the full dirac picture, e.g., 
* (IJ|KL) .ne. (IJ_|KL_) (see below)
*  pick (IJ_|KL_)
          JIN = -JIN
          LIN = -LIN
          IBLLEN = NUMINT(KR(IIN),KR(JIN),IT)
        end if
        if (IBLLEN.eq.0) then
*  pick (I_J|K_L)
          IIN = -IIN
          JIN = -JIN
          KIN = -KIN
          LIN = -LIN
          IBLLEN = NUMINT(KR(IIN),KR(JIN),IT)
        end if
        if (IBLLEN.eq.0) then
*  pick (I_J_|K_L_)
          JIN = -JIN
          LIN = -LIN
          IBLLEN = NUMINT(KR(IIN),KR(JIN),IT)
        end if
*
      else if (ISPINFREE.eq.0) then
*
**************************************************
* Full spin-dependent spinors as basis functions *
**************************************************
*
* Determine relevant block for first particle
* and corresponding subblock for second particle.
*
* Revised algorithm.
* 1) Count number of barred indices:
        NBAR = 0
        if (IIN.lt.0) NBAR = NBAR + 1
        if (JIN.lt.0) NBAR = NBAR + 1
        if (KIN.lt.0) NBAR = NBAR + 1
        if (LIN.lt.0) NBAR = NBAR + 1

*
* 2) Distinguish different MK cases:
* IT default 1 for first block referring to index pair (IJ|
        IT = 1
*************************************
        if (NBAR.eq.0) then
*************************************
*         (uu|uu) class.
*************************************
        else if (NBAR.eq.1) then
*************************************
*         (ub|uu) DMK = +1 or 
*         (bu|uu) DMK = -1 class.
          if (JIN.lt.0) then
            I1KSYM = -1
            I2KSYM = 1
            if (IDGTYP.eq.1) IT = 1  !!?? FIXME: check this
            if (IDGTYP.eq.2) IT = 1  !!?? FIXME: check this
            if (IDGTYP.eq.3) IT = 3
            if (abs(JIN).gt.IIN) then
              if (KRASYM) then
                IKF = -1
                call index_exch(IIN,JIN,KRFAC,IKF,I1KSYM)
              end if
            end if
          else
            I1KSYM = -1
            I2KSYM = 1
          end if
*************************************
        else if (NBAR.eq.2) then
*************************************
*         (uu|bb) DMK = 0 or
*         (ub|bu) DMK = 0 or
*         (ub|ub) DMK = +2 or
*         (bu|bu) DMK = -2 class.
          if (JIN.lt.0) then
            if (KIN.lt.0) then
*             (ub|bu) DMK = 0
              I1KSYM = -1
              I2KSYM = -1
              IT = 2
              if (abs(JIN).gt.IIN) then
                if (KRASYM) then
                  IKF = -1
                  call index_exch(IIN,JIN,KRFAC,IKF,I1KSYM)
                  call index_exch(KIN,LIN,KRFAC,IKF,I2KSYM)
                end if
              end if
            else
*             (ub|ub) DMK = +2
              I1KSYM = -1
              I2KSYM = -1
              if (IIN.lt.abs(JIN)) then
                if (KRASYM) then
                  IKF = -1
                  call index_exch(IIN,JIN,KRFAC,IKF,I1KSYM)
                  call index_exch(KIN,LIN,KRFAC,IKF,I2KSYM)
                end if
              end if
            end if
          else
            if (IIN.lt.0) then
*             (bu|bu) DMK = -2
              I1KSYM = -1
              I2KSYM = -1
              if (abs(IIN).gt.JIN) then
                if (KRASYM) then
                  IKF = -1
                  call index_exch(IIN,JIN,KRFAC,IKF,I1KSYM)
                  call index_exch(KIN,LIN,KRFAC,IKF,I2KSYM)
                end if
              end if
            else
*             (uu|bb) DMK = 0
              I1KSYM = 1
              I2KSYM = 1
              if (IDGTYP.eq.1.or.IDGTYP.eq.2) then
                if (KRASYM) then
                  IKF = -1
                  call index_exch(KIN,LIN,KRFAC,IKF,I2KSYM)
                end if
              else if (IDGTYP.eq.3) then
                IT = 2
              end if
            end if
          end if
*************************************
        else if (NBAR.eq.3) then
*************************************
*         (ub|bb) DMK = +1 or
*         (bu|bb) DMK = -1 class.
          I1KSYM = -1
          I2KSYM = 1
          if (IIN.lt.0) then
          else
            if (IDGTYP.eq.1) IT = 1  !!??  FIXME: check this
            if (IDGTYP.eq.2) IT = 1  !!??  FIXME: check this
            if (IDGTYP.eq.3) IT = 3
          end if
          if (KRASYM) then
            IKF = -1
            call index_exch(KIN,LIN,KRFAC,IKF,I2KSYM)
          end if
*************************************
        else if (NBAR.eq.4) then
*************************************
*         (bb|bb) class.
          I1KSYM = 1
          I2KSYM = 1
          if (KRASYM) then
            IKF = -1
            call index_exch(IIN,JIN,KRFAC,IKF,I1KSYM)
            call index_exch(KIN,LIN,KRFAC,IKF,I2KSYM)
          end if
        end if
*
* Find relevant block for first particle.
        IBLLEN = NUMINT(KR(IIN),KR(JIN),IT)
*
* Final step: Use complex conjugation symmetry.
        if (IBLLEN.eq.0) then
          if (NTEST.ge.5) 
     &        write(6,*) 'CAUTION ! Using c.c. symmetry on indices!'
          if (KIN.lt.0.and.LIN.lt.0) then
            IT = 2
          else
            IT = 1      !!??  FIXME : At least questionable !!
          end if
          if (IRI.eq.2) then
            IRIF = -1
          else
            IRIF = 1
          end if
          IKF = 1
          call index_exch(IIN,JIN,IRIFAC,IKF,1)
          call index_exch(KIN,LIN,IRIFAC,IKF,IRIF)
          IBLLEN = NUMINT(KR(IIN),KR(JIN),IT)
        end if
*
      end if
*
      if (NTEST.ge.5) then
        write(6,*) 'Kramers indices after symm. reordering:'
        write(6,'(A,4I6)') 'KR(NIN) ',KR(IIN),KR(JIN),KR(KIN),KR(LIN)
      end if
*
*  Fetch integral
      IBLKOFF = IOFF12(KR(IIN),KR(JIN),IT)
      if (NTEST.ge.5) then
        write(6,*) ' IT,IBLKOFF,IBLLEN,KIN,LIN :',
     &               IT,IBLKOFF,IBLLEN,KIN,LIN
      end if
      do INRUN=1,IBLLEN,1
        if (KR(INDK(INRUN+IBLKOFF)).eq.KR(KIN).
     &  and.KR(INDL(INRUN+IBLKOFF)).eq.KR(LIN))
     &  then
          YOUT = IRIFAC * KRFAC * XIN(INRUN+IBLKOFF)
          if (NTEST.ge.10) then
            write(6,'(A,2I4)') 'IRIFAC,KRFAC ',IRIFAC,KRFAC
            write(6,'(A,2I4,1F18.12)') 
     &          'IBLKOFF,INRUN,YOUT '
     &          ,IBLKOFF,INRUN,YOUT
          end if
          goto 10
        end if
      end do
*
*  Integral has not been found on existing symmetry (IJ| block
*  of input integrals. It has been deleted in DIRAC due to small
*  value. Insert ZERO in its place.
      YOUT = ZERO
      if (NTEST.ge.10) then
        write(6,'(10X,A,1F18.12)') 
     &      'Integral set to ',ZERO
      end if
   10 continue
*     ^ Generate next output index quadruple
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Determine 1- or 2- particle integral
* Fetch integral.
*
      subroutine getind_molf(YINT,INDEX,IRI,IDIM,INVREO,
     &                       KR,NASHT,DMFC,VMUUUU,IPRNT)
*
#include "implicit.h"
#include "ipoist8.inc"
*
      dimension VMUUUU(*),DMFC(*)
      dimension ICSAVE(4),INDEX(4),INVREO(-IDIM:IDIM),KR(-IDIM:IDIM)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
      if (NTEST.ge.5) then
        write(6,*)
        write(6,*) ' GETIND_MOLF : INDEX array'
        write(6,*) (INDEX(I),I=1,4,1)
      end if
*
      if (NTEST.ge.20) then
        write(6,*) '1-el. integral array in getind_molf:'
        call wrtmat(DMFC,IDIM*2,IDIM*2,IDIM*2,IDIM*2)
      end if
*
      NCT = 0
      do IC=1,4,1
        if (INDEX(IC).ne.0) then
          NCT = NCT + 1
          ICSAVE(NCT) = IC
        end if
      end do
*
      if (NCT.eq.2) then
        if (NTEST.ge.10) write(6,*) '1-electron integral.'
        I = KR(INVREO(INDEX(ICSAVE(1))))
        J = KR(INVREO(INDEX(ICSAVE(2))))
* Compute input address
        IIN =  (IRI-1)*NASHT*NASHT +
     &          (J-1) *NASHT +
     &           I
        if (NTEST.ge.5) then
          write(6,*) 'Index pair ',I,J,' address :',IIN
        end if
        YINT = DMFC(IIN)
*
* Sign for imaginary 1-electron integrals:
* (According to h_(IJ) = -h_(I_J_).
*  Reference are the unbarred spinors.)
        if (IRI.eq.2.and.INDEX(1).eq.0) YINT = - YINT
*
        if (NTEST.ge.2) write(6,'(2I3,2X,1F18.10)') I,J,YINT
*
      else if (NCT.eq.4) then
c       print*,'INDEX : ',(INDEX(IP),IP=1,4,1)
c       print*,'INVREO : ',(INVREO(IP),IP=-3,3,1)
c       print*,'KR : ',(KR(IP),IP=-3,3,1)
        I = KR(INVREO(INDEX(1)))
        J = KR(INVREO(INDEX(2)))
        K = KR(INVREO(INDEX(3)))
        L = KR(INVREO(INDEX(4)))
*
* Compute input address
        IIN = (IRI-1)*NASHT*NASHT*NASHT*NASHT +
     &         (L-1) *NASHT*NASHT*NASHT +
     &         (K-1) *NASHT*NASHT +
     &         (J-1) *NASHT +
     &          I
        if (NTEST.ge.5) then
          write(6,*) 'Index quad ',I,J,K,L,' address :',IIN
        end if
        YINT = VMUUUU(IIN)
        if (NTEST.ge.2) write(6,'(A,4I3,2X,1F18.10)')
     &                  'MOLF: ',I,J,K,L,YINT
c       stop
*
      else
        write(6,*) 'No QFT implemented YET in LUCIAREL.'
        write(6,*) 'you specified ',NCT,' indices.'
        call abend2('Quitting in getind_molf.')
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Determine 2- particle integral
* Fetch integral.
*
*=========================
* Timo Fleig, Oct 20, 2004
*=========================
*
      subroutine getind_nz3(YINT,INDEX,IRI,IDIM,INVREO,
     &                      KR,NASHT,NNASHX,H2AC,IPRNT)
*
#include "implicit.h"
#include "dgroup.h"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
#endif
#include "parluci.h"
*
      dimension H2AC(NASHT/2,NASHT/2,NNASHX,NZ,*)
      dimension ICSAVE(4),INDEX(4),INVREO(-IDIM:IDIM),KR(-IDIM:IDIM),
     &          INDEX_LOC(4)
*
      do IL = 1,4,1
        INDEX_LOC(IL) = INDEX(IL)
      end do
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRNT)
*
      if (IRI.eq.1) then
        MZ = 1
      else if (IRI.eq.2) then
        MZ = NZ
      end if
*
      if (NTEST.ge.2) then
        write(6,'(/A,6X,4I3)') '  INDEX array: ',(INDEX(I),I=1,4,1)
      end if
      if (NTEST.ge.5) then
        write(6,'(/A)') ' GETIND_NZ3 : '
        write(6,*) '  NZ,MZ : ',NZ,MZ
        write(6,*) '  No. of active orbitals : ',NASHT/2
        write(6,*) '  Compound index dimension: ',NNASHX
      end if
*
      if (NTEST.ge.200) then
        if( .not. lowsrt_ijkl ) then
          write(6,'(/A)')'2-el. integrals in getind_nz3, (NZ,3) format:'
          call prdnz3(H2AC(1,1,1,1,1),NASHT/2,NNASHX,NZ,IPQTOQ(1,0),6)
        end if
      end if
*
      NCT = 0
      do IC=1,4,1
        if (INDEX_LOC(IC).ne.0) then
          NCT = NCT + 1
          ICSAVE(NCT) = IC
        end if
      end do
*
      if (NCT.eq.4) then
*
* Interface to (NZ,3) format for 2-electron integrals
*
* Use Kramers symmetry/particle symmetry to generate a stored  
* index sequence on (NZ,3) format:
*  NZ3TP = 1     (pq|rs)
*  NZ3TP = 2     (pq_|rs_)
*  NZ3TP = 3     (p_q|rs_)
          IRIFAC = 1
          ISEQTP = 0
          if (INDEX_LOC(1).lt.0.and.INDEX_LOC(2).lt.0.and.
     &        INDEX_LOC(3).lt.0.and.INDEX_LOC(4).lt.0) then
            call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,-1,1)
            call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,-1,1)
            NZ3TP = 1
            ISEQTP = 1
          else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).gt.0.and.
     &             INDEX_LOC(3).gt.0.and.INDEX_LOC(4).gt.0) then
            NZ3TP = 1
            ISEQTP = 2
          else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).gt.0.and.
     &             INDEX_LOC(3).lt.0.and.INDEX_LOC(4).lt.0) then
            call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,-1,1)
            NZ3TP = 1
            ISEQTP = 3
          else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).lt.0.and.
     &             INDEX_LOC(3).lt.0.and.INDEX_LOC(4).gt.0) then
            call index_exch(INDEX_LOC(1),INDEX_LOC(3),1,1,1)
            call index_exch(INDEX_LOC(2),INDEX_LOC(4),1,1,1)
            INDEX_LOC(1) = -INDEX_LOC(1)
            INDEX_LOC(4) = -INDEX_LOC(4)
            NZ3TP = 3
            ISEQTP = 4
          else if (INDEX_LOC(1).gt.0.and.INDEX_LOC(2).lt.0.and.
     &             INDEX_LOC(3).gt.0.and.INDEX_LOC(4).lt.0) then
            INDEX_LOC(2) = -INDEX_LOC(2)
            INDEX_LOC(4) = -INDEX_LOC(4)
            NZ3TP = 2
            ISEQTP = 5
          else if (INDEX_LOC(1).lt.0.and.INDEX_LOC(2).gt.0.and.
     &             INDEX_LOC(3).lt.0.and.INDEX_LOC(4).gt.0) then
            if (IRI.eq.2) IRIFAC = -1
            call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,1,1)
            call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,1,1)
            INDEX_LOC(2) = -INDEX_LOC(2)
            INDEX_LOC(4) = -INDEX_LOC(4)
            NZ3TP = 2
            ISEQTP = 6
          end if
          if (NTEST.ge.2) write(6,'(A,1X,4I3)') 
     &       '  After 3 list step:', (INDEX_LOC(I),I=1,4,1)
*
* Find smallest index and its position in Mulliken notation
          call findlow(INDEX_LOC,4,ILOW,IPOS)
*
* Transfer smallest to fourth index by using 
*   particle symmetry 
*   Kramers symmetry
*   complex conjugation symmetry (IRIFAC2 = -1 if imaginary part)
*   CAUTION HERE:
*    The sign has to be derived from the original index sequence!
*    Thus the below sign arrangements may look queer.
*    They are nevertheless correct.
*
          IRIFAC2 = 1
          if (IPOS.eq.1) then 
            call index_exch(INDEX_LOC(1),INDEX_LOC(3),1,1,1)
            call index_exch(INDEX_LOC(2),INDEX_LOC(4),1,1,1)
            if (IRI.eq.2.and.NZ3TP.ne.2) IRIFAC2 = -1
            call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,1,1)
            call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,1,1)
          else if (IPOS.eq.2) then
            if (IRI.eq.2.and.NZ3TP.eq.3) IRIFAC2 = -1
            call index_exch(INDEX_LOC(1),INDEX_LOC(3),1,1,1)
            call index_exch(INDEX_LOC(2),INDEX_LOC(4),1,1,1)
          else if (IPOS.eq.3) then
            if (IRI.eq.2.and.NZ3TP.eq.1) IRIFAC2 = -1
            call index_exch(INDEX_LOC(1),INDEX_LOC(2),1,1,1)
            call index_exch(INDEX_LOC(3),INDEX_LOC(4),1,1,1)
          end if
          if (NTEST.ge.2) write(6,'(A,1X,4I3)') 
     &       '  After posit. step:', (INDEX_LOC(I),I=1,4,1)
*
* Determine address in (NZ,3) format
          NASHTH = NASHT/2
*
* Determine compound index KL: NDIM2
          NDIM2 = 0
          do IC = 1,INDEX_LOC(3),1
            if (IC.eq.INDEX_LOC(3)) then
              JCEND = INDEX_LOC(4)
            else if (IC.ne.INDEX_LOC(3)) then
              JCEND = IC
            end if
            do JC = 1,JCEND,1
              NDIM2 = NDIM2 + 1
            end do
          end do
*
          IIN = (NZ3TP-1)        *MZ*NNASHX*NASHTH*NASHTH +
     &          (MZ-1)           *NNASHX*NASHTH*NASHTH +
     &          (NDIM2-1)        *NASHTH*NASHTH +
     &          (INDEX_LOC(2)-1) *NASHTH +
     &           INDEX_LOC(1)
*
C         YINT = H2AC(IIN) * IRIFAC * IRIFAC2
          if( lowsrt_ijkl ) then
#if defined (VAR_MPI2)
            yint4 = 0.0D0
!           file offset is 0!
            ioffset = iin - 1
            call mpi_file_read_at(ijkl_pre,ioffset,yint4,
     &                            1,mpi_real8,istat,ierr)
!           write(6,*) 'iin ,ioffset, yint4',iin ,ioffset, yint4
            yint = yint4 * irifac * irifac2
#endif
          else
            YINT = H2AC(INDEX_LOC(1),INDEX_LOC(2),NDIM2,MZ,NZ3TP) 
     &                  * IRIFAC * IRIFAC2
          end if
!         ntest = 2
          if (NTEST.ge.2) write(6,'(A,5I3,2X,4I3,2X,1F18.10)') 
     &                          'MZ,NZ3TP,IRIFAC,IRIFAC2,IIN  integral',
     &                           MZ,NZ3TP,IRIFAC,IRIFAC2,IIN,
     &                           (INDEX_LOC(ID),ID=1,4,1),YINT
!         ntest = 0
*
      else
        write(6,*) 'No QFT implemented YET in LUCIAREL.'
        write(6,*) 'you specified ',NCT,' indices.'
        call abend2('Quitting in getind_nz3.')
      end if
*
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Slave (yes, slave. We are NOT politically correct in 
* FORTRAN programs!) routine for fetching desired integral
* from DIRAC list; both 1- and 2-electron
*
      subroutine getind_sigden(YINT,INDEX,IRI,ITOTDIM,
     &                         NUMINT,KR,INDK,INDL,INVREO,
     &                         IOFF12,XIN,ONEINT,IPRNT)
*
#include "implicit.h"
#include "ipoist8.inc"
*
      dimension INDEX(4),ICSAVE(4),INVREO(-ITOTDIM:ITOTDIM)
      dimension ONEINT(ITOTDIM*2,ITOTDIM*2)
*
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
      if (NTEST.ge.20) then
        write(6,*) '1-el. integral array in getind_sigden:'
        call wrtmat(ONEINT,ITOTDIM*2,ITOTDIM*2,ITOTDIM*2,ITOTDIM*2)
      end if
*
* 1- or 2-electron integral?
      NCT = 0
      do IC=1,4,1
        if (INDEX(IC).ne.0) then
          NCT = NCT + 1
          ICSAVE(NCT) = IC
        end if
      end do 
      if (NCT.eq.2) then
        if (NTEST.ge.10) then
          write(6,*) '1-electron integral.'
          write(6,'(A,2I4)') 'Index  pair  : ',INDEX(ICSAVE(1)),
     &                                         INDEX(ICSAVE(2))
          write(6,'(A,2I4)') 'reordered to : ',INVREO(INDEX(ICSAVE(1))),
     &                                         INVREO(INDEX(ICSAVE(2)))
        end if
        I = INVREO(INDEX(ICSAVE(1)))
        J = INVREO(INDEX(ICSAVE(2)))
        YINT = ONEINT(I,J)
*
* Sign for imaginary 1-electron integrals:
* (According to h_(IJ) = -h_(I_J_) )
        if (IRI.eq.2.and.INDEX(2).eq.0) YINT = - YINT
* FIXME : Why is it the unbarred integrals have to be multiplied??
*
      else if (NCT.eq.4) then
        if (NTEST.ge.10) then
          write(6,*) '2-electron integral.'
          write(6,'(A,4I4)') 'Index quadr. : ',INDEX(1),
     &                                         INDEX(2),
     &                                         INDEX(3),
     &                                         INDEX(4)
          write(6,'(A,4I4)') 'reordered to : ',INVREO(INDEX(1)),
     &                                         INVREO(INDEX(2)),
     &                                         INVREO(INDEX(3)),
     &                                         INVREO(INDEX(4))
        end if
        I = INVREO(INDEX(1))
        J = INVREO(INDEX(2))
        K = INVREO(INDEX(3))
        L = INVREO(INDEX(4))
*
        call getind_index(YINT,I,J,K,L,NUMINT,KR,INDK,INDL,
     &                    IOFF12,IRI,XIN,IPRNT)
*	
      else
        write(6,*) 'No QFT implemented YET in LUCIAREL.'
        write(6,*) 'you specified ',4-NCT,' indices.'
        write(6,*) 'INDEX array: ',(INDEX(I),I=1,4,1)
        call abend2('Quitting.')
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETINL(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM)
*
* Obtain integrals XINT(IK,JL) = (IJ!KL) for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*.ORBINP
C     PARAMETER(MXPORB=500)
C     PARAMETER(MXPIRR = 20)
C     PARAMETER(MXPOBS = 20)
C     PARAMETER(MXPR4T = 10)
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB)
#include "mxpdim.inc"
#include "orbinp.inc"
*.Output
      DIMENSION XINT(*)
*
      CALL QENTER('GETIN')
      NI = NTSOB(ITP,ISM)
      NJ = NTSOB(JTP,JSM)
      NK = NTSOB(KTP,KSM)
      NL = NTSOB(LTP,LSM)
*
      IF(IKSM.NE.0) THEN
        NIK = NI*(NI+1)/2
      ELSE
        NIK = NI*NK
      END IF
      IF(JLSM.NE.0) THEN
        NJL = NJ*(NJ+1)/2
      ELSE
        NJL = NJ*NL
      END IF
*
      I = 0
      K = 1
      DO 100 IK = 1, NIK
C       CALL NXTIJ(L,J,NL,NJ,JLSM,NONEW)
        CALL NXTIJ(I,K,NI,NK,IKSM,NONIK)
        IABS_T = I + IBTSOB(ITP,ISM)-1
        KABS_T = K + IBTSOB(KTP,KSM)-1
        J = 0
        L = 1
        DO 50 JL = 1, NJL
          CALL NXTIJ(J,L,NJ,NL,JLSM,NONJL)
          LABS_T = L + IBTSOB(LTP,LSM)-1
          JABS_T = J + IBTSOB(JTP,JSM)-1
C         XINT((IK-1)*NJL+JL) = GTIJKL_REL(IABS_T,JABS_T,KABS_T,LABS_T)
C         IF(IXCHNG.EQ.1)
C    &    XINT((IK-1)*NJL+JL) = XINT((IK-1)*NJL+JL)
C    &    - GTIJKL_REL(IABS_T,LABS_T,KABS_T,JABS_T)
          XINT((JL-1)*NIK+IK) = GTIJKL_REL(IABS_T,JABS_T,KABS_T,LABS_T)
          IF(IXCHNG.EQ.1)
     &    XINT((JL-1)*NIK+IK) = XINT((JL-1)*NIK+IK)
     &    - GTIJKL_REL(IABS_T,LABS_T,KABS_T,JABS_T)
   50   CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' 2 electron integral block for TS blocks '
        WRITE(6,'(1H ,4(A,I2,A,I2,A))')
     &  '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &  '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
        CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
      END IF
*
      CALL QEXIT('GETIN')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETINM(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,NORBMOL,NIRR_PN)
*
* Obtain integrals XINT(IK,JL) = (IJ!KL) for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "clunit.inc"
*.LUNDIO
      Parameter ( mxBatch = 106  )
      Parameter ( mxSyBlk = 666  )
      Common / LundIOR / LuTr2,lTr2Rec,iDAdr(mxBatch),nBatch(mxSyBlk)
      dimension NORBMOL(NIRR_PN)
*.Output
      DIMENSION XINT(*)
*.Local
      Parameter ( lBuf    = 9600 )
      Dimension Scr(lBuf)
      Logical iSymj,kSyml,iSymk,jSyml,ijSymkl,ikSymjl
      Logical ijPerm,klPerm,ijklPerm
*
      iOrb=NOBPTS(ITP,ISM)
      jOrb=NOBPTS(JTP,JSM)
      kOrb=NOBPTS(KTP,KSM)
      lOrb=NOBPTS(LTP,LSM)
      iOff=IOBPTS(ITP,ISM)
      jOff=IOBPTS(JTP,JSM)
      kOff=IOBPTS(KTP,KSM)
      lOff=IOBPTS(LTP,LSM)
*
      IAD13 = 0
      NLONGER = 2
      call findhigh(NORBMOL,NIRR_PN,MAXDMMO)
      write(6,*) 'Highest symmetry dimension is  ',MAXDMMO
      KBUFF = MAXDMMO**4
*
*     Collect Coulomb terms
*
      ijPerm=.false.
      klPerm=.false.
      ijklPerm=.false.
      If ( iSm.gt.jSm ) Then
        iSym=iSm
        jSym=jSm
      Else
        iSym=jSm
        jSym=iSm
        ijPerm=.true.
      End If
      ijBlk=jSym+iSym*(iSym-1)/2
      If ( kSm.gt.lSm ) Then
        kSym=kSm
        lSym=lSm
      Else
        kSym=lSm
        lSym=kSm
        klPerm=.true.
      End If
      klBlk=lSym+kSym*(kSym-1)/2
      If ( klBlk.gt.ijBlk ) Then
        iTemp=iSym
        iSym=kSym
        kSym=iTemp
        iTemp=jSym
        jSym=lSym
        lSym=iTemp
        iTemp=ijBlk
        ijBlk=klBlk
        klBlk=iTemp
        ijklPerm=.true.
      End If
C?    print *,' GetInt reports to you'
C?    print *,' iSym,jSym,kSym,lSym',iSym,jSym,kSym,lSym
      iSyBlk=klBlk+ijBlk*(ijBlk-1)/2
      iBatch=nBatch(iSyBlk)
C?    print *,' iSyBlk,iBatch',iSyBlk,iBatch
      iSymj=iSym.eq.jSym
      kSyml=kSym.eq.lSym
      iSymk=iSym.eq.kSym
      jSyml=jSym.eq.lSym
      ijSymkl=iSymj.and.kSyml
      ikSymjl=iSymk.and.jSyml
      itOrb=NTOOBS(iSym)
      jtOrb=NTOOBS(jSym)
      ktOrb=NTOOBS(kSym)
      ltOrb=NTOOBS(lSym)
C?    print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
      ijPairs=itOrb*jtOrb
      klPairs=ktOrb*ltOrb
      nInts=ijPairs*klPairs
      If ( ikSymjl ) Then
        If ( iSymj ) Then
          ijPairs=itOrb*(itOrb+1)/2
          klPairs=ijPairs
          nInts=ijPairs*(klPairs+1)/2
        Else
          ijPairs=itOrb*jtOrb
          klPairs=ijPairs
          nInts=ijPairs*(klPairs+1)/2
        End If
      Else If ( ijSymkl ) Then
        ijPairs=itOrb*(itOrb+1)/2
        klPairs=ktOrb*(ktOrb+1)/2
        nInts=ijPairs*klPairs
      End If
C?    print *,' ijPairs,klPairs',ijPairs,klPairs
      iRecOld=-1
      iInt=0
C?    print *,' Start loop over Coulomb integrals'
C?    print *,' ijkl,ij,kl,i,j,k,l'
      Do lJeppe=lOff,lOff+lOrb-1
        jMin=jOff
        If ( JLSM.ne.0 ) jMin=lJeppe
        Do jJeppe=jMin,jOff+jOrb-1
          Do kJeppe=kOff,kOff+kOrb-1
            iMin = iOff
            If(IKSM.ne.0) iMin = kJeppe
            Do iJeppe=iMin,iOff+iOrb-1
              If ( ijPerm ) Then
                ii=jJeppe
                jj=iJeppe
              Else
                ii=iJeppe
                jj=jJeppe
              End If
              If ( klPerm ) Then
                kk=lJeppe
                ll=kJeppe
              Else
                kk=kJeppe
                ll=lJeppe
              End If
              If ( ijklPerm ) Then
                i=kk
                j=ll
                k=ii
                l=jj
              Else
                i=ii
                j=jj
                k=kk
                l=ll
              End If
              If ( iSymj ) Then
                If ( i.gt.j ) Then
                  ij=j+i*(i-1)/2
                Else
                  ij=i+j*(j-1)/2
                End If
              Else
                ij=j+(i-1)*jtOrb
              End If
              ijOff=ij+(ij-1)*(ij-2)/2-1
              If ( kSyml ) Then
                If ( k.gt.l ) Then
                  kl=l+k*(k-1)/2
                Else
                  kl=k+l*(l-1)/2
                End If
              Else
                kl=l+(k-1)*ltOrb
              End If
              klOff=kl+(kl-1)*(kl-2)/2-1
              If ( ikSymjl ) Then
                If ( ij.gt.kl ) Then
                  ijkl=ij+(kl-1)*ijPairs-klOff
                Else
                  ijkl=kl+(ij-1)*klPairs-ijOff
                End If
              Else
                ijkl=ij+(kl-1)*ijPairs
              End If
C?            print '(7i6)',ijkl,ij,kl,i,j,k,l
              iRec=(ijkl-1)/lTr2Rec
              If ( iRec.eq.iRecOld ) then
                ijkl=ijkl-iRec*lTr2Rec
              Else
                iDisk=iDAdr(iBatch)
                Do iSkip=1,iRec
                End Do
                ijkl=ijkl-iRec*lTr2Rec
                iRecOld=iRec
              End If
              iInt=iInt+1
              XInt(iInt)=Scr(ijkl)
            End Do
          End Do
        End Do
      End Do
*
*     Collect Exchange terms
*
      If ( IXCHNG.ne.0 ) Then
        ijPerm=.false.
        klPerm=.false.
        ijklPerm=.false.
        If ( iSm.gt.lSm ) Then
          iSym=iSm
          jSym=lSm
        Else
          iSym=lSm
          jSym=iSm
          ijPerm =.True.
        End If
        ijBlk=jSym+iSym*(iSym-1)/2
        If ( kSm.gt.jSm ) Then
          kSym=kSm
          lSym=jSm
        Else
          kSym=jSm
          lSym=kSm
          klPerm =.True.
        End If
        klBlk=lSym+kSym*(kSym-1)/2
        If ( klBlk.gt.ijBlk ) Then
          iTemp=iSym
          iSym=kSym
          kSym=iTemp
          iTemp=jSym
          jSym=lSym
          lSym=iTemp
          iTemp=ijBlk
          ijBlk=klBlk
          klBlk=iTemp
          ijklPerm = .True.
        End If
        iSyBlk=klBlk+ijBlk*(ijBlk-1)/2
        iBatch=nBatch(iSyBlk)
        iSymj=iSym.eq.jSym
        kSyml=kSym.eq.lSym
        iSymk=iSym.eq.kSym
        jSyml=jSym.eq.lSym
        ijSymkl=iSymj.and.kSyml
        ikSymjl=iSymk.and.jSyml
        itOrb=NTOOBS(iSym)
        jtOrb=NTOOBS(jSym)
        ktOrb=NTOOBS(kSym)
        ltOrb=NTOOBS(lSym)
        ijPairs=itOrb*jtOrb
        klPairs=ktOrb*ltOrb
        nInts=ijPairs*klPairs
        If ( ikSymjl ) Then
          If ( iSymj ) Then
            ijPairs=itOrb*(itOrb+1)/2
            klPairs=ijPairs
            nInts=ijPairs*(klPairs+1)/2
          Else
            ijPairs=itOrb*jtOrb
            klPairs=ijPairs
            nInts=ijPairs*(klPairs+1)/2
          End If
        Else If ( ijSymkl ) Then
          ijPairs=itOrb*(itOrb+1)/2
          klPairs=ktOrb*(ktOrb+1)/2
          nInts=ijPairs*klPairs
        End If
        iRecOld=-1
        iInt=0
C?      print *,' iSym,jSym,kSym,lSym',iSym,jSym,kSym,lSym
C?      print *,' iSyBlk,iBatch',iSyBlk,iBatch
C?      print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
C?      print *,' ijPairs,klPairs',ijPairs,klPairs
C?      print *,' Start loop over Exchange integrals'
C?      print *,' ijkl,ij,kl,i,j,k,l'
        Do lJeppe=lOff,lOff+lOrb-1
          jMin=jOff
          If ( JLSM.ne.0 ) jMin=lJeppe
          Do jJeppe=jMin,jOff+jOrb-1
            Do kJeppe=kOff,kOff+kOrb-1
              iMin = iOff
              If(IKSM.ne.0) iMin = kJeppe
              Do iJeppe=iMin,iOff+iOrb-1
              If ( ijPerm ) Then
                ii=lJeppe
                ll=iJeppe
              Else
                ii=iJeppe
                ll=lJeppe
              End If
              If ( klPerm ) Then
                kk=jJeppe
                jj=kJeppe
              Else
                kk=kJeppe
                jj=jJeppe
              End If
              If ( ijklPerm ) Then
                i=kk
                l=jj
                k=ii
                j=ll
              Else
                i=ii
                j=jj
                k=kk
                l=ll
              End If
                If ( iSymj ) Then
                  If ( i.gt.l ) Then
                    ij=l+i*(i-1)/2
                  Else
                    ij=i+l*(l-1)/2
                  End If
                Else
                  ij=l+(i-1)*jtOrb
                End If
                ijOff=ij+(ij-1)*(ij-2)/2-1
                If ( kSyml ) Then
                  If ( k.gt.j ) Then
                    kl=j+k*(k-1)/2
                  Else
                    kl=k+j*(j-1)/2
                  End If
                Else
                  kl=j+(k-1)*ltOrb
                End If
                klOff=kl+(kl-1)*(kl-2)/2-1
                If ( ikSymjl ) Then
                  If ( ij.gt.kl ) Then
                    ijkl=ij+(kl-1)*ijPairs-klOff
                  Else
                    ijkl=kl+(ij-1)*klPairs-ijOff
                  End If
                Else
                  ijkl=ij+(kl-1)*ijPairs
                End If
C?              print '(7i6)',ijkl,ij,kl,i,j,k,l
                iRec=(ijkl-1)/lTr2Rec
                If ( iRec.eq.iRecOld ) then
                  ijkl=ijkl-iRec*lTr2Rec
                Else
                  iDisk=iDAdr(iBatch)
                  Do iSkip=1,iRec
                  End Do
                  ijkl=ijkl-iRec*lTr2Rec
                  iRecOld=iRec
                End If
                iInt=iInt+1
                XInt(iInt)=XInt(iInt)-Scr(ijkl)
              End Do
            End Do
          End Do
        End Do
      End If
*
      Return
      End
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETINT_DBG(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &                      KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &                      IRI,IXCHNG,IKSM,JLSM,ICOUL)
*
*
* Obtain block of double group integrals
*     Jeppe Olsen
*
*  New list structure incorporated
*     Timo Fleig
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. output
      DIMENSION XINT(*)
*. Scratch (should be able to hold one integral block)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
*
      INTEGER IREO(4)
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Godmorgen velkommen til GETINT_DBG '
        WRITE(6,*) ' Godmorgen velkommen til GETINT_DBG '
        WRITE(6,*) ' Godmorgen velkommen til GETINT_DBG '
        WRITE(6,*) ' Godmorgen velkommen til GETINT_DBG '
*
        WRITE(6,*) ' ITYP,ISYM,IBU',  ITYP,ISYM,IBU
        WRITE(6,*) ' JTYP,JSYM,JBU',  JTYP,JSYM,JBU
        WRITE(6,*) ' KTYP,KSYM,KBU',  KTYP,KSYM,KBU
        WRITE(6,*) ' LTYP,LSYM,LBU',  LTYP,LSYM,LBU
        WRITE(6,*)
        WRITE(6,*) ' JLSM = ', JLSM
      END IF
*
*
* The lists of integrals are
*
* 1 : (u u | u u)
*--------------------
* 2 : (u b | u u)
* 3 : (b u | u b)
* 4 : (b u | b u)
*--------------------
* 5 : (u b | u b)
* 6 : (b u | u u)
*--------------------
*
* What list of integrals should be accessed
*
      IF(IBU.NE.JBU.AND.KBU.NE.LBU) THEN
*. List 3 or list 4 or list 5
         IF(IBU.EQ.KBU) THEN
*. List 4 or 5, two possibilities : (u b ! u b) (c.c)
*                                   (b u ! b u) (c.c)
            IF(IBU.EQ.2) THEN
* ( U B ! U B) : complex conjugation
*. Output => input index
               IREO(1) = 2
               IREO(2) = 1
               IREO(3) = 3
               IREO(4) = 4
*
               SCALE = 1.0D0
C?             WRITE(6,*) ' GETINT_DBG_S 1 '
               CALL GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &              KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,3)),WORK(KPINT2D(3)),
     &              IREO,SCALE)
*
C                   GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
C    &                           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
C    &                           IXCHNG,IKSM,JLSM,ICOUL,
C    &                           XINTLIST,PINT2D,IREOINDX,SCALE)
            else if (IBU.eq.1) then
*. (B U ! B U) : complex conjugation
               IREO(1) = 1
               IREO(2) = 2
               IREO(3) = 4
               IREO(4) = 3
               IF(IRI.EQ.1) THEN
                  SCALE = 1.0D0
               ELSE
                  SCALE = -1.0D0
               END IF
C?             WRITE(6,*) ' GETINT_DBG_S 2 '
               CALL GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &              KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,3)),WORK(KPINT2D(3)),
     &              IREO,SCALE)
            END IF
*           ^ End if IBU = 2 or 1
         ELSE
*. List 4 : Two possibilities : (b u ! u b) (directly obtained)
*                               (u b ! b u) (c.c)
            IF(IBU.EQ.1) THEN
*     . Output => input index
               IREO(1) = 1
               IREO(2) = 2
               IREO(3) = 3
               IREO(4) = 4
*
               SCALE = 1.0D0
C?             WRITE(6,*) ' GETINT_DBG_S 3 '
               CALL GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &                       KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &                       IXCHNG,IKSM,JLSM,ICOUL,
     &                       WORK(KINT2D(IRI,3)),WORK(KPINT2D(3)),
     &                       IREO,SCALE)
            ELSE
*. (U B ! B U) : particle 1-2 interchange
               IREO(1) = 3
               IREO(2) = 4
               IREO(3) = 1
               IREO(4) = 2
*
               SCALE = 1.0D0
C?             WRITE(6,*) ' GETINT_DBG_S 4 '
               CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &                        KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &                        IXCHNG,IKSM,JLSM,ICOUL,
     &                        WORK(KINT2D(IRI,3)),WORK(KPINT2D(3)),
     &                        IREO,SCALE)
            END IF
*         ^ End of the two possibilities
         END IF
*       ^ End of switch list3/list4,5
      END IF
*
*
      IF(IBU.EQ.JBU.AND.KBU.EQ.LBU) THEN
*
* Derived from list 1 : ( U U ! U U )
         IF(JBU.EQ.KBU) THEN
*
* List 1 : ( U U ! U U) or ( B B ! B B)
*
            IF(IBU.EQ.2) THEN
*. ( U U ! U U ) case, simple
               IREO(1) = 1
               IREO(2) = 2
               IREO(3) = 3
               IREO(4) = 4
*
               SCALE = 1.0D0
C?             WRITE(6,*) ' GETINT_DBG_S 5 '
               CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &              KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,1)),WORK(KPINT2D(1)),
     &              IREO,SCALE)
            ELSE
* ( B B ! B B) case, use kramers symmetry and c.c
               CALL KRA_SYM(ISYM,IBU,ISYMK,IBUK)
               CALL KRA_SYM(JSYM,JBU,JSYMK,JBUK)
               CALL KRA_SYM(KSYM,KBU,KSYMK,KBUK)
               CALL KRA_SYM(LSYM,LBU,LSYMK,LBUK)
               IREO(1) = 1
               IREO(2) = 2
               IREO(3) = 3
               IREO(4) = 4
               IF(IRI.EQ.1) THEN
                  SCALE = 1.0D0
               ELSE
                  SCALE = -1.0D0
               END IF
C?             WRITE(6,*) ' GETINT_DBG_S 6 '
               CALL  GETINT_DBG_S(XINT,ITYP,ISYMK,IBUK,JTYP,JSYMK,JBUK,
     &              KTYP,KSYMK,KBUK,LTYP,LSYMK,LBUK,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,1)),WORK(KPINT2D(1)),
     &              IREO,SCALE)
            END IF
         ELSE
*. ( U U ! B B ) or ( B B ! U U )

            IF(IBU.EQ.2) THEN
* ( U U ! B1 B2 ) = ( U U ! U2 U1)
               CALL KRA_SYM(KSYM,KBU,KSYMK,KBUK)
               CALL KRA_SYM(LSYM,LBU,LSYMK,LBUK)
               IREO(1) = 1
               IREO(2) = 2
               IREO(3) = 4
               IREO(4) = 3
C?             WRITE(6,*) ' GETINT_DBG_S 7 '
               SCALE = 1.0D0
               CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &              KTYP,KSYMK,KBUK,LTYP,LSYMK,LBUK,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,1)),WORK(KPINT2D(1)),
     &              IREO,SCALE)
            ELSE IF (IBU.EQ.1) THEN
* ( B1 B2 ! U U) = (U2 U1 ! U U )
               CALL KRA_SYM(ISYM,IBU,ISYMK,IBUK)
               CALL KRA_SYM(JSYM,JBU,JSYMK,JBUK)
               IREO(1) = 2
               IREO(2) = 1
               IREO(3) = 3
               IREO(4) = 4
C?             WRITE(6,*) ' GETINT_DBG_S 8 '
               WRITE(6,*) ' IRI = ', IRI
               SCALE = 1.0D0
               WRITE(6,*) ' KINT2D(IRI,1)', KINT2D(IRI,1)
               CALL  GETINT_DBG_S(XINT,ITYP,ISYMK,IBUK,JTYP,JSYMK,JBUK,
     &              KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &              IXCHNG,IKSM,JLSM,ICOUL,
     &              WORK(KINT2D(IRI,1)),WORK(KPINT2D(1)),
     &              IREO,SCALE)
            END IF
*        ^ End of ( U U ! B B ) or ( B B ! U U ) switch
         END IF
*       ^ End of the cases for list 1,
      END IF
*     ^ End of the use of list 1
*
* List 2 :  (u b ! u u)
*
      IF((IBU.NE.JBU.AND.KBU.EQ.LBU).OR.
     &     (IBU.EQ.JBU.AND.KBU.NE.LBU)) THEN
*
* 8 cases
*         ( U B ! U U ) Simple
*         ( B U ! U U ) use c.c and Kramers symmetry
*         ( U B ! B B ) Use Kramers symmetry
*         ( B U ! B B ) use c.c
*
* The remaining four cases have different symmetries at particle 2,
* use particle 1-2 particle two symmetry to get these in place
*
         IF(IBU.EQ.2.AND.JBU.EQ.1.AND.KBU.EQ.2.AND.LBU.EQ.2) THEN
*. the simple one
            IREO(1) = 1
            IREO(2) = 2
            IREO(3) = 3
            IREO(4) = 4
            SCALE = 1.0D0
C?          WRITE(6,*) ' GETINT_DBG_S 9 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF(IBU.EQ.1.AND.JBU.EQ.2.AND.KBU.EQ.2.AND.LBU.EQ.2)THEN
* (B U ! U U ) => ( U B ! U U) by c.c
            IREO(1) = 2
            IREO(2) = 1
            IREO(3) = 4
            IREO(4) = 3
            IF(IRI.EQ.1) THEN
               SCALE = 1.0D0
            ELSE
               SCALE = -1.0D0
            END IF
C?          WRITE(6,*) ' GETINT_DBG_S 10 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF(IBU.EQ.2.AND.JBU.EQ.1.AND.KBU.EQ.1.AND.LBU.EQ.1) THEN
*. ( U  B1 ! B2 B3) => (U B1 ! U3 U2) by kramers
            IREO(1) = 1
            IREO(2) = 2
            IREO(3) = 4
            IREO(4) = 3
            CALL KRA_SYM(KSYM,KBU,KSYMK,KBUK)
            CALL KRA_SYM(LSYM,LBU,LSYMK,LBUK)
            SCALE = 1.0D0
C?          WRITE(6,*) ' GETINT_DBG_S 11 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYMK,KBUK,LTYP,LSYMK,LBUK,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF (IBU.EQ.1.AND.JBU.EQ.2.AND.KBU.EQ.1.AND.LBU.EQ.1) THEN
*. (B1 U ! B2 B3 ) => ( B1 U ! U3 U2 ) (kr) => (U B1 ! B3 B2) * (c.c)
            IREO(1) = 2
            IREO(2) = 1
            IREO(3) = 3
            IREO(4) = 4
            CALL KRA_SYM(KSYM,KBU,KSYMK,KBUK)
            CALL KRA_SYM(LSYM,LBU,LSYMK,LBUK)
            IF(IRI.EQ.1) THEN
               SCALE = 1.0D0
            ELSE
               SCALE = -1.0D0
            END IF
C?          WRITE(6,*) ' GETINT_DBG_S 12 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYMK,KBUK,LTYP,LSYMK,LBUK,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
*. and then the transpose cases
         ELSE IF(IBU.EQ.2.AND.JBU.EQ.2.AND.KBU.EQ.2.AND.LBU.EQ.1) THEN
* ( U U ! U B ) => ( U B ! U U ) by particle one -  two interchange
            IREO(1) = 3
            IREO(2) = 4
            IREO(3) = 1
            IREO(4) = 2
            SCALE = 1.0D0
C?          WRITE(6,*) ' GETINT_DBG_S 13 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF(IBU.EQ.2.AND.JBU.EQ.2.AND.KBU.EQ.1.AND.LBU.EQ.2) THEN
* (U U ! B U ) => (1-2) (B U ! U U ) => ( U B ! U U) by c.c
*
            IREO(1) = 4
            IREO(2) = 3
            IREO(3) = 2
            IREO(4) = 1
            IF(IRI.EQ.1) THEN
               SCALE = 1.0D0
            ELSE
               SCALE = -1.0D0
            END IF
C?          WRITE(6,*) ' GETINT_DBG_S 14 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYM,IBU,JTYP,JSYM,JBU,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF( IBU.EQ.1.AND.JBU.EQ.1.AND.KBU.EQ.2.AND.LBU.EQ.1) THEN
*. (B2 B3 ! U B1 ) =>(1-2) ( U  B1 ! B2 B3) => (U B1 ! U3 U2) by kramers
            CALL KRA_SYM(ISYM,IBU,ISYMK,IBUK)
            CALL KRA_SYM(JSYM,JBU,JSYMK,JBUK)
            SCALE = 1.0D0
            IREO(1) = 4
            IREO(2) = 3
            IREO(3) = 1
            IREO(4) = 2
C?          WRITE(6,*) ' GETINT_DBG_S 15 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYMK,IBUK,JTYP,JSYMK,JBUK,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         ELSE IF (IBU.EQ.1.AND.JBU.EQ.1.AND.KBU.EQ.1.AND.LBU.EQ.2) THEN
*. (B2 B3 ! B1 U) => (1-2) (B1 U ! B2 B3 ) => (Kr)
*  ( B1 U ! U3 U2 )  => (c.c) (U B1 ! U3 U2) *
            CALL KRA_SYM(ISYM,IBU,ISYMK,IBUK)
            CALL KRA_SYM(JSYM,JBU,JSYMK,JBUK)
            IREO(1) = 4
            IREO(2) = 3
            IREO(3) = 2
            IREO(4) = 1
            IF(IRI.EQ.1) THEN
               SCALE = 1.0D0
            ELSE
               SCALE = -1.0D0
            END IF
C?          WRITE(6,*) ' GETINT_DBG_S 16 '
            CALL  GETINT_DBG_S(XINT,ITYP,ISYMK,IBUK,JTYP,JSYMK,JBUK,
     &           KTYP,KSYM,KBU,LTYP,LSYM,LBU,
     &           IXCHNG,IKSM,JLSM,ICOUL,
     &           WORK(KINT2D(IRI,2)),WORK(KPINT2D(2)),
     &           IREO,SCALE)
         END IF
*      ^ End of the eight choices for List 2
      END IF
*     ^ End of list 2
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETINT_DBG_S(XINT,ITP,ISM,IBU,JTP,JSM,JBU,
     &                        KTP,KSM,KBU,LTP,LSM,LBU,
     &                        IXCHNG,IKSM,JLSM,ICOUL,
     &                        INTLST,IJKLOF,IREOINDX,SCALE)
*
* Obtain integrals
*
*     ICOUL = 0 :
*                  XINT(IK,JL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 1 :
*                  XINT(IJ,KL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
* Storing for ICOUL = 1 not working if IKSM or JLSM .ne. 0
*
* Version for integrals stored in INTLST
*
* Double group integral version. At the moment no
* permutation symmetry
*
* INDXREO gives mapping between input and output indeces
* output  index => input  index
*
* Input index : Index in complete integral list
* Output index : Index in outblock block
*
* all symmetries etc refers to data for the output block
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "lucinp.inc"
*. Reordering of indeces
      INTEGER IREOINDX(4)
*. INDXREO(I) gives the
*. Integral list
      REAL*8 INTLST(*)
      Dimension IJKLof(NsmOB,NsmOb,NsmOB)
*.Output
      DIMENSION XINT(*)
*. A bit of local scratch
      INTEGER INSM(4),IOUTSM(4),ISCR(4)
*
      NTEST = 000
*
      if (NTEST.ge.10) then
         WRITE(6,*) ' Godmorgen dette er GETINT_DBG_S '
         WRITE(6,*) ' IREOINDX : ', (IREOINDX(I),I=1,4)
         WRITE(6,*) ' ISM, JSM, KSM, LSM',ISM, JSM, KSM, LSM
         WRITE(6,*) ' IBU, JBU, KBU, LBU',IBU, JBU, KBU, LBU
         WRITE(6,*) ' SCALE = ', SCALE
      end if
*
      IF(IBU.EQ.2) THEN
         IORB = NOBPTS(ITP,ISM)
         ITORB = NTOOBS(ISM)
      ELSE
         IORB = NOBPTS2(ITP,ISM)
         ITORB = NTOOBS2(ISM)
      END IF
*
      IF(JBU.EQ.2) THEN
         JORB = NOBPTS(JTP,JSM)
         JTORB = NTOOBS(JSM)
      ELSE
         JORB = NOBPTS2(JTP,JSM)
         JTORB = NTOOBS2(JSM)
      END IF
*
      IF(KBU.EQ.2) THEN
         KORB = NOBPTS(KTP,KSM)
         KTORB = NTOOBS(KSM)
      ELSE
         KORB = NOBPTS2(KTP,KSM)
         KTORB = NTOOBS2(KSM)
      END IF
*
      IF(LBU.EQ.2) THEN
         LORB = NOBPTS(LTP,LSM)
         LTORB = NTOOBS(LSM)
      ELSE
         LORB = NOBPTS2(LTP,LSM)
         LTORB = NTOOBS2(LSM)
      END IF
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*)
     &   ' GETINT_DGB_S: IORB,JORB,KORB,LORB',IORB,JORB,KORB,LORB
         WRITE(6,*)
     &   '  ITORB, JTORB, KTORB, LTORB' ,
     &        ITORB, JTORB, KTORB, LTORB
      END IF
*. Offsets relative to start of all orbitals, symmetry ordered
      IF(IBU.EQ.2) THEN
         IOFF = IBSO(ISM)
         IOFF0 = IOFF
         DO IITP = 1, ITP -1
            IOFF = IOFF + NOBPTS(IITP,ISM)
         END DO
      ELSE
         IOFF = IBSO2(ISM)
         IOFF0 = IOFF
         DO IITP = 1, ITP -1
            IOFF = IOFF + NOBPTS2(IITP,ISM)
         END DO
      END IF
*
      IF(JBU.EQ.2) THEN
         JOFF = IBSO(JSM)
         JOFF0 = JOFF
         DO JJTP = 1, JTP -1
            JOFF = JOFF + NOBPTS(JJTP,JSM)
         END DO
      ELSE
         JOFF = IBSO2(JSM)
         JOFF0 = JOFF
         DO JJTP = 1, JTP -1
            JOFF = JOFF + NOBPTS2(JJTP,JSM)
         END DO
      END IF
*
      IF(KBU.EQ.2) THEN
         KOFF = IBSO(KSM)
         KOFF0 = KOFF
         DO KKTP = 1, KTP -1
            KOFF = KOFF + NOBPTS(KKTP,KSM)
         END DO
      ELSE
         KOFF = IBSO2(KSM)
         KOFF0 = KOFF
         DO KKTP = 1, KTP -1
            KOFF = KOFF + NOBPTS2(KKTP,KSM)
         END DO
      END IF
*
      IF(LBU.EQ.2) THEN
         LOFF = IBSO(LSM)
         LOFF0 = LOFF
         DO LLTP = 1, LTP -1
            LOFF = LOFF + NOBPTS(LLTP,LSM)
         END DO
      ELSE
         LOFF = IBSO2(LSM)
         LOFF0 = LOFF
         DO LLTP = 1, LTP -1
            LOFF = LOFF + NOBPTS2(LLTP,LSM)
         END DO
      END IF
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*)
     &  ' GETINT_DBG_S: IOFF,JOFF,KOFF,LOFF',
     &                  IOFF,JOFF,KOFF,LOFF
         WRITE(6,*)
     &  ' GETINT_DBG_S: IOFF0,JOFF0,KOFF0,LOFF0',
     &                  IOFF0,JOFF0,KOFF0,LOFF0
      END IF
*
*. Symmetry of input block
      ISCR(IREOINDX(1)) = ISM
      ISCR(IREOINDX(2)) = JSM
      ISCR(IREOINDX(3)) = KSM
      ISCR(IREOINDX(4)) = LSM
*
      ISMI = ISCR(1)
      JSMI = ISCR(2)
      KSMI = ISCR(3)
      LSMI = ISCR(4)
*
      if (NTEST.ge.100) then
        WRITE(6,*) ' ISMI = ', ISMI
        WRITE(6,*) ' JSMI = ', JSMI
        WRITE(6,*) ' KSMI = ', KSMI
        WRITE(6,*) ' LSMI = ', LSMI
      end if
*
*. Input dimensions
      ISCR(IREOINDX(1)) = ITORB
      ISCR(IREOINDX(2)) = JTORB
      ISCR(IREOINDX(3)) = KTORB
      ISCR(IREOINDX(4)) = LTORB
*
      NI_IN = ISCR(1)
      NJ_IN = ISCR(2)
      NK_IN = ISCR(3)
      NL_IN = ISCR(4)
*
      if (NTEST.ge.100) then
        WRITE(6,*) ' NI_IN = ', NI_IN
        WRITE(6,*) ' NJ_IN = ', NJ_IN
        WRITE(6,*) ' NK_IN = ', NK_IN
        WRITE(6,*) ' NL_IN = ', NL_IN
      end if
*
*
*     Collect Coulomb terms
*
*
      if (NTEST.ge.100) WRITE(6,*) ' Before IBLOFF '
*
      IBLOFF = IJKLOF(ISMI,JSMI,KSMI)
*
      if (NTEST.ge.100) then
         WRITE(6,*) ' IBLOFF = ', IBLOFF
         WRITE(6,*) ' Block of two-electron integrals in input lists '
         CALL WRTMAT(INTLST(IBLOFF),NI_IN*NJ_IN,NK_IN*NL_IN,
     &               NI_IN*NJ_IN,NK_IN*NL_IN )
      end if
*
      IINT = 0
*. Loop over indeces of output blocks
      Do L=LOFF,LOFF+LORB-1
         if (NTEST.ge.100) write(6,*) 'Doing L= ',L
         JMIN = JOFF
         If (JLSM.ne.0) JMIN = L
         Do J=JMIN,JOFF+JORB-1
            if (NTEST.ge.100) write(6,*) 'Doing J= ',J
*
            Do K=KOFF,KOFF+KORB-1
               if (NTEST.ge.100) write(6,*) 'Doing K= ',K
               IMIN = IOFF
               If (IKSM.ne.0) IMIN = K
               IF (ICOUL.EQ.1) THEN
*. Address before integral (1,j!k,l)
                IINT = (L-LOFF)*JORB*KORB*IORB + (K-KOFF)*JORB*IORB
     &                  + (J-JOFF)*IORB
               END IF
*
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) > (ksm,lsm)
               Do I=IMIN,IOFF+IORB-1
                  if (NTEST.ge.100) write(6,*) 'Doing I= ',I
*. 4 output indeces are now known, corresponding input indeces
                  ISCR(IREOINDX(1)) = I - IOFF0 + 1
                  ISCR(IREOINDX(2)) = J - JOFF0 + 1
                  ISCR(IREOINDX(3)) = K - KOFF0 + 1
                  ISCR(IREOINDX(4)) = L - LOFF0 + 1
                  I_IN = ISCR(1)
                  J_IN = ISCR(2)
                  K_IN = ISCR(3)
                  L_IN = ISCR(4)
*. (I J K L ) in input list
                  IJKL_IN = (L_IN-1)*NK_IN*NJ_IN*NI_IN
     &                 +    (K_IN-1)      *NJ_IN*NI_IN
     &                 +    (J_IN-1)            *NI_IN
     &                 +     I_IN
*
                  IJKL = IBLOFF-1+ IJKL_IN
                  IINT = IINT + 1
*
                  if (NTEST.ge.100) then
                     WRITE(6,*) ' IJKL_IN' , IJKL_IN
                     WRITE(6,*) ' IINT and IJKL ',IINT, IJKL
                  end if
*
                  XINT(IINT) = SCALE*INTLST(IJKL)
*
                  if (NTEST.ge.100) then
                     write(6,*) 'Processed INTLST(',IJKL,')
     &                          to XINT(',IINT,')'
                  end if
               End Do
*
            End Do
         End Do
      End Do
*
*     Collect Exchange terms
*
      If ( IXCHNG.ne.0 ) Then
*
*
         IBLOFF = IJKLOF(ISMI,LSMI,KSMI)
*
         IINT = 0
         Do L=LOFF,LOFF+LORB-1
            JMIN = JOFF
            If (JLSM.ne.0) JMIN = L
            Do J=JMIN,JOFF+JORB-1
               Do K=KOFF,KOFF+KORB-1
*
                  IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                   IINT = (L-LOFF)*JORB*KORB*IORB + (K-KOFF)*JORB*IORB
     &                     + (J-JOFF)*IORB
                  END IF
*
                  If(IKSM.ne.0) IMIN = K
                  Do I=IMIN,IOFF+IORB-1
*. 4 output indeces are now known, corresponding input indeces
                     ISCR(IREOINDX(1)) = I - IOFF0 + 1
                     ISCR(IREOINDX(2)) = J - JOFF0 + 1
                     ISCR(IREOINDX(3)) = K - KOFF0 + 1
                     ISCR(IREOINDX(4)) = L - LOFF0 + 1
                     I_IN = ISCR(1)
                     J_IN = ISCR(2)
                     K_IN = ISCR(3)
                     L_IN = ISCR(4)
*. (I L !K J) in output list
                     ILKJ0 = (J_IN-1)*NK_IN*NL_IN*NI_IN
     &                     + (K_IN-1)      *NL_IN*NI_IN
     &                     + (L_IN-1)            *NI_IN
     &                     +  I_IN
                     ILKJ = ILKJ0 + IBLOFF - 1
                     IINT = IINT + 1
                     XINT(IINT)=XINT(IINT)-SCALE*INTLST(ILKJ)
                  End Do
*
               End Do
            End Do
         End Do
      End If
*
      IF(NTEST.GE.100) THEN
        IF(ICOUL.EQ.1) THEN
          WRITE(6,*) ' output block in output format '
          CALL WRTMAT(XINT,IORB*JORB,KORB*LORB,IORB*JORB,KORB*LORB)
        END IF
      END IF
*
      Return
      End
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETINT_REL(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                      IXCHNG,IKSM,JLSM,ICOUL)
*
* Interface to DIRAC relativistic integrals
* (in spinfree formalism) implemented
*
*---------------------------------------------
* Timo Fleig, October 2000
*---------------------------------------------
*
* Outer routine for accessing integral block
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
*./CRUN/ : INTIMP used
#include "crun.inc"
*./ORBINP/  : NOBPTS used
*.Memory
#include "wrkspc.inc"
*
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cgas.inc"

C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1
#include "glbbas.inc"
*
      CALL QENTER('GETIN')
      NTEST = 000
*
      IF(INTIMP.EQ.1.OR.INTIMP.EQ.5.or.INTIMP.eq.6) THEN
*. Interface to MOTRA or SIRIUS integrals
*  or to DIRAC (spinfree) integrals
        IF(INCORE.EQ.0) THEN
        CALL GETINM(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &              IXCHNG,IKSM,JLSM,NORBMOL,NIRR_PN)
        ELSE IF (INCORE.NE.0) THEN
*. Old
COLD      CALL GETINC(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
COLD &                  IXCHNG,IKSM,JLSM,WORK(KINT2),
COLD &                  WORK(KPINT2),NSMOB)
*
COLD       WRITE(6,*) ' Integrals with old fetcher : '
COLD    NI = NOBPTS(ITP,ISM)
COLD    NK = NOBPTS(KTP,KSM)
COLD    IF(IKSM.EQ.0) THEN
COLD      NIK = NI * NK
COLD    ELSE
COLD      NIK = NI*(NI+1)/2
COLD    END IF
COLD    NJ = NOBPTS(JTP,JSM)
COLD    NL = NOBPTS(LTP,LSM)
COLD    IF(JLSM.EQ.0) THEN
COLD      NJL = NJ * NL
COLD    ELSE
COLD      NJL = NJ*(NJ+1)/2
COLD    END IF
COLD    WRITE(6,*) ' 2 electron integral block for TS blocks '
COLD    WRITE(6,*) ' Ixchng :', IXCHNG
COLD    WRITE(6,*) ' After GETINC '
COLD    WRITE(6,'(1H ,4(A,I2,A,I2,A))')
COLD &  '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
COLD &  '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
COLD    CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
C
*. New
        write(6,*) 'GETINT_REL: call GETINCN with KINT2'
        CALL GETINCN(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &               IXCHNG,IKSM,JLSM,WORK(KINT2),
     &               WORK(KPINT2),NSMSH,WORK(KINH1),ICOUL)
        END IF
      else if (INTIMP.EQ.2) THEN
C. Interface to LUCAS integrals
        CALL GETINL(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM)
      else if (INTIMP.eq.8) then
*
* Molecule Sweden interface
*
        call GETINCN(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &               IXCHNG,IKSM,JLSM,WORK(KINT2),
     &               WORK(KPINT2),NSMSH,WORK(KINH1),ICOUL)
      END IF
*
      IF(NTEST.NE.0) THEN
        NI = NGSSH(ISM,ITP)
        NK = NGSSH(KSM,KTP)
        IF(IKSM.EQ.0) THEN
          NIK = NI * NK
        ELSE
          NIK = NI*(NI+1)/2
        END IF
        NJ = NGSSH(JSM,JTP)
        NL = NGSSH(LSM,LTP)
        IF(JLSM.EQ.0) THEN
          NJL = NJ * NL
        ELSE
          NJL = NJ*(NJ+1)/2
        END IF
        WRITE(6,*) ' 2 electron integral block for TS blocks '
        WRITE(6,*) ' Ixchng :', IXCHNG
        WRITE(6,*) ' After GETINC '
        WRITE(6,'(1H ,4(A,I2,A,I2,A))')
     &  '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &  '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
        CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
      END IF
*
      CALL QEXIT('GETIN')
C     STOP ' Jeppe forced me to stop in GETINT '
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Fetch integral from integral list according to SIGDEN
* excitation class formalism - driver routine
*
      subroutine getint_sigden(YINT,INDEX,IRI,ITOTDIM,
     &                         NASHT,F1,F2,IPRNT)
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "crun.inc"
#include "symm.inc"
*
      logical NZ3FORM
      dimension INDEX(4)
*
* New (memory saving) (NZ,3) integral interface (if true)
      NZ3FORM = .TRUE.
*
      if (INTIMP.eq.6) then
* MOLFDIR format Kramers reduced
        if (IRI.eq.1) then
          IRIOFF1 = 0
          IRIOFF2 = 0
        else if (IRI.eq.2) then
          IRIOFF1 = (2*ITOTDIM)**2
          IRIOFF2 = (2*ITOTDIM)**4
        end if
*
C
        call getind_sigden(YINT,INDEX,IRI,ITOTDIM,
     &                     WORK(KNUMINT),WORK(KKR),WORK(KINDK),
     &                     WORK(KINDL),WORK(KINVREO),WORK(KIOFF12),
     &                     WORK(KDITWOE+IRIOFF2),WORK(KDIONEE+IRIOFF1),
     &                     IPRNT)
      else if (INTIMP.eq.7) then
        if (.not.NZ3FORM) then
*
* Define the KR array (which is usually read from the MOLFDIR file
* in standalone calculations!)
          call kkr_set(WORK(KKR),ITOTDIM,SYDI,NIRR_PN,INVSYM)
*
* MOLFDIR format full spinor list
          call getind_molf(YINT,INDEX,IRI,ITOTDIM,WORK(KINVREO),
     &                     WORK(KKR),NASHT,F1,F2,IPRNT)
*
        else if (NZ3FORM) then
!
!         this is the default integral scheme of today's LUCIAREL.
          nct = 0
          do ic = 1,4,1
            if(INDEX(ic).ne.0) nct = nct + 1
          end do
          if(nct .eq. 4)then
!
!           2e- ints in DIRAC (NZ,3) format
            call kkr_set(WORK(KKR),ITOTDIM,SYDI,NIRR_PN,INVSYM)
            NNASHX = ((NASHT/2) * ((NASHT/2)+1)) / 2
            call getind_nz3(YINT,INDEX,IRI,ITOTDIM,WORK(KINVREO),
     &                      WORK(KKR),NASHT,NNASHX,F2,IPRNT)
          else
!
!           1e- ints
            call getint_nts(YINT,INDEX,IRI,ITOTDIM,F1,IPRNT)
          end if
        end if
      end if
*
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETINTN(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,
     &                   LTP,LSMGET,
     &                   IXCHNG,IKSMGET,JLSMGET,ICOUL,IRI)
*
* Interface to DIRAC relativistic integrals
* (in spinfree and full spin-dependent formalism) implemented
*
*---------------------------------------------
* Timo Fleig, October 2000
*---------------------------------------------
*
* Outer routine for accessing integral block
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cgas.inc"
#include "glbbas.inc"
#include "symm.inc"
*
      CALL QENTER('GETIN')
      NTEST = 000
*
      if (INTIMP.eq.1.or.INTIMP.eq.5.or.INTIMP.eq.6) then
*. Interface to MOTRA or SIRIUS integrals
*  or to DIRAC (spinfree) integrals
        if (INCORE.eq.0) then
          if (INTIMP.eq.6) then
            write(6,*) 'Out of core integral handling'
            write(6,*) 'not implemented yet.'
            call abend('Quitting.')
          end if
          call getinm(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,LTP,LSMGET,
     &                IXCHNG,IKSMGET,JLSMGET,SYDI,NIRR_PN)
        else if (INCORE.ne.0) then
* Integrals are held in core
          if (INTIMP.eq.6) then
* Read two-electron integrals from DIRAC file MDCINT
            call getind(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,
     &                  LTP,LSMGET,
     &                  WORK(KINVREO),WORK(KKR),WORK(KNUMINT),
     &                  WORK(KINDK),WORK(KINDL),WORK(KIOFF12),
     &                  WORK(KDITWOE+(IRI-1)*IMOFFSET))
          else
*. Interface to MOTRA or SIRIUS integrals
            write(6,*) 'call to getincn with KINT2_2'
            call getincn(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,
     &                   LTP,LSMGET,
     &                   IXCHNG,IKSMGET,JLSMGET,WORK(KINT2),
     &                   WORK(KPINT2),NSMSH,WORK(KINH1),ICOUL)
          end if
        end if
*
      else if (INTIMP.EQ.2) THEN
C. Interface to LUCAS integrals
        CALL GETINL(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,LTP,LSMGET,
     &                  IXCHNG,IKSMGET,JLSMGET)
*
      else if (INTIMP.eq.8) then
*
* Molecule Sweden interface
*
            write(6,*) 'call to getincn with KINT2_3'
        call GETINCN(XINT,ITP,ISMGET,JTP,JSMGET,KTP,KSMGET,LTP,LSMGET,
     &               IXCHNG,IKSMGET,JLSMGET,WORK(KINT2),
     &               WORK(KPINT2),NSMSH,WORK(KINH1),ICOUL)
      END IF
*
      IF(NTEST.NE.0) THEN
        NI = NGSSH(ISMGET,ITP)
        NK = NGSSH(KSMGET,KTP)
        IF(IKSMGET.EQ.0) THEN
          NIK = NI * NK
        ELSE
          NIK = NI*(NI+1)/2
        END IF
        NJ = NGSSH(JSMGET,JTP)
        NL = NGSSH(LSMGET,LTP)
        IF(JLSMGET.EQ.0) THEN
          NJL = NJ * NL
        ELSE
          NJL = NJ*(NJ+1)/2
        END IF
        WRITE(6,*) ' 2 electron integral block for TS blocks '
        WRITE(6,*) ' Ixchng :', IXCHNG
        WRITE(6,*) ' After GETINC '
        WRITE(6,'(1H ,4(A,I2,A,I2,A))')
     &  '(',ITP,',',ISMGET,')','(',JTP,',',JSMGET,')',
     &  '(',KTP,',',KSMGET,')','(',LTP,',',LSMGET,')'
        CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
      END IF
*
      CALL QEXIT('GETIN')
*
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_PRPOPSYM_KRCI(WORK,LWORK,NTEST)
C***********************************************************************
C
C     calculate operator symmetry of all property operators NPROP_KRCI 
C     in list LPROP_KRCI in LUCIAREL sub-double group notation.
C     
C     symmetry of property operator IVOP according to
C      _           _              _      
C     | (IVOP) =  | (IBSP_KRCI) x | (IASP_KRCI)⁺
C
C     final symmetry stored in ISYMOPPRP_KRCI (krciprop.h).
C     
C     symmetry setup in LUCIAREL
C     ==========================
C
C             boson               fermion
C             =====               =======
C
C     C2   :  a , b               1E,  2E
C     C2h  :  ag, bg, au, bu      1Eg, 2Eg, 1Eu, 2Eu  
C
C     Written by S. Knecht - Sep 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "consts.h"
C
#include "dcbgen.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "krciprop.h"
#include "krmcluci_inf.h"
#include "mxpdim.inc"
#include "symm.inc"
      DIMENSION WORK(*)
      CHARACTER REPVOP_KRCI*2, FILELAB*8
C
#include "memint.h"
      MZ = MIN(NZ,2)
      CALL QENTER('GET_PRPOPSYM_KRCI')
C
C     allocate for VOP matrix
C
      CALL MEMGET('REAL',KVOP,4*N2ASHX*MZ,WORK,KFREE,LFREE)
C
      IASP_DC       = 0
      IBSP_DC       = 0
      IASP_KRCI_SYM = 0
      IBSP_KRCI_SYM = 0
C
C     loop over list of operators  
C
      DO 20 I = 1, NPROP_KRCI
C
        CALL DZERO(WORK(KVOP),  4*N2ASHX*MZ)
C
C       read VOP matrix from file
C
        INDXPR  = LPROP_KRCI(I)
        FILELAB = PRPNAM(INDXPR)(1:8)
C
        CALL KRCI_PRPFILE(LUKRM3,FILELAB,WORK(KVOP),
     &                    4*N2ASHX*MZ,1)
C   
C       find largest absolute element in VOP matrix
C
        INDEX_VOP   = IDAMAX(4*N2ASHX*MZ,WORK(KVOP),1)
        ABS_MAX_ELM = DABS(WORK(KVOP+INDEX_VOP-1))
C       ... if largest absolute element is zero, 
C           nothing to do for this operator
        IF( ABS_MAX_ELM .eq. D0 ) GOTO 20
C
        IF( NTEST .gt. 10 ) THEN
           WRITE(LUPRI,'(/A,A8)') 
     &     '  absolute largest element for operator ',FILELAB
           WRITE(LUPRI,*) ' found to be:', ABS_MAX_ELM
        END IF
C
C       which two spinors (in LUCIAREL symmetry notation) 
C       IASP_KRCI_SYM and IBSP_KRCI_SYM does this operator connect?
C
        CALL GET_IAIB_DIRAC(IASP_KRCI_SYM,IBSP_KRCI_SYM,ABS_MAX_ELM,
     &                      WORK(KVOP),2*NASHT,2*NASHT,2*NASHT,2*NASHT,
     &                      MZ,IPQTOQ(1,0),NTEST)
C            GET_IAIB_DIRAC(IASP_KRCI_SYM,IBSP_KRCI_SYM,ELM_ABSMAX,
C    &                      AMAT,NROW,NCOL,LRQ,LCQ,JZ,IQP,NTEST)
                  
C
C       compute symmetry of operator IVOP and store in ISYMOPPRP_KRCI
C
        ISYMOPPRP_KRCI(I)=IDBGMULT(INVELM(IASP_KRCI_SYM),IBSP_KRCI_SYM)
 20   CONTINUE
C
C     print section
C
      IF( NTEST .ge. 0 )THEN
         WRITE (LUPRI,'(/A,/A)') 
     &   '  Determination of operator representation in KRCI module',
     &   '  _______________________________________________________'
     
         WRITE (LUPRI,'(/A/A)')
     &   '  operator     boson rep in DIRAC     boson rep in KR-CI',
     &   '  ________     __________________     __________________'
C            ZDIPLEN             B1u                     Au
         DO I = 1, NPROP_KRCI
           INDXPR  = LPROP_KRCI(I)
           ISYMVOP = IPRPSYM(INDXPR)
           ISYMVOP_KRCI = ISYMOPPRP_KRCI(I)
           IF( DBGSYM(1:6) .eq. 'Cs    ' .or. 
     &         DBGSYM(1:6) .eq. 'C2    ')THEN
               IF( ISYMVOP_KRCI .eq. 1 )THEN
                  REPVOP_KRCI = 'A '
               ELSE
                  REPVOP_KRCI = 'B '
               END IF
           ELSE IF( DBGSYM(1:6) .eq. 'C2h   ' )THEN
               IF( ISYMVOP_KRCI .eq. 1 )THEN
                  REPVOP_KRCI = 'Ag'
               ELSE IF( ISYMVOP_KRCI .eq. 2 )THEN
                  REPVOP_KRCI = 'Bg'
               ELSE IF( ISYMVOP_KRCI .eq. 3 )THEN
                  REPVOP_KRCI = 'Au'
               ELSE 
                  REPVOP_KRCI = 'Bu'
               END IF
           END IF
           WRITE(LUPRI,'(A11,A15,A21)') 
     &           PRPNAM(INDXPR)(1:8), REP(ISYMVOP-1), REPVOP_KRCI
         END DO
      END IF
C
C     release memory
      CALL MEMREL('KRCIGETPRP',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('GET_PRPOPSYM_KRCI')
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      Subroutine Getobs2(ECORE,NAOS_ENV,NMOS_ENV)
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Open one-electron integral file
*.    Read in orbital information  and core energy
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*
      Implicit Real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "SysDef.inc"
*
      Parameter( mxOrb = 510 )
      Parameter( LuOne = 14  )
*
      DIMENSION NAOS_ENV(*),NMOS_ENV(*)
*
      COMMON/MOLOBSR/
     & IOList(64),iToc(64),nBas(8),nOrb(8),nFro(8),nDel(8),
     & Nsym
      Character*4 Name(2,mxOrb)
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*     open the transformed one-electron integral file                  *
*----------------------------------------------------------------------*
      LuX = LuOne
*----------------------------------------------------------------------*
*     Set up the scatter/gather list and read table of contents        *
*----------------------------------------------------------------------*
      iDisk=0
      Call GSList(IOLIST,8,iToc,64,ECORE,RtoI,nSym,1,
     &            nBas,8,nOrb,8,nFro,8,nDel,8,Name,2*4*mxOrb/ItoB)
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*
* Copy to external arrays
*
      CALL ICOPVE(NORB,NMOS_ENV,NSYM)
      CALL ICOPVE(NBAS,NAOS_ENV,NSYM)
*
      Return
      End
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION GMIJKL_REL(IORB,JORB,KORB,LORB,INTLST,IJKLOF)
*
* Obtain integral (IORB JORB ! KORB LORB) MOLCAS version
* Integrals assumed in core
*
* Version for integrals stored in INTLST
*
      IMPLICIT REAL*8(A-H,O-Z)
*.ORBINP
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB),
C    &              NOBPTS(6+MXPR4T,MXPOBS),IOBPTS(6+MXPR4T,MXPOBS),
C    &              ITOOBS(MXPOBS)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
*. Integral list
      Real * 8 Intlst(*)
      Dimension IJKLOF(NsmOB,NsmOb,NsmOB)
      Logical iSymj,kSyml,ISYMK,JSYML,ijSymkl,IKSYMJL
      Logical ijklPerm
*.
      NTEST = 000
*
*. The orbital list corresponds to type ordered indeces, reform to
*. symmetry ordering
*
      IABS = IREOTS(IORB)
      ISM = ISMFTO(IORB)
      IOFF = IBSO(ISM)
*
      JABS = IREOTS(JORB)
      JSM = ISMFTO(JORB)
      JOFF = IBSO(JSM)
*
      KABS = IREOTS(KORB)
      KSM = ISMFTO(KORB)
      KOFF = IBSO(KSM)
*
      LABS = IREOTS(LORB)
      LSM = ISMFTO(LORB)
      LOFF = IBSO(LSM)
*
      If( Ntest.ge. 100) THEN
        write(6,*) ' GMIJKL at your service '
        WRITE(6,*) ' IORB IABS ISM IOFF ',IORB,IABS,ISM,IOFF
        WRITE(6,*) ' JORB JABS JSM JOFF ',JORB,JABS,JSM,JOFF
        WRITE(6,*) ' KORB KABS KSM KOFF ',KORB,KABS,KSM,KOFF
        WRITE(6,*) ' LORB LABS LSM LOFF ',LORB,LABS,LSM,LOFF
      END IF
*
      If ( jSm.gt.iSm .or. ( iSm.eq.jSm .and. JABS.gt.IABS)) Then
        iSym=jSm
        jSym=iSm
        I = JABS - JOFF + 1
        J = IABS - IOFF + 1
      Else
        iSym=iSm
        jSym=jSm
        I = IABS - IOFF + 1
        J = JABS - JOFF + 1
      End If
      ijBlk=jSym+iSym*(iSym-1)/2
      If ( lSm.gt.kSm  .or. ( kSm.eq.lSm .and. LABS.gt.KABS)) Then
        kSym=lSm
        lSym=kSm
        K = LABS -LOFF + 1
        L = KABS - KOFF + 1
      Else
        kSym=kSm
        lSym=lSm
        K = KABS - KOFF + 1
        L = LABS -LOFF + 1
      End If
      klBlk=lSym+kSym*(kSym-1)/2
*
      ijklPerm=.false.
      If ( klBlk.gt.ijBlk ) Then
        iTemp=iSym
        iSym=kSym
        kSym=iTemp
        iTemp=jSym
        jSym=lSym
        lSym=iTemp
        iTemp=ijBlk
        ijBlk=klBlk
        klBlk=iTemp
        ijklPerm=.true.
*
        iTemp = i
        i = k
        k = itemp
        iTemp = j
        j = l
        l = iTemp
      End If
      If(Ntest .ge. 100 ) then
        write(6,*) ' i j k l ',i,j,k,l
        write(6,*) ' Isym,Jsym,Ksym,Lsym',Isym,Jsym,Ksym,Lsym
      End if
*
*  Define offset for given symmetry block
      IBLoff = IJKLof(Isym,Jsym,Ksym)
      If(ntest .ge. 100 )
     &WRITE(6,*) ' IBLoff Isym Jsym Ksym ', IBLoff,ISym,Jsym,Ksym
      iSymj=iSym.eq.jSym
      kSyml=kSym.eq.lSym
      iSymk=iSym.eq.kSym
      jSyml=jSym.eq.lSym
      ikSymjl=iSymk.and.jSyml
      ijSymkl=iSymj.and.kSyml
*
      itOrb=NTOOBS(iSym)
      jtOrb=NTOOBS(jSym)
      ktOrb=NTOOBS(kSym)
      ltOrb=NTOOBS(lSym)
C?    print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
      If ( iSymj ) Then
        ijPairs=itOrb*(itOrb+1)/2
        ij=j+i*(i-1)/2
      Else
        ijPairs=itOrb*jtOrb
        ij=j + (i-1)*jtOrb
      End if
*
      IF(KSYML ) THEN
        klPairs=ktOrb*(ktOrb+1)/2
        kl=l+k*(k-1)/2
      ELSE
        klPairs=ktOrb*ltOrb
        kl=l+(k-1)*ltOrb
      End If
C?    print *,' ijPairs,klPairs',ijPairs,klPairs
*
      If ( ikSymjl ) Then
        If ( ij.gt.kl ) Then
          klOff=kl+(kl-1)*(kl-2)/2-1
          ijkl=ij+(kl-1)*ijPairs-klOff
        Else
          ijOff=ij+(ij-1)*(ij-2)/2-1
          ijkl=kl+(ij-1)*klPairs-ijOff
        End If
      Else
        ijkl=ij+(kl-1)*ijPairs
      End If
      If( ntest .ge. 100 )
     & write(6,*) ' ijkl ', ijkl
*
      GMIJKL_REL = Intlst(iblOff-1+ijkl)
      If( ntest .ge. 100 )
     & write(6,*) ' GMIJKL_REL ', GMIJKL_REL
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION GTIJKL_REL(I,J,K,L)
*
* Obtain  integral (I J ! K L )
* where I,J,K and l refers to active orbitals in
* Type ordering
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"

*.GLobal pointers
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "crun.inc"
      IF(INTIMP .EQ. 2 ) THEN
*. LUCAS ordering
        I12S = 0
        I34S = 0
        I1234S = 1
CINA    GTIJKL_REL = GIJKLL(IREOTS(1+NINOB),WORK(KPINT2),WORK(KLSM2),
CINA &           WORK(KINT2),
CINA &           ISMFTO,IBSO,NACOB,NSMOB,NOCOBS,I,J,K,L)
        write(6,*) ' function GTIJKL_REL: call to GIJKLL with KINT2' 
        GTIJKL_REL = GIJKLL(IREOTS(1),WORK(KPINT2),WORK(KLSM2),
     &           WORK(KINT2),
     &           ISMFTO,IBSO,NACOB,NSMOB,NOCOBS,I,J,K,L)
       ELSE IF (INTIMP.EQ.1.OR.INTIMP.EQ.5) THEN
*. MOLCAS OR SIRIUS IMPORT ( I hope integrals are in core !! )
        write(6,*) ' function GTIJKL_REL 2: call to GIJKLL with KINT2' 
          GTIJKL_REL = GMIJKL_REL(I,J,K,L,WORK(KINT2),WORK(KPINT2))
       END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Routine connects LUCIAREL info with MOLECULE SWEDEN info
*  with respect to basis functions and related things.
*  (at least is planned to do so ...)

      subroutine gtmoswe(XINT_L,NTOTINT)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "wrkspc.inc"
*. Input
#include "cintfo.inc"
      dimension NORBMOL(MXNSPIR),NBASMOL(MXNSPIR)
      dimension NFRO(MXNSPIR)
*. Output
      DIMENSION XINT_L(*)
*
      IDUMMY = 1
      call memmar(KDUMMY,IDUMMY,'MARK  ', IDUMMY,'GTMOSW')
      call memmar(KLSWINT,NINT2,'ADDL  ',2,'SWEINT')
      write(6,*) 'SWEINT length NINT2  ',NINT2
*
*  Call integral fetching routine ...
*
      do ISYM = 1,MXNSPIR,1
         NORBCT = 0
         NBASCT = 0
         do IGAS = 1,NGAS,1
            NORBCT = NORBCT + NGSSH(ISYM,IGAS)
            NBASCT = NBASCT + NGSSH(ISYM,IGAS)
         end do
         NORBMOL(ISYM) = NORBCT
         NBASMOL(ISYM) = NBASCT
         NFRO(ISYM) = 0
      end do
*
*  The following info is superfluous, as it is provided by LUCIAREL's
*  input. Still, we retain the possibility of using a frozen core in
*  rdmole2 by keeping the arrays NBASMOL, NORBMOL, and NFRO for
*  future purposes.
*
C     LUONEM =
C     REWIND LUONEM
C     CALL SEARCH(LUONEM,'INDATACI')
C     READ (LUONEM) NIRR_PN,(NBASMOL(I),I=1,NIRR_PN),
C    &              (NORBMOL(I),I=1,NIRR_PN),
C    &              PNUCCO,(NFRO(I),I=1,NIRR_PN)
C     IF(IPRINT.GT.1)
C    &   WRITE(6,611) NIRR_PN,(NBASMOL(I),I=1,NIRR_PN),
C    &                (NORBMOL(I),I=1,NIRR_PN),
C    &                (NFRO(I),I=1,NIRR_PN),PNUCCO
611   FORMAT(' NIRR_PN :',I5/' NBASMOL(1:NIRR_PN),NORBMOL(1:NIRR_PN),
     &         NFRO(1:NIRR_PN) '/24I5/' MODIFIED CORE-ENERGY',D13.6)
*
      call rdmole2(WORK(KLSWINT),NORBMOL,NBASMOL,NFRO,NTOTINT,
     &             XINT_L)
CTF
*
      call memmar(KDUMMY,IDUMMY,'FLUSM ', IDUMMY,'GTMOSW')
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine integral_type(ISOB,NGAS,ITPSV,LUINGR)
*
#include "implicit.h"
#include "ipoist8.inc"
      character*11 CELEC2
      character*13 CELEC
      dimension ISOB(NGAS,4)
      ISUM = 0
      do J=1,4,1
        do I=1,NGAS,1
          ISUM = ISUM + ISOB(I,J)
        end do
      end do
      if (ISUM.eq.2) CELEC = 'One-electron '
      if (ISUM.eq.4) CELEC = 'Two-electron '
      NAC = 0
      NBC = 0
      NAA = 0
      NBA = 0
      do I=1,NGAS,1
        NAC = NAC + ISOB(I,1)
        NBC = NBC + ISOB(I,2)
        NAA = NAA + ISOB(I,3)
        NBA = NBA + ISOB(I,4)
      end do
      if (NAC.eq.NAA.and.NBC.eq.NBA) CELEC2 = 'spin-free  '
      if (NAC.ne.NAA.or.NBC.ne.NBA)  CELEC2 = 'spin-orbit '
      if (ISUM.ne.ITPSV) 
     &     write(LUINGR,'(A,A,A)') CELEC,CELEC2,'integrals.'
      ITPSV = ISUM
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION I2EAD2(ISM,IORB,JSM,JORB,KSM,KORB,LSM,LORB)
*
* Find adress of integral in LUCIA order
*
* Indeces IORB,JORB,KORB,LORB are relative to start of symmetry block of orbs.
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
*
      I2EAD2 = I2EAD2S(ISM,IORB,JSM,JORB,KSM,KORB,LSM,LORB,WORK(KPINT2))
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION I2EAD2S(ISM,IORB,JSM,JORB,KSM,KORB,LSM,LORB,IJKLOF)
*
* Obtain address of integral (IORB JORB ! KORB LORB) in MOLCAS order
* IORB JORB KORB LORB corresponds to SYMMETRY ordered indeces !!
* and are relative to start of oribitals with given symmetry
* Integrals assumed in core
*
*
* Constructed for the MOLECULE-LUCIA interface for dbg calculations
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cgas.inc"
#include "csm.inc"
*
      Dimension IJKLOF(NsmSH,NsmSH,NsmSH)
      Logical iSymj,kSyml,ISYMK,JSYML,ijSymkl,IKSYMJL
      Logical ijklPerm
*.
      NTEST = 000
*
      If( NTEST.ge. 100) THEN
        write(6,*) ' GMIJKL at your service '
        WRITE(6,*) ' IORB ISM ',IORB,ISM
        WRITE(6,*) ' JORB JSM ',JORB,JSM
        WRITE(6,*) ' KORB KSM ',KORB,KSM
        WRITE(6,*) ' LORB LSM ',LORB,LSM
      END IF
*
      IABS = 0
      JABS = 0
      KABS = 0
      LABS = 0
*
      If (jSm.gt.iSm .or. (iSm.eq.jSm .and. JABS.gt.IABS)) Then
        iSym=jSm
        jSym=iSm
        I = JORB
        J = IORB
      Else
        iSym=iSm
        jSym=jSm
        I = IORB
        J = JORB
      End If
      ijBlk=jSym+iSym*(iSym-1)/2
      If (lSm.gt.kSm  .or. (kSm.eq.lSm .and. LABS.gt.KABS)) Then
        kSym=lSm
        lSym=kSm
        K = LORB
        L = KORB
      Else
        kSym=kSm
        lSym=lSm
        K = KORB
        L = LORB
      End If
      klBlk=lSym+kSym*(kSym-1)/2
*
      ijklPerm=.false.
      If ( klBlk.gt.ijBlk ) Then
        iTemp=iSym
        iSym=kSym
        kSym=iTemp
        iTemp=jSym
        jSym=lSym
        lSym=iTemp
        iTemp=ijBlk
        ijBlk=klBlk
        klBlk=iTemp
        ijklPerm=.true.
*
        iTemp = i
        i = k
        k = itemp
        iTemp = j
        j = l
        l = iTemp
      End If
      If(Ntest .ge. 100 ) then
        write(6,*) ' i j k l ',i,j,k,l
        write(6,*) ' Isym,Jsym,Ksym,Lsym',Isym,Jsym,Ksym,Lsym
      End if
*
*  Define offset for given symmetry block
      IBLoff = IJKLof(Isym,Jsym,Ksym)
      If(ntest .ge. 100 )
     &WRITE(6,*) ' IBLoff Isym Jsym Ksym ', IBLoff,ISym,Jsym,Ksym
      iSymj=iSym.eq.jSym
      kSyml=kSym.eq.lSym
      iSymk=iSym.eq.kSym
      jSyml=jSym.eq.lSym
      ikSymjl=iSymk.and.jSyml
      ijSymkl=iSymj.and.kSyml
*
      itOrb=NGSSHT(iSym)
      jtOrb=NGSSHT(jSym)
      ktOrb=NGSSHT(kSym)
      ltOrb=NGSSHT(lSym)
C?    print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
      If ( iSymj ) Then
        ijPairs=itOrb*(itOrb+1)/2
        ij=j+i*(i-1)/2
      Else
        ijPairs=itOrb*jtOrb
        ij=j + (i-1)*jtOrb
      End if
*
      IF(KSYML ) THEN
        klPairs=ktOrb*(ktOrb+1)/2
        kl=l+k*(k-1)/2
      ELSE
        klPairs=ktOrb*ltOrb
        kl=l+(k-1)*ltOrb
      End If
C?    print *,' ijPairs,klPairs',ijPairs,klPairs
*
      If ( ikSymjl ) Then
        If ( ij.gt.kl ) Then
          klOff=kl+(kl-1)*(kl-2)/2-1
          ijkl=ij+(kl-1)*ijPairs-klOff
        Else
          ijOff=ij+(ij-1)*(ij-2)/2-1
          ijkl=kl+(ij-1)*klPairs-ijOff
        End If
      Else
        ijkl=ij+(kl-1)*ijPairs
      End If
      If( ntest .ge. 100 )
     & write(6,*) ' ijkl ', ijkl
*
      I2EAD2S = iblOff-1+ijkl
      If( NTEST .ge. 100 ) then
	write(6,*) 'i j k l ', i,j,k,l
	write(6,*) ' ibloff ijkl ',ibloff,ijkl
        write(6,*) ' I2EADS  = ', I2EADS
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE INTIM_REL(IPRTRA)
*
* Interface to external integrals
*
* If NOINT .ne. 0, only pointers are constructed
* Jeppe Olsen, Winter of 1991
*              Last revision : Jan 98, double group integrals added
* Note : In this version only pointers are constructed, the rest is 
*        realized in picasso.
*    Timo Fleig
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "glbbas.inc"
#include "clunit.inc"
#include "lucinp.inc"
#include "csm.inc"
#include "cintfo.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "input.inc"
      COMMON/CECORER/ECORE
*
      NTEST = 0000
*
*
*. : Pointers for symmetry blocks of integrals
*
      if (SIGDEN_ROUTE.eq.'TRADIT') then
        CALL INTPNT_PNG(WORK(KPINT1),WORK(KPINT2))
*
* Pointers for double group integrals
        IF(IDBG.EQ.1) CALL INTPNT_DBG
      end if
      IF (SIGDEN_ROUTE .ne.'SIGDEN') THEN
C
CSK     22-11-2007 not important for SIGDEN_ROUTE == SIGDEN
C
C       Pointer for orbital indices for symmetry blocked matrices
C
        CALL ORBINH1(WORK(KINH1),NTOOBS,NTOOB,NSMOB)
      END IF
C
      I_AM_TIMO_WORLD = 0
      if (INTIMP.eq.8) I_AM_TIMO_WORLD = 1
*
      IF (I_AM_TIMO_WORLD.EQ.0.AND.
     &       INTIMP.EQ.1.AND.NOINT.EQ.0) THEN
*
* ===============
* . MOLCAS format
* ===============
*
        WRITE(6,*) ' Integrals imported from MOLCAS files'
*.Initialize buffers, open
C       CALL MKLUNDIO
*. Load one-electron integrals
CTF     CALL GETH0(WORK(KINT1))
*. And two-electron integrals if desired
        IF(INCORE.EQ.1) THEN
          CALL INTIMM_REL(WORK(KINT2),NIRR_PN,SYDI,IPRTRA)
        END IF
*
* Check read integrals
        if (NTEST.ge.1000) then
          do IINT = 0,NINT2-1,1
            write(6,'(A4,I3,1X,F18.14)') 'No. ',IINT+1,WORK(KINT2+IINT)
          end do
        end if
*
* Formatted input of symmetry non-vanishing integrals
*
      ELSE IF (INTIMP.EQ.3.AND.NOINT.EQ.0 ) THEN
        WRITE(6,*) ' Integrals imported formatted (E22.15) '
*.1 : One-electron integrals
        REWIND LU2INT
        READ(LU2INT,'(E22.15)') (WORK(KINT1-1+INT1),INT1=1,NINT1)
*.2 : Two-electron integrals
        IF(INCORE.EQ.1.OR.EXTSPC.EQ.0) THEN
          READ(LU2INT,'(E22.15)') (WORK(KINT2-1+INT2),INT2=1,NINT2)
        END IF
*.3 : Core energy
        WRITE(6,*) ' Core energy will be read '
        READ(LU2INT,'(E22.15)') ECORE
*
*. Proceed in Normal MOLCAS way so
        INTIMP = 1
*
* ===============
* . SIRIUS format
* ===============
*
      ELSE IF (INTIMP .EQ. 5.AND.NOINT.EQ.0  ) THEN
        WRITE(6,*) ' Integrals imported from SIRIUS files'
*. Load one-electron integrals
CTF     CALL GETH0S(WORK(KINT1),NOCOB)
*. And two-electron integrals if desired
        IF(INCORE.EQ.1) THEN
          LSCR = 100000
          call memmar(KLSCR,LSCR,'ADDL  ',2,'KLSCR ')
CTF       CALL READMO(WORK(KLSCR),LSCR,WORK(KINT2))
        END IF
      else if (I_AM_TIMO_WORLD.eq.1.and.INTIMP.eq.8
     &         .and.NOINT.eq.0) then
*
* ===================
*  MOLECULE format
* ===================
*
         if (ARCHI.eq.'RS6') then
           write(6,*) 'Reading from file 13 in MOLECULE format'
           write(6,*) 'Reading from file 13 in MOLECULE format'
           write(6,*) 'Reading from file 13 in MOLECULE format'
           write(6,*) 'Reading from file 13 in MOLECULE format'
         else if (ARCHI.eq.'SGO') then
           write(6,*) 'Reading from ASCII files'
           write(6,*) 'Reading from ASCII files'
           write(6,*) 'Reading from ASCII files'
           write(6,*) 'Reading from ASCII files'
         end if
CTF
*  Interface to MOLECULE-SWEDEN integrals (Jan Almloef, 1986)
*  Patch from program tradpt written by Christel Marian in 1987.
CTF
         call gtmoswe(WORK(KINT2),NTOTINT)
         write(6,*) 'Total number of 2-el. ints
     &               that have been read:         ',NTOTINT
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
      SUBROUTINE INTIMM_REL(XINT,MAXSYM,IDIM,IPRT)
*
* Import all two electron integrals from MOTRA 2e-file
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"

#include "SysDef.inc"

*.ORBINP
#include "mxpdim.inc"
#include "orbinp.inc"
*.LUNDIO
      Parameter ( mxBatch = 106  )
      Parameter ( mxSyBlk = 666  )
      Common / LundIOR / LuTr2,lTr2Rec,iDAdr(mxBatch),nBatch(mxSyBlk)
*
      dimension IDIM(MAXSYM)
*.Output
      DIMENSION XINT(*)
*. For testing
      Ntestl = 0
      Ntest = max(NTESTL,IPRT)
*
      if (NTEST.ge.1) then
        write(6,*) '--------------------------------'
        write(6,*) ' Information from INTIMM        '
        write(6,*) '--------------------------------'
      end if
*
      Ioff = 1
      Do 101 Ism = 1, Maxsym
        Do 102 Jsm = 1,Ism
          Do 103 Ksm = 1, ISm
            If(Ism .eq. Ksm ) Then
             LsmMX = Jsm
            Else
             LsmMX = Ksm
            End if
            DO 104 Lsm = 1, LsmMX
            If ( ieor(iSm-1,jSm-1).ne.ieor(kSm-1,lSm-1) )  goto 104
*
              IJbl = Ism*(Ism-1)/2 + Jsm
              KLbl = Ksm*(Ksm-1)/2 + Lsm
              IJKLbl = IJbl*(IJbl-1)/2 + KLbl
*
              NIORB = IDIM(ISM)
              NJORB = IDIM(JSM)
              NKORB = IDIM(KSM)
              NLORB = IDIM(LSM)
*
              If(Ism.Eq.Jsm) Then
                nIJ = NIorb*(NIorb+1)/2
              Else
                nIJ = NIORB*NJORB
              End If
              If(Ksm.Eq.Lsm) Then
               nKL = nKorb*(nKorb+1)/2
              Else
                nKL = nKorb*nLorb
              End If
*
              If(Ism .Eq. Ksm .And. Jsm .Eq. Lsm ) Then
                nIJKL= nIJ*(nIJ+1)/2
              Else
                nIJKL= nIJ*nKL
              End if
*
              If ( Ntest .Ne. 0 ) then
                Write(6,*) ' Ism Jsm Ksm Lsm ', Ism,Jsm,Ksm,Lsm
                Write(6,*) ' Ioff, nIJKL ', Ioff,nIJKL
              End if
              Ibatch = nbatch(IJKLbl)
              iDisk=iDAdr(iBatch)
              If(Ntest.Ne.0) Write(6,*)
     &        ' IJKLbl iBatch iDisk ', IJKLbl,iBatch,iDisk
*. Loop over records
              nRec = nIJKL/lTr2Rec
              If(Nrec*lTr2Rec.Lt. nIJKL) nRec = nRec + 1
              Ioffo = Ioff
              Do 50 IRec = 1, Nrec
                If ( IRec .Ne. Nrec ) Then
                   Nintrc = lTr2Rec
                Else
                   Nintrc = nIJKL -(Nrec-1)*Ltr2Rec
                End if
                Ioff = Ioff + Nintrc
   50        Continue
*
             If( Ntest .Ne. 0 ) then
               Write(6,*) ' Integral block '
               If(Ism.Eq.Ksm .And. Jsm. Eq. Lsm ) Then
                 Call Prsym(Xint(Ioffo),nIJ)
               Else
                 Call Wrtmat(Xint(Ioffo),nIJ,nKL,nIJ,nKL)
               End if
             End if
*
  104       Continue
  103     Continue
  102   Continue
  101 Continue
*
      Return
      End
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Set KR array according to number of functions in fermion
* symmetries.
*
      subroutine kkr_set(KR,IDIM,ISYDI,NIRR,INVSYM)
*
#include "implicit.h"
#include "ipoist8.inc"
*
      dimension KR(-IDIM:IDIM),ISYDI(NIRR)
*
      NTESTL = 0
*
      ICT = 1
      IV = 1
      do IRR = 1,1+INVSYM,1
        IRRR = IRR
        if (INVSYM.eq.1.and.NIRR.gt.2.and.IRR.eq.2) IRRR = IRR + 1
        do I = 1,ISYDI(IRRR),1
          KR(IV) = ICT
          KR(-IV) = ICT + ISYDI(IRRR)
          ICT = ICT + 1
          IV = IV + 1
        end do
        ICT = ICT + ISYDI(IRR)
      end do
*
      if (NTESTL.ge.1) then
        write(6,*) ' KKR has set KR array :'
        do ID = 1,IDIM,1
          write(6,'(3X,2I4)') KR(-ID),KR(ID)
        end do
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  This little routine transfers input information
*  from LUCIA to PICASSO, so that PICASSO no longer
*  needs an input of his own ...
*
      subroutine luc2pic(ISOOP,IPRTRA)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cgas.inc"
#include "symm.inc"
#include "input.inc"
#include "crun.inc"
*
      NGASPI = NGAS
      if (NGAS.gt.MXPNGAS) then
        write(6,*) 'Too many GAS spaces.'
        write(6,*) 'Number of GAS spaces is ',NGAS
        write(6,*) 'Increase MXPNGAS in mxpdim.inc or reconsider.'
        Call Abend2( 'Quitting.' )
      end if
*
      NOBPIRR = 0
      IOFPGSM(1) = 0
      do IRR=1,NIRR_PN,1
         do IGAS=1,NGASPI,1
            NOBPIRR = NOBPIRR + NGSSH(IRR,IGAS)
         end do
         if (IRR.lt.NIRR_PN) IOFPGSM(IRR+1) = IOFPGSM(IRR) + NOBPIRR
         ITOTDM(IRR) = NOBPIRR
         SYDI(IRR) = ITOTDM(IRR)
         NOBPIRR = 0
      end do
*
* Determine largest symmetry block
* Determine triangle lengths
* ... and the sum of these
      J = 1
      IMAXDIM = J
      ITRIASUM = 0
      ITRIASUM_F = 0
      do I=1,NIRR_PN,1
         if (SYDI(I).gt.SYDI(J)) then
            IMAXDIM = I
            J = I
         end if
*
         ITRIA(I) = (ITOTDM(I)**2 + ITOTDM(I)) / 2
         ITRIA_F(I) = (NORB_F(I)**2 + NORB_F(I)) / 2
         ITRIASUM = ITRIASUM + ITRIA(I)
         ITRIASUM_F = ITRIASUM_F + ITRIA_F(I)
*
         do J=1,NGASPI,1
            NGSOBPI(I,J) = NGSSH(I,J)
         end do
      end do
*
* ISO1 = ISOOP in readin.
* '0' means inactive, '1' means active SO operator of this type.
* ISO2 not implemented yet.
      if (ISOOP.eq.0) IOPER = 1   ! Only sf one-electron operator
      if (ISOOP.eq.1) IOPER = 3   ! sf+SO one-electron operator
C     Only SO operator doesn't make any sense!
*
      IPRINT = IPRTRA
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Calculate total number of integrals on record
* of a symmetry reduced list (r/i, particle symmetry)
*
      subroutine numints_red_rel(IGLSUM,NSMBLK,IPRINT)
*
*------------------------------------
*  Timo Fleig, October 2000
*------------------------------------
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "symm.inc"
#include "parluci.h"
*
      dimension IDEKS(MXPORB),IQISYM(MXPIRR)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      if (NTEST.ge.1) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================='
        write(LUWRT,*) ' numints_red_rel speaking:'
        write(LUWRT,*) '=========================='
        if (NTEST.ge.5) then
          write(LUWRT,*) 'Point group multiplication table:'
          call iwrtma(IPGMULT,NIRR_PN,NIRR_PN,MXNSPIR,MXNSPIR)
          write(LUWRT,*) 'Symmetry dimension array:'
          call iwrtma(SYDI,NIRR_PN,1,NIRR_PN,1)
        end if
      end if
*
* Set up triangular array
      II = 0
      do I=1,MXPORB,1
        IDEKS(I) = II
        II = II + I
      end do
*
* Symmetry packing labels
      NOCC = 0
      do IS=1,NIRR_PN,1
        if (SYDI(IS).eq.0) goto 10
        NOCC = NOCC + 1
        IQISYM(NOCC) = IS
   10   continue
      end do
      if (NTEST.ge.5) then
        write(LUWRT,*) 'Symmetry packing labels:'
        call iwrtma(IQISYM,NOCC,1,NOCC,1)
      end if
*
* Number of symmetry blocks and total number of integrals
      INIT = 1
      do ISP=1,NIRR_PN,1
        if (NTEST.ge.10) write(LUWRT,*) 'ISP =',ISP
        NSP=IQISYM(ISP)
        if (NSP.eq.0) then
          write(LUWRT,*) 'Number of symm for P is  ',NSP
          write(LUWRT,*) 'Skipping two-el. int. read.'
          goto 333
        end if
        if (NTEST.ge.10) write(LUWRT,*) 'NSP =',NSP
        NOP=SYDI(NSP)
        do ISQ=1,ISP,1
          if (NTEST.ge.10) write(LUWRT,*) 'ISQ =',ISQ
          NSQ = IQISYM(ISQ)
          NOQ = SYDI(NSQ)
          ISPQ = IPGMULT(NSP,NSQ)
          if (NTEST.ge.10) write(LUWRT,*) 'NSQ =',NSQ
          do ISR=1,ISP,1
            if (NTEST.ge.10) write(LUWRT,*) 'ISR =',ISR
            NSR = IQISYM(ISR)
            NOR = SYDI(NSR)
            ISPQR = IPGMULT(ISPQ,NSR)
            ISSM=ISR
            if (ISR.EQ.ISP) ISSM=ISQ
            if (NTEST.ge.10) write(LUWRT,*) 'NSR =',NSR
            do 40 ISS=1,ISSM,1
              if (NTEST.ge.10) write(LUWRT,*) 'ISS =',ISS
              NSS=IQISYM(ISS)
              NOS=SYDI(NSS)
              if (NTEST.ge.10) then
                write(LUWRT,'(A,4I6)') 'NOP,NOQ,NOR,NOS',NOP,NOQ,NOR,NOS
                write(LUWRT,'(A,2I3)') 'ISPQ,ISPQR',ISPQ,ISPQR
              end if
              if (NTEST.ge.10) write(LUWRT,*) 'NSS =',NSS
              if(NSS.NE.ISPQR) goto 40
              if(NOS*NOR*NOQ*NOP.eq.0) goto 40
*
*  NUMBER OF INTEGRALS IN THIS SYMMETRYBLOCK
*
              if (NSP.EQ.NSQ) then
*
*  TYPES <PP|PP> AND <PP|QQ>
*
                NPQ = IDEKS(NOP+1)
                NRS = IDEKS(NOR+1)
              else
*
*  TYPES <PQ|PQ> AND <PQ|RS>
*
                NPQ = NOP*NOQ
                NRS = NOR*NOS
              end if
              if (NSP.EQ.NSR) then
*
*  TYPES <PP|PP> AND <PQ|PQ>
*
                NIGL = NPQ*(NPQ+1)/2
              else
*
*  TYPES <PP|QQ> AND <PQ|RS>
*
                NIGL = NPQ*NRS
              end if
              INIT = INIT + 1
              IGLSUM = IGLSUM + NIGL
              if (NTEST.ge.1) then
                 write(LUWRT,601) NSP,NSQ,NSR,NSS,NIGL
              end if
40          continue
          end do
        end do
* Number of non-vanishing symmetry blocks
        NSMBLK = INIT - 1
333     continue
      end do
601   FORMAT(' NUMBER OF INTEGRALS IN SYMMETRYBLOCK '
     &             ,4I2,' IS : ',I10)
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
      subroutine rd1int(ONEINT,FNONEAO,LUONEAO,NATOMS,IFOCKEX)
************************************************************************
*                                                                      *
*     Objective: Read the header of the one-electron integral file     *
*                Extract also symmetry and basis set information.      *
*                In addition read the overlap, the nuclear attraction  *
*                and kinetic integrals.                                *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*                                                                      *
*     modified for LUCIAREL adapted to MOLCAS (4.1, 4.3 running)       *
*                                                                      *
*     Timo Fleig, 16.05.99    (another sunny day in may; a sun-day)    *
*                              may-day ?
*                 25-27.05.99 (cont'd, more sunny days in may ...)     *
*                 30.11.99    (read nucl. rep. en. when INPORB ..)     *
*                             (A fucking rainy and windy november day) *
*                                                                      *
************************************************************************
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "symm.inc"
#include "onel.inc"
#include "input.inc"
*
      common/CECORER/ECORE
*
      character*6 FNONEAO
      character*8 OneLbl,ComLbl
      character*72 HEADER
*
      dimension NBAS_SCF(MXNSPIR)
      dimension COOR_SCF(3,MXNATOMS),ONEINT(ITRIASUM)
      dimension HONE_SCF(MXPORB*MXPORB)
*
      NTESTL = 0
      NTEST = max(IPRINT,NTESTL)
*
*  Initialize
      do IINT=1,ITRIASUM,1
        ONEINT(IINT) = 0.0
      end do
*
      Call qEnter('rd1int')
*----------------------------------------------------------------------*
*     Start procedure                                                  *
*     Open the one electron integral file                              *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      Call OpnOne(iRC,iOpt,FNONEAO,LUONEAO)
      If ( iRc.ne.0 ) Goto 990
*----------------------------------------------------------------------*
*     Read one-electron integral file header etc.                      *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='Title   '
      Call RdOne(iRc,iOpt,OneLbl,iComp,HEADER,iSyLbl)
C     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read no.of symm. species                                         *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='nSym    '
      Call RdOne(iRc,iOpt,OneLbl,iComp,NSYM_SCF,iSyLbl)
      if (NSYM_SCF.ne.NIRR_PN) then
        write(6,*) 'Number of symms on comfile does not'
        write(6,*) 'equal those from input.'
        Call Abend2( 'severe error, quitting.' )
      else
        write(6,*) 'Number of pg symmetries: ',NSYM_SCF
      end if
C     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read symm. oper per symm. species                                *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=0
c     iComp=0
c     iSyLbl=1
c     OneLbl='SymOp   '
c     Call RdOne(iRc,iOpt,OneLbl,iComp,IOPER_SCF,iSyLbl)
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read no. of basis functions per symm. species                    *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='nBas    '
      Call RdOne(iRc,iOpt,OneLbl,iComp,NBAS_SCF,iSyLbl)
      do ISYM=1,NIRR_PN,1
        write(6,*) 'Dimension of symm ',ISYM,' is ',NBAS_SCF(ISYM)
        if (SYDI(ISYM).ne.NBAS_SCF(ISYM)) then
          write(6,*) 'Dimensions of spatial symmetries corrupted'
          write(6,'(A5,I1,A4,I3)')
     &             'SYDI(',ISYM,') = ',SYDI(ISYM)
          write(6,'(A9,I1,A4,I3)')
     &             'NBAS_SCF(',ISYM,') = ',NBAS_SCF(ISYM)
          Call Abend2( 'fatal error in rd1int' )
        end if
      end do
C     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read no. of basis functions per symm. species                    *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=0
c     iComp=0
c     iSyLbl=1
c     OneLbl='AtomType'
c     Call RdOne(iRc,iOpt,OneLbl,iComp,BsLbl,iSyLbl)
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read no. of unique atoms in the system                           *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='nAtoms  '
      Call RdOne(iRc,iOpt,OneLbl,iComp,NATOMS,iSyLbl)
      write(6,*) 'Number of unique atoms:  ',NATOMS
      if (NATOMS.gt.MXNATOMS) then
        write(6,*) 'Number of atoms is too large.'
        write(6,*) 'Adjust MXNATOMS in mxpdim.inc or reconsider.'
        Call Abend2( 'Quitting.' )
      end if
C     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read atom labels                                                 *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=0
c     iComp=0
c     iSyLbl=1
c     OneLbl='AtomLbl '
c     Call RdOne(iRc,iOpt,OneLbl,iComp,AtLbl,iSyLbl)
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read coordinates of atoms                                        *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='Coord   '
      Call RdOne(iRc,iOpt,OneLbl,iComp,COOR_SCF,iSyLbl)
      do IATM=1,NATOMS,1
        write(6,'(A9,I3,A1,1X,3F10.6)') 'Atom no. ',IATM,':',
     &                                 (COOR_SCF(IC,IATM),IC=1,3,1)
      end do
C     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read nuclear repulsion energy                                    *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      iComp=0
      iSyLbl=1
      OneLbl='PotNuc  '
      Call RdOne(iRc,iOpt,OneLbl,iComp,POTNUC_SCF,iSyLbl)
      write(6,'(A,1X,F18.15)') 'Nuclear repulsion energy: ',POTNUC_SCF
C     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read overlap integrals                                           *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=6
c     iComp=1
c     iSyLbl=1
c     OneLbl='Mltpl  0'
c     Call RdOne(iRc,iOpt,OneLbl,iComp,HONE_SCF,iSyLbl)
c     write(6,*)
c     write(6,*) 'Overlap integrals:'
c     do IINT=1,ITRIASUM,1
c       write(6,'(A,X,I4,A,F18.14)')
c    &        'Symm orb integral ',IINT,':',HONE_SCF(IINT)
c     end do
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read core Hamiltonian                                            *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      OneLbl='OneHam  '
      Call RdOne(iRc,iOpt,OneLbl,iComp,HONE_SCF,iSyLbl)
      write(6,*)
      write(6,*) 'Reading core Hamiltonian'
      write(6,*) 'from MOLCAS file ONEINT ...'
      write(6,*)
*  Copy to SO integral storage array
      do IINT=1,ITRIASUM,1
        ONEINT(IINT) = ONEINT(IINT) + HONE_SCF(IINT)
      end do
*  Check
      if (NTEST.ge.1) then
        do IINT=1,ITRIASUM,1
          write(6,'(A,1X,I4,A,F20.12)')
     &          'Symm orb integral ',IINT,':',HONE_SCF(IINT)
        end do
      end if
C     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read nuclear attraction integrals                                *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=6
c     iComp=1
c     iSyLbl=1
c     OneLbl='Attract '
c     Call RdOne(iRc,iOpt,OneLbl,iComp,HONE_SCF,iSyLbl)
c     write(6,*)
c     write(6,*) 'Nuclear attraction integrals:'
c     do IINT=1,ITRIASUM,1
c       write(6,'(A,X,I4,A,F20.12)')
c    &        'Symm orb integral ',IINT,':',HONE_SCF(IINT)
c     end do
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read kinetic energy integrals                                    *
*----------------------------------------------------------------------*
c     iRc=-1
c     iOpt=6
c     iComp=1
c     iSyLbl=1
c     OneLbl='Kinetic '
c     Call RdOne(iRc,iOpt,OneLbl,iComp,HONE_SCF,iSyLbl)
c     write(6,*)
c     write(6,*) 'kinetic energy integrals:'
c     do IINT=1,ITRIASUM,1
c       write(6,'(A,X,I4,A,F20.12)')
c    &        'Symm orb integral ',IINT,':',HONE_SCF(IINT)
c     end do
c     If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     If this is a perturbative reaction field calculation then        *
*     modify the one-electron Hamiltonian by the reaction field and    *
*     the nuclear attraction by the cavity self-energy                 *
*----------------------------------------------------------------------*
c     If ( iRFpert.ne.0 ) then
c        nTemp=0
c        Do iSym=1,nSym
c           nTemp=nTemp+nBas(iSym)*(nBas(iSym)+1)/2
c        End Do
c        Call GetMem('RFFLD','Allo','Real',lTemp,nTemp)
c        iRc=-1
c        iOpt=0
c        Call OpnCom(iRc,iOpt,FnCom,LuCom)
c        If ( iRc.ne.0 ) Goto 993
c        iRc=-1
c        iOpt=0
c        ComLbl='ERFself'
c        Call RdCom(iRc,iOpt,ComLbl,ERFself,2)
c        PotNuc=PotNuc+ERFself
c        If ( iRc.ne.0 ) Goto 994
c        iRc=-1
c        iOpt=0
c        ComLbl='RFFLD'
c        Call RdCom(iRc,iOpt,ComLbl,Work(lTemp),2*nTemp)
c        Call Daxpy(nTemp,1.0D0,Work(lTemp),1,Hone,1)
c        If ( iRc.ne.0 ) Goto 994
c        iRc=-1
c        iOpt=0
c        Call ClsCom(iRc,iOpt)
c        If ( iRc.ne.0 ) Goto 995
c        Call GetMem('RFFLD','Free','Real',lTemp,nTemp)
c     End If
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=0
      call clsone(iRc,iOpt)
      If ( iRc.ne.0 ) Goto 992
*
c     nTot=0
c     nTot1=0
c     nTot2=0
c     n2max=0
c     Do iSym=1,nSym
c       iBas=nBas(iSym)
c       nTot=nTot+iBas
c       nTot1=nTot1+iBas*(iBas+1)/2
c       nTot2=nTot2+iBas*iBas
c       n2max=Max(n2max,iBas*iBas)
c     End Do
      Call qExit('Rd1Int')
*
      if (IFOCKEX.eq.0.and.ORBFILE) then
        ECORE = POTNUC_SCF
        write(6,*) 'Passing POTNUC to ECORE'
      end if
*
      write(6,*) 'file opened successfully, returning ...'
      Return
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
990   Call ErrMsg('Rd1Int',001)
      Call Quit(16)
C991   Call StkA(0,OneLbl)
      Call ErrMsg('Rd1Int',002)
      Call Quit(16)
992   Call ErrMsg('Rd1Int',003)
      Call Quit(16)
993   Call ErrMsg('Rd1Int',004)
      Call Quit(16)
C994   Call StkA(0,ComLbl)
      Call ErrMsg('Rd1Int',005)
      Call Quit(16)
995   Call ErrMsg('Rd1Int',006)
      Call Quit(16)
      End
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Import integrals from MOLCAS/SEWARD environment
*  For run with externally generated Fock matrix.
*  (Driver routine)
*
      subroutine rd_traone(ONEINT,LU2,FNONEMO,IPRINT)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "symm.inc"
*
      character*6 FNONEMO
      dimension ONEINT(ITRIASUM)
*
C     call geth0(ONEINT(1),LU2,NIRR_PN,SYDI(1),
C    &           NFRO_L(1),FNONEMO,IPRINT)
C Not required in DIRAC environment.
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Read NR orbitals from RASSCF file JOBIPH
*  or from formatted ASCII file
*
      subroutine rdcoef_molcas(COEF,CMOAO,HELPMAT,
     &                  ICHECK,ITOTDIM,ISYMDIM,ICRMXDIM,
     &                  ICR2DIM,ICR3DIM,
     &                  IDIM,NIRR_PN,NIRR_DG,FNIPH,LUIPH,
     &                  INTITER,IMAXDIM,IOUTLST,IDOUGRP,
     &                  IFOCKEX,LUMORB,LUMOUT,ITRA_TYPE,
     &                  LUMOUTR,LUMOUTI)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "coeff.inc"
#include "input.inc"
#include "SysDef.inc"
#include "files_r.inc"
*
*   Coefficient common block especially for Jeppe's routines
*
*     integer NBSPFSM,NUSPFSM,IBSPEXP
*     real*8 CSPEXP
*     dimension CSPEXP(MXNSPIPDMGRP**2)
*     dimension NBSPFSM(MXNSPIR),NUSPFSM(MXNSPIR)
*     dimension IBSPEXP(2,2),ISYMDIM(NIRR_PN)
*     common /COEFF/CSPEXP,NBSPFSM,NUSPFSM,IBSPEXP
*
      character*2 HEADER(72)
      character*6 FNIPH
      character*72 TITLE(18)
      logical okay
*
      common/CECORER/ECORE
*
      dimension IDIM(MXNSPIR)
      dimension COEF(ITOTDIM,ITOTDIM,2)
      dimension CMOAO(NIRR_PN*ICRMXDIM**2)
      dimension HELPMAT(ICR2DIM**2+ICR3DIM**2)
      dimension IOLIST(102),IPAD(16)
      dimension ISYMDIM(NIRR_PN)
      dimension NFRO_M(MXNSPIR),NISH_M(MXNSPIR),NASH_M(MXNSPIR),
     &          NDEL_M(MXNSPIR),NBAS_M(MXNSPIR),
     &          NAME_M(4*2*MXPORB/ITOB),IROOT_M(MXROOT),
     &          NRAS1_M(MXNSPIR),NRAS2_M(MXNSPIR),NRAS3_M(MXNSPIR)
*
      Call qEnter('rdcoef_molcas')
CTF
*  In VERY pathological cases, the field HELPMAT might be larger
*  than SYMDIM(IDIMGRP)**2*2, as it is dimensioned in oneltra.
*  Therefore, we check whether ICHECK has been destroyed or not:
      if (ICHECK.ne.987654)
     &    Call Abend2( 'Problem in rdcoef; HELPMAT too large' )
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
*  Either the coeffs are already on CSPAOSF, then skip the
*  following "instructions" or otherwise execute them; read
*  coeffs from file 19, construct matrix CSPAOSF.
*
*  determine total length required for CMOAO
      ILEN = 0
      do I=1,NIRR_PN,1
         ILEN = ILEN + IDIM(I)**2
      end do
C?    write(6,*) 'ILEN is    ',ILEN
*
 100  format(4F18.14)
      ZERO = 0.0D0
*
      if (INTITER.eq.1.and.IOUTLST.eq.1) then
*
        if (ORBFILE) then
*
          if (ITRA_TYPE.le.3) then
*----------------------------------------------------------------------*
*  Read coefficients from formatted ASCII file INPORB
*----------------------------------------------------------------------*
*
            IOFFS = 1
            open(unit=LUMORB,file='INPORB',
     &           status='OLD',form='FORMATTED')
*
              read(LUMORB,'(A)')
              do IRREP=1,NIRR_PN,1
                do I=1,IDIM(IRREP),1
                  read(LUMORB,'(A)')
                  read(LUMORB,100)
     &                (CMOAO(IN),IN=IOFFS,IOFFS+IDIM(IRREP)-1,1)
                  IOFFS = IOFFS + IDIM(IRREP)
                end do
              end do
            close (unit=LUMORB)
          else if (ITRA_TYPE.ge.4) then
*----------------------------------------------------------------------*
*  Read coefficients from formatted ASCII files LUMOUTR/I (opt. coeffs.)
*----------------------------------------------------------------------*
            open(unit=LUMOUTR,file=NLUMOUTR_MPI(1:LUMOUTRI_MPIL),
     &           status='OLD',form='FORMATTED')
              read(LUMOUTR,'(A)')
              do IRREP=1,IOUTLST-1,1
                do ICOL=1,ISYMDIM(IRREP),1
                  read(LUMOUTR,'(A)')
                  do IROW=1,ISYMDIM(IRREP),1
                    read(LUMOUTR,100)
                  end do
                end do
              end do
              do ICOL=1,ITOTDIM,1
                read(LUMOUTR,'(A)')
                read(LUMOUTR,100)
     &              (COEF(IROW,ICOL,1),IROW=1,ITOTDIM,1)
              end do
            close (unit=LUMOUTR)
*
            open(unit=LUMOUTI,file=NLUMOUTI_MPI(1:LUMOUTRI_MPIL),
     &           status='OLD',form='FORMATTED')
              read(LUMOUTI,'(A)')
              do IRREP=1,IOUTLST-1,1
                do ICOL=1,ISYMDIM(IRREP),1
                  read(LUMOUTI,'(A)')
                  do IROW=1,ISYMDIM(IRREP),1
                    read(LUMOUTI,100)
                  end do
                end do
              end do
              do ICOL=1,ITOTDIM,1
                read(LUMOUTI,'(A)')
                read(LUMOUTI,100)
     &              (COEF(IROW,ICOL,2),IROW=1,ITOTDIM,1)
              end do
            close (unit=LUMOUTI)
          end if
        else
*
*----------------------------------------------------------------------*
*     Read MO coefficients from JOBIPH generated by RASSCF             *
*----------------------------------------------------------------------*
*
          if (NTEST.ge.3) then
            write(6,*) 'Checking'
            write(6,*) 'MXNSPIR =',MXNSPIR
            write(6,*) 'MXPORB =',MXPORB
            write(6,*) 'MXTIT =',MXTIT
            write(6,*) 'MXROOT =',MXROOT
            write(6,*) 'RTOI = ',RTOI
            write(6,*) 'ITOB = ',ITOB
          end if
*
          Inquire (File=FNIPH,exist=okay)
          If ( okay ) Then
            IAD16 = 0
*
C           CALL GSLIST(IOLIST,23,NACTEL_M,1,ISPIN,1,NSYM,1,LSYM,1,
C    &             NFRO_M,MXNSPIR,NISH_M,MXNSPIR,NASH_M,MXNSPIR,
C    &             NDEL_M,MXNSPIR,NBAS_M,MXNSPIR,
C    &             NAME_M,4*2*MXPORB/ITOB,NCONF,1,HEADER,144/ItoB,
C    &             TITLE,4*18*MXTIT/ITOB,POTNUC,RtoI,
C    &             LROOTS,1,NROOTS,1,IROOT_M,MXROOT,NRAS1_M,MXNSPIR,
C    &             NRAS2_M,MXNSPIR,NRAS3_M,MXNSPIR,
C    &             NHOLE1,1,NELE3,1,IPT2,1,
C    &                                     x23,n23,x24,n24,
C    &             x25,n25,x26,n26,x27,n27,x28,n28,x29,n29,
C    &             x30,n30,x31,n31,x32,n32,x33,n33,x34,n34,
C    &             x35,n35,x36,n36,x37,n37,x38,n38,x39,n39,
C    &             x40,n40,x41,n41,x42,n42,x43,n43,x44,n44,
C    &             x45,n45,x46,n46,x47,n47,x48,n48,x49,n49)
*
            IAD16 = IPAD(1)
*
            IAD16 = IPAD(2)
*
            if (NTEST.ge.1) then
              write(6,*) 'Checking read coeffs from JOBIPH:'
              ICOFF = 1
              do ISYM=1,NIRR_PN,1
                write(6,'(A16,I3)') 'Symmetry number ',ISYM
                if (ISYM.ge.2) ICOFF = ICOFF + IDIM(ISYM-1)**2
                call wrtmat(CMOAO(ICOFF),IDIM(ISYM),IDIM(ISYM),
     &                            IDIM(ISYM),IDIM(ISYM))
              end do
            end if
          end if
        end if
*
        if (ITRA_TYPE.le.3) then
          call cspaosf_from_moao(CSPEXP,NBSPFSM(1),NUSPFSM(1),
     &                           IBSPEXP,CMOAO,IDIM(1),IDIM(1))
        end if
*
*  The coefficients have been written to the common block element
*  CSPEXP. The fetcher will get the appropriate blocks for this
*  transformation run and arrange them on COEF.
*
*  If ITRA_TYPE.ge.4 optimized spinors have been read in directly
*  from files and no further 'treatment' is required.
      end if
*
*  Erase output matrix COEF completely
      if (ITRA_TYPE.le.3) then
        do IR=1,2,1
          call setvec(COEF(1,1,IR),ZERO,ITOTDIM**2)
        end do
c       do I=1,ITOTDIM,1
c          do J=1,ITOTDIM,1
c             COEF(J,I,1) = 0.0
c             COEF(J,I,2) = 0.0
c          end do
c       end do
*
        if (IDOUGRP.eq.4) then
           ISPSM = IOUTLST + 4
        else if (IDOUGRP.ge.5.and.IDOUGRP.le.7) then
           ISPSM = IOUTLST + 2
        else if (IDOUGRP.eq.8) then
           ISPSM = 2
        end if
*
        do IREIM=1,2,1
           do IBU=1,2,1
              do IAB=1,2,1
                 call fetch_cspexp_blk(CSPEXP,IBSPEXP,IDIM(1),IDIM(1),
     &                                 IREIM,IBU,ISPSM,IAB,
     &                                 CMOAO,LROW,LCOL)
                 ILENBLK = LROW*LCOL
*
*  Determine whether alpha preceeds beta or vice versa
                 call evenodd(IEOLST,IOUTLST)
                 if (IEOLST.eq.1) IPERMAB = 1
                 if (IEOLST.eq.2) IPERMAB = 2
*
*  Determine active dimensions and offsets for vectomat
                 IDIMACT = 1
                 if (IDOUGRP.eq.7.and.IOUTLST.eq.2) IDIMACT = 2
                 if (IOUTLST.ge.3) IDIMACT = 3
                 if (IPERMAB.eq.1.and.IAB.eq.1) then
                    IROWOFF = 1
                    if (IBU.eq.1) ICOLOFF = IDIM(IDIMACT) + 1
                    if (IBU.eq.2) ICOLOFF = 1
                 end if
                 if (IPERMAB.eq.1.and.IAB.eq.2) then
                    IROWOFF = IDIM(IDIMACT) + 1
                    if (IBU.eq.1) ICOLOFF = IDIM(IDIMACT) + 1
                    if (IBU.eq.2) ICOLOFF = 1
                 end if
                 if (IPERMAB.eq.2.and.IAB.eq.1) then
                    IROWOFF = IDIM(IDIMACT) + 1
                    if (IBU.eq.2) ICOLOFF = IDIM(IDIMACT) + 1
                    if (IBU.eq.1) ICOLOFF = 1
                 end if
                 if (IPERMAB.eq.2.and.IAB.eq.2) then
                    IROWOFF = 1
                    if (IBU.eq.2) ICOLOFF = IDIM(IDIMACT) + 1
                    if (IBU.eq.1) ICOLOFF = 1
                 end if
*
                 if (NTEST.ge.10) then
                    write(6,*) 'Testing before vectomat ...'
                    write(6,*) 'IREIM =  ',IREIM
                    write(6,*) 'IBU =  ',IBU
                    write(6,*) 'IAB =  ',IAB
                    write(6,*) 'IROWOFF =  ',IROWOFF
                    write(6,*) 'ICOLOFF =  ',ICOLOFF
                 end if
*
                 call vectomat(CMOAO,COEF(1,1,IREIM),ILENBLK,
     &                         ITOTDIM,ITOTDIM,
     &                         LROW,LCOL,IROWOFF,ICOLOFF)
*
                 if (NTEST.gt.2) then
                    write(6,*) 'Checking coefficients ...'
                    call wrtmat(COEF(1,1,IREIM),ITOTDIM,ITOTDIM,
     &                          ITOTDIM,ITOTDIM)
                 end if
              end do
           end do
        end do
      end if
*     ^ End if type of transformation
*
      if (NTEST.gt.2) then
         write(6,*) 'Full coefficient matrices: real  imaginary'
         call wrtmat(COEF(1,1,1),ITOTDIM,ITOTDIM,ITOTDIM,ITOTDIM)
         call wrtmat(COEF(1,1,2),ITOTDIM,ITOTDIM,ITOTDIM,ITOTDIM)
      end if
*
      if (.not.ORBFILE) then
        if (IFOCKEX.eq.0) then
          ECORE = POTNUC
          write(6,*) 'Nuclear repulsion energy :  ',POTNUC
        end if
      end if
*
      call qexit('rdcoef_molcas')
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Get relevant block of coefficients from CSPAOSF(*) or
*  read coefficients from file 19
*
      subroutine rdcoef_molecu(ADJCO,COEF,CMOAO,HELPMAT,
     &                  ICHECK,ITOTDIM,ICRMXDIM,ICR2DIM,ICR3DIM,
     &                  IDIM,NIRR_PN,NIRR_DG,LU19,INTITER,
     &                  IMAXDIM,IOUTLST,IDOUGRP,ITRA_TYPE)
*------------------------------------------------------
* Modify according to ITRA_TYPE (like in rdcoef_molcas)
*------------------------------------------------------
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "coeff.inc"
#include "input.inc"
*
*   Coefficient common block especially for Jeppe's routines
*
*     integer NBSPFSM,NUSPFSM,IBSPEXP
*     real*8 CSPEXP
*     dimension CSPEXP(MXNSPIPDMGRP**2)
*     dimension NBSPFSM(MXNSPIR),NUSPFSM(MXNSPIR)
*     dimension IBSPEXP(2,2)
*     common /COEFF/CSPEXP,NBSPFSM,NUSPFSM,IBSPEXP
*
      dimension IDIM(MXNSPIR)
      dimension COEF(ITOTDIM,ITOTDIM,2)
      dimension ADJCO(ICRMXDIM,ICRMXDIM)
      dimension CMOAO(NIRR_PN*ICRMXDIM**2)
      dimension HELPMAT(ICR2DIM**2+ICR3DIM**2)
CTF
*  In VERY pathological cases, the field HELPMAT might be larger
*  than SYMDIM(IDIMGRP)**2*2, as it is dimensioned in oneltra.
*  Therefore, we check whether ICHECK has been destroyed or not:
      if (ICHECK.ne.987654)
     &    Call Abend2( 'Problem in rdcoef; HELPMAT too large' )
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
*  Either the coeffs are already on CSPAOSF, then skip the
*  following "instructions" or otherwise execute them; read
*  coeffs from file 19, construct matrix CSPAOSF.
*
*  determine total length required for CMOAO
      ILEN = 0
      do I=1,NIRR_PN,1
         ILEN = ILEN + IDIM(I)**2
      end do
C?    write(6,*) 'ILEN is    ',ILEN
*
*  Read irrepwise to array ADJCO, transfer to vector CMOAO, next irrep
      if (INTITER.eq.1.and.IOUTLST.eq.1) then
         if (IDOUGRP.eq.4) then
            IDIMREO4 = IDIM(4)
            IDIMREO3 = IDIM(3)
            IDIMREO2 = IDIM(2)
            IDIM(2) = IDIMREO3
            IDIM(3) = IDIMREO4
            IDIM(4) = IDIMREO2
         end if
         K = 1
         read(LU19,'(A)')
         do IRREP=1,NIRR_PN,1
            do I=1,IDIM(IRREP),1
               read(LU19,'(4F18.14)') (ADJCO(I,J),J=1,IDIM(IRREP),1)
            end do
            do I=1,IDIM(IRREP),1
               do J=1,IDIM(IRREP),1
                  CMOAO(K) = ADJCO(I,J)
                  K = K + 1
               end do
            end do
         end do
         if (IDOUGRP.eq.4) then
            IDIMREO4 = IDIM(4)
            IDIMREO3 = IDIM(3)
            IDIMREO2 = IDIM(2)
            IDIM(2) = IDIMREO4
            IDIM(3) = IDIMREO2
            IDIM(4) = IDIMREO3
         end if
*
         if (NTEST.ge.2) then
            do I=1,ILEN,1
               write(6,*) 'CMOAO(',I,') is  ',CMOAO(I)
            end do
         end if
*
*  Reorder, if IDOUGRP is C2H (4)
         if (IDOUGRP.eq.4) then
            IAUSTA = IDIM(1)**2 + 1
            IBUEND = IAUSTA + IDIM(3)**2 + IDIM(4)**2 - 1
            IBGSTA = IBUEND + 1
            IBGEND = IBGSTA + IDIM(2)**2 - 1
            IBGLEN = IDIM(2)**2
            if (NTEST.ge.3) then
               write(6,*) 'IAUSTA = ',IAUSTA
               write(6,*) 'IBUEND = ',IBUEND
               write(6,*) 'IBGSTA = ',IBGSTA
               write(6,*) 'IBGEND = ',IBGEND
               write(6,*) 'IBGLEN = ',IBGLEN
            end if
*
            J = 1
            do I=IBGSTA,IBGEND,1
               HELPMAT(J) = CMOAO(I)
               J = J + 1
            end do
            KLAST = ILEN
            do I=IBUEND,IAUSTA,-1
               CMOAO(KLAST) = CMOAO(I)
               KLAST = KLAST - 1
            end do
            IN = 1
            do I=IAUSTA,IAUSTA+IBGLEN-1,1
               CMOAO(I) = HELPMAT(IN)
               IN = IN + 1
            end do
*
            if (NTEST.ge.3) then
               write(6,*) 'Checking reordered MOAO coeffs in C2H:'
               do I=1,ILEN,1
                  write(6,*) 'CMOAO(',I,')= ',CMOAO(I)
               end do
            end if
         end if
*
         call cspaosf_from_moao(CSPEXP,NBSPFSM(1),NUSPFSM(1),
     &                          IBSPEXP,CMOAO,IDIM(1),IDIM(1))
*
*  The coefficients have been written to the common block element
*  CSPEXP. The fetcher will get the appropriate blocks for this
*  transformation run and arrange them on COEF.
      end if
*
*  Erase output matrix COEF completely
      do I=1,ITOTDIM,1
         do J=1,ITOTDIM,1
            COEF(J,I,1) = 0.0
            COEF(J,I,2) = 0.0
         end do
      end do
*
      if (IDOUGRP.eq.4) then
         ISPSM = IOUTLST + 4
      else if (IDOUGRP.eq.5.or.IDOUGRP.eq.6) then
         ISPSM = IOUTLST + 2
      end if
*
      do IREIM=1,2,1
         do IBU=1,2,1
            do IAB=1,2,1
               call fetch_cspexp_blk(CSPEXP,IBSPEXP,IDIM(1),IDIM(1),
     &                               IREIM,IBU,ISPSM,IAB,
     &                               CMOAO,LROW,LCOL)
               ILENBLK = LROW*LCOL
*
*  Determine whether alpha preceeds beta or vice versa
               call evenodd(IEOLST,IOUTLST)
               if (IEOLST.eq.1) IPERMAB = 1
               if (IEOLST.eq.2) IPERMAB = 2
*
*  Determine active dimensions and offsets for vectomat
               IDIMACT = 1
               if (IOUTLST.ge.3) IDIMACT = 3
               if (IDOUGRP.eq.1) IDIMACT = 1
               if (IPERMAB.eq.1.and.IAB.eq.1) then
                  IROWOFF = 1
                  if (IBU.eq.1) ICOLOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.2) ICOLOFF = 1
               end if
               if (IPERMAB.eq.1.and.IAB.eq.2) then
                  IROWOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.1) ICOLOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.2) ICOLOFF = 1
               end if
               if (IPERMAB.eq.2.and.IAB.eq.1) then
                  IROWOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.2) ICOLOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.1) ICOLOFF = 1
               end if
               if (IPERMAB.eq.2.and.IAB.eq.2) then
                  IROWOFF = 1
                  if (IBU.eq.2) ICOLOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.1) ICOLOFF = 1
               end if
*
               if (NTEST.ge.10) then
                  write(6,*) 'Testing before vectomat ...'
                  write(6,*) 'IREIM =  ',IREIM
                  write(6,*) 'IBU =  ',IBU
                  write(6,*) 'IAB =  ',IAB
                  write(6,*) 'IROWOFF =  ',IROWOFF
                  write(6,*) 'ICOLOFF =  ',ICOLOFF
               end if
*
               call vectomat(CMOAO,COEF(1,1,IREIM),ILENBLK,
     &                       ITOTDIM,ITOTDIM,
     &                       LROW,LCOL,IROWOFF,ICOLOFF)
               if (NTEST.gt.2) then
                  write(6,*) 'Checking coefficients ...'
                  call wrtmat(COEF(1,1,IREIM),ITOTDIM,ITOTDIM,
     &                        ITOTDIM,ITOTDIM)
               end if
            end do
         end do
      end do
*
      if (NTEST.gt.2) then
         write(6,*) 'Full coefficient matrices: real  imaginary'
         call wrtmat(COEF(1,1,1),ITOTDIM,ITOTDIM,ITOTDIM,ITOTDIM)
         call wrtmat(COEF(1,1,2),ITOTDIM,ITOTDIM,ITOTDIM,ITOTDIM)
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Import integrals from MOLCAS/SEWARD environment
*  (Driver routine)
*
      subroutine rdint_molcas(ONEINT,FNONEAO,LU2,NATOMS,IFOCKEX)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "symm.inc"
#include "onel.inc"
#include "input.inc"
*
      character*6 FNONEAO
*
      dimension ONEINT(ITRIASUM)
*
      call rd1int(ONEINT,FNONEAO,LU2,NATOMS,IFOCKEX)
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  read integrals from direct access file fort.2
*
      subroutine rdint_molecu(AOTRISM,HELPMAT,
     &                 LU2,LU9,LUSGI1,LENGTH,IHELEN)
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "symm.inc"
#include "onel.inc"
#include "input.inc"
*
      character*8 CAR(6)
      integer*4 IDUM,IDUMMY
*
* LUCIA common block for nuclear repulsion energy
      COMMON/CECORER/ECORE
*
      dimension AOTRISM(LENGTH)
      dimension TITLE(18)
      dimension IDUMMY(MXNSPIR)
      dimension HELPMAT(IHELEN)
*
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRINT)
*
* Calculate length of integral vector
      ILEN = 0
      do I=1,NIRR_PN,1
         ILEN = ILEN + ((ITOTDM(I)**2 + ITOTDM(I)) /2)
      end do
      write(6,*) 'Number of sf integrals is    ',ILEN
      if (LENGTH.ne.ILEN) then
         write(6,*) 'Problem in rdint'
         write(6,*) 'Dimensioning incorrect'
         write(6,*) 'ILEN is  ',ILEN
         write(6,*) 'Array length is  ',LENGTH
         Call Abend2( 'I will stop' )
      end if
*
* Refers to IBM architectures at the moment ...
* Read integrals and write out for porting to SGI
* ('O' stands for the type 'octane')
      if (ARCHI.ne.'SGO') then
         IADD2 = 0
         if (NTEST.ge.3) then
            write(6,*) 'IADD2 is    ',IADD2
            write(6,*) 'Checking address vector IAD2:'
            do I=1,3,1
               write(6,*) IAD2(I)
            end do
         end if
         IADD2 = IAD2(2)
*
* Initializing ...
         do I=1,LENGTH,1
            AOTRISM(I) = 0.0
         end do
*
* Read sf integrals from file ...
*
* Write ints to plain file for SGI run ...
         rewind LUSGI1
         if (NTEST.ge.3) then
            do I=1,ILEN,1
               write(LUSGI1,*) AOTRISM(I)
            end do
         end if
         rewind LUSGI1
      else if (ARCHI.eq.'SGO') then
         rewind LUSGI1
         do I=1,ILEN,1
            read(LUSGI1,*) AOTRISM(I)
         end do
      else
         write(6,*) 'rdint: Invalid architecture given. I will stop.'
         Call Abend('Quitting')
      end if
*
      if (NTEST.ge.2) then
         do I=1,ILEN,1
            write(6,*) 'AO INT ',I,AOTRISM(I)
         end do
      end if
*
*  If DOUGRP is C2H, then integrals have to be reordered into
*  double group ordering of symmetries, i.e.
*  (presuming that group generators Z and XY have been used)
*
*     MOLECULE      PICASSO     Order in PICASSO
*     --------      -------     ----------------
*       A_g         Gamma1+           1
*       A_u         Gamma1-           3
*       B_u         Gamma2-           4
*       B_g         Gamma2+           2
*
* ITOTDM are point group dimensions and given in LUCIAREL ordering!
*
      if (DOUGRP.eq.4) then
         IAUSTA = (ITOTDM(1)**2 + ITOTDM(1)) / 2 + 1
         IBUEND = IAUSTA + (ITOTDM(3)**2 + ITOTDM(3)) / 2
     &                   + (ITOTDM(4)**2 + ITOTDM(4)) / 2 - 1
         IBGSTA = IBUEND + 1
         IBGEND = IBGSTA + (ITOTDM(2)**2 + ITOTDM(2)) / 2 - 1
         if (IBGEND.ne.ILEN) then
           write(6,*) 'IBGEND is  ',IBGEND
           write(6,*) 'ILEN is  ',ILEN
           write(6,*) 'Problemo in reordering of C2h integrals.'
           Call Abend2( 'Quitting.' )
         end if
         IBGLEN = (ITOTDM(2)**2 + ITOTDM(2)) / 2
         if (NTEST.ge.3) then
            write(6,*) 'IAUSTA = ',IAUSTA
            write(6,*) 'IBUEND = ',IBUEND
            write(6,*) 'IBGSTA = ',IBGSTA
            write(6,*) 'IBGEND = ',IBGEND
            write(6,*) 'IBGLEN = ',IBGLEN
         end if
*
         J = 1
         do I=IBGSTA,IBGEND,1
            HELPMAT(J) = AOTRISM(I)
            J = J + 1
         end do
         KLAST = ILEN
         do I=IBUEND,IAUSTA,-1
            AOTRISM(KLAST) = AOTRISM(I)
            KLAST = KLAST - 1
         end do
         IN = 1
         do I=IAUSTA,IAUSTA+IBGLEN-1,1
            AOTRISM(I) = HELPMAT(IN)
            IN = IN + 1
         end do
*
         if (NTEST.ge.3) then
            write(6,*) 'Checking reordered AO integrals in C2H:'
            do I=1,ILEN,1
               write(6,*) 'AOTRISM(',I,')= ',AOTRISM(I)
            end do
         end if
      end if
*
*  Nuclear repulsion energy
*  Stored on common block 'onel.inc'
      rewind LU9
      read(LU9)
      read(LU9)
      read(LU9) TITLE,CAR,IDUM,(IDUMMY(I),I=1,IDUM),POTNUC
C?    write(6,*) IDUM
C?    do I=1,IDUM,1
C?       write(6,*) IDUMMY(I)
C?    end do
      write(6,*) 'Nuclear repulsion energy:  ',POTNUC
*
* Pass nuclear repulsion energy to LUCIA ...
      ECORE = POTNUC
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine rdmole2(SWEINT,NORBMOL,NBASMOL,NFRO,IGLSUM,XLUCINT)
*
* Sweden integral list        => Lucia integral list
* ^ Read from file to SWEINT     ^ Saved in array(sic) XLUCINT
*
*
*  Patch from tradpt (Christel Marian)
*  Timo Fleig
*
      implicit real*8 (A-H,O,P,R-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
#include "clunit.inc"
#include "parluci.h"
*
      parameter (LOCLEN=666)
*. Input
      dimension NORBMOL(NIRR_PN),NBASMOL(NIRR_PN)
      dimension SWEINT(*)
      dimension IOFFSET(LOCLEN)
      dimension NFRO(MXPIRR),QNTIL(MXPIRR),QMJ(MXPIRR),NIT(LOCLEN),
     &          QKJ(MXPIRR),QLJ(MXPIRR),QNJ(MXPIRR),QLG(MXPIRR),
     &          QISYM(MXPIRR),IDEKS(MXPORB)
*. Output
      DIMENSION XLUCINT(*)
*
      NTEST = 500
*
*  Set up triangular index array
      II = 0
      do I=1,MXPORB,1
        IDEKS(I) = II
        II = II + I
      end do
*
*  If DOUGRP is C2H, then some reordering has to be done.
*  First: Rename dimensions of symmetries
*         Go from LUCIAREL to MOLECULE ordering (Ag,Au,Bu,Bg)
*  Second: Reorder integrals according to LUCIAREL symmetry ordering
*          (done below)
*
      if (DOUGRP.eq.4) then
         IREOTEMP4 = NORBMOL(4)
         IREOTEMP3 = NORBMOL(3)
         IREOTEMP2 = NORBMOL(2)
         NORBMOL(2) = IREOTEMP3
         NORBMOL(3) = IREOTEMP4
         NORBMOL(4) = IREOTEMP2
*
         IREOTEMP4 = NBASMOL(4)
         IREOTEMP3 = NBASMOL(3)
         IREOTEMP2 = NBASMOL(2)
         NBASMOL(2) = IREOTEMP3
         NBASMOL(3) = IREOTEMP4
         NBASMOL(4) = IREOTEMP2
*
         IREOTEMP4 = NFRO(4)
         IREOTEMP3 = NFRO(3)
         IREOTEMP2 = NFRO(2)
         NFRO(2) = IREOTEMP3
         NFRO(3) = IREOTEMP4
         NFRO(4) = IREOTEMP2
      end if
*
*  CALCULATE VALUES FOR 1 ST OUTPUT-RECORD, PACK THEM TO OCCUPIED
*  IRREPS
*
      NOCC=0
      NBOX=0
      NCO=0
      LCT=0
      QNTIL(1)=0
      IORBS=NBASMOL(1)
      do 30 IS=1,NIRR_PN
        IORBS=IORBS+NBASMOL(IS)
        NCV=NFRO(IS)+NORBMOL(IS)
        if(NCV.EQ.0) goto 20
        NOCC=NOCC+1
*
*  # MOs - DELETED / OCCUPIED SYMMETRY
*
        QMJ(NOCC)=NCV
*
*  # FROZEN MO'S  / OCCUPIED SYMMETRY
*
        QKJ(NOCC)=NFRO(IS)
        NCO=NCO+NFRO(IS)
*
*  # VALENCE+VIRTUAL MO'S  / OCCUPIED SYMMETRY
*
        NO=NORBMOL(IS)
        QLJ(NOCC)=NO
        NBOX=NBOX+NO
*
*  # SOs / OCCUPIED SYMMETRY
*
        QNJ(NOCC)=NBASMOL(IS)
*
*  SYMMETRY-PACKING LABELS
*
        QISYM(NOCC)=IS
*
*  LAST MO IN EACH SYMMETRY
*
20      if(IS.NE.NIRR_PN) QNTIL(IS+1)=QNTIL(IS)+NO
30    CONTINUE
*
*  CALCULATE NUMBER OF INTEGRALS (NIGL) FOR EACH SYMMETRY-BLOCK
*  AND STORE 0-ADRESSES IN NIT
*
      NIT(1) = 0
      IOFFSET(1) = 1
      INIT = 1
      NPQSUM = 0
      IGLSUM = 0
*
*
* Calculate total number of integrals on record
* Determine old offsets for MOLECULE
      do ISP=1,NIRR_PN,1
         NSP=QISYM(ISP)
         if (NSP.eq.0) then
            write(LUWRT,*) 'Number of symm for P is  ',NSP
            write(LUWRT,*) 'Skipping two-el. int. read.'
            goto 333
         end if
         NOP=NORBMOL(NSP)
         do ISQ=1,ISP,1
            NSQ = QISYM(ISQ)
            NOQ = NORBMOL(NSQ)
            ISPQ = IPGMULT(ISP,ISQ)
            do ISR=1,ISP,1
               NSR = QISYM(ISR)
               NOR = NORBMOL(NSR)
               ISPQR = IPGMULT(ISPQ,ISR)
               ISSM=ISR
               if (ISR.EQ.ISP) ISSM=ISQ
               do 40 ISS=1,ISSM,1
                  NSS=QISYM(ISS)
                  NOS=NORBMOL(NSS)
                  if(NSS.NE.ISPQR) goto 40
                  if(NOS*NOR*NOQ*NOP.eq.0) goto 40
*
*  NUMBER OF INTEGRALS IN THIS SYMMETRYBLOCK
*
                  if (NSP.EQ.NSQ) then
*
*  TYPES <PP|PP> AND <PP|QQ>
*
                    NPQ = IDEKS(NOP+1)
                    NRS = IDEKS(NOR+1)
                  else
*
*  TYPES <PQ|PQ> AND <PQ|RS>
*
                    NPQ = NOP*NOQ
                    NRS = NOR*NOS
                  end if
                  if (NSP.EQ.NSR) then
*
*  TYPES <PP|PP> AND <PQ|PQ>
*
                    NIGL = NPQ*(NPQ+1)/2
                  else
*
*  TYPES <PP|QQ> AND <PQ|RS>
*
                    NIGL = NPQ*NRS
                  end if
                  NIT(INIT) = NIGL
                  INIT = INIT + 1
                  IOFFSET(INIT) = IOFFSET(INIT-1) + NIGL
                  IGLSUM = IGLSUM + NIGL
                  if (NTEST.ge.100) then
                     write(LUWRT,601) NSP,NSQ,NSR,NSS,NIGL
                  end if
40             continue
            end do
         end do
* Number of non-vanishing symmetry blocks
         NSMBLK = INIT - 1
         write(LUWRT,*) '# Symm blocks: ',NSMBLK
333      continue
      end do
601   FORMAT(' NUMBER OF INTEGRALS IN SYMMETRYBLOCK '
     &             ,4I2,' IS : ',I10)
*
*  Restore LUCIAREL ordering
*
      if (DOUGRP.eq.4) then
         IREOTEMP4 = NORBMOL(4)
         IREOTEMP3 = NORBMOL(3)
         IREOTEMP2 = NORBMOL(2)
         NORBMOL(2) = IREOTEMP4
         NORBMOL(3) = IREOTEMP2
         NORBMOL(4) = IREOTEMP3
*
         IREOTEMP4 = NBASMOL(4)
         IREOTEMP3 = NBASMOL(3)
         IREOTEMP2 = NBASMOL(2)
         NBASMOL(2) = IREOTEMP4
         NBASMOL(3) = IREOTEMP2
         NBASMOL(4) = IREOTEMP3
*
         IREOTEMP4 = NFRO(4)
         IREOTEMP3 = NFRO(3)
         IREOTEMP2 = NFRO(2)
         NFRO(2) = IREOTEMP4
         NFRO(3) = IREOTEMP2
         NFRO(4) = IREOTEMP3
*
* Do the same stuff for LUCIAREL offsets
*
         NOCC=0
         NBOX=0
         NCO=0
         LCT=0
         QNTIL(1)=0
         IORBS=NBASMOL(1)
         do 50 IS=1,NIRR_PN
           IORBS=IORBS+NBASMOL(IS)
           NCV=NFRO(IS)+NORBMOL(IS)
           if(NCV.EQ.0) goto 60
           NOCC=NOCC+1
*
*  # MOs - DELETED / OCCUPIED SYMMETRY
*
           QMJ(NOCC)=NCV
*
*  # FROZEN MO'S  / OCCUPIED SYMMETRY
*
           QKJ(NOCC)=NFRO(IS)
           NCO=NCO+NFRO(IS)
*
*  # VALENCE+VIRTUAL MO'S  / OCCUPIED SYMMETRY
*
           NO=NORBMOL(IS)
           QLJ(NOCC)=NO
           NBOX=NBOX+NO
*
*  # SOs / OCCUPIED SYMMETRY
*
           QNJ(NOCC)=NBASMOL(IS)
*
*  SYMMETRY-PACKING LABELS
*
           QISYM(NOCC)=IS
*
*  LAST MO IN EACH SYMMETRY
*
60         if(IS.NE.NIRR_PN) QNTIL(IS+1)=QNTIL(IS)+NO
50       CONTINUE
* Determine new offsets for LUCIAREL
         IOFFSET(INIT) = 1
         NPQSUM = 0
         IGLSUM_L = 0
         do ISP=1,NIRR_PN,1
            NSP=QISYM(ISP)
            if (NSP.eq.0) then
               write(6,*) 'Number of symm for P is  ',NSP
               write(6,*) 'Skipping two-el. int. read.'
               goto 334
            end if
            NOP=NORBMOL(NSP)
            write(6,*) 'ISP =  ',ISP
            write(6,*) 'NOP =   ',NOP
            do ISQ=1,ISP,1
               NSQ = QISYM(ISQ)
               NOQ = NORBMOL(NSQ)
               write(6,*) 'ISQ =  ',ISQ
               write(6,*) 'NOQ =   ',NOQ
CTF            ISPQ = IPGMULT(ISP,ISQ)
               ISPQ = IPGMULT(NSP,NSQ)
               do ISR=1,ISP,1
                  NSR = QISYM(ISR)
                  NOR = NORBMOL(NSR)
                  write(6,*) 'ISR =  ',ISR
                  write(6,*) 'NOR =   ',NOR
CTF               ISPQR = IPGMULT(ISPQ,ISR)
                  ISPQR = IPGMULT(ISPQ,NSR)
                  ISSM=ISR
                  if (ISR.EQ.ISP) ISSM=ISQ
                  do 80 ISS=1,ISSM,1
                     NSS=QISYM(ISS)
                     NOS=NORBMOL(NSS)
                     write(6,*) 'ISS =  ',ISS
                     write(6,*) 'NOS =   ',NOS
                     if (NSS.ne.ISPQR) then
                       NIGL = 0
                       goto 80
                     end if
                     if (NOS*NOR*NOQ*NOP.eq.0) then
                       NIGL = 0
                       goto 80
                     end if
*
*  NUMBER OF INTEGRALS IN THIS SYMMETRYBLOCK
*
                     if (NSP.EQ.NSQ) then
*
*  TYPES <PP|PP> AND <PP|QQ>
*
                       NPQ = IDEKS(NOP+1)
                       NRS = IDEKS(NOR+1)
                     else
*
*  TYPES <PQ|PQ> AND <PQ|RS>
*
                       NPQ = NOP*NOQ
                       NRS = NOR*NOS
                     end if
                     if (NSP.EQ.NSR) then
*
*  TYPES <PP|PP> AND <PQ|PQ>
*
                       NIGL = NPQ*(NPQ+1)/2
                     else
*
*  TYPES <PP|QQ> AND <PQ|RS>
*
                       NIGL = NPQ*NRS
                     end if
                     NIT(INIT) = NIGL
80                continue
                  IOFFSET(INIT+1) = IOFFSET(INIT) + NIGL
                  INIT = INIT + 1
                  IGLSUM_L = IGLSUM_L + NIGL
                  if (NTEST.ge.100) then
                     write(6,*) NSP,NSQ,NSR,NSS,NIGL
                     write(6,*) 'IGLSUM_L =  ',IGLSUM_L
                  end if
               end do
            end do
334         continue
         end do
      end if
*
*--------------------------------------------
* Read total integral record from file LU2INT
*--------------------------------------------
      IAD13 = 0
      LUSGITE = LUPRP
      if (ARCHI.eq.'RS6') then
*
*  Write for SGI read
         if (IPRINT.ge.2) then
            rewind LUSGITE
            write(6,*) 'writing transformed 2-el. integrals'
            do I=1,IGLSUM,1
               write(6,*) I,SWEINT(I)
               write(LUSGITE,*) SWEINT(I)
            end do
         end if
*
      else if (ARCHI.eq.'SGO') then
*
* I'm running on an SGI. Please give me a file ...
         rewind LUSGITE
         open(unit=LUSGITE,file='TWOEL',
     &        status='OLD',form='FORMATTED')
         do I=1,IGLSUM,1
            read(LUSGITE,*) SWEINT(I)
         end do
         close(unit=LUSGITE)
      else
         write(6,*) 'Architecture not implemented.'
         write(6,*) 'I am quitting!'
         Call Abend('Quitting')
      end if
*
*  Checking integrals previously read from LU 55 ...
      if (NTEST.ge.100) then
         write(6,*) 'Transformed spin-free 2-el. integrals:'
         do I=1,IGLSUM,1
            write(6,'(1G18.12)') SWEINT(I)
         end do
      end if
*
      if (DOUGRP.eq.4) then
* Do reordering of integrals if DOUGRP is C2H
         ILOCLEN = LOCLEN
         call reotwin_c2h(SWEINT(1),XLUCINT(1),
     &                    IOFFSET(1),IOFFSET(NSMBLK+1),NIT(1),
     &                    NIT(NSMBLK+1),IGLSUM,IGLSUM_L,
     &                    NSMBLK,NSMBLK_L,ILOCLEN)
      end if
*
* Now loop over symmetries again and pick relevant
* integral block from array SWEINT.
* Call reordering routine for internal LUCIA(REL)
* integral order.
*
*  P
*   Q
*    R
*     S
*
      INIT = 1
      IPOFF = 1
      do ISP=1,NIRR_PN,1
         NSP=QISYM(ISP)
         if (NSP.eq.0) then
            write(6,*) 'Number of symm for P is  ',NSP
            write(6,*) 'Skipping two-el. int. read.'
            goto 444
         end if
         NOP=NORBMOL(NSP)
         IQOFF = 1
         do ISQ=1,ISP,1
            NSQ = QISYM(ISQ)
            NOQ = NORBMOL(NSQ)
            ISPQ = IPGMULT(ISP,ISQ)
            IROFF = 1
            do ISR=1,ISP,1
               NSR = QISYM(ISR)
               NOR = NORBMOL(NSR)
               ISPQR = IPGMULT(ISPQ,ISR)
               ISSM=ISR
               if (ISR.EQ.ISP) ISSM=ISQ
               ISOFF = 1
               do 100 ISS=1,ISSM,1
                  NSS=QISYM(ISS)
                  NOS=NORBMOL(NSS)
                  if(NSS.NE.ISPQR) goto 49
                  if(NOS*NOR*NOQ*NOP.eq.0) goto 49
                  NINT = NIT(INIT)
*. offset in MOLECULE integrals to given symmetry block
                  IOFF = IOFFSET(INIT)
                  if (DOUGRP.eq.4) IOFF = IOFFSET(INIT+NSMBLK)
                  call reomoluc(SWEINT(IOFF),ISP,ISQ,ISR,ISS,
     &                          IPOFF,IQOFF,IROFF,ISOFF,
     &                          NOP,NOQ,NOR,NOS,
     &                          NIT(INIT+NSMBLK),XLUCINT)
*
                  INIT = INIT + 1
49                continue
                  ISOFF = ISOFF + NOS
100            continue
               IROFF = IROFF  + NOR
            end do
            IQOFF = IQOFF + NOQ
         end do
         IPOFF = IPOFF + NOP
444      continue
      end do
*
      write(6,*) 'IGLSUM (Number of ints) = ',IGLSUM
      INTNUM = IGLSUM
*
*  Check whether last integral is zero ...
      SWEABS = ABS(SWEINT(IGLSUM))
      CHECK = SWEABS - 1E-9
      if (CHECK.lt.0) then
        write(6,*) 'Problem with two-electron integral file'
        write(6,*) 'LASTABS is ',LASTABS
        write(6,*) 'CHECK is ',CHECK
        write(6,*) 'SWEINT(',IGLSUM,') is  ',SWEINT(IGLSUM)
        Call Abend2( 'Quitting in rdmole2' )
      end if
*
* Check integrals
      if (NTEST.ge.100) then
        write(6,*) 'Integrals on SWEINT after rdmole2: '
        write(6,*) 'Transformed spin-free 2-el. integrals:'
        do I=1,IGLSUM,1
          write(6,'(1G18.12)') SWEINT(I)
        end do
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine rdone(IDUM,IDUM1,DUM,IDUM2,IDUM3,IDUM4)
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Import one-electron integrals from DIRAC environment
*
      subroutine rdone_dirac_rel(ONEINT,LUONE,DIRACONE,MORB,NDIM,
     &                           ISPINFREE,IBOSYM,IDUMAR,DUMAR)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "input.inc"
*
      character*14 CHA(8)
      common/CECORER/ECORE
*
      character*7 DIRACONE
      logical BREIT
*
      dimension IDUMAR(*),IBOSYM(*),NSTR(2)
      dimension ONEINT((2*NDIM)**2*2),DUMAR(*)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      open(LUONE,FILE=DIRACONE,FORM='UNFORMATTED')
      read(LUONE) MORB,BREIT,ECORE_DIR,NFSYM,NZ
      read(LUONE) NSYMRP,(CHA(I),I=1,NSYMRP),(IDUMAR(I),I=1,NSYMRP),
     &                                         (NSTR(I),I=1,NFSYM)
      read(LUONE)
      read(LUONE)
      read(LUONE) (IDUMAR(ID),IDUMAR(ID),DUMAR(ID),ID=1,MORB,1),
     &            (IBOSYM(I),I=1,MORB,1)
      if (NTEST.ge.3) then
        write(6,*) 
        write(6,*) 'Testing read MRCONEE'
        write(6,*) 'MORB = ',MORB
        write(6,*) 'BREIT = ',BREIT
        write(6,*) 'ECORE_DIR = ',ECORE_DIR
        write(6,*) 'Number of spinors per fermion symmetry:'
        write(6,'(2I5)') (NSTR(I),I=1,NFSYM)
        write(6,*) 'Read boson symmetries of spinors:'
        do IB=1,MORB,1
          write(6,'(A11,I5,A3,I12)') 'Spinor no. ',IB,' : ',IBOSYM(IB)
        end do
      end if
*
      if (2*NDIM.ne.MORB) then
        write(6,*) 'Orbital dimensions of '
        write(6,*) 'DIRAC and input do not match.'
        write(6,*) 'DIRAC : ',MORB,' spinors'
        write(6,*) 'Input : ',NDIM,' orbitals'
        call abend2('Quitting.')
      end if
*
* Read the integrals
      NDIMSQ = MORB*MORB
      if (ISPINFREE.eq.1) then
        read(LUONE) (ONEINT(I),TDUM,I=1,MORB*MORB)
      else if (ISPINFREE.eq.0) then
        read(LUONE) (ONEINT(I),ONEINT(I+NDIMSQ),I=1,MORB*MORB)
*
* Set IBOSYM array so that U->alpha and B->beta spinors
*              To get correct Fermion symmetry labeling in
*              LUCIAREL we need to specify same(!) gerade boson
*              irrep for all gerade orbitals and same ungerade
*              boson irrep for all ungerade orbitals (note: which
*              gerade or ungerade boson irrep we choose doesn't
*              matter).
        IOFFFS = 0
        IRROFF = 0
        do IFSYM = 1,NFSYM,1
          do ISP = 1+IOFFFS,NSTR(IFSYM)+IOFFFS,1
            IBOSYM(ISP) = 0 + IRROFF
          end do
*
          if (NTEST.ge.3) then
            write(6,*) ' Reset IBOSYM array, ferm symm. :',IFSYM
            write(6,*) (IBOSYM(I),I=1+IOFFFS,NSTR(IFSYM)+IOFFFS,1)
          end if
*
          IOFFFS = NSTR(IFSYM)
          IRROFF = 2
          if (NZ.eq.1) IRROFF = 1
        end do
*
      end if
*
      if (NTEST.ge.5) then
        write(6,*)
        MORB_L = MORB
        write(6,*) 'Real One-electron ints from dirac:'
        call wrtmat(ONEINT(1),MORB_L,MORB_L,MORB_L,MORB_L)
        if (ISPINFREE.eq.0) then
          write(6,*) 'Imaginary One-electron ints from dirac:'
          call wrtmat(ONEINT(1+NDIMSQ),MORB_L,MORB_L,MORB_L,MORB_L)
        end if
      end if  
*
      ECORE = ECORE_DIR
* 
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Read spin-orbit integrals
*
      subroutine rdsoint_molcas(SOINXLINE,SOINYLINE,
     &                   SOINZLINE,NTRIAEL,
     &                   LU69,LU70,LU71,
     &                   LUSGI2,LUSGI3,LUSGI4)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
      character*8 CHARTEST
      character*30 INTTYPE
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
*
      parameter (ALPHA = 137.036)
*
      dimension SOINZLINE(NTRIAEL),SOINXLINE(NTRIAEL),
     &          SOINYLINE(NTRIAEL)
*
      NTESTL = 0
      if (NTESTL.eq.0) then
        write(6,*) 'No testing in rdsoint_molcas.'
      end if
      NTEST = max(NTESTL,IPRINT)
*
      if (ARCHI.eq.'SGO') then
         open(unit=LU69,file='AOPROPER_MF_SYM',
     &        status='OLD',form='UNFORMATTED')
         read(LU69)
         read(LU69) (SOINXLINE(I),I=1,NTRIAEL,1)
         read(LU69)
         read(LU69) (SOINYLINE(I),I=1,NTRIAEL,1)
         read(LU69)
         read(LU69) (SOINZLINE(I),I=1,NTRIAEL,1)
         close(unit=LU69)
*
* Change to atomic units:
         do I=1,NTRIAEL,1
           SOINXLINE(I) = SOINXLINE(I) * (1/(4*ALPHA**2))
           SOINYLINE(I) = SOINYLINE(I) * (1/(4*ALPHA**2))
           SOINZLINE(I) = SOINZLINE(I) * (1/(4*ALPHA**2))
         end do
*
         if (NTEST.ge.1) then
            write(6,*) 'writing SO - X integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINXLINE(I)
            end do
            write(6,*) 'writing SO - Y integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINYLINE(I)
            end do
            write(6,*) 'writing SO - Z integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINZLINE(I)
            end do
         end if
*
         write(6,*)
         write(6,*) 'SO-MF integrals read from file '
         write(6,*) 'AOPROPER_MF_SYM'
      else if (ARCHI.ne.'SGO') then
        write(6,*) 'Architecture set to  ',ARCHI
        write(6,*) 'I do not know this yet (under MOLCAS).'
        Call Abend2( 'Quitting.' )
      else
         write(6,*) 'rdsoint: Invalid architecture given. I will stop.'
         Call Abend('Quitting')
      end if
*
      write(6,*)
      write(6,'(A,I8)')
     &   'Number of SO integrals on full triangle: ',NTRIAEL
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Read spin-orbit integrals
*
      subroutine rdsoint_molcas5(SOINXLINE,SOINYLINE,
     &                           SOINZLINE,HELPMAT,IDIMLOC,
     &                           NTRIAEL,FNONEAO,LUONEAO)
*
*--------------------------------------------------------
* For compatibility with new storage of AMFI integrals
* under MOLCAS versions 5 and later.
*
*   Timo Fleig, summer 2000
*
* Generalized for use with external Fock matrices
*   Timo Fleig, 28.05.2001
*--------------------------------------------------------
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
      integer TRIANGLE
      character*6 FNONEAO
      character*8 CHARTEST,LABEL
      character*30 INTTYPE
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
*
      parameter (ALPHA = 137.036)
*
      dimension IDIMLOC(NIRR_PN)
      dimension SOINZLINE(NTRIAEL),SOINXLINE(NTRIAEL),
     &          SOINYLINE(NTRIAEL),HELPMAT(NTRIAEL)
      dimension TRIANGLE(2*MXNSPIPDMGRP)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
*  Define help array
      INRS = 0
      do IR=1,NIRR_PN,1
        INRS = INRS + IDIMLOC(IR)
      end do
      do IT=1,INRS,1
        TRIANGLE(IT) = 0
        do ITIN=1,IT,1
          TRIANGLE(IT) = TRIANGLE(IT) + ITIN
        end do
      end do
*
*************************************
*  Read from MOLCAS ONEINT file
*************************************
*
      Call qEnter('rdsoint')
*  Open file
      if (NTEST.ge.1) write(6,'(A14,I3)') 'OPENING LUONE ',LUONEAO
      iRc=-1
      iOpt=0
      Call OpnOne(iRC,iOpt,FNONEAO,LUONEAO)
      if (IRC.ne.0) goto 990
*
      do ICOMP=1,3,1
        IRC = -1
        IOPT = 1                ! get size NSIZ only
        LABEL = 'AMFI    '
        call rdone(IRC,IOPT,LABEL,ICOMP,NSIZ,ISYLBL)
        if (NTEST.ge.1) then
          write(6,*) 'Reading ',NSIZ,' SO integrals.'
        end if
C       if (IRC.ne.0) goto 991
*
        IOPT = 0
        call rdone(IRC,IOPT,LABEL,ICOMP,HELPMAT,ISYLBL)
C       if (IRC.ne.0) goto 991
        if (NTEST.ge.2) then
          write(6,*) 'Block of symmetry SO integrals read:'
          write(6,*) 'Symmetry label ',ISYLBL
          call wrtmat(HELPMAT,NSIZ,1,NSIZ,1)
        end if
*
        IN = 1
        NSKIP = 0
        IOUT = 1
        STRFAC = -1.0
* "Strange factor": Compared to the old AMFI, all integrals need to
* be multiplied by -1, in analogy to RASSI code comments. The
* reason is unknown!
* There remains a slight numerical difference (ca. 1*10^-11 E_H),
* but this is probably due to slightly different unit translation
* factors and will be ignored.
        if (ISYLBL.eq.1) then
          do ISYM=1,NIRR_PN,1
            do JSYM=1,ISYM-1,1
              NSKIP = NSKIP + IDIMLOC(JSYM)
            end do
            do IROW=1,IDIMLOC(ISYM),1
              IOUT = IOUT + NSKIP
              do ICOL=1,IROW,1
                if (ICOMP.eq.1) then
                  SOINXLINE(IOUT) = STRFAC*HELPMAT(IN)
                else if (ICOMP.eq.2) then
                  SOINYLINE(IOUT) = STRFAC*HELPMAT(IN)
                else if (ICOMP.eq.3) then
                  SOINZLINE(IOUT) = STRFAC*HELPMAT(IN)
                end if
                if (NTEST.ge.5) then
             write(6,'(A4,I6,A11,I6)') 'Int ',IN,' copied to ',IOUT
                end if
                IN = IN + 1
                IOUT = IOUT + 1
              end do
            end do
            NSKIP = 0
          end do
          if (IN.ne.NSIZ+1) then
            write(6,*) 'Number of distributed SO ints'
            write(6,*) 'does not match size of input array.'
            call abend2('Problem in rdsoint_molcas.')
          end if
        else if (ISYLBL.eq.2) then
          ICT = 0
          IDIMSUM = 0
          IOFF = 0
          do ISYM=1,NIRR_PN,2
            do JSYM=1,ISYM-1,1
              NSKIP = NSKIP + IDIMLOC(JSYM)
            end do
            if (NTEST.ge.5) write(6,*) 'NSKIP   ',NSKIP
            do KSYM=1,ISYM,1
              IDIMSUM = IDIMSUM + IDIMLOC(KSYM)
            end do
            if (NTEST.ge.5) write(6,*) 'IDIMSUM   ',IDIMSUM
            IOUT = IOUT + NSKIP
            if (NTEST.ge.5) write(6,*) 'IOUT before work ',IOUT
            if (IDIMLOC(ISYM+1).eq.0) goto 111
            do ICOL=1,IDIMLOC(ISYM),1
              do IROW=1,IDIMLOC(ISYM+1),1
                IOFF = TRIANGLE(IDIMSUM+IROW-1)
                if (NTEST.ge.5) write(6,*) 'IOFF  ',IOFF
                if (ICOMP.eq.1) then
                  SOINXLINE(IOFF+IOUT) = STRFAC*HELPMAT(IN)
                else if (ICOMP.eq.2) then
                  SOINYLINE(IOFF+IOUT) = STRFAC*HELPMAT(IN)
                end if
                if (NTEST.ge.5) then
          write(6,'(A4,I6,A11,I6)') 'Int ',IN,' copied to ',IOFF+IOUT
                end if
                IN = IN + 1
                ICT = ICT + 1
              end do
              ICT = 0
              IOUT = IOUT + 1
            end do
 111        continue
            IDIMSUM = 0
            NSKIP = 0
          end do
          if (IN.ne.NSIZ+1) then
            write(6,*) 'Number of distributed SO ints'
            write(6,*) 'does not match size of input array.'
            call abend2('Problem in rdsoint_molcas.')
          end if
        end if
      end do
*
      if (NTEST.ge.2) then
         write(6,*) 'writing SO - X integrals ... '
         do I=1,NTRIAEL,1
            write(6,*) I,SOINXLINE(I)
         end do
         write(6,*) 'writing SO - Y integrals ... '
         do I=1,NTRIAEL,1
            write(6,*) I,SOINYLINE(I)
         end do
         write(6,*) 'writing SO - Z integrals ... '
         do I=1,NTRIAEL,1
            write(6,*) I,SOINZLINE(I)
         end do
      end if
*
*  Close ONEINT file
      IRC = -1
      IOPT = 0
      call clsone(IRC,IOPT)
      If ( iRc.ne.0 ) Goto 992
      Call qExit('rdsoint')
*
      return
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
990   Call ErrMsg('RsoInt',001)
      Call Quit(16)
C991   Call StkA(0,OneLbl)
      Call ErrMsg('RsoInt',002)
      Call Quit(16)
992   Call ErrMsg('RsoInt',003)
      Call Quit(16)
993   Call ErrMsg('RsoInt',004)
      Call Quit(16)
C994   Call StkA(0,ComLbl)
      Call ErrMsg('RsoInt',005)
      Call Quit(16)
995   Call ErrMsg('RsoInt',006)
      Call Quit(16)
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Read spin-orbit integrals
*
      subroutine rdsoint_molecu(SOINXLINE,SOINYLINE,
     &                          SOINZLINE,
     &                          HELPMAT,NTRIAEL,
     &                          LU69,LU70,LU71,
     &                          LUSGI2,LUSGI3,LUSGI4)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
      character*8 CHARTEST
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
*
      dimension SOINZLINE(NTRIAEL),SOINXLINE(NTRIAEL),
     &          SOINYLINE(NTRIAEL)
      dimension HELPMAT(*)
      dimension IL1(4,4),IL2(10),IL3(10),IL4(10)
*
      if (ARCHI.ne.'SGO') then
         open(unit=LU69,file='SOMFINT_X',
     &        status='OLD',form='UNFORMATTED')
         read(LU69) (SOINXLINE(I),I=1,NTRIAEL,1)
         close(unit=LU69)
*
         if (IPRINT.ge.2) then
            rewind LUSGI2
            write(6,*) 'writing SO - X integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINXLINE(I)
               write(LUSGI2,*) SOINXLINE(I)
            end do
            rewind LUSGI2
         end if
*
         open(unit=LU70,file='SOMFINT_Y',
     &        status='OLD',form='UNFORMATTED')
         read(LU70) (SOINYLINE(I),I=1,NTRIAEL,1)
         close(unit=LU70)
*
         if (IPRINT.ge.2) then
            rewind LUSGI3
            write(6,*) 'writing SO - Y integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINYLINE(I)
               write(LUSGI3,*) SOINYLINE(I)
            end do
            rewind LUSGI3
         end if
*
         open(unit=LU71,file='SOMFINT_Z',
     &        status='OLD',form='UNFORMATTED')
         read(LU71) (SOINZLINE(I),I=1,NTRIAEL,1)
         close(unit=LU71)
*
         if (IPRINT.ge.2) then
            rewind LUSGI4
            write(6,*) 'writing SO - Z integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINZLINE(I)
               write(LUSGI4,*) SOINZLINE(I)
            end do
            rewind LUSGI4
         end if
      else if (ARCHI.eq.'SGO') then
         rewind LUSGI2
         rewind LUSGI3
         rewind LUSGI4
         do I=1,NTRIAEL,1
            read(LUSGI2,*) SOINXLINE(I)
            read(LUSGI3,*) SOINYLINE(I)
            read(LUSGI4,*) SOINZLINE(I)
         end do
         if (IPRINT.ge.2) then
            write(6,*) 'writing SO - X integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINXLINE(I)
            end do
            write(6,*) 'writing SO - Y integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINYLINE(I)
            end do
            write(6,*) 'writing SO - Z integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINZLINE(I)
            end do
         end if
      else
         write(6,*) 'rdsoint: Invalid architecture given. I will stop.'
         Call Abend('Quitting')
      end if
*
*  If DOUGRP is C2H, then integrals have to be reordered into
*  double group ordering of symmetries, i.e.
*  (presuming that group generators Z and XY have been used)
*
*     MOLECULE      PICASSO     Order in PICASSO
*     --------      -------     ----------------
*       A_g         Gamma1+           1
*       A_u         Gamma1-           3
*       B_u         Gamma2-           4
*       B_g         Gamma2+           2
*
*  This is too complicated, as the matrix is a full matrix,
*  not just symmetry non-vanishing blocks. So we proceed
*  differently: No reordering is done here. Instead, intsel
*  picks the non-vanishing 2-symmetry blocks appropriately.
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Read DIRAC-COULOMB two-electron integrals
*  (Later also BREIT integrals)
*  (GAUNT integrals can be included in COULOMB type ints.,
*   because they have the same symmetry.)
*
      subroutine rdtwo_dirac_rel(RKLR,RKLI,LUTWO,DIRACTWO,
     &                           NDIM,ISPINFREE,MORB,
     &                           KR,IOFF12,NUMINT,ITPLOC,INDK,INDL)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "input.inc"
*
      character*7 DIRACTWO
      character*10 DATEX,DATEXB,TIMEX*8,TIMEXB*8
      logical REALAR
c     logical BREIT
*
      dimension RKLR((2*NDIM)**4),RKLI((2*NDIM)**4)
      dimension KR(-NDIM:NDIM),INDK((2*NDIM)**4),INDL((2*NDIM)**4),
     &          IOFF12(MORB,MORB,4),NUMINT(MORB,MORB,4),
     &          ITPLOC(MORB,MORB)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      REALAR = .TRUE.
      TRESINT = 1.D-12
      NINT = 0
      KLS = 0
!radovan: m2 and m1 are not set
      KLF = M2
      MM = M1 * M2
*
      open (LUTWO,FILE=DIRACTWO,FORM='UNFORMATTED')
      rewind (LUTWO)
      read (LUTWO) datex,timex,nkr,
     & (kr(i),kr(-i),i=1,nkr)
*
      if (2*NKR.ne.MORB) then
        write(6,*) 'MRCONEE and MDCINT are incompatible.'
        stop 'Quitting.'
      end if
      write(6,*) 'Integrals have been calculated on the'
      write(6,'(A12,1X,A)') '            ',DATEX,' at ',TIMEX
      if (NTEST.ge.5) then
        write(6,*)
        write(6,*) 'Info from MDCINT file'
        write(6,*) 'NKR is ',NKR
        write(6,*) 'KR array:'
        do I=1,NKR,1
          write(6,'(2I8)') KR(-I),KR(I)
        end do
      end if
C
C     -------------------------
C     READ IN COULOMB INTEGRALS
C     -------------------------
C
* Initialize
      ISYMOFF = 0
      ITYPE = 1
      IKR_SAVE = 0
      JKR_SAVE = 0
      do M=1,MORB,1
        do N=1,MORB,1
          ITPLOC(M,N) = 0
          do IT=1,4,1
            NUMINT(M,N,IT) = 0
          end do
        end do
      end do
*
      rewind (LUTWO)
      read (LUTWO)
*
* Read a block of integrals
  10  continue
      read (LUTWO) IKR,JKR,NZ,
     &             (INDK(INZ+ISYMOFF),INDL(INZ+ISYMOFF),INZ=1,NZ,1),
     &             (RKLR(INZ+ISYMOFF),RKLI(INZ+ISYMOFF),INZ=1,NZ,1)
*
* End of file 
* (FIXME: No multiple files implemented, yet!)
      if (IKR.eq.0) then
        write(6,*) 'End of file MDCINT.'
        goto 20
      end if
*
      ITPLOC(KR(IKR),KR(JKR)) = ITPLOC(KR(IKR),KR(JKR)) + 1
      IOFF12(KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR))) = ISYMOFF
      NUMINT(KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR))) = NZ
*
      if (NTEST.ge.10) then
        write(6,'(A,2I6)') 'IKR,JKR : ',IKR,JKR
        write(6,'(A,3I4,A,I6)') 
     &        'IOFF12(',KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR)),') = ',
     &         IOFF12(KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR)))
        write(6,'(A,3I4,A,I6)')
     &        'NUMINT(',KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR)),') = ',
     &         NUMINT(KR(IKR),KR(JKR),ITPLOC(KR(IKR),KR(JKR)))
        write(6,'(A,1X,2I6)') 'Integral block: ',KR(IKR),KR(JKR)
        do INT=1,NZ,1
          write(6,'(A,2I8,1X,2F18.12)') 
     &      'Indices K,L ',KR(INDK(INT+ISYMOFF)),KR(INDL(INT+ISYMOFF)),
     &                     RKLR(INT+ISYMOFF),RKLI(INT+ISYMOFF)
        end do
      end if
*
      ISYMOFF = ISYMOFF + NZ
      IKR_SAVE = IKR
      JKR_SAVE = JKR
      goto 10
*
   20 if (NTEST.ge.5) then
        write(6,*)
        write(6,*) '===================='
        write(6,*) 'rdtwo_dirac speaking'
        write(6,*) '===================='
        write(6,*)
        write(6,*) ' Index pair, int type, offset, number of integrals:'
        do IKPAIR=1,-1,-2
          do II=IKPAIR,MORB*IKPAIR/2,IKPAIR
            do JKPAIR=1,-1,-2
              do IJ=JKPAIR,MORB*JKPAIR/2,JKPAIR
                do IT=1,ITPLOC(KR(II),KR(IJ)),1
                  if (NUMINT(KR(II),KR(IJ),IT).gt.0) then
                    write(6,'(1X,5I5)') KR(II),KR(IJ),IT,
     &                                  IOFF12(KR(II),KR(IJ),IT),
     &                                  NUMINT(KR(II),KR(IJ),IT)
                    do INT=1,NUMINT(KR(II),KR(IJ),IT),1
                      write(6,'(5X,2I4,2F18.10)') 
     &                 KR(INDK(INT+IOFF12(KR(II),KR(IJ),IT))),
     &                 KR(INDL(INT+IOFF12(KR(II),KR(IJ),IT))),
     &                 RKLR(INT+IOFF12(KR(II),KR(IJ),IT)),
     &                 RKLI(INT+IOFF12(KR(II),KR(IJ),IT))
                    end do
                  end if
                end do
              end do
            end do
          end do
        end do
      end if


c     DO INZ = 1, NZ
c        IF (ABS(RKLI(INZ)).GT.TRESINT) REALAR = .FALSE.
c     ENDDO
c     MINT = MINT + NZ
c     If (ikr.EQ.0) GoTo 40 ! End of File
c     i = indmo(kr(ikr))
c     j = indmo(kr(jkr))
c     itr = indmo(kr(-ikr))
c     jtr = indmo(kr(-jkr))
c     SignIJ = SIGN(1,ikr) * SIGN(1,jkr)
c     DO INZ = 1, NZ
c        kkr = indk(inz)
c        lkr = indl(inz)
c        k = indmo(kr(kkr))
c        l = indmo(kr(lkr))
c        ktr = indmo(kr(-kkr))
c        ltr = indmo(kr(-lkr))
c        SignKL = SIGN(1,kkr) * SIGN(1,lkr)
C--> Original integral
c        IF (I.NE.0.AND.J.NE.0.AND.K.NE.0.AND.L.NE.0) THEN
c        TWOR = RKLR(INZ)
c        TWOI = RKLI(INZ)
c        IJ = (I-1) * NORB + J
c        KL = (K-1) * NORB + L
c        JI = (J-1) * NORB + I
c        LK = (L-1) * NORB + K
c        IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
c           TIR(IJ+(KL-KLS-1)*M1) = TWOR
c           TII(IJ+(KL-KLS-1)*M1) = TWOI
c        ENDIF
c        IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
c           TIR(JI+(LK-KLS-1)*M1) = TWOR
c           TII(JI+(LK-KLS-1)*M1) = - TWOI
c        ENDIF
c        IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
c           TIR(KL+(IJ-KLS-1)*M1) = TWOR
c           TII(KL+(IJ-KLS-1)*M1) = TWOI
c        ENDIF
c        IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
c           TIR(LK+(JI-KLS-1)*M1) = TWOR
c           TII(LK+(JI-KLS-1)*M1) = - TWOI
c        ENDIF
c        ENDIF ! check on deleted spinors
c     ENDDO
C
C
c  40 IF (MINT.EQ.0) WRITE (6,1010) 'MDCINT'
C
c     IF (.NOT.BREIT) GOTO 61
C
C     -----------------------
C     READ IN BREIT INTEGRALS
C     -----------------------
C
c  50 read (mdbint,ERR=10011,END=10011) ikr,jkr,nz,
c    &             (indk(inz),indl(inz),inz=1,nz),
c    &             (rklr(inz),rkli(inz),inz=1,nz)
c     DO INZ = 1, NZ
c        IF (ABS(RKLI(INZ)).GT.TRESINT) REALAR = .FALSE.
c     ENDDO
c     MINTB = MINTB + NZ
c     If (ikr.EQ.0) GoTo 60 ! End of File
c     i = indmo(kr(ikr))
c     j = indmo(kr(jkr))
c     itr = indmo(kr(-ikr))
c     jtr = indmo(kr(-jkr))
c     SignIJ = - SIGN(1,ikr) * SIGN(1,jkr)
c     DO INZ = 1, NZ
c        kkr = indk(inz)
c        lkr = indl(inz)
c        k = indmo(kr(kkr))
c        l = indmo(kr(lkr))
c        ktr = indmo(kr(-kkr))
c        ltr = indmo(kr(-lkr))
c        SignKL = - SIGN(1,kkr) * SIGN(1,lkr)
c--> Original integral
c        IF (I.NE.0.AND.J.NE.0.AND.K.NE.0.AND.L.NE.0) THEN
c        TWOR = RKLR(INZ)
c        TWOI = RKLI(INZ)
c        IJ = (I-1) * NORB + J
c        KL = (K-1) * NORB + L
c        JI = (J-1) * NORB + I
c        LK = (L-1) * NORB + K
c        IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
c           TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
c           TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
c           TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
c           TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
c           TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
c           TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
c           TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
c           TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        ENDIF ! check on deleted spinors
c--> Time reversal left
c        IF (ITR.NE.0.AND.JTR.NE.0.AND.K.NE.0.AND.L.NE.0) THEN
c        TWOR = RKLR(INZ) * SignIJ
c        TWOI = RKLI(INZ) * SignIJ
c        IJ = (JTR-1) * NORB + ITR
c        KL = (K-1) * NORB + L
c        JI = (ITR-1) * NORB + JTR
c        LK = (L-1) * NORB + K
c        IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
c           TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
c           TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
c           TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
c           TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
c           TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
c           TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
c           TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
c           TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        ENDIF ! check on deleted spinors
c--> Time reversal right
c        IF (I.NE.0.AND.J.NE.0.AND.KTR.NE.0.AND.LTR.NE.0) THEN
c        TWOR = RKLR(INZ) * SignKL
c        TWOI = RKLI(INZ) * SignKL
c        IJ = (I-1) * NORB + J
c        KL = (LTR-1) * NORB + KTR
c        JI = (J-1) * NORB + I
c        LK = (KTR-1) * NORB + LTR
c        IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
c           TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
c           TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
c           TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
c           TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
c           TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
c           TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
c           TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
c           TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        ENDIF ! check on deleted spinors
c--> Time reversal both
c        IF (ITR.NE.0.AND.JTR.NE.0.AND.KTR.NE.0.AND.LTR.NE.0) THEN
c        TWOR = RKLR(INZ) * SignIJ * SignKL
c        TWOI = RKLI(INZ) * SignIJ * SignKL
c        IJ = (JTR-1) * NORB + ITR
c        KL = (LTR-1) * NORB + KTR
c        JI = (ITR-1) * NORB + JTR
c        LK = (KTR-1) * NORB + LTR
c        IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
c           TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
c           TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
c           TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
c           TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
c           TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
c           TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
c        ENDIF
c        IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
c           TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
c           TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
c        ENDIF
c        ENDIF ! check on deleted spinors
c     ENDDO
c
c     GOTO 50
c
c  60 IF (MINTB.EQ.0) WRITE (6,1010) 'MDBINT'
c
c     -----------------------------------------
c     ADD BREIT INTEGRALS AND COULOMB INTEGRALS
c     -----------------------------------------
c
c     DO IJKL = 1, (KLF-KLS) * M1
c        TIR(IJKL) = TIR(IJKL) + TIR(IJKL+MM)
c        TII(IJKL) = TII(IJKL) + TII(IJKL+MM)
c     ENDDO
c
   61 CONTINUE
C
C     ----------------------------------
C     WRITE 2-ELECTRON INTEGRALS TO FILE
C     ----------------------------------
C
c     IJKL = 0
c     DO 200 KL = KLS + 1, KLF
c     TIR (KL+(KL-KLS-1)*M1) = 0.5D0 * TIR(KL+(KL-KLS-1)*M1)
c     TII (KL+(KL-KLS-1)*M1) = 0.5D0 * TII(KL+(KL-KLS-1)*M1)
c     NKL = 0
c     DO 100 IJ = 1, M1
c        IJKL = IJKL + 1
c        IF (ABS(TIR(IJKL))+ABS(TII(IJKL)).GT.TRESINT) THEN
c           NKL = NKL + 1
c           INDK(NKL) = IJ
c           RKLR(NKL) = TIR(IJKL)
c           RKLI(NKL) = TII(IJKL)
c         ENDIF
c 100     CONTINUE
c     NINT = NINT + NKL
c     IF (REALAR) THEN
c        WRITE (MRCTWOE) KL,NKL,(INDK(IX),IX=1,NKL),(RKLR(IX),IX=1,NKL)
c     ELSE
c        WRITE (MRCTWOE) KL,NKL,(INDK(IX),IX=1,NKL),(RKLR(IX),IX=1,NKL),
c    &                   (RKLI(IX),IX=1,NKL)
c     ENDIF
c 200 CONTINUE
c     KLS = KLF
c     KLF = MIN0 (M1, KLF+M2)
c     IF (KLS.NE.M1) GOTO 10
c     IF (.NOT.BREIT) THEN
c        WRITE (6,1000) MINT,DATEX,TIMEX,NINT
c     ELSE
c        WRITE(6,1001) MINT,DATEX,TIMEX,MINTB,DATEXB,TIMEXB,NINT
c     ENDIF
c     IF (NINT.GT.0) THEN
c        WRITE (MRCTWOE) -NINT,0
c     ELSE
c        WRITE (MRCTWOE) -1,0
c     ENDIF
c     IF (REALAR) THEN
c        WRITE (6,1002)
c     ELSE
c        WRITE (6,1003)
c     ENDIF
C
c     CLOSE (MDCINT,STATUS='KEEP')
c     IF (BREIT) CLOSE (MDBINT,STATUS='KEEP')
      RETURN
 1000 FORMAT (
     &/I8,' Unique Coulomb integrals read from MDCINT (',A10,1X,A8,')'
     &/I8,' Coulomb integrals written to MRCTWOE')
 1001 FORMAT (
     &/I8,' Unique Coulomb integrals read from MDCINT (',A10,1X,A8,')'
     &/I8,' Unique Breit integrals read from MDBINT (',A10,1X,A8,')'
     &/I8,' Combined integrals written to MRCTWOE')
 1002 FORMAT (/' Real arithmetic will be used in EXPAND and TWOLIN')
 1003 FORMAT (/' Complex arithmetic will be used')
 1010 FORMAT (//' CAUTION : No non-zero integrals found on ',A//)
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Determine matching dimensions of boson functions,
* and reordering arrays for DIRAC environment.
*
      subroutine reo_info_dirac(ITOTDIM,NIRR_PN,SYMDIM,IBOSYM,SYDI,
     &                          DOUGRP,ISREO,INVREO,INVSYM,IDGTYP,
     &                          ISPINFREE,IPRNT)
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
*
      integer SYDI(NIRR_PN),SYMDIM(NIRR_PN),DOUGRP
      dimension ISREO(2*ITOTDIM),INVREO(-ITOTDIM:ITOTDIM),
     &          IBOSYM(2*ITOTDIM)
*
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
* Check for matching dimensions of boson symmetry functions
      if (ISPINFREE.eq.1) then
        do IRR=1,NIRR_PN,1
          IBCT = 0
          IOFFBS = 0
          do ISP=1,ITOTDIM,1
            if (IBOSYM(ISP+IOFFBS)+1.eq.IRR) IBCT = IBCT + 1
            if (INVSYM.eq.1.and.ISP.eq.SYMDIM(1)/2) IOFFBS = SYMDIM(1)/2
          end do
*
          if (IBCT.ne.SYDI(IRR)) then
            write(6,*)
            write(6,*) 'Problem with IBSPI array.'
            write(6,*) 'Boson symmetry dimensions do not match'
            write(6,*) 'with GAS matrix input.'
            write(6,*) 'Symmetry ',IRR
            write(6,*) 'No. of input functions ',SYDI(IRR)
            write(6,*) 'No. of functions on IBSPI array ',IBCT
            call abend('Quitting.')
          end if
        end do
        write(6,*)
        write(6,*) 'Boson symm. dimension check passed.'
*
      else if (ISPINFREE.eq.0) then
*
* Set IBOSYM array so that U->alpha and B->beta spinors
*              To get correct Fermion symmetry labeling in
*              LUCIAREL we need to specify same(!) gerade boson
*              irrep for all gerade orbitals and same ungerade
*              boson irrep for all ungerade orbitals (note: which
*              gerade or ungerade boson irrep we choose doesn't
*              matter).
        IOFFFS = 0
        IRROFF = 0
        do IFSYM = 1,NIRR_PN,1
          do ISP = 1+IOFFFS,2*SYMDIM(IFSYM)+IOFFFS,1
            IBOSYM(ISP) = 0 + IRROFF
          end do
*
          if (NTEST.ge.3) then
            write(6,*) ' Reset IBOSYM array, ferm symm. :',IFSYM
            write(6,*) (IBOSYM(I),I=1+IOFFFS,2*SYMDIM(IFSYM)+IOFFFS,1)
          end if
*
          IOFFFS = 2*SYMDIM(IFSYM)
          IRROFF = 2
          if (IDGTYP.eq.3) IRROFF = 1
        end do
      else
        write(6,*) 'ISPINFREE = ',ISPINFREE
        write(6,*) ' has illegal value.'
        stop 'Quitting.'
      end if
*
*  Create symmetry reordering array for this list
*  ... and the inverted array
      ICT = 1
      JCT = 1
      IOFFBS = 0
      do INSYM=1,NIRR_PN,1
        if (INSYM.le.2) then
          IDIMGRP = 1
        else
          IDIMGRP = 2
        end if
        IDIV = 1
        if (IDGTYP.eq.3) then
          IDIMGRP = INSYM
          IDIV = 2
        end if
*
* ... continue with offsets and skipping of even spinor
*     symmetries (these have redundant, i.e. Kramers-paired
*     classification
* For the moment, we hardwire the irrep order of C2h to the
*  generators XY and XYZ (necessarily in this order !!).
*  I.e., the irreps are classified
*  as Ag, Bg, Au, Bu. This should be made more flexible later on.
        ILEN = SYMDIM(IDIMGRP) / IDIV
        if (INSYM.ge.3) IOFFBS = 2*SYMDIM(1)
        if (DOUGRP.eq.4) then
*  INSYM_SEA corresponds to the irrep ordering generated by
*  the choice of symmetry generators in the integral program.
          if (INSYM.eq.1) INSYM_SEA = 1
          if (INSYM.eq.2) INSYM_SEA = 2
          if (INSYM.eq.3) INSYM_SEA = 3
          if (INSYM.eq.4) INSYM_SEA = 4
        else
          INSYM_SEA = INSYM
        end if
        if (INSYM.eq.3) ICT = 1
        if (NTEST.ge.10) then
          write(6,*) 'Current run       : ',INSYM
          write(6,*) 'Seeking symmetry  : ',INSYM_SEA
          write(6,*) 'Doing IND from ',1+IOFFBS, 'to',IOFFBS+ILEN
        end if
        do IND=1+IOFFBS,IOFFBS+ILEN,1
          if ((IBOSYM(IND)+1).eq.INSYM_SEA) then
            ISREO(IND) = ICT
            IND_REL = IND
            if (INSYM.ge.3) IND_REL = IND - SYMDIM(1)
            INVREO(JCT) = IND_REL
            INVREO(-JCT) = -IND_REL
            JCT = JCT + 1
            ICT = ICT + 1
          end if
        end do
      end do
*  ISREO(orbital number DIRAC) = boson symmetry ordered orbital number
*        for one-electron part
*  INVREO(boson symmetry ordered orbital number) = orbital number DIRAC
*        for two-electron part
*
*
      if (NTEST.ge.3) then
        write(6,*) 'Symmetry reordering array DIRAC:'
        do IGRP=1,NIRR_PN/2,1
          if (IGRP.eq.1) then
            write(6,*) 'First dimension group:'
            call iwrtma(ISREO(1),1,SYMDIM(IGRP),1,SYMDIM(IGRP))
          else if (IGRP.eq.2) then
            IFEROFF = 2*SYMDIM(1)
            write(6,*) 'Second dimension group:'
            call iwrtma(ISREO(1+IFEROFF),1,SYMDIM(IGRP),1,SYMDIM(IGRP))
            write(6,*) 'Ferm. index offset for second IDIMGRP: ',IFEROFF
          end if
        end do
        write(6,*) 'Inverse reordering array (One Kramers branch):'
        call iwrtma(INVREO(1),1,ITOTDIM,1,ITOTDIM)
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Reorder Two-electron integrals from MOLECULE symmetry blocks
*  ordering to LUCIAREL symmetry block ordering
*  (A_u and B_g have to be interchanged)
*
      subroutine reotwin_c2h(SWEINT,HELPVEC,MOLOFF,LUCOFF,
     &                       MOLNIT,LUCNIT,IGLSUM_M,IGLSUM_L,
     &                       NSMBLK,NSMBLK_L,ILOCLEN)
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
      dimension SWEINT(IGLSUM_M),HELPVEC(IGLSUM_M)
      dimension MOLOFF(NSMBLK),LUCOFF(NSMBLK_L)
      dimension MOLNIT(NSMBLK),LUCNIT(NSMBLK_L)
*
      NTESTL = 500
!radovan: iprint is not set
      NTEST = max(NTESTL,IPRINT)
*
      if (IGLSUM_M.ne.IGLSUM_L) then
         write(6,*) 'Problem in rdmole2'
         write(6,*) 'IGLSUM_M =  ',IGLSUM_M
         write(6,*) 'IGLSUM_L =  ',IGLSUM_L
         Call Abend2( 'Quitting.' )
      end if
*
      if (2*NSMBLK.gt.ILOCLEN) then
        write(6,*) 'Problem in rdmole2'
        write(6,*) '2*NSMBLK =  ',2*NSMBLK,'larger than '
        write(6,*) 'ILOCLEN =  ',ILOCLEN
        Call Abend2( 'Quitting.' )
      end if
*
      do ISMBLK=1,NSMBLK,1
         if (NTEST.ge.100) write(6,*) 'Doing ISMBLK =  ',ISMBLK
         INOFF = MOLOFF(ISMBLK)
         if (ISMBLK.eq.1) IOUTOFF = LUCOFF(1)
         if (ISMBLK.eq.2) IOUTOFF = LUCOFF(5)
         if (ISMBLK.eq.3) IOUTOFF = LUCOFF(7)
         if (ISMBLK.eq.4) IOUTOFF = LUCOFF(9)
         if (ISMBLK.eq.5) IOUTOFF = LUCOFF(11)
         if (ISMBLK.eq.6) IOUTOFF = LUCOFF(15)
         if (ISMBLK.eq.7) IOUTOFF = LUCOFF(16)
         if (ISMBLK.eq.8) IOUTOFF = LUCOFF(18)
         if (ISMBLK.eq.9) IOUTOFF = LUCOFF(19)
         if (ISMBLK.eq.10) IOUTOFF = LUCOFF(14)
         if (ISMBLK.eq.11) IOUTOFF = LUCOFF(2)
         if (ISMBLK.eq.12) IOUTOFF = LUCOFF(10)
         if (ISMBLK.eq.13) IOUTOFF = LUCOFF(6)
         if (ISMBLK.eq.14) IOUTOFF = LUCOFF(12)
         if (ISMBLK.eq.15) IOUTOFF = LUCOFF(13)
         if (ISMBLK.eq.16) IOUTOFF = LUCOFF(3)
         if (ISMBLK.eq.17) IOUTOFF = LUCOFF(8)
         if (ISMBLK.eq.18) IOUTOFF = LUCOFF(17)
         if (ISMBLK.eq.19) IOUTOFF = LUCOFF(4)
         if (NTEST.ge.100) write(6,*) 'INOFF =  ',INOFF
         if (NTEST.ge.100) write(6,*) 'IOUTOFF =  ',IOUTOFF
         do I=1,MOLNIT(ISMBLK),1
            HELPVEC(IOUTOFF+I-1) = SWEINT(INOFF+I-1)
         end do
         if (NTEST.ge.100) then
            write(6,*) 'transferred integrals:'
            do IWR=1,MOLNIT(ISMBLK),1
               write(6,*) SWEINT(INOFF+IWR-1)
            end do
         end if
      end do
*
      do ICOP=1,IGLSUM_M,1
         SWEINT(ICOP) = HELPVEC(ICOP)
      end do
*
      if (NTEST.ge.100) then
         write(6,*) 'Reordered integrals on SWEINT:'
         do I=1,IGLSUM_M,1
            write(6,*) I,SWEINT(I)
         end do
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine reomoluc(SYMINT_M,ISM,JSM,KSM,LSM,
     &                    IOFF,JOFF,KOFF,LOFF,
     &                    NIORB,NJORB,NKORB,NLORB,INLEN,
     &                    SYMINT_L)
*
*
* Transfer integrals of given symmetry block from Molecule to
* internal LUCIA order
*
*
*. Input :
* SYMINT_M : Input block of MOLECULE_SWEDEN INTEGRAL
*. Output
* SYMINT_L : Complete list of LUCIA integrals
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
*. Input
      dimension SYMINT_M(*)
*. Output
      DIMENSION SYMINT_L(*)
*
      NTEST = 0000
*
      if (NTEST.ge.1000) then
         write(6,*) 'This is reomoluc speaking'
         write(6,*) 'ISM = ',ISM
         write(6,*) '     NIORB = ',NIORB
         write(6,*) 'JSM = ',JSM
         write(6,*) '     NJORB = ',NJORB
         write(6,*) 'KSM = ',KSM
         write(6,*) '     NKORB = ',NKORB
         write(6,*) 'LSM = ',LSM
         write(6,*) '     NLORB = ',NLORB
         write(6,*) 'Input integrals:'
         write(6,*) 'Length of current block: ',INLEN
         do I=1,INLEN,1
            write(6,*) SYMINT_M(I)
         end do
      end if
*
*. Index restrictions
*
      IF(ISM.EQ.JSM) THEN
        IJ_RES = 1
      ELSE
        IJ_RES = 0
      END IF
*
      IF(KSM.EQ.LSM) THEN
       KL_RES = 1
      ELSE
       KL_RES = 0
      END IF
*
      IF(ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
       IJKL_RES = 1
      ELSE
       IJKL_RES = 0
      END IF
*
* K
*  L
*   I
*    J
*
      IJKL_M = 0
      DO K = KOFF, KOFF + NKORB - 1
        IF(KL_RES.EQ.1) THEN
          L_LAST = K
        ELSE
          L_LAST = LOFF + NLORB - 1
        END IF
        DO L = LOFF, L_LAST
          IF(IJKL_RES.EQ.1) THEN
            I_FIRST = K
          ELSE
            I_FIRST = IOFF
          END IF
          DO I = I_FIRST, IOFF + NIORB - 1
            IF(IJ_RES.EQ.1) THEN
              J_LAST = I
            ELSE
              J_LAST = JOFF + NJORB - 1
            END IF
            IF(IJKL_RES.EQ.1.AND.I.EQ.K) THEN
              J_FIRST = L
            ELSE
              J_FIRST = JOFF
            END IF
            DO J = J_FIRST, J_LAST
*. Well, go back to relative ordering
              IR = I - IOFF + 1
              JR = J - JOFF + 1
              KR = K - KOFF + 1
              LR = L - LOFF + 1
              if (NTEST.ge.1000) then
                write(6,'(A,4I4)') 'IR JR KR LR  ',IR,JR,KR,LR
                write(6,'(A,4I4)') 'I J K L  ',I,J,K,L
                write(6,'(A,4I4)') 'IOFF JOFF KOFF LOFF  ',
     &                              IOFF,JOFF,KOFF,LOFF
              end if
C             IJKL_L = I2EAD2(ISM,IR,JSM,JR,KSM,KR,LSM,LR)
              IJKL_M = IJKL_M + 1
              SYMINT_L(IJKL_L) = SYMINT_M(IJKL_M)
              if (NTEST.ge.1000) then
                 WRITE(6,'(A,6I4,E22.12)')
     &           ' I J K L IJKL_M IJKL_L VAL', I,J,K,L,
     &           IJKL_M,IJKL_L, SYMINT_M(IJKL_M)
              end if
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.1000) then
         write(6,*) ' First IJKL_M integrals of SYMINT '
         CALL WRTMAT(SYMINT_L,IJKL_M-1,1,IJKL_M-1,1)
      end if
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      subroutine sspir(ISOI)
*
      implicit real*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "symm.inc"
#include "coeff.inc"
#include "crun.inc"
*
*   Coefficient common block especially for Jeppe's routines
*
*     integer NBSPFSM,NUSPFSM,IBSPEXP
*     real*8 CSPEXP
*     dimension CSPEXP(MXNSPIPDMGRP**2)
*     dimension NBSPFSM(MXNSPIR),NUSPFSM(MXNSPIR)
*     dimension IBSPEXP(2,2)
*     common /COEFF/CSPEXP,NBSPFSM,NUSPFSM,IBSPEXP
*
*
* Number and symmetries of orbitals corresponding to a given shell
*
* =====
* Input
* =====
*
*   DOUGRP  : type of double group
*   NIRREP : Number of irreducible representations per double group
*   MXPIRR : Largest allowed number of shell irreps
*   MXPSBS : Largest allowed number of spinor symmetries
*   MXNSPIR : Maximum number of spatial irreps
*   MXNSPIPDMGRP: Maximum number of spinors per dimension group
*
* ======
* Output
* ======
*
*
* Jeppe Olsen , Winter of 1991
*
* modified: Timo Fleig, Fall of 1996
*
      integer ISSPIR,NSSPIR
*. Output
*
*  Scratch
      do IRREP=1,NIRR_DG,1
         NBSPFSM(IRREP) = 0
         NUSPFSM(IRREP) = 0
      end do
*
      if (DOUGRP.eq.1) then
*=======
* D2h(*)
*=======
        NSMOB = 0
      else if (DOUGRP.eq.2.or.DOUGRP.eq.3) then
         write(*,*) 'These cases have not yet been implemented.'
         Call Abend1( 5 )
*=======
* C2v(*)
*=======
      else if (DOUGRP.eq.4) then
*=======
* C2h(*)
*=======
         SYMDIM(1) = ITOTDM(1) + ITOTDM(2)
         SYMDIM(2) = ITOTDM(3) + ITOTDM(4)
         SYMDIM_F(1) = NORB_F(1) + NORB_F(2)
         SYMDIM_F(2) = NORB_F(3) + NORB_F(4)
         DDIM = (SYMDIM(1))*2 + (SYMDIM(2))*2
         NBSPFSM(5) = ITOTDM(2)
         NBSPFSM(6) = ITOTDM(1)
         NBSPFSM(7) = ITOTDM(4)
         NBSPFSM(8) = ITOTDM(3)
         NUSPFSM(5) = ITOTDM(1)
         NUSPFSM(6) = ITOTDM(2)
         NUSPFSM(7) = ITOTDM(3)
         NUSPFSM(8) = ITOTDM(4)
*
         do I=1,2,1
            ITRIADB(I) = (SYMDIM(I)**2 + SYMDIM(I)) / 2
            ITRIADB_F(I) = (SYMDIM_F(I)**2 + SYMDIM_F(I)) / 2
         end do
         do I=1,2,1
            ITRIADBSUM_F = ITRIADBSUM_F + ITRIADB_F(I)
         end do
      else if (DOUGRP.eq.5.or.DOUGRP.eq.6) then
*=======
* C2(*) and Cs(*)
*=======
         SYMDIM(1) = ITOTDM(1) + ITOTDM(2)
         SYMDIM_F(1) = NORB_F(1) + NORB_F(2)
         DDIM = (SYMDIM(1))*2
         NBSPFSM(3) = ITOTDM(2)
         NBSPFSM(4) = ITOTDM(1)
         NUSPFSM(3) = ITOTDM(1)
         NUSPFSM(4) = ITOTDM(2)
         ITRIADB(1) = (SYMDIM(1)**2 + SYMDIM(1)) / 2
         ITRIADB_F(1) = (SYMDIM_F(1)**2 + SYMDIM_F(1)) / 2
         ITRIADBSUM_F = ITRIADB_F(1)
      else if (DOUGRP.eq.7) then
*=======
* Ci(*)
*=======
         SYMDIM(1) = ITOTDM(1)*2
         SYMDIM(2) = ITOTDM(2)*2
         SYMDIM_F(1) = NORB_F(1)*2
         SYMDIM_F(2) = NORB_F(2)*2
         DDIM = SYMDIM(1) + SYMDIM(2)
         NBSPFSM(3) = ITOTDM(1)
         NBSPFSM(4) = ITOTDM(2)
         NUSPFSM(3) = ITOTDM(1)
         NUSPFSM(4) = ITOTDM(2)
         ITRIADB(1) = (SYMDIM(1)**2 + SYMDIM(1)) / 2
         ITRIADBSUM = ITRIADB(1)
         ITRIADB_F(1) = (SYMDIM_F(1)**2 + SYMDIM_F(1)) / 2
         ITRIADBSUM_F = ITRIADB_F(1)
      else if (DOUGRP.eq.8) then
*=======
* C1(*)
*=======
         SYMDIM(1) = ITOTDM(1)*2
         SYMDIM_F(1) = NORB_F(1)*2
         NBSPFSM(2) = ITOTDM(1)
         NUSPFSM(2) = ITOTDM(1)
         ITRIADB(1) = (SYMDIM(1)**2 + SYMDIM(1)) / 2
         ITRIADBSUM = ITRIADB(1)
         ITRIADB_F(1) = (SYMDIM_F(1)**2 + SYMDIM_F(1)) / 2
         ITRIADBSUM_F = ITRIADB_F(1)
      else
         write(*,*) ' Sorry  DOUGRP out of range , DOUGRP = ', DOUGRP
         write(*,*) ' SSPIR fatally wounded '
         Call Abend1( 5 )
      end if
*
*  Determine spinor info and symmetry within the double group
*
Cjan      call spiinf(LU1)
*
Cjan      write(*,*) 'This is sspir speaking:'
Cjan      do I=1,DDIM,1
Cjan         write(*,*) 'Spinor ',I,' type is   ',SPINOR(I)
Cjan         write(*,*) 'Spinor symmetry is        ',SPISYM(I)
Cjan      end do
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*  1-particle: molfdir full spinor format
*  2-particle: molfdir full spinor format
*
*--------------------------------------------------------------------*
*   November 12, 2001, Timo Fleig
*                - updated for correct one-particle function ordering
*                  Timo Fleig, June 2004
*                - Revised FSM/IUB ordering section             
*                  Timo Fleig, November 2006
*                - revised for proper NO-occupation number runs
*                  Stefan Knecht, October 2007
*                - CINC_NCT0 controls NCT = 0 option
*                  one-particle density run may fail otherwise
*                  Stefan Knecht, July 2008
*--------------------------------------------------------------------*
*
      subroutine trnsfdens(F1,F2,DENS,INDEX,LABEL,IHINDX,SIGN,IOP_REO,
     &                     SIGN_OPREO,ISF,NASHTL,IDEN_TYPE,IPRNT)
*
*  Input:  DENS
*  Output: F1,F2
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
#include "dgroup.h"
#include "dcborb.h"
#include "parluci.h"
#include "noccn_inf.inc"
*
      dimension F1(*),F2(*)
      dimension DENS(NALLINT),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,NALLINT),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),IOFFKF(2,MXPNGAS,2),
     &          IREOKF(MXPORB)
      integer*8 NINT_TOT, NINT_TP_I8, IRIOFF
      integer   NINT_TP
      LOGICAL   CINC_NCT0
*
      NTESTL = 00
*
      CINC_NCT0 = .TRUE.
      IF ( IDEN_TYPE .eq. 1 ) THEN
         CINC_NCT0 = .FALSE.
         NTESTL = 00
      END IF
*      
      NTEST = max(IPRNT,NTESTL)
*
      NINT_TOT = N1ELINT + N2ELINT
      XNINT_TP_I8 = 0.0D0
*
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHTL2 = NASHTL * 2
*
      if (NTEST.ge.2) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================================='
        write(LUWRT,*) ' Transferring densities to MOLFDIR format.'
        write(LUWRT,*) '=========================================='
      end if
*
*  Fermion symmetry/Kramers reordering (for output density indexing)
*  KRMC wants - IFSYM loop
*              - unbarred/barred loop
*               - GAS loop
* Determine offset array
      IADDKF = 1
      do IFSM = 1,NFSYM
        do IUBT = 1,2
          do KGAS = 1,NGAS_DC
            IOFFKF(IFSM,KGAS,IUBT) = IADDKF
            IADDKF = IADDKF + NGSH(IFSM,KGAS)
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        do IUBT = 1,2
          if (IUBT.eq.1) write(LUWRT,*) ' IOFFKF for unbarred.'
          if (IUBT.eq.2) write(LUWRT,*) ' IOFFKF for   barred.'
          call iwrtmamn(IOFFKF(1,1,IUBT),NFSYM,NGAS_DC,2,MXPNGAS,LUWRT)
        end do
      end if
*
* Determine reordering array
      II = 1
      do IUBT = 1,2
        do IFSM = 1,NFSYM
          do KGAS = 1,NGAS_DC
            do IRN = 1,NGSH(IFSM,KGAS)
              IREOKF(II) = IOFFKF(IFSM,KGAS,IUBT) + IRN - 1
              II = II + 1
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        write(LUWRT,*) ' IREOKF array'
        call iwrtmamn(IREOKF,1,II-1,1,MXPORB,LUWRT)
      end if
*
* Process densities
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
*
        if (NTEST.ge.4) then
          write(6,*)
          if (IRI.eq.1) write(LUWRT,*) ' Real part     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary part: '
          write(6,*) ' -------------------'
        end if
*
        do IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(LUWRT,*)
            write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
            write(LUWRT,*)        ' -------------------------'
            write(LUWRT,*)
            write(LUWRT,*) ' LABEL array (for integral indices!):  '
            write(LUWRT,*) ' Internally, LABEL is reord. for operators.'
            call iwrtmamn(LABEL(1,IXTP),1,4,1,4,LUWRT)
            write(LUWRT,*) 'Input LUCIAREL operator order, output dens.'
            call iwrtmamn(IOP_REO(1,IXTP),1,4,1,4,LUWRT)
          end if
*
          NINT_TP_I8 = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8 = REAL(NINT_TP_I8)
          NINT_TP_X = IGIVE_I_B(XNINT_TP_I8)
          NINT_TP = NINT_TP_X
csk       write(LUWRT,*) ' NINT_TP is for IRUN',NINT_TP,IRUN
          do INT = 1,NINT_TP,1
*
*           determine type of density element
*
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
csk         write(luwrt,*) 'NCT is',NCT,'for IRUN',IRUN
*
            if (NCT.eq.2) then
              I = IREOKF(INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHTL)
              J = IREOKF(INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHTL)
*
              IOUT =  (IRI-1)*NASHTL2*NASHTL2 +
     &                 (J-1) *NASHTL2 +
     &                  I
              F1(IOUT) = DENS(IRUN+IRIOFF)
              if (NTEST.ge.10) then
                write(LUWRT,'(2I3,1F18.14,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
              end if
*
            else if (NCT.eq.4) then
              if (SIGN_OPREO(IXTP).eq.ONE.or.
     &            SIGN_OPREO(IXTP).eq.ONEM) then
*
*
*  INDEX array refers to integral indexing.
*  The corresponding second quantization operators then result
*  from the index reordering as given in GET_HX_RELA. 
*  Retaining integral indexing for densities as they come out of
*  the program, and labelled by the INDEX array; this has been
*  generated correctly within DIST_INTS_SIGDEN for all excitation 
*  types. An overall sign for all 2-particle densities is accounted
*  for in GET_HX_RELA (here called SIGN_OPREO):
                IND_R(1) = IREOKF(INDEX_L(1)-(LABEL(1,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator I '
                  write(LUWRT,*) ' INDEX_L(1) ',INDEX_L(1)
                  write(LUWRT,*) ' LABEL(1,IXTP) ',LABEL(1,IXTP)
                end if
*
                IND_R(2) = IREOKF(INDEX_L(2)-(LABEL(2,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator L '
                  write(LUWRT,*) ' INDEX_L(2) ',INDEX_L(2)
                  write(LUWRT,*) ' LABEL(2,IXTP) ',LABEL(2,IXTP)
                end if
*
                IND_R(3) = IREOKF(INDEX_L(3)-(LABEL(3,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator K '
                  write(LUWRT,*) ' INDEX_L(3) ',INDEX_L(3)
                  write(LUWRT,*) ' LABEL(3,IXTP) ',LABEL(3,IXTP)
                end if
*
                IND_R(4) = IREOKF(INDEX_L(4)-(LABEL(4,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator J '
                  write(LUWRT,*) ' INDEX_L(4) ',INDEX_L(4)
                  write(LUWRT,*) ' LABEL(4,IXTP) ',LABEL(4,IXTP)
                end if
*
                IOUT = (IRI-1)*NASHTL2*NASHTL2*NASHTL2*NASHTL2 +
     &                  (IND_R(4)-1) *NASHTL2*NASHTL2*NASHTL2 +
     &                  (IND_R(3)-1) *NASHTL2*NASHTL2 +
     &                  (IND_R(2)-1) *NASHTL2 +
     &                   IND_R(1)
                F2(IOUT) = DENS(IRUN+IRIOFF) * SIGN(IXTP) 
     &                                       * SIGN_OPREO(IXTP)
                if (NTEST.ge.10) then
                  write(LUWRT,*)
                  write(LUWRT,'(4I3,1F18.14,A,I4)') 
     &                      IND_R(1),IND_R(2),IND_R(3),IND_R(4),
     &                      F2(IOUT),' to pos.',IOUT
                end if
*
* ... and copy to the entry for allowed index permutation
* (ab|cd) = -(ad|cb)  in terms of integral indexing
*
                I_R = IND_R(1)
                L_R = IND_R(2)
                K_R = IND_R(3)
                J_R = IND_R(4)
                SIGN_CP = -1 * ONE
                IOUT_CP = (IRI-1) *NASHTL2*NASHTL2*NASHTL2*NASHTL2 +
     &                    (L_R-1) *NASHTL2*NASHTL2*NASHTL2 +
     &                    (K_R-1) *NASHTL2*NASHTL2 +
     &                    (J_R-1) *NASHTL2 +
     &                     I_R
                F2(IOUT_CP) = F2(IOUT) * SIGN_CP
                if (NTEST.ge.10) then
                  write(LUWRT,'(4I3,1F18.14,A,I4)') 
     &                      I_R,J_R,K_R,L_R,F2(IOUT_CP),
     &                      ' to pos.',IOUT_CP
                end if
*
              end if
*             ^ Non-redundant density type (see get_hx_rela)
            else if( (NCT.eq. 0) .and. .NOT. CINC_NCT0 )then
csk               write(LUWRT,*) ' what to do with 4 zero inidces???'
csk               write(LUWRT,*) ' continuing ...'
            else 
              write(LUWRT,*) 'No QFT implemented YET in LUCIAREL.'
              write(LUWRT,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting in trnsfdens.')
            end if
            IRUN = IRUN + 1
          end do
*         ^ Loop over integrals in excitation type
        end do
*       ^ Loop over excitation types
      end do
*     ^ Loop over Real/Imaginary
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*
*  1-particle: molfdir full spinor format
*
*--------------------------------------------------------------------*
*   November 12, 2001, Timo Fleig
*                - updated for correct one-particle function ordering
*                  Timo Fleig, June 2004
*                - Revised FSM/IUB ordering section             
*                  Timo Fleig, November 2006
*                - revised for proper NO-occupation number runs
*                  Stefan Knecht, October 2007
*      
*                - only 1-particle ordering, 
*                  Stefan Knecht - October 2007
*
*--------------------------------------------------------------------*
*
      subroutine trnsfdens_no(F1,DENS,INDEX,LABEL,IHINDX,SIGN,IOP_REO,
     &                     SIGN_OPREO,ISF,NASHTL,IPRNT)
*
*  Input:  DENS
*  Output: F1
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
#include "dgroup.h"
#include "dcborb.h"
#include "parluci.h"
#include "noccn_inf.inc"
*
      dimension F1(*)
      dimension DENS(NALLINT),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,N1ELINT),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),IOFFKF(2,MXPNGAS,2),
     &          IREOKF(MXPORB)
      integer*8 NINT_TOT, NINT_TP_I8, IRIOFF
*
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
      NINT_TOT   = N1ELINT + N2ELINT
*
C
      CALL DZERO(F1,MIN(NZ,2)*(2*NASHT)**2 )
CSK   WRITE(LUWRT,*) ' F1 before transfer'
CSK   CALL WRTMATMN(F1,1,MIN(NZ,2)*(2*NASHT)**2,1,
CSK  &                   MIN(NZ,2)*(2*NASHT)**2,LUWRT)
C
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHTL2 = NASHTL * 2
*
      if (NTEST.ge.2) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================================='
        write(LUWRT,*) ' Transferring densities to MOLFDIR format.'
        write(LUWRT,*) '=========================================='
      end if
*
*  Fermion symmetry/Kramers reordering (for output density indexing)
*  KRMC wants - IFSYM loop
*              - unbarred/barred loop
*               - GAS loop
* Determine offset array
      IADDKF = 1
      do IFSM = 1,NFSYM
        do IUBT = 1,2
          do KGAS = 1,NGAS_DC
            IOFFKF(IFSM,KGAS,IUBT) = IADDKF
            IADDKF = IADDKF + NGSH(IFSM,KGAS)
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        do IUBT = 1,2
          if (IUBT.eq.1) write(LUWRT,*) ' IOFFKF for unbarred.'
          if (IUBT.eq.2) write(LUWRT,*) ' IOFFKF for   barred.'
          call iwrtmamn(IOFFKF(1,1,IUBT),NFSYM,NGAS_DC,2,MXPNGAS,LUWRT)
        end do
      end if
*
* Determine reordering array
      II = 1
      do IUBT = 1,2
        do IFSM = 1,NFSYM
          do KGAS = 1,NGAS_DC
            do IRN = 1,NGSH(IFSM,KGAS)
              IREOKF(II) = IOFFKF(IFSM,KGAS,IUBT) + IRN - 1
              II = II + 1
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        write(LUWRT,*) ' IREOKF array'
        call iwrtmamn(IREOKF,1,II-1,1,MXPORB,LUWRT)
      end if
*
* Process densities
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
*
        if (NTEST.ge.4) then
          write(6,*)
          if (IRI.eq.1) write(LUWRT,*) ' Real part     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary part: '
          write(6,*) ' -------------------'
        end if
*
*
        do IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(LUWRT,*)
            write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
            write(LUWRT,*)        ' -------------------------'
            write(LUWRT,*)
            write(LUWRT,*) ' LABEL array (for integral indices!):  '
            write(LUWRT,*) ' Internally, LABEL is reord. for operators.'
            call iwrtmamn(LABEL(1,IXTP),1,4,1,4,LUWRT)
            write(LUWRT,*) 'Input LUCIAREL operator order, output dens.'
            call iwrtmamn(IOP_REO(1,IXTP),1,4,1,4,LUWRT)
          end if
*
          NINT_TP_I8 = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8 = REAL(NINT_TP_I8)
          NINT_TP_X = IGIVE_I_B(XNINT_TP_I8)
          NINT_TP = NINT_TP_X
csk          write(LUWRT,*) ' NINT_TP is for IRUN and iRI',NINT_TP,IRUN,IRI
*
          IF( IRI. eq. 1 .and. IRUN .eq. (N1ELINT + 1) ) THEN
*
*           real part finished
*
            GOTO 222
*           
          ELSE IF( IRI. eq. 2 .and. IRUN .eq. ( N1ELINT + 1) ) THEN
*
*           well, we have finished the 1-particle part 
*
            GOTO 333
*
          END IF
          do INT = 1,NINT_TP,1
*
* Determine type of density element
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
csk            write(luwrt,*) 'NCT is',NCT,'for IRUN',IRUN
*
            if (NCT.eq.2) then
*
              I = IREOKF(INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHTL)
              J = IREOKF(INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHTL)
*
              IOUT =  (IRI-1)*NASHTL2*NASHTL2 +
     &                 (J-1) *NASHTL2 +
     &                  I
              F1(IOUT) = DENS(IRUN+IRIOFF)
              if (NTEST.ge.10) then
                write(LUWRT,'(2I3,1F18.14,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
                write(LUWRT,*) ' original position IRUN+IRIOFF',
     &                           IRUN+IRIOFF
              end if
*
            else if( NCT .eq. 0 .or. NCT .eq. 4)then
csk               write(LUWRT,*) ' what to do with 4 zero inidces???'
csk               write(LUWRT,*) ' continuing ...'
csk               write(LUWRT,*) ' what to do with 4 nonzero inidces???'
csk               write(LUWRT,*) ' continuing ...'
            else 
              write(LUWRT,*) 'No QFT implemented YET in LUCIAREL.'
              write(LUWRT,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting in trnsfdens_no.')
            end if
            IRUN = IRUN + 1
          end do
*         ^ Loop over integrals in excitation type
        end do
*       ^ Loop over excitation types
 222    CONTINUE
      end do
*     ^ loop over real/imaginary
*
 333  return
      end
#if defined (VAR_MPI2)
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*
*  1-particle: molfdir full spinor format
*
*--------------------------------------------------------------------*
*   November 12, 2001, Timo Fleig
*                - updated for correct one-particle function ordering
*                  Timo Fleig, June 2004
*                - Revised FSM/IUB ordering section             
*                  Timo Fleig, November 2006
*                - revised for proper NO-occupation number runs
*                  Stefan Knecht, October 2007
*      
*                - only 1-particle ordering, 
*                  Stefan Knecht - October 2007
*
*                - shared memory mode: Stefan Knecht - February 2008
*
*--------------------------------------------------------------------*
*
      subroutine trnsfdens_no_sm(F1,DENS,IT_TTPL,IT_TTOL,INDEX,
     &                           LABEL,IHINDX,SIGN,IOP_REO,
     &                           SIGN_OPREO,ISF,NASHTL,IPRNT)
*
*  Input:  DENS
*  Output: F1
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
#include "implicit.h"
#include "mpif.h"
#include "infpar.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
#include "krmc_shmem.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
#include "dgroup.h"
#include "dcborb.h"
#include "noccn_inf.inc"

*
      dimension F1(*), IT_TTPL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IT_TTOL(*)
      dimension DENS(NALLINT),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,N1ELINT),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),IOFFKF(2,MXPNGAS,2),
     &          IREOKF(MXPORB)
      INTEGER*8 NINT_TOT, NINT_TP_I8, IRIOFF, ITB_I8
      INTEGER*8 IRUN_XTP, IOFF_GET 
      INTEGER ID_GET
*
      IRUN_XTP = 0
      IOFF_GET = 0
      ID_GET   = 0
      ITB_I8   = 0
      NTESTL   = 00
      NTEST = max(IPRNT,NTESTL)
*
      NINT_TOT   = N1ELINT + N2ELINT
*
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHTL2 = NASHTL * 2
*
      if (NTEST.ge.2) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================================='
        write(LUWRT,*) ' Transferring densities to MOLFDIR format.'
        write(LUWRT,*) '=========================================='
      end if
*
*  Fermion symmetry/Kramers reordering (for output density indexing)
*  KRMC wants - IFSYM loop
*              - unbarred/barred loop
*               - GAS loop
* Determine offset array
      IADDKF = 1
      do IFSM = 1,NFSYM
        do IUBT = 1,2
          do KGAS = 1,NGAS_DC
            IOFFKF(IFSM,KGAS,IUBT) = IADDKF
            IADDKF = IADDKF + NGSH(IFSM,KGAS)
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        do IUBT = 1,2
          if (IUBT.eq.1) write(LUWRT,*) ' IOFFKF for unbarred.'
          if (IUBT.eq.2) write(LUWRT,*) ' IOFFKF for   barred.'
          call iwrtmamn(IOFFKF(1,1,IUBT),NFSYM,NGAS_DC,2,MXPNGAS,LUWRT)
        end do
      end if
*
* Determine reordering array
      II = 1
      do IUBT = 1,2
        do IFSM = 1,NFSYM
          do KGAS = 1,NGAS_DC
            do IRN = 1,NGSH(IFSM,KGAS)
              IREOKF(II) = IOFFKF(IFSM,KGAS,IUBT) + IRN - 1
              II = II + 1
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        write(LUWRT,*) ' IREOKF array'
        call iwrtmamn(IREOKF,1,II-1,1,MXPORB,LUWRT)
      end if
*
* Process densities
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
*
        if (NTEST.ge.4) then
          write(6,*)
          if (IRI.eq.1) write(LUWRT,*) ' Real part     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary part: '
          write(6,*) ' -------------------'
        end if
*
*
        do IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(LUWRT,*)
            write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
            write(LUWRT,*)        ' -------------------------'
            write(LUWRT,*)
            write(LUWRT,*) ' LABEL array (for integral indices!):  '
            write(LUWRT,*) ' Internally, LABEL is reord. for operators.'
            call iwrtmamn(LABEL(1,IXTP),1,4,1,4,LUWRT)
            write(LUWRT,*) 'Input LUCIAREL operator order, output dens.'
            call iwrtmamn(IOP_REO(1,IXTP),1,4,1,4,LUWRT)
          end if
*
          NINT_TP_I8 = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8 = REAL(NINT_TP_I8)
          NINT_TP_X = IGIVE_I_B(XNINT_TP_I8)
          NINT_TP = NINT_TP_X
csk          write(LUWRT,*) ' NINT_TP is for IRUN and iRI',NINT_TP,IRUN,IRI
*
          IF( IRI. eq. 1 .and. IRUN .eq. (N1ELINT + 1) ) THEN
*
*           real part finished
*
            GOTO 222
*           
          ELSE IF( IRI. eq. 2 .and. IRUN .eq. ( N1ELINT + 1) ) THEN
*
*           well, we have finished the 1-particle part 
*
            GOTO 333
*
          END IF
          IRUN_XTP = 0
          do INT = 1,NINT_TP,1
*
* Determine type of density element
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
csk            write(luwrt,*) 'NCT is',NCT,'for IRUN',IRUN
*
            if (NCT.eq.2) then
*
              I = IREOKF(INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHTL)
              J = IREOKF(INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHTL)
*
              IOUT =  (IRI-1)*NASHTL2*NASHTL2 +
     &                 (J-1) *NASHTL2 +
     &                  I
csk              F1(IOUT) = DENS(IRUN+IRIOFF)
              IXXTTP = IXTP
              IF( IRI .eq. 2 ) THEN
                IXXTTP = ITTP + NSPOBEX_TP
              END IF
              IF( IRUN_XTP .eq. 0 )THEN
C               process id and window offset (complex case included)
                ID_GET = IT_TTPL( IXXTTP )
C               check for availability
                IF( ID_GET .lt. 0) THEN
                  IRUN = IRUN + NINT_TP
                  GOTO 155
                END IF
C               window offset
                ITB_I8 = IT_TTOL( IXXTTP )
              END IF
              IOFF_GET = ITB_I8 + IRUN_XTP
              YDENS = 0.0D0
              IF( MYNEW_ID_SM .eq. ID_GET )THEN
                YDENS = DENS(IOFF_GET+1)
              ELSE
                CALL MPIXGET(YDENS,1,2,ID_GET,
     &                       IOFF_GET,1,2,MY_T_WIN)
              END IF
              F1(IOUT) = YDENS
              if (NTEST.ge.10) then
                write(LUWRT,'(2I3,1F18.14,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
                write(LUWRT,*) ' original position IRUN+IRIOFF',
     &                           IRUN+IRIOFF, IOFF_GET
              end if
*
            else if( NCT .eq. 0 .or. NCT .eq. 4)then
csk               write(LUWRT,*) ' what to do with 4 zero inidces???'
csk               write(LUWRT,*) ' continuing ...'
csk               write(LUWRT,*) ' what to do with 4 nonzero inidces???'
csk               write(LUWRT,*) ' continuing ...'
            else 
              write(LUWRT,*) 'No QFT implemented YET in LUCIAREL.'
              write(LUWRT,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting in trnsfdens_no_sm.')
            end if
            IRUN_XTP = IRUN_XTP + 1
            IRUN = IRUN + 1
          end do
*         ^ Loop over integrals in excitation type
 155      CONTINUE
        end do
*       ^ Loop over excitation types
 222    CONTINUE
      end do
*     ^ Loop over Real/Imaginary
*
 333  CONTINUE
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*
*  1-particle: molfdir full spinor format
*
*--------------------------------------------------------------------*
*   November 12, 2001, Timo Fleig
*                - updated for correct one-particle function ordering
*                  Timo Fleig, June 2004
*                - Revised FSM/IUB ordering section             
*                  Timo Fleig, November 2006
*                - revised for proper NO-occupation number runs
*                  Stefan Knecht, October 2007
*      
*                - only 1-particle ordering, 
*                  Stefan Knecht - October 2007
*
*--------------------------------------------------------------------*
*
      subroutine trnsfdens_no_split(F1,DENS,DENS_SCR,IT_TTPL,INDEX,
     &                              LABEL,IHINDX,SIGN,
     &                              IOP_REO,SIGN_OPREO,ISF,NASHTL,
     &                              IPRNT)
*
*  Input:  DENS
*  Output: F1
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
#include "dgroup.h"
#include "dcborb.h"
#include "parluci.h"
#include "noccn_inf.inc"
#include "krmc_shmem.h"
*
      dimension F1(*), DENS_SCR(*), IT_TTPL(*)
      dimension DENS(*),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,N1ELINT),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),IOFFKF(2,MXPNGAS,2),
     &          IREOKF(MXPORB)
      integer*8 NINT_TOT, ITB_BUFF, ITB_BUFF_ADD
*
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
      NINT_TOT     = N1ELINT + N2ELINT
      ITB_BUFF_ADD = 0
      ITB_BUFF     = 0
*
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHTL2 = NASHTL * 2
*
      if (NTEST.ge.2) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================================='
        write(LUWRT,*) ' Transferring densities to MOLFDIR format.'
        write(LUWRT,*) '=========================================='
      end if
*
*  Fermion symmetry/Kramers reordering (for output density indexing)
*  KRMC wants - IFSYM loop
*              - unbarred/barred loop
*               - GAS loop
* Determine offset array
      IADDKF = 1
      do IFSM = 1,NFSYM
        do IUBT = 1,2
          do KGAS = 1,NGAS_DC
            IOFFKF(IFSM,KGAS,IUBT) = IADDKF
            IADDKF = IADDKF + NGSH(IFSM,KGAS)
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        do IUBT = 1,2
          if (IUBT.eq.1) write(LUWRT,*) ' IOFFKF for unbarred.'
          if (IUBT.eq.2) write(LUWRT,*) ' IOFFKF for   barred.'
          call iwrtmamn(IOFFKF(1,1,IUBT),NFSYM,NGAS_DC,2,MXPNGAS,LUWRT)
        end do
      end if
*
* Determine reordering array
      II = 1
      do IUBT = 1,2
        do IFSM = 1,NFSYM
          do KGAS = 1,NGAS_DC
            do IRN = 1,NGSH(IFSM,KGAS)
              IREOKF(II) = IOFFKF(IFSM,KGAS,IUBT) + IRN - 1
              II = II + 1
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        write(LUWRT,*) ' IREOKF array'
        call iwrtmamn(IREOKF,1,II-1,1,MXPORB,LUWRT)
      end if
*
* Process densities
      do IRI = 1,IRILP,1
        IRUN = 1
        ITB_BUFF = 1
        do JRI = 1,IRI-1,1
          ITB_BUFF = ITB_BUFF + LEN_T_BUFF
        end do
csk     write(LUWRT,*) '  ITB_BUFF, LEN_T_BUFF',ITB_BUFF,LEN_T_BUFF

*
        if (NTEST.ge.4) then
          write(LUWRT,*)
          if (IRI.eq.1) write(LUWRT,*) ' Real part     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary part: '
          write(LUWRT,*) ' -------------------'
        end if
*
*
        do 111 IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(LUWRT,*)
            write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
            write(LUWRT,*)        ' -------------------------'
            write(LUWRT,*)
            write(LUWRT,*) ' LABEL array (for integral indices!):  '
            write(LUWRT,*) ' Internally, LABEL is reord. for operators.'
            call iwrtmamn(LABEL(1,IXTP),1,4,1,4,LUWRT)
            write(LUWRT,*) 'Input LUCIAREL operator order, output dens.'
            call iwrtmamn(IOP_REO(1,IXTP),1,4,1,4,LUWRT)
          end if
*
          IXXTTP = (( IRI - 1 ) * NSPOBEX_TP) + IXTP
C         active T block? (complex case included)
          ID_GET = IT_TTPL( IXXTTP )
csk       WRITE(6,*) ' type IXXTTP and ID_GET',IXXTTP,ID_GET,MYPROC
          ITB_BUFF_ADD = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8  = REAL(ITB_BUFF_ADD)
          NINT_TP      = IGIVE_I_B(XNINT_TP_I8)
          IF( ID_GET .lt. 0 ) THEN
            IRUN         = IRUN + NINT_TP
            ITB_BUFF_ADD = 0
          ELSE
C           receive or copy DENS block - at present: MASTER only
            IF( MYPROC .eq. MASTER )THEN
              CALL DZERO(DENS_SCR,NINT_TP)
              IF( ID_GET .eq. 0 )THEN
                CALL DCOPY(NINT_TP,DENS(ITB_BUFF),1,DENS_SCR,1)
csk             WRITE(LUWRT,*) ' density block at ITB_BUFF'
csk             CALL WRTMATMN(DENS_SCR,1,NINT_TP,1,NINT_TP,LUWRT)
              ELSE
                CALL MPI_RECV(DENS_SCR,NINT_TP,MPI_REAL8,ID_GET,821,
     &                        MPI_COMM_WORLD,ISTAT,IERR)
                ITB_BUFF_ADD = 0
              END IF
            ELSE
              IRUN         = IRUN + NINT_TP
              IF( MYPROC .eq. ID_GET )THEN
                CALL MPI_SEND(DENS(ITB_BUFF),NINT_TP,MPI_REAL8,
     &                        MASTER,821,MPI_COMM_WORLD,IERR)
                ITB_BUFF     = ITB_BUFF + ITB_BUFF_ADD
              ELSE
                ITB_BUFF_ADD = 0
              END IF
            END IF
          END IF
*
*
          IF( IRI. eq. 1 .and. IRUN .eq. (N1ELINT + 1) ) THEN
*
*           real part finished
*
            GOTO 222
*           
          ELSE IF( IRI. eq. 2 .and. IRUN .eq. ( N1ELINT + 1) ) THEN
*
*           well, we have finished the 1-particle part 
*
            GOTO 333
*
          END IF
C
          IF( MYPROC .ne. MASTER ) THEN
              GOTO 111
          ELSE
              IF( ID_GET .lt. 0 ) GOTO 111
          END IF
C
          IOFF_SET = 1
C
          do INT = 1,NINT_TP,1
*
*           Determine type of density element
*
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
            if (NCT.eq.2) then
*
              I = IREOKF(INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHTL)
              J = IREOKF(INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHTL)
*
              IOUT =  (IRI-1)*NASHTL2*NASHTL2 +
     &                 (J-1) *NASHTL2 +
     &                  I
              F1(IOUT) = DENS_SCR(IOFF_SET)
              if (NTEST.ge.10) then
                write(LUWRT,'(2I3,1F18.14,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
                write(LUWRT,*) ' original position IOFF_SET',
     &                           IOFF_SET
              end if
*
            else if( NCT .eq. 0 .or. NCT .eq. 4)then
csk               write(LUWRT,*) ' what to do with 4 zero inidces???'
csk               write(LUWRT,*) ' continuing ...'
csk               write(LUWRT,*) ' what to do with 4 nonzero inidces???'
csk               write(LUWRT,*) ' continuing ...'
            else 
              write(LUWRT,*) 'No QFT implemented YET in LUCIAREL.'
              write(LUWRT,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting in trnsfdens_no.')
            end if
            IRUN     = IRUN     + 1
            IOFF_SET = IOFF_SET + 1
          end do
*         ^ Loop over integrals in excitation type
          ITB_BUFF     = ITB_BUFF + ITB_BUFF_ADD
 111    CONTINUE
*       ^ Loop over excitation types
 222    CONTINUE
      end do
*     ^ Loop over Real/Imaginary
*
 333  CONTINUE
*
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Transfer density matrices from excitation class ordering to
* KRMC storage mode
*  1-particle: molfdir full spinor format
*  2-particle: molfdir full spinor format
*
*--------------------------------------------------------------------*
*   November 12, 2001, Timo Fleig
*                - updated for correct one-particle function ordering
*                  Timo Fleig, June 2004
*                - Revised FSM/IUB ordering section             
*                  Timo Fleig, November 2006
*                - T split version                             
*                  Stefan Knecht, May 2008
*--------------------------------------------------------------------*
*
      subroutine trnsfdens_split(F1,F2,DENS,DENS_SCR,IT_TTPL,
     &                           INDEX,LABEL,IHINDX,SIGN,IOP_REO,
     &                           SIGN_OPREO,ISF,NASHTL,IPRNT)
*
*  Input:  DENS
*  Output: F1,F2
*
* #############################################
*  LABEL is defined as 
*     1 for unbarred (alpha)
*    -1 for barred   (beta)
* #############################################
*
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
#include "dgroup.h"
#include "dcborb.h"
#include "parluci.h"
#include "noccn_inf.inc"
#include "krmc_shmem.h"
*
      dimension F1(*),F2(*), DENS_SCR(*), IT_TTPL(*)
      dimension DENS(*),SIGN(NSPOBEX_TP),SIGN_OPREO(NSPOBEX_TP)
      dimension INDEX(4,NALLINT),ICSAVE(4),INDEX_L(4),IND_R(4)
      dimension LABEL(4,NSPOBEX_TP),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),IOFFKF(2,MXPNGAS,2),
     &          IREOKF(MXPORB)
      integer*8 NINT_TOT, ITB_BUFF, ITB_BUFF_ADD
      integer   NINT_TP
*
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
      NINT_TOT     = N1ELINT + N2ELINT
      ITB_BUFF_ADD = 0
      ITB_BUFF     = 0
*
      if (ISF.eq.1.or.NZ.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      ONE = 1.0D0
      ONEM = -1.0D0
      NASHTL2 = NASHTL * 2
*
      if (NTEST.ge.2) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================================='
        write(LUWRT,*) ' Transferring densities to MOLFDIR format.'
        write(LUWRT,*) '=========================================='
      end if
*
*  Fermion symmetry/Kramers reordering (for output density indexing)
*  KRMC wants - IFSYM loop
*              - unbarred/barred loop
*               - GAS loop
* Determine offset array
      IADDKF = 1
      do IFSM = 1,NFSYM
        do IUBT = 1,2
          do KGAS = 1,NGAS_DC
            IOFFKF(IFSM,KGAS,IUBT) = IADDKF
            IADDKF = IADDKF + NGSH(IFSM,KGAS)
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        do IUBT = 1,2
          if (IUBT.eq.1) write(LUWRT,*) ' IOFFKF for unbarred.'
          if (IUBT.eq.2) write(LUWRT,*) ' IOFFKF for   barred.'
          call iwrtmamn(IOFFKF(1,1,IUBT),NFSYM,NGAS_DC,2,MXPNGAS,LUWRT)
        end do
      end if
*
* Determine reordering array
      II = 1
      do IUBT = 1,2
        do IFSM = 1,NFSYM
          do KGAS = 1,NGAS_DC
            do IRN = 1,NGSH(IFSM,KGAS)
              IREOKF(II) = IOFFKF(IFSM,KGAS,IUBT) + IRN - 1
              II = II + 1
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.5) then
        write(LUWRT,*) ' IREOKF array'
        call iwrtmamn(IREOKF,1,II-1,1,MXPORB,LUWRT)
      end if
*
* Process densities
      do IRI = 1,IRILP,1
        IRUN     = 1
        ITB_BUFF = 1
        do JRI = 1,IRI-1,1
          ITB_BUFF = ITB_BUFF + LEN_T_BUFF
        end do
*
        if (NTEST.ge.4) then
          write(LUWRT,*)
          if (IRI.eq.1) write(LUWRT,*) ' Real part     : '
          if (IRI.eq.2) write(LUWRT,*) ' Imaginary part: '
          write(LUWRT,*) ' -------------------'
        end if
*
        do 111 IXTP = 1,NSPOBEX_TP,1
*
          if (NTEST.ge.4) then
            write(LUWRT,*)
            write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
            write(LUWRT,*)        ' -------------------------'
            write(LUWRT,*)
            write(LUWRT,*) ' LABEL array (for integral indices!):  '
            write(LUWRT,*) ' Internally, LABEL is reord. for operators.'
            call iwrtmamn(LABEL(1,IXTP),1,4,1,4,LUWRT)
            write(LUWRT,*) 'Input LUCIAREL operator order, output dens.'
            call iwrtmamn(IOP_REO(1,IXTP),1,4,1,4,LUWRT)
          end if
*
          IXXTTP = (( IRI - 1 ) * NSPOBEX_TP) + IXTP
C         active T block? (complex case included)
          ID_GET = IT_TTPL( IXXTTP )
csk       WRITE(6,*) ' type IXXTTP and ID_GET',IXXTTP,ID_GET,MYPROC
          ITB_BUFF_ADD = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8  = REAL(ITB_BUFF_ADD)
          NINT_TP      = IGIVE_I_B(XNINT_TP_I8)
          IF( ID_GET .lt. 0 ) THEN
            IRUN         = IRUN + NINT_TP
            ITB_BUFF_ADD = 0
            GOTO 111
          ELSE
C           receive or copy DENS block - at present: MASTER only
            IF( MYPROC .eq. MASTER )THEN
              CALL DZERO(DENS_SCR,NINT_TP)
              IF( ID_GET .eq. 0 )THEN
                CALL DCOPY(NINT_TP,DENS(ITB_BUFF),1,DENS_SCR,1)
              ELSE
                CALL MPI_RECV(DENS_SCR,NINT_TP,MPI_REAL8,ID_GET,821,
     &                        MPI_COMM_WORLD,ISTAT,IERR)
                ITB_BUFF_ADD = 0
              END IF
            ELSE
              IRUN         = IRUN + NINT_TP
              IF( MYPROC .eq. ID_GET )THEN
                CALL MPI_SEND(DENS(ITB_BUFF),NINT_TP,MPI_REAL8,
     &                        MASTER,821,MPI_COMM_WORLD,IERR)
                ITB_BUFF     = ITB_BUFF + ITB_BUFF_ADD
              ELSE
                ITB_BUFF_ADD = 0
              END IF
            END IF
          END IF
C
C         MASTER enters loop over integrals
C
          IF( MYPROC .ne. MASTER ) GOTO 111
C
          IOFF_SET = 1
C
          do INT = 1,NINT_TP,1
*
*           Determine type of density element
*
            NCT = 0
            do IC=1,4,1
              if (INDEX(IC,IRUN).ne.0) then
                NCT = NCT + 1
                ICSAVE(NCT) = IC
                IUBFAC = 1
                if (INDEX(IC,IRUN).lt.0) IUBFAC = -1
                INDEX_L(IC) = INDEX(IC,IRUN) * IUBFAC
              end if
            end do
*
            if (NCT.eq.2) then
              I = IREOKF(INDEX_L(ICSAVE(1))  
     &            - (LABEL(IHINDX(1,IXTP),IXTP) - 1)/2 * NASHTL)
              J = IREOKF(INDEX_L(ICSAVE(2))  
     &            - (LABEL(IHINDX(2,IXTP),IXTP) - 1)/2 * NASHTL)
*
              IOUT =  (IRI-1)*NASHTL2*NASHTL2 +
     &                 (J-1) *NASHTL2 +
     &                  I
              F1(IOUT) = DENS_SCR(IOFF_SET)
              if (NTEST.ge.10) then
                write(LUWRT,'(2I3,1F18.14,A,I4)') 
     &                    I,J,F1(IOUT),' to pos.',IOUT
              end if
*
            else if (NCT.eq.4) then
              if (SIGN_OPREO(IXTP).eq.ONE.or.
     &            SIGN_OPREO(IXTP).eq.ONEM) then
*
*
*  INDEX array refers to integral indexing.
*  The corresponding second quantization operators then result
*  from the index reordering as given in GET_HX_RELA. 
*  Retaining integral indexing for densities as they come out of
*  the program, and labelled by the INDEX array; this has been
*  generated correctly within DIST_INTS_SIGDEN for all excitation 
*  types. An overall sign for all 2-particle densities is accounted
*  for in GET_HX_RELA (here called SIGN_OPREO):
                IND_R(1) = IREOKF(INDEX_L(1)-(LABEL(1,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator I '
                  write(LUWRT,*) ' INDEX_L(1) ',INDEX_L(1)
                  write(LUWRT,*) ' LABEL(1,IXTP) ',LABEL(1,IXTP)
                end if
*
                IND_R(2) = IREOKF(INDEX_L(2)-(LABEL(2,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator L '
                  write(LUWRT,*) ' INDEX_L(2) ',INDEX_L(2)
                  write(LUWRT,*) ' LABEL(2,IXTP) ',LABEL(2,IXTP)
                end if
*
                IND_R(3) = IREOKF(INDEX_L(3)-(LABEL(3,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator K '
                  write(LUWRT,*) ' INDEX_L(3) ',INDEX_L(3)
                  write(LUWRT,*) ' LABEL(3,IXTP) ',LABEL(3,IXTP)
                end if
*
                IND_R(4) = IREOKF(INDEX_L(4)-(LABEL(4,IXTP)-1)/2*NASHTL)
                if (NTEST.ge.50) then
                  write(LUWRT,*)
                  write(LUWRT,*) ' operator J '
                  write(LUWRT,*) ' INDEX_L(4) ',INDEX_L(4)
                  write(LUWRT,*) ' LABEL(4,IXTP) ',LABEL(4,IXTP)
                end if
*
                IOUT = (IRI-1)*NASHTL2*NASHTL2*NASHTL2*NASHTL2 +
     &                  (IND_R(4)-1) *NASHTL2*NASHTL2*NASHTL2 +
     &                  (IND_R(3)-1) *NASHTL2*NASHTL2 +
     &                  (IND_R(2)-1) *NASHTL2 +
     &                   IND_R(1)
                F2(IOUT) = DENS_SCR(IOFF_SET) * SIGN(IXTP) 
     &                                        * SIGN_OPREO(IXTP)
                if (NTEST.ge.10) then
                  write(LUWRT,*)
                  write(LUWRT,'(4I3,1F18.14,A,I4)') 
     &                      IND_R(1),IND_R(2),IND_R(3),IND_R(4),
     &                      F2(IOUT),' to pos.',IOUT
                end if
*
* ... and copy to the entry for allowed index permutation
* (ab|cd) = -(ad|cb)  in terms of integral indexing
*
                I_R = IND_R(1)
                L_R = IND_R(2)
                K_R = IND_R(3)
                J_R = IND_R(4)
                SIGN_CP = -1 * ONE
                IOUT_CP = (IRI-1) *NASHTL2*NASHTL2*NASHTL2*NASHTL2 +
     &                    (L_R-1) *NASHTL2*NASHTL2*NASHTL2 +
     &                    (K_R-1) *NASHTL2*NASHTL2 +
     &                    (J_R-1) *NASHTL2 +
     &                     I_R
                F2(IOUT_CP) = F2(IOUT) * SIGN_CP
                if (NTEST.ge.10) then
                  write(LUWRT,'(4I3,1F18.14,A,I4)') 
     &                      I_R,J_R,K_R,L_R,F2(IOUT_CP),
     &                      ' to pos.',IOUT_CP
                end if
*
              end if
*             ^ Non-redundant density type (see get_hx_rela)
            else 
              write(LUWRT,*) 'No QFT implemented YET in LUCIAREL.'
              write(LUWRT,*) 'you specified ',4-NCT,' indices.'
              call abend2('Quitting in trnsfdens.')
            end if
            IRUN = IRUN + 1
            IOFF_SET = IOFF_SET + 1
          end do
*         ^ Loop over integrals in excitation type
        ITB_BUFF     = ITB_BUFF + ITB_BUFF_ADD
 111    CONTINUE
*       ^ Loop over excitation types
      end do
*     ^ Loop over Real/Imaginary
*
      return
      end
#endif
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Write density matrix for current eigenvector in 
* SIGDEN ordering.
*
*   October 18, 2001, Timo Fleig
*
*   Revised for operator ordering, Mar 20, 2002
*
      subroutine wrtdens(DENS,INDEX,IHINDX,SIGN,NASHT,LABEL,IOP_REO,
     &                   SIGN_REO,ISF)
*
#include "implicit.h"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
#include "parluci.h"
*
      dimension DENS(NALLINT),SIGN(NSPOBEX_TP)
      dimension INDEX(4,NALLINT),IND_R(4),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP),LABEL(4,NSPOBEX_TP),
     &          SIGN_REO(NSPOBEX_TP)
      integer*8 NINT_TOT, NINT_TP_I8, IRIOFF 
*
      NINT_TOT = N1ELINT + N2ELINT
      if (ISF.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
        write(LUWRT,*) '  '
        write(LUWRT,*) '  '
        write(LUWRT,'(2X,A)') 
     &            'Writing densities in reordered form (for MCSCF):'
        write(LUWRT,*)'  - a+ a b b+  '
        write(LUWRT,'(2X,A)') 
     &            '(means sign multiplied to all 2-particle densities)'
        write(LUWRT,*)'  (done in GET_HX_RELA)'
        write(LUWRT,*) '  '
        if (IRI.eq.1) write(LUWRT,*) ' Real density elements     : '
        if (IRI.eq.2) write(LUWRT,*) ' Imaginary density elements: '
        write(LUWRT,*) ' -------------------------------'
        do IXTP = 1,NSPOBEX_TP,1
          write(LUWRT,*) '  '
          write(LUWRT,'(A,I4)') ' Excitation type ... ',IXTP
          write(LUWRT,*)        ' -------------------------'
          NINT_TP_I8 = NELM_CC(IXTP,NSPOBEX_TP,NINT_TOT)
          XNINT_TP_I8 = REAL(NINT_TP_I8)
          NINT_TP_X = IGIVE_I_B(XNINT_TP_I8)
          NINT_TP = NINT_TP_X
csk          write(LUWRT,*) ' NINT_TP is for IRUN',NINT_TP,IRUN
          do INT = 1,NINT_TP,1
*
            if (INDEX(1,IRUN).ne.0.and.INDEX(2,IRUN).ne.0) then
*  Two-particle density element.
*
*  INDEX array refers to integral indexing.
              IUBFAC = 1
              if (INDEX(1,IRUN).lt.0) IUBFAC = -1
              IND_R(1) = INDEX(1,IRUN) - (LABEL(1,IXTP) - 1)/2 * NASHT
*
              IUBFAC = 1
              if (INDEX(2,IRUN).lt.0) IUBFAC = -1
              IND_R(2) = INDEX(2,IRUN) - (LABEL(2,IXTP) - 1)/2 * NASHT
*
              IUBFAC = 1
              if (INDEX(3,IRUN).lt.0) IUBFAC = -1
              IND_R(3) = INDEX(3,IRUN) - (LABEL(3,IXTP) - 1)/2 * NASHT
*
              IUBFAC = 1
              if (INDEX(4,IRUN).lt.0) IUBFAC = -1
              IND_R(4) = INDEX(4,IRUN) - (LABEL(4,IXTP) - 1)/2 * NASHT
            else
*  One-particle density element. No modification.
              write(LUWRT,*) ' One-particle density element for IRUN',
     &                         IRUN
              IND_R(1) = INDEX(1,IRUN)
              IND_R(2) = INDEX(2,IRUN)
              IND_R(3) = INDEX(3,IRUN)
              IND_R(4) = INDEX(4,IRUN)
            end if
            write(LUWRT,'(1X,4I3,2X,1F18.10)') 
     &           (IND_R(I),I=1,4,1),DENS(IRUN+IRIOFF) * SIGN(IXTP)
     &                                                * SIGN_REO(IXTP)
            IRUN = IRUN + 1
          end do
        end do
      end do
c     stop 'AFTER WRTDENS.'
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
      subroutine LCR_RD_MRCMDC(TINTR,TINTI,INDK,INDL,RKLR,RKLI)
*
#include "implicit.h"
#include "mxpdim.inc"
#include "clunit.inc"
*
      character*4 REPNA(14),REPN(14)
      character*7 FNAM1
      character*6 FNAM2,FNAM3
      character*8 TIMEX
      character*10 DATEX,DATEXB,TIMEXB(8)
      character*10 NEXTFIL
      logical BREIT,SPINFR,REALAR
*
      dimension MULTB(MXPIRR*2,MXPIRR*2)
      dimension IRPMO(MXPORB),IRPAMO(MXPORB),EPS(MXPORB)
      dimension TINTR(*),TINTI(*)
      dimension KR(-MXPORB:MXPORB)
      dimension INDK(*),INDL(*)
      dimension RKLR(*),RKLI(*)
*
      FNAM1 = 'MRCONEE'
      FNAM2 = 'MDCINT'
*
* Read one-electron integrals
*
      open (LUSC96,FILE=FNAM1,FORM='UNFORMATTED')
       print*,'File ',FNAM1,' opened.'
      rewind (LUSC96)
      read (LUSC96) NORB,BREIT,ECORE,NFSYM,NZ,SPINFR
       print*,'NORB,BREIT,ECORE,NFSYM: ',NORB,BREIT,ECORE,NFSYM
       print*,' (NORB is the number of spinors.)'
       print*,'NZ,SPINFR: ',NZ,SPINFR
      REALAR = NZ.eq.1.or.SPINFR
       print*,'REALAR ',REALAR
c     IF (NORB.gt.MXPORB) call Q2JOB(3,'P.ORBSIN','MXPORB',NORB)
      read (LUSC96) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
       print*,'NSYMRP: ',NSYMRP
       print*,'REPN: ',(REPN(IRP),IRP=1,NSYMRP)
      read (LUSC96) NREP,(REPNA(IRP),IRP=1,2*NREP)
       print*,'NREP: ',NREP
       print*,'REPNA: ',(REPNA(IRP),IRP=1,2*NREP)
      read (LUSC96) ((MULTB(I,J),I=1,2*NREP),J=1,2*NREP)
       call iwrtma(MULTB(1,1),NREP*2,NREP*2,MXPIRR*2,MXPIRR*2)
      read (LUSC96) (IRPMO(IMO),IRPAMO(IMO),EPS(IMO),IMO=1,NORB)
c      print*,'IRPMO: ',(IRPMO(I),I=1,NORB)
c      print*,'IRPAMO: ',(IRPAMO(I),I=1,NORB)
c      print*,'EPS: ',(EPS(I),I=1,NORB)
      if (REALAR) then
         read (LUSC96) (TINTR(JI),TDUM,JI=1,NORB*NORB)
      else
         read (LUSC96) (TINTR(JI),TINTI(JI),JI=1,NORB*NORB)
      end if
       do IP=1,NORB**2
        print*,'TINTR(',IP,') = ',TINTR(IP)
       end do
      close (LUSC96)
*
* Read two-electron integrals
*
      open (LUSC97,file=FNAM2,form='UNFORMATTED')
       print*,'File ',FNAM2,' opened.'
      rewind (LUSC97)
      read (LUSC97,err=10000,end=10000) DATEX,TIMEX,NKR,
     &                                  (KR(I),KR(-I),I=1,NKR)
       print*,'DATEX,TIMEX,NKR ',DATEX,TIMEX,NKR
       print*,'KR(NKR)..KR(-NKR) ',(KR(I),I=-NKR,NKR)
      IF (2*NKR.NE.NORB)
     &   CALL QUIT('MRCONEE AND MDCINT ARE INCOMPATIBLE')
      IF (BREIT) THEN
        call quit ('BREIT interaction not yet implemented.')
c       open (MDBINT,file=FNAM(3),form='UNFORMATTED')
c       read (mdbint,err=10001,end=10001) datexb,timexb,nkrb
c       IF (NKR.NE.NKRB) CALL QUIT('MDCINT AND MDBINT ARE INCOMPATIBLE')
      END IF
      IF (REALAR) THEN
c        DO I = 1, M1
c           RKLI(I) = 0.0
c        ENDDO
      ENDIF

c     IF (REALAR) THEN
c        read (LUSC97,ERR=10010,END=10010) ikr,jkr,nz,
c    &                (indk(inz),indl(inz),inz=1,nz),
c    &                (rklr(inz),inz=1,nz)
c     ELSE
c        read (LUSC97,ERR=10010,END=10010) ikr,jkr,nz,
c    &                (indk(inz),indl(inz),inz=1,nz),
c    &                (rklr(inz),rkli(inz),inz=1,nz)
c     END IF
       print*,'IKR,JKR,NZ ',IKR,JKR,NZ
       print*,'INDK : ',(INDK(I),I=1,NZ)
       print*,'INDL : ',(INDL(I),I=1,NZ)
       print*,'RKLR : ',(RKLR(I),I=1,NZ)



      close (LUSC97)
       print*,'File ',FNAM2,' closed.'


      return
*
10000 CALL QUIT('ERROR READING HEADER OF MDCINT')
10001 CALL QUIT('ERROR READING HEADER OF MDBINT')
10010 CALL QUIT('ERROR READING INTEGRALS FROM MDCINT')
10011 CALL QUIT('ERROR READING INTEGRALS FROM MDBINT')
*
      end
