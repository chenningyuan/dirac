*  Distribute external Fock type one-electron integrals
*  to double group lists
*  Pick active spin-orbit integrals
*
      subroutine dist_1eints(HELPMAT,ONEINT,SPININTS,ITRALINE,
     &                       IMXDIM,IPRINT)
*
      implicit real*8 (A-H,O-Z)
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
        call memman(IDUM,IDUM,'MARK  ',IDUM,'DIST1E')
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
        call memman(KPICTRA,ITRALEN,'ADDL  ',2,'IPITRA')
        call setvec(WORK(KPICTRA),0.D0,ITRALEN)
        call memman(KIGSIND,NIRR_PN*SYMDIM(IDIMGRP),
     &              'ADDL  ',2,'IGSIND')
        call memman(KISMIND,SYMDIM(IDIMGRP),'ADDL  ',2,'ISMIND')
        call memman(KIGSOFF,2*NGASPI,'ADDL  ',2,'IGSOFF')
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
        call memman(IDUM,IDUM,'FLUSM ',IDUM,'DIST1E')
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
