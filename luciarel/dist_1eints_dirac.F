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
      call memman(IDUM,IDUM,'MARK  ',IDUM,'DISTDI')
      call memman(KIGSIND,NIRR_PN*SYMDIM(IDIMGRP),
     &              'ADDL  ',2,'IGSIND')
      call memman(KISMIND,SYMDIM(IDIMGRP),'ADDL  ',2,'ISMIND')
      call memman(KIGSOFF,2*NGASPI,'ADDL  ',2,'IGSOFF')
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
      call memman(IDUM,IDUM,'FLUSM ',IDUM,'DISTDI') 
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
