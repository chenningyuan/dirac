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
        if (IRC.ne.0) goto 991
*
        IOPT = 0
        call rdone(IRC,IOPT,LABEL,ICOMP,HELPMAT,ISYLBL)
        if (IRC.ne.0) goto 991
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
991   Call StkA(0,OneLbl)
      Call ErrMsg('RsoInt',002)
      Call Quit(16)
992   Call ErrMsg('RsoInt',003)
      Call Quit(16)
993   Call ErrMsg('RsoInt',004)
      Call Quit(16)
994   Call StkA(0,ComLbl)
      Call ErrMsg('RsoInt',005)
      Call Quit(16)
995   Call ErrMsg('RsoInt',006)
      Call Quit(16)
      end
