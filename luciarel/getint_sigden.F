* Fetch integral from integral list according to SIGDEN
* excitation class formalism - driver routine
*
      subroutine getint_sigden(YINT,INDEX,IRI,ITOTDIM,
     &                         NASHT,F1,F2,IPRNT)
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "crun.inc"
#include "symm.inc"
*
      dimension INDEX(4)
*
      if (INTIMP.eq.5) then
        write(6,*) 'MOLCAS sigden import under development.'
        write(6,*) 'Please be patient!'
        call abend2('Quitting.')
C FIXME: Write me
c       call getinm_sigden()
*
      else if (INTIMP.eq.6) then
* MOLFDIR format Kramers reduced
        if (IRI.eq.1) then
          IRIOFF1 = 0
          IRIOFF2 = 0
        else if (IRI.eq.2) then
          IRIOFF1 = (2*ITOTDIM)**2
          IRIOFF2 = (2*ITOTDIM)**4
        end if
*
        call getind_sigden(YINT,INDEX,IRI,ITOTDIM,
     &                     WORK(KNUMINT),WORK(KKR),WORK(KINDK),
     &                     WORK(KINDL),WORK(KINVREO),WORK(KIOFF12),
     &                     WORK(KDITWOE+IRIOFF2),WORK(KDIONEE+IRIOFF1),
     &                     IPRNT)
      else if (INTIMP.eq.7) then
*
* Define the KR array (which is usually read from the MOLFDIR file
* in standalone calculations!)
        call kkr_set(WORK(KKR),ITOTDIM,SYDI,NIRR_PN,INVSYM)
*
* MOLFDIR format full spinor list
        call getind_molf(YINT,INDEX,IRI,ITOTDIM,WORK(KINVREO),
     &                   WORK(KKR),NASHT,F1,F2,IPRNT)
      end if
*
      return
      end
