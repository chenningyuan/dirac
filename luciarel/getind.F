*
* Get 2-electron integrals from DIRAC environment structure
*
      subroutine getind(XOUT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  INVREO,KR,NUMINT,
     &                  INDK,INDL,IOFF12,
     &                  XIN)
*
      implicit real*8 (A-H,O-Z)
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
