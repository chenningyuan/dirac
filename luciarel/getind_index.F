* Get 2-electron integral with selected indices
      subroutine getind_index(YOUT,IIN,JIN,KIN,LIN,NUMINT,
     &                        KR,INDK,INDL,IOFF12,IRI,XIN,IPRNT)
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "symm.inc"
      include 'crun.inc'
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
              if (IDGTYP.eq.2) then
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
