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
*
#include "mxpdim.inc"
#include "symm.inc"
*
      dimension IDEKS(MXPORB),IQISYM(MXPIRR)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      if (NTEST.ge.1) then
        write(6,*)
        write(6,*) '=========================='
        write(6,*) ' numints_red speaking:    '
        write(6,*) '=========================='
        if (NTEST.ge.5) then
          write(6,*) 'Point group multiplication table:'
          call iwrtma(IPGMULT,NIRR_PN,NIRR_PN,MXNSPIR,MXNSPIR)
          write(6,*) 'Symmetry dimension array:'
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
        write(6,*) 'Symmetry packing labels:'
        call iwrtma(IQISYM,NOCC,1,NOCC,1)
      end if
*
* Number of symmetry blocks and total number of integrals
      INIT = 1
      do ISP=1,NIRR_PN,1
        if (NTEST.ge.10) write(6,*) 'ISP =',ISP
        NSP=IQISYM(ISP)
        if (NSP.eq.0) then
          write(6,*) 'Number of symm for P is  ',NSP
          write(6,*) 'Skipping two-el. int. read.'
          goto 333
        end if
        if (NTEST.ge.10) write(6,*) 'NSP =',NSP
        NOP=SYDI(NSP)
        do ISQ=1,ISP,1
          if (NTEST.ge.10) write(6,*) 'ISQ =',ISQ
          NSQ = IQISYM(ISQ)
          NOQ = SYDI(NSQ)
          ISPQ = IPGMULT(NSP,NSQ)
          if (NTEST.ge.10) write(6,*) 'NSQ =',NSQ
          do ISR=1,ISP,1
            if (NTEST.ge.10) write(6,*) 'ISR =',ISR
            NSR = IQISYM(ISR)
            NOR = SYDI(NSR)
            ISPQR = IPGMULT(ISPQ,NSR)
            ISSM=ISR
            if (ISR.EQ.ISP) ISSM=ISQ
            if (NTEST.ge.10) write(6,*) 'NSR =',NSR
            do 40 ISS=1,ISSM,1
              if (NTEST.ge.10) write(6,*) 'ISS =',ISS
              NSS=IQISYM(ISS)
              NOS=SYDI(NSS)
              if (NTEST.ge.10) then
                write(6,'(A,4I6)') 'NOP,NOQ,NOR,NOS',NOP,NOQ,NOR,NOS
                write(6,'(A,2I3)') 'ISPQ,ISPQR',ISPQ,ISPQR
              end if
              if (NTEST.ge.10) write(6,*) 'NSS =',NSS
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
                 write(6,601) NSP,NSQ,NSR,NSS,NIGL
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
