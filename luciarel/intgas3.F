*  Distribute integrals to partial lists according to
*  chosen GAS spaces
*
*  Order is symmetry - type at this stage
*  Within one block, the order is UU,UB,BU,BB
*  and within these blocks as in memory.
*
      subroutine intgas3(GASINT,
     &                   OINT,NGSOB,NGAS,NIRR_PN,INTLIST,
     &                   ITOTDIM,IDIM,ISMIND,IGSOFF,IGSIND,
     &                   IACSYM,IOUTLST,MXGSDIM,
     &                   ISPSF_MO,IDBG)
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "onel.inc"
#include "input.inc"
*
      dimension GASINT(4*MXGSDIM**2,NGAS,NGAS)
      dimension OINT(ITOTDIM,ITOTDIM),
     &          NGSOB(MXPIRR,MXPNGAS)
      dimension ISMIND(ITOTDIM),IGSIND(ITOTDIM,IACSYM:IACSYM+1)
      dimension IDIM(2),IGSOFF(2*NGAS),ISPSF_MO(MXNDGIRR,2)
      dimension IUBOFF(4,MXPNGAS,MXPNGAS),NFCS(4,MXPNGAS,MXPNGAS)
*
      NTESTL = 0
      NTEST = MAX(IPRINT,NTESTL)
*
      if (NTEST.ge.2) then
         write(6,*) 'This is intgas3.'
         write(6,*) 'Checking NGSOB matrix'
         do I=1,NIRR_PN,1
            write(6,*) (NGSOB(I,J),J=1,NGAS,1)
         end do
      end if
*
*  Calculate input offsets for start of GAS blocks
      JCT = 1
      IDGGASOFF(1,1) = 1
      ICOUNTER = 1
      IUB = 1
      do JGAS=2,NGAS,1
         if (IDBG.ge.4.and.IDBG.le.6) then
           IDGGASOFF(JGAS,IUB) = JCT + NGSOB(IACSYM+IUB-1,JGAS-1)
           JCT = JCT + NGSOB(IACSYM+IUB-1,JGAS-1)
         else if (IDBG.ge.7.and.IDBG.le.8) then
           IDGGASOFF(JGAS,IUB) = JCT + NGSOB(IOUTLST,JGAS-1)
           JCT = JCT + NGSOB(IOUTLST,JGAS-1)
         end if
         ICOUNTER = JGAS
      end do
      if (IDBG.ge.4.and.IDBG.le.6) then
        JCT = JCT + NGSOB(IACSYM+IUB-1,ICOUNTER)
      else if (IDBG.ge.7.and.IDBG.le.8) then
        JCT = JCT + NGSOB(IOUTLST,ICOUNTER)
      end if
      IDGGASOFF(1,2) = JCT
      IUB = 2
      do JGAS=2,NGAS,1
         if (IDBG.ge.4.and.IDBG.le.6) then
           IDGGASOFF(JGAS,IUB) = JCT + NGSOB(IACSYM+IUB-1,JGAS-1)
           JCT = JCT + NGSOB(IACSYM+IUB-1,JGAS-1)
         else if (IDBG.ge.7.and.IDBG.le.8) then
           IDGGASOFF(JGAS,IUB) = JCT + NGSOB(IOUTLST,JGAS-1)
           JCT = JCT + NGSOB(IOUTLST,JGAS-1)
         end if
      end do
* Check
      if (NTEST.ge.2) then
        write(6,*) 'Checking offset matrix IDGGASOFF(GAS,IUB) ...'
        call iwrtma(IDGGASOFF,NGAS,2,MXPNGAS,2)
      end if
*
*  Determine output offsets with respect to number of functions
*  in UB/GAS block
      IDBGSM = IOUTLST + NIRR_PN
      do IGAS=1,NGAS,1
        do JGAS = 1,NGAS,1
          IUBOFF(1,IGAS,JGAS) = 1
          do IIJUB = 1,4,1
            if (IIJUB.eq.1) then
             IUB = 1
             JUB = 1
            else if (IIJUB.eq.2) then
             IUB = 1
             JUB = 2
            else if (IIJUB.eq.3) then
             IUB = 2
             JUB = 1
            else
             IUB = 2
             JUB = 2
            end if
*
            IPGSM = ISPSF_MO(IDBGSM,IUB)
            JPGSM = ISPSF_MO(IDBGSM,JUB)
            NFCS(IIJUB,IGAS,JGAS) =
     &           NGSOB(IPGSM,IGAS) * NGSOB(JPGSM,JGAS)
*
            if (IIJUB.le.3) then
              IUBOFF(IIJUB+1,IGAS,JGAS) = IUBOFF(IIJUB,IGAS,JGAS)
     &                                    + NFCS(IIJUB,IGAS,JGAS)
            end if
          end do
        end do
      end do
*
      if (NTEST.ge.3) then
        do IGAS=1,NGAS,1
          do JGAS=1,NGAS,1
            do IIJUB=1,4,1
              write(6,*) 'IIJUB is  ',IIJUB
              write(6,'(A,X,I6,I3,I3,A,X,I6)')
     &                   'NFCS(',IIJUB,IGAS,JGAS,') = ',
     &                    NFCS(IIJUB,IGAS,JGAS)
              write(6,'(A,X,I6,I3,I3,A,X,I6)')
     &                   'IUBOFF(',IIJUB,IGAS,JGAS,') = ',
     &                    IUBOFF(IIJUB,IGAS,JGAS)
            end do
          end do
        end do
      end if
*
*  Distribute integrals
*
      call evenodd(IEOSYM,IOUTLST)
      IGSCT = 0
      IUBCT = 0
      if (IEOSYM.eq.1) then
*  Unbarred up left
        do JGAS=1,NGAS,1
           do IGAS=1,NGAS,1
              do IUB=1,2,1
                 do JUB=1,2,1
                    ISTART = IDGGASOFF(IGAS,IUB)
                    IEND = IDGGASOFF(IGAS,IUB) +
     &                      NGSOB(IOUTLST+IUB-1,IGAS) - 1
                    if (IDBG.ge.7.and.IDBG.le.8) then
                      IEND = IDGGASOFF(IGAS,IUB) +
     &                       NGSOB(IOUTLST,IGAS) - 1
                    end if
                    JSTART = IDGGASOFF(JGAS,JUB)
                    JEND = IDGGASOFF(JGAS,JUB) +
     &                      NGSOB(IOUTLST+JUB-1,JGAS) - 1
                    if (IDBG.ge.7.and.IDBG.le.8) then
                      JEND = IDGGASOFF(JGAS,JUB) +
     &                       NGSOB(IOUTLST,JGAS) - 1
                    end if
                    if (NTEST.ge.5) then
                      write(6,*) 'ISTART =  ',ISTART
                      write(6,*) 'IEND =  ',IEND
                      write(6,*) 'JSTART =  ',JSTART
                      write(6,*) 'JEND =  ',JEND
                    end if
*
                    if (IUB.eq.1.and.JUB.eq.1) NIJUB = 1
                    if (IUB.eq.1.and.JUB.eq.2) NIJUB = 2
                    if (IUB.eq.2.and.JUB.eq.1) NIJUB = 3
                    if (IUB.eq.2.and.JUB.eq.2) NIJUB = 4
                    IOUTOFF = IUBOFF(NIJUB,IGAS,JGAS)
                    if (NTEST.ge.5) then
                      write(6,*) 'Output offset is ',IOUTOFF
                    end if
*
                    do J=JSTART,JEND,1
                      do I=ISTART,IEND,1
                        GASINT(IOUTOFF+IUBCT,IGAS,JGAS)
     &                       = OINT(I,J)
                        if (NTEST.ge.5) then
                          write(6,*) 'Input integral ',OINT(I,J)
                          write(6,*) 'Output position ',IOUTOFF+IUBCT
                        end if
                        IUBCT = IUBCT + 1
                        IGSCT = IGSCT + 1
                      end do
                    end do
                    IPUBPGS(NIJUB,IGAS,JGAS,INTLIST)
     &                     = IUBCT
                    IUBCT = 0
                 end do
              end do
              if (NTEST.ge.5) then
                 do I=1,IGSCT,1
                    write(6,*) I,GASINT(I,IGAS,JGAS)
                 end do
                 write(6,*) '... end of block ...'
              end if
              IGSCT = 0
           end do
        end do
      else if (IEOSYM.eq.2) then
*  Barred up left: reorder!
        do JGAS=1,NGAS,1
           do IGAS=1,NGAS,1
              do IBU=2,1,-1
                 do JBU=2,1,-1
                    ISTART = IDGGASOFF(IGAS,IBU)
                    IEND = IDGGASOFF(IGAS,IBU) +
     &                      NGSOB(IACSYM+IBU-1,IGAS) - 1
                    if (IDBG.ge.7.and.IDBG.le.8) then
                      IEND = IDGGASOFF(IGAS,IBU) +
     &                       NGSOB(IOUTLST,IGAS) - 1
                    end if
                    JSTART = IDGGASOFF(JGAS,JBU)
                    JEND = IDGGASOFF(JGAS,JBU) +
     &                      NGSOB(IACSYM+JBU-1,JGAS) - 1
                    if (IDBG.ge.7.and.IDBG.le.8) then
                      JEND = IDGGASOFF(JGAS,JBU) +
     &                       NGSOB(IOUTLST,JGAS) - 1
                    end if
                    if (NTEST.ge.3) then
                      write(6,*) 'ISTART =  ',ISTART
                      write(6,*) 'IEND =  ',IEND
                      write(6,*) 'JSTART =  ',JSTART
                      write(6,*) 'JEND =  ',JEND
                    end if
*
                    if (IBU.eq.1.and.JBU.eq.1) NIJUB = 4
                    if (IBU.eq.1.and.JBU.eq.2) NIJUB = 3
                    if (IBU.eq.2.and.JBU.eq.1) NIJUB = 2
                    if (IBU.eq.2.and.JBU.eq.2) NIJUB = 1
                    IOUTOFF = IUBOFF(NIJUB,IGAS,JGAS)
                    if (NTEST.ge.3) then
                      write(6,*) 'Output offset is ',IOUTOFF
                    end if
*
                    do J=JSTART,JEND,1
                       do I=ISTART,IEND,1
                          GASINT(IOUTOFF+IUBCT,IGAS,JGAS)
     &                         = OINT(I,J)
                          if (NTEST.ge.5) then
                            write(6,*) 'Input integral ',OINT(I,J)
                            write(6,*) 'Output position ',IOUTOFF+IUBCT
                          end if
                          IUBCT = IUBCT + 1
                          IGSCT = IGSCT + 1
                       end do
                    end do
                    IPUBPGS(NIJUB,IGAS,JGAS,INTLIST)
     &                      = IUBCT
                    IUBCT = 0
                 end do
              end do
              if (NTEST.ge.2) then
                 do I=1,IGSCT,1
                    write(6,*) I,GASINT(I,IGAS,JGAS)
                 end do
                 write(6,*) '... end of block ...'
              end if
              IGSCT = 0
           end do
        end do
      end if
*
      if (NTEST.ge.2) then
         do IGAS=1,NGAS,1
            do JGAS=1,NGAS,1
               write(6,*) 'IPUBPGS for IGAS',IGAS,'and JGAS',JGAS
               call iwrtma(IPUBPGS(1,IGAS,JGAS,INTLIST),4,1,4,1)
            end do
         end do
      end if
*
      return
      end
