*  Generate index reordering matrices depending on
*  INDREO(Input index) = output index
*  IREOBACK(Output index) = input index
*  g functions and higher angular momenta not (yet) implemented !!
*
      subroutine gen_reo_bas(IDOUGRP,INDREO,IREOBACK,
     &                       NSAOBAS,NMAXBAS,NFRO,
     &                       IPRINT)
*
      implicit real*8 (A-H,O-Z)
*
      dimension INDREO(NMAXBAS),NSAOBAS(4),
     &          IREOBACK(NMAXBAS),NFRO(8)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      if (IDOUGRP.le.3) then
        write(6,*) 'Point group for index reordering not '
        write(6,*) 'implemented yet. '
        write(6,*) 'Double group label is ',IDOUGRP
        Call Abend2( 'Quitting.' )
      else if (IDOUGRP.eq.4) then
*---------------------------------
*  Group C_2h
*---------------------------------
*  s functions
        do IS=1,NSAOBAS(1),1
          INDREO(IS) = IS
          IREOBACK(IS) = IS
        end do
        ICT = NSAOBAS(1)
*  d+2,d0,d-2 functions
        do IDCOMP=1,3,1
          IDCT = IDCOMP
          ISUM = 0
          do I=1,IDCOMP-1,1
            ISUM = ISUM + NSAOBAS(3)
          end do
          do ID=ICT+1+ISUM,ICT+NSAOBAS(3)+ISUM,1
            INDREO(ID) = IDCT + ICT
            IREOBACK(IDCT+ICT) = ID
            IDCT = IDCT + 3
          end do
        end do
        ICT = ICT + 3*NSAOBAS(3)
        ICTSAVE = ICT
*
        IFROSUM = 0
        do ISYM=1,1,1
          IFROSUM = IFROSUM + NFRO(ISYM)
        end do
        do I=1,ICT-NFRO(1),1
          IREOBACK(I) = IREOBACK(I+NFRO(1)) - NFRO(1)
        end do
*
*  d+,d- functions
        do IDCOMP=1,2,1
          IDCT = IDCOMP
          ISUM = 0
          do I=1,IDCOMP-1,1
            ISUM = ISUM + NSAOBAS(3)
          end do
          do ID=ICT+1+ISUM,ICT+NSAOBAS(3)+ISUM,1
            INDREO(ID) = IDCT + ICT
            IREOBACK(IDCT+ICT) = ID
            IDCT = IDCT + 2
          end do
        end do
        ICT = ICT + 2*NSAOBAS(3)
        ICTSAVE = ICT
*
        IFROSUM = 0
        IFROPRE = 0
        do ISYM=1,2,1
          IFROSUM = IFROSUM + NFRO(ISYM)
        end do
        do ISYMM=1,1,1
          IFROPRE = IFROPRE + NFRO(ISYMM)
        end do
        do I=ICTSAVE-IFROPRE+1,ICT-IFROSUM,1
          IREOBACK(I) = IREOBACK(I+IFROPRE) - IFROPRE
        end do
*
*  p_0 functions
        do IP=ICT+1,ICT+NSAOBAS(2),1
          INDREO(IP) = IP
          IREOBACK(IP) = IP
        end do
        ICT = ICT + NSAOBAS(2)
*  f-2,f0,f+2 functions
        do IFCOMP=1,3,1
          IFCT = IFCOMP
          ISUM = 0
          do I=1,IFCOMP-1,1
            ISUM = ISUM + NSAOBAS(4)
          end do
          do IF=ICT+1+ISUM,ICT+NSAOBAS(4)+ISUM,1
            INDREO(IF) = IFCT + ICT
            IREOBACK(IFCT+ICT) = IF
            IFCT = IFCT + 3
          end do
        end do
        ICT = ICT + 3*NSAOBAS(4)
        ICTSAVE = ICT
*
        IFROSUM = 0
        IFROPRE = 0
        do ISYM=1,3,1
          IFROSUM = IFROSUM + NFRO(ISYM)
        end do
        do ISYMM=1,2,1
          IFROPRE = IFROPRE + NFRO(ISYMM)
        end do
        do I=ICTSAVE-IFROPRE+1,ICT-IFROSUM,1
          IREOBACK(I) = IREOBACK(I+IFROPRE) - IFROPRE
        end do
*
*  p+,p- functions
        do IPCOMP=1,2,1
          IPCT = IPCOMP
          ISUM = 0
          do I=1,IPCOMP-1,1
            ISUM = ISUM + NSAOBAS(2)
          end do
          do IP=ICT+1+ISUM,ICT+NSAOBAS(2)+ISUM,1
            INDREO(IP) = IPCT + ICT
            IREOBACK(IPCT+ICT) = IP
            IPCT = IPCT + 2
          end do
        end do
        ICT = ICT + 2*NSAOBAS(2)
*  f-3,f-1,f+1,f+3 functions
        do IFCOMP=1,4,1
          IFCT = IFCOMP
          ISUM = 0
          do I=1,IFCOMP-1,1
            ISUM = ISUM + NSAOBAS(4)
          end do
          do IF=ICT+1+ISUM,ICT+NSAOBAS(4)+ISUM,1
            INDREO(IF) = IFCT + ICT
            IREOBACK(IFCT+ICT) = IF
            IFCT = IFCT + 4
          end do
        end do
        ICT = ICT + 4*NSAOBAS(4)
*
        IFROSUM = 0
        IFROPRE = 0
        do ISYM=1,4,1
          IFROSUM = IFROSUM + NFRO(ISYM)
        end do
        do ISYMM=1,3,1
          IFROPRE = IFROPRE + NFRO(ISYMM)
        end do
        do I=ICTSAVE-IFROPRE+1,ICT-IFROSUM,1
          IREOBACK(I) = IREOBACK(I+IFROPRE) - IFROPRE
        end do
*
      else if (IDOUGRP.eq.5) then
*---------------------------------
*  Group C_2
*---------------------------------
*  s functions
        do IS=1,NSAOBAS(1),1
          INDREO(IS) = IS
          IREOBACK(IS) = IS
        end do
        ICT = NSAOBAS(1)
*  p_0 functions
        do IP=ICT+1,ICT+NSAOBAS(2),1
          INDREO(IP) = IP
          IREOBACK(IP) = IP
        end do
        ICT = ICT + NSAOBAS(2)
*  d+2,d0,d-2 functions
        do IDCOMP=1,3,1
          IDCT = IDCOMP
          ISUM = 0
          do I=1,IDCOMP-1,1
            ISUM = ISUM + NSAOBAS(3)
          end do
          do ID=ICT+1+ISUM,ICT+NSAOBAS(3)+ISUM,1
            INDREO(ID) = IDCT + ICT
            IREOBACK(IDCT+ICT) = ID
            IDCT = IDCT + 3
          end do
        end do
        ICT = ICT + 3*NSAOBAS(3)
*  f-2,f0,f+2 functions
        do IFCOMP=1,3,1
          IFCT = IFCOMP
          ISUM = 0
          do I=1,IFCOMP-1,1
            ISUM = ISUM + NSAOBAS(4)
          end do
          do IF=ICT+1+ISUM,ICT+NSAOBAS(4)+ISUM,1
            INDREO(IF) = IFCT + ICT
            IREOBACK(IFCT+ICT) = IF
            IFCT = IFCT + 3
          end do
        end do
        ICT = ICT + 3*NSAOBAS(4)
        ICTSAVE = ICT
*
        IFROSUM = 0
        do ISYM=1,1,1
          IFROSUM = IFROSUM + NFRO(ISYM)
        end do
        do I=1,ICT-NFRO(1),1
          IREOBACK(I) = IREOBACK(I+NFRO(1)) - NFRO(1)
        end do
*
*  p+,p- functions
        do IPCOMP=1,2,1
          IPCT = IPCOMP
          ISUM = 0
          do I=1,IPCOMP-1,1
            ISUM = ISUM + NSAOBAS(2)
          end do
          do IP=ICT+1+ISUM,ICT+NSAOBAS(2)+ISUM,1
            INDREO(IP) = IPCT + ICT
            IREOBACK(IPCT+ICT) = IP
            IPCT = IPCT + 2
          end do
        end do
        ICT = ICT + 2*NSAOBAS(2)
*  d+,d- functions
        do IDCOMP=1,2,1
          IDCT = IDCOMP
          ISUM = 0
          do I=1,IDCOMP-1,1
            ISUM = ISUM + NSAOBAS(3)
          end do
          do ID=ICT+1+ISUM,ICT+NSAOBAS(3)+ISUM,1
            INDREO(ID) = IDCT + ICT
            IREOBACK(IDCT+ICT) = ID
            IDCT = IDCT + 2
          end do
        end do
        ICT = ICT + 2*NSAOBAS(3)
*  f-3,f-1,f+1,f+3 functions
        do IFCOMP=1,4,1
          IFCT = IFCOMP
          ISUM = 0
          do I=1,IFCOMP-1,1
            ISUM = ISUM + NSAOBAS(4)
          end do
          do IF=ICT+1+ISUM,ICT+NSAOBAS(4)+ISUM,1
            INDREO(IF) = IFCT + ICT
            IREOBACK(IFCT+ICT) = IF
            IFCT = IFCT + 4
          end do
        end do
        ICT = ICT + 4*NSAOBAS(4)
*
        IFROSUM = 0
        do ISYM=1,2,1
          IFROSUM = IFROSUM + NFRO(ISYM)
        end do
        do I=ICTSAVE-NFRO(1)+1,ICT-IFROSUM,1
          IREOBACK(I) = IREOBACK(I+NFRO(1)) - NFRO(1)
        end do
*
      else if (IDOUGRP.eq.6) then
*---------------------------------
*  Group C_s
*---------------------------------
      else if (IDOUGRP.eq.7) then
*---------------------------------
*  Group C_i
*---------------------------------
      else if (IDOUGRP.eq.8) then
*---------------------------------
*  Group C_1
*---------------------------------

      end if
*
      if (NTEST.ge.1) then
        write(6,*)
        write(6,*) 'Reorder matrices for SAO indices: '
        write(6,*) 'Row: input index, Column: Irrep   '
        call iwrtma(INDREO,NMAXBAS,1,NMAXBAS,1)
        write(6,*)
        write(6,*) 'Reorder back matrices for SAO indices (fock): '
        write(6,*) 'Row: input index, Column: Irrep   '
        call iwrtma(IREOBACK,NMAXBAS,1,NMAXBAS,1)
      end if
*
      return
      end
