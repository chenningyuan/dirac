*  Arrange AO integral matrix for transformation.
*  Read integrals are on AOTRISM in packed triangle form,
*  all symmetries subsequently. This routine picks the
*  desired blocks and unpacks them to AOINT matrix.
*  This matrix finally is transformed by oneltra.
*
*!!!! CAUTION !!!
*!!!! CAUTION !!!
*
*  If spiorb = true, then HELPMAT2 contains the input
*  integrals. Otherwise, i.e. for spin-free integrals,
*  the input integrals are on AOTRISM
*  (This is terrible, but pragmatic ...)
*
*!!!! CAUTION !!!
*!!!! CAUTION !!!
*
      subroutine intsel(AOINT,AOTRISM,HELPMAT2,
     &                  PROD,IPRODDIM,
     &                  IPGIDIM,IDGIDIM,ICURRDIM,IDIMGRP,NIRR_PN,
     &                  ITRIA,ITRIASUM,ITRIADB,ITRIADBSUM,
     &                  SPIORB,IDBG,NTRIAEL,ENVIRO,COORD)
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "input.inc"
*
      logical SPIORB
      character*72 ENVIRO
      character*1 COORD
*
      dimension IPGIDIM(8),IDGIDIM(2),ITRIA(8),ITRIADB(4)
      dimension AOINT(ICURRDIM,ICURRDIM)
      dimension AOTRISM(ITRIASUM),
     &          HELPMAT2(NTRIAEL),PROD(IPRODDIM)
*
*  Spin-free or spin-orbit integrals:
*  The whole stuff only has to be done if we have symmetry
*  blocks of integrals which have to be picked and expanded.
*  In the SO case, this is obsolete.
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      ZERO = 0.0D0
*  Erase AOINT (just to make sure)
      do ICOL=1,IDGIDIM(IDIMGRP),1
        call setvec(AOINT(1,ICOL),ZERO,IDGIDIM(IDIMGRP))
      end do
*
*---------------------------------------------------------------
      if (SPIORB) goto 500
*
*  Processing spin-orbit integrals (HELPMAT2 contains input ints)
*---------------------------------------------------------------
*
      if (NTEST.gt.2) then
      write(6,*) 'Checking matrix AOTRISM ...'
         do I=1,ITRIASUM,1
            write(6,'(1G18.12)') AOTRISM(I)
         end do
      end if
*
*  Erase matrix HELPMAT2 (it has been used in rdcoef)
      do I=1,NTRIAEL,1
         HELPMAT2(I) = ZERO
      end do
*
*  Determine offset for relevant symmetry block(s)
      ISTART = 1
      if (IDIMGRP.eq.2) then
         do J=1,2,1
            ISTART = ISTART + ((IPGIDIM(J)**2 + IPGIDIM(J)) /2)
         end do
         if (IDBG.eq.7)
     &       ISTART = 1 + ((IPGIDIM(1)**2 + IPGIDIM(1)) /2)
      end if
C?    write(6,*) 'ISTART is    ',ISTART
*
*  Expand to full triangle matrix
      ISYM = 2*IDIMGRP - 1
      if (IDBG.eq.7) ISYM = IDIMGRP
      IEND = ISTART + ((IPGIDIM(ISYM)**2 + IPGIDIM(ISYM)) /2) - 1
C?    write(6,*) 'IEND is    ',IEND
      if (IDIMGRP.eq.2) then
         IENDEND = ITRIASUM
      else if (IDIMGRP.eq.1) then
         IENDEND = 0
         do J=1,2,1
            IENDEND = IENDEND + ((IPGIDIM(J)**2 + IPGIDIM(J)) /2)
         end do
         if (IDBG.eq.7)
     &       IENDEND = ((IPGIDIM(1)**2 + IPGIDIM(1)) /2)
      end if
      NSKIP = IPGIDIM(ISYM)
*
      if (IDBG.ge.4.and.IDBG.le.6) then
        J = 1
        do I=ISTART,IEND,1
          HELPMAT2(J) = AOTRISM(I)
          J = J + 1
        end do
*
        INSYM = IPGIDIM(ISYM+1)
        KOUT = J
        INP = IEND + 1
        do ILINE=1,INSYM,1
          KOUT = KOUT + NSKIP
          do I=INP,INP+ILINE-1,1
            HELPMAT2(KOUT) = AOTRISM(I)
            if (NTEST.gt.2) then
              write(6,*) 'HELPMAT2(KOUT) is    ',
     &                    HELPMAT2(KOUT),KOUT
              write(6,*) 'AOTRISM(I) is    ',AOTRISM(I),I
            end if
            KOUT = KOUT + 1
          end do
          INP = INP + ILINE
        end do
      else if (IDBG.ge.7.and.IDBG.le.8) then
        if (IDIMGRP.eq.1) then
          J = 1
          do I=ISTART,IEND,1
            HELPMAT2(J) = AOTRISM(I)
            J = J + 1
          end do
        else if (IDIMGRP.eq.2) then
          INSYM = IPGIDIM(2)
          KOUT = 1
          do I=ISTART,IEND,1
            HELPMAT2(KOUT) = AOTRISM(I)
            if (NTEST.gt.2) then
              write(6,*) 'HELPMAT2(KOUT) is    ',
     &                    HELPMAT2(KOUT),KOUT
              write(6,*) 'AOTRISM(I) is    ',AOTRISM(I),I
            end if
            KOUT = KOUT + 1
          end do
        end if
      end if
*
 100  continue
*
*  Check matrix HELPMAT2:
      if (NTEST.gt.2) then
      write(6,*) 'Checking matrix HELPMAT2:'
         do I=1,NTRIAEL,1
            write(6,'(1G18.12)') HELPMAT2(I)
         end do
      end if
*
*------------------------------------------
 500  continue
*------------------------------------------
*
*  Unfold to full matrix for transformation
      if (SPIORB) then
*
        if (NTEST.ge.5) then
          write(6,*) 'Checking input SO integrals in intsel:'
          call wrtmat(HELPMAT2,NTRIAEL,1,NTRIAEL,1)
        end if
*
*  Erase PROD (just to make sure)
        call setvec(PROD,ZERO,IPRODDIM)
*
        if (IDIMGRP.eq.2) then
           if (ENVIRO(1:6).eq.'MOLCAS') then
*  Pick second dimension group:
*  Check for pathological cases ...
*  (See LUCIAREL notes for derivation of these formulae!)
             ILENREQ = IPGIDIM(3)*IPGIDIM(4)
             if (ILENREQ.gt.IPRODDIM) then
                write(6,*) 'Problem in intsel (SO)'
                write(6,*) 'Required length:  ',ILENREQ
                write(6,*) 'Length usable: ',IPRODDIM
                Call Abend2( 'Quitting.' )
             end if
*
             if (IDBG.ge.4.and.IDBG.le.6) then
               J = 1
               INP = ITRIADB(1) + IDGIDIM(1) + 1
               do ISYM=1,IDGIDIM(2),1
                do I=INP,INP+ISYM-1,1
                     PROD(J) = HELPMAT2(I)
                     J = J + 1
                  end do
                  INP = INP + ISYM + IPGIDIM(1) + IPGIDIM(2)
               end do
               IGLOFAC = +1
               call unfold(AOINT(1,1),PROD(1),IDGIDIM(IDIMGRP),
     &                     -1,ITRIADB(IDIMGRP),IGLOFAC)
             else if (IDBG.eq.7) then
               J = 1
               INP = ITRIA(1) + IPGIDIM(1) + 1
               do ISYM=1,IPGIDIM(2),1
                 do I=INP,INP+ISYM-1,1
                   PROD(J) = HELPMAT2(I)
                   J = J + 1
                 end do
                 INP = INP + ISYM + IPGIDIM(1)
               end do
               IGLOFAC = +1
               call unfold(AOINT(1,1),PROD(1),IDGIDIM(2),
     &                     -1,ITRIA(2),IGLOFAC)
               if (COORD.eq.'Z') then
                 do ICOL=1,IPGIDIM(2),1
                   call copvec(AOINT(1,ICOL),
     &               AOINT(IPGIDIM(2)+1,IPGIDIM(2)+ICOL),IPGIDIM(2))
                   call setvec(AOINT(1,IPGIDIM(2)+ICOL),
     &                         ZERO,IPGIDIM(2))
                   call setvec(AOINT(IPGIDIM(2)+1,ICOL),
     &                         ZERO,IPGIDIM(2))
                 end do
               else if (COORD.eq.'X'.or.COORD.eq.'Y') then
                 do ICOL=1,IPGIDIM(2),1
                   call copvec(AOINT(1,ICOL),
     &               AOINT(1,IPGIDIM(2)+ICOL),IPGIDIM(2))
                   call copvec(AOINT(1,ICOL),
     &               AOINT(IPGIDIM(2)+1,ICOL),IPGIDIM(2))
                   call setvec(AOINT(1,ICOL),ZERO,IPGIDIM(2))
                   call setvec(AOINT(IPGIDIM(2)+1,IPGIDIM(2)+ICOL),
     &                         ZERO,IPGIDIM(2))
                 end do
               end if
             end if
*
           else if (ENVIRO(1:6).eq.'MOLECU') then
*  Pick blocks from MOLECULE ordering (Ag,Au,Bu,Bg)
             IAUBUDIM = IPGIDIM(3) + IPGIDIM(4)
             IAUSTA = ITRIA(1) + IPGIDIM(1) + 1
             do I=1,IAUBUDIM,1
               do II=1,I,1
                 PROD(J) = HELPMAT2(IAUSTA+II-1)
                 J = J + 1
               end do
               IAUSTA = IAUSTA + IPGIDIM(1) + I
             end do
           end if
*          ^ end ENVIRO
        else if (IDIMGRP.eq.1) then
           if (ENVIRO(1:6).eq.'MOLCAS') then
             if (IDBG.ge.4.and.IDBG.le.6) then
               IGLOFAC = +1
               call unfold(AOINT(1,1),HELPMAT2(1),IDGIDIM(IDIMGRP),
     &                     -1,ITRIADB(IDIMGRP),IGLOFAC)
             else if (IDBG.ge.7.and.IDBG.le.8) then
               IGLOFAC = +1
               if (COORD.eq.'Z') then
                 call vectopart(AOINT(1,1),HELPMAT2(1),IDGIDIM(1),
     &                          -1,ITRIA(1),IGLOFAC,IPGIDIM(1),1,1)
                 call vectopart(AOINT(1,1),HELPMAT2(1),IDGIDIM(1),
     &                          -1,ITRIA(1),IGLOFAC,IPGIDIM(1),
     &                          IPGIDIM(1)+1,IPGIDIM(1)+1)
               else if (COORD.eq.'X'.or.COORD.eq.'Y') then
                 call vectopart(AOINT(1,1),HELPMAT2(1),IDGIDIM(1),
     &                          -1,ITRIA(1),IGLOFAC,IPGIDIM(1),
     &                          IPGIDIM(1)+1,1)
                 call vectopart(AOINT(1,1),HELPMAT2(1),IDGIDIM(1),
     &                          -1,ITRIA(1),IGLOFAC,IPGIDIM(1),
     &                          1,IPGIDIM(1)+1)
               end if
             end if
*
           else if (ENVIRO(1:6).eq.'MOLECU') then
             I3TRIADIM = IPGIDIM(1) + IPGIDIM(3) + IPGIDIM(4)
             I3TRIA = (I3TRIADIM**2 + I3TRIADIM) / 2
             IBGSTA = I3TRIA + 1
             IBGDIM = IPGIDIM(2)
             IAGDIM = IPGIDIM(1)
             IAGEND = ITRIA(1)
             J = 1
             do I=1,IBGDIM,1
               do K=1,IAGDIM,1
                 PROD(J) = HELPMAT2(IBGSTA+K-1)
                 J = J + 1
               end do
               do II=1,I,1
                 PROD(J) = HELPMAT2(IBGSTA+I3TRIADIM+II-1)
                 J = J + 1
               end do
             end do
             K = 1
             do JJ=1,J-1,1
               HELPMAT2(IAGEND+K) = PROD(JJ)
               K = K + 1
             end do
             IGLOFAC = +1
             call unfold(AOINT(1,1),HELPMAT2(1),IDGIDIM(IDIMGRP),
     &                   -1,ITRIADB(IDIMGRP),IGLOFAC)
           end if
*          ^ end ENVIRO
        else
           write(6,*) 'Problem in intsel.'
           write(6,*) 'IDIMGRP is  ',IDIMGRP
           Call Abend2( 'Quitting.' )
        end if
*       ^ end IDIMGRP
*
*  IGLOFAC is the global factor for l integrals, as pure r x nabla
*  is calculated by genmfint/amfi. The l integrals should be
*  h/i = -ih, so assume imaginary integrals and multiply all
*  SO ints with -1: no, we're not doing that right now ...
*
      else
         IGLOFAC = 1
         call unfold(AOINT(1,1),HELPMAT2(1),IDGIDIM(IDIMGRP),1,
     &               ITRIA(IDIMGRP),IGLOFAC)
         if (IDBG.ge.7.and.IDBG.le.8) then
           do ICOL=1,IPGIDIM(IDIMGRP),1
             call copvec(AOINT(1,ICOL),
     &              AOINT(IPGIDIM(IDIMGRP)+1,IPGIDIM(IDIMGRP)+ICOL),
     &              IPGIDIM(IDIMGRP))
           end do
         end if
      end if
*     ^ end SPIORB
*
      if (NTEST.ge.3) then
         write(6,*) 'Checking AOINT in intsel:'
         if (SPIORB) then
            write(6,*) 'These are spin-orbit integrals'
         else
            write(6,*) 'These are spin-free integrals'
         end if
         call wrtmat(AOINT,IDGIDIM(IDIMGRP),IDGIDIM(IDIMGRP),
     &               IDGIDIM(IDIMGRP),IDGIDIM(IDIMGRP))
      end if
*
      return
      end
