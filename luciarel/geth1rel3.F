*  get desired block of one-electron integrals
*
*
*
*================================================
*================================================
*     CAUTION !!!!!
*
* At the moment, geth1rel is only designed for an
* initial CI calculation. As soon as the transformed
* integral matrices contain mixed off-diagonal and
* real/imaginary contributions, the summation has
* to be expanded to all lists !!!!
*================================================
*================================================
*
*
      subroutine geth1rel3(SPININTS,
     &                     H,ISM,ITYP,JSM,JTYP,IRISIG,IRICO,IPM1,
     &                     IJUB,IRETURN,SCALIN,INTSPACE,ENVIRO)
*
*-----------------------------------
*   Timo Fleig
*-----------------------------------
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
#include "onel.inc"
*
      character*72 ENVIRO
*
      dimension SPININTS(INTSPACE)
      dimension H(*)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      if (NTEST.ge.3) then
        write(6,*)
        write(6,*) 'Checking input in geth1rel3 '
        write(6,*) 'ISM       ',ISM
        write(6,*) 'JSM       ',JSM
        write(6,*) 'ITYP      ',JTYP
        write(6,*) 'IRISIG    ',IRISIG
        write(6,*) 'IRICO     ',IRICO
        write(6,*) 'IPM1      ',IPM1
        write(6,*) 'IJUB      ',IJUB
      end if
*
      ITRALINE = INTSPACE / 32
      ITRACHK = ITRALINE*32
      if (ITRACHK.ne.INTSPACE) then
         write(6,*) 'Problem with int space in geth1rel3.'
         Call Abend2( 'quitting' )
      end if
*
      IOFFTRA(1) = 1
      do I=2,32,1
         IOFFTRA(I) = IOFFTRA(I-1) + ITRALINE
      end do
*
      call findhigh(NGSOBPI,MXPIRR*MXPNGAS,MXGSDIM)
      if (NTEST.ge.5) then
         write(6,*) 'Highest GAS dimension is  ',MXGSDIM
      end if
*
      I_AM_OKAY = 1
      if (ISM.ne.JSM) then
         write(6,*) 'No integrals in this symmetry block'
         I_AM_OKAY = 0
         goto 100
      else
         if (NTEST.eq.5) then
            write(6,*) 'ISM equals JSM. I am fetching integrals'
         end if
*
*  determine real or imaginary integral parts
         SCALIN = 1.D0
         if (IRISIG.eq.1.and.IRICO.eq.1) then
            IRIINT = 1
         else if (IRISIG.eq.1.and.IRICO.eq.2) then
            IRIINT = 2
            SCALIN = -1.D0
         else if (IRISIG.eq.2.and.IRICO.eq.1) then
            IRIINT = 2
         else if (IRISIG.eq.2.and.IRICO.eq.2) then
            IRIINT = 1
         else
            write(6,*) 'IRI problem in geth1rel - returning'
            Call Abend()
         end if
*
         if (ISUMNEW.eq.1.and.ENVIRO(1:6).ne.'DIRAC ') then
            if (IOPER.eq.1) then
               continue
            else if (IOPER.eq.2) then
               do LIST=1,NIRR_DG,1
                  call evenodd(IEOLIST,LIST)
                  do I=1,ITRALINE,1
                     if (IEOLIST.eq.1) then
                        IRUN = I + IOFFTRA(LIST+8) - 1
                        IRUN1 = I + IOFFTRA(LIST+16) - 1
                        SPININTS(IRUN) =
     &                     SPININTS(IRUN1)
                     else if (IEOLIST.eq.2) then
                        IRUN = I + IOFFTRA(LIST+8) - 1
                        IRUN2 = I + IOFFTRA(LIST+24) - 1
                        SPININTS(IRUN) =
     &                     SPININTS(IRUN)  +
     &                     SPININTS(IRUN2)
                     end if
                  end do
               end do
            else if (IOPER.eq.3) then
               do LIST=1,NIRR_DG,1
                  call evenodd(IEOLIST,LIST)
                  do I=1,ITRALINE,1
                     if (IEOLIST.eq.1) then
                        IRUN = I + IOFFTRA(LIST) - 1
                        IRUN2 = I + IOFFTRA(LIST+16) - 1
                        SPININTS(IRUN) =
     &                        SPININTS(IRUN) +
     &                        SPININTS(IRUN2)
                     else if (IEOLIST.eq.2) then
                        IRUN = I + IOFFTRA(LIST) - 1
                        IRUN1 = I + IOFFTRA(LIST+8) - 1
                        IRUN3 = I + IOFFTRA(LIST+24) - 1
                        SPININTS(IRUN) =
     &                        SPININTS(IRUN1) +
     &                        SPININTS(IRUN3)
                     end if
                  end do
               end do
            end if
         end if
*
*  Determine integral list
         INTLIST = 2*(ISM-NIRR_PN) + IRIINT - 2
         if (IOPER.eq.2) INTLIST = INTLIST + 8
*  List offset
         LSTOFF = IOFFTRA(INTLIST)
*  Offset for ITYP (inner index)
         ITPOFF = 0
         do I=1,ITYP-1,1
            ITPOFF = ITPOFF + 4*MXGSDIM**2
         end do
*  Offset for JTYP (outer index)
         JTPOFF = 0
         do J=1,JTYP-1,1
            JTPOFF = JTPOFF + 4*MXGSDIM**2*NGASPI
         end do
*  Final offset
         IFINOFF = LSTOFF + JTPOFF + ITPOFF
*
         if (IPM1.ne.0.and.IJUB.ne.0) then
            write(6,*) 'Combination of IJUB and IPM1 impossible'
            write(6,*) 'I am quitting'
            Call Abend()
         end if
*
         if (NTEST.ge.5) then
           write(6,*)
           write(6,*) 'List and offsets in geth1rel3 '
           write(6,*) 'INTLIST    .... ',INTLIST
           write(6,*) 'LSTOFF     .... ',LSTOFF
           write(6,*) 'ITPOFF     .... ',ITPOFF
           write(6,*) 'JTPOFF     .... ',JTPOFF
           write(6,*) 'IFINOFF    .... ',IFINOFF
         end if
*
*  Determine relevant offsets and copy desired block to H
*
*  ISEARCH in the off-diagonal case depends on the symmetry of
*  the desired spinors. It has to be permuted, if the spinor
*  irrep is within the second half of the fermion irreps:
*
         if (IPM1.eq.1) then
            ISEARCH = 2
         else if (IPM1.eq.-1) then
            ISEARCH = 3
         else if (IPM1.eq.0) then
            if (IJUB.eq.1) then
               ISEARCH = 1
            else if (IJUB.eq.2) then
               ISEARCH = 4
            else
               write(6,*) 'IJUB out of range. Quitting.'
               Call Abend()
            end if
         else
            write(6,*) 'IPM1 out of range. Quitting.'
            Call Abend()
         end if
         ILEN = IPUBPGS(ISEARCH,ITYP,JTYP,INTLIST)
*
         ISTART = 1
         do I=1,ISEARCH-1,1
            ISTART = ISTART + IPUBPGS(I,ITYP,JTYP,INTLIST)
         end do
         IEND = ISTART + ILEN - 1
         K = 1
         ISTART = ISTART + IFINOFF - 1
         IEND = IEND + IFINOFF - 1
*
         if (NTEST.ge.10) then
           write(6,*) 'ISEARCH block  ',ISEARCH
           write(6,*) 'Block length   ',ILEN
           write(6,*) 'Start address  ',ISTART
           write(6,*) 'End address    ',IEND
         end if
*
* Now do the work ...
         do J=ISTART,IEND,1
            H(K) = SPININTS(J)
            K = K + 1
         end do
*
         if (NTEST.ge.10) then
            write(6,*) 'The integrals are:'
             do L=1,ILEN,1
                write(6,*) H(L)
             end do
            write(6,*) 'ISM is  ',ISM
            write(6,*) 'JSM is  ',JSM
            write(6,*) 'ITYP is  ',ITYP
            write(6,*) 'JTYP is  ',JTYP
            write(6,*) 'I have fetched integrals successfully'
            write(6,*) 'and am now returning'
         end if
         IRETURN = 1
      end if
100   continue
      if (I_AM_OKAY.eq.0) then
         write(6,*) 'I am returning without integrals'
         IRETURN = 0
      end if
*
* Run variable for summation of integrals. Only in the first call
* for geth1rel3, the summation should be carried out.
      ISUMNEW = ISUMNEW + 1
      if (NTEST.ge.5)
     &  write(6,*) 'ISUMNEW has current value  ',ISUMNEW
*
500   return
      end
