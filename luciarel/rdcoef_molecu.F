*  Get relevant block of coefficients from CSPAOSF(*) or
*  read coefficients from file 19
*
      subroutine rdcoef_molecu(ADJCO,COEF,CMOAO,HELPMAT,
     &                  ICHECK,ITOTDIM,ICRMXDIM,ICR2DIM,ICR3DIM,
     &                  IDIM,NIRR_PN,NIRR_DG,LU19,INTITER,
     &                  IMAXDIM,IOUTLST,IDOUGRP,ITRA_TYPE)
*------------------------------------------------------
* Modify according to ITRA_TYPE (like in rdcoef_molcas)
*------------------------------------------------------
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "coeff.inc"
#include "input.inc"
*
*   Coefficient common block especially for Jeppe's routines
*
*     integer NBSPFSM,NUSPFSM,IBSPEXP
*     real*8 CSPEXP
*     dimension CSPEXP(MXNSPIPDMGRP**2)
*     dimension NBSPFSM(MXNSPIR),NUSPFSM(MXNSPIR)
*     dimension IBSPEXP(2,2)
*     common /COEFF/CSPEXP,NBSPFSM,NUSPFSM,IBSPEXP
*
      dimension IDIM(MXNSPIR)
      dimension COEF(ITOTDIM,ITOTDIM,2)
      dimension ADJCO(ICRMXDIM,ICRMXDIM)
      dimension CMOAO(NIRR_PN*ICRMXDIM**2)
      dimension HELPMAT(ICR2DIM**2+ICR3DIM**2)
CTF
*  In VERY pathological cases, the field HELPMAT might be larger
*  than SYMDIM(IDIMGRP)**2*2, as it is dimensioned in oneltra.
*  Therefore, we check whether ICHECK has been destroyed or not:
      if (ICHECK.ne.987654)
     &    Call Abend2( 'Problem in rdcoef; HELPMAT too large' )
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
*  Either the coeffs are already on CSPAOSF, then skip the
*  following "instructions" or otherwise execute them; read
*  coeffs from file 19, construct matrix CSPAOSF.
*
*  determine total length required for CMOAO
      ILEN = 0
      do I=1,NIRR_PN,1
         ILEN = ILEN + IDIM(I)**2
      end do
C?    write(6,*) 'ILEN is    ',ILEN
*
*  Read irrepwise to array ADJCO, transfer to vector CMOAO, next irrep
      if (INTITER.eq.1.and.IOUTLST.eq.1) then
         if (IDOUGRP.eq.4) then
            IDIMREO4 = IDIM(4)
            IDIMREO3 = IDIM(3)
            IDIMREO2 = IDIM(2)
            IDIM(2) = IDIMREO3
            IDIM(3) = IDIMREO4
            IDIM(4) = IDIMREO2
         end if
         K = 1
         read(LU19,'(A)')
         do IRREP=1,NIRR_PN,1
            do I=1,IDIM(IRREP),1
               read(LU19,'(4F18.14)') (ADJCO(I,J),J=1,IDIM(IRREP),1)
            end do
            do I=1,IDIM(IRREP),1
               do J=1,IDIM(IRREP),1
                  CMOAO(K) = ADJCO(I,J)
                  K = K + 1
               end do
            end do
         end do
         if (IDOUGRP.eq.4) then
            IDIMREO4 = IDIM(4)
            IDIMREO3 = IDIM(3)
            IDIMREO2 = IDIM(2)
            IDIM(2) = IDIMREO4
            IDIM(3) = IDIMREO2
            IDIM(4) = IDIMREO3
         end if
*
         if (NTEST.ge.2) then
            do I=1,ILEN,1
               write(6,*) 'CMOAO(',I,') is  ',CMOAO(I)
            end do
         end if
*
*  Reorder, if IDOUGRP is C2H (4)
         if (IDOUGRP.eq.4) then
            IAUSTA = IDIM(1)**2 + 1
            IBUEND = IAUSTA + IDIM(3)**2 + IDIM(4)**2 - 1
            IBGSTA = IBUEND + 1
            IBGEND = IBGSTA + IDIM(2)**2 - 1
            IBGLEN = IDIM(2)**2
            if (NTEST.ge.3) then
               write(6,*) 'IAUSTA = ',IAUSTA
               write(6,*) 'IBUEND = ',IBUEND
               write(6,*) 'IBGSTA = ',IBGSTA
               write(6,*) 'IBGEND = ',IBGEND
               write(6,*) 'IBGLEN = ',IBGLEN
            end if
*
            J = 1
            do I=IBGSTA,IBGEND,1
               HELPMAT(J) = CMOAO(I)
               J = J + 1
            end do
            KLAST = ILEN
            do I=IBUEND,IAUSTA,-1
               CMOAO(KLAST) = CMOAO(I)
               KLAST = KLAST - 1
            end do
            IN = 1
            do I=IAUSTA,IAUSTA+IBGLEN-1,1
               CMOAO(I) = HELPMAT(IN)
               IN = IN + 1
            end do
*
            if (NTEST.ge.3) then
               write(6,*) 'Checking reordered MOAO coeffs in C2H:'
               do I=1,ILEN,1
                  write(6,*) 'CMOAO(',I,')= ',CMOAO(I)
               end do
            end if
         end if
*
         call cspaosf_from_moao(CSPEXP,NBSPFSM(1),NUSPFSM(1),
     &                          IBSPEXP,CMOAO,IDIM(1),IDIM(1))
*
*  The coefficients have been written to the common block element
*  CSPEXP. The fetcher will get the appropriate blocks for this
*  transformation run and arrange them on COEF.
      end if
*
*  Erase output matrix COEF completely
      do I=1,ITOTDIM,1
         do J=1,ITOTDIM,1
            COEF(J,I,1) = 0.0
            COEF(J,I,2) = 0.0
         end do
      end do
*
      if (IDOUGRP.eq.4) then
         ISPSM = IOUTLST + 4
      else if (IDOUGRP.eq.5.or.IDOUGRP.eq.6) then
         ISPSM = IOUTLST + 2
      end if
*
      do IREIM=1,2,1
         do IBU=1,2,1
            do IAB=1,2,1
               call fetch_cspexp_blk(CSPEXP,IBSPEXP,IDIM(1),IDIM(1),
     &                               IREIM,IBU,ISPSM,IAB,
     &                               CMOAO,LROW,LCOL)
               ILENBLK = LROW*LCOL
*
*  Determine whether alpha preceeds beta or vice versa
               call evenodd(IEOLST,IOUTLST)
               if (IEOLST.eq.1) IPERMAB = 1
               if (IEOLST.eq.2) IPERMAB = 2
*
*  Determine active dimensions and offsets for vectomat
               IDIMACT = 1
               if (IOUTLST.ge.3) IDIMACT = 3
               if (IDOUGRP.eq.1) IDIMACT = 1
               if (IPERMAB.eq.1.and.IAB.eq.1) then
                  IROWOFF = 1
                  if (IBU.eq.1) ICOLOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.2) ICOLOFF = 1
               end if
               if (IPERMAB.eq.1.and.IAB.eq.2) then
                  IROWOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.1) ICOLOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.2) ICOLOFF = 1
               end if
               if (IPERMAB.eq.2.and.IAB.eq.1) then
                  IROWOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.2) ICOLOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.1) ICOLOFF = 1
               end if
               if (IPERMAB.eq.2.and.IAB.eq.2) then
                  IROWOFF = 1
                  if (IBU.eq.2) ICOLOFF = IDIM(IDIMACT) + 1
                  if (IBU.eq.1) ICOLOFF = 1
               end if
*
               if (NTEST.ge.10) then
                  write(6,*) 'Testing before vectomat ...'
                  write(6,*) 'IREIM =  ',IREIM
                  write(6,*) 'IBU =  ',IBU
                  write(6,*) 'IAB =  ',IAB
                  write(6,*) 'IROWOFF =  ',IROWOFF
                  write(6,*) 'ICOLOFF =  ',ICOLOFF
               end if
*
               call vectomat(CMOAO,COEF(1,1,IREIM),ILENBLK,
     &                       ITOTDIM,ITOTDIM,
     &                       LROW,LCOL,IROWOFF,ICOLOFF)
               if (NTEST.gt.2) then
                  write(6,*) 'Checking coefficients ...'
                  call wrtmat(COEF(1,1,IREIM),ITOTDIM,ITOTDIM,
     &                        ITOTDIM,ITOTDIM)
               end if
            end do
         end do
      end do
*
      if (NTEST.gt.2) then
         write(6,*) 'Full coefficient matrices: real  imaginary'
         call wrtmat(COEF(1,1,1),ITOTDIM,ITOTDIM,ITOTDIM,ITOTDIM)
         call wrtmat(COEF(1,1,2),ITOTDIM,ITOTDIM,ITOTDIM,ITOTDIM)
      end if
*
      return
      end
