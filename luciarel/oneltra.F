*  This piece of code is supposed to transform any type
*  of one-electron integrals from a basis of atomic functions
*  blocked by spatial symmetry to a basis of molecular
*  spinors.
*
      subroutine oneltra(SPININTS,
     &                   AOTRISM,SOINXLINE,SOINYLINE,SOINZLINE,
     &                   AOINT,HELPMAT2,TRAINT,COEF,
     &                   ADJCO,PROD,
     &                   SOINTZ,SOINTX,SOINTY,
     &                   ISMIND,IGSOFF,IGSIND,
     &                   ITRALINE,IMXDIM,MXGSDIM,
     &                   IOUTLST,ICURRDIM,NTRIAEL,
     &                   INTITER,ITRA_TYPE,
     &                   INTSPACE,ENVIRO,IFOCKEX)
*
      implicit real*8 (A-H,O-Z)
*
*  for lblas, COORD for x,y,z spin-orbit
      character*1 transa,transb,COORD
      character*72 ENVIRO
*
      logical SPIORB
*
      dimension SPININTS(ITRALINE,32)
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
*
      dimension AOINT(ICURRDIM,ICURRDIM),
     &     COEF(ICURRDIM,ICURRDIM,2),
     &     ADJCO(ICURRDIM,ICURRDIM,2),
     &     PROD(ICURRDIM,ICURRDIM,2),
     &     AOTRISM(ITRIASUM),
     &     HELPMAT2((2*IMXDIM)**2),
     &     TRAINT(ICURRDIM,ICURRDIM,4),
     &     SOINTZ(ICURRDIM,ICURRDIM),
     &     SOINTX(ICURRDIM,ICURRDIM),
     &     SOINTY(ICURRDIM,ICURRDIM),
     &     SOINZLINE(NTRIAEL),
     &     SOINXLINE(NTRIAEL),
     &     SOINYLINE(NTRIAEL)
*
      dimension ISMIND(ICURRDIM),
     &     IGSOFF(2*NGASPI),
     &     IGSIND(ICURRDIM,NIRR_PN),
     &     IAD2(3)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      write(*,*) '****************************'
      write(*,*) '                            '
      write(*,*) ' This is oneltra in action  '
      write(*,*) '                            '
      write(*,*) '****************************'
*
*     We are transforming integrals with two general indices
*
 300  format(4G18.12)
*
*  transpose coefficients ...
      IMRECO = 1
      call adjoint(ADJCO(1,1,1),COEF(1,1,1),
     &     SYMDIM(IDIMGRP),
     &     SYMDIM(IDIMGRP),IMRECO)
*
*  ... and imaginary parts of coefficients
      IMRECO = -1
      call adjoint(ADJCO(1,1,2),COEF(1,1,2),
     &     SYMDIM(IDIMGRP),
     &     SYMDIM(IDIMGRP),IMRECO)
*
*  What type of operator?
      if (IOPER.eq.1.or.IOPER.eq.3) then
*
*  Select and arrange integrals
         SPIORB = .FALSE.
         IPRODDIM = (SYMDIM(IDIMGRP))**2
         call intsel(AOINT(1,1),AOTRISM(1),HELPMAT2(1),
     &               PROD(1,1,1),IPRODDIM,
     &               ITOTDM(1),SYMDIM(1),SYMDIM(IDIMGRP),
     &               IDIMGRP,NIRR_PN,
     &               ITRIA(1),ITRIASUM,ITRIADB(1),ITRIADBSUM,
     &               SPIORB,DOUGRP,NTRIAEL,ENVIRO,COORD)
*
 10      continue
*
         if (NTEST.ge.1) then
            write(*,*) 'Full AO integral matrix for this'
            write(*,*) 'transformation run reads:'
            write(*,*) 'AOINT'
            do I=1,SYMDIM(IDIMGRP),1
               write(*,300) (AOINT(I,J),J=1,SYMDIM(IDIMGRP))
            end do
         end if
*
*     HEAVY DUTY SECTION ...
*     Transform
*
*     Do integrals x real
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,AOINT,SYMDIM(IDIMGRP),
     &              COEF(1,1,1),SYMDIM(IDIMGRP),beta,PROD(1,1,1),
     &              SYMDIM(IDIMGRP))
*
*     Do integrals x imag.
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,AOINT,SYMDIM(IDIMGRP),
     &              COEF(1,1,2),SYMDIM(IDIMGRP),beta,PROD(1,1,2),
     &              SYMDIM(IDIMGRP))
*
         if (NTEST.gt.2) then
            write(*,*) 'The halftransformed integrals read (real):'
            write(*,*) 'PROD'
            do I=1,SYMDIM(IDIMGRP),1
               write(*,300) (PROD(I,J,1),J=1,SYMDIM(IDIMGRP))
            end do
            write(*,*) 'The halftransformed integrals read (imag):'
            write(*,*) 'PROD * i'
            do I=1,SYMDIM(IDIMGRP),1
               write(*,300) (PROD(I,J,2),J=1,SYMDIM(IDIMGRP),1)
            end do
         end if
*
*     Do second half transformation
*
*     First block: real x real
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,ADJCO(1,1,1),
     &              SYMDIM(IDIMGRP),PROD(1,1,1),SYMDIM(IDIMGRP),
     &              beta,TRAINT(1,1,1),SYMDIM(IDIMGRP))
*
*     Second block: imag. x imag.
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,ADJCO(1,1,2),
     &              SYMDIM(IDIMGRP),PROD(1,1,2),SYMDIM(IDIMGRP),
     &              beta,TRAINT(1,1,2),SYMDIM(IDIMGRP))
*
*     Include factor i*i in second real block:
*
         do I=1,SYMDIM(IDIMGRP),1
            do J=1,SYMDIM(IDIMGRP),1
               TRAINT(J,I,2) = - TRAINT(J,I,2)
            end do
         end do
*
*     Add real parts
*
         do I=1,SYMDIM(IDIMGRP),1
            do J=1,SYMDIM(IDIMGRP),1
               TRAINT(J,I,1) = TRAINT(J,I,1) + TRAINT(J,I,2)
            end do
         end do
*
*     Third block: real x imag.
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,ADJCO(1,1,1),
     &              SYMDIM(IDIMGRP),PROD(1,1,2),SYMDIM(IDIMGRP),
     &              beta,TRAINT(1,1,3),SYMDIM(IDIMGRP))
*
*     Fourth block: imag. x real
*
         transa = 'n'
         transb = 'n'
         alpha = 1.
         beta = 0.
         call dgemm(transa,transb,SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &              SYMDIM(IDIMGRP),alpha,ADJCO(1,1,2),
     &              SYMDIM(IDIMGRP),PROD(1,1,1),SYMDIM(IDIMGRP),
     &              beta,TRAINT(1,1,4),SYMDIM(IDIMGRP))
*     add
         do I=1,SYMDIM(IDIMGRP),1
            do J=1,SYMDIM(IDIMGRP),1
               TRAINT(J,I,2) = TRAINT(J,I,3) + TRAINT(J,I,4)
            end do
         end do
*
         if (NTEST.ge.1) then
            write(*,*) 'The fully transformed integrals read:'
            write(*,*) 'TRAINT (real)'
            do I=1,SYMDIM(IDIMGRP),1
               write(*,300) (TRAINT(I,J,1),J=1,SYMDIM(IDIMGRP),1)
            end do
            write(*,*) 'TRAINT * i (imag.)'
            do I=1,SYMDIM(IDIMGRP),1
               write(*,300) (TRAINT(I,J,2),J=1,SYMDIM(IDIMGRP),1)
            end do
         end if
*
*     Calculate appropriate integral output list and
*     write transformed integrals to list
*
         do K=1,2,1
            INTLIST = 2*IOUTLST - 2 + K
*
*  Fold and write to appropriate integral output list!
*
            if (NTEST.ge.3) then
             if (K.eq.1) then
                write(*,*) 'Real part of integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(*,*) 'are being passed to intgas'
             else
                write(*,*) 'Im. part of integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(*,*) 'are being passed to intgas'
             end if
            end if
*
            if (ITRA_TYPE.le.3) then
*  Distribute to GAS matrices
               IACSYM = 1
               if (IOUTLST.gt.2) IACSYM = 3
               call intgas3(SPININTS(1,INTLIST),
     &                     TRAINT(1,1,K),NGSOBPI(1,1),NGASPI,NIRR_PN,
     &                     INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &                     ISMIND,IGSOFF,IGSIND,IACSYM,
     &                     IOUTLST,MXGSDIM,ISPSF_MO,DOUGRP)
            else
*  Provide space for transformed integrals for MCSCF run
*  Store for MCSCF run ...
               if (K.eq.1) LSM = 1
               if (K.eq.2) LSM = -1
               call fold(SPININTS(1,INTLIST),TRAINT(1,1,K),
     &                   SYMDIM(IDIMGRP),LSM,ITRIADB(IOUTLST))
            end if
         end do
      end if
*
* =====================================
*  Now the spin-orbit integrals ...
* =====================================
*
      if (IOPER.eq.2.or.IOPER.eq.3) then
         SPIORB = .TRUE.
         COORD = 'Z'
         IPRODDIM = (SYMDIM(IDIMGRP))**2
         call intsel(SOINTZ,AOTRISM,SOINZLINE,
     &               PROD(1,1,1),IPRODDIM,
     &               ITOTDM(1),SYMDIM(1),SYMDIM(IDIMGRP),
     &               IDIMGRP,NIRR_PN,
     &               ITRIA(1),ITRIASUM,ITRIADB(1),ITRIADBSUM,
     &               SPIORB,DOUGRP,NTRIAEL,ENVIRO,COORD)
*
*  Call transformation routine
*
         call sotra(SOINTZ,TRAINT,ADJCO,COEF,PROD,
     &              IOUTLST,SYMDIM(IDIMGRP),
     &              ITOTDM(1),COORD,DOUGRP)
*
*     Write to appropriate output list
*
         do K=1,2,1
            INTLIST = 2*IOUTLST + K + 22
            if (NTEST.ge.3) then
             if (K.eq.1) then
                write(*,*) 'Real part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             else
                write(*,*) 'Im. part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             end if
            end if
            if (ITRA_TYPE.le.3) then
               IACSYM = 1
               if (IOUTLST.gt.2) IACSYM = 3
               call intgas3(SPININTS(1,INTLIST),
     &                     TRAINT(1,1,K),NGSOBPI(1,1),NGASPI,NIRR_PN,
     &                     INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &                     ISMIND,IGSOFF,IGSIND,IACSYM,
     &                     IOUTLST,MXGSDIM,ISPSF_MO,DOUGRP)
            else
*  Store for MCSCF run ...
               if (K.eq.1) LSM = 1
               if (K.eq.2) LSM = -1
               call fold(SPININTS(1,INTLIST),TRAINT(1,1,K),
     &                   SYMDIM(IDIMGRP),LSM,ITRIADB(IOUTLST))
            end if
         end do
*
         SPIORB = .TRUE.
         COORD = 'X'
         IPRODDIM = (SYMDIM(IDIMGRP))**2
         call intsel(SOINTX,AOTRISM,SOINXLINE,
     &               PROD(1,1,1),IPRODDIM,
     &               ITOTDM(1),SYMDIM(1),SYMDIM(IDIMGRP),
     &               IDIMGRP,NIRR_PN,
     &               ITRIA(1),ITRIASUM,ITRIADB(1),ITRIADBSUM,
     &               SPIORB,DOUGRP,NTRIAEL,ENVIRO,COORD)
*
         call sotra(SOINTX,TRAINT,ADJCO,COEF,PROD,
     &              IOUTLST,SYMDIM(IDIMGRP),
     &              ITOTDM(1),COORD,DOUGRP)
*
         do K=1,2,1
            INTLIST = 2*IOUTLST + K + 6
            if (NTEST.ge.3) then
             if (K.eq.1) then
                write(*,*) 'Real part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             else
                write(*,*) 'Im. part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             end if
            end if
            if (ITRA_TYPE.le.3) then
               IACSYM = 1
               if (IOUTLST.gt.2) IACSYM = 3
               call intgas3(SPININTS(1,INTLIST),
     &                     TRAINT(1,1,K),NGSOBPI(1,1),NGASPI,NIRR_PN,
     &                     INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &                     ISMIND,IGSOFF,IGSIND,IACSYM,
     &                     IOUTLST,MXGSDIM,ISPSF_MO,DOUGRP)
            else
*  Store for MCSCF run ...
               if (K.eq.1) LSM = 1
               if (K.eq.2) LSM = -1
               call fold(SPININTS(1,INTLIST),TRAINT(1,1,K),
     &                   SYMDIM(IDIMGRP),LSM,ITRIADB(IOUTLST))
            end if
         end do
*
         SPIORB = .TRUE.
         COORD = 'Y'
         IPRODDIM = (SYMDIM(IDIMGRP))**2
         call intsel(SOINTY,AOTRISM,SOINYLINE,
     &               PROD(1,1,1),IPRODDIM,
     &               ITOTDM(1),SYMDIM(1),SYMDIM(IDIMGRP),
     &               IDIMGRP,NIRR_PN,
     &               ITRIA(1),ITRIASUM,ITRIADB(1),ITRIADBSUM,
     &               SPIORB,DOUGRP,NTRIAEL,ENVIRO,COORD)
*
         call sotra(SOINTY,TRAINT,ADJCO,COEF,PROD,
     &              IOUTLST,SYMDIM(IDIMGRP),
     &              ITOTDM(1),COORD,DOUGRP)
*
         do K=1,2,1
            INTLIST = 2*IOUTLST + K + 14
            if (NTEST.ge.3) then
             if (K.eq.1) then
                write(*,*) 'Real part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             else
                write(*,*) 'Im. part of Integrals of list  '
                write(*,*) IOUTLST,'INTLIST is   ',INTLIST
                write(6,*) 'coordinate ',COORD
                write(*,*) 'are being passed to intgas'
             end if
            end if
            if (ITRA_TYPE.le.3) then
               IACSYM = 1
               if (IOUTLST.gt.2) IACSYM = 3
               call intgas3(SPININTS(1,INTLIST),
     &                     TRAINT(1,1,K),NGSOBPI(1,1),NGASPI,NIRR_PN,
     &                     INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &                     ISMIND,IGSOFF,IGSIND,IACSYM,
     &                     IOUTLST,MXGSDIM,ISPSF_MO,DOUGRP)
            else
*  Store for MCSCF run ...
               if (K.eq.1) LSM = 1
               if (K.eq.2) LSM = -1
               call fold(SPININTS(1,INTLIST),TRAINT(1,1,K),
     &                   SYMDIM(IDIMGRP),LSM,ITRIADB(IOUTLST))
            end if
         end do
      end if
*
      if (IOPER.le.0.or.IOPER.gt.3) then
         write(*,*) 'IOPER out of range.'
         write(*,*) 'IOPER =  ',IOPER
         Call Abend()
      end if
*
 1000 return
      end
