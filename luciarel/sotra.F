*  Transform spin-orbit integrals
*  l_z, l_x, l_y
*
      subroutine sotra(SOINT,TRAINT,ADJCO,COEF,PROD,
     &                 IOUTLST,IDIM,ITOTDM,
     &                 COORD,IDBG)
*
      implicit real*8 (A-H,O-Z)
      character*2 transa,transb
      character*1 COORD
*
#include "mxpdim.inc"
#include "input.inc"
*
      dimension SOINT(IDIM,IDIM),
     &          TRAINT(IDIM,IDIM,4),
     &          ADJCO(IDIM,IDIM,2),
     &          COEF(IDIM,IDIM,2),
     &          PROD(IDIM,IDIM,2)
*
      dimension ITOTDM(MXNSPIR)
*
      NTESTL = 00
      NTEST = MAX(NTESTL,IPRINT)
*
      if (NTEST.ge.10) then
        write(6,*) 'Checking input block in sotra:'
        call wrtmat(SOINT,IDIM,IDIM,IDIM,IDIM)
      end if
*
 300  format(4G18.12)
      ZERO = 0.0D0
*
*  Erase TRAINT
      do INUM=1,4,1
        do ICOL=1,IDIM,1
          call setvec(TRAINT(1,ICOL,INUM),ZERO,IDIM)
        end do
      end do
*
*     Write to symmetry matrix (see notes for structure)
*
*     include spin integration factor according to
*       <alpha| s(z) |alpha> = 1
*       <beta| s(z) |beta> = -1
*
*       <alpha| s(x) |beta> = 1
*       <beta| s(x) |alpha> = 1
*
*       <alpha| s(y) |beta> = -i
*       <beta| s(y) |alpha> = i
*
* This scheme also applies for quaternionic double groups.
* It is not relevant which spatial symmetries are coupled by
* the spin-orbit operator. In this transformation, the structure
* of the trafo coefficient matrices is always the same, i.e.
* it has the same U/B alpha/beta structure. So the l integral
* matrices are set up in the same way, no matter which group
* is considered. sf and SO(z) SAO integrals over the l operator
* remain in the diagonal blocks, SO(x,y) integrals off-diagonal.
* The difference in coupled spatial symmetries is already taken
* care of in the atomic l integrals, and the corresponding
* expansion coefficients are copied into the transformation
* matrix accordingly.
*
* Caution: When the spinor symmetries G4 (or G6, G8 with respect
*          to C_2v^*) are being treated, the coefficient matrices
*          are rearranged, so the sign has to be applied in the
*          upper instead of the lower blocks!
*
      IACDIM = 1
      if (IOUTLST.gt.2) IACDIM = 3
      if (IDBG.eq.7.and.IOUTLST.eq.2) IACDIM = 2
*
      call evenodd(IEOD,IOUTLST)
      if (COORD.eq.'Y'.or.COORD.eq.'Z') then
        if (IEOD.eq.1) then
          do K=1,2,1
            do I=1,IDIM,1
              do J=ITOTDM(IACDIM)+1,IDIM,1
                TRAINT(J,I,K+2) = - COEF(J,I,K)
              end do
            end do
            do I=1,IDIM,1
              do J=1,ITOTDM(IACDIM),1
                TRAINT(J,I,K+2) = COEF(J,I,K)
              end do
            end do
          end do
        else if (IEOD.eq.2) then
          do K=1,2,1
            do I=1,IDIM,1
              do J=1,ITOTDM(IACDIM),1
                TRAINT(J,I,K+2) = - COEF(J,I,K)
              end do
            end do
            do I=1,IDIM,1
              do J=ITOTDM(IACDIM)+1,IDIM,1
                TRAINT(J,I,K+2) = COEF(J,I,K)
              end do
            end do
          end do
        end if
      else if (COORD.eq.'X') then
        do K=1,2,1
          do I=1,IDIM,1
            do J=1,IDIM,1
              TRAINT(J,I,K+2) = COEF(J,I,K)
            end do
          end do
        end do
      end if
*
*     Transform
*
*     Do integrals x real part
*     Result is imaginary
*
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,SOINT(1,1),
     &           IDIM,TRAINT(1,1,3),IDIM,
     &           beta,PROD(1,1,1),IDIM)
*     Do integrals x imag.
*     Result is (-1) real!!
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,SOINT(1,1),
     &           IDIM,TRAINT(1,1,4),IDIM,
     &           beta,PROD(1,1,2),IDIM)
*     include factor (-1)
      do I=1,IDIM,1
         do J=1,IDIM,1
            PROD(J,I,2) = - PROD(J,I,2)
         end do
      end do
*
*     Transform second index
*
*     Real x real
*
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,ADJCO(1,1,1),
     &           IDIM,PROD(1,1,2),IDIM,
     &           beta,TRAINT(1,1,1),IDIM)
*     imag. x imag.
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,ADJCO(1,1,2),
     &           IDIM,PROD(1,1,1),IDIM,
     &           beta,TRAINT(1,1,2),IDIM)
*
      do I=1,IDIM,1
         do J=1,IDIM,1
            TRAINT(J,I,2) = - TRAINT(J,I,2)
         end do
      end do
*
      do I=1,IDIM,1
         do J=1,IDIM,1
            TRAINT(J,I,1) = TRAINT(J,I,1) + TRAINT(J,I,2)
         end do
      end do
*     Real x imag.
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,ADJCO(1,1,1),
     &           IDIM,PROD(1,1,1),IDIM,
     &           beta,TRAINT(1,1,3),IDIM)
*     Imag. x real
      transa = 'n'
      transb = 'n'
      alpha = 1.
      beta = ZERO
      call dgemm(transa,transb,IDIM,IDIM,
     &           IDIM,alpha,ADJCO(1,1,2),
     &           IDIM,PROD(1,1,2),IDIM,
     &           beta,TRAINT(1,1,4),IDIM)
*
      do I=1,IDIM,1
         do J=1,IDIM,1
            TRAINT(J,I,2) = TRAINT(J,I,3) + TRAINT(J,I,4)
         end do
      end do
*
      if (NTEST.ge.3) then
        if (COORD.eq.'Y') then
          write(6,*) 'SO-Y before factor:'
          call wrtmat(TRAINT(1,1,2),IDIM,IDIM,
     &                IDIM,IDIM)
        end if
      end if
*
*  Include factor of i for the SO(y) operator from spin integration
*  These are the imaginary parts, so include factor (-1)
      if (COORD.eq.'Y') then
         do I=1,IDIM,1
            do J=1,IDIM,1
               TRAINT(J,I,3) = - TRAINT(J,I,2)
            end do
         end do
*
*  now interchange real and imaginary parts
         do I=1,IDIM,1
            do J=1,IDIM,1
               TRAINT(J,I,2) = TRAINT(J,I,1)
            end do
         end do
         do I=1,IDIM,1
            do J=1,IDIM,1
               TRAINT(J,I,1) = TRAINT(J,I,3)
            end do
         end do
      end if
*
      if (NTEST.ge.1) then
         write(*,*) 'The fully transformed SO ints. read:'
         write(*,*) 'TRAINT (real)  ',COORD
         call wrtmat(TRAINT(1,1,1),IDIM,IDIM,
     &               IDIM,IDIM)
         write(*,*) 'TRAINT (imag.)  ',COORD
         call wrtmat(TRAINT(1,1,2),IDIM,IDIM,
     &               IDIM,IDIM)
      end if
      return
      end
