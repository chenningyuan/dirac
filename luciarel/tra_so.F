* Transform spin-orbit integrals only
*
      subroutine tra_so(SPININTS,HE_F,
     &                  AOTRISM_F,SOINPUT_F,
     &                  HELPMAT2,TRAINT,TRAINT_F,
     &                  COEF_F,ADJCO_F,PROD,
     &                  PROD_F,SOINT_F,
     &                  ISMIND,IGSOFF,IGSIND,
     &                  ITRALINE,IMXDIM,IMXDIM_F,MXGSDIM,
     &                  IOUTLST,ICURRDIM,ICURRDIM_F,NTRIAEL,
     &                  INTITER,ITRA_TYPE,
     &                  INTSPACE,ENVIRO,COORD,
     &                  NORB_F,
     &                  NTRSM_FULL,IFOCKEX)
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
      dimension ISMIND(ICURRDIM),
     &     IGSOFF(2*NGASPI),
     &     IGSIND(ICURRDIM,NIRR_PN),
     &     IAD2(3),NORB_F(NIRR_PN),
     &     NOSYM(2),NOSYM_F(2)
*
      dimension COEF_F(ICURRDIM_F,ICURRDIM_F,2),
     &     HE_F((4*IMXDIM_F)**2),
     &     ADJCO_F(ICURRDIM_F,ICURRDIM_F,2),
     &     PROD(ICURRDIM,ICURRDIM,2),
     &     PROD_F(ICURRDIM_F,ICURRDIM_F,2),
     &     AOTRISM_F(ITRIASUM_F),
     &     HELPMAT2((2*IMXDIM)**2),
     &     TRAINT(ICURRDIM,ICURRDIM,4),
     &     TRAINT_F(ICURRDIM_F,ICURRDIM_F,4),
     &     SOINT_F(ICURRDIM_F,ICURRDIM_F),
     &     SOINPUT_F(NTRSM_FULL)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      write(*,*) '****************************'
      write(*,*) '                            '
      write(*,*) ' This is tra_so in action   '
      write(*,*) ' Using external Fock matrix '
      write(*,*) '                            '
      write(*,*) '****************************'
*
*     We are transforming integrals with two general indices
*
      if (NTEST.ge.10) then
        write(6,*) 'Coefficients in tra_so after rdcoef_m'
        do K=1,2,1
          if (K.eq.1) write(6,*) 'REAL'
          if (K.eq.2) write(6,*) 'IMAGINARY'
          call wrtmat(COEF_F(1,1,K),SYMDIM_F(IDIMGRP),
     &                SYMDIM_F(IDIMGRP), SYMDIM_F(IDIMGRP),
     &                SYMDIM_F(IDIMGRP))
        end do
      end if
*
*  transpose coefficients ...
      IMRECO = 1
      call adjoint(ADJCO_F(1,1,1),COEF_F(1,1,1),
     &     SYMDIM_F(IDIMGRP),
     &     SYMDIM_F(IDIMGRP),IMRECO)
*
*  ... and imaginary parts of coefficients
      IMRECO = -1
      call adjoint(ADJCO_F(1,1,2),COEF_F(1,1,2),
     &     SYMDIM_F(IDIMGRP),
     &     SYMDIM_F(IDIMGRP),IMRECO)
*
* =====================================
*  Do spin-orbit integrals ...
* =====================================
*
      SPIORB = .TRUE.
      IHELEN = (4*IMXDIM_F)**2
      call intsel(SOINT_F,AOTRISM_F,SOINPUT_F,
     &            HE_F,IHELEN,
     &            NORB_F(1),SYMDIM_F(1),SYMDIM_F(IDIMGRP),
     &            IDIMGRP,NIRR_PN,
     &            ITRIA_F(1),ITRIASUM_F,
     &            ITRIADB_F(1),ITRIADBSUM_F,
     &            SPIORB,DOUGRP,NTRSM_FULL,ENVIRO,
     &            COORD)
*
*  Call transformation routine
*
      call sotra(SOINT_F,TRAINT_F,ADJCO_F,
     &           COEF_F,PROD_F,
     &           IOUTLST,
     &           SYMDIM_F(IDIMGRP),NORB_F(1),
     &           COORD,DOUGRP)
*
*     Write to appropriate output list
*
      if (COORD.eq.'X') then
        LISTIND = 6
        if (NTEST.ge.5) then
          do K=1,2,1
            write(6,*)
            write(6,'(A,A,A)') 'Writing SO-',COORD,' integrals:'
            if (K.eq.1) write(6,*) 'Real transformed ints'
            if (K.eq.2) write(6,*) 'Imaginary transformed ints'
            call wrtmat(TRAINT_F(1,1,K),SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),SYMDIM_F(IDIMGRP))
          end do
        end if
      else if (COORD.eq.'Y') then
        LISTIND = 14
        if (NTEST.ge.5) then
          do K=1,2,1
            write(6,*)
            write(6,'(A,A,A)') 'Writing SO-',COORD,' integrals:'
            if (K.eq.1) write(6,*) 'Real transformed ints'
            if (K.eq.2) write(6,*) 'Imaginary transformed ints'
            call wrtmat(TRAINT_F(1,1,K),SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),SYMDIM_F(IDIMGRP))
          end do
        end if
      else if (COORD.eq.'Z') then
        LISTIND = 22
        if (NTEST.ge.5) then
          do K=1,2,1
            write(6,*)
            write(6,'(A,A,A)') 'Writing SO-',COORD,' integrals:'
            if (K.eq.1) write(6,*) 'Real transformed ints'
            if (K.eq.2) write(6,*) 'Imaginary transformed ints'
            call wrtmat(TRAINT_F(1,1,K),SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),
     &                  SYMDIM_F(IDIMGRP),SYMDIM_F(IDIMGRP))
          end do
        end if
      end if
*
      do K=1,2,1
*
*  Now pick integrals over active/virtual orbitals
        if (DOUGRP.eq.7.or.DOUGRP.eq.8) then
          NOSYM_F(1) = NORB_F(IOUTLST)
          NOSYM_F(2) = NORB_F(IOUTLST)
          NOSYM(1) = ITOTDM(IOUTLST)
          NOSYM(2) = ITOTDM(IOUTLST)
          call pick_active(2,NOSYM_F,NOSYM,ADJCO_F(1,1,1),
     &                     ADJCO_F(1,1,2),TRAINT_F(1,1,K),
     &                     TRAINT(1,1,K),1,1,1,1,
     &                     SYMDIM_F(IDIMGRP),SYMDIM(IDIMGRP),
     &                     2,2,NTEST)
        else
          ISYMOFF = 1
          if (DOUGRP.eq.4.and.IOUTLST.ge.3) ISYMOFF = 3
          call pick_active(2,NORB_F(ISYMOFF),ITOTDM(ISYMOFF),
     &                     ADJCO_F(1,1,1),
     &                     ADJCO_F(1,1,2),TRAINT_F(1,1,K),
     &                     TRAINT(1,1,K),1,1,1,1,
     &                     SYMDIM_F(IDIMGRP),SYMDIM(IDIMGRP),
     &                     2,2,NTEST)
        end if
*
        if (NTEST.ge.1) then
          write(6,*) 'Spin-orbit integrals after pick_active:'
          if (K.eq.1) then
           write(6,'(A5,I3,A11,A2,A5)')
     &                'List ',IOUTLST,',coordinate ',COORD,', real:'
           call wrtmat(TRAINT(1,1,1),SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &                               SYMDIM(IDIMGRP),SYMDIM(IDIMGRP))
          else if (K.eq.2) then
           write(6,'(A5,I3,A11,A2,A5)')
     &                'List ',IOUTLST,',coordinate ',COORD,',imag.:'
           call wrtmat(TRAINT(1,1,2),SYMDIM(IDIMGRP),SYMDIM(IDIMGRP),
     &                               SYMDIM(IDIMGRP),SYMDIM(IDIMGRP))
          end if
        end if
*
        INTLIST = 2*IOUTLST + K + LISTIND
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
     &                 TRAINT(1,1,K),NGSOBPI(1,1),NGASPI,NIRR_PN,
     &                 INTLIST,SYMDIM(IDIMGRP),ITOTDM(IACSYM),
     &                 ISMIND,IGSOFF,IGSIND,IACSYM,
     &                 IOUTLST,MXGSDIM,ISPSF_MO,DOUGRP)
        else
*  Store for MCSCF run ...
          if (K.eq.1) LSM = 1
          if (K.eq.2) LSM = -1
          call fold(SPININTS(1,INTLIST),TRAINT(1,1,K),
     &              SYMDIM(IDIMGRP),LSM,ITRIADB(IOUTLST))
        end if
      end do
*
* Check transformed active integral matrix for
* hermiticity/antihermiticity:
      do K=1,2,1
        if (K.eq.1) KF = 1
        if (K.eq.2) KF = -1
        XMXDIFF = 1.0D-11
        DO I = 1,SYMDIM(IDIMGRP)
          DO J = 1, I
            if (NTEST.ge.20) then
              write(6,'(A,1F18.10)')
     &          'DIFF = ',ABS(TRAINT(I,J,K) - KF * TRAINT(J,I,K))
              write(6,'(5X,2F18.10)') TRAINT(I,J,K),TRAINT(J,I,K)
            end if
            IF(ABS(TRAINT(I,J,K) - KF * TRAINT(J,I,K)).
     &             GT.XMXDIFF) THEN
              write(6,*) 'Transformed SO ints not (anti)hermitian.'
              if (K.eq.1) then
                write(6,*) 'Real integral I,J,TRAINT(I,J) '
                write(6,*) I,J,TRAINT(I,J,1)
                write(6,*) 'does not match'
                write(6,*) 'Real integral J,I,TRAINT(J,I) '
                write(6,*) J,I,TRAINT(J,I,1)
                stop 'Problem in tra_so.'
              else if (K.eq.2) then
                write(6,*) 'Imag integral I,J,TRAINT(I,J) '
                write(6,*) I,J,TRAINT(I,J,1)
                write(6,*) 'does not match'
                write(6,*) 'Imag integral J,I,TRAINT(J,I) '
                write(6,*) J,I,-TRAINT(J,I,1)
                stop 'Problem in tra_so.'
              end if
            END IF
          END DO
        END DO
      end do
*
C     if (IOUTLST.eq.2) stop 'Timo forced me to stop in tra_so.'
*
 1000 return
      end
