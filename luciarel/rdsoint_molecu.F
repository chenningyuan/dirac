*  Read spin-orbit integrals
*
      subroutine rdsoint_molecu(SOINXLINE,SOINYLINE,
     &                   SOINZLINE,
     &                   HELPMAT,NTRIAEL,
     &                   LU69,LU70,LU71,
     &                   LUSGI2,LUSGI3,LUSGI4)
*
      implicit real*8 (A-H,O-Z)
*
      character*8 CHARTEST
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
*
      dimension SOINZLINE(NTRIAEL),SOINXLINE(NTRIAEL),
     &          SOINYLINE(NTRIAEL)
      dimension HELPMAT(*)
      dimension IL1(4,4),IL2(10),IL3(10),IL4(10)
*
      if (ARCHI.ne.'SGO') then
         open(unit=LU69,file='SOMFINT_X',
     &        status='OLD',form='UNFORMATTED')
         read(LU69) (SOINXLINE(I),I=1,NTRIAEL,1)
         close(unit=LU69)
*
         if (IPRINT.ge.2) then
            rewind LUSGI2
            write(6,*) 'writing SO - X integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINXLINE(I)
               write(LUSGI2,*) SOINXLINE(I)
            end do
            rewind LUSGI2
         end if
*
         open(unit=LU70,file='SOMFINT_Y',
     &        status='OLD',form='UNFORMATTED')
         read(LU70) (SOINYLINE(I),I=1,NTRIAEL,1)
         close(unit=LU70)
*
         if (IPRINT.ge.2) then
            rewind LUSGI3
            write(6,*) 'writing SO - Y integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINYLINE(I)
               write(LUSGI3,*) SOINYLINE(I)
            end do
            rewind LUSGI3
         end if
*
         open(unit=LU71,file='SOMFINT_Z',
     &        status='OLD',form='UNFORMATTED')
         read(LU71) (SOINZLINE(I),I=1,NTRIAEL,1)
         close(unit=LU71)
*
         if (IPRINT.ge.2) then
            rewind LUSGI4
            write(6,*) 'writing SO - Z integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINZLINE(I)
               write(LUSGI4,*) SOINZLINE(I)
            end do
            rewind LUSGI4
         end if
      else if (ARCHI.eq.'SGO') then
         rewind LUSGI2
         rewind LUSGI3
         rewind LUSGI4
         do I=1,NTRIAEL,1
            read(LUSGI2,*) SOINXLINE(I)
            read(LUSGI3,*) SOINYLINE(I)
            read(LUSGI4,*) SOINZLINE(I)
         end do
         if (IPRINT.ge.2) then
            write(6,*) 'writing SO - X integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINXLINE(I)
            end do
            write(6,*) 'writing SO - Y integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINYLINE(I)
            end do
            write(6,*) 'writing SO - Z integrals ... '
            do I=1,NTRIAEL,1
               write(6,*) I,SOINZLINE(I)
            end do
         end if
      else
         write(6,*) 'rdsoint: Invalid architecture given. I will stop.'
         Call Abend()
      end if
*
*  If DOUGRP is C2H, then integrals have to be reordered into
*  double group ordering of symmetries, i.e.
*  (presuming that group generators Z and XY have been used)
*
*     MOLECULE      PICASSO     Order in PICASSO
*     --------      -------     ----------------
*       A_g         Gamma1+           1
*       A_u         Gamma1-           3
*       B_u         Gamma2-           4
*       B_g         Gamma2+           2
*
*  This is too complicated, as the matrix is a full matrix,
*  not just symmetry non-vanishing blocks. So we proceed
*  differently: No reordering is done here. Instead, intsel
*  picks the non-vanishing 2-symmetry blocks appropriately.
*
      return
      end
