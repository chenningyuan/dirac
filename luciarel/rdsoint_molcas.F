*  Read spin-orbit integrals
*
      subroutine rdsoint_molcas(SOINXLINE,SOINYLINE,
     &                   SOINZLINE,NTRIAEL,
     &                   LU69,LU70,LU71,
     &                   LUSGI2,LUSGI3,LUSGI4)
*
      implicit real*8 (A-H,O-Z)
*
      character*8 CHARTEST
      character*30 INTTYPE
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
*
      parameter (ALPHA = 137.036)
*
      dimension SOINZLINE(NTRIAEL),SOINXLINE(NTRIAEL),
     &          SOINYLINE(NTRIAEL)
*
      NTESTL = 0
      if (NTESTL.eq.0) then
        write(6,*) 'No testing in rdsoint_molcas.'
      end if
      NTEST = max(NTESTL,IPRINT)
*
      if (ARCHI.eq.'SGO') then
         open(unit=LU69,file='AOPROPER_MF_SYM',
     &        status='OLD',form='UNFORMATTED')
         read(LU69)
         read(LU69) (SOINXLINE(I),I=1,NTRIAEL,1)
         read(LU69)
         read(LU69) (SOINYLINE(I),I=1,NTRIAEL,1)
         read(LU69)
         read(LU69) (SOINZLINE(I),I=1,NTRIAEL,1)
         close(unit=LU69)
*
* Change to atomic units:
         do I=1,NTRIAEL,1
           SOINXLINE(I) = SOINXLINE(I) * (1/(4*ALPHA**2))
           SOINYLINE(I) = SOINYLINE(I) * (1/(4*ALPHA**2))
           SOINZLINE(I) = SOINZLINE(I) * (1/(4*ALPHA**2))
         end do
*
         if (NTEST.ge.1) then
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
*
         write(6,*)
         write(6,*) 'SO-MF integrals read from file '
         write(6,*) 'AOPROPER_MF_SYM'
      else if (ARCHI.ne.'SGO') then
        write(6,*) 'Architecture set to  ',ARCHI
        write(6,*) 'I do not know this yet (under MOLCAS).'
        Call Abend2( 'Quitting.' )
      else
         write(6,*) 'rdsoint: Invalid architecture given. I will stop.'
         Call Abend()
      end if
*
      write(6,*)
      write(6,'(A,I8)')
     &   'Number of SO integrals on full triangle: ',NTRIAEL
*
      return
      end
