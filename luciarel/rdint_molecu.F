*  read integrals from direct access file fort.2
*
      subroutine rdint_molecu(AOTRISM,HELPMAT,
     &                 LU2,LU9,LUSGI1,LENGTH,IHELEN)
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "symm.inc"
#include "onel.inc"
#include "input.inc"
*
      character*8 CAR(6)
      integer*4 IDUM,IDUMMY
*
* LUCIA common block for nuclear repulsion energy
      COMMON/CECORE/ECORE
*
      dimension AOTRISM(LENGTH)
      dimension TITLE(18)
      dimension IDUMMY(MXNSPIR)
      dimension HELPMAT(IHELEN)
*
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRINT)
*
* Calculate length of integral vector
      ILEN = 0
      do I=1,NIRR_PN,1
         ILEN = ILEN + ((ITOTDM(I)**2 + ITOTDM(I)) /2)
      end do
      write(6,*) 'Number of sf integrals is    ',ILEN
      if (LENGTH.ne.ILEN) then
         write(6,*) 'Problem in rdint'
         write(6,*) 'Dimensioning incorrect'
         write(6,*) 'ILEN is  ',ILEN
         write(6,*) 'Array length is  ',LENGTH
         Call Abend2( 'I will stop' )
      end if
*
* Refers to IBM architectures at the moment ...
* Read integrals and write out for porting to SGI
* ('O' stands for the type 'octane')
      if (ARCHI.ne.'SGO') then
         IADD2 = 0
         call dafile(LU2,2,IAD2(1),3,IADD2)
         if (NTEST.ge.3) then
            write(6,*) 'IADD2 is    ',IADD2
            write(6,*) 'Checking address vector IAD2:'
            do I=1,3,1
               write(6,*) IAD2(I)
            end do
         end if
         IADD2 = IAD2(2)
*
* Initializing ...
         do I=1,LENGTH,1
            AOTRISM(I) = 0.0
         end do
*
* Read sf integrals from file ...
         call dafile(LU2,2,AOTRISM(1),2*ILEN,IADD2)
*
* Write ints to plain file for SGI run ...
         rewind LUSGI1
         if (NTEST.ge.3) then
            do I=1,ILEN,1
               write(LUSGI1,*) AOTRISM(I)
            end do
         end if
         rewind LUSGI1
      else if (ARCHI.eq.'SGO') then
         rewind LUSGI1
         do I=1,ILEN,1
            read(LUSGI1,*) AOTRISM(I)
         end do
      else
         write(6,*) 'rdint: Invalid architecture given. I will stop.'
         Call Abend()
      end if
*
      if (NTEST.ge.2) then
         do I=1,ILEN,1
            write(6,*) 'AO INT ',I,AOTRISM(I)
         end do
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
* ITOTDM are point group dimensions and given in LUCIAREL ordering!
*
      if (DOUGRP.eq.4) then
         IAUSTA = (ITOTDM(1)**2 + ITOTDM(1)) / 2 + 1
         IBUEND = IAUSTA + (ITOTDM(3)**2 + ITOTDM(3)) / 2
     &                   + (ITOTDM(4)**2 + ITOTDM(4)) / 2 - 1
         IBGSTA = IBUEND + 1
         IBGEND = IBGSTA + (ITOTDM(2)**2 + ITOTDM(2)) / 2 - 1
         if (IBGEND.ne.ILEN) then
           write(6,*) 'IBGEND is  ',IBGEND
           write(6,*) 'ILEN is  ',ILEN
           write(6,*) 'Problemo in reordering of C2h integrals.'
           Call Abend2( 'Quitting.' )
         end if
         IBGLEN = (ITOTDM(2)**2 + ITOTDM(2)) / 2
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
            HELPMAT(J) = AOTRISM(I)
            J = J + 1
         end do
         KLAST = ILEN
         do I=IBUEND,IAUSTA,-1
            AOTRISM(KLAST) = AOTRISM(I)
            KLAST = KLAST - 1
         end do
         IN = 1
         do I=IAUSTA,IAUSTA+IBGLEN-1,1
            AOTRISM(I) = HELPMAT(IN)
            IN = IN + 1
         end do
*
         if (NTEST.ge.3) then
            write(6,*) 'Checking reordered AO integrals in C2H:'
            do I=1,ILEN,1
               write(6,*) 'AOTRISM(',I,')= ',AOTRISM(I)
            end do
         end if
      end if
*
*  Nuclear repulsion energy
*  Stored on common block 'onel.inc'
      rewind LU9
      read(LU9)
      read(LU9)
      read(LU9) TITLE,CAR,IDUM,(IDUMMY(I),I=1,IDUM),POTNUC
C?    write(6,*) IDUM
C?    do I=1,IDUM,1
C?       write(6,*) IDUMMY(I)
C?    end do
      write(6,*) 'Nuclear repulsion energy:  ',POTNUC
*
* Pass nuclear repulsion energy to LUCIA ...
      ECORE = POTNUC
      return
      end
