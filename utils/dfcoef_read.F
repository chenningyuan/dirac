!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

c     program to dump data from a DFCOEF file.
c     Usage: dfcoef_read FILE
c     where FILE is the path to a DFCOEF file.
C     By Ulf Ekstrom 2009.

c     this is a cut down version of REACMO (dirac/dirgp.F)
c     -Ulf
      SUBROUTINE DUMPCMO(IUNIT)
      implicit none
      CHARACTER TEXT*74
      integer iunit,idim,i,j,k,nsym,pos1, totorbs,
     &     ftell, marker_size
      double precision toterg,eigs
      DIMENSION IDIM(3,2),eigs(100000)
C
      REWIND IUNIT
C
C     Read title line
C     ===============
C
      READ (IUNIT,END=10,ERR=20) TEXT,NSYM,
     &       ((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG
      print *,'TEXT ',TEXT
      print *,'NSYM',NSYM
      do i=1,NSYM
         print *,'In fermion symmetry ',i
         print *,'  Positronic orbitals', IDIM(1,I)
         print *,'  Electronic orbitals', IDIM(2,I)
         print *,'  AO basis function  ', IDIM(3,I)
      enddo
      print *,'TOTERG',TOTERG
      
      totorbs = 0
      do i=1,NSYM
         totorbs = totorbs + IDIM(1,I)+IDIM(2,I)
      enddo

      pos1 = ftell(iunit)
      if (pos1.ne.94+12*NSYM) then
         print *,'Sorry, unknown on-disk format, I cannot'//
     &        ' guess NZ value.'
         stop 'Unknown binary format'
      endif
c     Number of bytes per "marker" surrounding each file section
c     Typically 4 (or 8 in 64 bit formats). Will be used to guess NZ
c     for reading mo coefficients.
      marker_size = 4

c     for now we just skip down to the eigenvalues
      print *,'Skipping mo coefficients'
      READ (IUNIT,END=10,ERR=20)

      read (iunit) (eigs(i),i=1,totorbs)
      print *,'Eigenvalues:'
      k = 1
      do i=1,nsym       
         if (idim(1,i).gt.0) then
            print *,' Positronic (sym =',i,')'
            do j=1,idim(1,i)
               print *,j,eigs(k)
               k = k + 1
            enddo
         endif
         print *,' Electronic (sym =',i,')'
         do j=1,idim(2,i)
            print *,j,eigs(k)
            k = k + 1
         enddo
      enddo
      return
 10   continue
      print *,'End of file while looking for data.'
 20   continue
      print *,'Error reading DFCOEF file'
      return
      END

      PROGRAM DFCOEF_READ
      implicit none
      INTEGER   NARG,JARG,LUINP,LUPRI
      CHARACTER*80 ARGC
      PARAMETER (LUINP = 1, LUPRI = 6)
      NARG = command_argument_count()
      IF (NARG .NE. 1) THEN
         WRITE(6,'(/3A/)') 'Usage: dfcoef_read FILE'
         stop
      END IF
      CALL get_command_argument(1,ARGC)
      OPEN(LUINP,FILE=ARGC,STATUS='UNKNOWN',FORM='UNFORMATTED')
      call dumpcmo(luinp)
      END
