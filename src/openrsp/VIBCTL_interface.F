      subroutine VIBCTL_ifc(nc,w,ALPHA,GPRIME,THETA,DIPGRAD,
     &                      dALPHAdR,dGPRIMEdR,dTHETAdR,wrk,lwrk)
      !ajt jan09 Reduced from nw to 1 frequency (thus several calls instead)
      !          Removed negative signs on all the tensors
#include "implicit.h"
#include "mxcent.h"
#include "cbilnr.h"
#include "cbivib.h"
#include "abainf.h"
#include "moldip.h"
      dimension ALPHA(3,3),GPRIME(3,3),THETA(3,6)
      dimension DIPGRAD(3,nc),dALPHAdR(3,3,nc)
      dimension dGPRIMEdR(3,3,nc),dTHETAdR(3,6,nc),wrk(lwrk)
      NFRVAL = 1 !NFRVAL in cbilnr
      FRVAL(1) = w  !FRVAL in cbilnr
      FRVAL(2:MXFR) = 0
C Apr. 22, 2009, Bin Gao
C * change lunit to -1
      lunit = -1
      call GPOPEN(lunit,'DALTON.WLK',' ','NEW',
     &                  'UNFORMATTED',(0),.false.)
      write (lunit)
      write (lunit) !three blank records
      write (lunit)
      write (lunit) NFRVAL,FRVAL,
     &                ALPHA(:,:)  ,(0d0,i=1,3*3*(MXFR-1)),           !alpha(x,y,w)
     &            (dALPHAdR(:,:,i),(0d0,k=1,3*3*(MXFR-1)),i=1,nc),   !dalpha/dg(x,y,w,c(1:MXCOOR))
     &                             (0d0,i=1,3*3*MXFR*(MXCOOR-nc)),
     &               GPRIME(:,:)  ,(0d0,i=1,3*3*(MXFR-1)),           !gprime(x,z,w)
     &           (dGPRIMEdR(:,:,i),(0d0,k=1,3*3*(MXFR-1)),i=1,nc),   !dgprime/dg(x,z,w,c(1:MXCOOR))
     &                             (0d0,i=1,3*3*MXFR*(MXCOOR-nc)),
     &               (0d0,i=1,3*3),(0d0,i=1,3*3*(MXFR-1)),           !gprime_lon(x,z,w)
     &              ((0d0,k=1,3*3),(0d0,k=1,3*3*(MXFR-1)),i=1,nc),   !dgprime_lon/dg(x,z,w,c(1:MXCOOR))
     &                             (0d0,i=1,3*3*MXFR*(MXCOOR-nc)),
     &                THETA(:,1:2),THETA(:,4),THETA(:,2:3),THETA(:,5),
     &                THETA(:,4:6),(0d0,i=1,9*3*(MXFR-1)),           !theta(x,xy,w)
     &          (dTHETAdR(:,1:2,i),dTHETAdR(:,4,i),dTHETAdR(:,2:3,i),
     &           dTHETAdR(:,5,i),dTHETAdR(:,4:6,i),                  !dtheta/dg(x,xy,w,c(1:MXCOOR))
     &                             (0d0,k=1,9*3*(MXFR-1)),i=1,nc),
     &                             (0d0,i=1,9*3*MXFR*(MXCOOR-nc)),
     &                             (0d0,i=1,9*3*MXFR)
      call GPCLOSE(lunit,'KEEP')
      call VIBINI()
      NUMHES = .false.
      HESFIL = .true. !DALTON.HES must be there!
      NINTCM = 0
      RAMAN  = any(dALPHAdR(:,:,:)/=0)
      VROA   = (RAMAN.and.any(dGPRIMEdR(:,:,:)/=0).and.
     &                    any(dTHETAdR(:,:,:)/=0))
      DIPDER = any(DIPGRAD(:,:)/=0)
      DIPFLT(:,1:nc) = DIPGRAD
      DIPFLT(:,nc+1:MXCOOR) = 0
      DOSYM(1) = .true.
      IPRINT = 6
      call VIBCTL(wrk,(lwrk)) !(lwrk) makes a local copy
      end


      subroutine HBTinterface(nc,DIPGRAD,wrk,lwrk)
#include "implicit.h"
#include "mxcent.h"
#include "cbilnr.h"
#include "cbivib.h"
#include "abainf.h"
#include "moldip.h"
      dimension DIPGRAD(3,nc)
      dimension wrk(*)

      call VIBINI()
      NUMHES = .false.
      HESFIL = .true. !DALTON.HES must be there!
      NINTCM = 0
      RAMAN  = .false.
      DIPDER = any(DIPGRAD(:,:)/=0)
      DIPFLT(:,1:nc) = -DIPGRAD; DIPFLT(:,nc+1:)=0
C      DIPDER = .true.
      DOSYM(1) = .true.
      IPRINT = 6
      i=lwrk; call VIBCTL(wrk,i)
      end


      subroutine NUCLEI_ifc(n,Z,IS,G)
      !ajt Used in prop_contribs.f90
#include "mxcent.h"
#include "nuclei.h"
      integer n,IS(n); real*8 Z(n),G(3,n)
      if (n/=NATOMS) stop 'NUCLEI_ifc: n/=NATOMS'
      Z = CHARGE(1:n)
#ifndef PRG_DIRAC
      IS = ISOTOP(1:n)
#else
!radovan: isotop array is not there in DIRAC
!         use workaround
      IS = 1
#endif
      G = CORD(:,1:n)
      end


      subroutine DIPNUC_ifc(n,F,FG)
      !ajt Used in prop_contribs.f90
#include "mxcent.h"
#include "dipole.h"
      integer :: n
      real(8) :: F(3),FG(3,3*n)
      call DIPNUC((/0d0/),(/0d0/),0,.true.)
      F = DIPMN
      FG = DDIPN(:,1:3*n)
      end


#ifndef PRG_DIRAC
      subroutine QDRNUC_ifc(Q)
      !ajt Used in prop_contribs.f90
#include "mxcent.h"
#include "nuclei.h"
#include "quadru.h"
#include "priunit.h"
      real(8) :: Q(6)
      call NUCQDR(CORD(:,1:NUCDEP),(/0d0,0d0,0d0/),LUPRI,0)
      Q = (/QDRNUC(1,1),QDRNUC(1:2,2),QDRNUC(1:3,3)/)
      end subroutine
#endif /* ifndef PRG_DIRAC */


      subroutine GRADNN_ifc(n,G)
      !ajt Used in prop_contribs.f90
#include "mxcent.h"
#ifndef PRG_DIRAC
#include "cbinuc.h"
#endif /* ifndef PRG_DIRAC */
#include "energy.h"
#include "priunit.h"
      integer n
      real(8) G(3*n)
      IPRINT = 0; MAXDIF = 1
      call NUCREP((0d0),(0d0),(0d0))
      G(1:3*n) = GRADNN(1:3*n)
      end subroutine


      subroutine HESSNN_ifc(na,H)
C     ajt aug09 Used in prop_integrals.f90
#include "mxcent.h"
#ifndef PRG_DIRAC
#include "cbinuc.h"
#endif /* ifndef PRG_DIRAC */
#include "energy.h"
#include "priunit.h"
      integer na, i, j
      real(8) H(3*na,3*na)
#ifndef PRG_DIRAC
      real(8) HESSNN(MXCOOR,MXCOOR)
      IPRINT = 0
      MAXDIF = 2
C     second and third arg only used when IPRINT > 1
      call NUCREP(HESSNN,(0d0),(0d0))
C     ajt This might come out with halved diagonal
      do j = 1, 3*na
         do i = 1, 3*na
            H(i,j) = HESSNN(max(i,j),min(i,j))
         end do
      end do
#endif /* ifndef PRG_DIRAC */
      end subroutine


      subroutine WRITE_HOTFCHT_SPECTRA(STRING,NATOMS)
#include "hotfccom.h"
      CHARACTER STRING*(*)
      LOGICAL OptionA,ANGSTROM
      INTEGER IDUMMY,LUINFO,IPOS,IPOS2,numbertypes,atomtypes,J
      CHARACTER(80) LINE
      CHARACTER(2) NAME
      CHARACTER(1) T1,ID3
      CHARACTER(15) T15
      DOUBLE PRECISION X,Y,Z,AC !AngstromConvert
      
      IF(HOTFC)WRITE(IHOTFC,'(A)')' '
      WRITE(IHOTFC,'(A,I5)')'NumberOfAtoms =',NATOMS
      WRITE(IHOTFC,'(A,A)')'SpectraType =',STRING
      WRITE(IHOTFC,'(A)')'UpperLimit/cm^-1 = 20000'
      WRITE(IHOTFC,'(A)')'LowerLimit/cm^-1 = 0'
      WRITE(IHOTFC,'(A)')'NumberOfTemperatures = 1'
      WRITE(IHOTFC,'(A)')'Temperatures/K = 0'
      INQUIRE(file='EXCMOLECULE.INP',EXIST=OptionA) 
      IF(OptionA)THEN
         AC = 0.52917720830
         !READ EXCITED STATE GEOMETRY
         LUINFO=-1
         CALL GPOPEN (LUINFO,'EXCMOLECULE.INP','OLD',' ','FORMATTED',
     &        IDUMMY,.FALSE.)
         READ (LUINFO,'(a80)') LINE
         IF (LINE(1:5) .EQ. 'BASIS') THEN
            READ (LUINFO,'(a80)') LINE
         ENDIF      
         READ (LUINFO,'(a80)') LINE
         READ (LUINFO,'(a80)') LINE
         READ (LUINFO,'(a80)') LINE
         
         ANGSTROM=.FALSE.
         IPOS = INDEX(LINE,'Angstrom')
         IF (IPOS .NE. 0) THEN
            ANGSTROM=.TRUE.
         ENDIF
         IPOS = INDEX(LINE,'Atomtypes')
         IF (IPOS .NE. 0) THEN
            IPOS2 = INDEX(LINE(IPOS:80),'=')
            IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 10)) THEN
               CALL QUIT('Incorrect input for # atomtypes')
            ELSE
               READ (LINE((IPOS+IPOS2):80),*) Atomtypes
            ENDIF
         ELSE                   !ASSUME OLD INPUT
!            READ (LINE,'(BN,A1,I4,I3,A2,10A1,D10.2,6I5)',IOSTAT=ios) CRT,&
!            & Atomtypes,MolecularCharge,SYMTXT,((KASYM(I,J),I=1,3),J=1,3),&
!            & ID3,IntegralThreshold
            READ (LINE,'(BN,A1,I4,A15,D10.2)',IOSTAT=ios) T1
     &           Atomtypes,T15,ID3
            IF(ID3 .EQ. 'A' .OR. ID3 .EQ. 'X') THEN
               Angstrom = .TRUE.
            ELSE
               Angstrom = .FALSE.
            ENDIF
         ENDIF
         
!         IF(ANGSTROM)THEN
!            WRITE(IHOTFC,'(A)')'ReferenceStructureOfFinalState/Ang='
!         ELSE
            WRITE(IHOTFC,'(A)')' '
            WRITE(IHOTFC,'(A)')'ReferenceStructureOfFinalState/a_0='
!         ENDIF
         numbertypes=1
         DO 
            IF(numbertypes>Atomtypes) EXIT
            READ (LUINFO,'(a80)') LINE
            IPOS = INDEX(LINE,'Atoms')
            IF (IPOS .NE. 0) THEN
               IPOS2 = INDEX(LINE(IPOS:),'=')
               IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
                  CALL QUIT('Incorrect input for # of atoms')
               ELSE
                  READ (LINE((IPOS+IPOS2):),*) nAtoms
                  numbertypes=numbertypes+1
               ENDIF
               DO J=1,nAtoms
                  READ (LUINFO,*) NAME,X,Y,Z
                  IF(ANGSTROM)THEN
                     IF(HOTFC)WRITE(IHOTFC,'(3F15.10)') X/AC,Y/AC,Z/AC
                  ELSE
                     IF(HOTFC)WRITE(IHOTFC,'(3F15.10)') X,Y,Z
                  ENDIF
               ENDDO
            ELSE                !OLD INPUT
               READ (LINE,'(BN,F10.0,I5)',IOSTAT=ios) AtomicCharge,
     &              nAtoms
               IF(IOS .NE. 0) THEN
                  CALL QUIT('Error in number of atoms')
               ENDIF
               numbertypes=numbertypes+1
               DO J=1,nAtoms
                  READ (LUINFO,*) NAME,X,Y,Z
                  IF(ANGSTROM)THEN
                     IF(HOTFC)WRITE(IHOTFC,'(3F15.10)') X/AC,Y/AC,Z/AC
                  ELSE
                     IF(HOTFC)WRITE(IHOTFC,'(3F15.10)') X,Y,Z
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         call gpclose(luinfo,'KEEP')
      ENDIF

      END SUBROUTINE

      SUBROUTINE WRITE_HBT_DIPOL(X,Y,Z)
#include "hotfccom.h"
      DOUBLE PRECISION X,Y,Z,AC !AngstromConvert
      WRITE(IHOTFC,'(A)')' '
      WRITE(IHOTFC,'(A)')'ElectronicTransitionDipoleMoment ='      
!      WRITE(IHOTFC,'(3F15.10)') X, Y, Z
      WRITE(IHOTFC,'(3F15.10)') 0.D0, 0.D0, 0.D0
      END SUBROUTINE
