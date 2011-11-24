      SUBROUTINE RECP_LNK_READCP(LUINFO,NUC1,NONT0)
!     Read ECP Input 
      USE RECP_NTR
      IMPLICIT NONE
#include "priunit.h"
#include "mxcent.h"
#include "molinp.h"
#include "nuclei.h"
#include "recpval.h"
#include "argoscom.h"
!     global variables
      INTEGER LUINFO,NUC1,NONT0
      REAL*8  Q,Q1
!     local variables
      INTEGER I,J,K,IOS,RECPIN1,LUBAS,INTQ0,ISTART
      LOGICAL ECPLIB
      CHARACTER*3 ECPKEYWORD,SOREP_CHR
      CHARACTER*30 ECPNAM,ECPTMP
      CHARACTER*80 ELINE

      RECPIN1 = 51
!     Check RECP
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(BN,A3)',IOSTAT=IOS) ECPKEYWORD
      IF (IOS.NE.0) GOTO 8000

      IF (ECPKEYWORD.NE.'ECP') THEN
!        Discard the line we just read.
         NMLINE = NMLINE - 1
         GOTO 9000
      ENDIF

!     Set RECP_SET/ Open 'RECP_INPUT2' file
!     =====================================
      IF (.NOT.ALLOCATED(RECP_SET)) THEN
         ALLOCATE(RECP_SET(2))
         RECP_SET(1) = 0   ! Number of RECP 
         
!        Open file 'RECP_INPUT2'
         OPEN (UNIT=RECPIN1,FILE='RECP_INPUT2',
     &         FORM='FORMATTED',STATUS='NEW')
         WRITE(RECPIN1,'(A)') '%RECP'
      ENDIF

!     Allocate variables for RECP parameters
      CALL RECP_NTRA_READCP(MXCENT)

!     Check types of RECP read
!     ========================
      ECPLIB = .FALSE.
      IF (MLINE(NMLINE)(4:6).EQ.'LIB') THEN

!        Read From library
!        -----------------
!        CALL GETENV ('BASDIR', BASDIR)
         ECPLIB = .TRUE.
         READ(MLINE(NMLINE),'(7X,A30)') ECPNAM
         ISTART = INDEX(ECPNAM,' ')
         ECPTMP = ECPNAM(1:ISTART - 1)   ! ecp name for search
         INTQ0  = NINT(CHARGE(ECPNUC+1)) ! charge index for search
         LUBAS  = 20                     ! file unit for open
         CALL RECP_LNK_FINDEL(ECPTMP,INTQ0,LUBAS)
         CALL RECP_LNK_RDECPLINE(ELINE,LUBAS)  
         READ (ELINE,'(3X,I3,I2,I2)',IOSTAT=IOS) 
     &        RECP_CORE(NUC1), AREP_ANG(NUC1), SOREP_ANG(NUC1)
         IF (IOS.NE.0) GOTO 8000
!        CALL RECP_LNK_FINPOS(NUC1,CORE0,AREP0,SOREP0,LUBAS)
      ELSE 
!        Explicit RECP parameter
!        -----------------------
!        Read #CORE, #AREP
         AREP_ANG(NUC1) = 0
         READ (MLINE(NMLINE),'(3X,I3,I2)',IOSTAT=IOS)
     &        RECP_CORE(NUC1), AREP_ANG(NUC1)
         IF (IOS.NE.0) GOTO 8000
!        Read #SOREP
         READ (MLINE(NMLINE),'(8X,A3)',IOSTAT=IOS) SOREP_CHR
         IF (IOS.NE.0) GOTO 8000
         IF (SOREP_CHR.EQ.'   ') THEN  ! Skip SOREP (blank)
            SOREP_ANG(NUC1) = 0
         ELSE
            READ (MLINE(NMLINE),'(8X,I2)',IOSTAT=IOS)SOREP_ANG(NUC1)
            IF (IOS.NE.0) GOTO 8000
         ENDIF
      ENDIF

!     Count the number of RECP   
      RECP_SET(1) = RECP_SET(1) + 1

!     Reset charge value
!     ==================
      Q1 = CHARGE(ECPNUC+1)
!     print *,'ECPNUC',ECPNUC,'Q',Q1,'CHARGE',CHARGE(ECPNUC+1)
!     print *,'CORE',RECP_CORE(NUC1)
      DO I = 1, NONT0
         ECPNUC = ECPNUC + 1
         IF (ECPNUC.GT.MXCENT) GOTO 8300
         CHARGE(ECPNUC) = Q1 - RECP_CORE(NUC1)*1.0D0
         IF (CHARGE(ECPNUC).LT.0) GOTO 8310
!        GNUEXP(ECPNUC) = 6.66D+37
      ENDDO
      WRITE(LUPRI,'(2X,A,I3,A)') 
     &     '* This atomic center has RECP with ',RECP_CORE(NUC1),
     &     ' core electrons.'
      WRITE(LUPRI,'(2X,A,F8.3)') 
     &     '  The charge value is changed as',CHARGE(ECPNUC)

!     Write RECP information into file
!     ================================
      WRITE (RECPIN1,'(4(I3))') 
     &      RECP_CORE(NUC1),(AREP_ANG(NUC1)-1),SOREP_ANG(NUC1),NUC1

      IF (.NOT.ALLOCATED(RECP_MCRS)) THEN 
         ALLOCATE(RECP_MCRS(NONTYP))
         DO I = 1, NONTYP 
            RECP_MCRS(I) = 0
         ENDDO
      ENDIF
      RECP_MCRS(NUC1) = RECP_SET(1)
!     ---------
!     Read AREP
!     ---------
      IF (AREP_ANG(NUC1).EQ.0) GOTO 9000
      DO I = 1, AREP_ANG(NUC1)
!        Read the number of element in each block 
         IF (ECPLIB) THEN
            CALL RECP_LNK_RDECPLINE(ELINE,LUBAS)
            READ (ELINE,'(BN,I3)',IOSTAT=IOS) AREP_BLK(NUC1,I)
         ELSE
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),'(BN,I3)',IOSTAT=IOS) AREP_BLK(NUC1,I)
         ENDIF
         IF (IOS.NE.0) GOTO 8000
         WRITE(RECPIN1,*) AREP_BLK(NUC1,I)

!        Read in each block 
         DO J = 1, AREP_BLK(NUC1,I)
            IF (ECPLIB) THEN
               CALL RECP_LNK_RDECPLINE(ELINE,LUBAS)
               READ (ELINE,*,IOSTAT=IOS)
     &              AREP_R(NUC1,I,J),AREP_E(NUC1,I,J),AREP_C(NUC1,I,J)
            ELSE
               CALL RDLINE(LUINFO)
               READ (MLINE(NMLINE),*,IOSTAT=IOS)
     &              AREP_R(NUC1,I,J),AREP_E(NUC1,I,J),AREP_C(NUC1,I,J)
            ENDIF
            WRITE(RECPIN1,*)
     &           AREP_R(NUC1,I,J),AREP_E(NUC1,I,J),AREP_C(NUC1,I,J)
            IF (IOS.NE.0) GOTO 8000
         ENDDO
      ENDDO

!     ---------------------
!     Read SOREP (if exist)
!     ---------------------
      IF (SOREP_ANG(NUC1).EQ.0) GOTO 9000 
      DO I = 1, SOREP_ANG(NUC1)
!        Read in each block         
         IF (ECPLIB) THEN
            CALL RECP_LNK_RDECPLINE(ELINE,LUBAS)
            READ (ELINE,'(BN,I3)',IOSTAT=IOS) SOREP_BLK(NUC1,I)
         ELSE
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),'(BN,I3)',IOSTAT=IOS) SOREP_BLK(NUC1,I)
         ENDIF
         IF (IOS.NE.0) GOTO 8000 
         WRITE(RECPIN1,*) SOREP_BLK(NUC1,I)

!        Read in each block         
         DO J = 1, SOREP_BLK(NUC1,I)
            IF (ECPLIB) THEN
               CALL RECP_LNK_RDECPLINE(ELINE,LUBAS)
               READ (ELINE,*,IOSTAT=IOS)
     &           SOREP_R(NUC1,I,J),SOREP_E(NUC1,I,J),SOREP_C(NUC1,I,J)
            ELSE
               CALL RDLINE(LUINFO)
               READ (MLINE(NMLINE),*,IOSTAT=IOS)
     &           SOREP_R(NUC1,I,J),SOREP_E(NUC1,I,J),SOREP_C(NUC1,I,J)
            ENDIF
            WRITE(RECPIN1,*)
     &           SOREP_R(NUC1,I,J),SOREP_E(NUC1,I,J),SOREP_C(NUC1,I,J)
            IF (IOS.NE.0) GOTO 8000 
         ENDDO
      ENDDO
      GOTO 9000

! Error message
 8000 CONTINUE
      WRITE(LUPRI,*) 'Error in reading line'
      WRITE(LUPRI,*) 'Line',NMLINE,': ',MLINE(NMLINE)
      CALL QUIT('Error in reading line')
 8300 CONTINUE     
      WRITE(LUPRI,*) 'ECPNUC exceeded MXCENT'
      CALL QUIT('Error in reading line')
 8310 CONTINUE     
      WRITE(LUPRI,*) 'The nuclear charge value is less than 0'
      WRITE(LUPRI,*) 'Total nuclear charge   : ',Q1
      WRITE(LUPRI,*) 'Number of core elecron : ',RECP_CORE(NUC1)
      CALL QUIT('Wrong charge or core-electron number')
 9000 CONTINUE
      IF (ALLOCATED(RECP_SET)) THEN
         IF (NUC1.EQ.NONTYP) THEN
            WRITE(RECPIN1,'(A)') '%END RECP'
            CALL RECP_NTRD_READCP  ! deallocate
         ENDIF
      ENDIF
      END

      SUBROUTINE RECP_LNK_CHECKMOL(DOCART)
      USE RECP_NTR
      IMPLICIT NONE
#include "priunit.h"
      LOGICAL DOCART

      IF (.NOT.ALLOCATED(RECP_SET)) GOTO 9000

!     Check if Cartesian gaussians is used 
      IF (.NOT.DOCART) THEN
         WRITE(LUPRI,'(A)') ' >>> ERROR in RECP_LNK_CHECKMOL'
         WRITE(LUPRI,'(A)') ' ECP calculation is only allowed for the'
         WRITE(LUPRI,'(A)') ' Cartesian Gaussians in this version'
         CALL QUIT('Error in setting Cartesian Gaussians') 
      ENDIF

 9000 CONTINUE
      END

      SUBROUTINE RECP_LNK_IRREP(MAXREP)
      USE RECP_NTR
      IMPLICIT NONE
#include "priunit.h"
#include "argoscom.h"
#include "argoscomch.h"
#include "pgroup.h"
#include "inc_geninp.h"
      INTEGER I,J,RECPIN1,MAXREP,ND(8)
      RECPIN1 = 51

      IF (.NOT.ALLOCATED(RECP_SET)) GOTO 9000

      WRITE(RECPIN1,'(A)') '%IRREP'

!     General options
!     ===============
!     WRITE(RECPIN1,'(6I3)') (ngenp),NUCIND,blk,(ncons),blk,(nmcrs) 

!     Point group symmetry
!     ====================
      DO I = 0, MAXREP
         ND(I) = 1
      ENDDO
      WRITE(RECPIN1,1000) (MAXREP+1), (ND(I),itypag(I),I=0,MAXREP)

!     P1 X P2 = P3
!     ============
      DO I = 0, nstag-1
         DO J = 0, nstag-1
            IF (itypag(I).EQ.IRTAB(IGRPAG*10 + J)) IRREP(J+1) = (I+1)
         ENDDO
      ENDDO

      IF (ndptag.EQ.0) THEN
!        --------------------
!        C1,Ci,Cs,C2 symmetry
!        --------------------
         WRITE(RECPIN1,1001) ndptag

      ELSEIF (ndptag.EQ.1) THEN
!        -------------------
!        D2,C2v,C2h symmetry
!        -------------------
         WRITE(RECPIN1,1001) ndptag
         WRITE(RECPIN1,1002) IRREP(2),IRREP(3),IRREP(4)

      ELSEIF (ndptag.EQ.7) THEN
!        ------------
!        D2h symmetry
!        ------------
         WRITE(RECPIN1,1001) ndptag
         WRITE(RECPIN1,1002) IRREP(2),IRREP(3),IRREP(4) !B1g*B2g=B3g 2*3=4
         WRITE(RECPIN1,1002) IRREP(2),IRREP(5),IRREP(6) !B1g*Au =B1u 2*5=6
         WRITE(RECPIN1,1002) IRREP(2),IRREP(7),IRREP(8) !B1g*B2u=B3u 2*7=8
         WRITE(RECPIN1,1002) IRREP(3),IRREP(5),IRREP(7) !B2g*Au =B2u 3*5=7
         WRITE(RECPIN1,1002) IRREP(3),IRREP(6),IRREP(8) !B2g*B1u=B3u 3*6=8
         WRITE(RECPIN1,1002) IRREP(4),IRREP(5),IRREP(8) !B3g*Au =B3u 4*5=8
         WRITE(RECPIN1,1002) IRREP(4),IRREP(6),IRREP(7) !B3g*B1u=B2u 4*6=7
      ENDIF

!1000 FORMAT(I3,8(I2,1X,A3))
 1000 FORMAT(I3,12(I3,A3))
 1001 FORMAT(I3)
 1002 FORMAT(3(2X,I1))
      WRITE(RECPIN1,'(A)') '%END IRREP'
 9000 CONTINUE
      END

      SUBROUTINE RECP_LNK_RDGEO(SYMTXT,NSYMOP,KASYM,blk,
     &           IQM,JCO,KATOM,KANG)
!     Write geometry data for RECP integral
      USE RECP_NTR
!     IMPLICIT NONE
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
!     INTEGER      KATOM,KANG
      DIMENSION    IQM(KATOM),JCO(KANG,KATOM)
      CHARACTER*1  KASYM(3,3)
      CHARACTER*11 SYMTXT
      REAL*8       XYZ3(3,8,MXATOM,MXATOM)
      INTEGER      MJ2(MXATOM,MXATOM),ISEND(8,MXATOM,8,MXATOM),NFUNC
      INTEGER      I,J,K,N,ICENT,RECPIN1,NSYMOP,IGENQ(3),NGENP,blk
      RECPIN1 = 51

      IF (.NOT.ALLOCATED(RECP_SET)) GOTO 9000

!     RECP integral is capable less than D2h symmetry
      IF (SYMTXT(1:2).EQ.'  ') GOTO 8000 
      IF (SYMTXT(1:1).EQ.'D' ) GOTO 8000
      IF (SYMTXT(1:1).EQ.'C' ) GOTO 8000

      CALL RECP_LNK_PATM(ISEND,XYZ3,MJ2,NSYMOP,KASYM)
      CALL RECP_LNK_PGEN(ISEND,NSYMOP,NGENP,MJ2)

!     Write geometry data into 'RECPIN1' file
!     ----------------------------------
      WRITE(RECPIN1,'(A)') '%GEOMETRY'
      ICENT = 1
      IBLOCK = 0
!     DO I = 1,NONTYP ! atomic type
!        WRITE(RECPIN1,*) ICENT, NONT(I), CHARGE(ICENT)
!        DO N = 1,NONT(I) ! symmetry-independent center
!           WRITE(RECPIN1,*) NAMN(ICENT),(CORD(J,ICENT),J=1,3)
!           ICENT = ICENT + 1
!        ENDDO
!     ENDDO
      DO I = 1,NONTYP ! atomic type
!        WRITE(RECPIN1,*) ICENT, NONT(I), CHARGE(ICENT)
         DO N = 1,NONT(I) ! symmetry-independent center
            NFUNC = 0
            DO J = 1, IQM(I)
               NFUNC = NFUNC + JCO(J,I)
            ENDDO
            WRITE(RECPIN1,'(A3,2I3,F14.10)') NAMN(ICENT),NFUNC,
     &           (MJ2(I,N)-1),CHARGE(ICENT)
            DO J = 1, MJ2(I,N)-1
               WRITE(RECPIN1,'(A3,3F14.10)') 
     &              NAMN(ICENT),(XYZ3(K,J,N,I),K=1,3)
            ENDDO

            ICENT = ICENT + 1

!           Write interchange
            IGENQ(1)=2
            IGENQ(2)=3
            IGENQ(3)=5
            IF (MJ2(I,N).GT.2) THEN
               DO J=1,NGENP
                  WRITE (RECPIN1,'(26I3)')
     &            (ISEND(K,N,IGENQ(J),I),K=1,MJ2(I,N)-1)
               ENDDO
            ENDIF

!           Match basis-block and AO2SO-matrix
            DO J = 1, IQM(I)
               IBLOCK = IBLOCK + 1
               DO K = 1, JCO(J,I)
                  WRITE(RECPIN1,'(2(I3))') NCONS0(K,J,I),IBLOCK
               ENDDO
            ENDDO

!           Index for mcrs
            WRITE(RECPIN1,'(I3)') RECP_MCRS(I)
         ENDDO
      ENDDO
      WRITE(RECPIN1,'(A)') '%END GEOMETRY'
      IF (ALLOCATED(RECP_MCRS)) DEALLOCATE(RECP_MCRS)
      DEALLOCATE ( NCONS0 ) 

!     General options
!     ===============
      WRITE(RECPIN1,'(A)') '%GENERAL'
!     Number of RECP parameters : nmcrs
!     Sum of all JCOs           : ncons
!     WRITE(RECPIN1,'(6I3)') NGENP,NUCIND,blk,(ncons),blk,(nmcrs) 
      WRITE(RECPIN1,'(6I3)') NGENP,NUCIND,blk,
     &     RECP_SET(2),blk,RECP_SET(1)
      WRITE(RECPIN1,'(A)') '%END GENERAL'

      CLOSE(UNIT=RECPIN1)
      GOTO 9000
 8000 CONTINUE
      WRITE(*,*) 'Maximum symmetry of RECP integral is D2h' 
      CALL QUIT('Error in geometry reading')
 9000 CONTINUE
      END


      SUBROUTINE RECP_LNK_PATM(ISEND,XYZ3,MJ2,NSYMOP,KASYM)
      IMPLICIT NONE
#include "mxcent.h"
#include "nuclei.h"
      CHARACTER*1 KASYM(3,3)
      INTEGER   LGEN(3,3),NIREP,LTRANS(3,8),ICLAT,IREP1,ICENT
      REAL*8    DISTSQ,NCLOSE,sqdif,tolsq
      INTEGER   I,J,K,J2,J3,NSYMOP,  MJ2(MXATOM,MXATOM),
     &          ISEND(8,MXATOM,8,MXATOM)
      REAL*8    XYZ1(3,8,MXATOM), XYZ2(3,8,MXATOM,MXATOM),
     &          XYZ3(3,8,MXATOM,MXATOM),XYZ4(3)

      tolsq=1.0d-12

!     DO I = 1, 3
!        DO J = 1, 3
!           print *, 'KASYM',I,J,KASYM(J,I)
!        ENDDO
!     ENDDO
      
!     read rotation 
!     -------------      
      DO I = 1, NSYMOP
         DO J = 1, 3
            LGEN(J,I) = 1
         ENDDO
         DO J = 1, 3
            IF(KASYM(J,I).EQ.'X') LGEN(1,I)=-1
            IF(KASYM(J,I).EQ.'Y') LGEN(2,I)=-1
            IF(KASYM(J,I).EQ.'Z') LGEN(3,I)=-1
         ENDDO
!        print *,'LGEN',LGEN(1,i),LGEN(2,i),LGEN(3,i)
      ENDDO
 
!     generate symetry point group operations
!     ---------------------------------------
      DO j=1,3
         LTRANS(j,1)=1
         IF (NSYMOP.GE.1) THEN
            LTRANS(j,2)=LGEN(j,1)
            IF (NSYMOP.GE.2) THEN
               LTRANS(j,3)=LGEN(j,2)
               LTRANS(j,4)=LGEN(j,2)*LGEN(j,1)
               IF (NSYMOP.GE.3) THEN
                  LTRANS(j,5)=LGEN(j,3)
                  LTRANS(j,6)=LGEN(j,3)*LGEN(j,1)
                  LTRANS(j,7)=LGEN(j,3)*LGEN(j,2)
                  LTRANS(j,8)=LGEN(j,3)*LGEN(j,2)*LGEN(j,1)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
 
!     set coordinate
!     --------------
      ICENT=1
      DO I=1,NONTYP
         DO J=1,NONT(I)
            DO K=1,3
               XYZ1(K,J,I)=CORD(k,ICENT) ! CORD from DIRAC
               IF(NSYMOP.eq.0)THEN
                 XYZ3(k,1,j,i)=XYZ1(K,J,I)
                 MJ2(i,j)=2
               ENDIF
            ENDDO
            ICENT=ICENT+1
         ENDDO
      ENDDO
      IF(NSYMOP.eq.0) GOTO 999
 
!     operation on to origin coordinate
!     ---------------------------------
      NIREP = 2 ** NSYMOP
 
      do i = 1, NONTYP 
         do j = 1, NONT(I)   
            MJ2(i,j)=1       
            do j2=1,NIREP
               do k=1,3
                  XYZ2(k,j2,j,i)=LTRANS(k,j2)*XYZ1(k,j,i)
                  XYZ4(k)=XYZ2(k,j2,j,i)
               enddo
!              check 
               do J3=1,j2-1
                  DISTSQ=sqdif(XYZ4(1),XYZ2(1,J3,j,i),3)
                  if (DISTSQ.lt.tolsq) then
                     goto 50
                  endif
               enddo
               do k=1,3
                  XYZ3(k,MJ2(i,j),j,i)=XYZ2(k,j2,j,i)
               enddo
               MJ2(i,j)=MJ2(i,j)+1
 50            continue
            enddo
         enddo
      enddo
!
!  operation  
!  ---------
      do i=1,NONTYP
         do IREP1=1,NIREP
            do j=1,NONT(i)
               do j2=1,MJ2(i,j)-1  ! natom(nc)
                  do k=1,3
                     XYZ4(k)=LTRANS(k,IREP1)*XYZ3(k,j2,j,i)
!                    print *,'43',XYZ4(k),LTRANS(k,IREP1),XYZ3(k,j2,j,i)
                  enddo
                  
                  NCLOSE=10000 ! closest atom
                  do j3=1,MJ2(i,j)-1
                     DISTSQ=sqdif(XYZ4(1),XYZ3(1,j3,j,i),3)
!                    print *,'DISTSQ',DISTSQ
                     if (DISTSQ.lt.NCLOSE) then
                        NCLOSE=DISTSQ
                        ICLAT=j3    ! closest atom so far
                     endif
                  enddo
                  ISEND(j2,j,IREP1,i)=ICLAT
!                 print *,'ISEND',i,IREP1,j,j2,ISEND(j2,j,IREP1,i)
               enddo
            enddo
         enddo
      enddo
 999  continue
      END

      FUNCTION sqdif( a, b, n )
!     returns the sum of the squares of the differences between
!     the elements of the vectors a(:) and b(:)
      IMPLICIT NONE
      INTEGER n,it
      REAL*8 a(n), b(n), sqdif
      sqdif = 0d0
      DO it = 1, n
         sqdif = sqdif + (a(it) - b(it))**2
      ENDDO
      RETURN
      END

      SUBROUTINE RECP_LNK_PGEN(ISEND,NSYMOP,NGENP,MJ2)
      IMPLICIT NONE
#include "mxcent.h"
#include "nuclei.h"
      INTEGER   NGENP,IGENP(8),IGENQ(3),NSYMOP
      INTEGER   MJ2(MXATOM,MXATOM), ISEND(8,MXATOM,8,MXATOM)
      INTEGER   I,J,IGEN,JGEN,KGEN,IATOM,IATOM2

      NGENP=0
      IGENQ(1)=2
      IGENQ(2)=3
      IGENQ(3)=5
      DO JGEN=1,NSYMOP
         IGEN=IGENQ(JGEN)
         DO I=1,NUCIND         !NUCIND=nsag
            DO J=1,NONT(I)
               DO IATOM=1,MJ2(I,J)-1          !MJ2(I,J)-1 = nc(I)
                  IF (ISEND(IATOM,J,IGEN,I).NE.IATOM) GOTO 550
               ENDDO
            ENDDO
         ENDDO
         GOTO 559

 550     CONTINUE               ! this is not the identity
         DO KGEN=1,NGENP
            DO I=1,NUCIND
               DO J=1,NONT(I)
                  DO IATOM=1,MJ2(I,J)-1
                     IF (ISEND(IATOM,J,IGEN       ,I) .NE.
     &                   ISEND(IATOM,J,IGENP(KGEN),I))
     &               GOTO 551
                  ENDDO
               ENDDO
            ENDDO
            GOTO 559
 551        CONTINUE            ! this is not generator KGEN
         ENDDO

         IF (NGENP.GE.2) THEN
            DO I=1,NUCIND
               DO J=1,NONT(I)
                  DO IATOM=1,MJ2(I,J)-1
                     IATOM2 = ISEND(IATOM,J,IGENP(2),I)
                     IF (ISEND(IATOM ,J,IGEN    ,I) .NE.
     &                   ISEND(IATOM2,J,IGENP(1),I)) 
     &                  GOTO 552
                  ENDDO
               ENDDO
            ENDDO
            GOTO 559
         ENDIF

 552     CONTINUE     ! this is not a product of previous generators
         NGENP=NGENP+1
         IGENP(NGENP)=IGENQ(JGEN)
 559     CONTINUE
      ENDDO
      END

      SUBROUTINE RECP_LNK_RDORB(IQM,JCO,NUC,NRC,
     &           ALPHA,CPRIM,CPRIMU,KATOM,KANG,KBLOCK,KPRIM)
!     Write Basis-set information for RECP integral
!     From SUBROUTINE ORBPRO 
      USE RECP_NTR
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
!     INTEGER KATOM,KANG,KBLOCK,KPRIM
      DIMENSION IQM(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          ALPHA(KPRIM,KBLOCK),
     &          CPRIM(KPRIM,KPRIM,KBLOCK),CPRIMU(KPRIM,KPRIM,KBLOCK)
      INTEGER I,J,K,L,M,N,RECPIN1,KBCH,MXIQM,MXJCO
!     Skip if 'RECP_SET' is not setted
      IF (.NOT.ALLOCATED(RECP_SET)) GOTO 9000


      RECPIN1 = 51

!     For allocate NCONS0
      MXIQM = 0
      MXJCO = 0
      DO I = 1,NONTYP
         MXIQM = MAX(MXIQM,IQM(I)) 
         DO J = 1,IQM(I)
            MXJCO = MAX(MXJCO,JCO(J,I))
         ENDDO
      ENDDO
      ALLOCATE ( NCONS0(MXJCO,MXIQM,NONTYP) ) 

      WRITE(RECPIN1,'(A)') '%BASISSET'
      WRITE(RECPIN1,'(A6,I3)') 'NONTYP',NONTYP

      KBCH    = 0
      KBCH3   = 0
      RECP_SET(2) = 0
!     print *,'RECP_LNK_RDORB'
      DO I = 1,NONTYP             ! atomic type
         WRITE(RECPIN1,'(A6,I3)') 'IQM   ',IQM(I)
         DO J = 1,IQM(I)       ! angular momentum value
            WRITE(RECPIN1,'(A6,I3)') 'JCO   ',JCO(J,I)
            KBCH3= KBCH3+ 1
            DO K = 1, JCO(J,I) ! blocks of given atom and L-value
               KBCH = KBCH + 1
               NCONS0(K,J,I) = KBCH
!              ----------------------------------------------------
               WRITE(RECPIN1,'(5I3)') NUC(KBCH),NRC(KBCH),J
!              print *,NUC(KBCH), K, NRC(KBCH)
               DO L = 1, NUC(KBCH)
                  WRITE(RECPIN1,*) ALPHA(L,KBCH)
                  DO M = 1, NRC(KBCH)
                     WRITE(RECPIN1,*) CPRIMU(L,M,KBCH)
                  ENDDO
!                 print*,ALPHA(L,KBCH),
!    &                  (CPRIMU(L,M,KBCH),M=1,NRC(KBCH))
               ENDDO
!              ----------------------------------------------------
            ENDDO
         ENDDO
      ENDDO
      RECP_SET(2) = KBCH
      WRITE(RECPIN1,'(A)') '%END BASISSET'

 9000 CONTINUE
      END


      SUBROUTINE RECP_LNK_CTRAN(blk)
      USE RECP_NTR
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "pincom.h"
#include "aosotr.h"
#include "sphtrm.h"
#include "huckel.h"
#include "cbirea.h"
#include "argoscom.h"
#include "argoscomch.h"
! be careful with the variables in 'ibtfun.h' 
#include "ibtfun.h"
      INTEGER soef, blk, I0,J0,K0,L0,BLKCNT,RECPIN1,NHKTA0,irowmax0
!
!     Run over Nuclear Centres
!
      RECPIN1 = 51
      IF (.NOT.ALLOCATED(RECP_SET)) GOTO 9000

      NHKTA0   = 0
      irowmax0 = 0
      DO I0 = 1, NUCIND
         DO J0 = 1, KMAX
            NHKTA0 = MAX(NHKTA0,NHKT(J0))
            irowmax0 = MAX (irowmax0, (nucdeg(I0)*KHKT(J0)) )
         ENDDO
      ENDDO

!     print *,'alloc,CTRANARG',NUCIND,NHKTA0,irowmax0,irowmax0
      ALLOCATE(CTRANARG(NUCIND,NHKTA0,irowmax0,irowmax0))
      DO I0 = 1, NUCIND
         DO J0 = 1, NHKTA0
            DO K0 = 1, irowmax0
               DO L0 = 1, irowmax0
                  CTRANARG(I0,J0,K0,L0) = 0.0d0
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 

      blk = 0
      do 401 ICENT = 1, NUCIND   ! run over nuclear center
         nucdega = nucdeg(ICENT)
         DO 411 IA = 1, KMAX     ! run over shells
            ila = 0
            icol = 1
            MULA   = ISTBAO(IA)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            NUMCFA = NUMCF(IA)
            IC     = LCLASS(IA)
            ncenta = NCENT(IA)
            irowmax(icent,NHKTA) = nucdega * KHKTA
!YCP        IF (IA-1.GE.1) THEN ! arrays boundaries fix
            if (ncenta .eq. ICENT) then
               if ((ncenta.ne.NCENT(IA-1)) .or. 
     &             (NHKTA.ne.NHKT(IA-1))  ) then
                  blk = blk +1
                  IF (NUMCFA .EQ. 1) IBLOCK = IBLOCK + 1              
                  IOFF = NAORBD
                  DO 421 NA = 1, KHKTA      ! run over components
                     ila = 0
                     DO 431 LA = 0, MAXREP  ! run over symmetry operator
                        IVARB = IBTXOR(LA,ISYMAO(NHKTA,NA))
!
!              If orbital contributes to this representation:
!
!                (IVARB is a bitstring of basic operations with
!                1 in the positions where the basic operations
!                has a different parity for the irrep and the
!                AO(when centered in origo). These basic
!                operations with different parities must then
!                not be part of the stabilizer of the center
!                for the AO).
                        IF (IBTAND(MULA,IVARB) .EQ. 0) THEN
                           ila = ila + 1
                           JKB = 0
 
!                          loop over symmetry dependent centers
                           DO 441 KB = 0, MAXREP
                              IF (IBTAND(KB,MULA) .EQ. 0) THEN
                                 JKB              = JKB + 1
                                 irow = NA + KHKTA*(ila-1)
                                 icol = NA + KHKTA*(JKB-1)
                                 CTRANARG(icent,NHKTA,irow,icol) 
     &                            = PT(IBTAND(KB,IVARB))
                           soef = nint(CTRANARG(ICENT,NHKTA,irow,icol))
                                 ireparg(icent,NHKTA,irow) = LA + 1
!                                 write(lupri,'(a,i3,i3,i3,i3,a,i3)')
!    &                   'CTRANARG(',icent,NHKTA,irow,icol,') = ',soef
                              MAXANG(icent) = NHKTA
                              END IF
  441                      CONTINUE
                        ELSE
                        END IF
  431                CONTINUE
  421             continue
               endif
            endif
!YCP        ENDIF
  411    CONTINUE
  401 CONTINUE

!     --------------------------------------
!     Write AO to SO transformation matrices
!     --------------------------------------
      
      WRITE(RECPIN1,'(A)') '%AOTOSO'
      DO I0 = 1, NUCIND
         DO J0 = 1, MAXANG(I0)
            WRITE (RECPIN1,'(20i3)') irowmax(I0,J0),
     &      (ireparg(I0,J0,K0),K0=1,irowmax(I0,J0))
         ENDDO
      ENDDO

      BLKCNT = 0
      DO I0 = 1, NUCIND
         DO J0 = 1, MAXANG(I0)
            BLKCNT = BLKCNT + 1
            WRITE(RECPIN1,'(20i3)')
     &           irowmax(I0,J0),irowmax(I0,J0),BLKCNT
            DO K0 = 1, irowmax(I0,J0)
               WRITE(RECPIN1,'(20i3)')
     &         (nint(CTRANARG(I0,J0,K0,L0)),L0=1,irowmax(I0,J0))
            ENDDO
         ENDDO   
      ENDDO   
      WRITE(RECPIN1,'(A)') '%END AOTOSO'
      DEALLOCATE (CTRANARG)
 9000 CONTINUE
      END

      SUBROUTINE RECP_LNK_FINDEL(BASNAM,INTQ0,LUBAS)
!     From FINDEL
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basis file               *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*     EXST - Logical variable that helps inquire if     *
C*            there exists a file with that name.        *
C*********************************************************
      CHARACTER*(*) BASNAM
      CHARACTER*200 STRING !long because contains the full dir-path
      CHARACTER SIGN
      LOGICAL EXST, EXST2, USE_BASNAM
      INTEGER INTQ0
#include "inftap.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
C     The environment variable is a : separated string, e.g.:
C     /programs/Dirac/basis:/programs/Dalton/basis:.
C
C     We loop over these directories in search for the basis set and
c     element given.
C
C     Depending on the machine, the string has trailing blanks
C     (ASCII char=32) or null elements (ASCII=0).
      ISTART = 1
      IBLANK=INDEX(BASDIR,CHAR(32))
      INULL =INDEX(BASDIR,CHAR(0))
      ISTOP=MAX(IBLANK,INULL)
      IF(ISTOP.EQ.0) THEN
C       If string is empty add current directory 
        BASDIR(1:1) = '.'
        ISTOP = 1
      ELSE
        ISTOP=ISTOP-1
      ENDIF
C     
 10   CONTINUE
C     
C     Find next colon in BASDIR
C     
      ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:ISTOP), ':' )
C     
C     In there are no colons, we are the last directory in BASDIR
C     
      IF ( ICOLON .EQ. ISTART - 1 ) ICOLON = ISTOP + 1
      IEND = ICOLON - 1
      STRING = BASDIR(ISTART:IEND)//' '
C
C     Add trailing slash if not present.
C
      IEND = INDEX ( STRING(1:), ' ') - 1
!     stefan: make check with idb indicates a possible access penalty if
!     BASDIR is not set...
      if(iend.lt.1)then
        write(lupri,'(/a)') '  *** error in FINDEL: environment '//
     &        'variable BASDIR is not set'
        call quit('*** error in FINDEL: BASDIR is not set')
      end if
      IF ( STRING(IEND:IEND) .NE. '/' ) THEN
         IEND = IEND + 1
         STRING(IEND:IEND) = '/'
      END IF
      STRING = STRING(1:IEND)//BASNAM
      IEND = INDEX(STRING(1:),' ') - 1
Cluuk Activate when debugging basis set search (should not be on by
Cdefault)
CLuuk WRITE(LUPRI,'(A/,A1,A,A1)')
CLuuk&     ' Searching basis : trying file: ','"',STRING(1:IEND),'"'
C
C     Inquire if "/path/basis-set" exists.
C
      INQUIRE (FILE = STRING, EXIST = EXST)
CMI  ... on certain systems (like PC-linux,OpenMPI-gfortran) string
C     beginning with '/' is NOT the proper file name ! therefore
C     continue in inquiring:
      EXST2=EXST
      USE_BASNAM=.FALSE.
      IF (.NOT.EXST2) THEN
        INQUIRE (FILE = BASNAM, EXIST = EXST)
        IF (EXST) THEN
          WRITE(LUPRI,"(/2A/2A)")
     &     "Info: used proper ECP name : ",BASNAM,
     &     "Info: instead of           : ",STRING
          USE_BASNAM=.TRUE.
        ENDIF
      ENDIF
      IF (EXST) THEN
        IF (.NOT.USE_BASNAM) THEN
         OPEN(LUBAS, FILE=STRING)
        ELSE
         OPEN(LUBAS, FILE=BASNAM)
        ENDIF
C     
C
C        Searching the file for the element.
C
 20      CONTINUE
         READ(LUBAS,'(A88)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         READ (STRING, '(A1)',IOSTAT=IOS) SIGN
         IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
            READ (STRING, 110,IOSTAT=IOS) SIGN, NUCEL
            IF (IOS.NE.0) THEN
               CALL QUIT('FINDEL: Error in reading NUCEL')
            ENDIF
 110        FORMAT (BN, A1, I4)
            IF (INTQ0.EQ. NUCEL) THEN
               RETURN
            ELSE
               GOTO 20
            END IF
         ELSE
            GOTO 20
         END IF
      END IF
      ISTART = ICOLON + 1
      IF(ISTART.GT.ISTOP) GOTO 300
      GOTO 10
 300  CONTINUE
      LSTR = 1
      DO J = LEN(BASNAM), 1, -1
         IF ( BASNAM(J:J) .NE. ' ' ) THEN
            LSTR = J
            GOTO 301
         END IF
      END DO
 301  CONTINUE
      WRITE (LUPRI,'(/,3A)')
     +     'Basis "',BASNAM(1:LSTR),'" doesn''t exist'
      WRITE(LUPRI,'(2X,A,A)') 'ECP name, STRING=',STRING
      WRITE(LUPRI,'(2X,A,A)') 'ECP name, BASNAM=',BASNAM
      IF (BASNAM(1:6) .EQ. 'Huckel') THEN
         WRITE (LUPRI,'(A)')
     +        'NB! Basis Huckel must be available to use Huckel guess!'
      END IF
      CALL QUIT('Non-existing ECP')
C
C Error messages 
C
 2000 CONTINUE
      WRITE (LUPRI,'(/2A/A,I5)')
     &   'FINDEL: Error when reading from ECP file ',STRING(1:IEND),
     &   ' IOSTAT =',IOERR
      CALL QUIT('I/O error in FINDEL in HERBAS')
 200  CONTINUE
      WRITE (LUPRI,'(/I3,2A)') INTQ0,
     &     ' is an unsupported element for ECP ',BASNAM
      CALL QUIT('Unsupported ECP element')
      END

      SUBROUTINE RECP_LNK_FINPOS(NUC1,CORE0,AREP0,SOREP0,LUBAS)
!     From FINPOS
      IMPLICIT NONE
#include "priunit.h"
      INTEGER NUC1,CORE0,AREP0,SOREP0,LUBAS
      INTEGER NUC0,IOS
      CHARACTER SIGN
      CHARACTER*88 STRING

!     Read / Find nuclear charge index
 1000 CONTINUE
      READ (LUBAS,'(A88)',ERR=8100,END=8110) STRING
      READ (STRING,'(A1)') SIGN
      IF ((SIGN.EQ.'a') .OR. (SIGN.EQ.'A')) THEN
         READ (STRING,'(2X,I3)') NUC0
         IF (NUC0.NE.NUC1) GOTO 1000
!        Read / Skip if #/!/$ 
 1100    CONTINUE
         READ (LUBAS,'(A88)',ERR=8100,END=8110) STRING
         READ (STRING,'(A1)') SIGN
         IF ((SIGN.EQ.'#').OR.(SIGN.EQ.'!').OR.(SIGN.EQ.'$')) GOTO 1100
!        Read / Find CORE0,AREP0,SOREP0
         READ (STRING,'(3X,I3,I2,I2)',IOSTAT=IOS) CORE0,AREP0,SOREP0 
         IF (IOS.NE.0) GOTO 8100
         GOTO 9000
      ELSE
         GOTO 1000
      ENDIF
 
 8100 CONTINUE
      WRITE (LUPRI,'(/A)') 'Error in reading following line'
      WRITE (LUPRI,'(A)') STRING
      CALL QUIT('Error in reading ECP file')
 8110 CONTINUE
      WRITE (LUPRI,'(/A)') 'Error in reading file, your ECP '//
     &     'is not complete.'
      CALL QUIT('Incomplete ECP file')

 9000 CONTINUE
      END

      SUBROUTINE RECP_LNK_RDECPLINE(ELINE,LUBAS)
      IMPLICIT NONE
#include "priunit.h"
      INTEGER LUBAS
      CHARACTER SIGN
      CHARACTER*80 ELINE

 1000 CONTINUE
      READ (LUBAS,'(A80)',ERR=8100,END=8110) ELINE
      READ (ELINE,'(A1)') SIGN 
      IF ((SIGN.EQ.'#').OR.(SIGN.EQ.'!').OR.(SIGN.EQ.'$')) GOTO 1000
      GOTO 9000
!     Error message
 8100 CONTINUE
      WRITE (LUPRI,'(/A)') 'Error in reading following line'
      WRITE (LUPRI,'(A)') ELINE
      CALL QUIT('Error in reading ECP file')
 8110 CONTINUE
      WRITE (LUPRI,'(/A)') 'Error in reading file, your ECP '//
     &     'is not complete.'
      CALL QUIT('Incomplete ECP file')
 9000 CONTINUE
      END
