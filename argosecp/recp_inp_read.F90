MODULE RECP_INP_READ
CONTAINS

SUBROUTINE RECP_INP2_READGEN(RECPIN1,ngen,ns,naords,ncons,ngcs,ncrs)
! Read general information, Y. C. Park
  IMPLICIT NONE
  INTEGER RECPIN1,ngen,ns,naords,ncons,ngcs,ncrs

  CALL RECP_INP_READMENU(RECPIN1,'%GENERAL')
  READ (RECPIN1,*) ngen, ns, naords, ncons, ngcs, ncrs

! if ( ns .gt. msu ) then
!   call bummer('change msup (two places) to ',ns,2)
! elseif ( naords .gt. kaordp ) then
!   call bummer('change kaordp (one place) to ',naords,2)
! elseif ( ncons .gt. mconsu ) then
!   call bummer('change mconsp (two places) to ',ncons,2)
! elseif ( ngcs .gt. mgcsu ) then
!   call bummer('change mgcsup (one place) to ',ngcs,2)
! endif
END SUBROUTINE RECP_INP2_READGEN

SUBROUTINE RECP_INP2_READIRREP(RECPIN1,idp,nd,ityp)
  IMPLICIT NONE
  INTEGER        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  COMMON /PARMI/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  INTEGER idp(mstu,mstu,*),nd(*),NDPT,DPT1,IPRD
  INTEGER I,J,K,RECPIN1,IOS
  CHARACTER*3 ityp(*)
  CHARACTER*80 INPLINE
  CALL RECP_INP_READMENU(RECPIN1,'%IRREP')

! Set idp = 0
  IF (nst.GT.mstu) THEN
     WRITE (*,*) 'change mstup (two places) to ',nst
     CALL QUIT('Size of idp is small')
  ENDIF
  DO I = 1, mstu
     DO J = 1, mstu
        DO K = 1, mstu  ! nst->mstup
           idp(I,J,K) = 0
        ENDDO
     ENDDO
  ENDDO

! Read number of irreps, degeneracy, and character labels.
! READ (RECPIN1,'(I3,8(I2,1X,A3))') nst,(nd(I),ityp(I),I=1,nst)
  CALL RECP_INP_READLINE(RECPIN1,INPLINE)
  READ (INPLINE,'(I3,12(I3,A3))',IOSTAT=IOS) nst,(nd(I),ityp(I),I=1,nst)
  IF (IOS.NE.0) THEN
     WRITE(*,'(A80)') INPLINE
     CALL QUIT('Error in reading line') 
  ENDIF
! Assign the totally symmetric irrep binary products.
  DO I = 1,nst
    idp(I,I,1) = 1
    idp(I,1,I) = 1
    idp(1,I,I) = 1
  ENDDO

! Read the number of product triples.
  READ(RECPIN1,'(I3)') NDPT

! read in the totally symmetric irrep product triples
  DO DPT1 = 1, NDPT
     READ (RECPIN1,'(3(2X,I1))') I, J, K
     idp(I,J,K) = 1
     idp(I,K,J) = 1
     idp(J,I,K) = 1
     idp(J,K,I) = 1
     idp(K,I,J) = 1
     idp(K,J,I) = 1
  ENDDO
END SUBROUTINE RECP_INP2_READIRREP

SUBROUTINE RECP_INP2_READAO2SO(RECPIN1,naords,ngcs,maords,la,nir,CX_TEMP)
! Read in AO to SO transformation matrices 
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_mxvalue.h"
  INTEGER naords, ngcs, maords(*),la(mrup,*),nir(*)
  INTEGER AORDS1 ! ith-AO reduction set
  INTEGER IRREP1,IRREP2,IRREP_MAX
  INTEGER SO1,SO2,AO1,AO2,CX_TEMP,RECPIN1
  INTEGER, ALLOCATABLE :: TMP_AO2SO(:)
  REAL*8  AO2SO_TEMP1, AO2SO_TEMP2
  CALL RECP_INP_READMENU(RECPIN1,'%AOTOSO')

! Read number(IRREP2) and list(la) of irrep in ith-AO reduction set 
  IRREP_MAX = 0
  DO AORDS1 = 1, naords
    READ (RECPIN1,'(20i3)') IRREP2,(la(IRREP1,AORDS1), IRREP1=1,IRREP2)
    IF (IRREP2.GT.mrup) THEN
      WRITE (*,*) 'Size of la is small'
      WRITE (*,*) 'change mrup (one place) to ',IRREP2
      CALL QUIT('Size of la is small')
    ENDIF

    IF (IRREP2.GE.IRREP_MAX)  IRREP_MAX=IRREP2
    nir(AORDS1) = IRREP2
  ENDDO

! Allocate 
  CX_TEMP = IRREP_MAX*IRREP_MAX*naords
  CX_TEMP = CX_TEMP*CX_TEMP
  ALLOCATE( IPT_CX_TEMP(CX_TEMP) )
  ALLOCATE( IPT_AO2SO(IRREP_MAX,IRREP_MAX,naords) )
  ALLOCATE( TMP_AO2SO(IRREP_MAX) )
! ngcs : Number of transformation matrices relating AOs to SOs to be read in.
!        Usually this will be the same as the number of AO reduction sets.
  DO AORDS1 = 1,ngcs
     READ (RECPIN1,'(20I3)') SO2, AO2, maords(AORDS1)
!    if ( SO2 .gt. mcsup ) then
!       call bummer('change mcsup (one place) to ',SO2,2)
!    elseif ( AO2 .gt. mctup ) then
!       call bummer('change mctup (one place) to ',AO2,2)
!    endif

     DO SO1 = 1,SO2
!       READ (RECPIN1,'(20I3)') (IPT_AO2SO(AO1,SO1,AORDS1), AO1=1,AO2)
        READ (RECPIN1,'(20I3)') (TMP_AO2SO(AO1), AO1=1,AO2)
        DO AO1 = 1, AO2
           IPT_AO2SO(AO1,SO1,AORDS1) = TMP_AO2SO(AO1)*1.0d0
        ENDDO
 
        DO AO1 = 1,AO2
           AO2SO_TEMP1 = IPT_AO2SO(AO1,SO1,AORDS1)
           AO2SO_TEMP2 = ABS(AO2SO_TEMP1)
           IF (AO2SO_TEMP2.NE.0.0d0 .AND. AO2SO_TEMP2.NE.1.0d0 )  &
              IPT_AO2SO(AO1,SO1,AORDS1) = SIGN(SQRT(AO2SO_TEMP2),AO2SO_TEMP1)
!          PRINT *,'IPT_AO2SO(',AO1,SO1,AORDS1,')',IPT_AO2SO(AO1,SO1,AORDS1)
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE( TMP_AO2SO )
END SUBROUTINE RECP_INP2_READAO2SO

SUBROUTINE RECP_INP2_READORB(RECPIN1,zet,eta,lmnp1,nrcr,ncon)
  IMPLICIT NONE
#include "inc_mxvalue.h"
  REAL*8  zet(mconup,*), eta(mrcrup,mconup,*)
  INTEGER lmnp1(*),nrcr(*),ncon(*),RECPIN1
! local variables
  INTEGER NONTYP0,NONT0,IQM0,JCO0
  INTEGER I,J,K,L,M1,M2,N,BLOCK1,EXP1,COEFF1,EXP2,COEFF2

! Find starting point to read
  CALL RECP_INP_READMENU(RECPIN1,'%BASISSET')

  BLOCK1 = 0
  READ (RECPIN1,'(BN,6X,I3)') NONTYP0 
  DO I = 1, NONTYP0
!    READ (RECPIN1,'(BN,6X,I3)') NONT0
!    DO N = 1,NONT0 
        READ (RECPIN1,'(BN,6X,I3)') IQM0
        DO J = 1,IQM0
           READ (RECPIN1,'(BN,6X,I3)') JCO0
           DO K = 1, JCO0
              BLOCK1 = BLOCK1 + 1

!             Read numbers to read in each block
!             ==================================
              READ (RECPIN1,'(3I3)') EXP2,COEFF2,lmnp1(BLOCK1)
              IF (COEFF2.EQ.0) COEFF2 = 1

!             Check maximum value
              IF (EXP2  .GT.mconup) GOTO 8110 
              IF (COEFF2.GT.mrcrup) GOTO 8120

!             Save to global variables
              nrcr(BLOCK1) = COEFF2
              ncon(BLOCK1) = EXP2

!             Read exponent and coefficient in each block
!             ===========================================
              DO EXP1 = 1,EXP2
                 READ (RECPIN1,*) zet(EXP1,BLOCK1) 
                 DO COEFF1=1,COEFF2
                    READ (RECPIN1,*) eta(COEFF1,EXP1,BLOCK1)
                 ENDDO
              ENDDO 
           ENDDO 
        ENDDO
!    ENDDO
  ENDDO

  GOTO 9000

! Error message
! -------------
  8110 CONTINUE
        WRITE(*,*) 'Number of exp. is exceed mconup',EXP2
        CALL QUIT('Increase mconup')
  8120 CONTINUE
        WRITE(*,*) 'Number of coeff. is exceed mrcrup',COEFF2
        CALL QUIT('Increase mrcrup')
  9000 CONTINUE
END SUBROUTINE RECP_INP2_READORB

SUBROUTINE RECP_INP2_READCSO(RECPIN1,lcr,lls,nkcrl,nkcru,nklsl,  &
           nklsu,ncr,ncrs,ncru,lproju,ns,mcrs,zcr,ccr)
  IMPLICIT NONE
#include "inc_mxvalue.h"
  INTEGER lcr(*),lls(*),nkcrl(6,*),nkcru(6,*),nklsl(4,*),nklsu(4,*),ncr(*)
  INTEGER ncrs,ncru,lproju,ns,mcrs(*),RECPIN1
  REAL*8  zcr(*),ccr(*)
  INTEGER I,J,K,RECPCORE
  INTEGER AREP2,SOREP2,BLOCK1,kcru,kcrl,MCRS1
! Skip if there is no RECP
  IF (ncrs.EQ.0) GOTO 9000
  CALL RECP_INP_READMENU(RECPIN1,'%RECP')

  ncru   = -2
  lproju = 0
  kcru   = 0 

! DO I = 1, ns
!    mcrs(I) = 0
! ENDDO

  DO I = 1, ncrs
!    
     READ (RECPIN1,'(4(I3))')RECPCORE,AREP2,SOREP2,MCRS1
     lproju = max(lproju, AREP2-1, SOREP2)

     IF (AREP2 .GT.6) GOTO 8110
     IF (SOREP2.GT.4) GOTO 8120
     lcr(I) = AREP2
     lls(I) = SOREP2

!    Match Geometry and RECP
!    mcrs(MCRS1) = I

!    Read AREP 
!    ---------
     DO J = 1,(AREP2+1)
!    DO J = 1,AREP2
        READ (RECPIN1,*) BLOCK1
        kcrl = kcru + 1
        kcru = kcru + BLOCK1
        nkcrl(J,I) = kcrl
        nkcru(J,I) = kcru

        DO K = kcrl,kcru
           READ (RECPIN1,*) ncr(K), zcr(K), ccr(K)
           IF (ncr(K).LT.0) GOTO 8410
           ncru = max(ncru,ncr(K))
        ENDDO
     ENDDO

!    Read SOREP 
!    ----------
     DO J = 1,SOREP2
        READ (RECPIN1,*) BLOCK1
        kcrl = kcru + 1
        kcru = kcru + BLOCK1
        nklsl(J,I) = kcrl
        nklsu(J,I) = kcru

        DO K = kcrl,kcru
           READ (RECPIN1,*) ncr(K), zcr(K), ccr(K)
           ccr(K)=ccr(K)*0.5
           IF (ncr(K).LT.-2) GOTO 8510
           ncru = max(ncru,ncr(K))
        ENDDO
     ENDDO
  ENDDO

  IF (kcru.GT.mcrup) GOTO 8900
  GOTO 9000

  8110 CONTINUE
       WRITE(*,*)'AREP2 too large ',(AREP2+1)
       CALL QUIT('Error in reading AREP')
  8120 CONTINUE
       WRITE(*,*)'SOREP2 too large ',SOREP2
       CALL QUIT('Error in reading SOREP')
  8410 CONTINUE
       WRITE(*,*)'ncr(kcr) too small ',ncr(K)
       CALL QUIT('Error in reading AREP')
  8510 CONTINUE
       WRITE(*,*)'ncr(kcr) too small ',ncr(K)
       CALL QUIT('Error in reading SOREP')
  8900 CONTINUE
       WRITE(*,*)'change mcrup (one place) to ',kcru
       CALL QUIT('Error in reading RECP')
  9000 CONTINUE
END SUBROUTINE RECP_INP2_READCSO

SUBROUTINE RECP_INP2_READGEO(RECPIN1,ns,ng,NBASIS,NUC2,ica,X,Y,Z,CHG1,NUCNAME, &
           nst,nso,nbft,mcons,mgcs,nt,ntl,ntu,la,lb,nrcr,nir,lmn1u,lmnp1,  &
           maords,ncrs,mcrs)
  IMPLICIT NONE
#include "inc_mxvalue.h"
  INTEGER RECPIN1,ns,ng,NBASIS(*),NUC2(*),ica(mcup,msup,*) 
  REAL*8  X(mcup,*),Y(mcup,*),Z(mcup,*),CHG1(*)
  CHARACTER*3 NUCNAME(*)
  INTEGER nst,nso(*),nbft,mcons(*),mgcs(*),nt(*),ntl(*),ntu(*),la(mrup,*),lb(*)
  INTEGER nrcr(*),nir(*),lmn1u,lmnp1(*),maords(*),ncrs,mcrs(*)
! local variables
  INTEGER I,J,J2,K,M1,M2,BASIS1,AO2SO1,la2
  CHARACTER*3 NUCNAME0

  CALL RECP_INP_READMENU(RECPIN1,'%GEOMETRY')

  DO I = 1, nst
     nso(i) = 0
  ENDDO

  DO I = 1, ns
     mcrs(I) = 0
  ENDDO

  J2   = 0
  nbft = 0 ! total number of symmetry orbital basis functions

  DO I = 1, ns
!    Read nuclear information
!    ========================
!                     mtype      nf        nc      chg
     READ (RECPIN1,'(A3,2I3,F14.10)') NUCNAME(I),NBASIS(I),NUC2(I),CHG1(I) 

!    Error if the number of nuclear is more than mcu
     IF (NUC2(I).GT.mcup) GOTO 8100

!    Read nuclear geometry
     DO J = 1,NUC2(I)
       READ (RECPIN1,'(A3,3F14.10)') NUCNAME0, X(J,I), Y(J,I), Z(J,I)
     ENDDO
!    Read 
     IF (NUC2(I).GT.1) THEN
        DO J = 2, ng
           READ (RECPIN1,*) (ica(K,I,J),K=1,NUC2(I) )
        ENDDO
      ENDIF

!    Match basis-set block to l-quantum number 
!    =========================================
     DO J = 1, NBASIS(I)
        READ (RECPIN1,'(2(I3))') BASIS1,AO2SO1
!       Save to global variables
        J2 = J2 + 1
        mcons(J2)= BASIS1
        mgcs(J2) = AO2SO1

!       lmn1u   : maximum value of lmnp1
        lmn1u   = max(lmn1u,lmnp1(BASIS1))
!       nt      : (upper/lower : ntu/ntl)
        nt(J2)  = (lmnp1(BASIS1)+1)*( (lmnp1(BASIS1)+1) - 1)/2
        ntu(J2) = (nt(J2)*(lmnp1(BASIS1)+2))/3
        ntl(J2) = ntu(J2)-nt(J2)+1

        DO M1 = 1,nrcr(BASIS1)
           DO M2 = 1,nir(maords(AO2SO1))
              nbft = nbft+1
              la2  = la(M2,maords(AO2SO1))
              nso(la2) = nso(la2)+1
              lb(nbft) = nso(la2)
           ENDDO
        ENDDO
     ENDDO
!    Read
     IF (ncrs.GT.0) READ (RECPIN1,'(I3)') mcrs(I)
  ENDDO

  IF (J2.GT.msfup) GOTO 8510
  IF (nbft.GT.msfrup) GOTO 8520
  GOTO 9000

! Error message
  8100 CONTINUE
       WRITE(*,*)'change mcup (one place) to ',NUC2(I)
  8510 CONTINUE
        WRITE(*,*) 'change msfup to ',J2
        CALL QUIT('Increase msfup')
  8520 CONTINUE
        WRITE(*,*) 'change msfrup (one place) to ',nbft
        CALL QUIT('Increase msfrup')
  9000 CONTINUE
END SUBROUTINE RECP_INP2_READGEO

SUBROUTINE RECP_INP_READMENU(RECPINP,MENUNAME)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: MENUNAME
  CHARACTER*80 INPLINE
  INTEGER :: IOS,RECPINP,MENUSIZE
  MENUSIZE = LEN_TRIM(MENUNAME)
  REWIND RECPINP

  1000 CONTINUE
  CALL RECP_INP_READLINE(RECPINP,INPLINE) 
  IF (INPLINE(1:MENUSIZE).NE.MENUNAME) GOTO 1000
END SUBROUTINE RECP_INP_READMENU

SUBROUTINE RECP_INP_READLINE(RECPINP,INPLINE)
  IMPLICIT NONE
  CHARACTER*80 INPLINE
  INTEGER :: IOS,RECPINP
  READ (RECPINP,'(A80)',IOSTAT=IOS) INPLINE
  IF (IOS.NE.0) THEN
     WRITE(*,*) 'Error in reading line'
     WRITE(*,*) INPLINE
     CALL QUIT('Error in reading line')
  ENDIF
END SUBROUTINE RECP_INP_READLINE

END MODULE RECP_INP_READ
