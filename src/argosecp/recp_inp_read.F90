MODULE RECP_INP_READ
CONTAINS

SUBROUTINE RECP_INP2_READGEN(ngen,ns,naords,ncons,ngcs,ncrs)
! Read general information, Y. C. Park
  USE RECP_NTR
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER I,ngen,ns,naords,ncons,ngcs,ncrs

  CALL OUTPUT_DBG_TITLE('INP_READGEN : Input parameter')

  ngen  = RECPIN_GENVAL(1)
  ns    = RECPIN_GENVAL(2)
  naords= RECPIN_GENVAL(3)
  ngcs  = RECPIN_GENVAL(3) !naords=ngcs
  ncrs  = RECP_SET(1) ! Number of RECP parameters
  ncons = RECP_SET(2) ! Sum of all JCOs

  CALL OUTPUT_DBG_IVAL('ngen',ngen,'I3',2)
  CALL OUTPUT_DBG_IVAL('ns',ns,'I3',2)
  CALL OUTPUT_DBG_IVAL('naords',naords,'I3',2)
  CALL OUTPUT_DBG_IVAL('Number of RECP parameters (ncrs)',ncrs,'I3',2)
  CALL OUTPUT_DBG_IVAL('Number of basis block (ncons)',ncons,'I3',2)
! DEALLOCATE ( RECP_SET )

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

SUBROUTINE RECP_INP2_READIRREP(idp,nd,ityp)
  USE RECP_NTR
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER        mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  COMMON /PARMI/ mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
  INTEGER idp(mstu,mstu,*),nd(*),DPT1,IPRD
  INTEGER I,J,K,IOS
  CHARACTER*3 ityp(*)
  CHARACTER*80 INPLINE
  INTEGER DPT_TABLE(21)
  DATA DPT_TABLE / 2,3,4, 2,5,6, 2,7,8, 3,5,7, 3,6,8, 4,5,8, 4,6,7 /

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

! Set number of irreps, degeneracy, and character labels.
  nst = RECP_MAXREP(1)
  DO I = 1, RECP_MAXREP(1)
     nd(I)   = 1
     ityp(I) = RECP_ITYP(I)
  ENDDO 

  IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A)') ' * Irreps and Character tables'
  IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(I2,8A3)') nst,(ityp(I),I=1,nst)

! Assign the totally symmetric irrep binary products.
  DO I = 1,nst
    idp(I,I,1) = 1
    idp(I,1,I) = 1
    idp(1,I,I) = 1
  ENDDO

! Read the number of product triples.
  IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,I3)') ' * NDPT : ',RECP_NDPT(1)
! NDPT = RECP_NDPT(1)

! read in the totally symmetric irrep product triples
  DO DPT1 = 1, RECP_NDPT(1)
     I = RECP_IRREP( DPT_TABLE( 3*(DPT1-1)+1 ) )
     J = RECP_IRREP( DPT_TABLE( 3*(DPT1-1)+2 ) )
     K = RECP_IRREP( DPT_TABLE( 3*(DPT1-1)+3 ) )
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(I3,A,3I3)') DPT1,':',I,J,K

     idp(I,J,K) = 1
     idp(I,K,J) = 1
     idp(J,I,K) = 1
     idp(J,K,I) = 1
     idp(K,I,J) = 1
     idp(K,J,I) = 1
  ENDDO
  CALL RECP_NTRD_IRREP
END SUBROUTINE RECP_INP2_READIRREP

SUBROUTINE RECP_INP2_READAO2SO(naords,ngcs,maords,la,nir,CX0)
! Read in AO to SO transformation matrices 
  USE RECP_IPT
  USE RECP_NTR
  IMPLICIT NONE
#include "inc_mxvalue.h"
#include "inc_print.h"
  INTEGER naords, ngcs, maords(*),la(mrup,*),nir(*)
  INTEGER AORDS1 ! ith-AO reduction set
  INTEGER IRREP1,IRREP_MAX
  INTEGER SO1,AO1,AO2,CX0
  REAL*8  AO2SO_TEMP1, AO2SO_TEMP2

! Read number(nir) and list(la) of irrep in ith-AO reduction set 
! ==============================================================
  IRREP_MAX = 0
  DO AORDS1 = 1, naords
!    Save to nir,la     
     nir(AORDS1) = RECPIN_MXROW(AORDS1)
     DO IRREP1 = 1,RECPIN_MXROW(AORDS1)
        la(IRREP1,AORDS1) = RECPIN_LA1(IRREP1,AORDS1)
     ENDDO
     
     IF (RECPIN_MXROW(AORDS1).GT.mrup) THEN
        WRITE (*,*) 'Size of la is small'
        WRITE (*,*) 'change mrup (one place) to ',RECPIN_MXROW(AORDS1)
        CALL QUIT('Size of la is small')
     ENDIF

     IF (RECPIN_MXROW(AORDS1).GE.IRREP_MAX) IRREP_MAX=RECPIN_MXROW(AORDS1)
  ENDDO

! Allocate 
  CX0 = IRREP_MAX*IRREP_MAX*naords
  CX0 = CX0*CX0
  ALLOCATE( IPT_CX0(CX0) )
  ALLOCATE( IPT_AO2SO(IRREP_MAX,IRREP_MAX,naords) )

! ngcs : Number of transformation matrices relating AOs to SOs to be read in.
!        Usually this will be the same as the number of AO reduction sets.
  DO AORDS1 = 1,ngcs
     maords(AORDS1) = AORDS1
     DO SO1 = 1,RECPIN_MXROW(AORDS1)
        DO AO1 = 1,RECPIN_MXROW(AORDS1)
           IPT_AO2SO(AO1,SO1,AORDS1) = RECPIN_AO2SO(AO1,SO1,AORDS1)*1.0d0
        ENDDO
 
        DO AO1 = 1,RECPIN_MXROW(AORDS1)
           AO2SO_TEMP1 = IPT_AO2SO(AO1,SO1,AORDS1)
           AO2SO_TEMP2 = ABS(AO2SO_TEMP1)
           IF (AO2SO_TEMP2.NE.0.0d0 .AND. AO2SO_TEMP2.NE.1.0d0 )  &
              IPT_AO2SO(AO1,SO1,AORDS1) = SIGN(SQRT(AO2SO_TEMP2),AO2SO_TEMP1)
        ENDDO
     ENDDO

     IF (RECP_DBG.GE.2) THEN 
        WRITE(RECP_OUT,'(A,A,I3)') & 
           ' * AO to SO matrix saved (IPT_AO2SO)',' #',AORDS1 
        
        DO SO1 = 1,RECPIN_MXROW(AORDS1)
           WRITE(RECP_OUT,'(3X,50(F5.1))')  &
              (IPT_AO2SO(AO1,SO1,AORDS1), AO1=1,RECPIN_MXROW(AORDS1))
        ENDDO
     ENDIF
  ENDDO

  CALL RECP_NTRD_AOTOSO

END SUBROUTINE RECP_INP2_READAO2SO

SUBROUTINE RECP_INP2_READORB(zet,eta,lmnp1,nrcr,ncon)
  USE RECP_NTR
  IMPLICIT NONE
#include "inc_mxvalue.h"
#include "inc_print.h"
  REAL*8  zet(mconup,*), eta(mrcrup,mconup,*)
  INTEGER lmnp1(*),nrcr(*),ncon(*)
! local variables
  INTEGER NONTYP0,NONT0,IQM0,JCO0
  INTEGER I,J,K,L,M1,M2,N,BLOCK1,EXP1,COEFF1

  BLOCK1 = 0
  IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(/,A)') ' * Basis-set read'

  DO I = 1, RECPIN_NONTYP(1)
     DO J = 1,RECPIN_IQM(I)
        DO K = 1, RECPIN_JCO(J,I)
!          Read numbers to read in each block
!          ==================================
           BLOCK1 = BLOCK1 + 1
           IF (RECPIN_NRC(BLOCK1).EQ.0) RECPIN_NRC(BLOCK1) = 1

!          Check maximum value
           IF (RECPIN_NUC(BLOCK1).GT.mconup) GOTO 8110 
           IF (RECPIN_NRC(BLOCK1).GT.mrcrup) GOTO 8120

!          Save to global variables
           ncon(BLOCK1) = RECPIN_NUC(BLOCK1)  ! EXP2 
           nrcr(BLOCK1) = RECPIN_NRC(BLOCK1)  ! COEFF2
           lmnp1(BLOCK1)= RECPIN_LMNP(BLOCK1)
           IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,3I3)') ' ncon, nrcr, lmnp1 : ',&
              ncon(BLOCK1),nrcr(BLOCK1),lmnp1(BLOCK1)

!          Read exponent and coefficient in each block
!          ===========================================
           DO EXP1 = 1,RECPIN_NUC(BLOCK1)
              zet(EXP1,BLOCK1) = RECPIN_ALPHA(EXP1,BLOCK1)
              DO COEFF1=1,RECPIN_NRC(BLOCK1)
                 eta(COEFF1,EXP1,BLOCK1) = RECPIN_CPRIMU(EXP1,COEFF1,BLOCK1) 
              ENDDO
!             Print
              IF (RECP_DBG.GE.2) THEN
                 IF (RECPIN_NRC(BLOCK1).LE.3) THEN
                    WRITE(RECP_OUT,'(4F15.7)') zet(EXP1,BLOCK1), &
                         (eta(COEFF1,EXP1,BLOCK1),COEFF1=1,RECPIN_NRC(BLOCK1)) 
                 ELSEIF (RECPIN_NRC(BLOCK1).LE.6) THEN
                    WRITE(RECP_OUT,'(4F15.7)') zet(EXP1,BLOCK1), &
                         (eta(COEFF1,EXP1,BLOCK1),COEFF1=1,3) 
                    WRITE(RECP_OUT,'(15X,3F15.7)') &
                         (eta(COEFF1,EXP1,BLOCK1),COEFF1=4,RECPIN_NRC(BLOCK1)) 
                 ELSE
                    WRITE(RECP_OUT,'(4F15.7)') zet(EXP1,BLOCK1), &
                         (eta(COEFF1,EXP1,BLOCK1),COEFF1=1,3) 
                    WRITE(RECP_OUT,'(15X,3F15.7)') &
                         (eta(COEFF1,EXP1,BLOCK1),COEFF1=4,6) 
                    WRITE(RECP_OUT,'(A)') ' Only 6 coefficients were shown.' 
                 ENDIF 
              ENDIF
           ENDDO 

        ENDDO 
     ENDDO
  ENDDO

  CALL RECP_NTRD_BASIS
  GOTO 9000

! Error message
! -------------
  8110 CONTINUE
        WRITE(*,*) 'Number of exp. is exceed mconup',RECPIN_NUC(BLOCK1)
        CALL QUIT('Increase mconup')
  8120 CONTINUE
        WRITE(*,*) 'Number of coeff. is exceed mrcrup',RECPIN_NRC(BLOCK1)
        CALL QUIT('Increase mrcrup')
  9000 CONTINUE
END SUBROUTINE RECP_INP2_READORB


SUBROUTINE RECP_INP2_READCSO(lcr,lls,nkcrl,nkcru,nklsl,  &
           nklsu,ncr,ncrs,ncru,lproju,ns,zcr,ccr)
  USE RECP_NTR
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "mxcent.h"
#include "inc_mxvalue.h"
#include "inc_print.h"
#include "recpval.h"
#include "argoscom.h"
  INTEGER lcr(*),lls(*),nkcrl(6,*),nkcru(6,*),nklsl(4,*),nklsu(4,*),ncr(*)
  INTEGER ncrs,ncru,lproju,ns
  REAL*8  zcr(*),ccr(*)
  INTEGER I1,I,J,K,INDEX1,RECPCORE
  INTEGER AREP2,SOREP2,kcru,kcrl

! Skip if there is no RECP
  IF (ncrs.EQ.0) GOTO 9000

  CALL OUTPUT_DBG_TITLE('INP_READCSO : RECP parameter')
  ncru   = -2
  lproju = 0
  kcru   = 0 

  I1 = 1
  DO I = 1, ns

     RECPCORE = RECP_CORE(I)
     AREP2    = (AREP_ANG(I)-1)
     SOREP2   = SOREP_ANG(I)

!    'ns' is sum of NONT(I)
     IF (ECPNONT(I).EQ.0) GOTO 1000

!    Print information
     IF (RECP_DBG.GE.2) THEN
        WRITE (RECP_OUT,'(/,3X,5(A,I3,X))') 'No.',I, 'NONT',ECPNONT(I), & 
           'Core',RECPCORE,'AREP block',(AREP2+1),'SOREP block',SOREP2
     ENDIF

!    If Atom doesn't have the RECP parameter
     IF (AREP_ANG(I).EQ.0) THEN
        IF (RECP_DBG.GE.2) & 
           WRITE (RECP_OUT,'(3X,A)') 'This nuclear center has no core potential'
        GOTO 1000 
     ENDIF

     lproju = max(lproju, (AREP2-1), SOREP2)

     IF (AREP2 .GT.6) GOTO 8110
     IF (SOREP2.GT.4) GOTO 8120
     lcr(I1) = AREP2
     lls(I1) = SOREP2
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,I3,A,I3)')  &
        ' * Number of AREP BLOCK-1 : lcr(',I1,') = ', lcr(I1)
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,I3,A,I3)')  &
        ' * Number of SOREP BLOCK : lls(',I1,') = ', lls(I1)

!    ---------
!    Read AREP 
!    ---------
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(/,A)') ' * AREP parts'
     DO J = 1,(AREP2+1)
        IF (RECP_DBG.GE.2) WRITE(RECP_OUT,*) AREP_BLK(I,J) 
        kcrl = kcru + 1
        kcru = kcru + AREP_BLK(I,J)
        nkcrl(J,I1) = kcrl
        nkcru(J,I1) = kcru

        DO K = kcrl,kcru
           ncr(K) = RECPIN_R(K)
           zcr(K) = RECPIN_E(K)
           ccr(K) = RECPIN_C(K)
           IF (RECP_DBG.GE.2) &
              WRITE(RECP_OUT,'(I3,A,I2,2F14.7)') K,':',ncr(K),zcr(K),ccr(K)
           IF (ncr(K).LT.0) GOTO 8410
           ncru = max(ncru,ncr(K))

        ENDDO
     ENDDO

!    ----------
!    Read SOREP 
!    ----------
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(/,A)') ' * SO parts'
     DO J = 1,SOREP2
        IF (RECP_DBG.GE.2) WRITE(RECP_OUT,*) SOREP_BLK(I,J)
        kcrl = kcru + 1
        kcru = kcru + SOREP_BLK(I,J) 
        nklsl(J,I1) = kcrl
        nklsu(J,I1) = kcru

        DO K = kcrl,kcru
           ncr(K) = RECPIN_R(K)
           zcr(K) = RECPIN_E(K)
           ccr(K) = RECPIN_C(K)
           IF (RECP_DBG.GE.2) &
              WRITE(RECP_OUT,'(I3,A,I2,2F14.7)') K,':',ncr(K),zcr(K),ccr(K)
           ccr(K)=ccr(K)*0.5
           IF (ncr(K).LT.-2) GOTO 8510
           ncru = max(ncru,ncr(K))
        ENDDO
     ENDDO
     I1 = I1 + 1

     1000 CONTINUE

  ENDDO

  CALL RECP_NTRD_READCP  ! deallocate
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

SUBROUTINE RECP_INP2_READGEO(ns,ng,mcu,msu,NBASIS,NUC2,ica,X,Y,Z,CHG1,NUCNAME, &
           nst,nso,nbft,mcons,mgcs,nt,ntl,ntu,la,lb,nrcr,nir,lmn1u,lmnp1,  &
           maords,ncrs,mcrs)
  USE RECP_NTR
  USE RECP_FUNCTION
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "inc_mxvalue.h"
#include "inc_print.h"
  INTEGER ns,ng,mcu,msu,NBASIS(*),NUC2(*),ica(mcup,msup,*) 
  REAL*8  X(mcup,*),Y(mcup,*),Z(mcup,*),CHG1(*)
  CHARACTER*3 NUCNAME(*)
  INTEGER nst,nso(*),nbft,mcons(*),mgcs(*),nt(*),ntl(*),ntu(*),la(mrup,*),lb(*)
  INTEGER nrcr(*),nir(*),lmn1u,lmnp1(*),maords(*),ncrs,mcrs(*)
! local variables
  INTEGER I,J,J2,K,M1,M2,BASIS1,la2

  CALL OUTPUT_LOC('RECP_INP2_READGEO','E')
  CALL RECP_SETZERO_I1(nso,nst)
  CALL RECP_SETZERO_I1(mcrs,ns)

  J2   = 0
  nbft = 0 ! total number of symmetry orbital basis functions
  lmn1u= 0
  CALL OUTPUT_DBG_TITLE('INP_READGEO : nuclear information')
  DO I = 1, ns
!    Read nuclear information
!    ========================
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(/,A,I2,A,I2,A)') &
        ' * Read nuclear information (ns :',I,' of ',ns,')'

     NUCNAME(I) = RECPIN_NAMN(I)     ! mtype
     NBASIS(I)  = RECPIN_BASISBLK(I) ! nf
     NUC2(I)    = RECPIN_NUC2(I)     ! nc
     CHG1(I)    = RECPIN_CHARGE(I)   ! chg

     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,A4)') &
        ' * Atom name : ',NUCNAME(I)
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,I3)') &
        ' * Number of basis-set for this atom : ',NBASIS(I)
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,I3)') &
        ' * Number of atoms for this type (duplication included) : ',NUC2(I)
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,F15.7)') &
        ' * Actual Charge of this atom : ',CHG1(I)

!    Error if the number of nuclear is more than mcu
     IF (NUC2(I).GT.mcup) GOTO 8100

!    Read nuclear geometry
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A)') &
        ' * Read Nuclear Geometry'
     DO J = 1,NUC2(I)
!      NUCNAME0 =  
       X(J,I)   = RECPIN_GEO(1,J,I)
       Y(J,I)   = RECPIN_GEO(2,J,I)
       Z(J,I)   = RECPIN_GEO(3,J,I)
       IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(3X,A,I3,A,3F15.7)') &
          'No.',J,' Cord: ',X(J,I),Y(J,I),Z(J,I)
     ENDDO
!    Read interchange 
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A)') &
        ' * Read duplication from nuclear symmetry (nuc1, nuc2, ...)'
     IF (NUC2(I).GT.1) THEN
        DO J = 2, ng
           DO K = 1, NUC2(I)
!             ica(K,I,J) = RECPIN_ICA(K,I,J)  FIXED
              ica(K,I,J) = RECPIN_ICA(K,I,J-1)
           ENDDO
           IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(3X,10I3)') (ica(K,I,J),K=1,NUC2(I) )
        ENDDO
     ENDIF

!    Match basis-set block to l-quantum number 
!    =========================================
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A)') &
        ' * Match basis function and Angular Momentum number'

     DO J = 1, NBASIS(I)
!       Save to global variables
        J2 = J2 + 1
        BASIS1   = RECPIN_BLKMATCH(1,J2)
        mcons(J2)= BASIS1                !Basis block
        mgcs(J2) = RECPIN_BLKMATCH(2,J2) !L-value(=J2)
        IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(3X,3(A,I3))') &
           'J2',J2,' Basis block ',mcons(J2),' L-value',mgcs(J2) 
            

!       lmn1u   : maximum value of lmnp1
        lmn1u   = max(lmn1u,lmnp1(BASIS1))
!       nt      : (upper/lower : ntu/ntl)
        nt(J2)  = (lmnp1(BASIS1)+1)*( (lmnp1(BASIS1)+1) - 1)/2
        ntu(J2) = (nt(J2)*(lmnp1(BASIS1)+2))/3
        ntl(J2) = ntu(J2)-nt(J2)+1

        DO M1 = 1,nrcr(BASIS1)
           DO M2 = 1,nir( maords(mgcs(J2)) )
              nbft = nbft+1
              la2  = la( M2,maords(mgcs(J2)) )
              nso(la2) = nso(la2)+1
              lb(nbft) = nso(la2)
           ENDDO
        ENDDO
     ENDDO
!    Read ECP number
     IF (ncrs.GT.0) mcrs(I) = RECP_MCRS(I)
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,I3)') & 
        ' * RECP No. to read :',RECP_MCRS(I)
  ENDDO


  IF (ALLOCATED(RECP_MCRS)) DEALLOCATE(RECP_MCRS)
  CALL RECP_NTRD_GEO

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



SUBROUTINE RECP_INP_MENU
! Read menu file
  USE RECP_OUTPUT
  IMPLICIT NONE
#include "priunit.h"
#include "inc_print.h"
  INTEGER I,IOS,NTABLE
  PARAMETER (NTABLE=2)
  CHARACTER WORD0*7,WORD1*7,TABLE(NTABLE)*7
  DATA WORD0 /'*ECP   '/
  DATA TABLE /'.PRINT ','.DEBUG '/

! Original parameter
  PRINT_LEVEL = 1
  RECP_DBG    = 0
  RECP_OUT    = 6 ! DIRAC output
  
! Read MENU file (*ECP)
! ---------------------
  REWIND (LUCMD,IOSTAT=IOS)
  1000 CONTINUE
  READ (LUCMD,'(A7)',END=9900,ERR=9000) WORD1
  IF (WORD1.EQ.WORD0) THEN
     GOTO 2000
  ELSE
     GOTO 1000
  ENDIF
  CALL OUTPUT_DBG_CHAR1('No *ECP is detected in MENU file',1) 
  GOTO 9900 

! Check options in *ECP
! ---------------------
  2000 CONTINUE
  READ (LUCMD,'(A7)',END=9900,ERR=9000) WORD1 
  IF (WORD1(1:1).EQ.'!') GOTO 2000
  IF (WORD1(1:1).EQ.'#') GOTO 2000
  IF (WORD1(1:1).EQ.'*') GOTO 9900
  IF (WORD1(1:1).EQ.'.') THEN
     DO I = 1, NTABLE
        IF (WORD1.EQ.TABLE(I)) THEN
           GOTO(3010,3020),I
        ENDIF
     ENDDO
  ENDIF
  GOTO 9100

! Reading vaues : .PRINT  
! ----------------------
  3010 CONTINUE 
  READ (LUCMD,'(A7)',ERR=9000) WORD1
  IF (WORD1(1:1).EQ.'!') GOTO 3010
  IF (WORD1(1:1).EQ.'#') GOTO 3010
  READ (WORD1,*,ERR=3011) PRINT_LEVEL 
  WRITE (RECP_OUT,'(/,A)') ' * Print level for RECP : ',PRINT_LEVEL
  GOTO 2000
! Error message
  3011 CONTINUE
  WRITE(RECP_OUT,'(/,A,I7)') ' .PRINT value error : ',PRINT_LEVEL
  GOTO 9000

! Reading vaues : .DEBUG 
! ----------------------
  3020 CONTINUE
  READ (LUCMD,'(A7)',ERR=9000) WORD1
  IF (WORD1(1:1).EQ.'!') GOTO 3020
  IF (WORD1(1:1).EQ.'#') GOTO 3020
  READ (WORD1,*,ERR=3021) RECP_DBG
  GOTO 2000
! Error message
  3021 CONTINUE
  WRITE(RECP_OUT,'(/,A,I7)') ' .DEBUG value error : ',RECP_DBG
  GOTO 9000

! Error message
! -------------
  9000 CONTINUE 
     WRITE(RECP_OUT,'(/,2A)') ' Error in reading menu file : ',WORD1
     CALL QUIT('Error in reading line')
  9100 CONTINUE
     WRITE(RECP_OUT,'(/,3A)') ' Keyword ',WORD1,' is not recognized in ECP'
     CALL QUIT('Wrong input')
  9900 CONTINUE 
END SUBROUTINE RECP_INP_MENU

END MODULE RECP_INP_READ
