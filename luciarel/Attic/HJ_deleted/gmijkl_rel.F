      FUNCTION GMIJKL_REL(IORB,JORB,KORB,LORB,INTLST,IJKLOF)
*
* Obtain integral (IORB JORB ! KORB LORB) MOLCAS version
* Integrals assumed in core
*
* Version for integrals stored in INTLST
*
      IMPLICIT REAL*8(A-H,O-Z)
*.ORBINP
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB),
C    &              NOBPTS(6+MXPR4T,MXPOBS),IOBPTS(6+MXPR4T,MXPOBS),
C    &              ITOOBS(MXPOBS)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
*. Integral list
      Real * 8 Intlst(*)
      Dimension IJKLOF(NsmOB,NsmOb,NsmOB)
      Logical iSymj,kSyml,ISYMK,JSYML,ijSymkl,IKSYMJL
      Logical ijklPerm
*.
      NTEST = 000
*
*. The orbital list corresponds to type ordered indeces, reform to
*. symmetry ordering
*
      IABS = IREOTS(IORB)
      ISM = ISMFTO(IORB)
      IOFF = IBSO(ISM)
*
      JABS = IREOTS(JORB)
      JSM = ISMFTO(JORB)
      JOFF = IBSO(JSM)
*
      KABS = IREOTS(KORB)
      KSM = ISMFTO(KORB)
      KOFF = IBSO(KSM)
*
      LABS = IREOTS(LORB)
      LSM = ISMFTO(LORB)
      LOFF = IBSO(LSM)
*
      If( Ntest.ge. 100) THEN
        write(6,*) ' GMIJKL at your service '
        WRITE(6,*) ' IORB IABS ISM IOFF ',IORB,IABS,ISM,IOFF
        WRITE(6,*) ' JORB JABS JSM JOFF ',JORB,JABS,JSM,JOFF
        WRITE(6,*) ' KORB KABS KSM KOFF ',KORB,KABS,KSM,KOFF
        WRITE(6,*) ' LORB LABS LSM LOFF ',LORB,LABS,LSM,LOFF
      END IF
*
      If ( jSm.gt.iSm .or. ( iSm.eq.jSm .and. JABS.gt.IABS)) Then
        iSym=jSm
        jSym=iSm
        I = JABS - JOFF + 1
        J = IABS - IOFF + 1
      Else
        iSym=iSm
        jSym=jSm
        I = IABS - IOFF + 1
        J = JABS - JOFF + 1
      End If
      ijBlk=jSym+iSym*(iSym-1)/2
      If ( lSm.gt.kSm  .or. ( kSm.eq.lSm .and. LABS.gt.KABS)) Then
        kSym=lSm
        lSym=kSm
        K = LABS -LOFF + 1
        L = KABS - KOFF + 1
      Else
        kSym=kSm
        lSym=lSm
        K = KABS - KOFF + 1
        L = LABS -LOFF + 1
      End If
      klBlk=lSym+kSym*(kSym-1)/2
*
      ijklPerm=.false.
      If ( klBlk.gt.ijBlk ) Then
        iTemp=iSym
        iSym=kSym
        kSym=iTemp
        iTemp=jSym
        jSym=lSym
        lSym=iTemp
        iTemp=ijBlk
        ijBlk=klBlk
        klBlk=iTemp
        ijklPerm=.true.
*
        iTemp = i
        i = k
        k = itemp
        iTemp = j
        j = l
        l = iTemp
      End If
      If(Ntest .ge. 100 ) then
        write(6,*) ' i j k l ',i,j,k,l
        write(6,*) ' Isym,Jsym,Ksym,Lsym',Isym,Jsym,Ksym,Lsym
      End if
*
*  Define offset for given symmetry block
      IBLoff = IJKLof(Isym,Jsym,Ksym)
      If(ntest .ge. 100 )
     &WRITE(6,*) ' IBLoff Isym Jsym Ksym ', IBLoff,ISym,Jsym,Ksym
      iSymj=iSym.eq.jSym
      kSyml=kSym.eq.lSym
      iSymk=iSym.eq.kSym
      jSyml=jSym.eq.lSym
      ikSymjl=iSymk.and.jSyml
      ijSymkl=iSymj.and.kSyml
*
      itOrb=NTOOBS(iSym)
      jtOrb=NTOOBS(jSym)
      ktOrb=NTOOBS(kSym)
      ltOrb=NTOOBS(lSym)
C?    print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
      If ( iSymj ) Then
        ijPairs=itOrb*(itOrb+1)/2
        ij=j+i*(i-1)/2
      Else
        ijPairs=itOrb*jtOrb
        ij=j + (i-1)*jtOrb
      End if
*
      IF(KSYML ) THEN
        klPairs=ktOrb*(ktOrb+1)/2
        kl=l+k*(k-1)/2
      ELSE
        klPairs=ktOrb*ltOrb
        kl=l+(k-1)*ltOrb
      End If
C?    print *,' ijPairs,klPairs',ijPairs,klPairs
*
      If ( ikSymjl ) Then
        If ( ij.gt.kl ) Then
          klOff=kl+(kl-1)*(kl-2)/2-1
          ijkl=ij+(kl-1)*ijPairs-klOff
        Else
          ijOff=ij+(ij-1)*(ij-2)/2-1
          ijkl=kl+(ij-1)*klPairs-ijOff
        End If
      Else
        ijkl=ij+(kl-1)*ijPairs
      End If
      If( ntest .ge. 100 )
     & write(6,*) ' ijkl ', ijkl
*
      GMIJKL_REL = Intlst(iblOff-1+ijkl)
      If( ntest .ge. 100 )
     & write(6,*) ' GMIJKL_REL ', GMIJKL_REL
*
      RETURN
      END
