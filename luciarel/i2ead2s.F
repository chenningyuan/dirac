      FUNCTION I2EAD2S(ISM,IORB,JSM,JORB,KSM,KORB,LSM,LORB,IJKLOF)
*
* Obtain address of integral (IORB JORB ! KORB LORB) in MOLCAS order
* IORB JORB KORB LORB corresponds to SYMMETRY ordered indeces !!
* and are relative to start of oribitals with given symmetry
* Integrals assumed in core
*
*
* Constructed for the MOLECULE-LUCIA interface for dbg calculations
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cgas.inc"
#include "csm.inc"
*
      Dimension IJKLOF(NsmSH,NsmSH,NsmSH)
      Logical iSymj,kSyml,ISYMK,JSYML,ijSymkl,IKSYMJL
      Logical ijklPerm
*.
      NTEST = 000
*
      If( NTEST.ge. 100) THEN
        write(6,*) ' GMIJKL at your service '
        WRITE(6,*) ' IORB ISM ',IORB,ISM
        WRITE(6,*) ' JORB JSM ',JORB,JSM
        WRITE(6,*) ' KORB KSM ',KORB,KSM
        WRITE(6,*) ' LORB LSM ',LORB,LSM
      END IF
*
      IABS = 0
      JABS = 0
      KABS = 0
      LABS = 0
*
      If (jSm.gt.iSm .or. (iSm.eq.jSm .and. JABS.gt.IABS)) Then
        iSym=jSm
        jSym=iSm
        I = JORB
        J = IORB
      Else
        iSym=iSm
        jSym=jSm
        I = IORB
        J = JORB
      End If
      ijBlk=jSym+iSym*(iSym-1)/2
      If (lSm.gt.kSm  .or. (kSm.eq.lSm .and. LABS.gt.KABS)) Then
        kSym=lSm
        lSym=kSm
        K = LORB
        L = KORB
      Else
        kSym=kSm
        lSym=lSm
        K = KORB
        L = LORB
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
      itOrb=NGSSHT(iSym)
      jtOrb=NGSSHT(jSym)
      ktOrb=NGSSHT(kSym)
      ltOrb=NGSSHT(lSym)
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
      I2EAD2S = iblOff-1+ijkl
      If( NTEST .ge. 100 ) then
	write(6,*) 'i j k l ', i,j,k,l
	write(6,*) ' ibloff ijkl ',ibloff,ijkl
        write(6,*) ' I2EADS  = ', I2EADS
      END IF
*
      RETURN
      END
