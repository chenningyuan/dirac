      subroutine rdmole2(SWEINT,NORBMOL,NBASMOL,NFRO,IGLSUM,XLUCINT)
*
* Sweden integral list        => Lucia integral list
* ^ Read from file to SWEINT     ^ Saved in array(sic) XLUCINT
*
*
*  Patch from tradpt (Christel Marian)
*  Timo Fleig
*
      implicit real*8 (A-H,O,P,R-Z)
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
#include "clunit.inc"
*
      parameter (LOCLEN=666)
*. Input
      dimension NORBMOL(NIRR_PN),NBASMOL(NIRR_PN)
      dimension SWEINT(*)
      dimension IOFFSET(LOCLEN)
      dimension NFRO(MXPIRR),QNTIL(MXPIRR),QMJ(MXPIRR),NIT(LOCLEN),
     &          QKJ(MXPIRR),QLJ(MXPIRR),QNJ(MXPIRR),QLG(MXPIRR),
     &          QISYM(MXPIRR),IDEKS(MXPORB)
*. Output
      DIMENSION XLUCINT(*)
*
      NTEST = 500
*
*  Set up triangular index array
      II = 0
      do I=1,MXPORB,1
        IDEKS(I) = II
        II = II + I
      end do
*
*  If DOUGRP is C2H, then some reordering has to be done.
*  First: Rename dimensions of symmetries
*         Go from LUCIAREL to MOLECULE ordering (Ag,Au,Bu,Bg)
*  Second: Reorder integrals according to LUCIAREL symmetry ordering
*          (done below)
*
      if (DOUGRP.eq.4) then
         IREOTEMP4 = NORBMOL(4)
         IREOTEMP3 = NORBMOL(3)
         IREOTEMP2 = NORBMOL(2)
         NORBMOL(2) = IREOTEMP3
         NORBMOL(3) = IREOTEMP4
         NORBMOL(4) = IREOTEMP2
*
         IREOTEMP4 = NBASMOL(4)
         IREOTEMP3 = NBASMOL(3)
         IREOTEMP2 = NBASMOL(2)
         NBASMOL(2) = IREOTEMP3
         NBASMOL(3) = IREOTEMP4
         NBASMOL(4) = IREOTEMP2
*
         IREOTEMP4 = NFRO(4)
         IREOTEMP3 = NFRO(3)
         IREOTEMP2 = NFRO(2)
         NFRO(2) = IREOTEMP3
         NFRO(3) = IREOTEMP4
         NFRO(4) = IREOTEMP2
      end if
*
*  CALCULATE VALUES FOR 1 ST OUTPUT-RECORD, PACK THEM TO OCCUPIED
*  IRREPS
*
      NOCC=0
      NBOX=0
      NCO=0
      LCT=0
      QNTIL(1)=0
      IORBS=NBASMOL(1)
      do 30 IS=1,NIRR_PN
        IORBS=IORBS+NBASMOL(IS)
        NCV=NFRO(IS)+NORBMOL(IS)
        if(NCV.EQ.0) goto 20
        NOCC=NOCC+1
*
*  # MOs - DELETED / OCCUPIED SYMMETRY
*
        QMJ(NOCC)=NCV
*
*  # FROZEN MO'S  / OCCUPIED SYMMETRY
*
        QKJ(NOCC)=NFRO(IS)
        NCO=NCO+NFRO(IS)
*
*  # VALENCE+VIRTUAL MO'S  / OCCUPIED SYMMETRY
*
        NO=NORBMOL(IS)
        QLJ(NOCC)=NO
        NBOX=NBOX+NO
*
*  # SOs / OCCUPIED SYMMETRY
*
        QNJ(NOCC)=NBASMOL(IS)
*
*  SYMMETRY-PACKING LABELS
*
        QISYM(NOCC)=IS
*
*  LAST MO IN EACH SYMMETRY
*
20      if(IS.NE.NIRR_PN) QNTIL(IS+1)=QNTIL(IS)+NO
30    CONTINUE
*
*  CALCULATE NUMBER OF INTEGRALS (NIGL) FOR EACH SYMMETRY-BLOCK
*  AND STORE 0-ADRESSES IN NIT
*
      NIT(1) = 0
      IOFFSET(1) = 1
      INIT = 1
      NPQSUM = 0
      IGLSUM = 0
*
*
* Calculate total number of integrals on record
* Determine old offsets for MOLECULE
      do ISP=1,NIRR_PN,1
         NSP=QISYM(ISP)
         if (NSP.eq.0) then
            write(6,*) 'Number of symm for P is  ',NSP
            write(6,*) 'Skipping two-el. int. read.'
            goto 333
         end if
         NOP=NORBMOL(NSP)
         do ISQ=1,ISP,1
            NSQ = QISYM(ISQ)
            NOQ = NORBMOL(NSQ)
            ISPQ = IPGMULT(ISP,ISQ)
            do ISR=1,ISP,1
               NSR = QISYM(ISR)
               NOR = NORBMOL(NSR)
               ISPQR = IPGMULT(ISPQ,ISR)
               ISSM=ISR
               if (ISR.EQ.ISP) ISSM=ISQ
               do 40 ISS=1,ISSM,1
                  NSS=QISYM(ISS)
                  NOS=NORBMOL(NSS)
                  if(NSS.NE.ISPQR) goto 40
                  if(NOS*NOR*NOQ*NOP.eq.0) goto 40
*
*  NUMBER OF INTEGRALS IN THIS SYMMETRYBLOCK
*
                  if (NSP.EQ.NSQ) then
*
*  TYPES <PP|PP> AND <PP|QQ>
*
                    NPQ = IDEKS(NOP+1)
                    NRS = IDEKS(NOR+1)
                  else
*
*  TYPES <PQ|PQ> AND <PQ|RS>
*
                    NPQ = NOP*NOQ
                    NRS = NOR*NOS
                  end if
                  if (NSP.EQ.NSR) then
*
*  TYPES <PP|PP> AND <PQ|PQ>
*
                    NIGL = NPQ*(NPQ+1)/2
                  else
*
*  TYPES <PP|QQ> AND <PQ|RS>
*
                    NIGL = NPQ*NRS
                  end if
                  NIT(INIT) = NIGL
                  INIT = INIT + 1
                  IOFFSET(INIT) = IOFFSET(INIT-1) + NIGL
                  IGLSUM = IGLSUM + NIGL
                  if (NTEST.ge.100) then
                     write(6,601) NSP,NSQ,NSR,NSS,NIGL
                  end if
40             continue
            end do
         end do
* Number of non-vanishing symmetry blocks
         NSMBLK = INIT - 1
         write(6,*) '# Symm blocks: ',NSMBLK
333      continue
      end do
601   FORMAT(' NUMBER OF INTEGRALS IN SYMMETRYBLOCK '
     &             ,4I2,' IS : ',I10)
*
*  Restore LUCIAREL ordering
*
      if (DOUGRP.eq.4) then
         IREOTEMP4 = NORBMOL(4)
         IREOTEMP3 = NORBMOL(3)
         IREOTEMP2 = NORBMOL(2)
         NORBMOL(2) = IREOTEMP4
         NORBMOL(3) = IREOTEMP2
         NORBMOL(4) = IREOTEMP3
*
         IREOTEMP4 = NBASMOL(4)
         IREOTEMP3 = NBASMOL(3)
         IREOTEMP2 = NBASMOL(2)
         NBASMOL(2) = IREOTEMP4
         NBASMOL(3) = IREOTEMP2
         NBASMOL(4) = IREOTEMP3
*
         IREOTEMP4 = NFRO(4)
         IREOTEMP3 = NFRO(3)
         IREOTEMP2 = NFRO(2)
         NFRO(2) = IREOTEMP4
         NFRO(3) = IREOTEMP2
         NFRO(4) = IREOTEMP3
*
* Do the same stuff for LUCIAREL offsets
*
         NOCC=0
         NBOX=0
         NCO=0
         LCT=0
         QNTIL(1)=0
         IORBS=NBASMOL(1)
         do 50 IS=1,NIRR_PN
           IORBS=IORBS+NBASMOL(IS)
           NCV=NFRO(IS)+NORBMOL(IS)
           if(NCV.EQ.0) goto 60
           NOCC=NOCC+1
*
*  # MOs - DELETED / OCCUPIED SYMMETRY
*
           QMJ(NOCC)=NCV
*
*  # FROZEN MO'S  / OCCUPIED SYMMETRY
*
           QKJ(NOCC)=NFRO(IS)
           NCO=NCO+NFRO(IS)
*
*  # VALENCE+VIRTUAL MO'S  / OCCUPIED SYMMETRY
*
           NO=NORBMOL(IS)
           QLJ(NOCC)=NO
           NBOX=NBOX+NO
*
*  # SOs / OCCUPIED SYMMETRY
*
           QNJ(NOCC)=NBASMOL(IS)
*
*  SYMMETRY-PACKING LABELS
*
           QISYM(NOCC)=IS
*
*  LAST MO IN EACH SYMMETRY
*
60         if(IS.NE.NIRR_PN) QNTIL(IS+1)=QNTIL(IS)+NO
50       CONTINUE
* Determine new offsets for LUCIAREL
         IOFFSET(INIT) = 1
         NPQSUM = 0
         IGLSUM_L = 0
         do ISP=1,NIRR_PN,1
            NSP=QISYM(ISP)
            if (NSP.eq.0) then
               write(6,*) 'Number of symm for P is  ',NSP
               write(6,*) 'Skipping two-el. int. read.'
               goto 334
            end if
            NOP=NORBMOL(NSP)
            write(6,*) 'ISP =  ',ISP
            write(6,*) 'NOP =   ',NOP
            do ISQ=1,ISP,1
               NSQ = QISYM(ISQ)
               NOQ = NORBMOL(NSQ)
               write(6,*) 'ISQ =  ',ISQ
               write(6,*) 'NOQ =   ',NOQ
CTF            ISPQ = IPGMULT(ISP,ISQ)
               ISPQ = IPGMULT(NSP,NSQ)
               do ISR=1,ISP,1
                  NSR = QISYM(ISR)
                  NOR = NORBMOL(NSR)
                  write(6,*) 'ISR =  ',ISR
                  write(6,*) 'NOR =   ',NOR
CTF               ISPQR = IPGMULT(ISPQ,ISR)
                  ISPQR = IPGMULT(ISPQ,NSR)
                  ISSM=ISR
                  if (ISR.EQ.ISP) ISSM=ISQ
                  do 80 ISS=1,ISSM,1
                     NSS=QISYM(ISS)
                     NOS=NORBMOL(NSS)
                     write(6,*) 'ISS =  ',ISS
                     write(6,*) 'NOS =   ',NOS
                     if (NSS.ne.ISPQR) then
                       NIGL = 0
                       goto 80
                     end if
                     if (NOS*NOR*NOQ*NOP.eq.0) then
                       NIGL = 0
                       goto 80
                     end if
*
*  NUMBER OF INTEGRALS IN THIS SYMMETRYBLOCK
*
                     if (NSP.EQ.NSQ) then
*
*  TYPES <PP|PP> AND <PP|QQ>
*
                       NPQ = IDEKS(NOP+1)
                       NRS = IDEKS(NOR+1)
                     else
*
*  TYPES <PQ|PQ> AND <PQ|RS>
*
                       NPQ = NOP*NOQ
                       NRS = NOR*NOS
                     end if
                     if (NSP.EQ.NSR) then
*
*  TYPES <PP|PP> AND <PQ|PQ>
*
                       NIGL = NPQ*(NPQ+1)/2
                     else
*
*  TYPES <PP|QQ> AND <PQ|RS>
*
                       NIGL = NPQ*NRS
                     end if
                     NIT(INIT) = NIGL
80                continue
                  IOFFSET(INIT+1) = IOFFSET(INIT) + NIGL
                  INIT = INIT + 1
                  IGLSUM_L = IGLSUM_L + NIGL
                  if (NTEST.ge.100) then
                     write(6,*) NSP,NSQ,NSR,NSS,NIGL
                     write(6,*) 'IGLSUM_L =  ',IGLSUM_L
                  end if
               end do
            end do
334         continue
         end do
      end if
*
*--------------------------------------------
* Read total integral record from file LU2INT
*--------------------------------------------
      IAD13 = 0
      LUSGITE = LUPRP
      if (ARCHI.eq.'RS6') then
         call dafile(LU2INT,2,SWEINT(1),2*IGLSUM,IAD13)
*
*  Write for SGI read
         if (IPRINT.ge.2) then
            rewind LUSGITE
            write(6,*) 'writing transformed 2-el. integrals'
            do I=1,IGLSUM,1
               write(6,*) I,SWEINT(I)
               write(LUSGITE,*) SWEINT(I)
            end do
         end if
*
      else if (ARCHI.eq.'SGO') then
*
* I'm running on an SGI. Please give me a file ...
         rewind LUSGITE
         open(unit=LUSGITE,file='TWOEL',
     &        status='OLD',form='FORMATTED')
         do I=1,IGLSUM,1
            read(LUSGITE,*) SWEINT(I)
         end do
         close(unit=LUSGITE)
      else
         write(6,*) 'Architecture not implemented.'
         write(6,*) 'I am quitting!'
         Call Abend()
      end if
*
*  Checking integrals previously read from LU 55 ...
      if (NTEST.ge.100) then
         write(6,*) 'Transformed spin-free 2-el. integrals:'
         do I=1,IGLSUM,1
            write(6,'(1G18.12)') SWEINT(I)
         end do
      end if
*
      if (DOUGRP.eq.4) then
* Do reordering of integrals if DOUGRP is C2H
         ILOCLEN = LOCLEN
         call reotwin_c2h(SWEINT(1),XLUCINT(1),
     &                    IOFFSET(1),IOFFSET(NSMBLK+1),NIT(1),
     &                    NIT(NSMBLK+1),IGLSUM,IGLSUM_L,
     &                    NSMBLK,NSMBLK_L,ILOCLEN)
      end if
*
* Now loop over symmetries again and pick relevant
* integral block from array SWEINT.
* Call reordering routine for internal LUCIA(REL)
* integral order.
*
*  P
*   Q
*    R
*     S
*
      INIT = 1
      IPOFF = 1
      do ISP=1,NIRR_PN,1
         NSP=QISYM(ISP)
         if (NSP.eq.0) then
            write(6,*) 'Number of symm for P is  ',NSP
            write(6,*) 'Skipping two-el. int. read.'
            goto 444
         end if
         NOP=NORBMOL(NSP)
         IQOFF = 1
         do ISQ=1,ISP,1
            NSQ = QISYM(ISQ)
            NOQ = NORBMOL(NSQ)
            ISPQ = IPGMULT(ISP,ISQ)
            IROFF = 1
            do ISR=1,ISP,1
               NSR = QISYM(ISR)
               NOR = NORBMOL(NSR)
               ISPQR = IPGMULT(ISPQ,ISR)
               ISSM=ISR
               if (ISR.EQ.ISP) ISSM=ISQ
               ISOFF = 1
               do 100 ISS=1,ISSM,1
                  NSS=QISYM(ISS)
                  NOS=NORBMOL(NSS)
                  if(NSS.NE.ISPQR) goto 49
                  if(NOS*NOR*NOQ*NOP.eq.0) goto 49
                  NINT = NIT(INIT)
*. offset in MOLECULE integrals to given symmetry block
                  IOFF = IOFFSET(INIT)
                  if (DOUGRP.eq.4) IOFF = IOFFSET(INIT+NSMBLK)
                  call reomoluc(SWEINT(IOFF),ISP,ISQ,ISR,ISS,
     &                          IPOFF,IQOFF,IROFF,ISOFF,
     &                          NOP,NOQ,NOR,NOS,
     &                          NIT(INIT+NSMBLK),XLUCINT)
*
                  INIT = INIT + 1
49                continue
                  ISOFF = ISOFF + NOS
100            continue
               IROFF = IROFF  + NOR
            end do
            IQOFF = IQOFF + NOQ
         end do
         IPOFF = IPOFF + NOP
444      continue
      end do
*
      write(6,*) 'IGLSUM (Number of ints) = ',IGLSUM
      INTNUM = IGLSUM
*
*  Check whether last integral is zero ...
      SWEABS = ABS(SWEINT(IGLSUM))
      CHECK = SWEABS - 1E-9
      if (CHECK.lt.0) then
        write(6,*) 'Problem with two-electron integral file'
        write(6,*) 'LASTABS is ',LASTABS
        write(6,*) 'CHECK is ',CHECK
        write(6,*) 'SWEINT(',IGLSUM,') is  ',SWEINT(IGLSUM)
        Call Abend2( 'Quitting in rdmole2' )
      end if
*
* Check integrals
      if (NTEST.ge.100) then
        write(6,*) 'Integrals on SWEINT after rdmole2: '
        write(6,*) 'Transformed spin-free 2-el. integrals:'
        do I=1,IGLSUM,1
          write(6,'(1G18.12)') SWEINT(I)
        end do
      end if
*
      return
      end
