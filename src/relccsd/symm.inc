!
!
!     The essential double group symmetry information.
!
!     MXREP is the number of fermion or boson irreps in the Abelian
!     subgroup. It is doubled at some places (MULTB) where we need
!     both kind of irreps at the same time. Usually loops run
!     only over either fermion or boson irreps at a time.
!
!
!     NAOPI,NAVPI - orbitals per Abelian symmetry assigned for the IH-FSCC Pi space
!         MXREP is from symm.inc
!

      INTEGER :: MXREP
      PARAMETER (MXREP=32)

      INTEGER :: NZ, NSYMRP, NREP, MULTB
      COMMON /SYMI/ NZ,NSYMRP,NREP,MULTB(MXREP*2,MXREP*2,2)

      COMMON /SYMC/ REPN(16),REPNA(MXREP*2) 
      CHARACTER*4 REPN,REPNA
      COMMON /SYML/ SPFR
      LOGICAL SPFR

      INTEGER :: NDIMT1, NDIMT2, NFOO, NFVO, NFVV, NV1, NV2, NV3, NV4
      INTEGER :: NDIMX, NDIMY, NORBT
      COMMON /DIMS4/ NDIMT1,NDIMT2,NFOO,NFVO,NFVV,NV1,NV2,NV3,NV4       &
     &              ,NDIMX,NDIMY,NORBT
      INTEGER*8 NV5,NV6
      COMMON /DIMS8/ NV5,NV6

!     Sub dimensions of symmetry packed arrays.
      INTEGER,TARGET :: NE, NO, NV, NC
      INTEGER,TARGET :: MOO, MVO, MOV, MVV
      INTEGER,TARGET :: NOO, NOV, NVO, NVV
      INTEGER,TARGET :: NOOT, NVVT
      INTEGER,TARGET :: NOOOT, NVOOT, NVVOT, NVVVT, NOVO
      INTEGER,TARGET :: NVOO, NOVVT, NOOVT
      COMMON /SYMD/ NE(MXREP),NO(MXREP),NV(MXREP),NC(MXREP),            &
     & MOO(MXREP),MVO(MXREP),MOV(MXREP),MVV(MXREP),                     &
     & NOO(MXREP),NOV(MXREP),NVO(MXREP),NVV(MXREP),                     &
     & NOOT(MXREP),NVVT(MXREP),                                         &
     & NOOOT(MXREP),NVOOT(MXREP),NVVOT(MXREP),NVVVT(MXREP),NOVO(MXREP), &
     & NVOO(MXREP),NOVVT(MXREP),NOOVT(MXREP)

!     Pointers for symmetry packed arrays.
      INTEGER*8,TARGET :: JVVVO,IVOVVT,IVVVVTT,KVVVOT,KT3VVVOT                 
      COMMON /SYMP8/ JVVVO(MXREP+1),IVOVVT(MXREP+1),                    &
     & IVVVVTT(MXREP+1),KVVVOT(MXREP+1),KT3VVVOT(MXREP+1)               

      INTEGER :: IO, IV, IOO
      INTEGER :: IVO, IVV, JVOXX
      INTEGER,TARGET :: JOOOO, JVOOO, JVVOO, JVOVO
      INTEGER,TARGET :: J2VOVO, J2VOOO
      INTEGER,TARGET :: J2VOVV, J2OOOO, J2VVOO, J2VVVO
      INTEGER :: IOOOOTT, IVOOOT, IVVOOTT, IVOVO
      INTEGER :: IVVOO, IVVOOT
      INTEGER,TARGET :: KVVOOT, KOOOVT, KVOOV, KVOOO
      INTEGER,TARGET :: KOOOOT, KOOVVT, KOOVOT
      INTEGER,TARGET :: LOOOOT, LVVOOT, LOVOOT, LVOOOT
      INTEGER,TARGET :: LOOVO, LVOVO, LOVVOT, LVOOVT
      INTEGER :: IIOOT, IIOV, IIVO
      INTEGER :: IIVV, IIVVT
      INTEGER,TARGET :: JJOO, JJOV
      INTEGER,TARGET :: JJVO, JJVV
      INTEGER,TARGET :: KKOOOT, KKVVOT, KKVVVT
      INTEGER,TARGET :: KKVOO, KKOOVT
      INTEGER,TARGET :: LLOOOT, LLVOOT, LLOVO
      INTEGER,TARGET :: LOOVVT, LVOVVT
      INTEGER,TARGET :: LLOVVT, LLVVVT
      INTEGER :: IVVVT, JVVVT, IVVVTMX
      INTEGER :: NT3VVT, NT3VVVT, LLT3VVVT
      COMMON /SYMP4/ IO(MXREP+1),IV(MXREP+1),IOO(MXREP+1),              &
     &               IVO(MXREP+1),IVV(MXREP+1),JVOXX(MXREP+1),          &
     & JOOOO(MXREP+1),JVOOO(MXREP+1),JVVOO(MXREP+1),JVOVO(MXREP+1),     &
     & J2VOVO(MXREP+1),J2VOOO(MXREP+1),                                 &
     & J2VOVV(MXREP+1),J2OOOO(MXREP+1),J2VVOO(MXREP+1),J2VVVO(MXREP+1), &
     & IOOOOTT(MXREP+1),IVOOOT(MXREP+1),IVVOOTT(MXREP+1),IVOVO(MXREP+1),&
     & IVVOO(MXREP+1),IVVOOT(MXREP+1),                                  &
     & KVVOOT(MXREP+1),KOOOVT(MXREP+1),KVOOV(MXREP+1),KVOOO(MXREP+1),   &
     & KOOOOT(MXREP+1),KOOVVT(MXREP+1),KOOVOT(MXREP+1),                 &
     & LOOOOT(MXREP+1),LVVOOT(MXREP+1),LOVOOT(MXREP+1),LVOOOT(MXREP+1), &
     & LOOVO(MXREP+1),LVOVO(MXREP+1),LOVVOT(MXREP+1), LVOOVT(MXREP+1),  &
     & IIOOT(MXREP,MXREP),IIOV(MXREP,MXREP),IIVO(MXREP,MXREP),          &
     & IIVV(MXREP,MXREP),IIVVT(MXREP,MXREP),                            &
     & JJOO(MXREP,MXREP),JJOV(MXREP,MXREP),                             &
     & JJVO(MXREP,MXREP),JJVV(MXREP,MXREP),                             &
     & KKOOOT(MXREP,MXREP),KKVVOT(MXREP,MXREP),                         &
     & KKVOO(MXREP,MXREP),KKOOVT(MXREP,MXREP),KKVVVT(MXREP,MXREP),      &
     & LLOOOT(MXREP,MXREP),LLVOOT(MXREP,MXREP),LLOVO(MXREP,MXREP),      &
     & LOOVVT(MXREP+1),LVOVVT(MXREP+1),                                 &
     & LLOVVT(MXREP,MXREP),LLVVVT(MXREP,MXREP),                         &
     & IVVVT(MXREP+1,MXREP),JVVVT(MXREP+1,MXREP),IVVVTMX,               &
     & NT3VVT(MXREP),NT3VVVT(MXREP),LLT3VVVT(MXREP,MXREP)   

!     For Fock space : we now have 4 classes of spinors
      INTEGER :: NIO, NIV, NAO, NAV
      INTEGER :: NSO, NSV, NAOPI, NAVPI
      INTEGER :: NIVIO, NAVAO, NAVSO, NSOAV, NAOV
      INTEGER :: MAOAV
      INTEGER :: NIOIOT, NAOAOT, NAVAVT, NIVIVT
      INTEGER :: NIVIVIOT, NIOIOIVT, NSVSVT, NSOSOT
      INTEGER :: NSVSVOT, NVSOSOT, NAOSOSOT, NSVSVAVT
      INTEGER :: LTSHFS
      COMMON /SYMFSD/ NIO(MXREP),NIV(MXREP),NAO(MXREP),NAV(MXREP),      &
     &                NSO(MXREP),NSV(MXREP),                            &
     &                NAOPI(MXREP),NAVPI(MXREP),                        &
     & NIVIO(MXREP),NAVAO(MXREP),NAVSO(MXREP),NSOAV(MXREP),NAOV(MXREP), &
     & MAOAV(MXREP),                                                    &
     & NIOIOT(MXREP),NAOAOT(MXREP),NAVAVT(MXREP),NIVIVT(MXREP),         &
     & NIVIVIOT(MXREP),NIOIOIVT(MXREP),NSVSVT(MXREP),NSOSOT(MXREP),     &
     & NSVSVOT(MXREP),NVSOSOT(MXREP),NAOSOSOT(MXREP),NSVSVAVT(MXREP),   &
     & LTSHFS

      INTEGER :: IAVAV, IAOAO, IAVSO
      INTEGER :: KIVIVIOAVT, KIOIOIVAOT, KVVOAVT
      INTEGER :: KSVSVAVAVT, LAOAOSOSOT
      INTEGER :: LAOVOOT, IISOAV, KKIVIVIOT
      INTEGER :: JAOAVAOAV, JAOAVOV
      INTEGER :: JJAOAV
      INTEGER :: KKIOIOIVT, KKSVSVOT
      INTEGER :: KKSVSVAVT
      INTEGER :: LLVSOSOT, LLAOSOSOT
      INTEGER :: NFS_HE, IFS_HE, IFS_HE2
      COMMON /SYMFSP/ IAVAV(MXREP+1),IAOAO(MXREP+1),IAVSO(MXREP+1),     &
     & KIVIVIOAVT(MXREP+1),KIOIOIVAOT(MXREP+1),KVVOAVT(MXREP+1),        &
     & KSVSVAVAVT(MXREP+1),LAOAOSOSOT(MXREP+1),                         &
     & LAOVOOT(MXREP+1),IISOAV(MXREP,MXREP),KKIVIVIOT(MXREP,MXREP),     &
     & JAOAVAOAV(MXREP+1),JAOAVOV(MXREP+1),                             &
     & JJAOAV(MXREP,MXREP),                                             &
     & KKIOIOIVT(MXREP,MXREP),KKSVSVOT(MXREP,MXREP),                    &
     & KKSVSVAVT(MXREP,MXREP),                                          &
     & LLVSOSOT(MXREP,MXREP),LLAOSOSOT(MXREP,MXREP),                    &
     & NFS_HE(MXREP,6),IFS_HE(MXREP+1,6),IFS_HE2(MXREP+1,6) 
