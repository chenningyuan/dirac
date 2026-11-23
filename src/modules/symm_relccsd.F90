!-----------------------------------------
! Module with variables previously defined
! at src/relccsd/symm.inc include file
!-----------------------------------------
!  The essential double group symmetry information.
!
!  MXREP is the number of fermion or boson irreps in the Abelian
!  subgroup. It is doubled at some places (MULTB) where we need
!  both kind of irreps at the same time. Usually loops run
!  only over either fermion or boson irreps at a time.
!
!  NAOPI,NAVPI - orbitals per Abelian symmetry assigned for the IH-FSCC Pi space
!
!  NE..NVVVT2 (/SYMD/): sub dimensions of symmetry packed arrays;
!  /SYMP8/, /SYMP4/: pointers for symmetry packed arrays;
!  /SYMFSD/, /SYMFSP/: for Fock space (4 classes of spinors).
!  The TARGET attributes of the header are kept (pointers point to these arrays).
   module symm_relccsd_mod
   implicit none
   integer, parameter :: MXREP = 32
!  /SYMI/
   integer :: NSYMRP, NREP, MULTB(MXREP*2,MXREP*2,3)
!  /SYMC/
   character(len=4) :: REPN(16), REPNA(MXREP*2)
!  /SYML/
   logical :: SPFR
!  /DIMS4/
   integer :: NDIMT1, NDIMT2, NFOO, NFVO, NFVV, NV1, NV2, NV3, NV4, NDIMX, NDIMY, NORBT
!  /DIMS8/
   integer*8 :: NV5, NV6
!  /SYMD/
   integer, target :: NE(MXREP), NO(MXREP), NV(MXREP), NC(MXREP), NCONT(MXREP), MOO(MXREP), &
                      MVO(MXREP), MOV(MXREP), MVV(MXREP), POO(MXREP), PVO(MXREP), POV(MXREP), &
                      PVV(MXREP), NOO(MXREP), NOV(MXREP), NVO(MXREP), NVV(MXREP), NOOT(MXREP), &
                      NVVT(MXREP), NOOOT(MXREP), NVOOT(MXREP), NVVOT(MXREP), NVVVT(MXREP), &
                      NOVO(MXREP), NVOO(MXREP), NOVVT(MXREP), NOOVT(MXREP), NVOV(MXREP), &
                      NVVO(MXREP), NVOV2(MXREP), NOVO2(MXREP), NOVV(MXREP), NOOV(MXREP), &
                      NOOOT2(MXREP), NVVVT2(MXREP)
!  /SYMP8/
   integer*8, target :: JVVVO(MXREP+1), IVOVVT(MXREP+1), IVVVVTT(MXREP+1), KVVVOT(MXREP+1), &
                        KT3VVVOT(MXREP+1), KT3VVVVT(MXREP+1)
!  /SYMP4/
   integer :: IO(MXREP+1), IV(MXREP+1), IOO(MXREP+1), IVO(MXREP+1), IVV(MXREP+1), JVOXX(MXREP+1), &
              IVVVT(MXREP+1,MXREP), JVVVT(MXREP+1,MXREP), IVVVTMX, NT3VVT(MXREP), NT3VVVT(MXREP), &
              LLT3VVVT(MXREP,MXREP)
   integer, target :: JOOOO(MXREP+1), JVOOO(MXREP+1), JVVOO(MXREP+1), JVOVO(MXREP+1), &
                      JVVVOI4(MXREP+1), J2VOVO(MXREP+1), J2VOOO(MXREP+1), IOVVO(MXREP+1), &
                      J2VOVV(MXREP+1), J2OOOO(MXREP+1), J2VVOO(MXREP+1), J2VVVO(MXREP+1), &
                      IOOOOTT(MXREP+1), IVOOOT(MXREP+1), IVVOOTT(MXREP+1), IVOVO(MXREP+1), &
                      IVOOV(MXREP+1), IVVOO(MXREP+1), IVVOOT(MXREP+1), KVOVO(MXREP+1), &
                      KVVOOT(MXREP+1), KOOOVT(MXREP+1), KVOOV(MXREP+1), KVOOO(MXREP+1), &
                      KOOOOT(MXREP+1), KOOVVT(MXREP+1), KOOVOT(MXREP+1), KVOVV(MXREP+1), &
                      KOVOO(MXREP+1), KVVOVT(MXREP+1), IOVVVT(MXREP+1), IOVOOT(MXREP+1), &
                      LOOOOT(MXREP+1), LVVOOT(MXREP+1), LOVOOT(MXREP+1), LVOOOT(MXREP+1), &
                      LOOVO(MXREP+1), LVOVO(MXREP+1), LOVVOT(MXREP+1), LVOOVT(MXREP+1), &
                      IIOOT(MXREP,MXREP), IIOV(MXREP,MXREP), IIVO(MXREP,MXREP), IIVV(MXREP,MXREP), &
                      IIVVT(MXREP,MXREP), LOVVO(MXREP+1), JJOO(MXREP,MXREP), JJOV(MXREP,MXREP), &
                      LVVVO(MXREP+1), JJVO(MXREP,MXREP), JJVV(MXREP,MXREP), J2VOOV(MXREP+1), &
                      KKOOOT(MXREP,MXREP), KKVVOT(MXREP,MXREP), KKVOV(MXREP,MXREP), &
                      KKVOO(MXREP,MXREP), KKOOVT(MXREP,MXREP), KKVVVT(MXREP,MXREP), &
                      LLOOOT(MXREP,MXREP), LLVOOT(MXREP,MXREP), LLOVO(MXREP,MXREP), &
                      LOOVVT(MXREP+1), LVOVVT(MXREP+1), JVOOV(MXREP+1), LLOVVT(MXREP,MXREP), &
                      LLVVVT(MXREP,MXREP), LLVVO(MXREP,MXREP), LOVOV(MXREP+1), LLOOV(MXREP,MXREP), &
                      LLVOV(MXREP,MXREP), KKOVO(MXREP,MXREP), KKOVV(MXREP,MXREP), IOOOO(MXREP+1), &
                      IIOO(MXREP,MXREP), JOOOV(MXREP+1), JOOVO(MXREP+1), JOOVV(MXREP+1), &
                      JOVOO(MXREP+1), JOVOV(MXREP+1), JOVVV(MXREP+1), JVOVV(MXREP+1), &
                      JOVVO(MXREP+1), JVVOV(MXREP+1), J2OOOV(MXREP+1), J2OOVO(MXREP+1), &
                      J2OOVV(MXREP+1), J2OVOO(MXREP+1), J2OVOV(MXREP+1), J2OVVV(MXREP+1), &
                      J2OVVO(MXREP+1), J2VVOV(MXREP+1), KOVVOT(MXREP+1), KOVOV(MXREP+1)
!  /SYMFSD/
   integer :: NIO(MXREP), NIV(MXREP), NAO(MXREP), NAV(MXREP), NSO(MXREP), NSV(MXREP), NAOPI(MXREP), &
              NAVPI(MXREP), NIVIO(MXREP), NAVAO(MXREP), NAVSO(MXREP), NSOAV(MXREP), NAOV(MXREP), &
              MAOAV(MXREP), NIOIOT(MXREP), NAOAOT(MXREP), NAVAVT(MXREP), NIVIVT(MXREP), &
              NIVIVIOT(MXREP), NIOIOIVT(MXREP), NSVSVT(MXREP), NSOSOT(MXREP), NSVSVOT(MXREP), &
              NVSOSOT(MXREP), NAOSOSOT(MXREP), NSVSVAVT(MXREP), LTSHFS
!  /SYMFSP/
   integer :: IAVAV(MXREP+1), IAOAO(MXREP+1), IAVSO(MXREP+1), KIVIVIOAVT(MXREP+1), &
              KIOIOIVAOT(MXREP+1), KVVOAVT(MXREP+1), KSVSVAVAVT(MXREP+1), LAOAOSOSOT(MXREP+1), &
              LAOVOOT(MXREP+1), IISOAV(MXREP,MXREP), KKIVIVIOT(MXREP,MXREP), JAOAVAOAV(MXREP+1), &
              JAOAVOV(MXREP+1), JJAOAV(MXREP,MXREP), KKIOIOIVT(MXREP,MXREP), KKSVSVOT(MXREP,MXREP), &
              KKSVSVAVT(MXREP,MXREP), LLVSOSOT(MXREP,MXREP), LLAOSOSOT(MXREP,MXREP), &
              NFS_HE(MXREP,6), IFS_HE(MXREP+1,6), IFS_HE2(MXREP+1,6)
   end module symm_relccsd_mod
