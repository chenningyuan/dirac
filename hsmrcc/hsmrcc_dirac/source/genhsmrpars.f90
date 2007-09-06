subroutine genHSMRpars

   use Parameters
   use OrbMod
   use IOModule
   implicit none
   integer i
   integer N_CHOOSE_R

   NHMO = NCMO + NVMO
   NPMO = NEMO + NVMO
   NELE = NCMO + NHOLE
   NTMO = NCMO + NVMO + NEMO
   NVMO = NHOLE + NPTCL

   NHMOSQ = ((NHMO*(NHMO-1))/2)
   NPMOSQ = ((NPMO*(NPMO-1))/2)
   NVMOSQ = ((NVMO*(NVMO-1))/2)

   NSING  = NHMO*NPMO
   NDBLE  = NHMOSQ*NPMOSQ
   NTSIZE = NSING + NDBLE

   NMAINORB = 0
   do i = 1, NVMO
      if (IsMainOrb(i)) NMAINORB = NMAINORB + 1
   enddo

   NINTMORB = 0
   do i = 1, NVMO
      if (IsIntmOrb(i)) NINTMORB = NINTMORB + 1
   enddo

   if ((NMAINORB+NINTMORB).ne.NVMO) then
      write(UOUT,*) "unaccounted orbitals.. errors somewhere"
      call flush(UOUT)
      call abort
   endif

   NINTMORB_CORE = 0
   NINTMORB_VIRT = 0
   do i = 1, NVMO
      if (IsIntmCore(i)) NINTMORB_CORE = NINTMORB_CORE + 1
      if (IsIntmVirt(i)) NINTMORB_VIRT = NINTMORB_VIRT + 1
   enddo

   if ((NINTMORB_CORE+NINTMORB_VIRT).ne.NINTMORB) then
      write(UOUT,*) "something wrong with intemediate orbital space"
      call flush(UOUT)
      call abort
   endif

   if (NHOLE.LT.NINTMORB_CORE) then
      write(UOUT,*) "invalid number of active electrons"
      call flush(UOUT)
      call abort
   endif

   NSTATE = N_CHOOSE_R(NVMO,NHOLE)
   NSTATE_MAIN = N_CHOOSE_R(NMAINORB,(NHOLE-NINTMORB_CORE))
   NSTATE_INTM = NSTATE - NSTATE_MAIN

   MAXEXC = MIN(NHOLE,NPTCL)

end subroutine GenHSMRpars
