!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine GenHSMRpars
!-----------------------------------------------------------------------
!
!   Generate set of parameters for HSMRCC run.
!
!   called from InitHSMRBase (source/inihsbase.F90)
!
!-----------------------------------------------------------------------

   use Parameters ! most of parameters are stored here ...
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
      write(UOUT,*) "GenHSMRpars: unaccounted orbitals.. errors somewhere"
      call flushout
      call quit('GenHSMRpars: unaccounted orbitals.. errors somewhere')
   endif

   NINTMORB_CORE = 0
   NINTMORB_VIRT = 0
   do i = 1, NVMO
      if (IsIntmCore(i)) NINTMORB_CORE = NINTMORB_CORE + 1
      if (IsIntmVirt(i)) NINTMORB_VIRT = NINTMORB_VIRT + 1
   enddo

   if ((NINTMORB_CORE+NINTMORB_VIRT).ne.NINTMORB) then
      write(UOUT,*) "GenHSMRpars: something wrong with intemediate orbital space"
      call flushout
      call quit('GenHSMRpars: something wrong with intemediate orbital space')
   endif

   if (NHOLE.LT.NINTMORB_CORE) then
      write(UOUT,*) "GenHSMRpars: invalid number of active electrons"
      call flushout
      call quit('GenHSMRpars: invalid number of active electrons')
   endif

   NSTATE = N_CHOOSE_R(NVMO,NHOLE)
   NSTATE_MAIN = N_CHOOSE_R(NMAINORB,(NHOLE-NINTMORB_CORE))
   NSTATE_INTM = NSTATE - NSTATE_MAIN

   MAXEXC = MIN(NHOLE,NPTCL)

   ! ... control print-out
   if ( get_iprint() >= 10 ) then
     call HEADER('Output from GenHSMRpars',-1)
     
     write(UOUT,'(2x,a,i7)') 'NSTATE=',NSTATE
     write(UOUT,'(2x,a,i7)') 'NSTATE_MAIN=',NSTATE_MAIN
     write(UOUT,'(2x,a,i7)') 'NSTATE_INTM=',NSTATE_INTM
     write(UOUT,'(2x,a,i7)') 'NAX_EXC=',MAXEXC

     write(UOUT,'(2x,a,i7)') 'NMAINORB=',NMAINORB
     write(UOUT,'(2x,a,i7)') 'NINTMORB=',NINTMORB

     write(UOUT,'(2x,a,i7)') 'NINTMORB_CORE =',NINTMORB_CORE
     write(UOUT,'(2x,a,i7)') 'NINTMORB_VIRT =',NINTMORB_VIRT

   endif


end subroutine GenHSMRpars
