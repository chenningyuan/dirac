!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine GenHSMRpars
!-----------------------------------------------------------------------
!
!  Generate set of parameters for the HSMRCC run.
!
!  called from InitHSMRBase (source/inihsbase.F90)
!
!  NOTE: Needs NCMO, NVMO, NHOLE, NEMO, NPCTL ... variables must be ready !
!
!  needs functions: IsMainOrb, IsIntmOrb, IsIntmCore, IsIntmVirt ...
!
!-----------------------------------------------------------------------
   use Parameters ! most of parameters are stored here ...
   use OrbMod
   use IOModule
   use Options
   implicit none
   integer i
   integer N_CHOOSE_R

   if ( get_iprint() >= 10 ) then

     call HEADER('Output from GenHSMRpars',-1)

     call print_where

     call print_basic_params

     call print_orbmod_functions

   endif

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

!------------------------------------------------------------------
!    ... Needs: IsMainOrb, IsIntmOrb, IsIntmCore, IsIntmVirt  ...
!------------------------------------------------------------------

!  ... main space, Pm, NMAINORB
   NMAINORB = 0
   do i = 1, NVMO
     if (IsMainOrb(i)) NMAINORB = NMAINORB + 1
   enddo

   if (NMAINORB==0) then
    write(UOUT,'(2x,a)') 'IsMainOrb(i), i=1,NVMO gives all false !'
    call quit('GenHSMRpars: NMAINORB .eq. 0 !')
   endif

!  ... intermediate space, Pi, make ready NINTMORB
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
   if ( get_iprint() >= 7 ) then
     call print_gms_parameters
   endif

end subroutine GenHSMRpars

