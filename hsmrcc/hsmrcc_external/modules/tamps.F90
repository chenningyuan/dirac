module TAmps

   use Parameters
   use IOModule
   implicit none
   double precision, allocatable         :: TSING(:)
   double precision, allocatable, target :: TOLD(:), TNEW(:)

contains

subroutine AllocTAMP

   implicit none

   if (.not.allocated(TSING)) allocate(TSING(1:NSING))
   if (.not.allocated(TOLD))  allocate(TOLD(0:NTSIZE))
   if (.not.allocated(TNEW))  allocate(TNEW(0:NTSIZE))

   if (get_iprint() > 8) then
     write(UOUT,*) '    TSING allocated with size=',NSING
     write(UOUT,*) 'TOLD,TNEW allocated with size=',NTSIZE
   endif

   return
end subroutine AllocTAMP

subroutine DeAllocTAMP

   implicit none

   if (allocated(TSING)) deallocate(TSING)
   if (allocated(TOLD)) deallocate(TOLD)
   if (allocated(TNEW)) deallocate(TNEW)

   return

end subroutine DeAllocTAMP

subroutine SaveTAMP

   implicit none
   integer :: vac

   do vac = 1, NSTATE
      read (UTNEW,rec=vac) TNEW
      write(UTAMP,rec=vac) TNEW
   enddo

   return

end subroutine SaveTAMP

end module TAmps
