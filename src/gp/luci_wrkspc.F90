module luci_wrkspc
    
  implicit none
  real(8),allocatable ::  work(:)

contains

  subroutine allocate_wrkspc
     use memory_allocator
     integer :: lsize
     call legacy_lwork_get(lsize)
     if (.not.allocated(work)) allocate(work(lsize))
!    print*,"Allocated ",lsize," real words for use by the LUCI modules"
  end


  subroutine delete_wrkspc
     deallocate(work)
!    print*,"Deallocated memory for the LUCI modules"
  end

end module
