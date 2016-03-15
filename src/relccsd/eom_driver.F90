  module eom_driver

!this module calculates the sigma vectors for IP(2h-1p), EA(1h-2p) and EE(2h-2p)  

  use contraction
  use intermediates_1b_2b

!    public eom_ip  
!    public eom_ea
    public eom_ee

    private

   contains

   subroutine eom_driver()


    call getampt(t1,t2)

    call IM_container()
   
    call davidson 


   end subroutine


   end module 
