

module fde_types

   use fde_mag_cfg

   implicit none

   public

   type fde_grid
      integer :: id
      integer :: npoints
      real(kind=8), pointer :: r(:,:)    => null()
      real(kind=8), pointer :: w(:)      => null()
   end type

   type grid_function 
      integer :: id
      integer :: npoints
      integer :: nspin
      real(kind=8), pointer :: n(:)          => null()
      real(kind=8), pointer :: n_b(:,:)      => null()
      real(kind=8), pointer :: s_b(:,:,:)    => null()
      real(kind=8), pointer :: gn(:,:)       => null()
      real(kind=8), pointer :: gn_b(:,:,:)   => null()
      real(kind=8), pointer :: gs_b(:,:,:,:) => null()
      real(kind=8), pointer :: elpot(:)      => null()
      real(kind=8), pointer :: nucpot(:)     => null()
      real(kind=8), pointer :: hn(:,:)       => null()
   end type

   public new_fde_grid
   public del_fde_grid
!gosia: unused   public set_fde_grid_size

   public set_fde_grid_nr_elements


   public new_grid_function
   public del_grid_function
   public set_grid_function_nspin

   contains


      subroutine new_grid_function(gf,nr_elements)
         type(grid_function) :: gf
         integer :: nr_elements, prior_nr_elements

         gf%npoints = nr_elements
         gf%nspin   = 1
         if (.not.associated(gf%n)) then 
            allocate(gf%n(nr_elements))
         else
            prior_nr_elements = size(gf%n)
            print *, 'Warning ! gf%n already allocated !'
            if (prior_nr_elements .ne. nr_elements) then
               print *, prior_nr_elements,' elements then,',nr_elements,'now !'
               stop
            endif 
         endif

         if (.not.associated(gf%gn)) then
            allocate(gf%gn(3,nr_elements))
         else
            prior_nr_elements = size(gf%gn,2)
            print *, 'Warning ! gf%gn already allocated !'
            if (prior_nr_elements .ne. nr_elements) then
               print *, prior_nr_elements,' elements then,',nr_elements,'now !'
               stop
            endif
         endif

         if (.not.associated(gf%hn)) then
            allocate(gf%hn(6,nr_elements))
         else
            prior_nr_elements = size(gf%hn,2)
            print *, 'Warning ! gf%hn already allocated !'
            if (prior_nr_elements .ne. nr_elements) then
               print *, prior_nr_elements,' elements then,',nr_elements,'now !'
               stop
            endif
         endif

         if (.not.associated(gf%elpot)) then
            allocate(gf%elpot(nr_elements))
         else
            prior_nr_elements = size(gf%elpot)
            print *, 'Warning ! gf%elpot already allocated !'
            if (prior_nr_elements .ne. nr_elements) then
               print *, prior_nr_elements,' elements then,',nr_elements,'now !'
               stop
            endif
         endif

         if (.not.associated(gf%nucpot)) then
            allocate(gf%nucpot(nr_elements))
         else
            prior_nr_elements = size(gf%nucpot)
            print *, 'Warning ! gf%nucpot already allocated !'
            if (prior_nr_elements .ne. nr_elements) then
               print *, prior_nr_elements,' elements then,',nr_elements,'now !'
               stop
            endif
         endif

         if (fde_rsp_mag_lao_import .or. fde_rsp_mag_lao_export) then
            if (.not. associated(gf%n_b)) then
               allocate(gf%n_b(3,nr_elements))
            else
               prior_nr_elements = size(gf%n_b)
               print *, 'Warning ! gf%n_b already allocated !'
               if (prior_nr_elements .ne. nr_elements) then
                  print *, prior_nr_elements,' elements then,',nr_elements,'now !'
                  stop
               endif
            end if
            if (.not. associated(gf%gn_b)) then
               allocate(gf%gn_b(3,3,nr_elements))
            else
               prior_nr_elements = size(gf%gn_b)
               print *, 'Warning ! gf%gn_b already allocated !'
               if (prior_nr_elements .ne. nr_elements) then
                  print *, prior_nr_elements,' elements then,',nr_elements,'now !'
                  stop
               endif
            end if
            if (.not. fde_rsp_mag_lao_no_s) then
               if (.not. associated(gf%s_b)) then
                  allocate(gf%s_b(3,3,nr_elements))
               else
                  prior_nr_elements = size(gf%s_b)
                  print *, 'Warning ! gf%s_b already allocated !'
                  if (prior_nr_elements .ne. nr_elements) then
                     print *, prior_nr_elements,' elements then,',nr_elements,'now !'
                     stop
                  endif
               end if
               if (.not. associated(gf%gs_b)) then
                  allocate(gf%gs_b(3,3,3,nr_elements))
               else
                  prior_nr_elements = size(gf%gs_b)
                  print *, 'Warning ! gf%gs_b already allocated !'
                  if (prior_nr_elements .ne. nr_elements) then
                     print *, prior_nr_elements,' elements then,',nr_elements,'now !'
                     stop
                  endif
               end if
            end if
         end if
         
         gf%n  = 0.0d0
         gf%gn = 0.0d0
         gf%hn = 0.0d0
         gf%elpot = 0.0d0
         gf%nucpot = 0.0d0

         if (fde_rsp_mag_lao_import .or. fde_rsp_mag_lao_export) then
            gf%n_b  = 0.0d0
            gf%gn_b = 0.0d0
            if (.not. fde_rsp_mag_lao_no_s) then
               gf%s_b  = 0.0d0
               gf%gs_b = 0.0d0
            end if
         end if

      end subroutine new_grid_function


      subroutine del_grid_function(gf)
         type(grid_function) :: gf
         integer :: nr_elements

         gf%npoints = 0
         if (associated(gf%n)) &
            deallocate(gf%n)

         if (associated(gf%gn)) &
            deallocate(gf%gn)

         if (associated(gf%hn)) &
            deallocate(gf%hn)

         if (associated(gf%elpot)) &
            deallocate(gf%elpot)

         if (associated(gf%nucpot)) &
            deallocate(gf%nucpot)

         if (associated(gf%n_b)) &
            deallocate(gf%n_b)

         if (associated(gf%gn_b)) &
            deallocate(gf%gn_b)

         if (associated(gf%s_b)) &
            deallocate(gf%s_b)

         if (associated(gf%gs_b)) &
            deallocate(gf%gs_b)

      end subroutine del_grid_function


      subroutine set_grid_function_nspin(gf,nspin)
         type(grid_function) :: gf
         integer :: nspin 

         gf%nspin = nspin
      end subroutine set_grid_function_nspin
      

     subroutine new_fde_grid(g,nr_elements)
       type(fde_grid) :: g
       integer :: nr_elements

       g%npoints = nr_elements

       if (.not.associated(g%r)) then
          allocate(g%r(3,nr_elements))
       else
          prior_nr_elements = size(g%r,2)
          print *, 'Warning ! g%r already allocated !'
          if (prior_nr_elements .ne. nr_elements) then
             print *, prior_nr_elements,' elements then,',nr_elements,'now !'
             stop
          endif
       endif

       if (.not.associated(g%w)) then
          allocate(g%w(nr_elements))
       else
          prior_nr_elements = size(g%w)
          print *, 'Warning ! gf%nucpot already allocated !'
          if (prior_nr_elements .ne. nr_elements) then
             print *, prior_nr_elements,' elements then,',nr_elements,'now !'
             stop
          endif
       endif

     end subroutine new_fde_grid


     subroutine del_fde_grid(g)
        type(fde_grid) :: g

        g%npoints = 0
        if (associated(g%r)) &
           deallocate(g%r)

        if (associated(g%w)) &
           deallocate(g%w)

     end subroutine del_fde_grid

end
