module talsh_eom

  use talsh
  use tensor_algebra
  use exacorr_datatypes
  use intermediates
  use talsh_ao_to_mo, only : tensor_norm2
  use exacorr_global, only : debug_exa
  use exacorr_utils,  only : print_date
  use talsh_davidson

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE=(-1.D0,0.D0),ONE=(1.D0,0.D0),ONE_HALF=(0.5D0,0.D0)
  real(8),    parameter :: ONE_QUARTER=0.25D0

  private

  public talsh_eom_driver

  contains

     subroutine talsh_eom_trial_create(workspace, flavor)
        type(talsh_davidson_workspace), intent(inout) :: workspace
        character(len=2), intent(in) :: flavor

        integer :: nocc, nvir, i, j, a, b, ierr, iroot, k(1), N, N1, N2 
        real(8), allocatable :: diagonal(:)
        integer :: i1, j1, a1

        complex(8), pointer, contiguous:: hdiag1_hp_tens(:,:), hdiag2_hp_tens(:,:,:,:)
        complex(8), pointer, contiguous:: hdiag1_h_p_tens(:),  hdiag2_h_p_tens(:,:,:)
        complex(8), pointer, contiguous:: trial1_hp_tens(:,:), trial2_hp_tens(:,:,:,:)
        complex(8), pointer, contiguous:: trial1_h_p_tens(:),  trial2_h_p_tens(:,:,:)
        type(C_PTR):: hdiag1_body_p, hdiag2_body_p
        type(C_PTR):: trial1_body_p, trial2_body_p

        logical :: enable_double_manifold = .true.
        logical :: debug_trial = .true.

        ierr = talsh_tensor_get_body_access(workspace%Hdiag1,hdiag1_body_p,C8,0,DEV_HOST)
        ierr = talsh_tensor_get_body_access(workspace%Hdiag2,hdiag2_body_p,C8,0,DEV_HOST)

! first we copy the diagonal to a local working variable ("diagonal"), in order to sort it later with minloc() 
        select case(flavor)
           case("IP")
              call c_f_pointer(hdiag1_body_p,hdiag1_h_p_tens,workspace%tens1_dims(1:workspace%tens1_rank))
              call c_f_pointer(hdiag2_body_p,hdiag2_h_p_tens,workspace%tens2_dims(1:workspace%tens2_rank))

              if (enable_double_manifold) then
                 allocate(diagonal(workspace%N(1) + workspace%N(2)))
              else
                 allocate(diagonal(workspace%N(1)))
              end if
              diagonal = 0.0d0

              do i = 1, workspace%tens1_dims(1) 
                 diagonal(i) = real(hdiag1_h_p_tens(i))
              end do

! aspg: have to still figure out how to handle the single and double manifolds handled separately
! and find which are the workspace%nroots lowest roots more generally (i.e. when some of the trials
! belong to single and some to double
!
              if (enable_double_manifold) then
                  b = workspace%N(1)
                  do i = 1, workspace%tens2_dims(3)
                     do j = 1, workspace%tens2_dims(2)
                        do a = 1, workspace%tens2_dims(1)
                           b = b + 1
                           diagonal(b) = realpart(hdiag2_h_p_tens(a,j,i))
                        end do
                     end do
                  end do
               end if

              if (debug_trial) then
                 do i = 1, size(diagonal,1)
                     write (6,*) 'H_II(',i,'):',diagonal(i)
                 end do
              end if

        end select

! now set up the initial trial vectors, as unit matrices (or the equivalent for higher-dimensional tensors)

        if (debug_trial) then 
           select case(flavor)
              case("IP")
                 do iroot = 1, workspace%nroots
   
                    ierr = talsh_tensor_get_body_access(workspace%trial1(iroot),trial1_body_p,C8,0,DEV_HOST)
                    ierr = talsh_tensor_get_body_access(workspace%trial2(iroot),trial2_body_p,C8,0,DEV_HOST)
                    call c_f_pointer(trial1_body_p,trial1_h_p_tens,workspace%tens1_dims(1:workspace%tens1_rank))
                    call c_f_pointer(trial2_body_p,trial2_h_p_tens,workspace%tens2_dims(1:workspace%tens2_rank))

                    write (6,*) 'debug >>> t1_h_p_tens, root',iroot,' >>>'
!                   call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
                    do i = 1, workspace%tens1_dims(1) 
                        write (6,*) i, trial1_h_p_tens(i)
                    end do
                    write (6,*) 'debug <<< t1_h_p_tens <<<'

                    write (6,*) 'debug >>> t2_h_p_tens, root',iroot,' >>>'
                    call talsh_tensor_print_info(workspace%trial2(iroot))

!                   call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
                    b = 0
                    do i = 1, workspace%tens2_dims(3)
                       do j = 1, workspace%tens2_dims(2)
                          do a = 1, workspace%tens2_dims(1)
                             b = b + 1
                             write (6,*) b, a, j, i, trial2_h_p_tens(a,j,i)
                          end do
                       end do
                    end do
                    write (6,*) 'debug <<< t2_h_p_tens <<<'
                 end do
           end select

           write (6,*) 'debug >>> real part of hdiag1_h_p_tens >>>'
           write (6,*) diagonal
           write (6,*) 'debug <<< real part of hdiag1_h_p_tens <<<'
 
        end if

! forcing full excitation manifold 
        if (enable_double_manifold) then
              workspace%nroots = workspace%N(1) + workspace%N(2)
        else
! forcing only the single part now
           if (workspace%nroots.gt.workspace%N(1)) then
              workspace%nroots = workspace%N(1)
              write(6,*) 'WARNING!! resetting number of roots to :',workspace%N(1)
           end if 
        end if

        write(6,*) 'roots',workspace%nroots
        j = 0
        do iroot = 1, workspace%nroots

!          k = minloc(diagonal, MASK=(diag_H .ge. in_window_HbarSS_threshold))
           if (size(diagonal).gt.2) then
              k = minloc(diagonal)
              i = k(1)
              if (debug_trial) write(6,*) 'debug, roots',iroot,diagonal(i),i
              diagonal(i) = huge(diagonal)
              if (debug_trial) write(6,*) 'debug, making diagonal huge',diagonal(i)
           else 
              j = j + 1
              i = j
           end if

!! we may have diagonal elements which do not correspond to an element within the windows, so we check that
!! this is the case
!          if (index_in_window(i).eqv..true.) then
           if (.true.) then
              select case(flavor)
                 case("IP")

                    ierr = talsh_tensor_get_body_access(workspace%trial1(iroot),trial1_body_p,C8,0,DEV_HOST)
                    ierr = talsh_tensor_get_body_access(workspace%trial2(iroot),trial2_body_p,C8,0,DEV_HOST)
                    call c_f_pointer(trial1_body_p,trial1_h_p_tens,workspace%tens1_dims(1:workspace%tens1_rank))
                    call c_f_pointer(trial2_body_p,trial2_h_p_tens,workspace%tens2_dims(1:workspace%tens2_rank))

! aspg, come back and see how to handle this
                    if (enable_double_manifold) then ! we have to go into the 2h1p section to select trial vectors
                       if (j .le. workspace%N(1)) then 
                          trial1_h_p_tens(j) = ONE
                       else
! translate j into the spinor indexes:w
                          print *,'j:',j,i1,j1,a1
                          call from_linear_to_indexes_oov(j,i1,j1,a1,workspace%tens2_dims(3),workspace%tens2_dims(1))
                          trial2_h_p_tens(a1,j1,i1) = ONE
                       end if
                    else
! aspg, for now, we only set the tens1 part of the trial vectors
                       trial1_h_p_tens(j) = ONE
                    end if
              end select
           end if
        end do

        deallocate (diagonal)

! debug information
        if (debug_trial) then
           select case(flavor)
              case("IP")
                 do iroot = 1, workspace%nroots
   
                    ierr = talsh_tensor_get_body_access(workspace%trial1(iroot),trial1_body_p,C8,0,DEV_HOST)
                    ierr = talsh_tensor_get_body_access(workspace%trial2(iroot),trial2_body_p,C8,0,DEV_HOST)
                    call c_f_pointer(trial1_body_p,trial1_h_p_tens,workspace%tens1_dims(1:workspace%tens1_rank))
                    call c_f_pointer(trial2_body_p,trial2_h_p_tens,workspace%tens2_dims(1:workspace%tens2_rank))
   
                    write (6,*) 'debug >>> trial1_h_p_tens, root',iroot,' >>>'
!                   call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
                    do i = 1, workspace%tens1_dims(1)
                        write (6,*) i, trial1_h_p_tens(i)
                    end do
                    write (6,*) 'debug <<< trial1_h_p_tens <<<'
   
                    write (6,*) 'debug >>> trial2_h_p_tens, root',iroot,' >>>'
!                   call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
                    b = 0
                    do i = 1, workspace%tens2_dims(3)
                       do j = 1, workspace%tens2_dims(2)
                          do a = 1, workspace%tens2_dims(1)
                             b = b + 1
                             write (6,*) b, a, j, i, trial2_h_p_tens(a,j,i)
                       end do
                       end do
                    end do
                 end do
           end select
        end if

     end subroutine

    subroutine talsh_eom_Hdiag_create(workspace, Fmi_tensor, Fae_tensor, flavor)
      type(talsh_davidson_workspace), intent(inout) :: workspace
      type(talsh_tens_t), intent(inout) :: Fmi_tensor, Fae_tensor
      character(len=2), intent(in) :: flavor
      integer :: nocc, nvir, i, j, a, b, ierr 
      integer :: Fmi_dims(2), Fae_dims(2), Fmi_rank, Fae_rank
      real(8) :: diagonal
      complex(8), pointer, contiguous:: hdiag1_hp_tens(:,:), hdiag2_hp_tens(:,:,:,:)
      complex(8), pointer, contiguous:: hdiag1_h_p_tens(:),  hdiag2_h_p_tens(:,:,:)
      complex(8), pointer, contiguous:: Fmi_p_tens(:,:), Fae_p_tens(:,:) 
      type(C_PTR):: hdiag1_body_p, hdiag2_body_p
      type(C_PTR):: Fmi_body_p, Fae_body_p
      logical :: debug_hdiag_create = .false.


      diagonal = ZERO

      ierr = talsh_tensor_dimensions(Fmi_tensor,Fmi_rank,Fmi_dims)
      ierr = talsh_tensor_get_body_access(Fmi_tensor,Fmi_body_p,C8,0,DEV_HOST)
      call c_f_pointer(Fmi_body_p,fmi_p_tens,Fmi_dims(1:2))

      ierr = talsh_tensor_dimensions(Fae_tensor,Fae_rank,Fae_dims)
      ierr = talsh_tensor_get_body_access(Fae_tensor,Fae_body_p,C8,0,DEV_HOST)
      call c_f_pointer(Fae_body_p,fae_p_tens,Fae_dims(1:2))

      ierr = talsh_tensor_get_body_access(workspace%Hdiag1,hdiag1_body_p,C8,0,DEV_HOST)
      ierr = talsh_tensor_get_body_access(workspace%Hdiag2,hdiag2_body_p,C8,0,DEV_HOST)

! transfer the diagonals of 

      select case(flavor)
         case("IP")
            call c_f_pointer(hdiag1_body_p,hdiag1_h_p_tens,workspace%tens1_dims(1:workspace%tens1_rank))
            call c_f_pointer(hdiag2_body_p,hdiag2_h_p_tens,workspace%tens2_dims(1:workspace%tens2_rank))

            do i = 1,  workspace%tens1_dims(1)
               diagonal = -fmi_p_tens(i,i)
               if (debug_hdiag_create) write (6,*) 'f_ii:',fmi_p_tens(i,i),'diag:',diagonal
               hdiag1_h_p_tens(i) = diagonal
            end do
            if (debug_hdiag_create) then
               write (6,*) 'debug >>> hdiag1_h_p_tens >>>'
               write (6,*) hdiag1_h_p_tens(1:nocc)
               write (6,*) 'debug <<< hdiag1_h_p_tens <<<'
            end if

            do i = 1, workspace%tens2_dims(3)
               do j = 1, workspace%tens2_dims(2)
                  do a = 1, workspace%tens2_dims(1)
                     diagonal = fae_p_tens(a,a) - fmi_p_tens(j,j) - fmi_p_tens(i,i)
                     if (debug_hdiag_create) &
                     write (6,*) 'f_ii:',fmi_p_tens(i,i),'f_jj:',fmi_p_tens(j,j),'f_aa',fae_p_tens(a,a),'diag:',diagonal
                     hdiag2_h_p_tens(a,j,i) = diagonal 
                  end do
               end do
            end do
            if (debug_hdiag_create) then
               write (6,*) 'debug >>> hdiag2_h_p_tens >>>'
               write (6,*) hdiag2_h_p_tens(1:nvir,1:nocc,1:nocc)
!              write (6,*) 'debug <<< hdiag2_h_p_tens <<<'
            end if

         case("EA")
            call c_f_pointer(hdiag1_body_p,hdiag1_h_p_tens,workspace%tens1_dims(1:workspace%tens1_rank))
            call c_f_pointer(hdiag2_body_p,hdiag2_h_p_tens,workspace%tens2_dims(1:workspace%tens2_rank))

            do a = 1, workspace%tens1_dims(1)
               diagonal = fae_p_tens(a,a)
               hdiag1_h_p_tens(a) = diagonal
            end do

            do i = 1, workspace%tens2_dims(3)
               do b = 1, workspace%tens2_dims(2)
                  do a = 1, workspace%tens2_dims(1)
                     diagonal = fae_p_tens(a,a) + fae_p_tens(b,b) - fmi_p_tens(i,i) 
                     hdiag2_h_p_tens(a,b,i) = diagonal
                  end do
               end do
            end do

         case("EE")
            call c_f_pointer(hdiag1_body_p,hdiag1_hp_tens,workspace%tens1_dims(1:workspace%tens1_rank))
            call c_f_pointer(hdiag2_body_p,hdiag2_hp_tens,workspace%tens2_dims(1:workspace%tens2_rank))

            do i = 1, workspace%tens1_dims(2)
               do a = 1, workspace%tens1_dims(1)
                  diagonal = fae_p_tens(a,a) - fmi_p_tens(i,i)
                  hdiag1_hp_tens(a,i) = diagonal
               end do
            end do

            do j = 1,  workspace%tens2_dims(4)
               do i = 1,  workspace%tens2_dims(3)
                  do b = 1,  workspace%tens2_dims(2)
                     do a = 1,   workspace%tens2_dims(1)
                        diagonal = fae_p_tens(a,a) + fae_p_tens(b,b) - fmi_p_tens(j,j) - fmi_p_tens(i,i)
                        hdiag2_hp_tens(a,b,i,j) = diagonal
                     end do
                  end do
               end do
            end do

       end select

    end subroutine

    subroutine talsh_eom_driver (exa_input,eps_occ,eps_vir, t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor,t1_tensor)

!     input variables
      type(exacc_input),     intent(in   )           :: exa_input
      type(talsh_tens_t),    intent(inout)           :: t2_tensor
      type(talsh_intg_tens), intent(inout)           :: int_t
      type(talsh_tens_t),    intent(inout)           :: foo_tensor, fov_tensor, fvv_tensor
      type(talsh_tens_t),    intent(inout), optional :: t1_tensor, l2_tensor
      real(8),               intent(in)              :: eps_occ(:), eps_vir(:)
!     eom tensors
      type(talsh_eom_tens) :: eom_t
      type(talsh_davidson_workspace) :: workspace_t
      type(talsh_projectors) :: projector_t


      write(*,*) ""
      write(*,'(A)') "*******************************************************************************"
      write(*,'(A)') "***************************** Equation of Motion ******************************"
      write(*,'(A)') "*******************************************************************************"
      write(*,*) ""

      call talsh_davidson_setup(unit_output=6, &
                                solve_right=.true.,&
                                solve_left=.false.,&
                                symmetric=.false., &
                                complex_mode=.true.)

      call talsh_davidson_workspace_create(workspace_t, exa_input%nocc, exa_input%nvir, exa_input%eom_flavour, 2)
!     call talsh_tensor_print_info(workspace_t%Hdiag1)
!     call talsh_tensor_print_info(workspace_t%Hdiag2)
!     call talsh_eom_Hdiag_create(workspace_t, foo_tensor, fvv_tensor, exa_input%eom_flavour)
!     call talsh_tensor_print_info(workspace_t%Hdiag1)
!     call talsh_tensor_print_info(workspace_t%Hdiag2)
!     call talsh_eom_trial_create(workspace_t, exa_input%eom_flavour)

!************************************
!     EOM calculation: IP, EA or EE *
!************************************

      if (.not.exa_input%ccd) then
        select case(exa_input%eom_flavour)
          case("IP")
             print*, "foo_tensor norm", tensor_norm2(foo_tensor)
             print*, "fov_tensor norm", tensor_norm2(fov_tensor)
             print*, "fvv_tensor norm", tensor_norm2(fvv_tensor)

            call eom_ip(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
          case("EA")
            call eom_ea(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
          case("EE")
            if (.not.exa_input%eom_properties) then
        call eom_ee(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor=t1_tensor)
            else
        call eom_ee(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor,t1_tensor)
            end if
          case default
            call quit ('Error: No valid EOM-flavour chosen')
        end select
      else
        select case(exa_input%eom_flavour)
          case("IP")
            call eom_ip(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
          case("EA")
            call eom_ea(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
          case("EE")
            if (.not.exa_input%eom_properties) then
              call eom_ee(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
            else
              call eom_ee(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor)
            end if
        case default
          call quit ('Error: No valid EOM-flavour chosen')
        end select
      end if

      call talsh_davidson_workspace_destroy(workspace_t)

!     make some noise so that we know we are leaving
      print*, ""
      call print_date("Leaving eom_driver routine")
      return

    end subroutine talsh_eom_driver

    subroutine eom_ip(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
       use talsh_eom_sigma, only : sigma_right_ip, sigma_left_ip

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_davidson_workspace), intent(inout) :: workspace_t
      type(talsh_projectors), intent(inout) :: projector_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_tens_t),    intent(inout) :: foo_tensor, fov_tensor, fvv_tensor
      type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!     tensor dims
      integer :: nocc, nvir
      integer, dimension(2) :: oo_dims, ov_dims, vv_dims
      integer, dimension(3) :: voo_dims, oov_dims
      integer, dimension(4) :: oooo_dims, ooov_dims, ovoo_dims, vovo_dims
!     CCD switch
      logical :: CCD
!     error code
      integer :: ierr

      CCD = exa_input%ccd

!************************
!     construct tensors *
!************************

      nocc = exa_input%nocc
      nvir = exa_input%nvir

!     set tensor dims
      oo_dims = (/nocc,nocc/)
      ov_dims = (/nocc,nvir/)
      vv_dims = (/nvir,nvir/)
      voo_dims = (/nvir,nocc,nocc/)
      oooo_dims = (/nocc,nocc,nocc,nocc/)
      ooov_dims = (/nocc,nocc,nocc,nvir/)
      ovoo_dims = (/nocc,nvir,nocc,nocc/)
      vovo_dims = (/nvir,nocc,nvir,nocc/)
      if (exa_input%eom_properties) oov_dims = (/nocc,nocc,nvir/) ! nocc 1h-states, nocc*nocc*nvir 2h1p states

!     create all fixed intermediates
      ierr=talsh_tensor_construct(eom_t%fbar_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_ov,C8,ov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_vv,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_oooo,C8,oooo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_ooov,C8,ooov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_ovoo,C8,ovoo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vovo,C8,vovo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)

      call get_wbar_vovo(eom_t%wbar_vovo,t2_tensor,int_t)
      if (.not.CCD) then
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv,t1_tensor)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_w_oooo (eom_t%w_oooo,int_t,t2_tensor,t1_tensor)
        call get_w_ooov (eom_t%w_ooov,int_t,t1_tensor)
        call get_w_ovoo (eom_t%w_ovoo,int_t,eom_t%fbar_ov,eom_t%w_oooo,eom_t%wbar_vovo,t2_tensor,t1_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor,t1_tensor)
      else
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor)
        call get_w_oooo (eom_t%w_oooo,int_t,t2_tensor)
        call get_w_ooov (eom_t%w_ooov,int_t)
        call get_w_ovoo (eom_t%w_ovoo,int_t,eom_t%fbar_ov,eom_t%w_oooo,eom_t%wbar_vovo,t2_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor)
      end if

      print*, "t1        tensor norm", tensor_norm2(t1_tensor)
      print*, "t2        tensor norm", tensor_norm2(t2_tensor)
      print*, "fbar_oo   tensor norm", tensor_norm2(eom_t%fbar_oo)
      print*, "fbar_ov   tensor norm", tensor_norm2(eom_t%fbar_ov)
      print*, "fbar_vv   tensor norm", tensor_norm2(eom_t%fbar_vv)
      print*, "w_oooo    tensor norm", tensor_norm2(eom_t%w_oooo)
      print*, "w_ooov    tensor norm", tensor_norm2(eom_t%w_ooov)
      print*, "w_ovoo    tensor norm", tensor_norm2(eom_t%w_ovoo)
      print*, "w_vooo    tensor norm", tensor_norm2(eom_t%w_vovo)
      print*, "wbar_vovo tensor norm", tensor_norm2(eom_t%wbar_vovo)

      ierr=talsh_tensor_destruct(eom_t%wbar_vovo)

      call talsh_tensor_print_info(workspace_t%Hdiag1)
      call talsh_tensor_print_info(workspace_t%Hdiag2)
      call talsh_eom_Hdiag_create(workspace_t, eom_t%Fbar_oo, eom_t%Fbar_vv, exa_input%eom_flavour)
      call talsh_tensor_print_info(workspace_t%Hdiag1)
      call talsh_tensor_print_info(workspace_t%Hdiag2)
      call talsh_eom_trial_create(workspace_t, exa_input%eom_flavour)


!***********************************
!     calculate right sigma vector *
!***********************************

      call talsh_davidson_driver(exa_input,int_t,eom_t,t2_tensor, workspace_t, projector_t)
      call print_date("Calculated EOM-IP RHS")

!************************************************
!     calculate left eigenvector for properties *
!************************************************

!      if (exa_input%eom_properties) then
!        call sigma_left_ip(exa_input,int_t,eom_t,t2_tensor)
!        call print_date("Calculated EOM-IP LHS")
!      end if

    end subroutine eom_ip

    subroutine eom_ea(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
       use talsh_eom_sigma, only : sigma_right_ea, sigma_left_ea

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_davidson_workspace), intent(inout) :: workspace_t
      type(talsh_projectors), intent(inout) :: projector_t
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_tens_t),    intent(inout) :: foo_tensor, fov_tensor, fvv_tensor
      type(talsh_tens_t),    intent(inout), optional :: t1_tensor
!     tensor dims
      integer :: nocc, nvir
      integer, dimension(2) :: oo_dims, ov_dims, vv_dims
      integer, dimension(3) :: vvo_dims, ovv_dims
      integer, dimension(4) :: vovo_dims, vovv_dims, vvvo_dims, vvvv_dims
!     CCD switch
      logical :: CCD
!     error code
      integer :: ierr

      CCD = exa_input%ccd

!************************
!     construct tensors *
!************************

      nocc = exa_input%nocc
      nvir = exa_input%nvir

!     set tensor dims
      oo_dims = (/nocc,nocc/)
      ov_dims = (/nocc,nvir/)
      vv_dims = (/nvir,nvir/)
      vvo_dims = (/nvir,nvir,nocc/)
      vovo_dims = (/nvir,nocc,nvir,nocc/)
      vovv_dims = (/nvir,nocc,nvir,nvir/)
      vvvo_dims = (/nvir,nvir,nvir,nocc/)
      vvvv_dims = (/nvir,nvir,nvir,nvir/)
      if (exa_input%eom_properties) ovv_dims = (/nocc,nvir,nvir/)

!     create all fixed intermediates
      ierr=talsh_tensor_construct(eom_t%fbar_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_ov,C8,ov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_vv,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vovo,C8,vovo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vovv,C8,vovv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vvvo,C8,vvvo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vvvv,C8,vvvv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)

      call get_wbar_vovo(eom_t%wbar_vovo,t2_tensor,int_t)
      if (.not.CCD) then
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv,t1_tensor)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor,t1_tensor)
        call get_w_vovv (eom_t%w_vovv,int_t,t1_tensor)
        call get_w_vvvv (eom_t%w_vvvv,int_t,t2_tensor,t1_tensor)
        call get_w_vvvo (eom_t%w_vvvo,int_t,eom_t%w_vvvv,eom_t%wbar_vovo,eom_t%fbar_ov,t2_tensor,t1_tensor)
      else
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor)
        call get_w_vovv (eom_t%w_vovv,int_t)
        call get_w_vvvv (eom_t%w_vvvv,int_t,t2_tensor)
        call get_w_vvvo (eom_t%w_vvvo,int_t,eom_t%w_vvvv,eom_t%wbar_vovo,eom_t%fbar_ov,t2_tensor)
      end if
      ierr=talsh_tensor_destruct(eom_t%wbar_vovo)

!***********************************
!     calculate right sigma vector *
!***********************************

      call talsh_tensor_print_info(workspace_t%Hdiag1)
      call talsh_tensor_print_info(workspace_t%Hdiag2)
      call talsh_eom_Hdiag_create(workspace_t, eom_t%Fbar_oo, eom_t%Fbar_vv, exa_input%eom_flavour)
      call talsh_tensor_print_info(workspace_t%Hdiag1)
      call talsh_tensor_print_info(workspace_t%Hdiag2)
      call talsh_eom_trial_create(workspace_t, exa_input%eom_flavour)

      call talsh_davidson_driver(exa_input,int_t,eom_t,t2_tensor, workspace_t, projector_t)
      call print_date("Calculated EOM-EA RHS")

!************************************************
!     calculate left eigenvector for properties *
!************************************************

!      if (exa_input%eom_properties) then
!        call sigma_left_ea(exa_input,int_t,eom_t,t2_tensor)
!        call print_date("Calculated EOM-EA LHS")
!      end if

!**************
!     cleanup *
!**************

      ierr=talsh_tensor_destruct(eom_t%fbar_oo)
      ierr=talsh_tensor_destruct(eom_t%fbar_ov)
      ierr=talsh_tensor_destruct(eom_t%fbar_vv)
      ierr=talsh_tensor_destruct(eom_t%w_vovo)
      ierr=talsh_tensor_destruct(eom_t%w_vovv)
      ierr=talsh_tensor_destruct(eom_t%w_vvvo)
      ierr=talsh_tensor_destruct(eom_t%w_vvvv)

    end subroutine eom_ea

    subroutine eom_ee(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor,t1_tensor)
       use talsh_eom_sigma, only : sigma_right_ee, sigma_left_ee

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
      type(talsh_davidson_workspace), intent(inout) :: workspace_t
      type(talsh_projectors), intent(inout) :: projector_t
      type(talsh_eom_tens),  intent(inout) :: eom_t
      type(talsh_tens_t),    intent(inout) :: t2_tensor
      type(talsh_intg_tens), intent(inout) :: int_t
      type(talsh_tens_t),    intent(inout) :: foo_tensor, fov_tensor, fvv_tensor
      type(talsh_tens_t),    intent(inout), optional :: l2_tensor,t1_tensor
!     tensor dims
      integer :: nocc, nvir
      integer, dimension(2) :: oo_dims, ov_dims, vo_dims, vv_dims
      integer, dimension(4) :: oooo_dims, ooov_dims, ovoo_dims, oovv_dims, vvoo_dims, vovo_dims, vovv_dims, vvvo_dims, vvvv_dims
!     CCD switch
      logical :: CCD
!     error code
      integer :: ierr

      CCD = exa_input%ccd

!************************
!     construct tensors *
!************************

      nocc = exa_input%nocc
      nvir = exa_input%nvir

!     set tensor dims
      oo_dims = (/nocc,nocc/)
      ov_dims = (/nocc,nvir/)
      vo_dims = (/nvir,nocc/)
      vv_dims = (/nvir,nvir/)
      oooo_dims = (/nocc,nocc,nocc,nocc/)
      ooov_dims = (/nocc,nocc,nocc,nvir/)
      ovoo_dims = (/nocc,nvir,nocc,nocc/)
      if (exa_input%eom_properties) oovv_dims = (/nocc,nocc,nvir,nvir/)
      vvoo_dims = (/nvir,nvir,nocc,nocc/)
      vovo_dims = (/nvir,nocc,nvir,nocc/)
      vovv_dims = (/nvir,nocc,nvir,nvir/)
      vvvo_dims = (/nvir,nvir,nvir,nocc/)
      vvvv_dims = (/nvir,nvir,nvir,nvir/)

!     create all fixed intermediates
      ierr=talsh_tensor_construct(eom_t%fbar_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_ov,C8,ov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%fbar_vv,C8,vv_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%goo,C8,oo_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(eom_t%gvv,C8,vv_dims,init_val=ZERO)
      end if
      ierr=talsh_tensor_construct(eom_t%w_oooo,C8,oooo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_ooov,C8,ooov_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_ovoo,C8,ovoo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vovo,C8,vovo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vovv,C8,vovv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vvvo,C8,vvvo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%w_vvvv,C8,vvvv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(eom_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)

      call get_wbar_vovo(eom_t%wbar_vovo,t2_tensor,int_t)
      if (.not.CCD) then
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv,t1_tensor)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor,t1_tensor)
        call get_w_oooo (eom_t%w_oooo,int_t,t2_tensor,t1_tensor)
        call get_w_ooov (eom_t%w_ooov,int_t,t1_tensor)
        call get_w_ovoo (eom_t%w_ovoo,int_t,eom_t%fbar_ov,eom_t%w_oooo,eom_t%wbar_vovo,t2_tensor,t1_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor,t1_tensor)
        call get_w_vovv (eom_t%w_vovv,int_t,t1_tensor)
        call get_w_vvvv (eom_t%w_vvvv,int_t,t2_tensor,t1_tensor)
        call get_w_vvvo (eom_t%w_vvvo,int_t,eom_t%w_vvvv,eom_t%wbar_vovo,eom_t%fbar_ov,t2_tensor,t1_tensor)
      else
        call get_fbar_oo(eom_t%fbar_oo,int_t,foo_tensor,fov_tensor,t2_tensor)
        call get_fbar_ov(eom_t%fbar_ov,fov_tensor,CCD,int_t%oovv)
        call get_fbar_vv(eom_t%fbar_vv,int_t,fvv_tensor,fov_tensor,t2_tensor)
        call get_w_oooo (eom_t%w_oooo,int_t,t2_tensor)
        call get_w_ooov (eom_t%w_ooov,int_t)
        call get_w_ovoo (eom_t%w_ovoo,int_t,eom_t%fbar_ov,eom_t%w_oooo,eom_t%wbar_vovo,t2_tensor)
        call get_w_vovo (eom_t%w_vovo,int_t,t2_tensor)
        call get_w_vovv (eom_t%w_vovv,int_t)
        call get_w_vvvv (eom_t%w_vvvv,int_t,t2_tensor)
        call get_w_vvvo (eom_t%w_vvvo,int_t,eom_t%w_vvvv,eom_t%wbar_vovo,eom_t%fbar_ov,t2_tensor)
      end if
      ierr=talsh_tensor_destruct(eom_t%wbar_vovo)

!***********************************
!     calculate right sigma vector *
!***********************************

      call talsh_tensor_print_info(workspace_t%Hdiag1)
      call talsh_tensor_print_info(workspace_t%Hdiag2)
      call talsh_eom_Hdiag_create(workspace_t, eom_t%Fbar_oo, eom_t%Fbar_vv, exa_input%eom_flavour)
      call talsh_tensor_print_info(workspace_t%Hdiag1)
      call talsh_tensor_print_info(workspace_t%Hdiag2)
      call talsh_eom_trial_create(workspace_t, exa_input%eom_flavour)

      call talsh_davidson_driver(exa_input,int_t,eom_t,t2_tensor, workspace_t, projector_t)
      call print_date("Calculated EOM-EE RHS")

!************************************************
!     calculate left eigenvector for properties *
!************************************************

!     create mixed type intermediates
      if (exa_input%eom_properties) then
        call get_goo_lambda(eom_t%goo,l2_tensor,t2_tensor)
        call get_gvv_lambda(eom_t%gvv,l2_tensor,t2_tensor)
      end if

!     if (exa_input%eom_properties) then
!       if (.not.CCD) then
!         call sigma_left_ee(exa_input,int_t,eom_t,t1_tensor)
!       else
!         call sigma_left_ee(exa_input,int_t,eom_t)
!       end if
!       call print_date("Calculated EOM-EE LHS")
!     end if

!**************
!     cleanup *
!**************

      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(eom_t%goo)
        ierr=talsh_tensor_destruct(eom_t%gvv)
      end if
      ierr=talsh_tensor_destruct(eom_t%fbar_oo)
      ierr=talsh_tensor_destruct(eom_t%fbar_ov)
      ierr=talsh_tensor_destruct(eom_t%fbar_vv)
      ierr=talsh_tensor_destruct(eom_t%w_oooo)
      ierr=talsh_tensor_destruct(eom_t%w_ooov)
      ierr=talsh_tensor_destruct(eom_t%w_ovoo)
      ierr=talsh_tensor_destruct(eom_t%w_vovo)
      ierr=talsh_tensor_destruct(eom_t%w_vovv)
      ierr=talsh_tensor_destruct(eom_t%w_vvvo)
      ierr=talsh_tensor_destruct(eom_t%w_vvvv)

    end subroutine eom_ee

    subroutine from_linear_to_indexes_oov(b_desired,i,j,a,dim_o,dim_v)
       integer, intent(in)  :: b_desired, dim_o, dim_v
       integer, intent(out) :: i, j, a
       integer :: b
       
       b = 0
       do i = 1, dim_o
          do j = 1, dim_o
             do a = 1, dim_v
                b = b + 1
                print *, 'i:',i,'j:',j,'a:',a,'   b:',b
                if  (b .eq. b_desired) return
             end do
          end do
       end do
    end subroutine



end module
