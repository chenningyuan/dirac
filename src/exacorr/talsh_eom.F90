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

        integer :: nocc, nvir, i, j, a, b, ierr, iroot, k(1) 
        real(8), allocatable :: diagonal(:)

        integer :: hdiag1_dims(2), hdiag2_dims(4), hdiag1_rank, hdiag2_rank
        complex(8), pointer, contiguous:: hdiag1_hp_tens(:,:), hdiag2_hp_tens(:,:,:,:)
        complex(8), pointer, contiguous:: hdiag1_h_p_tens(:),  hdiag2_h_p_tens(:,:,:)
        type(C_PTR):: hdiag1_body_p, hdiag2_body_p

        integer :: trial1_dims(2), trial2_dims(4), trial1_rank, trial2_rank
        complex(8), pointer, contiguous:: trial1_hp_tens(:,:,:), trial2_hp_tens(:,:,:,:,:)
        complex(8), pointer, contiguous:: trial1_h_p_tens(:,:),  trial2_h_p_tens(:,:,:,:)
        type(C_PTR):: trial1_body_p, trial2_body_p

        ierr = talsh_tensor_dimensions(workspace%Hdiag1,hdiag1_rank,hdiag1_dims)
        write (6,*) 'hdiag1_rank',hdiag1_rank
        write (6,*) 'hdiag1_dims',hdiag1_dims
        write (6,*) 'hdiag1_body_p',hdiag1_body_p
        ierr = talsh_tensor_get_body_access(workspace%Hdiag1,hdiag1_body_p,C8,0,DEV_HOST)
        write (6,*) 'hdiag1_body_p',hdiag1_body_p

!       ierr = talsh_tensor_dimensions(workspace%Hdiag2,hdiag2_rank,hdiag2_dims)
!       ierr = talsh_tensor_get_body_access(workspace%Hdiag2,hdiag2_body_p,C8,0,DEV_HOST)

        select case(flavor)
           case("IP")
              call c_f_pointer(hdiag1_body_p,hdiag1_h_p_tens,hdiag1_dims(1:hdiag1_rank))
              call c_f_pointer(hdiag2_body_p,hdiag2_h_p_tens,hdiag2_dims(1:hdiag2_rank))

              write (6,*) 'foo',hdiag1_dims(1)
              write (6,*) 'bar',hdiag1_h_p_tens(1:hdiag1_dims(1))

              allocate(diagonal(hdiag1_dims(1)))
            
              do i = 1, hdiag1_dims(1)
                 diagonal(i) = realpart(hdiag1_h_p_tens(i))
              end do

        end select

        write (6,*) 'debug >>> hdiag1_h_p_tens >>>'
        write (6,*) diagonal
        write (6,*) 'debug <<< hdiag1_h_p_tens <<<'

        iroot = 1
        do while (iroot.le.workspace%nroots)
!          k = minloc(diagonal, MASK=(diag_H .ge. in_window_HbarSS_threshold))
           k = minloc(diagonal)
           i = k(1)
           write(6,*) 'debug, roots',iroot,diagonal(i),i
           diagonal(i) = huge(diagonal)
           write(6,*) 'debug, making diagonal huge',diagonal(i)

           ierr = talsh_tensor_dimensions(workspace%trial1(iroot),trial1_rank,trial1_dims)
           write (6,*) 'trial1_rank',trial1_rank
           write (6,*) 'trial1_dims',trial1_dims
           write (6,*) 'trial1_body_p',trial1_body_p
           ierr = talsh_tensor_get_body_access(workspace%trial1(iroot),trial1_body_p,C8,0,DEV_HOST)
           write (6,*) 'trial1_body_p',trial1_body_p
!       ierr = talsh_tensor_dimensions(workspace%trial2,trial2_rank,trial2_dims)
!       ierr = talsh_tensor_get_body_access(workspace%trial2,trial2_body_p,C8,0,DEV_HOST)

!! we may have diagonal elements which do not correspond to an element within the windows, so we check that
!! this is the case
!          if (index_in_window(i).eqv..true.) then
           if (.true.) then
              select case(flavor)
                 case("IP")

                    call c_f_pointer(trial1_body_p,trial1_h_p_tens,trial1_dims(1:trial1_rank))
                    call c_f_pointer(trial2_body_p,trial2_h_p_tens,trial2_dims(1:trial2_rank))

                    write(6,*)'determinant in window',I,trial1_rank,trial1_dims(1:trial1_rank)
                    trial1_h_p_tens(i,1) = ONE
                    iroot = iroot + 1
              end select
           end if
        write(6,*)'   ff',trial1_h_p_tens
        end do


!       select case(flavor)
!          case("IP")
!             write (6,*) 'debug >>> trial1_h_p_tens >>>'
!             do i = 1, trial1_dims(1)
!                write (6,*) trial1_h_p_tens(i,1:trial1_dims(3))
!             end do
!             write (6,*) 'debug <<< trial1_h_p_tens <<<'
!       end select


        deallocate (diagonal)

     end subroutine

    subroutine talsh_eom_Hdiag_create(workspace, eps_occ, eps_vir, flavor)
      type(talsh_davidson_workspace), intent(inout) :: workspace
      real(8), intent(in) :: eps_occ(:), eps_vir(:)
      character(len=2), intent(in) :: flavor
      integer :: nocc, nvir, i, j, a, b, ierr 
      integer :: hdiag1_dims(2), hdiag2_dims(4), tens1_rank, tens2_rank
      real(8) :: diagonal
      complex(8), pointer, contiguous:: hdiag1_hp_tens(:,:), hdiag2_hp_tens(:,:,:,:)
      complex(8), pointer, contiguous:: hdiag1_h_p_tens(:),  hdiag2_h_p_tens(:,:,:)
      type(C_PTR):: Hdiag1_body_p, Hdiag2_body_p


      diagonal = ZERO

      nocc = size(eps_occ,1)
      nvir = size(eps_vir,1)

      ierr = talsh_tensor_dimensions(workspace%Hdiag1,tens1_rank,hdiag1_dims)
      ierr = talsh_tensor_get_body_access(workspace%Hdiag1,Hdiag1_body_p,C8,0,DEV_HOST)

      ierr = talsh_tensor_dimensions(workspace%Hdiag2,tens2_rank,hdiag2_dims)
      ierr = talsh_tensor_get_body_access(workspace%Hdiag2,Hdiag2_body_p,C8,0,DEV_HOST)

      print *,'debug, in Hdiag_create: no=',nocc,'nv=',nvir
      print *,'debug, in Hdiag_create: tens1_rank',tens1_rank,'hdiag1_dims',hdiag1_dims
      print *,'debug, in Hdiag_create: tens2_rank',tens2_rank,'hdiag2_dims',hdiag2_dims

      select case(flavor)
         case("IP")
            call c_f_pointer(Hdiag1_body_p,hdiag1_h_p_tens,hdiag1_dims(1:1))
            call c_f_pointer(Hdiag2_body_p,hdiag2_h_p_tens,hdiag2_dims(1:3))

     write (6,*) 'hdiag1_dims',hdiag1_dims 
     write (6,*) 'hdiag1_dims(1)',hdiag1_dims(1) 
     write (6,*) 'hdiag1_dims(2)',hdiag1_dims(2) 
     write (6,*) 'hdiag2_dims',hdiag2_dims 
     write (6,*) 'hdiag2_dims(1)',hdiag2_dims(1) 
     write (6,*) 'hdiag2_dims(2)',hdiag2_dims(2) 
     write (6,*) 'hdiag2_dims(3)',hdiag2_dims(3) 
     write (6,*) 'hdiag2_dims(4)',hdiag2_dims(4) 

            do i = 1, nocc 
               diagonal = -eps_occ(i)
               hdiag1_h_p_tens(i) = diagonal
            end do

     write (6,*) 'debug >>> hdiag1_h_p_tens >>>'
     write (6,*) hdiag1_h_p_tens(1:nocc)
     write (6,*) 'debug <<< hdiag1_h_p_tens <<<'

            do i = 1, nocc 
               do j = 1, nocc 
                  do a = 1, nvir 
                     diagonal = eps_vir(a) - eps_occ(j) - eps_occ(i)
                     hdiag2_h_p_tens(a,j,i) = diagonal 
                  end do
               end do
            end do

     write (6,*) 'debug >>> hdiag2_h_p_tens >>>'
     write (6,*) hdiag2_h_p_tens(1:nvir,1:nocc,1:nocc)
!    write (6,*) 'debug <<< hdiag2_h_p_tens <<<'

         case("EA")
            call c_f_pointer(Hdiag1_body_p,hdiag1_h_p_tens,hdiag1_dims(1:1))
            call c_f_pointer(Hdiag2_body_p,hdiag2_h_p_tens,hdiag2_dims(1:3))
            do a = 1, nvir 
               diagonal = eps_vir(a)
               hdiag1_h_p_tens(a) = diagonal
            end do

            do a = 1, nvir 
               do b = 1, nvir 
                  do i = 1, nocc 
                     diagonal = eps_vir(a) + eps_vir(b) - eps_occ(i)
                     hdiag2_h_p_tens(a,b,i) = diagonal
                  end do
               end do
            end do

         case("EE")
            call c_f_pointer(Hdiag1_body_p,hdiag1_hp_tens,hdiag1_dims(1:2))
            call c_f_pointer(Hdiag2_body_p,hdiag2_hp_tens,hdiag2_dims(1:4))
            do i = 1, nocc
               do a = 1, nvir
                  diagonal = eps_vir(a) - eps_occ(i)
                  hdiag1_hp_tens(a,i) = diagonal
               end do
            end do

            do i = 1, nocc
               do j = 1, nocc
                  do a = 1, nvir
                     do b = 1, nvir
                        diagonal = eps_vir(b) + eps_vir(a) - eps_occ(j) - eps_occ(i)
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
      call talsh_tensor_print_info(workspace_t%Hdiag1)
      call talsh_tensor_print_info(workspace_t%Hdiag2)
      call talsh_eom_Hdiag_create(workspace_t, eps_occ, eps_vir, exa_input%eom_flavour)
      call talsh_tensor_print_info(workspace_t%Hdiag1)
      call talsh_tensor_print_info(workspace_t%Hdiag2)
      call talsh_eom_trial_create(workspace_t, exa_input%eom_flavour)

!************************************
!     EOM calculation: IP, EA or EE *
!************************************

      if (.not.exa_input%ccd) then
        select case(exa_input%eom_flavour)
          case("IP")
            call eom_ip(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
          case("EA")
            call eom_ea(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
          case("EE")
            if (.not.exa_input%eom_properties) then
              call eom_ee(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor=t1_tensor)
            else
              call eom_ee(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor,t1_tensor)
            end if
          case default
            call quit ('Error: No valid EOM-flavour chosen')
        end select
      else
        select case(exa_input%eom_flavour)
          case("IP")
            call eom_ip(exa_input,workspace_t,projector_t,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
          case("EA")
            call eom_ea(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
          case("EE")
            if (.not.exa_input%eom_properties) then
              call eom_ee(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor)
            else
              call eom_ee(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor)
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
      if (exa_input%eom_properties) oov_dims = (/nocc,nocc,nvir/)

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
      ierr=talsh_tensor_destruct(eom_t%wbar_vovo)

!     create all coefficient and sigma tensors
      ierr=talsh_tensor_construct(eom_t%r1_coeff,C8,(/nocc/),init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r1_sigma,C8,(/nocc/),init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l1_coeff,C8,(/nocc/),init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l1_sigma,C8,(/nocc/),init_val=ZERO)
      end if
      ierr=talsh_tensor_construct(eom_t%r2_coeff,C8,voo_dims,init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r2_sigma,C8,voo_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l2_coeff,C8,oov_dims,init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l2_sigma,C8,oov_dims,init_val=ZERO)
      end if

!***********************************
!     calculate right sigma vector *
!***********************************

      call talsh_davidson_driver(exa_input,int_t,eom_t,t2_tensor, workspace_t, projector_t)
      call print_date("Calculated EOM-IP RHS")

!************************************************
!     calculate left eigenvector for properties *
!************************************************

      if (exa_input%eom_properties) then
        call sigma_left_ip(exa_input,int_t,eom_t,t2_tensor)
        call print_date("Calculated EOM-IP LHS")
      end if

!**************
!     cleanup *
!**************

      ierr=talsh_tensor_destruct(eom_t%r1_coeff)
      ierr=talsh_tensor_destruct(eom_t%r1_sigma)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(eom_t%l1_coeff)
        ierr=talsh_tensor_destruct(eom_t%l1_sigma)
      end if
      ierr=talsh_tensor_destruct(eom_t%r2_coeff)
      ierr=talsh_tensor_destruct(eom_t%r2_sigma)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(eom_t%l2_coeff)
        ierr=talsh_tensor_destruct(eom_t%l2_sigma)
      end if
      ierr=talsh_tensor_destruct(eom_t%fbar_oo)
      ierr=talsh_tensor_destruct(eom_t%fbar_ov)
      ierr=talsh_tensor_destruct(eom_t%fbar_vv)
      ierr=talsh_tensor_destruct(eom_t%w_oooo)
      ierr=talsh_tensor_destruct(eom_t%w_ooov)
      ierr=talsh_tensor_destruct(eom_t%w_ovoo)
      ierr=talsh_tensor_destruct(eom_t%w_vovo)

    end subroutine eom_ip

    subroutine eom_ea(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,t1_tensor)
       use talsh_eom_sigma, only : sigma_right_ea, sigma_left_ea

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
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

!     create all coefficient and sigma tensors
      ierr=talsh_tensor_construct(eom_t%r1_coeff,C8,(/nvir/),init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r1_sigma,C8,(/nvir/),init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l1_coeff,C8,(/nvir/),init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l1_sigma,C8,(/nvir/),init_val=ZERO)
      end if
      ierr=talsh_tensor_construct(eom_t%r2_coeff,C8,vvo_dims,init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r2_sigma,C8,vvo_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l2_coeff,C8,ovv_dims,init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l2_sigma,C8,ovv_dims,init_val=ZERO)
      end if

!***********************************
!     calculate right sigma vector *
!***********************************

      call sigma_right_ea(exa_input,int_t,eom_t,t2_tensor)
      call print_date("Calculated EOM-EA RHS")

!************************************************
!     calculate left eigenvector for properties *
!************************************************

      if (exa_input%eom_properties) then
        call sigma_left_ea(exa_input,int_t,eom_t,t2_tensor)
        call print_date("Calculated EOM-EA LHS")
      end if

!**************
!     cleanup *
!**************

      ierr=talsh_tensor_destruct(eom_t%r1_coeff)
      ierr=talsh_tensor_destruct(eom_t%r1_sigma)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(eom_t%l1_coeff)
        ierr=talsh_tensor_destruct(eom_t%l1_sigma)
      end if
      ierr=talsh_tensor_destruct(eom_t%r2_coeff)
      ierr=talsh_tensor_destruct(eom_t%r2_sigma)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(eom_t%l2_coeff)
        ierr=talsh_tensor_destruct(eom_t%l2_sigma)
      end if
      ierr=talsh_tensor_destruct(eom_t%fbar_oo)
      ierr=talsh_tensor_destruct(eom_t%fbar_ov)
      ierr=talsh_tensor_destruct(eom_t%fbar_vv)
      ierr=talsh_tensor_destruct(eom_t%w_vovo)
      ierr=talsh_tensor_destruct(eom_t%w_vovv)
      ierr=talsh_tensor_destruct(eom_t%w_vvvo)
      ierr=talsh_tensor_destruct(eom_t%w_vvvv)

    end subroutine eom_ea

    subroutine eom_ee(exa_input,eom_t,t2_tensor,int_t,foo_tensor,fov_tensor,fvv_tensor,l2_tensor,t1_tensor)
       use talsh_eom_sigma, only : sigma_right_ee, sigma_left_ee

!     input variables
      type(exacc_input),     intent(in   ) :: exa_input
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

!     create all coefficient and sigma tensors
      ierr=talsh_tensor_construct(eom_t%r1_coeff,C8,vo_dims,init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r1_sigma,C8,vo_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l1_coeff,C8,ov_dims,init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l1_sigma,C8,ov_dims,init_val=ZERO)
      end if
      ierr=talsh_tensor_construct(eom_t%r2_coeff,C8,vvoo_dims,init_val=ONE)
      ierr=talsh_tensor_construct(eom_t%r2_sigma,C8,vvoo_dims,init_val=ZERO)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_construct(eom_t%l2_coeff,C8,oovv_dims,init_val=ONE)
        ierr=talsh_tensor_construct(eom_t%l2_sigma,C8,oovv_dims,init_val=ZERO)
      end if

!***********************************
!     calculate right sigma vector *
!***********************************

      call sigma_right_ee(exa_input,int_t,eom_t,t2_tensor)
      call print_date("Calculated EOM-EE RHS")

!************************************************
!     calculate left eigenvector for properties *
!************************************************

!     create mixed type intermediates
      if (exa_input%eom_properties) then
        call get_goo_lambda(eom_t%goo,l2_tensor,t2_tensor)
        call get_gvv_lambda(eom_t%gvv,l2_tensor,t2_tensor)
      end if

      if (exa_input%eom_properties) then
        if (.not.CCD) then
          call sigma_left_ee(exa_input,int_t,eom_t,t1_tensor)
        else
          call sigma_left_ee(exa_input,int_t,eom_t)
        end if
        call print_date("Calculated EOM-EE LHS")
      end if

!**************
!     cleanup *
!**************

      ierr=talsh_tensor_destruct(eom_t%r1_coeff)
      ierr=talsh_tensor_destruct(eom_t%r1_sigma)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(eom_t%l1_coeff)
        ierr=talsh_tensor_destruct(eom_t%l1_sigma)
      end if
      ierr=talsh_tensor_destruct(eom_t%r2_coeff)
      ierr=talsh_tensor_destruct(eom_t%r2_sigma)
      if (exa_input%eom_properties) then
        ierr=talsh_tensor_destruct(eom_t%l2_coeff)
        ierr=talsh_tensor_destruct(eom_t%l2_sigma)
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

end module
