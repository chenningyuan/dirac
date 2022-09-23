module talsh_trial_eomccsd

  use talsh
  use tensor_algebra
  use exacorr_datatypes
  use intermediates
  use talsh_common_routines, only : tensor_norm2, print_tensor
  use exacorr_utils,  only : print_date
  use talsh_davidson

  implicit none

  private

  complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE=(-1.D0,0.D0),ONE=(1.D0,0.D0),ONE_HALF=(0.5D0,0.D0)
  real(8),    parameter :: ONE_QUARTER=0.25D0

  public talsh_eomccsd_ip_trial_create 
  public talsh_eomccsd_ee_trial_create 

  public talsh_eomccsd_ip_Hdiag_create
  public talsh_eomccsd_ea_Hdiag_create
  public talsh_eomccsd_ee_Hdiag_create

  contains

    subroutine from_linear_to_indexes_pqr(b_desired, p, dim_p, q, dim_q, r, dim_r)
       integer, intent(in)  :: b_desired, dim_p, dim_q, dim_r
       integer, intent(out) :: p, q, r
       integer :: b

       b = 0
       do r = 1, dim_r
          do q = 1, dim_q
             do p = 1, dim_p
                b = b + 1
!               print *, 'i:',i,'j:',j,'a:',a,'   b:',b
                if  (b .eq. b_desired) return
             end do
          end do
       end do
    end subroutine


     subroutine talsh_eomccsd_ip_trial_create(workspace)
        type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace

        integer :: nocc, nvir, i, j, a, b, X, ierr, iroot, k(1), N, N1, N2 
        real(8), allocatable :: diagonal(:)
        integer :: i1, j1, a1, b1

        complex(8), pointer, contiguous:: hdiag1_h_p_tens(:,:),  hdiag2_h_p_tens(:,:,:,:)
        complex(8), pointer, contiguous:: trial1_h_p_tens(:,:),  trial2_h_p_tens(:,:,:,:)
        type(C_PTR):: hdiag1_body_p, hdiag2_body_p
        type(C_PTR):: trial1_body_p, trial2_body_p

        logical :: enable_double_manifold = .true.
        logical :: debug_trial = .false.

        ierr = talsh_tensor_get_body_access(workspace%Hdiag1,hdiag1_body_p,C8,0,DEV_HOST)
        ierr = talsh_tensor_get_body_access(workspace%Hdiag2,hdiag2_body_p,C8,0,DEV_HOST)

! first we copy the diagonal to a local working variable ("diagonal"), in order to sort it later with minloc() 
        call c_f_pointer(hdiag1_body_p,hdiag1_h_p_tens,workspace%tens1_dims(1:workspace%tens1_rank))
        call c_f_pointer(hdiag2_body_p,hdiag2_h_p_tens,workspace%tens2_dims(1:workspace%tens2_rank))

        if (workspace%nroots.le.workspace%N(1)) enable_double_manifold = .false.

        if (enable_double_manifold) then
           allocate(diagonal(workspace%N(1) + workspace%N(2)))
        else
           allocate(diagonal(workspace%N(1)))
        end if
        diagonal = 0.0d0

        b = 0
        do i = 1, workspace%tens1_dims(2) 
           do X = 1, workspace%tens1_dims(1) 
              b = b + 1
              diagonal(b) = real(hdiag1_h_p_tens(X,i))
           end do
        end do

! aspg: have to still figure out how to handle the single and double manifolds handled separately
! and find which are the workspace%nroots lowest roots more generally (i.e. when some of the trials
! belong to single and some to double
!
        if (enable_double_manifold) then
            b = workspace%N(1)
            do j = 1, workspace%tens2_dims(4)
               do i = 1, workspace%tens2_dims(3)
                  do X = 1, workspace%tens2_dims(2)
                     do a = 1, workspace%tens2_dims(1)
                        b = b + 1
                        if (i.lt.j) then
                           diagonal(b) = real(hdiag2_h_p_tens(a,X,j,i))
                        else 
                           diagonal(b) = 0.1*huge(diagonal)
                        end if
                     end do
                  end do
               end do
            end do
         end if

        if (debug_trial) then
           do i = 1, size(diagonal,1)
               write (6,*) 'H_II(',i,'):',diagonal(i)
           end do

        end if

! now set up the initial trial vectors, as unit matrices (or the equivalent for higher-dimensional tensors)

        if (debug_trial) then 
           do iroot = 1, workspace%nroots
   
              ierr = talsh_tensor_get_body_access(workspace%trial1(iroot),trial1_body_p,C8,0,DEV_HOST)
              ierr = talsh_tensor_get_body_access(workspace%trial2(iroot),trial2_body_p,C8,0,DEV_HOST)
              call c_f_pointer(trial1_body_p,trial1_h_p_tens,workspace%tens1_dims(1:workspace%tens1_rank))
              call c_f_pointer(trial2_body_p,trial2_h_p_tens,workspace%tens2_dims(1:workspace%tens2_rank))

              write (6,*) 'debug >>> t1_h_p_tens, root',iroot,' >>>'
!             call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
              do i = 1, workspace%tens1_dims(2) 
                 do X = 1, workspace%tens1_dims(1) 
                    write (6,*) X, i, trial1_h_p_tens(X,i)
                 end do
              end do
              write (6,*) 'debug <<< t1_h_p_tens <<<'

              write (6,*) 'debug >>> t2_h_p_tens, root',iroot,' >>>'
              call talsh_tensor_print_info(workspace%trial2(iroot))

!             call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
              do j = 1, workspace%tens2_dims(4)
                 do i = 1, workspace%tens2_dims(3)
                    do X = 1, workspace%tens2_dims(2)
                       do a = 1, workspace%tens2_dims(1)
                          b = b + 1
                          write (6,*) b, a, X, i, j, trial2_h_p_tens(a,X,i,j)
                       end do
                    end do
                 end do
              end do
              write (6,*) 'debug <<< t2_h_p_tens <<<'
           end do

           write (6,*) 'debug >>> real part of hdiag1_h_p_tens >>>'
           write (6,*) diagonal
           write (6,*) 'debug <<< real part of hdiag1_h_p_tens <<<'
 
        end if

! forcing full excitation manifold 
      !   if (enable_double_manifold) then
      !         workspace%nroots = workspace%N(1) + workspace%N(2)
      !          !  workspace%nroots = workspace%N(1) + (workspace%N(2)-workspace%tens2_dims(1)*workspace%tens2_dims(3))/2
      !   else
! forcing only the single part now
      !      if (workspace%nroots.gt.workspace%N(1)) then
      !         workspace%nroots = workspace%N(1)
      !         write(6,*) 'WARNING!! resetting number of roots to :',workspace%N(1)
      !      end if 
      !   end if

        write(6,*) 'roots',workspace%nroots
        j = 0
        X = 1
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

                    ierr = talsh_tensor_get_body_access(workspace%trial1(iroot),trial1_body_p,C8,0,DEV_HOST)
                    ierr = talsh_tensor_get_body_access(workspace%trial2(iroot),trial2_body_p,C8,0,DEV_HOST)
                    call c_f_pointer(trial1_body_p,trial1_h_p_tens,workspace%tens1_dims(1:workspace%tens1_rank))
                    call c_f_pointer(trial2_body_p,trial2_h_p_tens,workspace%tens2_dims(1:workspace%tens2_rank))

!                   print *,'j:',j,' i:',i,i1,j1,a1
! aspg, come back and see how to handle this
                    if (enable_double_manifold) then ! we have to go into the 2h1p section to select trial vectors
!                      print *,'a'
                       if (i .le. workspace%N(1)) then 
                          trial1_h_p_tens(X,i) = ONE
                       else
!                      print *,'b'
! translate j into the spinor indexes
                          call from_linear_to_indexes_pqrs(i-workspace%N(1), &
                                                           a1,workspace%tens2_dims(1), &
                                                           b1,workspace%tens2_dims(2),  &
                                                           i1,workspace%tens2_dims(3), &
                                                           j1,workspace%tens2_dims(4))
                          if (j1.gt.i1) then
                             trial2_h_p_tens(a1,b1,j1,i1) = ONE*(2.0**(-0.5))
                           !   trial2_h_p_tens(a1,b1,j1,i1) = ONE
                        !   else if (j1.lt.i1) then
                             trial2_h_p_tens(a1,b1,i1,j1) = MINUS_ONE*(2.0**(-0.5)) !XY
                           !   trial2_h_p_tens(a1,b1,i1,j1) = MINUS_ONE
                        !   else
                        !    trial2_h_p_tens(a1,b1,j1,i1) =ZERO
                          end if
                       end if
                    else
!                      print *,'c'
! aspg, for now, we only set the tens1 part of the trial vectors
                       trial1_h_p_tens(X,i) = ONE
                    end if
           end if
        end do

        deallocate (diagonal)

! debug information
        if (debug_trial) then
                 do iroot = 1, workspace%nroots
   
                    ierr = talsh_tensor_get_body_access(workspace%trial1(iroot),trial1_body_p,C8,0,DEV_HOST)
                    ierr = talsh_tensor_get_body_access(workspace%trial2(iroot),trial2_body_p,C8,0,DEV_HOST)
                    call c_f_pointer(trial1_body_p,trial1_h_p_tens,workspace%tens1_dims(1:workspace%tens1_rank))
                    call c_f_pointer(trial2_body_p,trial2_h_p_tens,workspace%tens2_dims(1:workspace%tens2_rank))
   
                    b = 0
                    write (6,*) 'debug >>> trial1_h_p_tens, root',iroot,' >>>'
!                   call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
                    do i = 1, workspace%tens1_dims(2)
                       do X = 1, workspace%tens1_dims(1)
                          b = b + 1
                          write (6,*) b, X, i, trial1_h_p_tens(X,i)
                       end do
                    end do
                    write (6,*) 'debug <<< trial1_h_p_tens <<<'
   
                    write (6,*) 'debug >>> trial2_h_p_tens, root',iroot,' >>>'
!                   call davidson_print_vectors_complex(trial1_h_p_tens,N1,1) !write (6,*) trial1_h_p_tens(1:trial1_dims(1))
                    do j = 1, workspace%tens2_dims(4)
                       do i = 1, workspace%tens2_dims(3)
                          do X = 1, workspace%tens2_dims(2)
                             do a = 1, workspace%tens2_dims(1)
                                b = b + 1
                                write (6,*) b, a, X, i, j, trial2_h_p_tens(a,X,i,j)
                             end do
                          end do
                       end do
                    end do
                 end do
        end if

     end subroutine

    subroutine talsh_eomccsd_ip_Hdiag_create(workspace, Fmi_tensor, Fae_tensor)
      type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace
      type(talsh_tens_t), intent(inout) :: Fmi_tensor, Fae_tensor
      integer :: nocc, nvir, ncont, i, j, a, X, ierr 
      integer :: Fmi_dims(2), Fae_dims(2), Fmi_rank, Fae_rank
      real(8) :: diagonal
      complex(8), pointer, contiguous:: hdiag1_tens(:,:),  hdiag2_tens(:,:,:,:)
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

      call c_f_pointer(hdiag1_body_p,hdiag1_tens,workspace%tens1_dims(1:workspace%tens1_rank))
      call c_f_pointer(hdiag2_body_p,hdiag2_tens,workspace%tens2_dims(1:workspace%tens2_rank))

      nvir  = workspace%tens2_dims(1) 
      ncont = workspace%tens2_dims(2)
      nocc  = workspace%tens2_dims(4)
      X = 1 ! continuum
 
      if (debug_hdiag_create) print *,'nocc',nocc,'nvir',nvir,'ncont',ncont

      do i = 1,  nocc
         diagonal = -fmi_p_tens(i,i)
         if (debug_hdiag_create) write (6,*) 'f_ii:',fmi_p_tens(i,i),'diag:',diagonal
         hdiag1_tens(X,i) = diagonal
      end do
      if (debug_hdiag_create) then
         write (6,*) 'debug >>> hdiag1_h_p_tens >>>'
!        write (6,*) hdiag1_tens(1:ncont,1:nocc)
         call print_tensor(workspace%Hdiag1, 1.0d-18, "workspace%Hdiag1")
         write (6,*) 'debug <<< hdiag1_h_p_tens <<<'
      end if

      do j = 1, nocc
         do i = 1, nocc
            do a = 1, nvir
               diagonal = fae_p_tens(a,a) - fmi_p_tens(j,j) - fmi_p_tens(i,i)
               if (debug_hdiag_create) &
               write (6,*) 'f_ii:',fmi_p_tens(i,i),'f_jj:',fmi_p_tens(j,j),'f_aa',fae_p_tens(a,a),'diag:',diagonal
               hdiag2_tens(a,X,i,j) = diagonal 
            end do
         end do
      end do
      if (debug_hdiag_create) then
         write (6,*) 'debug >>> hdiag2_h_p_tens >>>'
         call print_tensor(workspace%Hdiag2, 1.0d-18, "workspace%Hdiag2")
!        write (6,*) hdiag2_tens(1:nvir,1:ncont,1:nocc,1:nocc)
         write (6,*) 'debug <<< hdiag2_h_p_tens <<<'
      end if

    end subroutine


    subroutine talsh_eomccsd_ea_Hdiag_create(workspace, Fmi_tensor, Fae_tensor)
      type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace
      type(talsh_tens_t), intent(inout) :: Fmi_tensor, Fae_tensor
      integer :: nocc, nvir, i, a, b, ierr 
      integer :: Fmi_dims(2), Fae_dims(2), Fmi_rank, Fae_rank
      real(8) :: diagonal
      complex(8), pointer, contiguous:: hdiag1_h_p_tens(:),  hdiag2_h_p_tens(:,:,:)
      complex(8), pointer, contiguous:: Fmi_p_tens(:,:), Fae_p_tens(:,:) 
      type(C_PTR):: hdiag1_body_p, hdiag2_body_p
      type(C_PTR):: Fmi_body_p, Fae_body_p

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

    end subroutine

    subroutine talsh_eomccsd_ee_Hdiag_create(workspace, Fmi_tensor, Fae_tensor)
      type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace
      type(talsh_tens_t), intent(inout) :: Fmi_tensor, Fae_tensor
      integer :: nocc, nvir, i, j, a, b, ierr 
      integer :: Fmi_dims(2), Fae_dims(2), Fmi_rank, Fae_rank
      real(8) :: diagonal
      complex(8), pointer, contiguous:: hdiag1_hp_tens(:,:), hdiag2_hp_tens(:,:,:,:)
      complex(8), pointer, contiguous:: hdiag1_h_p_tens(:),  hdiag2_h_p_tens(:,:,:)
      complex(8), pointer, contiguous:: Fmi_p_tens(:,:), Fae_p_tens(:,:) 
      type(C_PTR):: hdiag1_body_p, hdiag2_body_p
      type(C_PTR):: Fmi_body_p, Fae_body_p

      diagonal = ZERO

      nvir  = workspace%tens1_dims(1)
      nocc  = workspace%tens1_dims(2)

      ierr = talsh_tensor_dimensions(Fmi_tensor,Fmi_rank,Fmi_dims)
      ierr = talsh_tensor_get_body_access(Fmi_tensor,Fmi_body_p,C8,0,DEV_HOST)
      call c_f_pointer(Fmi_body_p,fmi_p_tens,Fmi_dims(1:2))

      ierr = talsh_tensor_dimensions(Fae_tensor,Fae_rank,Fae_dims)
      ierr = talsh_tensor_get_body_access(Fae_tensor,Fae_body_p,C8,0,DEV_HOST)
      call c_f_pointer(Fae_body_p,fae_p_tens,Fae_dims(1:2))

      ierr = talsh_tensor_get_body_access(workspace%Hdiag1,hdiag1_body_p,C8,0,DEV_HOST)
      ierr = talsh_tensor_get_body_access(workspace%Hdiag2,hdiag2_body_p,C8,0,DEV_HOST)

      call c_f_pointer(hdiag1_body_p,hdiag1_hp_tens,workspace%tens1_dims(1:workspace%tens1_rank))
      call c_f_pointer(hdiag2_body_p,hdiag2_hp_tens,workspace%tens2_dims(1:workspace%tens2_rank))

      do i = 1, nocc
         do a = 1, nvir
            diagonal = fae_p_tens(a,a) - fmi_p_tens(i,i)
            hdiag1_hp_tens(a,i) = diagonal
!           print *,'--> a:',(a+nocc),' i:',i,' hdiag1_hp_tens(a,i):',hdiag1_hp_tens(a,i)
!           print *,'    faa:',fae_p_tens(a,a)
!           print *,'    fii:',fmi_p_tens(i,i)
         end do
      end do

      do j = 1, nocc
         do i = 1, nocc
            do b = 1, nvir
               do a = 1, nvir
                  diagonal = fae_p_tens(a,a) + fae_p_tens(b,b) - fmi_p_tens(j,j) - fmi_p_tens(i,i)
                if ((a.eq.b).or.(i.eq.j)) then 
                  ! if (i.lt.j.and.a.lt.b) then
! matching indexes of virtuals and/or occupied correspond to unphysical excited determinants and thus spurious matrix elements
! so we just project those out by making them huge
                     hdiag2_hp_tens(a,b,i,j) = huge(diagonal)
                  else
                     hdiag2_hp_tens(a,b,i,j) = diagonal
                  end if
!                 print *,'--> a:',(a+nocc),' b:',(b+nocc),&
!    &                    ' i:',i,' j:',j,'  hdiag2_hp_tens(a,b,i,j)',diagonal
!                 print *,'    faa:',fae_p_tens(a,a)
!                 print *,'    fbb:',fae_p_tens(b,b)
!                 print *,'    fii:',fmi_p_tens(i,i)
!                 print *,'    fjj:',fmi_p_tens(j,j)
               end do
            end do
         end do
      end do

      call print_tensor(workspace%Hdiag1, 1.0d-18, "Hdiag1")
      call print_tensor(workspace%Hdiag2, 1.0d-18, "Hdiag2")

    end subroutine


    subroutine from_linear_to_indexes_pq(k_desired, p, dim_p, q, dim_q)
       integer, intent(in)  :: k_desired, dim_p, dim_q
       integer, intent(out) :: p, q
       integer :: k

       k = 0
       do q = 1, dim_q
          do p = 1, dim_p
             k = k + 1
             if  (k .eq. k_desired) return
          end do
       end do
    end subroutine

    subroutine from_linear_to_indexes_pqrs(k_desired, p, dim_p, q, dim_q, r, dim_r, s, dim_s)
       integer, intent(in)  :: k_desired, dim_p, dim_q, dim_r, dim_s
       integer, intent(out) :: p, q, r, s
       integer :: k

       k = 0
       do s = 1, dim_s
          do r = 1, dim_r
             do q = 1, dim_q
                do p = 1, dim_p
                   k = k + 1
                   if  (k .eq. k_desired) return
                end do 
             end do 
          end do
       end do
    end subroutine 

     subroutine talsh_eomccsd_ee_trial_create(workspace)
        type(talsh_davidson_workspace_tens_t), intent(inout) :: workspace

        integer :: nocc, nvir, i, j, a, b, ierr, iroot, k(1), N, N1, N2 
        real(8), allocatable :: diagonal(:)
        integer :: i1, j1, a1, b1
        integer :: w

        complex(8), pointer, contiguous:: hdiag1_hp_tens(:,:),  hdiag2_hp_tens(:,:,:,:)
        complex(8), pointer, contiguous:: trial1_hp_tens(:,:),  trial2_hp_tens(:,:,:,:)
        type(C_PTR):: hdiag1_body_p, hdiag2_body_p
        type(C_PTR):: trial1_body_p, trial2_body_p

        logical :: enable_double_manifold = .false.
        logical :: debug_trial = .false.

        ierr = talsh_tensor_get_body_access(workspace%Hdiag1,hdiag1_body_p,C8,0,DEV_HOST)
        ierr = talsh_tensor_get_body_access(workspace%Hdiag2,hdiag2_body_p,C8,0,DEV_HOST)

! first we copy the diagonal to a local working variable ("diagonal"), in order to sort it later with minloc() 
        call c_f_pointer(hdiag1_body_p,hdiag1_hp_tens,workspace%tens1_dims(1:workspace%tens1_rank))
        call c_f_pointer(hdiag2_body_p,hdiag2_hp_tens,workspace%tens2_dims(1:workspace%tens2_rank))

        if (workspace%nroots.gt.workspace%N(1)) then
           enable_double_manifold = .true.
        end if

        if (enable_double_manifold) then
           allocate(diagonal(workspace%N(1) + workspace%N(2)))
        else
           allocate(diagonal(workspace%N(1)))
        end if
        diagonal = 0.0d0

        w = 0
        do i = 1, workspace%tens1_dims(2) 
           do a = 1, workspace%tens1_dims(1) 
              w = w + 1
              diagonal(w) = real(hdiag1_hp_tens(a,i))
           end do
        end do

! aspg: have to still figure out how to handle the single and double manifolds handled separately
! and find which are the workspace%nroots lowest roots more generally (i.e. when some of the trials
! belong to single and some to double
!
        if (enable_double_manifold) then
           do j = 1, workspace%tens2_dims(4)
              do i = 1, workspace%tens2_dims(3)
                 do b = 1, workspace%tens2_dims(2)
                    do a = 1, workspace%tens2_dims(1)
                       w = w + 1
! aspg, once we rethink antisymmetrization we should be able to just pick out the elements of Hdiag instead of redoing a construction
! of the diagonal here.
!                      diagonal(w) = real(hdiag2_hp_tens(a,b,i,j))
! the code below works nicely for when we want to only operate on the unique (triangular) index combinations
! and in the trial vector generation swap indices around while at the same time ensuring antisymetrization of elements of the tensor 
                       if (i.lt.j.and.a.lt.b) then
                         diagonal(w) = real(hdiag2_hp_tens(a,b,i,j))
                       else 
                         diagonal(w) = huge(diagonal)
                       end if
                    end do
                 end do
              end do
           end do
        end if

!       if (debug_trial) then
        if (.true.) then
           do i = 1, size(diagonal,1)
               write (6,*) 'H_II(',i,'):',diagonal(i)
           end do
        end if

! now set up the initial trial vectors, as unit matrices (or the equivalent for higher-dimensional tensors)

        if (debug_trial) then 
           do iroot = 1, workspace%nroots
              write (6,*) 'debug >>> t1_h_p_tens, root',iroot,' >>>'
              call print_tensor(workspace%trial1(iroot), 0.0d-18, "trial1")
              write (6,*) 'debug <<< t1_hp_tens <<<'

              write (6,*) 'debug >>> t2_hp_tens, root',iroot,' >>>'
              call print_tensor(workspace%trial2(iroot), 0.0d-18, "trial2")
              write (6,*) 'debug <<< t2_hp_tens <<<'
           end do

        end if

        write(6,*) 'debug, number of roots',workspace%nroots
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

                    ierr = talsh_tensor_get_body_access(workspace%trial1(iroot),trial1_body_p,C8,0,DEV_HOST)
                    ierr = talsh_tensor_get_body_access(workspace%trial2(iroot),trial2_body_p,C8,0,DEV_HOST)
                    call c_f_pointer(trial1_body_p,trial1_hp_tens,workspace%tens1_dims(1:workspace%tens1_rank))
                    call c_f_pointer(trial2_body_p,trial2_hp_tens,workspace%tens2_dims(1:workspace%tens2_rank))

                    call from_linear_to_indexes_pq(i,a1,workspace%tens1_dims(1),i1,workspace%tens1_dims(2))
! aspg, come back and see how to handle this
                    if (enable_double_manifold) then ! we have to go into the 2h1p section to select trial vectors
                       if (i.le.workspace%N(1)) then 
                        ! print *,'w:',a,' a:',i,a1,i1
                          trial1_hp_tens(a1,i1) = ONE
                       else
! translate j into the spinor indexes
                          call from_linear_to_indexes_pqrs(i-workspace%N(1),&
                                                           a1,workspace%tens2_dims(1), &
                                                           b1,workspace%tens2_dims(2), &
                                                           i1,workspace%tens2_dims(3), &
                                                           j1,workspace%tens2_dims(4))
!                         print *,'w:',i,' i:',i,a1,b1,i1,j1
                        !   trial2_hp_tens(a1,b1,i1,j1) = ONE

! below we mean to pick out a unique element of |abij>, where a>b and i>j and generate the
! permutations and antisymmetrize them
! see comments above concerning the diagonal for the trial vector generation.
                          if (j1.gt.i1.and.b1.gt.a1) then
                              ! trial2_hp_tens(b1,a1,j1,i1) = ONE*(2.0**(-1.0))
                              ! trial2_hp_tens(b1,a1,i1,j1) = MINUS_ONE*(2.0**(-1.0)) !XY
                              ! trial2_hp_tens(a1,b1,j1,i1) = MINUS_ONE*(2.0**(-1.0))
                              ! trial2_hp_tens(a1,b1,i1,j1) = ONE*(2.0**(-1.0)) !XY
                              trial2_hp_tens(b1,a1,j1,i1) = ONE
                              trial2_hp_tens(b1,a1,i1,j1) = MINUS_ONE
                              trial2_hp_tens(a1,b1,j1,i1) = MINUS_ONE
                              trial2_hp_tens(a1,b1,i1,j1) = ONE
                           !   trial2_hp_tens(a1,b1,i1,j1) = ONE
                          end if


                       end if
                    else
!                      print *,'c'
! aspg, for now, we only set the tens1 part of the trial vectors
                     print *, "a:", a1, "i1:", i1
                       trial1_hp_tens(a1,i1) = ONE
                    end if
           end if
        end do

        deallocate (diagonal)

! debug information
        if (debug_trial) then
           do iroot = 1, workspace%nroots
              write (6,*) 'debug >>> t1_h_p_tens, root',iroot,' >>>'
              call print_tensor(workspace%trial1(iroot), 0.0d-18, "trial1")
              write (6,*) 'debug <<< t1_hp_tens <<<'

              write (6,*) 'debug >>> t2_hp_tens, root',iroot,' >>>'
              call print_tensor(workspace%trial2(iroot), 0.0d-18, "trial2")
              write (6,*) 'debug <<< t2_hp_tens <<<'
           end do
        end if

     end subroutine


end module
