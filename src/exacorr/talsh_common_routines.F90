module talsh_common_routines

!This module contains subroutines both needed in solving t equations as well as Lambda equations.
!get_orbital_energies,diis,scale_with_denominators

    use talsh
    use tensor_algebra

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)
    private

    interface scale_with_denominators
        module procedure scale_with_denominators
    end interface scale_with_denominators

    public get_orbital_energies
    public scale_with_denominators
    public print_tensor

    contains

      subroutine scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor,t2_tensor,t3_tensor,eps_ijk)

       type(talsh_tens_t), optional, intent(inout)  :: t1_tensor, t2_tensor, t3_tensor
       real(8), intent(in)                          :: eps_occ(:),eps_vir(:)
       real(8), intent(in), optional                :: eps_ijk

       real(8)  :: denominator
       complex(8), pointer, contiguous:: t1_tens(:,:), t2_tens(:,:,:,:), t3_tens(:,:,:)
       type(C_PTR):: body_p

       integer:: t1_dims(1:2), t2_dims(1:4), t3_dims(1:3)

       integer   :: ierr, tens_rank, a, b, c, i, j

       integer, optional, intent(in) :: nocc

       if (present(t1_tensor)) then

          ierr = talsh_tensor_dimensions(t1_tensor,tens_rank,t1_dims)
          if (ierr.ne.0 .or. tens_rank.ne.2) stop 'program error: t1 tensor corrupted'
          ierr=talsh_tensor_get_body_access(t1_tensor,body_p,C8,0,DEV_HOST)
          call c_f_pointer(body_p,t1_tens,t1_dims(1:2)) ! to use <t1_tens> as a regular Fortran 2d array

          !if (t1_dims(2) .eq. nocc*2) then
          if (t1_dims(2) .eq. nocc) then
           !t amplitudes
           do i = 1, t1_dims(2)
             do a = 1, t1_dims(1)
                denominator = eps_occ(i) - eps_vir(a)
                t1_tens(a,i) = t1_tens(a,i) / denominator
             end do
           end do
          else
           !Lambda
           do i = 1, t1_dims(1)
             do a = 1, t1_dims(2)
                denominator = eps_occ(i) - eps_vir(a)
                t1_tens(i,a) = t1_tens(i,a) / denominator
             end do
           end do
          end if

       end if

       if (present(t2_tensor)) then

          ierr = talsh_tensor_dimensions(t2_tensor,tens_rank,t2_dims)
          if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: t2 tensor corrupted'
          ierr=talsh_tensor_get_body_access(t2_tensor,body_p,C8,0,DEV_HOST)
          call c_f_pointer(body_p,t2_tens,t2_dims(1:4)) ! to use <t2_tens> as a regular Fortran 4d array

          !if (t2_dims(3) .eq. nocc*2) then
          if (t2_dims(3) .eq. nocc) then
           !t amplitudes
           do j = 1, t2_dims(4)
             do i = 1, t2_dims(3)
               do b = 1, t2_dims(2)
                 do a = 1, t2_dims(1)
                    denominator = eps_occ(i) + eps_occ(j) - eps_vir(a) - eps_vir(b)
                    t2_tens(a,b,i,j) = t2_tens(a,b,i,j) / denominator
                 end do
               end do
             end do
           end do
          else
           !Lambda
           do i = 1, t2_dims(1)
             do j = 1, t2_dims(2)
               do a = 1, t2_dims(3)
                 do b = 1, t2_dims(4)
                    denominator = eps_occ(i) + eps_occ(j) - eps_vir(a) - eps_vir(b)
                    t2_tens(i,j,a,b) = t2_tens(i,j,a,b) / denominator
                 end do
               end do
             end do
           end do
          end if

       end if

       if (present(t3_tensor)) then

          ierr = talsh_tensor_dimensions(t3_tensor,tens_rank,t3_dims)
          if (ierr.ne.0 .or. tens_rank.ne.3) stop 'program error: t3 tensor corrupted'
          ierr=talsh_tensor_get_body_access(t3_tensor,body_p,C8,0,DEV_HOST)
          call c_f_pointer(body_p,t3_tens,t3_dims(1:3)) ! to use <t3_tens> as a regular Fortran 3d array

          do c = 1, t3_dims(3)
            do b = 1, t3_dims(2)
              do a = 1, t3_dims(1)
                 denominator = eps_ijk - eps_vir(a) - eps_vir(b) - eps_vir(c)
                 t3_tens(a,b,c) = t3_tens(a,b,c) / denominator
              end do
            end do
          end do

       end if

      end subroutine scale_with_denominators

      subroutine get_orbital_energies (nmo, mo_list, eps)

!This subroutine reads the desired subset of mo energies from DIRAC

!      Written by Lucas Visscher, winter 2016/2017 (but in Goa, India, temperature about 25 C)

       use exacorr_mo
       use exacorr_global

       integer, intent(in ) :: nmo          ! the length of the mo basis
       integer, intent(in ) :: mo_list(:)   ! and their indices
       real(8), intent(out) :: eps(:)
       type(cmo)            :: cspinor

       call get_mo_coefficients (cspinor,mo_list,nmo)
       eps = cspinor%energy
       call dealloc_mo(cspinor)

      end subroutine get_orbital_energies

      subroutine print_tensor(h_tensor, rank, acc, t_name)

        integer, intent(in)                  :: rank
        type(talsh_tens_t), intent(inout)    :: h_tensor
        real(8), intent(in)                  :: acc
        character(LEN=*), intent(in)         :: t_name

        integer  ::  dims2(1:2)
        integer  ::  dims4(1:4)
        integer  ::  ierr
        integer  ::  rrank
        integer  ::  i, j, k, l
        type(C_PTR)         :: body_p
        complex(8), pointer :: tens2(:, :)
        complex(8), pointer :: tens4(:, :, :, :)
        
        if (rank.eq.2) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'

          ierr=talsh_tensor_get_body_access(h_tensor,body_p,C8,0,DEV_HOST)
          call c_f_pointer(body_p,tens2,dims2)

          print *, "print 2D tensor: ", t_name
          do i = 1, dims2(1)
            do j = 1, dims2(2)
              if(abs(real(tens2(i,j),8)).gt.acc .or. abs(dimag(tens2(i,j))).gt.acc) then
                print *, i, j, tens2(i,j)
              end if
            end do
          end do
          print *, " end print 2D tensor ", t_name

        else if (rank.eq.4) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'

          ierr=talsh_tensor_get_body_access(h_tensor,body_p,C8,0,DEV_HOST)
          call c_f_pointer(body_p,tens4,dims4)

          print *, "print 4D tensor: ", t_name
          do i = 1, dims4(1)
            do j = 1, dims4(2)
              do k = 1, dims4(3)
                do l = 1, dims4(4)
                  if(abs(real(tens4(i,j,k,l),8)).gt.acc .or. abs(dimag(tens4(i,j,k,l))).gt.acc) then
                    print *, i, j, k, l, tens4(i, j, k, l)
                  end if
                end do
              end do
            end do
          end do
          print *, " end print 4D tensor: ", t_name

        else
          stop 'error in print_tensor: only rank 2 and 4 implemented'
        end if

      end subroutine print_tensor

end module talsh_common_routines
