module talsh_common_routines

!This module contains subroutines both needed in solving t equations as well as Lambda equations.
!get_orbital_energies,diis,scale_with_denominators, laplace multiplications

    use talsh
    use tensor_algebra
    use, intrinsic:: ISO_C_BINDING

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)
    private

    interface scale_with_denominators
        module procedure scale_with_denominators
    end interface scale_with_denominators

    public get_orbital_energies
    public scale_with_denominators
    public print_tensor
    public get_element
    public Pair_correlation_energies
    public talsh_init_R8
    public talsh_init_C8
    public tensor_norm2
    public laplace_mult
    public laplace_testing

    contains

      subroutine scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor,t2_tensor,t3_tensor,eps_ijk)

       type(talsh_tens_t), optional, intent(inout)  :: t1_tensor, t2_tensor, t3_tensor
       real(8), intent(in)                          :: eps_occ(:),eps_vir(:)
       real(8), intent(in), optional                :: eps_ijk

       real(8)  :: denominator
       complex(8), pointer, contiguous:: t1_tens(:,:), t2_tens(:,:,:,:), t3_tens(:,:,:)
       type(C_PTR):: body_p

       integer(INTD) :: t1_dims(1:2), t2_dims(1:4), t3_dims(1:3)

       integer(INTD) :: ierr, tens_rank
       integer       :: a, b, c, i, j

       integer, optional, intent(in) :: nocc

       if (present(t1_tensor)) then

          ierr = talsh_tensor_dimensions(t1_tensor,tens_rank,t1_dims)
          if (ierr.ne.0 .or. tens_rank.ne.2) stop 'program error: t1 tensor corrupted'
          ierr=talsh_tensor_get_body_access(t1_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
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
          ierr=talsh_tensor_get_body_access(t2_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
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
          ierr=talsh_tensor_get_body_access(t3_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
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

      subroutine print_tensor(h_tensor, acc, t_name)

        !routine for printing talsh tensors

        implicit none

        type(talsh_tens_t), intent(inout)    :: h_tensor
        real(8), intent(in)                  :: acc
        character(LEN=*), intent(in)         :: t_name

        integer(INTD)  ::  dims1(1)
        integer(INTD)  ::  dims2(1:2)
        integer(INTD)  ::  dims3(1:3)
        integer(INTD)  ::  dims4(1:4)
        integer(C_INT) ::  ierr
        integer(INTD)  ::  rank, rrank
        integer        ::  i, j, k, l
        type(C_PTR)         :: body_p
        complex(8), pointer :: tens1(:)
        complex(8), pointer :: tens2(:, :)
        complex(8), pointer :: tens3(:, :, :)
        complex(8), pointer :: tens4(:, :, :, :)
        real(8), pointer    :: tensR1(:)
        real(8), pointer    :: tensR2(:, :)
        real(8), pointer    :: tensR3(:, :, :)
        real(8), pointer    :: tensR4(:, :, :, :)
        integer             :: DataKind(1),nData

        ierr=talsh_tensor_data_kind(h_tensor,nData,DataKind)
        if (ierr.ne.0) stop 'error in getting DataKind'

        rank=talsh_tensor_rank(h_tensor)
        if (rank.eq.0) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'

          print *, "print 0D tensor: ", t_name, " (shape:", dims1,")"
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          select case (DataKind(1))
          case (C8)
            call c_f_pointer(body_p,tens1,dims1)
          case (R8)
             call c_f_pointer(body_p,tensR1,dims1)
           case default
            stop "wrong datakind in print_tensor"
          end select

          do i = 1, dims1(1)
            if (DataKind(1)==C8) then
              if(abs(tens1(i)).gt.acc) then
                print *, i, tens1(i)
              end if
            else
              if(abs(tensR1(i)).gt.acc) then
                print *, i, tensR1(i)
              end if
            end if
          end do
          print *, " end print 0D tensor ", t_name
        else if (rank.eq.1) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'

          print *, "print 1D tensor: ", t_name, " (shape:", dims1,")"
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          select case (DataKind(1))
          case (C8)
            call c_f_pointer(body_p,tens1,dims1)
          case (R8)
             call c_f_pointer(body_p,tensR1,dims1)
           case default
            stop "wrong datakind in print_tensor"
          end select
          
          do i = 1, dims1(1)
            if (DataKind(1)==C8) then
              if(abs(tens1(i)).gt.acc) then
                print *, i, tens1(i)
              end if
            else
              if(abs(tensR1(i)).gt.acc) then
                print *, i, tensR1(i)
              end if
            end if
          end do
          print *, " end print 1D tensor ", t_name
        else if (rank.eq.2) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'

          print *, "print 2D tensor: ", t_name, " (shape:", dims2,")"
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          select case (DataKind(1))
          case (C8)
            call c_f_pointer(body_p,tens2,dims2)
          case (R8)
             call c_f_pointer(body_p,tensR2,dims2)
           case default
            stop "wrong datakind in print_tensor"
          end select

          do i = 1, dims2(1)
            do j = 1, dims2(2)
              if (DataKind(1)==C8) then
                if(abs(tens2(i,j)).gt.acc) then
                  print *, i, j, tens2(i,j)
                end if
              else
                if(abs(tensR2(i,j)).gt.acc) then
                  print *, i, j, tensR2(i,j)
                end if
              end if
            end do
          end do
          print *, " end print 2D tensor ", t_name

        else if (rank.eq.3) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'

          print *, "print 3D tensor: ", t_name, " (shape:", dims3,")"
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          select case (DataKind(1))
          case (C8)
            call c_f_pointer(body_p,tens3,dims3)
          case (R8)
             call c_f_pointer(body_p,tensR3,dims3)
           case default
            stop "wrong datakind in print_tensor"
          end select

          do i = 1, dims3(1)
            do j = 1, dims3(2)
              do k = 1, dims3(3)
                if (DataKind(1)==C8) then
                  if(abs(tens3(i,j,k)).gt.acc) then
                    print *, i, j, k, tens3(i,j,k)
                  end if
                else
                  if(abs(tensR3(i,j,k)).gt.acc) then
                    print *, i, j, k, tensR3(i,j,k)
                  end if
                end if
              end do
            end do
          end do
          print *, " end print 3D tensor: ", t_name

        else if (rank.eq.4) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'

          print *, "print 4D tensor: ", t_name, " (shape:", dims4, ")"
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          select case (DataKind(1))
          case (C8)
            call c_f_pointer(body_p,tens4,dims4)
          case (R8)
             call c_f_pointer(body_p,tensR4,dims4)
           case default
            stop "wrong datakind in print_tensor"
          end select

          do i = 1, dims4(1)
            do j = 1, dims4(2)
              do k = 1, dims4(3)
                do l = 1, dims4(4)
                  if (DataKind(1)==C8) then
                    if(abs(tens4(i,j,k,l)).gt.acc) then
                      print *, i, j, k, l, tens4(i,j,k,l)
                    end if
                  else
                    if(abs(tensR4(i,j,k,l)).gt.acc) then
                      print *, i, j, k, l, tensR4(i,j,k,l)
                    end if
                  end if
                end do
              end do
            end do
          end do
          print *, " end print 4D tensor: ", t_name

        else
          stop 'error in print_tensor: only ranks 0 to 4 implemented'
        end if

      end subroutine print_tensor
      
      function get_element(h_tensor, dims) result(val)

        !routine to get an element 

        implicit none

        type(talsh_tens_t), intent(inout)    :: h_tensor
        integer(INTD), intent(in)            :: dims(:)
        complex(8)                           :: val

        integer(C_INT) ::  ierr
        integer(INTD)  ::  rank, rrank, ldims
        integer        ::  i, j, k, l
        integer(INTD)  ::  dims1(1)
        integer(INTD)  ::  dims2(1:2)
        integer(INTD)  ::  dims3(1:3)
        integer(INTD)  ::  dims4(1:4)
        type(C_PTR)         :: body_p
        complex(8), pointer :: tens1(:)
        complex(8), pointer :: tens2(:, :)
        complex(8), pointer :: tens3(:, :, :)
        complex(8), pointer :: tens4(:, :, :, :)
        integer             :: DataKind(1)

        DataKind(1)=C8

        rank=talsh_tensor_rank(h_tensor)
        ldims=size(dims)
        if (rank.ne.ldims) stop 'error in get_element: wrong number of indices'

        if (rank.eq.1) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in get_element: wrong rank'

          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens1,dims1)
          i=dims(1)
          val=tens1(i)

        else if (rank.eq.2) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'

          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens2,dims2)
          i=dims(1)
          j=dims(2)
          val=tens2(i,j)

        else if (rank.eq.3) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'

          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens3,dims3)
          i=dims(1)
          j=dims(2)
          k=dims(3)
          val=tens3(i,j,k)

        else if (rank.eq.4) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'

          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens4,dims4)
          i=dims(1)
          j=dims(2)
          k=dims(3)
          l=dims(4)
          val=tens4(i,j,k,l)

        else
          stop 'error in get_element: only ranks 0 to 4 implemented'
        end if

      end function get_element


     subroutine Pair_correlation_energies( t2_amp, int_oovv, nocc, nvir )

        implicit none
        type(talsh_tens_t)  :: t2_amp, int_oovv
        integer(INTD)       :: ierr
        integer             :: nocc
        integer             :: nvir
        integer             :: i, j, a, b
        type(C_PTR)         :: bodya_p, bodyb_p
        integer(INTD)       :: rranka, rrankb
        integer(INTD)       :: dimsa(1:4), dimsb(1:4)
        complex(8), pointer :: tens_int_oovv(:,:,:,:)
        complex(8), pointer :: tens_t2(:,:,:,:)
        complex(8), pointer :: tens_E_pair_corr(:,:)

100 format(*(ES10.2)) ! A scientific format to visualize the possible vanishing values of the matrix E_ij

        allocate(tens_E_pair_corr(nocc,nocc))

        ! Getting the tensor encoding the amplitudes T2(vir_a, vir_b, occ_i, occ_j)
        ierr = talsh_tensor_dimensions(t2_amp,rranka,dimsa)
        ierr = talsh_tensor_get_body_access(t2_amp,bodya_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(bodya_p,tens_t2,dimsa)

        ! Getting the tensor encoding the elements < occ_i, occ_j || vir_a, vir_b >
        ierr = talsh_tensor_dimensions(int_oovv,rrankb,dimsb)
        ierr = talsh_tensor_get_body_access(int_oovv,bodyb_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(bodyb_p,tens_int_oovv,dimsb)
        tens_E_pair_corr = (0.d0,0.d0)
        do i=1,nocc
          do j=1,nocc
            do a=1,nvir
              do b=1,nvir
                tens_E_pair_corr(i,j) = tens_E_pair_corr(i,j)               &
                + 0.5d0 * tens_t2(a,b,i,j) * tens_int_oovv(i,j,a,b)
              enddo
            enddo
          enddo
        enddo

        write(*,*)
        write(*,*) '======================================='
        write(*,*) '  Matrix of pair correlation energies'
        write(*,*) '============================================='
        do i=1,nocc
          write(*,100) real(tens_E_pair_corr(i,:))
        enddo
        write(*,*) "Sum of pair correlation energies : ", sum(real(tens_E_pair_corr))
        write(*,*) "Correlation energy               : ", sum(real(tens_E_pair_corr))/2.D0

        write(*,*) '============================================='
        write(*,*) '======================================='
        write(*,*)

        deallocate(tens_E_pair_corr)

      end subroutine Pair_correlation_energies

      subroutine talsh_init_C8(h_tensor, val)

        !routine for initializing a tensor

        implicit none

        type(talsh_tens_t), intent(inout)    :: h_tensor
        complex(8), intent(in)               :: val

        integer(INTD)  ::  dims1(1)
        integer(INTD)  ::  dims2(1:2)
        integer(INTD)  ::  dims3(1:3)
        integer(INTD)  ::  dims4(1:4)
        integer(INTD)  ::  dims5(1:5)
        integer(INTD)  ::  dims6(1:6)
        integer(C_INT) ::  ierr
        integer(INTD)  ::  rank, rrank
        type(C_PTR)         :: body_p
        complex(8), pointer :: tens1(:)
        complex(8), pointer :: tens2(:, :)
        complex(8), pointer :: tens3(:, :, :)
        complex(8), pointer :: tens4(:, :, :, :)
        complex(8), pointer :: tens5(:, :, :, :, :)
        complex(8), pointer :: tens6(:, :, :, :, :, :)
        integer             :: DataKind(1),nData

        ierr=talsh_tensor_data_kind(h_tensor,nData,DataKind)
        if (ierr.ne.0) stop 'error in getting DataKind'
        if (DataKind(1).ne.C8) stop 'wrong DataKind'

        rank=talsh_tensor_rank(h_tensor)
        if (rank.eq.0) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens1,dims1)
          tens1=val

        else if (rank.eq.1) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens1,dims1)
          tens1=val
        
        else if (rank.eq.2) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens2,dims2)
          tens2=val

        else if (rank.eq.3) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens3,dims3)
          tens3=val

        else if (rank.eq.4) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens4,dims4)
          tens4=val

        else if (rank.eq.5) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims5)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens5,dims5)
          tens5=val

        else if (rank.eq.6) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims6)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens6,dims6)
          tens6=val

        else
          stop 'error in talsh_init: only ranks 0 to 6 implemented'
        end if

      end subroutine talsh_init_C8


      subroutine talsh_init_R8(h_tensor, val)

        !routine for initializing a tensor

        implicit none

        type(talsh_tens_t), intent(inout)    :: h_tensor
        real(8), intent(in)               :: val

        integer(INTD)  ::  dims1(1)
        integer(INTD)  ::  dims2(1:2)
        integer(INTD)  ::  dims3(1:3)
        integer(INTD)  ::  dims4(1:4)
        integer(INTD)  ::  dims5(1:5)
        integer(INTD)  ::  dims6(1:6)
        integer(C_INT) ::  ierr
        integer(INTD)  ::  rank, rrank
        type(C_PTR)         :: body_p
        real(8), pointer :: tens1(:)
        real(8), pointer :: tens2(:, :)
        real(8), pointer :: tens3(:, :, :)
        real(8), pointer :: tens4(:, :, :, :)
        real(8), pointer :: tens5(:, :, :, :, :)
        real(8), pointer :: tens6(:, :, :, :, :, :)
        integer             :: DataKind(1),nData

        ierr=talsh_tensor_data_kind(h_tensor,nData,DataKind)
        if (ierr.ne.0) stop 'error in getting DataKind'
        if (DataKind(1).ne.R8) stop 'wrong DataKind'

        rank=talsh_tensor_rank(h_tensor)
        if (rank.eq.0) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens1,dims1)
          tens1=val

        else if (rank.eq.1) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens1,dims1)
          tens1=val
        
        else if (rank.eq.2) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens2,dims2)
          tens2=val

        else if (rank.eq.3) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens3,dims3)
          tens3=val

        else if (rank.eq.4) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens4,dims4)
          tens4=val

        else if (rank.eq.5) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims5)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens5,dims5)
          tens5=val

        else if (rank.eq.6) then
          ierr = talsh_tensor_dimensions(h_tensor,rrank,dims6)
          if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr=talsh_tensor_get_body_access(h_tensor,body_p,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens6,dims6)
          tens6=val

        else
          stop 'error in talsh_init: only ranks 0 to 6 implemented'
        end if

      end subroutine talsh_init_R8

      function tensor_norm2(a_tensor) result(a_norm2)

         !compute norm2(complete contraction of a tensor with itself) for debugging purposed

         type(talsh_tens_t), intent(inout):: a_tensor

         integer(INTD) :: a_dims(1:4)
         integer(INTD) :: ierr, a_rank

!        scalars (need to be defined as tensor types)
         type(talsh_tens_t)  :: result_tensor
         integer(C_INT)      :: result_dims(1)
         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         real(8) :: a_norm2

!        Initialize scalars that are to be used as tensors in contractions
         result_dims(1) = 1
         ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,result_dims)

         ierr = talsh_tensor_dimensions(a_tensor,a_rank,a_dims)
         select case (a_rank)
         case (0)
            ierr=talsh_tensor_contract("R()+=A+()*A()",result_tensor,a_tensor,a_tensor)
         case (1)
            ierr=talsh_tensor_contract("R()+=A+(p)*A(p)",result_tensor,a_tensor,a_tensor)
         case (2)
            ierr=talsh_tensor_contract("R()+=A+(p,q)*A(p,q)",result_tensor,a_tensor,a_tensor)
         case (3)
            ierr=talsh_tensor_contract("R()+=A+(p,q,r)*A(p,q,r)",result_tensor,a_tensor,a_tensor)
         case (4)
            ierr=talsh_tensor_contract("R()+=A+(p,q,r,s)*A(p,q,r,s)",result_tensor,a_tensor,a_tensor)
         case default
            call quit ('wrong dimension in tensor_norm2')
         end select

         a_norm2 = real(result_tens(1),8)
         ierr=talsh_tensor_destruct(result_tensor)

     end function tensor_norm2

     subroutine laplace_mult (tensor,factor,orb1,orb2,orb3,orb4,orb5,orb6)
      
      type(talsh_tens_t), intent(inout) :: tensor
      real(8), intent(in)               :: factor
      real(8), intent(in)               :: orb1(:),orb2(:),orb3(:)
      real(8), intent(in), optional     :: orb4(:),orb5(:),orb6(:)
      
      integer(INTD)       :: ierr, rank
      integer(INTD)       :: dims3(1:3),dims4(1:4),dims6(1:6)
      integer(INTD)       :: p, q, r, s, t, u
      type(C_PTR)         :: body_p
      complex(8), pointer :: tens3(:, :, :)
      complex(8), pointer :: tens4(:, :, :, :)
      complex(8), pointer :: tens6(:, :, :, :, :, :)
      

      rank=talsh_tensor_rank(tensor)
      if (rank.eq.3) then
         ierr = talsh_tensor_dimensions(tensor,rank,dims3)
         if (size(orb1).ne.dims3(1)) stop 'error in laplace_mult: vector 1 wrong size'
         if (size(orb2).ne.dims3(2)) stop 'error in laplace_mult: vector 2 wrong size'
         if (size(orb3).ne.dims3(3)) stop 'error in laplace_mult: vector 3 wrong size'

         ierr=talsh_tensor_get_body_access(tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,tens3,dims3)

         do r = 1, dims3(3)
           do q = 1, dims3(2)
             do p = 1, dims3(1)
               tens3(p,q,r) = tens3(p,q,r) * EXP(-(orb1(p)+orb2(q)+orb3(r))*factor)
             end do
           end do
         end do

      else if (rank.eq.4) then
         ierr = talsh_tensor_dimensions(tensor,rank,dims4)
         if (size(orb1).ne.dims4(1)) stop 'error in laplace_mult: vector 1 wrong size'
         if (size(orb2).ne.dims4(2)) stop 'error in laplace_mult: vector 2 wrong size'
         if (size(orb3).ne.dims4(3)) stop 'error in laplace_mult: vector 3 wrong size'
         if (size(orb4).ne.dims4(4)) stop 'error in laplace_mult: vector 4 wrong size'

         ierr=talsh_tensor_get_body_access(tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,tens4,dims4)

         do s = 1, dims4(4)
            do r = 1, dims4(3)
               do q = 1, dims4(2)
                  do p = 1, dims4(1)
                     tens4(p,q,r,s) = tens4(p,q,r,s) * EXP(-(orb1(p)+orb2(q)+orb3(r)+orb4(s))*factor)
                  end do
               end do
            end do
         end do

      else if (rank.eq.6) then
         ierr = talsh_tensor_dimensions(tensor,rank,dims6)
         if (size(orb1).ne.dims6(1)) stop 'error in laplace_mult: vector 1 wrong size'
         if (size(orb2).ne.dims6(2)) stop 'error in laplace_mult: vector 2 wrong size'
         if (size(orb3).ne.dims6(3)) stop 'error in laplace_mult: vector 3 wrong size'
         if (size(orb4).ne.dims6(4)) stop 'error in laplace_mult: vector 4 wrong size'
         if (size(orb5).ne.dims6(5)) stop 'error in laplace_mult: vector 5 wrong size'
         if (size(orb6).ne.dims6(6)) stop 'error in laplace_mult: vector 6 wrong size'

         ierr=talsh_tensor_get_body_access(tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,tens6,dims6)

         do u = 1, dims6(6)
           do t = 1, dims6(5)
             do s = 1, dims6(4)
               do r = 1, dims6(3)
                 do q = 1, dims6(2)
                   do p = 1, dims6(1)
                     tens6(p,q,r,s,t,u) = tens6(p,q,r,s,t,u) * &
                     EXP(-(orb1(p)+orb2(q)+orb3(r)+orb4(s)+orb5(t)+orb6(u))*factor)
                   end do
                 end do
               end do
             end do 
           end do
         end do

      else 
        stop 'error in laplace_mult: tensor is not of rank 4/6'
      end if

    end subroutine laplace_mult


    subroutine laplace_testing(test_org, eps_occ, eps_vir)

      type(talsh_tens_t)  , intent(inout)  :: test_org
      real(8), intent(in)                  :: eps_occ(:),eps_vir(:)

      integer(INTD)       :: ierr, rank
      integer(INTD)       :: dims4(1:4), dims6(1:6)
      integer(INTD)       :: nocc, nvir, a, b, c, i, j, k
      type(C_PTR)         :: body_p
      complex(8), pointer :: tens4(:, :, :, :)
      complex(8), pointer :: tens6(:, :, :, :, :, :)
      real(8)             :: test_sum, test_diff

      rank=talsh_tensor_rank(test_org)
      if (rank.eq.4) then
         ierr = talsh_tensor_dimensions(test_org,rank,dims4)
         if (ierr.ne.0 .or. rank.ne.4) stop 'error: tensor corrupted'

         if (dims4(1).eq.dims4(2)) then
           nvir = dims4(1)
         else
           stop 'error: asymmetric tensor in laplace_testing'
         end if
         if (dims4(3).eq.dims4(4)) then
           nocc = dims4(3)
         else
           stop 'error: asymmetric tensor in laplace_testing'
         end if

         ierr=talsh_tensor_get_body_access(test_org,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,tens4,dims4)


         test_sum=0.0D0
         test_diff=0.0D0
         do j = 1, nocc
           do i = 1, nocc
             do b = 1, nvir
               do a = 1, nvir
                test_sum=test_sum+tens4(a,b,i,j)
                test_diff=test_diff+tens4(a,b,i,j) - &
                     1.0D0/(eps_vir(a)+eps_vir(b)-eps_occ(i)-eps_occ(j))
               end do
             end do
           end do
         end do

         write(*,*) 'result of testing all elements:'
         write(*,*) 'test sum   = ',test_sum
         write(*,*) 'test error = ',test_diff

      else if (rank.eq.6) then
         ierr = talsh_tensor_dimensions(test_org,rank,dims6)
         if (ierr.ne.0 .or. rank.ne.6) stop 'error: tensor corrupted'

         if (dims6(1).eq.dims6(2) .and. dims6(1).eq.dims6(3)) then
           nvir = dims6(1)
         else
           stop 'error: asymmetric tensor in laplace_testing'
         end if
         if (dims6(4).eq.dims6(5) .and. dims6(4).eq.dims6(6)) then
           nocc = dims6(4)
         else
           stop 'error: asymmetric tensor in laplace_testing'
         end if

         ierr=talsh_tensor_get_body_access(test_org,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,tens6,dims6)

         test_sum=0.0D0
         test_diff=0.0D0
         do c = 1, nvir
          do b = 1, nvir
           do a = 1, nvir
            do k = 1, nocc
             do j = 1, nocc
              do i = 1, nocc
               test_sum=test_sum+tens6(a,b,c,i,j,k)
               test_diff=test_diff+tens6(a,b,c,i,j,k) &
                   - 1.0D0/(eps_vir(a)+eps_vir(b)+eps_vir(c)-eps_occ(i)-eps_occ(j)-eps_occ(k))
              end do
             end do   
            end do      
           end do
          end do
         end do

         write(*,*) 'result of testing all elements:'
         write(*,*) 'test sum   = ',test_sum
         write(*,*) 'test error = ',test_diff

      else
        stop 'error in laplace_testing: only rank 4 or 6'
      end if

    end subroutine laplace_testing


end module talsh_common_routines
