module eti_common_routines

!This module contains subroutines both needed in solving t equations as well as Lambda equations.
!get_orbital_energies,diis,scale_with_denominators, laplace multiplications

    use exacorr_tensor_interface
    use tensor_algebra
    use, intrinsic:: ISO_C_BINDING
    use exacorr_tensor_interface
    use exacorr_scaling_parameters

    implicit none

    private

    interface scale_with_denominators
        module procedure scale_with_denominators
    end interface scale_with_denominators

    public get_orbital_energies
    public scale_with_denominators
    public print_tensor
    public write_2dtensor_as_matrix
    public get_element
    public Pair_correlation_energies
    public eti_init_R8
    public eti_init_C8
    public tensor_norm2
    public laplace_mult
    public laplace_testing
    public print_tensor_elements

    contains

      subroutine scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor,t2_tensor,t3_tensor, &
                                          full_t3_tensor,full_t4_tensor,full_t5_tensor,eps_ijk)

       type(tensor_t), optional, intent(inout)  :: t1_tensor, t2_tensor, t3_tensor
       type(tensor_t), optional, intent(inout)  :: full_t3_tensor
       type(tensor_t), optional, intent(inout)  :: full_t4_tensor
       type(tensor_t), optional, intent(inout)  :: full_t5_tensor
       real(8), intent(in)                          :: eps_occ(:),eps_vir(:)
       real(8), intent(in), optional                :: eps_ijk

       real(8)  :: denominator
       complex(8), pointer, contiguous:: t1_tens(:,:), t2_tens(:,:,:,:), t3_tens(:,:,:)
       complex(8), pointer, contiguous:: full_t3_tens(:,:,:,:,:,:)
#ifndef VAR_NVFORTRAN 
       complex(8), pointer, contiguous:: full_t4_tens(:,:,:,:,:,:,:,:)
       complex(8), pointer, contiguous:: full_t5_tens(:,:,:,:,:,:,:,:,:,:)
#endif
       type(C_PTR):: body_p

       integer(INTD) :: t1_dims(1:2), t2_dims(1:4), t3_dims(1:3)
       integer(INTD) :: full_t3_dims(1:6)
       integer(INTD) :: full_t4_dims(1:8)
       integer(INTD) :: full_t5_dims(1:10)
       integer(INTD) :: ierr, tens_rank
       integer       :: a, b, c, d, e, i, j, k, l, m

       integer, optional, intent(in) :: nocc

       if (present(t1_tensor)) then

          ierr = get_nmodes(t1_tensor,tens_rank)
          ierr = get_extents_static(t1_tensor,t1_dims)
          if (ierr.ne.0 .or. tens_rank.ne.2) stop 'program error: t1 tensor corrupted'
          ierr=direct_data_pointer(t1_tensor,body_p)
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

          ierr = get_nmodes(t2_tensor,tens_rank)
          ierr = get_extents_static(t2_tensor,t2_dims)
          if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: t2 tensor corrupted'
          ierr=direct_data_pointer(t2_tensor,body_p)
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

       if (present(full_t3_tensor)) then

          ierr = get_nmodes(full_t3_tensor,tens_rank)
          ierr = get_extents_static(full_t3_tensor,full_t3_dims)
          if (ierr.ne.0 .or. tens_rank.ne.6) stop 'program error: full t3 tensor corrupted' 
          ierr=direct_data_pointer(full_t3_tensor,body_p)
          call c_f_pointer(body_p,full_t3_tens,full_t3_dims(1:6)) 

          if (full_t3_dims(4) .eq. nocc) then
            !t amplitudes
            do k = 1, full_t3_dims(6)
              do j = 1, full_t3_dims(5)
                do i = 1, full_t3_dims(4)
                  do c = 1, full_t3_dims(3)
                    do b = 1, full_t3_dims(2)
                      do a = 1, full_t3_dims(1)
                        denominator = eps_occ(i) + eps_occ(j) + eps_occ(k) 
                        denominator = denominator - eps_vir(a) - eps_vir(b) - eps_vir(c)
                        full_t3_tens(a,b,c,i,j,k) = full_t3_tens(a,b,c,i,j,k) / denominator
                      end do
                    end do
                  end do
                end do
              end do
            end do

          else
            !Lambda
            do i = 1, full_t3_dims(1)
              do j = 1, full_t3_dims(2)
                do k = 1, full_t3_dims(3)
                  do a = 1, full_t3_dims(4)
                    do b = 1, full_t3_dims(5)
                      do c = 1, full_t3_dims(6)
                        denominator = eps_occ(i) + eps_occ(j) + eps_occ(k) 
                        denominator = denominator - eps_vir(a) - eps_vir(b) - eps_vir(c)
                        full_t3_tens(i,j,k,a,b,c) = full_t3_tens(i,j,k,a,b,c) / denominator
                      end do
                    end do
                  end do
                end do
              end do
            end do

          end if

       end if

#ifndef VAR_NVFORTRAN
       if (present(full_t4_tensor)) then

          ierr = get_nmodes(full_t4_tensor,tens_rank)
          ierr = get_extents_static(full_t4_tensor,full_t4_dims)
          if (ierr.ne.0 .or. tens_rank.ne.8) stop 'program error: full t4 tensor corrupted'
          ierr=direct_data_pointer(full_t4_tensor,body_p)
          call c_f_pointer(body_p,full_t4_tens,full_t4_dims(1:8))

          if (full_t4_dims(5) .eq. nocc) then
            !t amplitudes
            do l = 1, full_t4_dims(8)
              do k = 1, full_t4_dims(7)
                do j = 1, full_t4_dims(6)
                  do i = 1, full_t4_dims(5)
                    do d = 1, full_t4_dims(4)
                      do c = 1, full_t4_dims(3)
                        do b = 1, full_t4_dims(2)
                          do a = 1, full_t4_dims(1)
                            denominator = eps_occ(i) + eps_occ(j) + eps_occ(k) + eps_occ(l) 
                            denominator = denominator - eps_vir(a) - eps_vir(b) - eps_vir(c) - eps_vir(d)
                            full_t4_tens(a,b,c,d,i,j,k,l) = full_t4_tens(a,b,c,d,i,j,k,l) / denominator
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

          else
            !Lambda
            do i = 1, full_t4_dims(1)
              do j = 1, full_t4_dims(2)
                do k = 1, full_t4_dims(3)
                  do l = 1, full_t4_dims(4)
                    do a = 1, full_t4_dims(5)
                      do b = 1, full_t4_dims(6)
                        do c = 1, full_t4_dims(7)
                          do d = 1, full_t4_dims(8)
                            denominator = eps_occ(i) + eps_occ(j) + eps_occ(k) + eps_occ(l) 
                            denominator = denominator - eps_vir(a) - eps_vir(b) - eps_vir(c) - eps_vir(d)
                            full_t4_tens(i,j,k,l,a,b,c,d) = full_t4_tens(i,j,k,l,a,b,c,d) / denominator
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

          end if

       end if

       if (present(full_t5_tensor)) then

          ierr = get_nmodes(full_t5_tensor,tens_rank)
          ierr = get_extents_static(full_t5_tensor,full_t5_dims)
          if (ierr.ne.0 .or. tens_rank.ne.10) stop 'program error: full t5 tensor corrupted'
          ierr=direct_data_pointer(full_t5_tensor,body_p)
          call c_f_pointer(body_p,full_t5_tens,full_t5_dims(1:10))

          if (full_t5_dims(6) .eq. nocc) then
            !t amplitudes
            do m = 1, full_t5_dims(10)
              do l = 1, full_t5_dims(9)
                do k = 1, full_t5_dims(8)
                  do j = 1, full_t5_dims(7)
                    do i = 1, full_t5_dims(6)
                      do e = 1, full_t5_dims(5)
                        do d = 1, full_t5_dims(4)
                          do c = 1, full_t5_dims(3)
                            do b = 1, full_t5_dims(2)
                              do a = 1, full_t5_dims(1)
                                denominator = eps_occ(i) + eps_occ(j) + eps_occ(k) + eps_occ(l) + eps_occ(m)
                                denominator = denominator &
                                            - eps_vir(a) - eps_vir(b) - eps_vir(c) - eps_vir(d) - eps_vir(e)
                                full_t5_tens(a,b,c,d,e,i,j,k,l,m) = full_t5_tens(a,b,c,d,e,i,j,k,l,m) / denominator
                              end do
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

          else
            !Lambda
            do i = 1, full_t5_dims(1)
              do j = 1, full_t5_dims(2)
                do k = 1, full_t5_dims(3)
                  do l = 1, full_t5_dims(4)
                    do m = 1, full_t5_dims(5)
                      do a = 1, full_t5_dims(6)
                        do b = 1, full_t5_dims(7)
                          do c = 1, full_t5_dims(8)
                            do d = 1, full_t5_dims(9)
                              do e = 1, full_t5_dims(10)
                                denominator = eps_occ(i) + eps_occ(j) + eps_occ(k) + eps_occ(l) + eps_occ(m)
                                denominator = denominator & 
                                            - eps_vir(a) - eps_vir(b) - eps_vir(c) - eps_vir(d) - eps_vir(e)
                                full_t5_tens(i,j,k,l,m,a,b,c,d,e) = full_t5_tens(i,j,k,l,m,a,b,c,d,e) / denominator
                              end do
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

          end if

       end if ! full_t5_tensor
#else
       if (present(full_t4_tensor) .or. present(full_t5_tensor)) then
          stop "T4/T5 unsupported due to NVFORTRAN restrictions on array dimensions"
       end if
#endif

       if (present(t3_tensor)) then

          ierr = get_nmodes(t3_tensor,tens_rank)
          ierr = get_extents_static(t3_tensor,t3_dims)
          if (ierr.ne.0 .or. tens_rank.ne.3) stop 'program error: t3 tensor corrupted'
          ierr=direct_data_pointer(t3_tensor,body_p)
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

        !routine for printing tensors

        implicit none

        type(tensor_t), intent(inout)    :: h_tensor
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

        ierr=get_data_kind_tensor_t(h_tensor, DataKind(1))
        if (ierr.ne.0) stop 'error in getting DataKind'

        ierr = get_nmodes(h_tensor,rank)
        if (rank.eq.0) then
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor, dims1)

          print *, "print 0D tensor: ", t_name, " (shape:", dims1,")"
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
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
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor, dims1)

          print *, "print 1D tensor: ", t_name, " (shape:", dims1,")"
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
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
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
!        if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor, dims2)

          print *, "print 2D tensor: ", t_name, " (shape:", dims2,")"
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
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
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor, dims3)

          print *, "print 3D tensor: ", t_name, " (shape:", dims3,")"
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
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
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor, dims4)

          print *, "print 4D tensor: ", t_name, " (shape:", dims4, ")"
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
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

      !------------------------------------------------------------------------------------------------------------------------------------
      subroutine write_2dtensor_as_matrix(nvir,dm_vv,file_out,mode)
! Routine for writing complex density matrix into the file'DM_complex'. 
! Written by Xiang Yuan
! Modified by Johann Pototschnig and shifted to talsh_common_routines
! andre: todo, ! generalize as a reshape to a 2d array to write to file ? 

       use exacorr_utils

! input variables
         integer,intent(in)                  :: nvir
         type(tensor_t), intent(inout)   :: dm_vv
         character(len=*),intent(in)         :: file_out
         integer, intent(in), optional       :: mode ! 1:binary, 2:ascii

! body pointers
         type(C_PTR)                     :: body_p2
         complex(8), pointer, contiguous :: dm_block(:,:)
         complex(8), allocatable         :: dm_mo_fort(:,:)
! error code
         integer :: ierr
         integer :: mp2dm

         allocate (dm_mo_fort(nvir,nvir))
         dm_mo_fort = ZERO
         
         ierr=direct_data_pointer(dm_vv,body_p2)
         call c_f_pointer(body_p2,dm_block,(/nvir,nvir/)) ! to use <dm_mo_block> as a regular Fortran 2d array
         dm_mo_fort(1:nvir,1:nvir) = dm_block
          
         call get_free_fileunit(mp2dm)
         if (present(mode)) then
           if (mode==2) then
              open (mp2dm,file=file_out,status='REPLACE',FORM='FORMATTED',access='SEQUENTIAL')
              ! write dimension of the matrix 
              write(mp2dm,*) nvir
              ! write matrix element
              do ierr=1,nvir
                write(mp2dm,*) dm_mo_fort(ierr,1:nvir)
              end do
            else
               call quit('error in choice of mode for writing mp2nos')
            end if
         else 
           open (mp2dm,file=file_out,status='REPLACE',FORM='UNFORMATTED',access='SEQUENTIAL')
           ! write dimension of the matrix 
           write(mp2dm) nvir
           ! write matrix element 
           write(mp2dm) dm_mo_fort
         end if

         close(mp2dm,status='KEEP')

        end subroutine write_2dtensor_as_matrix
      
      function get_element(h_tensor, dims) result(val)

        !routine to get an element 

        implicit none

        type(tensor_t), intent(inout)    :: h_tensor
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

        ierr = get_nmodes(h_tensor,rank)
        ldims=size(dims)
        if (rank.ne.ldims) stop 'error in get_element: wrong number of indices'

        if (rank.eq.1) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in get_element: wrong rank'
          ierr = get_extents_static(h_tensor,dims1)

          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens1,dims1)
          i=dims(1)
          val=tens1(i)

        else if (rank.eq.2) then
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims2)

          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens2,dims2)
          i=dims(1)
          j=dims(2)
          val=tens2(i,j)

        else if (rank.eq.3) then
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims3)

          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens3,dims3)
          i=dims(1)
          j=dims(2)
          k=dims(3)
          val=tens3(i,j,k)

        else if (rank.eq.4) then
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims4)

          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
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
        type(tensor_t)  :: t2_amp, int_oovv
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
        ierr = get_nmodes(t2_amp,rranka)
        ierr = get_extents_static(t2_amp,dimsa)
        ierr = direct_data_pointer(t2_amp,bodya_p)
        call c_f_pointer(bodya_p,tens_t2,dimsa)

        ! Getting the tensor encoding the elements < occ_i, occ_j || vir_a, vir_b >
        ierr = get_nmodes(int_oovv,rrankb)
        ierr = get_extents_static(int_oovv,dimsb)
        ierr = direct_data_pointer(int_oovv,bodyb_p)
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

      subroutine eti_init_C8(h_tensor, val)

        !routine for initializing a tensor

        implicit none

        type(tensor_t), intent(inout)    :: h_tensor
        complex(8), intent(in)               :: val

        integer(INTD)  ::  dims1(1)
        integer(INTD)  ::  dims2(1:2)
        integer(INTD)  ::  dims3(1:3)
        integer(INTD)  ::  dims4(1:4)
        integer(INTD)  ::  dims5(1:5)
        integer(INTD)  ::  dims6(1:6)
        integer(INTD)  ::  dims7(1:7)
        integer(INTD)  ::  dims8(1:8)
        integer(INTD)  ::  dims9(1:9)
        integer(INTD)  ::  dims10(1:10)
        integer(C_INT) ::  ierr
        integer(INTD)  ::  rank, rrank
        type(C_PTR)         :: body_p
        complex(8), pointer :: tens1(:)
        complex(8), pointer :: tens2(:, :)
        complex(8), pointer :: tens3(:, :, :)
        complex(8), pointer :: tens4(:, :, :, :)
        complex(8), pointer :: tens5(:, :, :, :, :)
        complex(8), pointer :: tens6(:, :, :, :, :, :)
        complex(8), pointer :: tens7(:, :, :, :, :, :, :)
#ifndef VAR_NVFORTRAN
        complex(8), pointer :: tens8(:, :, :, :, :, :, :, :)
        complex(8), pointer :: tens9(:, :, :, :, :, :, :, :, :)
        complex(8), pointer :: tens10(:, :, :, :, :,:, :, :, :, :)
#endif
        integer             :: DataKind(1),nData

        ierr=get_data_kind_tensor_t(h_tensor, DataKind(1))
        if (ierr.ne.0) stop 'error in getting DataKind'
        if (DataKind(1).ne.C8) stop 'wrong DataKind'

        ierr = get_nmodes(h_tensor,rank)
        if (rank.eq.0) then
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr = get_extents_static(h_tensor,dims1)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens1,dims1)
          tens1=val

        else if (rank.eq.1) then
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr = get_extents_static(h_tensor,dims1)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens1,dims1)
          tens1=val
        
        else if (rank.eq.2) then
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr = get_extents_static(h_tensor,dims2)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens2,dims2)
          tens2=val

        else if (rank.eq.3) then
! does this make sense ?
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr = get_extents_static(h_tensor,dims3)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens3,dims3)
          tens3=val

        else if (rank.eq.4) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims4)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens4,dims4)
          tens4=val

        else if (rank.eq.5) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims5)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr = get_extents_static(h_tensor,dims5)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens5,dims5)
          tens5=val

        else if (rank.eq.6) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims6)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims6)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens6,dims6)
          tens6=val

        else if (rank.eq.7) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims7)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims7)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens7,dims7)
          tens7=val

#ifndef VAR_NVFORTRAN
        else if (rank.eq.8) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims8)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims8)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens8,dims8)
          tens8=val

        else if (rank.eq.9) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims9)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims9)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens9,dims9)
          tens9=val

        else if (rank.eq.10) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims10)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims10)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens10,dims10)
          tens10=val
#endif 
        else
          stop 'error in eti_init: only ranks 0 to 10 implemented (0 to 7 for NVFORTRAN)'
        end if

      end subroutine eti_init_C8


      subroutine eti_init_R8(h_tensor, val)

        !routine for initializing a tensor

        implicit none

        type(tensor_t), intent(inout)    :: h_tensor
        real(8), intent(in)               :: val

        integer(INTD)  ::  dims1(1)
        integer(INTD)  ::  dims2(1:2)
        integer(INTD)  ::  dims3(1:3)
        integer(INTD)  ::  dims4(1:4)
        integer(INTD)  ::  dims5(1:5)
        integer(INTD)  ::  dims6(1:6)
        integer(INTD)  ::  dims7(1:7)
        integer(INTD)  ::  dims8(1:8)
        integer(INTD)  ::  dims9(1:9)
        integer(INTD)  ::  dims10(1:10)
        integer(C_INT) ::  ierr
        integer(INTD)  ::  rank, rrank
        type(C_PTR)         :: body_p
        real(8), pointer :: tens1(:)
        real(8), pointer :: tens2(:, :)
        real(8), pointer :: tens3(:, :, :)
        real(8), pointer :: tens4(:, :, :, :)
        real(8), pointer :: tens5(:, :, :, :, :)
        real(8), pointer :: tens6(:, :, :, :, :, :)
        real(8), pointer :: tens7(:, :, :, :, :, :, :)
#ifndef VAR_NVFORTRAN
        real(8), pointer :: tens8(:, :, :, :, :, :, :, :)
        real(8), pointer :: tens9(:, :, :, :, :, :, :, :, :)
        real(8), pointer :: tens10(:, :, :, :, :,:, :, :, :, :)
#endif
        integer             :: DataKind(1),nData

        ierr=get_data_kind_tensor_t(h_tensor, DataKind(1))
        if (ierr.ne.0) stop 'error in getting DataKind'
        if (DataKind(1).ne.R8) stop 'wrong DataKind'

        ierr = get_nmodes(h_tensor,rank)
        if (rank.eq.0) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr = get_extents_static(h_tensor,dims1)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens1,dims1)
          tens1=val

        else if (rank.eq.1) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims1)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr = get_extents_static(h_tensor,dims1)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens1,dims1)
          tens1=val
        
        else if (rank.eq.2) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr = get_extents_static(h_tensor,dims2)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens2,dims2)
          tens2=val

        else if (rank.eq.3) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr = get_extents_static(h_tensor,dims3)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens3,dims3)
          tens3=val

        else if (rank.eq.4) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims4)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens4,dims4)
          tens4=val

        else if (rank.eq.5) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims5)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in talsh_init: wrong rank'
          ierr = get_extents_static(h_tensor,dims5)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens5,dims5)
          tens5=val

        else if (rank.eq.6) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims6)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims6)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens6,dims6)
          tens6=val

        else if (rank.eq.7) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims7)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims7)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens7,dims7)
          tens7=val

#ifndef VAR_NVFORTRAN
        else if (rank.eq.8) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims8)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims8)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens8,dims8)
          tens8=val

        else if (rank.eq.9) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims9)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims9)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens9,dims9)
          tens9=val

        else if (rank.eq.10) then
!         ierr = talsh_tensor_dimensions(h_tensor,rrank,dims10)
!         if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in print_tensor: wrong rank'
          ierr = get_extents_static(h_tensor,dims10)
          ierr=direct_data_pointer(h_tensor,body_p) ! ,DataKind(1),int(0,C_INT),DEV_HOST)
          call c_f_pointer(body_p,tens10,dims10)
          tens10=val

#endif
        else
          stop 'error in eti_init: only ranks 0 to 10 implemented (0 to 7 for NVFORTRAN)'
        end if

      end subroutine eti_init_R8

      function tensor_norm2(a_tensor) result(a_norm2)

         !compute norm2(complete contraction of a tensor with itself) for debugging purposed

         type(tensor_t), intent(inout):: a_tensor

         integer(INTD), allocatable :: a_dims(:)
         integer(INTD) :: ierr, a_rank

!        scalars (need to be defined as tensor types)
         type(tensor_t)  :: result_tensor
         integer(C_INT)      :: result_dims(1)
         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         real(8) :: a_norm2

!        Initialize scalars that are to be used as tensors in contractions
         result_dims(1) = 1
         ierr=construct_tensor_t(result_tensor,C8,result_dims(1:0),init_val=ZERO)

         ierr = get_nmodes(a_tensor,a_rank)
         ierr = get_extents(a_tensor,a_dims)
         select case (a_rank)
         case (0)
            ierr=t_einsum("R()+=A+()*A()",result_tensor,a_tensor,a_tensor)
         case (1)
            ierr=t_einsum("R()+=A+(p)*A(p)",result_tensor,a_tensor,a_tensor)
         case (2)
            ierr=t_einsum("R()+=A+(p,q)*A(p,q)",result_tensor,a_tensor,a_tensor)
         case (3)
            ierr=t_einsum("R()+=A+(p,q,r)*A(p,q,r)",result_tensor,a_tensor,a_tensor)
         case (4)
            ierr=t_einsum("R()+=A+(p,q,r,s)*A(p,q,r,s)",result_tensor,a_tensor,a_tensor)
         case (5)
            ierr=t_einsum("R()+=A+(p,q,r,s,t)*A(p,q,r,s,t)",result_tensor,a_tensor,a_tensor)
         case (6)
            ierr=t_einsum("R()+=A+(p,q,r,s,t,u)*A(p,q,r,s,t,u)",result_tensor,a_tensor,a_tensor)
         case (7)
            ierr=t_einsum("R()+=A+(p,q,r,s,t,u,v)*A(p,q,r,s,t,u,v)",result_tensor,a_tensor,a_tensor)
         case (8)
            ierr=t_einsum("R()+=A+(p,q,r,s,t,u,v,x)*A(p,q,r,s,t,u,v,x)",result_tensor,a_tensor,a_tensor)
         case (9)
            ierr=t_einsum("R()+=A+(a,p,q,r,s,t,u,v,x)*A(a,p,q,r,s,t,u,v,x)",result_tensor,a_tensor,a_tensor)
         case (10)
            ierr=t_einsum("R()+=A+(a,b,p,q,r,s,t,u,v,x)*A(a,b,p,q,r,s,t,u,v,x)",result_tensor,a_tensor,a_tensor)
         case default
            call quit ('wrong dimension in tensor_norm2')
         end select

         ierr=direct_data_pointer(result_tensor,body_p)
         call c_f_pointer(body_p,result_tens,result_dims)
         a_norm2 = real(result_tens(1),8)
         ierr=destroy_tensor_t(result_tensor)

     end function tensor_norm2


     subroutine laplace_mult (tensor,factor,orb1,orb2,orb3,orb4,orb5,orb6)
      
      type(tensor_t), intent(inout) :: tensor
      real(8), intent(in)               :: factor
      real(8), intent(in)               :: orb1(:),orb2(:)
      real(8), intent(in), optional     :: orb3(:),orb4(:),orb5(:),orb6(:)
      
      integer(INTD)       :: ierr, rank
      integer(INTD)       :: dims2(1:2),dims3(1:3),dims4(1:4),dims6(1:6)
      integer(INTD)       :: p, q, r, s, t, u
      type(C_PTR)         :: body_p
      complex(8), pointer :: tens2(:, :)
      complex(8), pointer :: tens3(:, :, :)
      complex(8), pointer :: tens4(:, :, :, :)
      complex(8), pointer :: tens6(:, :, :, :, :, :)

      ierr = get_nmodes(tensor,rank)
      if (rank.eq.2) then
         ierr = get_extents_static(tensor,dims2)
         if (size(orb1).ne.dims2(1)) stop 'error in laplace_mult: vector 1 wrong size'
         if (size(orb2).ne.dims2(2)) stop 'error in laplace_mult: vector 2 wrong size'

         ierr=direct_data_pointer(tensor,body_p)
         call c_f_pointer(body_p,tens2,dims2)

         do r = 1, dims2(2)
           do p = 1, dims2(1)
               tens2(p,r) = tens2(p,r) * DEXP(-(orb1(p)+orb2(r))*factor)
           end do
         end do

      else if (rank.eq.3) then
         ierr = get_extents_static(tensor,dims3)
         if (size(orb1).ne.dims3(1)) stop 'error in laplace_mult: vector 1 wrong size'
         if (size(orb2).ne.dims3(2)) stop 'error in laplace_mult: vector 2 wrong size'
         if (size(orb3).ne.dims3(3)) stop 'error in laplace_mult: vector 3 wrong size'

         ierr=direct_data_pointer(tensor,body_p)
         call c_f_pointer(body_p,tens3,dims3)

         do r = 1, dims3(3)
           do q = 1, dims3(2)
             do p = 1, dims3(1)
               tens3(p,q,r) = tens3(p,q,r) * DEXP(-(orb1(p)+orb2(q)+orb3(r))*factor)
             end do
           end do
         end do

      else if (rank.eq.4) then
         ierr = get_extents_static(tensor,dims4)
         if (size(orb1).ne.dims4(1)) stop 'error in laplace_mult: vector 1 wrong size'
         if (size(orb2).ne.dims4(2)) stop 'error in laplace_mult: vector 2 wrong size'
         if (size(orb3).ne.dims4(3)) stop 'error in laplace_mult: vector 3 wrong size'
         if (size(orb4).ne.dims4(4)) stop 'error in laplace_mult: vector 4 wrong size'

         ierr=direct_data_pointer(tensor,body_p)
         call c_f_pointer(body_p,tens4,dims4)

         do s = 1, dims4(4)
            do r = 1, dims4(3)
               do q = 1, dims4(2)
                  do p = 1, dims4(1)
                     tens4(p,q,r,s) = tens4(p,q,r,s) * DEXP(-(orb1(p)+orb2(q)+orb3(r)+orb4(s))*factor)
                  end do
               end do
            end do
         end do

      else if (rank.eq.6) then
         ierr = get_extents_static(tensor,dims6)
         if (size(orb1).ne.dims6(1)) stop 'error in laplace_mult: vector 1 wrong size'
         if (size(orb2).ne.dims6(2)) stop 'error in laplace_mult: vector 2 wrong size'
         if (size(orb3).ne.dims6(3)) stop 'error in laplace_mult: vector 3 wrong size'
         if (size(orb4).ne.dims6(4)) stop 'error in laplace_mult: vector 4 wrong size'
         if (size(orb5).ne.dims6(5)) stop 'error in laplace_mult: vector 5 wrong size'
         if (size(orb6).ne.dims6(6)) stop 'error in laplace_mult: vector 6 wrong size'

         ierr=direct_data_pointer(tensor,body_p)
         call c_f_pointer(body_p,tens6,dims6)

         do u = 1, dims6(6)
           do t = 1, dims6(5)
             do s = 1, dims6(4)
               do r = 1, dims6(3)
                 do q = 1, dims6(2)
                   do p = 1, dims6(1)
                     tens6(p,q,r,s,t,u) = tens6(p,q,r,s,t,u) * &
                     DEXP(-(orb1(p)+orb2(q)+orb3(r)+orb4(s)+orb5(t)+orb6(u))*factor)
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

      type(tensor_t)  , intent(inout)  :: test_org
      real(8), intent(in)                  :: eps_occ(:),eps_vir(:)

      integer(INTD)       :: ierr, rank
      integer(INTD)       :: dims4(1:4), dims6(1:6)
      integer(INTD)       :: nocc, nvir, a, b, c, i, j, k
      type(C_PTR)         :: body_p
      complex(8), pointer :: tens4(:, :, :, :)
      complex(8), pointer :: tens6(:, :, :, :, :, :)
      real(8)             :: test_sum, test_diff

      ierr = get_nmodes(test_org,rank)
      if (rank.eq.4) then
         ierr = get_extents_static(test_org,dims4)
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

         ierr=direct_data_pointer(test_org,body_p)
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
         ierr = get_extents_static(test_org,dims6)
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

         ierr=direct_data_pointer(test_org,body_p)
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

    subroutine print_tensor_elements(tensor, relccsd_order)
       implicit none

       type(tensor_t), intent(inout)    :: tensor
       logical, intent(in), optional :: relccsd_order

       integer(INTD)  ::  rank
       integer(INTD)  ::  dims1(1)
       integer(INTD)  ::  dims2(1:2)
       integer(INTD)  ::  dims3(1:3)
       integer(INTD)  ::  dims4(1:4)
       integer(INTD)  ::  ierr
       integer  ::  i, j, k, l, b
       integer, pointer :: index_i_exacorr_to_relccsd(:)
       integer, pointer :: index_j_exacorr_to_relccsd(:)
       integer :: ku, kb, num_kr_pairs_i, num_kr_pairs_j, num_spinors_i, num_spinors_j
       type(C_PTR)         :: body
       complex(8), pointer :: tens1(:)
       complex(8), pointer :: tens2(:,:)
       complex(8) :: value_to_print
       logical :: print_in_relccsd_order

       ierr = get_nmodes(tensor,rank)

       if (rank.eq.1) then
         ierr = get_extents_static(tensor,dims1)
         ierr=direct_data_pointer(tensor,body)
         call c_f_pointer(body,tens1,dims1)

         b = 0
         do i = 1, dims1(1)
            b = b + 1
            value_to_print = tens1(i)
            write(6,'(3X,2I4,2x,2E28.14)') b, value_to_print
         end do

       else if (rank.eq.2) then
         ierr = get_extents_static(tensor,dims2)
         ierr=direct_data_pointer(tensor,body)
         call c_f_pointer(body,tens2,dims2)

         if (present(relccsd_order)) then
            print_in_relccsd_order = relccsd_order
         else
            print_in_relccsd_order = .false.
         end if

         if (print_in_relccsd_order) then
            num_spinors_i = dims2(1)
            num_spinors_j = dims2(2)

            num_kr_pairs_i= num_spinors_i / 2
            num_kr_pairs_j= num_spinors_j / 2

            allocate(index_i_exacorr_to_relccsd(num_spinors_i))
            allocate(index_j_exacorr_to_relccsd(num_spinors_j))
            index_i_exacorr_to_relccsd = 0
            index_j_exacorr_to_relccsd = 0

            ku = 1
            kb = 1
            do i = 1, num_spinors_i
               if (mod(i,2).eq.0) then
                  index_i_exacorr_to_relccsd(num_kr_pairs_i + kb) = i
                  kb = kb + 1
               else
                  index_i_exacorr_to_relccsd(ku) = i
                  ku = ku + 1
               end if
            end do

            ku = 1
            kb = 1
            do j = 1, num_spinors_j
               if (mod(j,2).eq.0) then
                  index_j_exacorr_to_relccsd(num_kr_pairs_j + kb) = j
                  kb = kb + 1
               else
                  index_j_exacorr_to_relccsd(ku) = j
                  ku = ku + 1
               end if
            end do
            print *,'index_j_exacorr_to_relccsd:',index_j_exacorr_to_relccsd

            b = 0
            do j = 1, dims2(2)
               do i = 1, dims2(1)
                  b = b + 1
                  value_to_print = tens2(index_i_exacorr_to_relccsd(i),index_j_exacorr_to_relccsd(j))
                  write(6,'(3X,I4,2x,2F28.10,4I4)') b, value_to_print, i, j, &
          &            index_i_exacorr_to_relccsd(i), index_j_exacorr_to_relccsd(j)
               end do
            end do

            deallocate(index_i_exacorr_to_relccsd)
            deallocate(index_j_exacorr_to_relccsd)

          else
            b = 0
            do j = 1, dims2(2)
               do i = 1, dims2(1)
                  b = b + 1
                  value_to_print = tens2(i,j)
                  write(6,'(3X,I4,2x,2F28.10,2I4)') b, value_to_print, i, j
               end do
            end do
          end if

       else
          call quit('unsupported tensor rank in print_tensor_elements')
       end if

    end subroutine




end module eti_common_routines
