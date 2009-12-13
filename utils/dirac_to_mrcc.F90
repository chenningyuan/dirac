module dirac_to_mrcc
  implicit none

  integer, parameter     :: filenumber_1e = 21
  integer, parameter     :: filenumber_2e = 22
  integer, parameter     :: filenumber_55 = 55
  logical, parameter     :: generate_full_list = .true.

  type SpinorInformation

      integer     :: irrep
      integer     :: abelian_irrep
      integer     :: occupation
      integer     :: index
      real(8)     :: energy

  end type SpinorInformation

  integer                 :: number_of_spinors, number_of_electrons, number_of_irreps, number_of_abelian_irreps
  real(8)                 :: core_energy
  integer, allocatable    :: kramer_to_spinor(:)
  integer, allocatable    :: multiplication_table(:,:)
  integer                 :: irrep_occupation(128)
  type(SpinorInformation), allocatable :: spinor(:)

  public initialize, write_fort55
  private process_1e, process_2e, irrep_reordered

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize ()

     integer       ::  number_of_kramerspairs, i, j, index, n
     character(10) ::  date_of_generation
     character(8)  ::  time_of_generation
     logical       ::  breit_included
     character(14) ::  representation_name

     open (filenumber_1e, file='MRCONEE', Form='UNFORMATTED')
     read (filenumber_1e) number_of_spinors, &
                          breit_included,    &
                          core_energy
     read (filenumber_1e) number_of_irreps,(representation_name,i=1,number_of_irreps), &
                          (irrep_occupation(i),i=1,number_of_irreps)
     read (filenumber_1e) number_of_abelian_irreps

!    count total number of electrons
     number_of_electrons = 0
     do i = 1, number_of_irreps
        number_of_electrons = number_of_electrons + irrep_occupation(i)
     end do

     allocate (multiplication_table(2*number_of_abelian_irreps,2*number_of_abelian_irreps))
     allocate (spinor(number_of_spinors))
     allocate (kramer_to_spinor(-number_of_spinors/2:number_of_spinors/2))

     read (filenumber_1e) ((multiplication_table(i,j),i=1,2*number_of_abelian_irreps),j=1,2*number_of_abelian_irreps)
     read (filenumber_1e) (spinor(i)%irrep,spinor(i)%abelian_irrep,spinor(i)%energy,i=1,number_of_spinors)

!    reorder spinors such that we have all occupied first
     index = 0
     do i = 1, number_of_irreps
        n = 0
        do j = 1, number_of_spinors
           if (spinor(j)%irrep == i) then
              n = n + 1
              if (n <= irrep_occupation(i)) then
                 index = index + 1
                 spinor(j)%index = index
                 spinor(j)%occupation = 1
              end if
           end if
        end do
     end do
     do i = 1, number_of_irreps
        n = 0
        do j = 1, number_of_spinors
           if (spinor(j)%irrep == i) then
              n = n + 1
              if (n > irrep_occupation(i)) then
                 index = index + 1
                 spinor(j)%index = index
                 spinor(j)%occupation = 0
              end if
           end if
        end do
     end do

     write (*,*) " Initialized reading from MRCONEE"
     write (*,*) " Core energy: ", core_energy
     write (*,*) " Breit interaction: ", breit_included

     open (filenumber_2e, file='MDCINT', Form='UNFORMATTED')
     read (filenumber_2e) date_of_generation,  &
                          time_of_generation,  &
                          number_of_kramerspairs, &
                          (kramer_to_spinor(i),kramer_to_spinor(-i),i=1,number_of_spinors/2) 

     if (2 * number_of_kramerspairs /= number_of_spinors ) then
        write (*,*) number_of_spinors, number_of_kramerspairs
        stop 'inconsistent MRCONEE and MDCINT files'
     end if 
     write (*,*) " Initialized reading from MDCINT"
     write (*,*) " MDCINT created: ",date_of_generation,"  ",time_of_generation
     
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine process_1e

     integer                :: i, j
     real(8), allocatable   :: integral(:,:,:)

     allocate (integral(number_of_spinors,number_of_spinors,2))
     read (filenumber_1e) ((integral(i,j,1),integral(i,j,2),i=1,number_of_spinors),j=1,number_of_spinors)

     do i = 1, number_of_spinors
        do j = 1, number_of_spinors
           write (filenumber_55,'(E28.20,4i4)') integral(i,j,1), &
                 spinor(i)%index,                              &
                 spinor(j)%index,                              &
                 0,                                            &
                 0                                              
        end do
     end do
     deallocate (integral)

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_2e
  
     integer                :: nonzero, ikr, jkr, inz
     integer, allocatable   :: indk(:), indl(:)
     real(8), allocatable   :: integral(:)

     allocate (integral(number_of_spinors**2))
     allocate (indk(number_of_spinors**2))
     allocate (indl(number_of_spinors**2))

     do
        read (filenumber_2e) ikr, jkr, nonzero, (indk(inz), indl(inz), inz=1, nonzero), (integral(inz), inz=1, nonzero)
        if (ikr == 0) exit
        do inz = 1, nonzero
           write (filenumber_55,'(E28.20,4i4)') integral(inz), &
                 spinor(kramer_to_spinor(ikr))%index,                 &
                 spinor(kramer_to_spinor(jkr))%index,                 &
                 spinor(kramer_to_spinor(indk(inz)))%index,           &
                 spinor(kramer_to_spinor(indl(inz)))%index             
!          make also kramers-related integral if desired
!          note that we assume real integrals, in which case the kr integral is identical
           if (generate_full_list) &
              write (filenumber_55,'(E28.20,4i4)') integral(inz), &
                    spinor(kramer_to_spinor(-ikr))%index,                 &
                    spinor(kramer_to_spinor(-jkr))%index,                 &
                    spinor(kramer_to_spinor(-indk(inz)))%index,           &
                    spinor(kramer_to_spinor(-indl(inz)))%index             
       end do
    end do
  
    deallocate(integral)
    deallocate(indk)
    deallocate(indl)
    
  
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function irrep_reordered (irrep)
     
     integer, intent(in) :: irrep

! Reorder the irreps such that boson irreps are given first and fermion irreps follow

     if (irrep > number_of_abelian_irreps) then
        irrep_reordered = irrep - number_of_abelian_irreps
     else
        irrep_reordered = irrep + number_of_abelian_irreps
     end if 
  
  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_fort55

  integer               :: i, j
  integer, allocatable  :: reordered_multiplication_table(:,:)

  allocate (reordered_multiplication_table(2*number_of_abelian_irreps,2*number_of_abelian_irreps))

! Reorder to list the boson irreps first
  do j = 1, 2 * number_of_abelian_irreps
     do i = 1, 2 * number_of_abelian_irreps
        reordered_multiplication_table(irrep_reordered(i), irrep_reordered(j)) = irrep_reordered(multiplication_table(i,j))
     end do
  end do
  spinor(:)%abelian_irrep = spinor(:)%abelian_irrep + number_of_abelian_irreps

  open (filenumber_55, file='fort.55', Form='FORMATTED')
  write (filenumber_55,'(2i6)') number_of_spinors, number_of_electrons
  write (filenumber_55,'(8i6)') (spinor(i)%abelian_irrep,i=1,number_of_spinors)
  write (filenumber_55,'(2i6)') -3
  write (filenumber_55,'(2i6)') 2 * number_of_abelian_irreps
  do j = 1, 2 * number_of_abelian_irreps
     write (filenumber_55,'(8i6)') (reordered_multiplication_table(i,j),i=1,2*number_of_abelian_irreps)
  end do 
  call process_2e
  call process_1e
  write (filenumber_55,'(E28.20,4i4)') core_energy,0,0,0,0

  deallocate (multiplication_table)
  deallocate (reordered_multiplication_table)

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module

  program make_mrcc_interface
  
  use dirac_to_mrcc

  call initialize
  call write_fort55

  end
