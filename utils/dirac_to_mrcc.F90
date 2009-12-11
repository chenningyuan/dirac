module dirac_to_mrcc
  implicit none

  integer                :: number_of_spinors, number_of_electrons, number_of_irreps, number_of_abelian_irreps
  integer, parameter     :: filenumber_1e = 21
  integer, parameter     :: filenumber_2e = 22
  integer, parameter     :: filenumber_55 = 55
  real(8)                ::  core_energy

  type SpinorInformation

      integer     :: irrep
      integer     :: abelian_irrep
      integer     :: number
      real(8)     :: energy

   end type SpinorInformation

  public initialize, write_fort55

  private read_general_info, process_1e, process_2e, spinor_reordered, irrep_reordered

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize ()

     integer       ::  number_of_kramerspairs
     character(10) ::  date_of_generation
     character(8)  ::  time_of_generation
     logical       ::  breit_included

     open (filenumber_1e, file='MRCONEE', Form='UNFORMATTED')
     read (filenumber_1e) number_of_spinors, &
                          breit_included,    &
                          core_energy
     read (filenumber_1e) number_of_irreps
     read (filenumber_1e) number_of_abelian_irreps
     write (*,*) " Initialized reading from MRCONEE"
     write (*,*) " Core energy: ", core_energy
     write (*,*) " Breit interaction: ", breit_included

     open (filenumber_2e, file='MDCINT', Form='UNFORMATTED')
     read (filenumber_2e) date_of_generation,  &
                          time_of_generation,  &
                          number_of_kramerspairs

     if (2 * number_of_kramerspairs /= number_of_spinors ) then
        write (*,*) number_of_spinors, number_of_kramerspairs
        stop 'inconsistent MRCONEE and MDCINT files'
     end if 
     write (*,*) " Initialized reading from MDCINT"
     write (*,*) " MDCINT created: ",date_of_generation,"  ",time_of_generation
     
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_general_info (multiplication_table,spinor)
     
     integer                             :: i, j
     integer, intent(out)                :: multiplication_table(2*number_of_abelian_irreps, 2*number_of_abelian_irreps)
     type(SpinorInformation)             :: spinor(number_of_spinors)

     read (filenumber_1e) ((multiplication_table(i,j),i=1,2*number_of_abelian_irreps),j=1,2*number_of_abelian_irreps)
     read (filenumber_1e) (spinor(i)%irrep,spinor(i)%abelian_irrep,spinor(i)%energy,i=1,number_of_spinors)


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
                 i,                                            &
                 j,                                            &
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
                 ikr,                                          &
                 jkr,                                          &
                 indk(inz),                                    &
                 indl(inz)                                      
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

  integer function spinor_reordered (spinor)

     integer, intent(in) :: spinor

! Reorder the spinors such that all occupied spinors are given first

     spinor_reordered = spinor
  
  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_fort55

  integer               :: i, j
  integer, allocatable  :: multiplication_table(:,:), reordered_multiplication_table(:,:)
  type(SpinorInformation), allocatable :: spinor(:)

  allocate (spinor(number_of_spinors))
  allocate (multiplication_table(2*number_of_irreps,2*number_of_irreps))
  allocate (reordered_multiplication_table(2*number_of_irreps,2*number_of_irreps))

  call read_general_info (multiplication_table,spinor)
  do j = 1, 2 * number_of_abelian_irreps
     do i = 1, 2 * number_of_abelian_irreps
        reordered_multiplication_table(irrep_reordered(i), irrep_reordered(j)) = irrep_reordered(multiplication_table(i,j))
     end do
  end do

  open (filenumber_55, file='fort.55', Form='FORMATTED')
  write (filenumber_55,'(2i6)') number_of_spinors, number_of_electrons
  write (filenumber_55,'(8i6)') (spinor(i)%abelian_irrep,i=1,number_of_spinors)
  write (filenumber_55,'(2i6)') -3
  write (filenumber_55,'(2i6)') number_of_abelian_irreps
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
