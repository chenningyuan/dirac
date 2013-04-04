module pcm_write
   
   implicit none
      
   character(*), parameter, public :: pcm_file_name = 'PCM_mep_asc'
   integer,      parameter, public :: pcm_file_unit = 800 
   logical                , public :: pcm_file_exists = .false.
   logical                , public :: pcm_file_open = .false.

   public report_after_pcm_input
   public pcm_write_file

   private
  
   contains
   
   subroutine report_after_pcm_input()
      use pcmmod_cfg
!
! A subroutine to print the relevant settings of both DIRAC and PCMSolver.
!

      write(*, *) ' ===== Polarizable Continuum Model calculation set-up ====='
      write(*, *) '* Polarizable Continuum Model using PCMSolver external module:'
      write(*, *) '  1: Converged potentials and charges at tesserae representative points written on file.'
      
      if (pcmmod_doss) then
         write(*, *) '  2: Calculate the SS block of the electrostatic potential matrix.'
      else
         write(*, *) '  2: Skip the SS block of the electrostatic potential matrix.', &
        ' Use a Simple Coulombic Correction.'
      end if
      
      if (pcmmod_separate) then
         write(*, *) '  3: Separate potentials and apparent charges in nuclear and electronic.'
      else
         write(*, *) '  3: Use total potentials and apparent charges.'
      end if
      
      if (pcmmod_dospf) then
         write(*, *) '  4: Perform elimination of spin-orbit dependent terms in the PCM operator.'
      end if
      
      if (pcmmod_print > 5 .and. pcmmod_print < 10) then
         write(*, *) '  5: Print potentials at tesserae representative points.'
      else if (pcmmod_print > 10) then
         write(*, *) '  5: Print potentials and charges at tesserae representative points.'
      else
         write(*, *) '  5: Do not print potentials and charges.'
      end if     
    

! Should print info for the cavity and the solver here.
      
      write(*, *) '* References: '
      write(*, *) '  - J. Tomasi, B. Mennucci and R. Cammi:',                &
      '   "Quantum Mechanical Continuum Solvation Models", Chem. Rev., 105 (2005) 2999',  & 
      '  - PCMSolver, an API for the Polarizable Continuum Model electrostatic problem',   &    
      '    L. Frediani et al. '

   end subroutine

   subroutine pcm_write_file(nr_points, potentials, charges)

! Passed variables
      integer, intent(in) :: nr_points
      real(8), intent(in) :: potentials(nr_points)
      real(8), intent(in) :: charges(nr_points)
 
! Local variables
      real(8) :: tot_chg
      integer :: ipoint

      inquire(file = pcm_file_name, exist = pcm_file_exists)
      if (pcm_file_exists) then
         open(pcm_file_unit, &
              file = pcm_file_name, &
              status = 'old', &
              form = 'formatted', &
              access = 'sequential')
         close(pcm_file_unit, status = 'delete')
      end if
      open(pcm_file_unit, &
           file = pcm_file_name, &
           status = 'new', &
           form = 'formatted', &
           access = 'sequential')
      rewind(pcm_file_unit)

      tot_chg = 0.0d0
      do ipoint = 1, nr_points
         tot_chg = tot_chg + charges(ipoint)
         write(pcm_file_unit, *) '@point', ipoint, "MEP =", potentials(ipoint), "ASC =", charges(ipoint)
      end do
      write(pcm_file_unit, *) 'Number of tesserae', nr_points, 'Sum of apparent charges', tot_chg
    
   end subroutine 
   
end module
