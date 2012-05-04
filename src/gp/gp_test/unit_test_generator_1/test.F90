program test

  use unit_testing
  use unit_test_generator
  implicit none

  integer :: ier
  integer :: nxy   
  integer :: ix    
  integer :: iy    
  integer :: irepij
  integer :: iprint
  integer :: nz    
  integer :: lupri 
  integer :: norbt 
  integer :: nasht 
  integer :: nnashx
  integer, allocatable, dimension(:,:) :: ipqtoq
  real*8,  allocatable, dimension(:,:,:,:) :: h2xy  
  real*8,  allocatable, dimension(:,:,:,:) :: h2yx  
  real*8,  allocatable, dimension(:,:,:,:,:) :: pv    
  real*8,  allocatable, dimension(:,:,:) :: fq    

  allocate( ipqtoq(4,8) , stat=ier ); call errall(ier,'test','ipqtoq')
  allocate( h2xy  (37,4,2,3) , stat=ier ); call errall(ier,'test','h2xy  ')
  allocate( h2yx  (37,4,2,3) , stat=ier ); call errall(ier,'test','h2yx  ')
  allocate( pv    (4,4,10,2,3) , stat=ier ); call errall(ier,'test','pv    ')
  allocate( fq    (37,4,2) , stat=ier ); call errall(ier,'test','fq    ')

  open(60,file='INPUT',status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'test: ERROR while creating INPUT file'
    stop
  endif

  read(60,*)
  read(60,*) nxy   
  read(60,*)
  read(60,*) ix    
  read(60,*)
  read(60,*) iy    
  read(60,*)
  read(60,*) irepij
  read(60,*)
  read(60,*) iprint
  read(60,*)
  read(60,*) nz    
  read(60,*)
  read(60,*) lupri 
  read(60,*)
  read(60,*) norbt 
  read(60,*)
  read(60,*) nasht 
  read(60,*)
  read(60,*) nnashx

  read(60,*)
  read(60,*) ipqtoq

  read(60,*)
  read(60,*) h2xy  

  read(60,*)
  read(60,*) h2yx  

  read(60,*)
  read(60,*) pv    

  read(60,*)
  read(60,*) fq    


  close(60,status='keep')

  !------------------------------------------------------------------------------!
  !             Generate test with input data writen to the file                 !
  !------------------------------------------------------------------------------!
  ! if raddfq90.inp_4_3_1 already exist delete it
  open(60,file='raddfq90.inp_4_3_1',status='unknown',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'test: ERROR while creating raddfq90.inp_4_3_1 file'
    stop
  endif
  close(60,status='delete')

  call set_suffix(nxy)
  call set_suffix(ix)
  call set_suffix(iy)
  call add_input_file ('raddfq90.inp')
  call set_driver_routine('raddfq90')
  call add_module('mcscf_routines')

  call add_opt_variable('key   ','1','integer')

  call add_inp_variable(nxy   ,'nxy   ')
  call add_inp_variable(ix    ,'ix    ')
  call add_inp_variable(iy    ,'iy    ')
  call add_inp_variable(irepij,'irepij')
  call add_inp_variable(iprint,'iprint')
  call add_inp_variable(nz    ,'nz    ')
  call add_inp_variable(lupri ,'lupri ')
  call add_inp_variable(norbt ,'norbt ')
  call add_inp_variable(nasht ,'nasht ')
  call add_inp_variable(nnashx,'nnashx')
  call add_inp_variable(ipqtoq,'ipqtoq')
  call add_inp_variable(h2xy  ,'h2xy  ')
  call add_inp_variable(h2yx  ,'h2yx  ')
  call add_inp_variable(pv    ,'pv    ')

  call add_inp_variable(fq    ,'fq    ')
  call add_out_variable(fq    ,'fq    ',compare='file-file',file='raddfq90.out')

  call generate_test()

  ! compare input file and test
  call compare_text_files(60, 'raddfq90.inp_4_3_1', 'raddfq90.inp_ref')
  call compare_text_files(60, 'test.F90_4_3_1', 'test.F90_ref')

  deallocate( ipqtoq )
  deallocate( h2xy   )
  deallocate( h2yx   )
  deallocate( pv     )
  deallocate( fq     )

  write(6,*) ' test result: =====  OK  ====='

contains

  subroutine errall(ier,chsub,charr)
  implicit none

  integer      :: ier
  character(*) :: chsub, charr

  if(ier /= 0)then
    write(6,*)
    write(6,'(a,a,a,a)')' Error in allocation of ',charr,' in subroutine ',chsub
    write(6,'(a,i4)'   )' iostat =',ier
    stop 'ERROR in allocation (see output form more information)'
  endif

  end subroutine

end program
