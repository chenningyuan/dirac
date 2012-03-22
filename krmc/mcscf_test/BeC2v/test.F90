program test

  use unit_testing
  implicit none

  integer :: ier
  integer :: key
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
  real*8,  allocatable, dimension(:) :: fq    
  real*8,  allocatable, dimension(:) :: h2xy  
  real*8,  allocatable, dimension(:) :: h2yx  
  real*8,  allocatable, dimension(:) :: pv    

  open(60,file='raddfq90.inp',status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'test: ERROR while creating raddfq90.inp file'
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

  allocate( ipqtoq(4,8) , stat=ier ); call errall(ier,'test','ipqtoq')
  allocate( fq    (norbt*nasht*nz) , stat=ier ); call errall(ier,'test','fq    ')
  allocate( h2xy  (norbt*nasht*nz*3) , stat=ier ); call errall(ier,'test','h2xy  ')
  allocate( h2yx  (norbt*nasht*nz*3) , stat=ier ); call errall(ier,'test','h2yx  ')
  allocate( pv    (nasht*nasht*nnashx*nz*3) , stat=ier ); call errall(ier,'test','pv    ')

  read(60,*)
  read(60,*) ipqtoq
  read(60,*)
  read(60,*) h2xy  
  read(60,*)
  read(60,*) h2yx  
  read(60,*)
  read(60,*) pv    

  close(60,status='keep')

  key = 1
  call raddfq90( &
          fq     &
         ,h2xy   &
         ,h2yx   &
         ,pv     &
         ,nxy    &
         ,ix     &
         ,iy     &
         ,irepij &
         ,iprint &
         ,nz     &
         ,ipqtoq &
         ,lupri  &
         ,norbt  &
         ,nasht  &
         ,nnashx &
         ,key    &
         )

  call compare_files(60,'raddfq90.out','raddfq90.out_ref',1.0d-14)

  deallocate( ipqtoq )
  deallocate( fq     )
  deallocate( h2xy   )
  deallocate( h2yx   )
  deallocate( pv     )

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
