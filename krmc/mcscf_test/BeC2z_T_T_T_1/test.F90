program test

  use unit_testing
  use mcscf_routines
  implicit none

  integer :: ier
  logical :: psymm 
  logical :: ksymm 
  logical :: tsymm 
  logical :: spinfr
  logical :: levyle
  integer :: nasht 
  integer :: nnashx
  integer :: nz    
  integer :: nzinci
  integer :: nm5   
  integer :: nbsym 
  integer :: nfsym 
  integer :: iprint
  integer :: lupri 
  integer, allocatable, dimension(:) :: nash  
  integer, allocatable, dimension(:,:) :: ipqtoq
  real*8,  allocatable, dimension(:,:,:,:,:) :: vmuuuu
  real*8,  allocatable, dimension(:,:,:,:,:) :: vduuuu

  open(60,file='m2dnz390.inp_T_T_T_1',status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'test: ERROR while creating m2dnz390.inp_T_T_T_1 file'
    stop
  endif

  read(60,*)
  read(60,*) psymm 
  read(60,*)
  read(60,*) ksymm 
  read(60,*)
  read(60,*) tsymm 
  read(60,*)
  read(60,*) spinfr
  read(60,*)
  read(60,*) levyle
  read(60,*)
  read(60,*) nasht 
  read(60,*)
  read(60,*) nnashx
  read(60,*)
  read(60,*) nz    
  read(60,*)
  read(60,*) nzinci
  read(60,*)
  read(60,*) nm5   
  read(60,*)
  read(60,*) nbsym 
  read(60,*)
  read(60,*) nfsym 
  read(60,*)
  read(60,*) iprint
  read(60,*)
  read(60,*) lupri 

  allocate( nash  (2) , stat=ier ); call errall(ier,'test','nash  ')
  read(60,*)
  read(60,*) nash  

  allocate( ipqtoq(4,8) , stat=ier ); call errall(ier,'test','ipqtoq')
  read(60,*)
  read(60,*) ipqtoq

  allocate( vmuuuu(8,8,8,8,2) , stat=ier ); call errall(ier,'test','vmuuuu')
  read(60,*)
  read(60,*) vmuuuu

  allocate( vduuuu(4,4,10,2,3) , stat=ier ); call errall(ier,'test','vduuuu')

  close(60,status='keep')

  call m2dnz390( &
          psymm  = psymm  &
         ,ksymm  = ksymm  &
         ,tsymm  = tsymm  &
         ,spinfr = spinfr &
         ,levyle = levyle &
         ,nasht  = nasht  &
         ,nnashx = nnashx &
         ,nz     = nz     &
         ,nzinci = nzinci &
         ,nm5    = nm5    &
         ,nbsym  = nbsym  &
         ,nfsym  = nfsym  &
         ,iprint = iprint &
         ,lupri  = lupri  &
         ,nash   = nash   &
         ,ipqtoq = ipqtoq &
         ,vmuuuu = vmuuuu &
         ,vduuuu = vduuuu &
         )

  call compare_vector(60,'m2dnz390.out_T_T_T_1',vduuuu,1.0d-15)

  deallocate( nash   )
  deallocate( ipqtoq )
  deallocate( vmuuuu )
  deallocate( vduuuu )

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
