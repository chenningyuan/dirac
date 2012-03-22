program test

  use unit_testing
  implicit none

  integer   :: ier, key_unit_test
  character :: chtmp*6

  integer :: nxy, ix, iy, irepij, iprint
  integer :: nz, ipqtoq(4,0:7), lupri, norbt, nasht, nnashx
  real*8, allocatable :: h2xy(:,:,:,:), h2yx(:,:,:,:)
  real*8, allocatable :: pv(:,:,:,:,:)
  real*8, allocatable :: fq(:,:,:)

  ! open file
  open(60,file='raddfq90.uinp',status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'test_raddfq90: ERROR while creating raddfq90.uinp file'
    stop
  endif

  rewind 60

  ! write input data to the file
  read(60,*) chtmp,nxy
  read(60,*) chtmp,ix
  read(60,*) chtmp,iy
  read(60,*) chtmp,irepij
  read(60,*) chtmp,iprint
  read(60,*) chtmp,nz
  read(60,*) chtmp,lupri
  read(60,*) chtmp,norbt
  read(60,*) chtmp,nasht
  read(60,*) chtmp,nnashx

  allocate( h2xy(norbt,nasht,nz,3) , stat=ier ); call errall(ier,'test','h2xy')
  allocate( h2yx(norbt,nasht,nz,3) , stat=ier ); call errall(ier,'test','h2yx')
  allocate( pv  (nasht,nasht,nnashx,nz,3) , stat=ier ); call errall(ier,'test','pv')
  allocate( fq  (norbt,nasht,nz) , stat=ier ); call errall(ier,'test','fq')
  fq = 0.0d0

  read(60,*)
  read(60,*) ipqtoq
  read(60,*)
  read(60,*) h2xy
  read(60,*)
  read(60,*) h2yx
  read(60,*)
  read(60,*) pv

  ! close file
  close(60,status='keep')

  key_unit_test = 1
  call raddfq90(fq,h2xy,h2yx,pv,nxy,ix,iy,irepij,iprint, &
                nz,ipqtoq,lupri,norbt,nasht,nnashx,key_unit_test)

  deallocate( h2xy )
  deallocate( h2yx )
  deallocate( pv   )
  deallocate( fq   )

  call compare_files(60,'raddfq90.uout','raddfq90.uref',1.0d-14)

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

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
