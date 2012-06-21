!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine get_dirac_ints(nort,eps,nn,ind,v,nint)
!------------------------------------------------------------------------------
!
!  provide data (integrals,orbital energies) for high-sect CCSD from DIRAC
!  
!  nort - number of orbitals  
!  
!  nint - number of 2-el.integrals
!  
!  called from readv(ccsd3.F) -  HighSectFSCC method
!
!------------------------------------------------------------------------------
#if defined MOD_HSFS || HSFS_DYN_ALLOC
use print_level
use DIRAC_HSFS

#if defined HSFS_DYN_ALLOC
use hsfs_parameters
#else
#include "param_ccsd.h"
#endif

#include "implicit.h"
#include "priunit.h"
   real*8, intent(inout) :: eps(*), v(*)
   integer, intent(inout) :: ind(*),nint, nort
   integer, dimension(8) :: irreps
   logical :: unread_file

   ! ... to be replaced by module 
   ipack(i1,i2,i3,i4)=ior(ior(ishft(i1,24),ishft(i2,16)),ior(ishft(i3,8) ,i4))

   call header('Output from "get_dirac_ints"',-1)

!  call set_print_level(iprint_hsfs)
!  write(LUPRI,'(2x,a,i5)') 'getiprint=',getiprint()
   iprint_hsfs = getiprint()
!  write(LUPRI,'(2x,a,i5)') 'iprint_hsfs=',iprint_hsfs

   call open_dirac_hsfs_file; rewind(luhsfs)
   write(LUPRI,*) luhsfs,' file open'

   read(luhsfs,iostat=ios) nirp,nort
   if (ios.ne.0) then
     write(LUPRI,*) 'broken reading, nirp=',nirp,' nort=',nort
     call quit('error in luhsfs reading 1!')
   endif

   write(LUPRI,'(2x,a,i2,a,i4)') 'read nirp=',nirp,' nort=',nort

   read(luhsfs,iostat=ios) (irreps(i),i=1,nirp)
   if (ios.ne.0) call quit('error in luhsfs reading 2!')
   write(LUPRI,'(2x,a,8i4)') 'read irreps:',(irreps(i),i=1,nirp)

   write(LUPRI,'(2x,a,i5)') 'going to read orbital energies, nort=',nort

   do i=1, nort
    read(luhsfs,iostat=ios) eps(i)
    if (ios.ne.0) call quit('error in luhsfs reading 3 !')

    if (iprint_hsfs >= 3) then
      write(LUPRI,*) i,eps(i)
    endif

   enddo 

   read(luhsfs,iostat=ios) nkr
   if (ios.ne.0) call quit('error in luhsfs nint reading!')

   write(LUPRI,*) 'nkr=',nkr

   write(LUPRI,'(2x,a)') 'Going to read integrals & indexes:'

   nintx=0
   unread_file=.true.
   do while (unread_file) 
    nintx = nintx + 1
    if (nintx > ndint) then
      write(LUPRI,*) 'nintx,ndint:',nintx,ndint
      call quit('get_dirac_ints: number of integrals exceeds allocated arrays !')
    endif

    read(luhsfs,iostat=ios) i1,i2,i3,i4,v(nintx)
    if (ios.ne.0) call quit('error in luhsfs i1,i2,i3,i4,v(i) reading !')
    ind(nintx)=ipack(i1,i2,i3,i4) ! possible source of overflow !
    
    unread_file=.not.( (i1.eq.0).and.(i2.eq.0)) 
    if (iprint_hsfs >= 8) then
      write(LUPRI,'(i5,a,4i4,1x,a,i11,a,f18.12)') nintx,'.',i1,i2,i3,i4,'pck in:',ind(nintx),' int:',v(nintx)
    endif

   enddo
    
   nint=nintx-1

   write(LUPRI,*) 'read all integrals&indexes from the file, their number,  nint=',nint

   call close_dirac_hsfs_file

!   call quit('get_dirac_ints')

#endif
end subroutine get_dirac_ints
