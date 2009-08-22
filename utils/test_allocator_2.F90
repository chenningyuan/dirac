!
!      Copyright (c) 2008 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.sdu.dk
!
!
! $Id: test_allocator.F90 7732 2008-09-16 20:30:19Z bast $
!
!
! $Id: test_allocator.F90 7732 2008-09-16 20:30:19Z bast $
!
! memory allocation subsystem
!
! initially conceived to be used in the dirac and dalton packages
! as a replacement to the memget/memrel routines and related
!
! design/discussions    : ulf ekstrom, andre gomes, jetze sikkema
! inital implementation : andre gomes
!                         vu amsterdam, winter/spring 2008
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test program for memory allocator for dirac/dalton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program alloc_test

         use memory_allocator

         integer, parameter :: kreal = 8, kcomplex = 8
         integer, parameter :: charArrayLength = 132

         integer                :: v2size, v1size, p_mem2,p_mem3,i,f1,f2
         integer                :: maxsize, minimal = 20, verbose = 0
         double precision, allocatable, target :: vdatadb(:)
         real(kind=kreal), allocatable :: matrix(:,:), test(:), vdatarl(:)
         real(kind=kreal), pointer             :: p_mem, p_mem42(:)
       
         character(len=charArrayLength) :: kbytesAvailable

         f1= -10
         f2=  10
         maxsize   = 1500
         minimal   = 16
         nstep     = 10
         direction = 1

         v2size = 10

!        write (*,*) 'memory address 0a: ',c_loc(vdatarl)
!        write (*,*) 'memory address 0b: ',c_loc(vdatadb)

         call allocator_init
!         call allocator_setGroupName('main')

         call alloc(test,(1024*1024/kind(test)), id="test")
!        do i=1, 11 
!        test(i) = 1.2D0
!        write (*,*) 'test(',i,'):',test(i)
!        test(i) = 0.0
!        enddo
         call dealloc(test)

         write (*,*) ' before call_allocation'
         call test_trace
         write (*,*) ' after  call_allocation'
         write (*,*) ' '
         call alloc(test,0, d1=f1, d2=f2, id="test")
!        do i=f1, f2 
!!       test(i) = 1.22D0
!        write (*,*) 'test(',i,'):',test(i)
!        enddo
         call dealloc(test)

         call allocator_report(0)
  
      end program alloc_test

      subroutine test_trace
         call call_allocatio1
      end subroutine test_trace

      subroutine call_allocatio1
         write (*,*) 'level 1'
         call call_trace_and_allo2
      end subroutine call_allocatio1

      subroutine call_trace_and_allo2
         use memory_allocator

         integer :: charArrayLength,depth
         character(len=132) :: trace_result
         real(kind=8), allocatable :: bla(:)

         charArrayLength=132

         write (trace_result,'(A)') 'main'
  
         call alloc(bla,(1024*1024/kind(bla)))

         call get_traceback_info (4,trace_result,1,depth)

         write (*,*) 'depth is ',depth
         if (trace_result(1:4).eq.'main') then
            write (*,*) 'here we go: depth is ',depth
         else
            write (*,'(A)') trace_result(1:5)
         endif
         call dealloc(bla)

      end subroutine call_trace_and_allo2

