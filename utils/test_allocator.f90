!
! $Id$
!
!
! $Id$
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

         integer, parameter :: kchar = 1, kint = 8, kreal = 8, kcomplex = 8

         integer(kind=kint)     :: v2size
         integer                :: v1size, p_mem2,p_mem3
!         integer                :: v2size, v1size, p_mem2,p_mem3
         double precision, allocatable, target :: vdatadb(:), vdatarl(:)
         double precision, pointer             :: p_mem
         v2size = 10

!        write (*,*) 'memory address 0a: ',c_loc(vdatarl)
!        write (*,*) 'memory address 0b: ',c_loc(vdatadb)

         call allocator_init
!         call allocator_setGroupName('main')

         call alloc(vdatadb,v2size)
         call alloc(vdatarl,v2size)

         p_mem => vdatadb(0)

         call get_memory_address(vdatadb, p_mem2)
         write (*,*) 'memory address 2: [',p_mem2,']'
         call get_memory_address(vdatarl, p_mem3)
         write (*,*) 'memory address 3: [',p_mem3,']'

!        write (*,*) 'memory address 4: ',c_loc(vdatarl)
!        write (*,*) 'memory address 5: ',c_loc(vdatadb)

         call compare_memory(p_mem2,p_mem3)

         write (*,*) 'outside rbla'
         do i=1,v2size
            vdatarl(i) = 0.001 * rand()
            write (*,*) 'prep data ',vdatarl(i)
         enddo

         write (*,*) 'before rbla'
         call rbla(vdatarl,v2size)
         write (*,*) 'before dbla'
         call dbla(vdatarl,vdatadb,3)

         do i=1,v2size
            write (*,*) 'd(',i,')=',vdatadb(i),' r(',i,')=',vdatarl(i)
            write (*,*) 'kind(d) =',kind(vdatadb(i))
            write (*,*) 'kind(r) =',kind(vdatarl(i))
         enddo

         call dealloc(vdatadb)
         call dealloc(vdatarl)

         call get_memory_address(vdatadb, p_mem2)
         write (*,*) 'memory address 2: [',p_mem2,']'
         call get_memory_address(vdatarl, p_mem3)
         write (*,*) 'memory address 3: [',p_mem3,']'
         call compare_memory(p_mem2,p_mem3)
!        write (*,*) 'memory address 4: ',c_loc(vdatarl)
!        write (*,*) 'memory address 5: ',c_loc(vdatadb)
         write (*,*) 'is 4 allocated? ',allocated(vdatarl)
         write (*,*) 'is 5 allocated? ',allocated(vdatadb)

!         call allocator_PrintStats()

         v2size = 100000000
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size)
         write (*,*) 'is it allocated? ',allocated(vdatadb)
         call dealloc(vdatadb)
         call dealloc(vdatadb)
          
         v2size = 1000000000
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size)
         write (*,*) 'is it allocated? ',allocated(vdatadb)
         call dealloc(vdatadb)

         v2size = 2207784000
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size)
         write (*,*) 'is it allocated? ',allocated(vdatadb)
         call dealloc(vdatadb)

         v2size = 22077845000
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size)
         write (*,*) 'is it allocated? ',allocated(vdatadb)
         call dealloc(vdatadb)

         v2size = 100000000000
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size)
         write (*,*) 'is it allocated? ',allocated(vdatadb)
         call dealloc(vdatadb)


         v2size = 100000000000
         write (*,*) 'large allocation test, n =',v2size
         call alloc(vdatadb,v2size)
         write (*,*) 'is it allocated? ',allocated(vdatadb)
         call dealloc(vdatadb)


      end program alloc_test

      subroutine compare_memory(location_1, location_2)
           integer           :: location_1, location_2 

           write (*,*) 'a=',location_1,'  b=',location_2
           if (location_1 .eq. location_2) then
              write (*,*) 'addresses are the same'
           else
              write (*,*) 'addresses are different'
           endif

      end subroutine compare_memory

      subroutine rbla(data, size)
           character(len=1)  :: dummy
           integer           :: size
           real              :: data(*)

!           call allocator_setGroupName('rbla')
           write (*,*) 'inside rbla 1'
           do i=1,size
              write (*,*) ' data ',data(i)
           enddo
           write (*,*) 'inside rbla 2'
           do i=1,size
              data(i) = 2.1 * rand()
              write (*,*) ' new data ',data(i)
           enddo
           write (*,*) 'leaving rbla'
!           call get_traceback_info (printPrettyTrace,dummy,1)
!           call get_traceback_info (printAllocCaller,dummy,1)
      end subroutine rbla

      subroutine dbla(data1, data2, size)
           integer           :: size
           double precision  :: data2(*)
           real              :: data1(*)

           write (*,*) 'inside dbla'
           do i=1,size
              data2(i) = data1(i) / 1000
           enddo
           write (*,*) 'leaving dbla'
      end subroutine dbla

