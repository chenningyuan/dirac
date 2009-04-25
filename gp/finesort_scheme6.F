      subroutine init_metadata_file(i_file)
      
         use memory_allocator
#include "implicit.h"
#include "bucket_sort.h"
#include "infpar.h"
         CHARACTER*10 NAME*5,FNODE,LFNAM

         IRECL  = IRECLEN(0,6,0)

         WRITE (NAME,'(A4,I1)') 'META',I_FILE
         FNODE = LFNAM(NAME)
         LMFIL = I_MD_UNIT+1

         OPEN (UNIT=LMFIL,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL,
     &         iostat=istat)

         if (istat.ne.0) call quit('Error opening metadata file!')

         irec_md_global_lastwritten = 0
c        call alloc(irec_md_lastwritten,N_BUFFER)
         do i=1,N_BUFFER
            irec_md_lastwritten(i) = 0
         enddo
      
      end subroutine init_metadata_file

c
c driver for fine-sorting in the 2HT step in scheme 6
c
       subroutine sort_buffers(thisnode,i_file)
       use memory_allocator
#include "implicit.h"
#include "bucket_sort.h"
#include "infpar.h"

c       integer, intent(in) :: i_file
       integer :: i_file, thisnode
       integer :: data_chunksize, buffer, idim3

       call print_resort_banner()

       data_chunksize  = NGBFSZ 

c aspg, hardcoded! change this later
       i_metadata_file = 1

c we find the maximum buffer size in terms of written records, and
c use that as scratch space
       max_records_per_buffer = -1

       do buffer = 1, N_BUFFER 
          if (thisnode.eq.node_for_buffer(buffer)) then
c            write (*,*) '   buffer ',buffer,' r_in_buf: ',
c     &                  records_in_buffer(buffer)
             if (max_records_per_buffer.lt.
     &           records_in_buffer(buffer)) then
                max_records_per_buffer = records_in_buffer(buffer) 
             endif 
          endif
       enddo

c aspg, debug
c      write (*,*) ' maximum records per buffer',max_records_per_buffer
c aspg, debug

       if (max_records_per_buffer.eq.0) then
          call quit('All zero-sized buffers is not possible!')
       endif
       

c more scratch allocation
       imdsz = (KLpairsInBuffer+1)
       isiz0 = max_records_per_buffer*data_chunksize
       isiz1 = isiz0*KLpairsInBuffer
       isiz2 = isiz1*2

c aspg, debug
       write (*,*) ' '
       write (*,*) ' size of data chunk (units) ',data_chunksize    
       write (*,*) ' metadata  storage unit     ',imdsz
       write (*,*) ' temporary storage unit 0, size: ',isiz0
       write (*,*) ' temporary storage unit 1, size: ',isiz1
       write (*,*) ' temporary storage unit 2, size: ',isiz2
       write (*,*) ' '
c aspg, debug


       call init_metadata_file(i_metadata_file)  ! one record per buffer, so 
                                                 ! record should be the dimension
                                                 ! of metadata+space to hold
                                                 ! records in the linked list?

c       write (*,*) '       ++ Sorting buffer '

       do buffer = 1, N_BUFFER 
c to implement
         if (thisnode.eq.node_for_buffer(buffer)) then
c        if buffer is in my node do work

 
c aspg, debug
c         write (*,*) 'calling read_and_resort'
c         write (*,*) '   buffer   ',buffer
cc        write (*,*) '   i_file   ',i_file
c         write (*,*) '   kl pairs ',KLpairsInBuffer
c aspg, debug

         call read_and_resort(buffer,
     &                        i_file, i_metadata_file,
     &                        data_chunksize,
     &                        KLpairsInBuffer,imdsz,isiz1)

         write (6,'(A31I8A9I5)') '     + Done fine-sorting buffer',
     &                       buffer,' at node ',thisnode

         call flush(6)
         endif

      enddo

      call print_resort_stats()

      call timer('sort_buffers',TIMSTR,TIMEND)

      end 

      subroutine print_resort_banner()
         write (*,*) ' '
         write (*,*) '+---> Fine-sorting 1HT integral file(s)'
         write (*,*) ' '
      end

      subroutine print_resort_stats()
         write (*,*) ' '
         write (*,*) '    + Fine-sort statistics coming soon'
         write (*,*) ' '
         write (*,*) ' '
      end
      
      subroutine read_and_resort(iactive_buffer,i_file,
     &                           i_metadata_file,
     &                           data_chunksize,
     &                           kl_in_buffer,
     &                           imdsz,isiz1)
      use memory_allocator
#include "implicit.h"
#include "bucket_sort.h"
#include "infpar.h"
       integer, intent(in) :: i_file, i_metadata_file
       integer :: iactive_buffer
       integer :: kl_in_buffer
       integer :: total_entries
       integer :: active_record
       integer :: mdata_idx, indexb
       REAL*8  :: RGBUF(NGBFSZ)
       DIMENSION IGBUF(2,NGBFSZ)
       integer, allocatable :: metadata(:,:), b(:,:), b2(:,:), count(:)
       real*8,  allocatable :: a(:), a2(:)
       integer, allocatable :: records_in_use(:,:)  ! keep track of which records are used in this buffer
c     common/totalreads/itotalreads

       iactive=iactive_buffer

       LGFIL    = I_BF_UNIT+I_FILE
       IREC     = LGREC(iactive_buffer,I_FILE)

       if (irec.eq.0) return ! if there is nothing in this buffer we leave

       LMFIL   = I_MD_UNIT+1
       imdrec  = irec_md_lastwritten(iactive_buffer)

       call alloc(records_in_use,3,max_records_per_buffer)
       call izero(records_in_use,3*max_records_per_buffer)

       call alloc(metadata,3,imdsz)

c aspg
c for using the "suck all data into memory" sorting strategy
c this could be an issue if the size of the buffer (number of kl pairs
c times the number of integrals with these kl pairs) becomes too
c big, because we may run out of memory... but the upside is that
c one avoids reading and writing several times...
       call alloc(a,isiz1)
       call izero(a,isiz1)

       call alloc(b,2,isiz1)
       call izero(b,2*isiz1)

c furthermore, for the key-indexed counting sorting algorithm we need the extra storage below
       call alloc(count,imdsz+1)
       count(1)=1
       forall (i=2:imdsz+1) count(i) = 0
c we don't initialize a2, b2 to zero because there will be no accumulation done in them; 
c rather, they will be assigned the values in a,b
       call alloc(a2,isiz1)
       call alloc(b2,2,isiz1)
c aspg

c
c using the "suck all data into memory" sorting strategy: the whole buffer goes 
c in memory and is subsequently sorted, and broken down into the original record
c layout on file.
c
      call put_buffer_in_memory(lgfil,irec,active_record,
     &                          records_in_use,igbuf,b,rgbuf,a,
     &                          total_entries,NGBFSZ,count,kl_in_buffer)
c aspg debug
c      write (*,*) 'buffer before sort'
c      call print_buffer(a,b,total_entries)
c aspg debug     

C selection sort is the simplest sorting applicable here, but
c scales quadratically wrt the size of "total_entries", so gets 
c inefficient very quickly
c
c       call selection_sort_in_memory(a,b,total_entries)
c
c key-indexed is linear in all operations, so hopefully it will show
c good performance in practice... but perhaps important here is to try
c to have moltra setup in such a way that the buffers are not toooo big,
c like if there are too many KL pairs for each buffer...
c
       call keyindexcount_sort_in_memory(a, b, total_entries,
     &                                   a2,b2,count,kl_in_buffer)
c aspg debug
c      write (*,*) 'buffer after sort'
c      call print_buffer(a,b,total_entries)
c aspg debug
c
c and now we want to break down A and B in such a way that they fit the
c layout of the read records, and we then overwrite those in the file
c
      indexb = 0
      mdata_idx = 1 

      call initialize_metadata(metadata,imdsz)
      
c     metadata(1,ii) = b(2,1)              ! first column in the buffer 
c     metadata(2,ii) = records_in_use(2,1)
c     metadata(3,ii) = records_in_use(3,1)
c     metadata(3,ii) = 0

      do j = 1, active_record
         n     = records_in_use(1,j)
         irec  = records_in_use(2,j)
         irecn = records_in_use(3,j)

c         write (*,*) ' layout:  n =',n,' irec =',irec,' irecn = ',irecn

         do k = 1, n            ! number of elements in the record
            indexb = indexb+1   ! indexb tracks total number of elements
            
            IGBUF(1,k) = b(1,indexb)
            IGBUF(2,k) = b(2,indexb)
            RGBUF(k)   = a(indexb)
c we change value of icol, so we have to mark which record we are still in
C and if there is more then one column in this record, we mark its beginning as well
            call set_metadata(irec,b,indexb,metadata)
         enddo

c         call print_metadata(metadata,mdata_idx)
         
c overwrite the unsorted buffer with the sorted one
         write(LGFIL,REC=irec) irecn,n,rgbuf,igbuf
         
         if ( irecn .eq. 0 ) then
            call save_metadata(lmfil,iactive,metadata,
     &                         irec_md_global_lastwritten,
     &                         irec_md_lastwritten, kl_in_buffer)

            deallocate(a)
            deallocate(b)
            deallocate(metadata)
            deallocate(records_in_use)

            return
         endif
         
      enddo
c aspg, strange bug! see comments above...
      iactive_buffer=iactive
      
      deallocate(a)
      deallocate(b)
      deallocate(metadata)
      deallocate(records_in_use)

      end subroutine read_and_resort

      subroutine set_metadata(irec,b,indexb,metadata)
         integer, intent(out) :: metadata(3,*)
         integer, intent(in)  :: irec, indexb, b(2,*)
         integer :: iprev

         iprev = indexb - 1
         
         md_idx_prev    = b(2,iprev) 
         md_idx_current = b(2,indexb)

         if (indexb.eq.1) then
            metadata(2,md_idx_current) = irec
         elseif ( indexb.gt.1.and.(b(2,iprev).ne.b(2,indexb))) then
            metadata(3,md_idx_prev)    = irec
            metadata(2,md_idx_current) = irec
         endif
         metadata(3,md_idx_current) = irec

      end subroutine set_metadata


      subroutine print_metadata(metadata,mdata_idx)
         integer, intent(in) :: mdata_idx, metadata(3,*)
         integer :: i

         write(6,*) 'metadata information stored so far'
         do i = 1, mdata_idx
            write(6,'(4I12)') i,(metadata(j,i),j=1,3)
         enddo

      end subroutine print_metadata


      subroutine initialize_metadata(metadata_array,metadata_size)
         integer, intent(in) :: metadata_size
         integer :: i, metadata_array(3,*)
         
         do i = 1, metadata_size
            metadata_array(1,i) = i
            metadata_array(2,i) = 0
            metadata_array(3,i) = 0
         enddo
         
      end subroutine initialize_metadata


      subroutine save_metadata(lmfil,iactive_buffer,metadata_array,
     &                         irec_md_global_lastwritten,
     &                         irec_md_lastwritten, KLpairsInBuffer)
         integer, intent(in) :: lmfil,iactive_buffer, KLpairsInBuffer
         integer :: l, imdrec, imdrecn
         integer :: metadata_array(3,*)
         integer :: irec_md_global_lastwritten, irec_md_lastwritten(*)

c         write(6,*) 'saving metadata information'
         
         do l = 1, KLpairsInBuffer

            if (l.eq.1) then
               imdrec  = 0
            else
               imdrec  = irec_md_global_lastwritten
            endif
               imdrecn = irec_md_global_lastwritten+1

            if ( metadata_array(1,l).ne.l ) then
               write(6,*) 'inconsistent kl? ',l,
     &                     metadata_array(1,l),
     &                     metadata_array(2,l),
     &                     metadata_array(3,l)

c               metadata_array(1,l) = l
            endif

            write(LMFIL,REC=imdrecn) imdrec,iactive_buffer,l,
     &                               metadata_array(1,l),
     &                               metadata_array(2,l),
     &                               metadata_array(3,l)

            irec_md_lastwritten(iactive_buffer) = imdrecn
            irec_md_global_lastwritten   = imdrecn

         enddo

      end subroutine save_metadata

      subroutine put_buffer_in_memory(lgfil,start_record,active_record,
     &                                records_in_use,igbuf,b,rgbuf,a,
     &                                total_entries,max_size,count,
     &                                KLpairsInBuffer)
         integer, intent(in) :: max_size, lgfil
         integer :: records_in_use(3,*), count(KLpairsInBuffer)
         integer, intent(out) :: total_entries
         integer :: igbuf(2,max_size), b(2,*)
         real*8  :: rgbuf(max_size), a(*)
         integer :: irec, start_record, active_record

         total_entries = 0
         active_record = 0

         irec = start_record
C      first we want to have all the data from the buffer in memory...
    1    if (irec.eq.0) then
            goto 2
         endif

           READ (LGFIL,REC=IREC) IRECN,N,RGBUF,IGBUF
           active_record = active_record+1

C     we first keep track of how much data each read has, becase later we will write back
c     to the same locations..

           records_in_use(1,active_record) = N
           records_in_use(2,active_record) = IREC
           records_in_use(3,active_record) = IRECN

           IREC = IRECN
C
C     first we put everyone for this buffer in memory 
C
           DO I = 1, N
              a(i+total_entries)   = RGBUF(I)
              b(1,i+total_entries) = IGBUF(1,I)
              b(2,i+total_entries) = IGBUF(2,I)    ! this is the ICOL equivalent
              count(IGBUF(2,I) + 1) = count(IGBUF(2,I) + 1) + 1  ! for key-indexed counting sort
           ENDDO
           total_entries = total_entries + N
C
C     We read in reverse order, if we have the first block we're done
C
         goto 1

C     so here we read everything, and kept the structure of the buffer in terms
c     of the records written... now what we want to reorder the contents of temp
c     variables A and B

    2 continue

c for key-indexed counting sort, this will set the final offset
c      do j = 2, KLpairsInBuffer+1
c         write (6,*) ' ---> # elements with IGBUF(2) =',(j-1),
c     &               ' is ',count(j)
c      enddo

      do j = 2, KLpairsInBuffer
         count(j) = count(j) + count(j-1) 
      enddo

c aspg, debug
      if (.false.) then
      write (*,*) ' ---> records in use for this buffer '
      do iriu = 1, active_record
         write (*,*) '      ',iriu,
     &                  '  ',records_in_use(1,iriu),
     &                  '  ',records_in_use(2,iriu),
     &                  '  ',records_in_use(3,iriu) 
      enddo

      
      do i = 1, KLpairsInBuffer
         write (6,*) ' ---> offsets for IGBUF(2) =',i,' is ',count(i)
      enddo
c aspg, debug
      call flush(6)
      endif

      end subroutine put_buffer_in_memory 





C soubroutines to sort buffer in-memory

       subroutine keyindexcount_sort_in_memory(a,b,total_entries,
     &               a2,b2,count,KLpairsInBuffer)
          integer :: total_entries, i, j, itmp, min, KLpairsInBuffer
          integer :: b(2,total_entries)
          integer :: b2(2,total_entries), count(KLpairsInBuffer)
          real*8  :: a(total_entries),   a2(total_entries)

c below we do the final two steps of key-indexed counting
          do i = 1, total_entries
             a2(   count(b(2,i) ))   = a(i)
             b2(1, count(b(2,i) ))   = b(1,i)
             b2(2, count(b(2,i) ))   = b(2,i)
             count( b(2,i) ) = count( b(2,i) ) + 1
          enddo

          do i = 1, total_entries
             a(i)   = a2(i) 
             b(1,i) = b2(1,i)
             b(2,i) = b2(2,i)
          enddo


       end subroutine keyindexcount_sort_in_memory


       subroutine selection_sort_in_memory(a,b,total_entries)
          integer :: total_entries, i, j, itmp, min
          integer :: b(2,total_entries)
          real*8  :: rtmp
          real*8  :: a(total_entries)

c below is a selection sort implementation...
          do i = 1, total_entries-1
             min = i
             do j = i+1, total_entries
                if ( b(2,j) .lt. b(2,min) ) min = j
             enddo   

             itmp      = b(2,i)
             b(2,i)    = b(2,min)
             b(2,min) = itmp

             itmp      = b(1,i)
             b(1,i)    = b(1,min)
             b(1,min) = itmp

             rtmp    = a(i)
             a(i)    = a(min)
             a(min) = rtmp

          enddo
          
       end subroutine selection_sort_in_memory
      
      subroutine print_buffer(a,b,total_entries)
#include "priunit.h"
         integer :: total_entries, i
         integer :: b(2,total_entries)
         real*8  :: a(total_entries)

         do i = 1, total_entries
            write (LUPRI,'(2X,3I1,F12.8)') i,b(1,i),b(2,i),a(i)
         enddo
      
      end subroutine print_buffer



      subroutine init_finesort_commons(kl_in_buf,NumBuffers,BuffNodes)
#include "implicit.h"
#include "bucket_sort.h"

         integer :: kl_in_buf, NumBuffers, BuffNodes(NumBuffers)
         integer :: ibuf

         KLpairsInBuffer = KL_IN_BUF

         do ibuf=1,N_BUFFER
            records_in_buffer(ibuf) = 0
            node_for_buffer(ibuf)   = BuffNodes(ibuf)
         enddo

      end subroutine init_finesort_commons


