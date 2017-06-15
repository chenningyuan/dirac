#if defined (VAR_MPI)

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE mpi_slave_server()

      use interface_to_mpi

      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
! This is a server module that executes all slave node tasks required in
! the parallel POLPRP run. The services are controlled by a corresponding
! broadcast code. The codes are contained in the polprp_servercodes.h
! file. A termination signal triggers exit from the server loop.
!
!--------------- Includes--------------------------------------
!
#include  "polprp_servercodes.h"
#include  "../relccsd/symm.inc"
#include  "../relccsd/ccpar.inc"
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      Integer                           :: taskid = 0   ! defines the server's task
      REAL*8,allocatable,dimension(:)   :: buf
      logical                           :: done
      INTEGER                           :: irp,istart,nint,mint
      INTEGER                           :: ialloc,locoff,istat
      INTEGER                           :: efflen
!
!---------------Interface area --------------------------------------
!
      interface

        SUBROUTINE mpi_vovv_complete(ra1,i8)
          REAL*8, dimension(:)             :: ra1
          INTEGER*8                        :: i8
        END SUBROUTINE

        SUBROUTINE mpi_vvvv_complete(ra1,i8)
          REAL*8, dimension(:)             :: ra1
          INTEGER*8                        :: i8
        END SUBROUTINE

      end interface
!
!---------------Executable code--------------------------------------
!

!----------------------------------------------------------------
      DO       ! ---------- INFINITE SERVER MAIN LOOP -----------
!----------------------------------------------------------------

        call interface_mpi_BCAST(taskid,1,MASTER,global_communicator)
        write(*,*) 'Node',MYPROC,' has received taskid',taskid

        select case(taskid)


!-----
!-----  TASK: say hello
!-----


          case(SERVER_HELLO)
          write(*,*) 'Node',MYPROC,' confirms existence.'


!-----
!-----  TASK: terminate server loop
!-----


          case(SERVER_DONE)
          write(*,*) 'Node',MYPROC,' terminates service.'
          exit   ! terminates infinite do loop


!-----
!-----  TASK: provide local share of VVVV integrals in irep IRP
!-----


          case(SERVER_VVVV)
          call interface_mpi_recv(irp,1,MASTER,MSGN,global_communicator)
          write(*,*) 'Node',MYPROC,' obtained irp:',irp
          if(irp.lt.1.or.irp.gt.nrep) then
            write(*,*)  '**** Warning on node',MYPROC,':'
            write(*,*)  '**** irep number out of bounds,'
            write(*,*)  '**** returning with zero batch length.'
            efflen = 0
            return
          endif
          efflen = NVVT(irp) * NVVT(irp) !efflen is max number of communicated integrals
                                         !in this irep
          if(efflen.gt.0) then
            allocate(buf(efflen*RCW),stat = ialloc)
            if(ialloc.ne.0)
     &          call quit('Allocation problem in mpi_slave_server!')
            buf = 0.0d0

            locoff = 0
            efflen = 0  !efflen will be recalculated depending on batch length
            istart = 0

            MINT = NVVT(IRP)  !maximum number of possible batches
            CALL GETVVVV (irp,istart,nint,done,buf,mint)

            write(*,*) 'Recovered ',nint,' batches on',MYPROC
!
! NINT contains the actual number of batches on *this* node.
!
            if (.not.done) then
               write(*,*) 'irrep',irp,' unfinished on node',myproc
               call quit('Error reading VVVV ints in mpi_slave_server.')
            endif

            efflen = NVVT(irp) * NINT
          endif

          call interface_mpi_send(efflen*rcw,1,MASTER,MSGN,
     &                            global_communicator)
          locoff = IDIST(3,1,irp)
          call interface_mpi_send(locoff*rcw,1,MASTER,MSGN,
     &                            global_communicator)

          if(efflen.gt.0) then
            call interface_mpi_send(buf,efflen*rcw,MASTER,MSGN,
     &                            global_communicator)
            deallocate(buf)
          endif


!-----
!-----  TASK: provide local share of VOVV integrals in irep IRP
!-----


          case(SERVER_VOVV)
          call interface_mpi_recv(irp,1,MASTER,MSGN,global_communicator)
          write(*,*) 'Node',MYPROC,' obtained irp:',irp
          if(irp.lt.1.or.irp.gt.nrep) then
            write(*,*)  '**** Warning on node',MYPROC,':'
            write(*,*)  '**** irep number out of bounds,'
            write(*,*)  '**** returning with zero batch length.'
            efflen = 0
            return
          endif
          efflen = NVO(irp) * NVVT(irp)  !efflen is max number of communicated integrals
                                       !in this irep
          if(efflen.gt.0) then
            allocate(buf(efflen*RCW),stat = ialloc)
            if(ialloc.ne.0)
     &        call quit('Allocation problem in mpi_slave_server!')
            buf = 0.0d0

            locoff = 0
            efflen = 0  !efflen will be recalculated depending on batch length
            istart = 0

            MINT = NVVT(IRP)  !maximum number of possible batches
            CALL GETVOVV (irp,istart,nint,done,buf,mint)

            write(*,*) 'Recovered ',nint,' batches on',MYPROC
!
! NINT contains the actual number of batches on *this* node.
!
            if (.not.done) then
              write(*,*) 'irrep',irp,' unfinished on node',myproc
              call quit('Error reading VOVV ints in mpi_slave_server')
            endif

            efflen = NVO(irp) * NINT
          endif
!
          call interface_mpi_send(efflen*rcw,1,MASTER,MSGN,
     &                            global_communicator)
          locoff = IDIST(3,2,irp)
          call interface_mpi_send(locoff*rcw,1,MASTER,MSGN,
     &                            global_communicator)

          if(efflen.gt.0) then
            call interface_mpi_send(buf,efflen*rcw,MASTER,MSGN,
     &                              global_communicator)
            deallocate(buf)
          endif


        end select

!----------------------------------------------------------------
      ENDDO    ! ----------- END OF SERVER MAIN LOOP ------------
!----------------------------------------------------------------


      return
      END SUBROUTINE  ! mpi_slave_server
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#endif
