#if defined (VAR_MPI)

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE mpi_slave_server(iobase,intbuf)

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
      Integer                           :: iobase,intbuf

      Integer                           :: taskid = 0   ! defines the server's task
      REAL*8,allocatable,dimension(:)   :: buf
      logical                           :: done
      INTEGER                           :: irp,istart,nint,mint
      INTEGER                           :: ialloc,locoff,istat
      INTEGER                           :: efflen
   
      Integer                               :: i,ixx,jdummy
      REAL*8,allocatable,dimension(:)       :: bufr,bufi
      Complex*16,allocatable,dimension(:)   :: cbuf
      Integer,allocatable,dimension(:)      :: ioi,ioj
      Integer                               :: nbufscur,nact
      Character*10                          :: locfn
      Character*6                           :: nmbase='PP_ELE'
     
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


!----------------------------------------------------------------
!-----  TASK: say hello
!----------------------------------------------------------------


          case(SERVER_HELLO)
          write(*,*) 'Node',MYPROC,' confirms existence.'


!----------------------------------------------------------------
!-----  TASK: terminate server loop
!----------------------------------------------------------------


          case(SERVER_DONE)
          write(*,*) 'Node',MYPROC,' terminates service.'
          exit   ! terminates infinite do loop


!----------------------------------------------------------------
!-----  TASK: provide local share of VVVV integrals in irep IRP
!----------------------------------------------------------------


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


!----------------------------------------------------------------
!-----  TASK: provide local share of VOVV integrals in irep IRP
!----------------------------------------------------------------


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


!----------------------------------------------------------------
!-----  TASK: store fraction of ADC matrix on local node
!----------------------------------------------------------------


          case(SERVER_MATDIS)
!
! create local file name
!
          IF (MYPROC .LT. 10) THEN
            WRITE (locfn,'(A6,A1,I1)') nmbase,'.',MYPROC
          ELSEIF (MYPROC .LT. 100) THEN
            WRITE (locfn,'(A6,A1,I2)') nmbase,'.',MYPROC
          ELSEIF (MYPROC .LT. 1000) THEN
            WRITE (locfn,'(A6,A1,I3)') nmbase,'.',MYPROC
           ELSE
          CALL QUIT("NMPROC.GT.1000! EXTEND mpi_xmat_dist routine!")
          ENDIF
          write(*,*) 'Node',myproc,' opened ',locfn,' for writing.'
          open(unit=iobase,file=locfn, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='new')
!
! create local buffers, freed afterwards
!
          if(carith) then
            allocate(cbuf(intbuf))
            allocate(bufr(intbuf))
            allocate(bufi(intbuf))
          else
            allocate(buf(intbuf))
          endif
          allocate(ioi(intbuf))
          allocate(ioj(intbuf))

          jdummy = 0
!
! enter loop for receiving chunks of ADC matrix
!
          call interface_mpi_recv(nbufscur,1,MASTER,MSGN,
     &                           global_communicator)
          write(*,*) 'Slave',MYPROC,' will receive',nbufscur,
     &        ' batches from master.'
          write(iobase,ERR=889) nbufscur
          do i=1,nbufscur
            if(carith) then
              call interface_mpi_recv(bufr,intbuf,MASTER,MSGN,
     &                                global_communicator)
              call interface_mpi_recv(bufi,intbuf,MASTER,MSGN,
     &                                global_communicator)
              do ixx=1,intbuf
                cbuf(ixx) = dcmplx(bufr(ixx),bufi(ixx))
              enddo
            else
              call interface_mpi_recv(buf,intbuf,MASTER,MSGN,
     &                                global_communicator)
            endif
            call interface_mpi_recv(ioi,intbuf,MASTER,MSGN,
     &                              global_communicator)
            call interface_mpi_recv(ioj,intbuf,MASTER,MSGN,
     &                              global_communicator)
            call interface_mpi_recv(nact,1,MASTER,MSGN,
     &                              global_communicator)
            if(carith) then
              write(iobase,ERR=889) (cbuf(ixx),ixx=1,intbuf),
     &                              (ioi(ixx),ixx=1,intbuf),
     &                              (ioj(ixx),ixx=1,intbuf),
     &                              nact,jdummy
            else
              write(iobase,ERR=889) (buf(ixx),ixx=1,intbuf),
     &                              (ioi(ixx),ixx=1,intbuf),
     &                              (ioj(ixx),ixx=1,intbuf),
     &                              nact,jdummy
            endif
          enddo
          flush(iobase)
          close(iobase)
          write(*,*) 'Node',myproc,' has closed local ADC file.'
!
! deallocating local buffers
!
          if(carith) then
            deallocate(cbuf,bufr,bufi)
          else
            deallocate(buf)
          endif
          deallocate(ioi,ioj)

        end select

!----------------------------------------------------------------
      ENDDO    ! ----------- END OF SERVER MAIN LOOP ------------
!----------------------------------------------------------------


      return
!
!  error exit for write routine.
!
 889  call quit('problem on node with writing partial ADC matrix')
      END SUBROUTINE  ! mpi_slave_server
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#endif
