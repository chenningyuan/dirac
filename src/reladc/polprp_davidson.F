!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE mpi_davidson_master_ext_r(iw,
     &                                     iobase,
     &                                     desrep,
     &                                     intbuf,
     &                                     nmbase,
     &                                     filecnf,
     &                                     irecl,
     &                                     ladc,
     &                                     nmain)
!
!  the real Davidson diagonalizer written for serial and parallel purpose.
!  Parallel run: only master executes this, the slaves are occupied by performing
!                the partial matrix multiplication
!  Serial run:   Everything is done here.
!
!  In order to reduce the jungle of include files, we explicitly transfer the
!  relvant matrix parameters as arguments to the main routine.
!  In case of a parallel run the
!  main routine is responsible for communicating the relevant
!  parameters to the helper processes.
!
!  No separate serial Davidson version is kept
!
!  Important: In case of a parallel run the number of buffers of the partial
!  ADC matrices is the first integer in the file. This is -1 in the serial run since
!  the final size can only be determined at the end of the construction.
!
#if defined (VAR_MPI)
      use interface_to_mpi
#endif
      use polprp_cfg   ! make input variables and davidson control variables
                       ! accessible.

      IMPLICIT NONE
!
!  formal subroutine parameters
!
      INTEGER                      :: iw
      INTEGER                      :: iobase
      INTEGER                      :: desrep
      INTEGER                      :: intbuf
      CHARACTER*6                  :: nmbase
      CHARACTER*6                  :: filecnf
      INTEGER                      :: irecl
      INTEGER                      :: ladc
      INTEGER                      :: nmain
!
!  include files (will later be replaced by a suitable polprp_par).
!
#if defined (VAR_MPI)

#include  "polprp_servercodes.h"
#include  "../relccsd/ccpar.inc"

#endif
!
! *** control variables for Davidson
!
      INTEGER         :: nroots
      INTEGER         :: maxdavsp
      INTEGER         :: maxdavit
      REAL*8          :: convthr
!
!---------------Local variables----------------------------------------
!
!*** scalars (logical/integer/real/complex)

      logical         :: notallconv
      logical         :: proxflag
      logical         :: orthoviolflag

      integer         :: i,j,k,ialloc
      integer         :: nsiz,jvec,ibold_maxrec
      integer         :: iroot,locdim
      integer         :: lwork,info
      integer         :: idavit,iobase1,io_bnew,io_bold
      integer         :: ida_iostat

      real*8                                :: qmn,a1,ortho
      real*8,parameter                      :: autoev = 27.2113957D0
      real*8                                :: cput1,cput2
      real*8                                :: time1,time2
      real*8                                :: ttime,a

!*** arrays

      real*8, allocatable, dimension(:,:)   :: evint
      real*8, allocatable, dimension(:,:)   :: bstart
      real*8, allocatable, dimension(:,:)   :: aij
      real*8, allocatable, dimension(:,:)   :: aijbase
      real*8, allocatable, dimension(:)     :: lambda,reserr
      real*8, allocatable, dimension(:)     :: work
      real*8, allocatable, dimension(:)     :: qm,adiag,xi
      real*8, allocatable, dimension(:)     :: dmp1,aux
      integer, allocatable, dimension(:)    :: ievec
      logical, allocatable, dimension(:)    :: notconverged

!*** characters

      character*11                          :: davstartname
      character*11                          :: davevcname
      character*8                           :: dsroot ='DAVSTART'
      character*8                           :: dsroot2='ADCEVECS'
      character*8                           :: dsroot3='ADCSTATE'
      character*8                           :: dsroot4='DAVIMEDI'
      character*60                          :: listout
      character*70                          :: field
      character*10                          :: fileadc
!
!---------------Interface area --------------------------------------
!
      interface

        SUBROUTINE xmvmul_sp_r(i1,i2,i3,ra1,ra2)
          INTEGER                          :: i1,i2,i3
          real*8, dimension(:)             :: ra1,ra2
        END SUBROUTINE

        SUBROUTINE qmooccalc_r(i1,ra1,ra2,i2,i3,i4,i5)
          INTEGER                      :: i1,i2,i3,i4,i5
          real*8, dimension(:)         :: ra1
          real*8, dimension(:,:)       :: ra2
        END SUBROUTINE

        SUBROUTINE reortho_ext_r(i1,i2,i3)
          integer                          :: i1,i2,i3
        END SUBROUTINE

      end interface
!
!---------------Executable code----------------------------------------
!
      Call PST('Entering real (S/P) Davidson diagonalizer+')
!
! initialze davidson control variables
!
      nroots   =   polprp_davroots   !number of requested Davidson roots
      maxdavsp =   polprp_davmaxsp   !maximum dimension of Davidson space
      maxdavit =   polprp_davmaxit   !maximum number of Dav. iterations
      convthr  =   polprp_davconv    !convergence of Davidson eigenvectors
!
! we are on master. Therefore the local adc filename is already fixed
! 
#if defined (VAR_MPI)
      WRITE (fileadc,'(A6,A1,I1)') nmbase,'.',MASTER
#else
      WRITE (fileadc,'(A6)') nmbase
#endif
!
!  write out matrix and Davidson parameters (master only)
!
      WRITE(iw,*)
      WRITE(iw,*) ' --------------------------------------'
      WRITE(iw,*) ' ------ Davidson parameters -----------'
      WRITE(iw,*) ' --------------------------------------'
      WRITE(iw,*)
#if defined (VAR_MPI)
      WRITE (iw,*) '*** parallel run on',NMPROC,' nodes.'
#else
      WRITE (iw,*) '*** serial run'
#endif
      WRITE (iw,*) '*** all Davidson vectors are externalized.'
      WRITE(iw,*)
      WRITE(iw,*) 'Number of Davidson roots:         ',nroots
      WRITE(iw,*) 'Maximum Davidson subspace:        ',maxdavsp
      WRITE(iw,*) 'Number of Macro iterations:       ',maxdavit
      WRITE(iw,'(A,E14.5)') ' Eigenvalue convergence up to:   ',convthr
      WRITE(iw,*)
      WRITE(iw,*) ' ------------------------------------'
      WRITE(iw,*) ' ------ Matrix parameters -----------'
      WRITE(iw,*) ' ------------------------------------'
      WRITE(iw,*)
      WRITE(iw,*) 'Write buffer length is:           ',intbuf
      WRITE(iw,*) 'Final states have symmetry:       ',desrep
      WRITE(iw,*) 'ADC matrix dimension (# of rows): ',ladc
      WRITE(iw,*) 'Packed partial matrix stored in:  ',fileadc
      WRITE(iw,*) 'Configuration data stored in:     ',filecnf
      WRITE(iw,*) 'Record length of Config file:     ',irecl
      WRITE(iw,*) 'Size of main space:               ',nmain
      WRITE(iw,*)
!
!   allocate work space arrays
!
      lwork = 3*maxdavsp
      allocate(aij(maxdavsp,maxdavsp))     !Dav. overlap matrix (DOM)+eigenvectors
         aij = 0.0d0
      allocate(aijbase(maxdavsp,maxdavsp)) !Dav. overlap matrix base, calcd in macro itrs
         aijbase = 0.0d0
      allocate(lambda(maxdavsp))      !Eigenvalues of DOM
         lambda = 0.0d0
      allocate(reserr(maxdavsp))      !Residual errors of eigenvalues
         reserr = 0.0d0
      allocate(work(lwork))           !workspace for DSYEV
         work = 0.0d0
      allocate(qm(ladc))              !intermediate vector, also used as auxvector
         qm = 0.0d0
      allocate(adiag(ladc))           !vector of ADC diagonal elements
         adiag = 0.0d0
      allocate(xi(ladc))              !intermediate vector
         xi = 0.0d0
      allocate(ievec(ladc))
         ievec = 0
      allocate(dmp1(ladc))
         dmp1 = 0.0d0
      allocate(aux(ladc))
         aux = 0.0d0
      allocate(notconverged(maxdavsp))

      allocate(evint(ladc,nroots),stat=ialloc)
      if(ialloc.ne.0) then
        call quit('Alloc. problem (evint) in davidson_sp_r!')
      endif
      evint   = 0.0d0
!__________________________________________________________
!|
!| open file for start vectors and read them in .
!| select the first 1,...,nroots vectors from the file
!| which are the lowest ones in the main space.
!|
!| this is only done by MASTER in case of a parallel run!
!|
!     IF(desrep.GT.9) THEN
!       WRITE(davstartname,'(A8,A1,I2)') dsroot,'.',desrep
!     ELSE
!       WRITE(davstartname,'(A8,A2,I1)') dsroot,'.0',desrep
!     ENDIF

!     OPEN(iobase,file=davstartname,access='sequential',
!    &     status='unknown',form='unformatted')
!     READ(iobase) nsiz
!     WRITE(*,*) 'Found',nsiz,' Start vectors on file ',davstartname

!     If(2*nroots.gt.maxdavsp) then
!        Call Quit('Increase Davidson space!')
!     Endif

!     If(nroots.gt.nsiz) then
!       Write(*,*) 'You request more states than there',
!    &             ' are in the Main space.'
!       Write(*,*) 'Adapting requested roots',nroots,' to',nsiz
!       nroots = nsiz
!     Endif

!     allocate(bstart(nsiz,nroots))
!     bstart = 0.0d0

!     do i=1,nroots
!       Read(iobase) (bstart(j,i),j=1,nsiz)
!     enddo
!     write(iw,*) nroots,' start vectors read from ',davstartname
!     write(iw,*) 'Length of start vectors:',nsiz

!     Close(iobase,status='delete')

! check orthogonality of the start vectors

!     ortho = 0.0d0
!     Do i=1,nroots
!     Do j=1,nroots
!       a1 = dot_product(bstart(:,i),bstart(:,j))
!       ortho = ortho + abs(a1)
!     Enddo
!     Enddo
!     ortho = ortho - dble(nroots)
!     WRITE(*,*) 'Deviation from orthogonality of the',nroots,
!    &   ' Davidson start vectors:',ortho
!|
!|  end reading initial vectors.
!|___________________________________________________________

!__________________________________________________________
!|
!|   
!|   read matrix diagonal and open direct access files
!|   for ADC matrix, the old and the new Davidson vectors.
!|   File names for vector storage are hardcoded since they
!|   are temporary and nowhere else needed.
!|
!|   Write start vectors into file for old B-vectors. Then start
!|   macro iteration loop.
!|
      WRITE(iw,*)
      WRITE(iw,*) ' ------------------------------------'
      WRITE(iw,*) ' ------ File information  -----------'
      WRITE(iw,*) ' ------------------------------------'
      WRITE(iw,*)
      open(unit=iobase,file='ADCDGTMP', access='SEQUENTIAL',
     &     form='UNFORMATTED',status='unknown')
      read(iobase,ERR=888) (adiag(i),i=1,ladc)
      close(iobase,status='DELETE')
      write(iw,*) 'Diagonal elements read into memory.'

      Open(iobase,FILE=fileadc,FORM='UNFORMATTED',STATUS='UNKNOWN')
      write(iw,*) 'ADC file ',fileadc,' opened.'

      io_bnew = iobase + 2
      open(io_bnew,file='DAVXBN',access='direct',
     &     recl=ladc*8,status='unknown')

      io_bold = iobase + 3
      ibold_maxrec = 1
      open(io_bold,file='DAVXBO',access='direct',
     &     recl=ladc*8,status='unknown')
!
! now make start vectors the initial Davidson vectors for the first 
! macro iteration
! 
!     Do i=1,nroots
!       qm(1:nsiz) = bstart(1:nsiz,i)
!       qm(nsiz+1:ladc) = 0.0d0
!       write(io_bold,rec=i,iostat=ida_iostat) qm
!       if(ida_iostat.ne.0) call
!    &     quit('Error writing direct access file for bold')
!       if(i.gt.ibold_maxrec) ibold_maxrec = i
!     Enddo

!     deallocate(bstart)
!
! do the Davidson restart check
!
      WRITE(iw,*)
      WRITE(iw,*) ' ------------------------------------'
      WRITE(iw,*) ' ------ Restart information ---------'
      WRITE(iw,*) ' ------------------------------------'
      WRITE(iw,*)
      call restart_check(desrep,nroots,maxdavsp,ladc,
     &                   dsroot,dsroot2,io_bold)

      ibold_maxrec= nroots
!|
!|
!|___________________________________________________________
!

      WRITE(*,*) '  '
      orthoviolflag = .false.
!
! set timers to zero
!
      cput1 = 0.0d0
      cput2 = 0.0d0
      time1 = 0.0d0
      time2 = 0.0d0
!
! start time measurement
!
      call cpuused(cput1)

!__________________________________________________________
!|
!|
!|             DAVIDSON  MACRO ITERATION LOOP
!|
!*********************************************************
      Do idavit = 1,maxdavit
!*********************************************************


      Write(*,*) '  ****   MACRO iteration',idavit
      Write(*,*)
!
!      STEP  A.
!
      locdim = nroots    ! we start with the smallest dav space
      notconverged = .true.  !all vectors are at first not converged

      WRITE(*,*)
      WRITE(*,'(3X,A,4X,A,7X,A,7X,A,5X,A)') 'root','dim',
     &           'exc.energy','error','status'
      WRITE(*,'(70A1)') ('-',i=1,70)
!_____________________________________________________
!|
!| perform the initial block multiplication.
!| always necessary in a macro iteration.
!| each node performs the partial matrix vector mult.
!| due to MPI we do multiplication vector by vector
!|
!| at this stage of the code locdim is the number of 'old' vectors
!|

#if defined (VAR_MPI)
      do jvec=1,locdim
        call interface_mpi_bcast(SERVER_MATMUL,1,MASTER,
     &                           global_communicator)
        call interface_mpi_BCAST(ladc,1,MASTER,global_communicator)
!
! old:
!       xi(:) = bold(:,jvec)
 
! new:
        if(jvec.gt.ibold_maxrec) then
          call quit('bold record number exceeds maximum!')
        endif
        read(io_bold,rec=jvec,iostat=ida_iostat) xi
        if(ida_iostat.ne.0) call
     &     quit('Error reading direct access file for bold')

        call interface_mpi_bcast(xi,ladc,MASTER,global_communicator)
! master must do some work, too:
        call xmvmul_sp_r(iobase,intbuf,ladc,xi,qm) !qm is aux vector!
! master's result is stored in aux array qm
! important: MASTER participates identically in the reduce operation
! all the results are then collected in the array xi
        xi = 0.0d0
        call interface_mpi_reduce(qm,xi,ladc,
     &                 op_MPI_SUM,MASTER,global_communicator)
! as soon as the current new vector is available (in xi)
! we calculate its contribution to the aijbase-array

        do i=jvec,locdim
          if(i.gt.ibold_maxrec) then
            call quit('bold record number exceeds maximum!')
          endif
          read(io_bold,rec=i,iostat=ida_iostat) dmp1

          aijbase(i,jvec) = dot_product(dmp1(:),xi(:))
          aijbase(jvec,i) = aijbase(i,jvec)
        enddo

        write(io_bnew,rec=jvec) xi
      enddo
#else
      do jvec=1,locdim
        if(jvec.gt.ibold_maxrec) then
          call quit('bold record number exceeds maximum!')
        endif
        read(io_bold,rec=jvec,iostat=ida_iostat) dmp1
        call xmvmul_sp_r(iobase,intbuf,ladc,dmp1,xi)

        do i=jvec,locdim
          if(i.gt.ibold_maxrec) then
            call quit('bold record number exceeds maximum!')
          endif
          read(io_bold,rec=i,iostat=ida_iostat) dmp1
          aijbase(i,jvec) = dot_product(dmp1(:),xi(:))
          aijbase(jvec,i) = aijbase(i,jvec)
        enddo

        write(io_bnew,rec=jvec) xi

      enddo
#endif

!|
!| 
!|         END OF initial block multiplication.
!| 
!|____________________________________________________
    
!_______________________________________________________________
!|
!|         Davidson micro iteration loop
!|
      DO      ! Loop terminated by an exit statement
!
!  form  b_i*A*b_j and diagonalize. Matrix is replaced by eigenvectors.
!  Compiler knows length of matrix columns
!
        Do j=1,locdim
        Do i=1,locdim
          aij(i,j) = aijbase(i,j)
        Enddo
        Enddo

        CALL DSYEV('V','L',locdim,aij,maxdavsp,lambda,work,lwork,info)
        If(INFO.NE.0) Then
          Call Quit('Error in DSYEV (real Davidson!)')
!       else
!         write(*,*) 'Eigenvalues of aij matrix:'
!         do i=1,locdim
!           write(*,*) i,lambda(i)
!         enddo
        Endif
!
!  STEP B + C.  Form residual vector q_m for each eigenvalue k
!
!***************************** loop over the roots
        Do k = 1,nroots   ! loop over sought eigenvalues lambda_k
!*****************************

          proxflag = .false.

          Call QMOOCCALC_R(io_bnew,qm,aij,ladc,maxdavsp,locdim,k)

          evint(:,k) = qm(:)   !if converged, then we construct the normalized eigenvector.
          a1 = dot_product(evint(:,k),evint(:,k))
          a1 = sqrt(a1)
          evint(:,k) = evint(:,k)/a1

          Do i = 1,locdim
            if(i.gt.ibold_maxrec) then
              call quit('bold record number exceeds maximum!')
            endif
            read(io_bold,rec=i,iostat=ida_iostat) aux
            qm(:) = qm(:) - aij(i,k)*lambda(k)*aux(:)
          Enddo
          qmn = sqrt(dot_product(qm,qm))
          reserr(k) = qmn
         
          if(qmn.le.convthr) then
            notconverged(k) = .false.
          else
            notconverged(k) = .true.
          endif
          write(listout,'(A3,I3,A3,I5,6X,F10.7,6X,F10.7)') '   ',
     &           k,'   ',locdim,lambda(k),qmn
          if(notconverged(k).eqv..false.) then
            listout(47:48) =  ' *'
          else
            listout(47:48) =  ' o'
          endif
!
!  STEP D. form new preconditioner Xi
!  here one has to take care of the division!!!
!  It may happen that lambda(k) is approx. eq. to adiag(j)!
!
          Do j=1,ladc
            a1=lambda(k) - adiag(j)
            if(dabs(a1).lt.1.0e-03) then
              proxflag = .true.
              listout(49:51) =  ' PX'
              xi(j) = 1.0d0
            else
              xi(j)=qm(j)/a1
            endif
          Enddo
!         write(*,*) '   @Norm of preconditioner xi:',
!    &           sqrt(dot_product(xi(:),xi(:)))
!
!  STEP E. form  vector d_(M+1) for the kth eigenvalue
!  orthogonalize to *all* previous b_i directions from xi
!  Importance of new vector is decided on its norm !
!
          dmp1 = xi
          Do i = 1,locdim + k - 1
            if(i.gt.ibold_maxrec) then
              call quit('bold record number exceeds maximum!')
            endif
            read(io_bold,rec=i,iostat=ida_iostat) aux
            a1=dot_product(xi,aux)
            dmp1(:) = dmp1(:) - a1*aux(:)
          Enddo
          a1 = sqrt(dot_product(dmp1,dmp1))
!         write(*,*) '   @Norm of preconditioner d(M+1):',a1
!
!  STEP F. form  vector b_(M+1): b_(M+1)=d_(M+1)/||d_(M+1)||
!  this vector has to be included in the array of old vectors!
!  important. This is useful only if the vector b_(M+1) adds
!  new information (should be linearly independent)
!
          dmp1(:)=dmp1(:)/a1
          if(locdim+k.gt.ibold_maxrec) then
            ibold_maxrec = locdim + k
!           write(*,*) 'New record',ibold_maxrec,' written.'
          endif
! old:
!         bold(:,locdim+k)=dmp1(:)
! new:
          write(io_bold,rec=locdim+k,iostat=ida_iostat) dmp1(:)
          xi = dmp1
!
!  intermediate ortho check of new b vector.
!
          ortho=0.0d0
! old:
!         xi = bold(:,locdim+k)
! new:
          read(io_bold,rec=locdim+k,iostat=ida_iostat) xi
          Do i=1,locdim+k-1
            read(io_bold,rec=i,iostat=ida_iostat) aux
            a1 = dot_product(xi,aux)
            ortho = ortho + a1
          Enddo
          If(abs(ortho).lt.1.0E-06) then
            listout(52:54) = ' OK'
          else
            listout(52:54) = ' VI'
            orthoviolflag = .true.
          endif
!
!  STEP G.  Form A*b_(M+1)
!  this vector has to be included in the array of new vectors!
!
#if defined (VAR_MPI)

          call interface_mpi_bcast(SERVER_MATMUL,1,MASTER,
     &                             global_communicator)
          call interface_mpi_bcast(ladc,1,MASTER,global_communicator)
!         new b vector needs to be in xi and is still there (see above)
          call interface_mpi_bcast(xi,ladc,MASTER,global_communicator)
!  master's work (slaves work simultaneously)
          call xmvmul_sp_r(iobase,intbuf,ladc,xi,qm) !qm is aux vector!
          call interface_mpi_reduce(qm,xi,ladc,   !xi contains the final result
     &                 op_MPI_SUM,MASTER,global_communicator)
! final result should be in qm!
          qm = xi
#else

          call xmvmul_sp_r(iobase,intbuf,ladc,xi,qm) !qm contains result.

#endif
          WRITE(io_bnew,rec=locdim+k) qm(:)
!
!  STEP H.  update Aij to A_i,M+1,  i=1,...,M+1
!
          do i=1,locdim+k
            read(io_bold,rec=i,iostat=ida_iostat) aux
            aijbase(i,locdim+k) = dot_product(aux(:),qm(:))
            aijbase(locdim+k,i) = aijbase(i,locdim+k)
          enddo

          Write(*,*)  listout
!
! check orthogonality violation of newly added vectors.
! If it can not be maintained we have to exit in order to
! avoid intruder states.
!
          If(orthoviolflag.eqv..true.) EXIT   !leave this do loop


!*****************************  end loop over the roots
        Enddo  ! (k) was runnng index
!*****************************
        Write(*,*) '    '        !spacer

!
! check vector convergence in microiterations
!
        notallconv = .false.
        do k=1,nroots
          notallconv=notallconv.or.notconverged(k)
        enddo

        if(notallconv.eqv..false.) then
          write(*,*) 'All vectors converged. Leaving',
     &               ' microiteration loop.'
          exit
        endif

        If(orthoviolflag.eqv..true.) exit
      
        locdim = locdim + nroots
!
! check for macroiteration, transfer new start vector set.
! Attention!! It was observed that orthogonality slowly gets
! lost when macroiterations are entered. We therefore
! reorthogonalize the new start vector set.
!
        if( (locdim + nroots).gt.maxdavsp) then
          write(*,*) 'Davidson space exhausted.'
          write(*,*) 'Generating new (orthogonalized) start vectors.'
          write(*,*) 'For the next macroiteration.'
          do k=1,nroots
            if(k.gt.ibold_maxrec) then
!             write(*,*) 'New record',k,' written.'
              ibold_maxrec = k
            endif
            write(io_bold,rec=k,iostat=ida_iostat) evint(:,k)
          enddo
          call reortho_ext_r(io_bold,ladc,nroots)
          exit
        endif

      ENDDO
!|
!|         End Davidson micro iteration loop
!|
!|______________________________________________________________
!
!  we have to check convergence again to distinguish between next
!  macroiteration or end of davidson
!
      if(notallconv.eqv..false.) exit

      If(orthoviolflag.eqv..true.) exit

      ENDDO   !idavit loop
!|
!|         End Davidson MACRO iteration loop
!|
!|______________________________________________________________

      write(*,*) 'Leaving Davidson macro iteration loop.'

      If(orthoviolflag.eqv..true.) then
        write(*,*) 'Orthogonality of new vectors can not'
        write(*,*) 'be maintained any longer. Exiting...'
      Endif

      if(notallconv.eqv..false.) then
        write(*,*) 'All vectors converged to desired threshold'
      else
        write(*,*) '**WARNING: Not all vectors converged to desired'
        write(*,*) '**         Threshold. Check residual errors!'
      endif
!
!   close (local) ADC matrix file and delete it
!
      close(iobase,status='delete')
      write(*,*) 'Master node has deleted file ',fileadc
!
!   open configuration file for state analysis and writeout
!
      iobase1 = iobase + 1
      OPEN(iobase1,FILE=FILECNF,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')

!__________________  write out ** LONG ** results (binary)
!|
!| The format is number of states, main block size,
!| length of ADC matrix, then for each state:
!| state number, excitation energy (eV),
!| pole strength,  coefficients in the generic order.
!|
!|
      IF(desrep.GT.9) THEN
        WRITE(davevcname,'(A8,A1,I2)') dsroot2,'.',desrep
      ELSE
        WRITE(davevcname,'(A8,A2,I1)') dsroot2,'.0',desrep
      ENDIF
      OPEN(iobase,file=davevcname,access='sequential',
     &     status='unknown',form='unformatted')
!
      write(iobase) nroots
      write(iobase) nmain
      write(iobase) ladc
      Do K = 1,nroots
        write(iobase) K
        write(iobase) lambda(K)
        a1 = 0.0d0
        do j=1,nmain
          a1 = a1 + evint(j,k)*evint(j,k)
        enddo
        write(iobase) a1
        write(iobase) (evint(j,k),j=1,ladc)
      Enddo   !nroots
!
      CLOSE(iobase)
!|
!|
!|_______________________________________________
!
!_____ write out ** SHORT ** results (overview,formatted,sorted,main)
!|
!|
!|
      IF(desrep.GT.9) THEN
        WRITE(davevcname,'(A8,A1,I2)') dsroot3,'.',desrep
      ELSE
        WRITE(davevcname,'(A8,A2,I1)') dsroot3,'.0',desrep
      ENDIF
      OPEN(iobase,file=davevcname,access='sequential',
     &     status='unknown',form='formatted')

      Do K = 1,nroots
        a1 = 0.0d0
        do j=1,nmain
          a1 = a1 + evint(j,k)*evint(j,k)
        enddo

        Write(iobase,'(A14,X,I4)') 'Excited state:',K
        Write(iobase,'(A19)') '-----------------'
        Write(iobase,'(A12,F10.6,A4,F10.6,A8,F10.6,A5,F10.6,A5,I4)') 
     &   'Exc. energy:',
     &   lambda(K),' au ',lambda(K)*autoev,
     &   ' eV  PS:',a1,' Err:',reserr(k),' Sym:',desrep
        Write(iobase,'(A28,F15.6)') 'Norm^2 of double exc. contr:',
     &         1.0d0 - a1
        Write(iobase,'(A)') ' '
!
! determine 20 most relevant configurations
! reuse xi-array for modulus calculation and sort
! the corresponding eigenvector with respect to coeffs.
!
        do j=1,ladc;ievec(j)=j;enddo
        xi=-abs(evint(:,k))
        Call DISORT(ladc,xi,ievec)  !sorts from smallest to largest
        Do j=1,min(20,ladc)
          Read(iobase1,rec=ievec(j)) field
          a1 = evint(ievec(j),k)
          Write(iobase,'(A70,F15.6,3X,F15.6)') field,a1,a1*a1
        Enddo
        Write(iobase,'(A17)') '-----------------'
        Write(iobase,'(A)') ' '

      Enddo
      CLOSE(iobase)
!|
!|
!|_______________________________________________
!
!_____ write out ** MEDIUM ** results (formatted,main)
!|
!|
      if(desrep.GT.9) THEN
        write(davevcname,'(A8,A1,I2)') dsroot4,'.',desrep
      else
        write(davevcname,'(A8,A2,I1)') dsroot4,'.0',desrep
      endif
      open(iobase,file=davevcname,access='sequential',
     &     status='unknown',form='formatted')
      write(iobase,'(A)') 'Generic h/p eigenvector components'
      write(iobase,'(A,I5)') 'for symmetry ',desrep
      write(iobase,'(A)') ' '
      Do K = 1,nroots
        a1 = 0.0d0
        do j=1,nmain
          a1 = a1 + evint(j,k)*evint(j,k)
        enddo
        Write(iobase,'(A14,X,I4)') 'Excited state:',K
        Write(iobase,'(A19)') '-----------------'
        Write(iobase,'(A12,F10.6,A4,F10.6,A8,F10.6,A5,F10.6,A5,I4)') 
     &   'Exc. energy:',
     &   lambda(K),' au ',lambda(K)*autoev,
     &   ' eV  PS:',a1,' Err:',reserr(k),' Sym:',desrep
        do j=1,nmain
          read(iobase1,rec=j) field
          a1 = evint(j,k)
          Write(iobase,'(A70,F15.6,3X,F15.6)') field,a1,a1*a1
        enddo
        
        write(iobase,'(A)') ' '
      Enddo
      close(iobase)
!
!|
!|
!|
!|_______________________________________________
!

      CLOSE(iobase1)
!
!  deallocate all intermediate arrays
!
      deallocate(aij)
      deallocate(aijbase)
      deallocate(lambda)
      deallocate(reserr)
      deallocate(work)
      deallocate(qm)
      deallocate(adiag)
      deallocate(xi)
      deallocate(ievec)
      deallocate(dmp1)
      deallocate(aux)
      deallocate(notconverged)

      deallocate(evint)
!
!  delete random access files
!
      close(io_bnew,status='delete')
      close(io_bold,status='delete')
!
!  measure time
!
      call cpuused(cput2)
      ttime= cput2-cput1
      Write(*,*) '++++++++++++'
      Write(*,*) '++++++ time in davidson: ',ttime
      Write(*,*) '++++++++++++'

      RETURN

888   call quit('Error reading the diagonal matrix file')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE mpi_davidson_master_ext_c(iw,
     &                                     iobase,
     &                                     desrep,
     &                                     intbuf,
     &                                     nmbase,
     &                                     filecnf,
     &                                     irecl,
     &                                     ladc,
     &                                     nmain)
!
!  the complex Davidson diagonalizer written for serial and parallel purpose.
!  Parallel run: only master executes this, the slaves are occupied by performing
!                the partial matrix multiplication
!  Serial run:   Everything is done here.
!
!  In order to reduce the jungle of include files, we explicitly transfer the
!  relvant matrix parameters as arguments to the main routine.
!  In case of a parallel run the
!  main routine is responsible for communicating the relevant
!  parameters to the helper processes.
!
!  No separate serial Davidson version is kept
!
!  Important: In case of a parallel run the number of buffers of the partial
!  ADC matrices is the first integer in the file. This is -1 in the serial run since
!  the final size can only be determined at the end of the construction.
!
#if defined (VAR_MPI)
      use interface_to_mpi
#endif
      use polprp_cfg   ! make input variables and davidson control variables
                       ! accessible.

      IMPLICIT NONE
!
!  formal subroutine parameters
!
      INTEGER                      :: iw
      INTEGER                      :: iobase
      INTEGER                      :: desrep
      INTEGER                      :: intbuf
      CHARACTER*6                  :: nmbase
      CHARACTER*6                  :: filecnf
      INTEGER                      :: irecl
      INTEGER                      :: ladc
      INTEGER                      :: nmain
!
!  include files (will be replaced by a suitable polprp_par).
!
#if defined (VAR_MPI)

#include  "polprp_servercodes.h"
#include  "../relccsd/ccpar.inc"

#endif
!
! *** control variables for Davidson, initialized in input module

      INTEGER         :: nroots
      INTEGER         :: maxdavsp
      INTEGER         :: maxdavit
      REAL*8          :: convthr
!
!---------------Local variables----------------------------------------
!
!*** scalars (logical/integer/real/complex)

      logical                :: notallconv
      logical                :: proxflag
      logical                :: orthoviolflag
                             
      integer                :: i,j,k,ialloc
      integer                :: nsiz,jvec,ibold_maxrec
      integer                :: iroot,locdim
      integer                :: lwork,info
      integer                :: idavit,iobase1,io_bnew,io_bold
      integer                :: ida_iostat

      real*8                 :: qmn,a1,ortho !(real auxiliary variables)
      real*8,parameter       :: autoev = 27.2113957D0
      real*8                 :: cput1,cput2
      real*8                 :: time1,time2
      real*8                 :: ttime,a
      real*8,parameter       :: rzero=0.0d0

      complex*16,parameter   :: czero=(0.0d0,0.0d0)
      complex*16             :: cax  !(complex auxiliary variables)

!*** arrays

      complex*16, allocatable, dimension(:,:)   :: evint
      complex*16, allocatable, dimension(:,:)   :: bstart
      complex*16, allocatable, dimension(:,:)   :: aij
      complex*16, allocatable, dimension(:,:)   :: aijbase
      real*8, allocatable, dimension(:)         :: lambda, reserr
      complex*16, allocatable, dimension(:)     :: work
      real*8, allocatable, dimension(:)         :: rwork
      complex*16, allocatable, dimension(:)     :: qm, adiag, xi
      complex*16, allocatable, dimension(:)     :: dmp1, aux
      integer, allocatable, dimension(:)        :: ievec
      logical, allocatable, dimension(:)        :: notconverged
      real*8, allocatable, dimension(:)         :: mpiaux
      real*8, allocatable, dimension(:)         :: mpires1, mpires2

!*** characters

      character*11                 :: davstartname
      character*11                 :: davevcname
      character*8                  :: dsroot ='DAVSTART'
      character*8                  :: dsroot2='ADCEVECS'
      character*8                  :: dsroot3='ADCSTATE'
      character*8                  :: dsroot4='DAVIMEDI'
      character*60                 :: listout
      character*70                 :: field
      character*10                 :: fileadc
!
!---------------Interface area --------------------------------------
!
      interface

        FUNCTION pp_zdotc(i1,ca1,ca2)
          double complex                   :: pp_zdotc
          INTEGER                          :: i1
          double complex, dimension(:)     :: ca1,ca2
        END FUNCTION

        FUNCTION rzvecnorm(i1,ca1)
          real*8                           :: rzvecnorm
          INTEGER                          :: i1
          complex*16, dimension(:)         :: ca1
        END FUNCTION

        SUBROUTINE xmvmul_sp_c(i1,i2,i3,ca1,ca2)
          INTEGER                          :: i1,i2,i3
          complex*16, dimension(:)         :: ca1,ca2
        END SUBROUTINE

        SUBROUTINE qmooccalc_c(i1,ca1,ca2,i2,i3,i4,i5)
          INTEGER                          :: i1,i2,i3,i4,i5
          complex*16, dimension(:)         :: ca1
          complex*16, dimension(:,:)       :: ca2
        END SUBROUTINE

        SUBROUTINE reortho_ext_c(i1,i2,i3)
          integer                          :: i1,i2,i3
        END SUBROUTINE

      end interface
!
!---------------Executable code----------------------------------------
!
      Call PST('Entering complex (S/P) Davidson diagonalizer+')
!
! initialze davidson control variables
!
      nroots   =   polprp_davroots   !number of requested Davidson roots
      maxdavsp =   polprp_davmaxsp   !maximum dimension of Davidson space
      maxdavit =   polprp_davmaxit   !maximum number of Dav. iterations
      convthr  =   polprp_davconv    !convergence of Davidson eigenvectors
!
! we are on master. Therefore the local adc filename is already fixed
! 
#if defined (VAR_MPI)
      WRITE (fileadc,'(A6,A1,I1)') nmbase,'.',MASTER
#else
      WRITE (fileadc,'(A6)') nmbase
#endif
!
!  write out matrix and Davidson parameters (master only)
!
      WRITE(iw,*)
      WRITE(iw,*) ' --------------------------------------'
      WRITE(iw,*) ' ------ Davidson parameters -----------'
      WRITE(iw,*) ' --------------------------------------'
      WRITE(iw,*)
#if defined (VAR_MPI)
      WRITE (iw,*) '*** parallel run on',NMPROC,' nodes.'
#else
      WRITE (iw,*) '*** serial run'
#endif
      WRITE (iw,*) '*** all Davidson vectors are externalized.'
      WRITE(iw,*)
      WRITE(iw,*) 'Number of Davidson roots:         ',nroots
      WRITE(iw,*) 'Maximum Davidson subspace:        ',maxdavsp
      WRITE(iw,*) 'Number of Macro iterations:       ',maxdavit
      WRITE(iw,'(A,E14.5)') ' Eigenvalue convergence up to:   ',convthr
      WRITE(iw,*)
      WRITE(iw,*) ' ------------------------------------'
      WRITE(iw,*) ' ------ Matrix parameters -----------'
      WRITE(iw,*) ' ------------------------------------'
      WRITE(iw,*)
      WRITE(iw,*) 'Write buffer length is:           ',intbuf
      WRITE(iw,*) 'Final states have symmetry:       ',desrep
      WRITE(iw,*) 'ADC matrix dimension (# of rows): ',ladc
      WRITE(iw,*) 'Packed partial matrix stored in:  ',fileadc
      WRITE(iw,*) 'Configuration data stored in:     ',filecnf
      WRITE(iw,*) 'Record length of Config file:     ',irecl
      WRITE(iw,*) 'Size of main space:               ',nmain
      WRITE(iw,*)
!
!   allocate and clear work space arrays
!
      lwork = 3*maxdavsp
      allocate(aij(maxdavsp,maxdavsp))     !Dav. overlap matrix (DOM)+eigenvectors
           aij = czero
      allocate(aijbase(maxdavsp,maxdavsp)) !Dav. overlap matrix base, calcd in macro itrs
           aijbase = czero
      allocate(lambda(maxdavsp))      !Eigenvalues of (hermitian) DOM (real)
           lambda = rzero
      allocate(reserr(maxdavsp))      !Residual errors of eigenvalues
           reserr = rzero
      allocate(work(lwork))           !workspace for ZHEEV
           work = czero
      allocate(rwork(lwork))          !workspace for ZHEEV
           rwork = rzero
      allocate(qm(ladc))              !intermediate vector, also used as auxvector
           qm = czero
      allocate(adiag(ladc))           !vector of ADC diagonal elements
           adiag = czero
      allocate(xi(ladc))              !intermediate vector
           xi = czero
      allocate(dmp1(ladc))
           dmp1 = czero
      allocate(aux(ladc))
           aux = czero
      allocate(ievec(ladc))
           ievec = 0
      allocate(notconverged(maxdavsp))
           notconverged = .true.
      allocate(mpiaux(ladc))       !helper array for MPI complex number communication
           mpiaux = rzero
      allocate(mpires1(ladc))       !helper array for MPI complex number communication
           mpires1 = rzero
      allocate(mpires2(ladc))       !helper array for MPI complex number communication
           mpires2 = rzero
 

      allocate(evint(ladc,nroots),stat=ialloc)
      if(ialloc.ne.0) then
        call quit('Alloc. problem (evint) in mpi_davidson_master_c!')
      endif
      evint = czero
!
! set all timers to zero
!
      cput1 = 0.0d0
      cput2 = 0.0d0
      time1 = 0.0d0
      time2 = 0.0d0
!
! start time measurement here
!
      call cpuused(cput1)
!__________________________________________________________
!|
!| open file for start vectors and read them in .
!| select the first 1,...,nroots vectors from the file
!| which are the lowest ones in the main space.
!|
!| this is only done by MASTER in case of a parallel run!
!|
      IF(desrep.GT.9) THEN
        WRITE(davstartname,'(A8,A1,I2)') dsroot,'.',desrep
      ELSE
        WRITE(davstartname,'(A8,A2,I1)') dsroot,'.0',desrep
      ENDIF

      OPEN(iobase,file=davstartname,access='sequential',
     &     status='unknown',form='unformatted')
      READ(iobase) nsiz
      WRITE(*,*) 'Found',nsiz,' Start vectors on file ',davstartname

      If(2*nroots.gt.maxdavsp) then
         Call Quit('Increase Davidson space!')
      Endif

      If(nroots.gt.nsiz) then
        Write(*,*) 'You request more states than there',
     &             ' are in the Main space.'
        Write(*,*) 'Adapting requested roots',nroots,' to',nsiz
        nroots = nsiz
      Endif

      allocate(bstart(nsiz,nroots))
      bstart = czero

      do i=1,nroots
        Read(iobase) (bstart(j,i),j=1,nsiz)
      enddo
      write(iw,*) nroots,' start vectors read from ',davstartname
      write(iw,*) 'Length of start vectors:',nsiz

      Close(iobase,status='delete')
!
! check convergence of the start vectors
! 
      ortho = 0.0d0
      Do i=1,nroots
      Do j=1,nroots
        cax = pp_zdotc(nsiz,bstart(:,i),bstart(:,j))
        ortho = ortho + abs(cax)
      Enddo
      Enddo
      ortho = ortho - dble(nroots)
      WRITE(*,*) 'Deviation from orthogonality of the',nroots,
     &   ' Davidson start vectors:',ortho
!|
!|  end reading initial vectors.
!|___________________________________________________________

      WRITE(*,*) '  '
      orthoviolflag = .false.
!__________________________________________________________
!|
!|   
!|   read matrix diagonal and open direct access files
!|   for ADC matrix, the old and the new Davidson vectors.
!|   File names for vector storage are hardcoded since they
!|   are temporary and nowhere else needed.
!|
!|   Write start vectors into file for old B-vectors. Then start
!|   macro iteration loop.
!|
      open(unit=iobase,file='ADCDGTMP', access='SEQUENTIAL',
     &     form='UNFORMATTED',status='unknown')
      read(iobase,ERR=888) (adiag(i),i=1,ladc)
      close(iobase,status='DELETE')
      write(iw,*) 'Diagonal elements read into memory.'

      Open(iobase,FILE=fileadc,FORM='UNFORMATTED',STATUS='UNKNOWN')
      write(iw,*) 'ADC file ',fileadc,' opened.'

      io_bnew = iobase + 2
      open(io_bnew,file='DAVXBN',access='direct',
     &     recl=ladc*16,status='unknown')

      io_bold = iobase + 3
      ibold_maxrec = 1
      open(io_bold,file='DAVXBO',access='direct',
     &     recl=ladc*16,status='unknown')
!
! now make start vectors the initial Davidson vectors for the first
! macro iteration
! 
      Do i=1,nroots
        qm(1:nsiz) = bstart(1:nsiz,i)
        qm(nsiz+1:ladc) = czero
        write(io_bold,rec=i,iostat=ida_iostat) qm
        if(ida_iostat.ne.0) call
     &     quit('Error writing direct access file for bold')
        if(i.gt.ibold_maxrec) ibold_maxrec = i
      Enddo

      deallocate(bstart)
!|
!|
!|___________________________________________________________

!__________________________________________________________
!|
!|
!|             DAVIDSON  MACRO ITERATION LOOP
!|
!*********************************************************
      Do idavit = 1,maxdavit
!*********************************************************

      Write(*,*) '  ****   MACRO iteration',idavit
      Write(*,*)
!
!      STEP  A.
!
      locdim = nroots    ! we start with the smallest dav space
      notconverged = .true.  !all vectors are at first not converged

      WRITE(*,*)
      WRITE(*,'(3X,A,4X,A,7X,A,7X,A,5X,A)') 'root','dim',
     &           'exc.energy','error','status'
      WRITE(*,'(70A1)') ('-',i=1,70)

!_____________________________________________________
!|
!| perform the initial block multiplication.
!| always necessary in a macro iteration.
!| each node performs the partial matrix vector mult.
!| due to MPI we do multiplication vector by vector
!|
!| at this stage of the code locdim is the number of 'old' vectors
!|

#if defined (VAR_MPI)
      do jvec=1,locdim
        call interface_mpi_bcast(SERVER_MATMUL,1,MASTER,
     &                           global_communicator)
        call interface_mpi_BCAST(ladc,1,MASTER,global_communicator)

! old:
!       xi(:) = bold(:,jvec)
! new:
        if(jvec.gt.ibold_maxrec) then
          call quit('bold record number exceeds maximum!')
        endif
        read(io_bold,rec=jvec,iostat=ida_iostat) xi
        if(ida_iostat.ne.0) call
     &     quit('Error reading direct access file for bold')


        mpiaux = real(xi)  ! communicate real part of xi to array mpiaux
        call interface_mpi_bcast(mpiaux,ladc,MASTER,global_communicator)
        mpiaux = aimag(xi)  ! communicate imaginary part of xi to array mpiaux
        call interface_mpi_bcast(mpiaux,ladc,MASTER,global_communicator)

        call xmvmul_sp_c(iobase,intbuf,ladc,xi,qm)

        mpires1 = rzero
        mpiaux  = real(qm)
        call interface_mpi_reduce(mpiaux,mpires1,ladc,
     &                 op_MPI_SUM,MASTER,global_communicator)
        mpires2 = rzero
        mpiaux  = aimag(qm)
        call interface_mpi_reduce(mpiaux,mpires2,ladc,
     &                 op_MPI_SUM,MASTER,global_communicator)
        do i=1,ladc
          xi(i) = dcmplx(mpires1(i), mpires2(i))
        enddo
        
        do i=jvec,locdim
          if(i.gt.ibold_maxrec) then
            call quit('bold record number exceeds maximum!')
          endif
          read(io_bold,rec=i,iostat=ida_iostat) dmp1

          aijbase(i,jvec) = pp_zdotc(ladc,dmp1(:),xi(:))
          aijbase(jvec,i) = dconjg(aijbase(i,jvec))
        enddo

        write(io_bnew,rec=jvec) xi
      enddo
#else
      do jvec=1,locdim
        if(jvec.gt.ibold_maxrec) then
          call quit('bold record number exceeds maximum!')
        endif
        read(io_bold,rec=jvec,iostat=ida_iostat) dmp1
        call xmvmul_sp_c(iobase,intbuf,ladc,dmp1,xi)

        do i=jvec,locdim
          if(i.gt.ibold_maxrec) then
            call quit('bold record number exceeds maximum!')
          endif
          read(io_bold,rec=i,iostat=ida_iostat) dmp1
          aijbase(i,jvec) = pp_zdotc(ladc,dmp1(:),xi(:))
          aijbase(jvec,i) = dconjg(aijbase(i,jvec))
        enddo

        write(io_bnew,rec=jvec) xi

      enddo
#endif
!|
!| 
!|         END OF initial block multiplication.
!| 
!|____________________________________________________
    
!_______________________________________________________________
!|
!|         Davidson micro iteration loop
!|
      DO      ! MICRO iteration loop. Terminated by an exit statement
!
!  (b_i)^* * A * b_j is formed and is diagonalized. Matrix is replaced by eigenvectors.
!  AIJ should be hermitian
!
        Do j=1,locdim
        Do i=1,locdim
          aij(i,j) = aijbase(i,j)
        Enddo
        Enddo
!
!  original declaration
!
!  zheev (JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO)
!

        CALL ZHEEV('V','L',locdim,aij,maxdavsp,lambda,
     &             work,lwork,rwork,info)
        if(info.ne.0) then
          write(*,*) 'ZHEEV problem in mpi_davidson_c! code:',info
          call quit('numerical problem')
!       else
!         write(*,*) 'Eigenvalues of aij matrix:'
!         do i=1,locdim
!           write(*,*) i,lambda(i)
!         enddo
        endif
!
!  STEP B + C.  Form residual vector q_m for each eigenvalue k
!
!***************************** loop over the roots
        Do k = 1,nroots   ! loop over sought eigenvalues lambda_k
!*****************************

          proxflag = .false.

          Call QMOOCCALC_C(io_bnew,qm,aij,ladc,maxdavsp,locdim,k)

          evint(:,k) = qm(:)   !if converged, then we construct the normalized eigenvector.
          a1 = rzvecnorm(ladc,evint(:,k))   !real number!
          evint(:,k) = evint(:,k)/dcmplx(a1,0.0d0)

          Do i = 1,locdim
            if(i.gt.ibold_maxrec) then
              call quit('bold record number exceeds maximum!')
            endif
            read(io_bold,rec=i,iostat=ida_iostat) aux
            qm(:) = qm(:) - aij(i,k)*dcmplx(lambda(k),0.0d0)*aux(:)
          Enddo
          qmn = rzvecnorm(ladc,qm)
          reserr(k) = qmn
         
          if(qmn.le.convthr) then
            notconverged(k) = .false.
          else
            notconverged(k) = .true.
          endif
          write(listout,'(A3,I3,A3,I5,6X,F10.7,6X,F10.7)') '   ',
     &           k,'   ',locdim,lambda(k),qmn
          if(notconverged(k).eqv..false.) then
            listout(47:48) =  ' *'
          else
            listout(47:48) =  ' o'
          endif
!
!  STEP D. form new preconditioner Xi
!  here one has to take care of the division!!!
!  It may happen that lambda(k) is approx. eq. to adiag(j)!
!
          Do j=1,ladc
            a1=lambda(k) - real(adiag(j))
            if(dabs(a1).lt.1.0e-03) then
              proxflag = .true.
              listout(49:51) =  ' PX'
              xi(j) = (1.0d0,0.0d0)
            else
              xi(j)=qm(j)/dcmplx(a1,0.0d0)
            endif
          Enddo
!
!  STEP E. form  vector d_(M+1) for the kth eigenvalue
!  orthogonalize to *all* previous b_i directions from xi
!  Importance of new vector is decided on its norm !
!
          dmp1 = xi
          Do i = 1,locdim + k - 1
            if(i.gt.ibold_maxrec) then
              call quit('bold record number exceeds maximum!')
            endif
            read(io_bold,rec=i,iostat=ida_iostat) aux
            cax=pp_zdotc(ladc,aux,xi)
            dmp1(:) = dmp1(:) - cax*aux(:)
          Enddo
!
!  STEP F. form  vector b_(M+1): b_(M+1)=d_(M+1)/||d_(M+1)||
!  this vector has to be included in the array of old vectors!
!  important. This is useful only if the vector b_(M+1) adds
!  new information (should be linearly independent)
!
          a1 = rzvecnorm(ladc,dmp1)
          dmp1(:)=dmp1(:)/dcmplx(a1,0.0d0)
          if(locdim+k.gt.ibold_maxrec) then
            ibold_maxrec = locdim + k
!           write(*,*) 'New record',ibold_maxrec,' written.'
          endif
!old:
!         bold(:,locdim+k)=dmp1(:)
!new:
          write(io_bold,rec=locdim+k,iostat=ida_iostat) dmp1(:)
!
!  intermediate ortho check of new b vector.
!
          ortho=0.0d0
          read(io_bold,rec=locdim+k,iostat=ida_iostat) xi
          Do i=1,locdim+k-1
            read(io_bold,rec=i,iostat=ida_iostat) aux
            cax = pp_zdotc(ladc,xi,aux)
            ortho = ortho + abs(cax)
          Enddo
          If(abs(ortho).lt.1.0E-06) then
            listout(52:54) = ' OK'
          else
            listout(52:54) = ' VI'
            orthoviolflag = .true.
          endif
!
!  STEP G.  Form A*b_(M+1)
!  this vector has to be included in the array of new vectors!
!

#if defined (VAR_MPI)
          call interface_mpi_bcast(SERVER_MATMUL,1,MASTER,
     &                             global_communicator)
          call interface_mpi_bcast(ladc,1,MASTER,global_communicator)
!         xi(:) = bold(:,locdim+k)  !xi still contains this (see above)
          mpiaux = real(xi)  ! communicate real part of xi to array mpiaux
          call interface_mpi_bcast(mpiaux,ladc,MASTER,
     &                             global_communicator)
          mpiaux = aimag(xi)  ! communicate imaginary part of xi to array mpiaux
          call interface_mpi_bcast(mpiaux,ladc,MASTER,
     &                             global_communicator)

          call xmvmul_sp_c(iobase,intbuf,ladc,xi,qm)

          mpires1 = rzero
          mpiaux  = real(qm)
          call interface_mpi_reduce(mpiaux,mpires1,ladc,
     &                   op_MPI_SUM,MASTER,global_communicator)
          mpires2 = rzero
          mpiaux  = aimag(qm)
          call interface_mpi_reduce(mpiaux,mpires2,ladc,
     &                   op_MPI_SUM,MASTER,global_communicator)
          do i=1,ladc
            xi(i) = dcmplx(mpires1(i), mpires2(i))
          enddo
! final result should be in qm!
          qm = xi
#else

          call xmvmul_sp_c(iobase,intbuf,ladc,xi,qm)

#endif
          WRITE(io_bnew,rec=locdim+k) qm(:)
!
!  STEP H.  update Aij to A_i,M+1,  i=1,...,M+1
!
          do i=1,locdim+k
            read(io_bold,rec=i,iostat=ida_iostat) aux
            aijbase(i,locdim+k) = pp_zdotc(ladc,aux(:),qm(:))
            aijbase(locdim+k,i) = dconjg(aijbase(i,locdim+k))
          enddo

          Write(*,*)  listout
!
! check orthogonality violation of newly added vectors.
! If it can not be maintained we have to exit in order to
! avoid intruder states.
!
          If(orthoviolflag.eqv..true.) EXIT   !leave this do loop


!*****************************  end loop over the roots
        Enddo  ! (k) was runnng index
!*****************************
        Write(*,*) '    '        !spacer

!
! check vector convergence in microiterations
!
        notallconv = .false.
        do k=1,nroots
          notallconv=notallconv.or.notconverged(k)
        enddo

        if(notallconv.eqv..false.) then
          write(*,*) 'All vectors converged. Leaving',
     &               ' microiteration loop.'
          exit
        endif

        If(orthoviolflag.eqv..true.) exit
      
        locdim = locdim + nroots
!
! check for macroiteration, transfer new start vector set.
! Attention!! It was observed that orthogonality slowly gets
! lost when macroiterations are entered. We therefore
! reorthogonalize the new start vector set.
!
        if( (locdim + nroots).gt.maxdavsp) then
          write(*,*) 'Davidson space exhausted.'
          write(*,*) 'Generating new (orthogonalized) start vectors.'
          write(*,*) 'For the next macroiteration.'
          do k=1,nroots
            if(k.gt.ibold_maxrec) then
!             write(*,*) 'New record',k,' written.'
              ibold_maxrec = k
            endif
            write(io_bold,rec=k,iostat=ida_iostat) evint(:,k)
          enddo
          call reortho_ext_c(io_bold,ladc,nroots)
          exit
        endif

      ENDDO
!|
!|         End Davidson micro iteration loop
!|
!|______________________________________________________________
!
!  we have to check convergence again to distinguish between next
!  macroiteration or end of davidson
!
      if(notallconv.eqv..false.) exit

      If(orthoviolflag.eqv..true.) exit

      ENDDO   !idavit loop

      write(*,*) 'Leaving Davidson macro iteration loop.'

      If(orthoviolflag.eqv..true.) then
        write(*,*) 'Orthogonality of new vectors can not'
        write(*,*) 'be maintained any longer. Exiting...'
      Endif

      if(notallconv.eqv..false.) then
        write(*,*) 'All vectors converged to desired threshold'
      else
        write(*,*) '**WARNING: Not all vectors converged to desired'
        write(*,*) '**         Threshold. Check residual errors!'
      endif
!|
!|         End Davidson MACRO iteration loop
!|
!|______________________________________________________________
!
!   close (local) ADC matrix file and delete it
!
      close(iobase,status='delete')
      write(*,*) 'Master node has deleted file ',fileadc
!
!   open configuration file
!
      iobase1 = iobase + 1
      OPEN(iobase1,FILE=FILECNF,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
!
!__________________ loop over the roots and write out **LONG** results
!|
!| In contrast to earlier versions we will no longer write out the ASCII
!| numbers since huge unused files are generated, otherwise. Instead we need
!| the data as binary (format='unformatted').
!|
!| The format is number of states, main block size,
!| length of ADC matrix, then for each state:
!| state number, excitation energy (eV),
!| pole strength,  coefficients in the generic h/p,2h/2p order.
!|
      IF(desrep.GT.9) THEN
        WRITE(davevcname,'(A8,A1,I2)') dsroot2,'.',desrep
      ELSE
        WRITE(davevcname,'(A8,A2,I1)') dsroot2,'.0',desrep
      ENDIF
      OPEN(iobase,file=davevcname,access='sequential',
     &     status='unknown',form='unformatted')
!
      write(iobase) nroots
      write(iobase) nmain
      write(iobase) ladc
      Do K = 1,nroots
        write(iobase) K
        write(iobase) lambda(K)
        a1 = 0.0d0
        do j=1,nmain
          a1 = a1 + real(dconjg(evint(j,k))*evint(j,k))
        enddo
        write(iobase) a1
        write(iobase) (evint(j,k),j=1,ladc)
      Enddo   !nroots
!
      CLOSE(iobase)
!|
!|
!|_______________________________________________
!
!__________________ loop over the roots and write out ** SHORT ** results
!|
!|
!|
      IF(desrep.GT.9) THEN
        WRITE(davevcname,'(A8,A1,I2)') dsroot3,'.',desrep
      ELSE
        WRITE(davevcname,'(A8,A2,I1)') dsroot3,'.0',desrep
      ENDIF
      OPEN(iobase,file=davevcname,access='sequential',
     &     status='unknown',form='formatted')

      Do K = 1,nroots
!
! determine pole strength
!
        a1 = 0.0d0
        do j=1,nmain
          a1 = a1 + real(dconjg(evint(j,k))*evint(j,k))
        enddo

        Write(iobase,'(A14,X,I4)') 'Excited state:',K
        Write(iobase,'(A19)') '-----------------'
        Write(iobase,'(A12,F10.6,A4,F10.6,A8,F10.6,A5,F10.6,A5,I4)') 
     &   'Exc. energy:',
     &   lambda(K),' au ',lambda(K)*autoev,
     &   ' eV  PS:',a1,' Err:',reserr(k),' Sym:',desrep
        Write(iobase,'(A28,F15.6)') 'Norm^2 of double exc. contr:',
     &         1.0d0 - a1
        Write(iobase,'(A)') ' '
!
! determine 20 most relevant configurations
! reuse xi-array for modulus calculation and sort
! the corresponding eigenvector with respect to coeffs.
! mpiaux serves as auxiliary here
!
        do j=1,ladc;ievec(j)=j;enddo
        mpiaux=-abs(evint(:,k))
        Call DISORT(ladc,mpiaux,ievec)  !sorts from smallest to largest
        Do j=1,min(20,ladc)
          Read(iobase1,rec=ievec(j)) field
          cax = evint(ievec(j),k)
          Write(iobase,'(A70,2F15.6,3X,F15.6)') field,cax,
     &                  real(dconjg(cax)*cax)
        Enddo
        Write(iobase,'(A17)') '-----------------'
        Write(iobase,'(A)') ' '

      Enddo
      CLOSE(iobase)
!|
!|
!|_______________________________________________
!
!________________________________________________
!|
!|
!|  write out results as MEDIUM output (symmetry specific)
!|
      if(desrep.GT.9) THEN
        write(davevcname,'(A8,A1,I2)') dsroot4,'.',desrep
      else
        write(davevcname,'(A8,A2,I1)') dsroot4,'.0',desrep
      endif
      open(iobase,file=davevcname,access='sequential',
     &     status='unknown',form='formatted')
      write(iobase,'(A)') 'Generic h/p eigenvector components'
      write(iobase,'(A,I5)') 'for symmetry ',desrep
      write(iobase,'(A)') ' '
      Do K = 1,nroots
        a1 = 0.0d0
        do j=1,nmain
          a1 = a1 + real(dconjg(evint(j,k))*evint(j,k))
        enddo
        Write(iobase,'(A14,X,I4)') 'Excited state:',K
        Write(iobase,'(A19)') '-----------------'
        Write(iobase,'(A12,F10.6,A4,F10.6,A8,F10.6,A5,F10.6,A5,I4)') 
     &   'Exc. energy:',
     &   lambda(K),' au ',lambda(K)*autoev,
     &   ' eV  PS:',a1,' Err:',reserr(k),' Sym:',desrep
        do j=1,nmain
          read(iobase1,rec=j) field
          cax = evint(j,k)
          Write(iobase,'(A70,2F15.6,3X,F15.6)') field,cax,
     &                  real(dconjg(cax)*cax)
        enddo
        write(iobase,'(A)') ' '
      Enddo
      close(iobase)
!
!|
!|
!|
!|_______________________________________________
!

      CLOSE(iobase1)
!
!  deallocate all intermediate arrays
!
      deallocate(aij)
      deallocate(aijbase)
      deallocate(lambda)
      deallocate(reserr)
      deallocate(work)
      deallocate(rwork)
      deallocate(qm)
      deallocate(adiag)
      deallocate(xi)
      deallocate(dmp1)
      deallocate(aux)
      deallocate(ievec)
      deallocate(notconverged)
      deallocate(mpiaux,mpires1,mpires2)

      deallocate(evint)
!
!  delete direct access file
!
      close(io_bnew,status='delete')
      close(io_bold,status='delete')
!
!  measure time
!
      call cpuused(cput2)
      ttime= cput2-cput1
      Write(*,*) '++++++++++++'
      Write(*,*) '++++++ time in davidson: ',ttime
      Write(*,*) '++++++++++++'

      RETURN

 888  call quit('Error reading the diagonal matrix file')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&
!&&
!&&   routines needed in serial and parallel 
!&&
!&&
!&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      real*8 function rzvecnorm(n,z)
      integer                     :: n
      complex*16, dimension (:)   :: z

      integer          :: i
      complex*16       :: sum,x

      sum = (0.0d0,0.0d0)
      do i=1,n
        x = z(i)
        sum = sum + dconjg(x) * x
      enddo
      if(abs(aimag(sum)).gt.1.0E-10) then
        write(*,*) 'Too large imaginary part in rzvecnorm!'
        call quit('numerical problem')
      endif
      rzvecnorm = sqrt(real(sum))

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE REORTHO_EXT_R(io_bold,ladc,nroots)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      INTEGER                        :: io_bold,ladc,nroots
!
!------------------------ local variables -------------------------
!
      REAL*8,  allocatable, dimension(:)     :: x,ax
      REAL*8                                 :: a1
      INTEGER                                :: i,k,ida_iostat

!
!---------------Executable ----------------------------------------
!
! reorthogonalize nvec vectors stored columnwise in the array a
!
      allocate(x(ladc))
      allocate(ax(ladc))

      Do k=1,nroots
        read(io_bold,rec=k,iostat=ida_iostat) x
        Do i=1,k-1
          read(io_bold,rec=i,iostat=ida_iostat) ax
          a1 = dot_product(x(:),ax(:))
          x = x - a1*ax(:)
        Enddo
        a1 = sqrt(dot_product(x,x))
        ax  =  x/a1
        write(io_bold,rec=k,iostat=ida_iostat) ax
      Enddo

      deallocate(x)
      deallocate(ax)
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE REORTHO_EXT_C(io_bold,ladc,nroots)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      INTEGER                        :: io_bold,ladc,nroots
!
!------------------------ local variables -------------------------
!
      complex*16, allocatable, dimension(:)   :: x,ax
      complex*16                              :: ca5
      real*8                                  :: rn
      INTEGER                                 :: i,k,ida_iostat
!
!------------------------ interfaces -------------------------------
!
      interface

        FUNCTION pp_zdotc(i1,ca1,ca2)
          double complex                   :: pp_zdotc
          INTEGER                          :: i1
          double complex, dimension(:)     :: ca1,ca2
        END FUNCTION

        FUNCTION rzvecnorm(i1,ca1)
          real*8                           :: rzvecnorm
          integer                          :: i1
          double complex, dimension(:)     :: ca1
        END FUNCTION

      end interface
!
!---------------Executable ----------------------------------------
!
! reorthogonalize nroots vectors stored columnwise in the array v
!
      allocate(x(ladc))
      allocate(ax(ladc))

      Do k=1,nroots
        read(io_bold,rec=k,iostat=ida_iostat) x
        Do i=1,k-1
          read(io_bold,rec=i,iostat=ida_iostat) ax
          ca5 = pp_zdotc(ladc,x(:),ax(:))
          x = x - ca5*ax(:)
        Enddo
        rn = rzvecnorm(ladc,x)
        ax = x/dcmplx(rn,0.0d0)
        write(io_bold,rec=k,iostat=ida_iostat) ax
      Enddo

      deallocate(x)
      deallocate(ax)
      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      double complex function pp_zdotc(n,zx,zy)

      integer                           :: n
      double complex, dimension(:)      :: zx,zy
!
!   local variables
!
      double complex                    :: a
      integer                           :: i
      double complex,parameter          :: czero = (0.0d0,0.0d0)
!
!   executable code
!
      a = czero

      do i=1,n
        a = a + dconjg(zx(i)) * zy(i)
      enddo

      pp_zdotc = a

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE xmvmul_sp_r(iobase,intbuf,ladc,PO,PN)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase,intbuf,ladc
      real*8,dimension(:)            :: PO,PN
!
!------------------------ local variables -------------------------
!
      REAL*8,  allocatable, dimension(:)      :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      REAL*8                                  :: A
      integer                                 :: nbufsloc
      integer                                 :: ioresult

!
!---------------Description--------------------------------------------
!
!  Perform the real symmetric multiplication of the partial or complete
!  ADC matrix with PO and store the result in the vector PN
!
!---------------Executable code--------------------------------------
!
!  for the matrix elements we note that only the lower triangular
!  part *AND* the diagonal of the ADC-matrix is stored. 
!  the if takes care of not double counting the diagonal.

      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
!
!  nbufsloc is an indicator if the ADC matrix is complete or partial
!  In serial mode it is complete and has -1 at the beginning. Here we need the
!  iostat mechanism to determine the file end
!
!  In parallel mode nbufsloc contains the actual number of buffers. This is known due
!  to the matrix distribution procedure.
!
!  By this we avoid an awkward handling of the number of buffers as arguments, or common blocks...
!
      pn = 0.0d0

      rewind(iobase)
      read(iobase) nbufsloc

      If(nbufsloc.lt.0) then  !any negative number indicates serial
!___________________________________________________________
!|
!|   serial case, exact buffer number is not known
!|   we read stream blocks until end of file is reached.

      Do

        READ(iobase,iostat=ioresult) (BUF(IXX),IXX=1,INTBUF),
     &                               (IOI(IXX),IXX=1,INTBUF),
     &                               (IOJ(IXX),IXX=1,INTBUF),
     &                               NACT,JDUMMY
        if (ioresult > 0) then
           call quit('ADC matrix read error (real). Terminating')
        endif
        if (ioresult < 0) then  ! EOF detected
           exit    !leave infinite do loop
        endif

! useful data in buffer

        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
        ENDDO

      Enddo
!|
!|   
!|__________________________________________________________
      else
!___________________________________________________________
!|
!|   parallel case, exact buffer number is known
!|
        DO I = 1,nbufsloc
          READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                         (IOI(IXX),IXX=1,INTBUF),
     &                         (IOJ(IXX),IXX=1,INTBUF),
     &                          NACT,JDUMMY
          DO K = 1, NACT
            IROW = IOI(K)
            ICOL = IOJ(K)
            A    = BUF(K)

            PN(IROW) = PN(IROW) + PO(ICOL) * A
            IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
          ENDDO
        ENDDO
!|
!|   
!|__________________________________________________________

      endif

      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 889  CALL QUIT('Matrix read error in xmvmul_sp_r')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE xmvmul_sp_c(iobase,intbuf,ladc,PO,PN)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase,intbuf,ladc
      complex*16,dimension(:)        :: PO,PN
!
!------------------------ local variables -------------------------
!
      complex*16, allocatable, dimension(:)   :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      complex*16                              :: A
      integer                                 :: nbufsloc
      integer                                 :: ioresult

!
!---------------Description--------------------------------------------
!
!  Perform the hermitian multiplication of the external partial 
!  ADC matrix with PO and store the result in the vector PN
!
!---------------Executable code--------------------------------------
!
!  for the matrix elements we note that only the lower triangular
!  part *AND* the diagonal of the ADC-matrix is stored. 
!  the if takes care of not double counting the diagonal.

      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
!
!  nbufsloc is an indicator if the ADC matrix is complete or partial
!  In serial mode it is complete and has -1 at the beginning. Here we need the
!  iostat mechanism to determine the file end
!
!  In parallel mode nbufsloc contains the actual number of buffers. This is known due
!  to the matrix distribution procedure.
!
!  By this we avoid an awkward handling of the number of buffers as arguments, or common blocks...
!
      pn = (0.0d0,0.0d0)

      rewind(iobase)
      read(iobase) nbufsloc

      If(nbufsloc.lt.0) then  !any negative number indicates serial
!___________________________________________________________
!|
!|   serial case, exact buffer number is not known
!|   we read stream blocks until end of file is reached.

        Do

          READ(iobase,iostat=ioresult) (BUF(IXX),IXX=1,INTBUF),
     &                                 (IOI(IXX),IXX=1,INTBUF),
     &                                 (IOJ(IXX),IXX=1,INTBUF),
     &                                 NACT,JDUMMY
          if (ioresult > 0) then
             call quit('ADC matrix read error (complex). Terminating')
          endif
          if (ioresult < 0) then  ! EOF detected
             exit    !leave infinite do loop
          endif

! useful data in buffer

          DO K = 1, NACT
            IROW = IOI(K)
            ICOL = IOJ(K)
            A    = BUF(K)

             PN(IROW) = PN(IROW) + PO(ICOL) * A
            IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW)*dconjg(A)
          ENDDO

        Enddo
!|
!|   
!|__________________________________________________________
      else
!___________________________________________________________
!|
!|   parallel case, exact buffer number is known
!|
        DO I = 1,nbufsloc
          READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                         (IOI(IXX),IXX=1,INTBUF),
     &                         (IOJ(IXX),IXX=1,INTBUF),
     &                          NACT,JDUMMY
          DO K = 1, NACT
            IROW = IOI(K)
            ICOL = IOJ(K)
            A    = BUF(K)

            PN(IROW) = PN(IROW) + PO(ICOL) * A
            IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW)*dconjg(A)
          ENDDO
        ENDDO
!|
!|   
!|__________________________________________________________

      endif

      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 889  CALL QUIT('Matrix read error in xmvmul_sp_c')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XDIAG_INI_R(iobase,intbuf,nbufs,ladc)

!  Fetch the diagonal entries of the sparse matrix A, needed for the
!  Davidson iterations.
!  This routine is only called in the serial case and obtains NBUFS.
!  It generates a local file handle and writes out diagonal.

      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase
      integer                        :: intbuf,nbufs,ladc
!
!------------------------ local variables -------------------------
!
      Character*8,parameter                   :: diagfile='ADCDGTMP'
      integer                                 :: serialindicator
      REAL*8,  allocatable, dimension(:)      :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      REAL*8,  allocatable, dimension(:)      :: adiag
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      REAL*8                                  :: A
      INTEGER                                 :: io_out

!
!--------------- Interfaces --------------------------------------------
!
      interface

        Function get_file_unit()
          integer                      :: get_file_unit
        END Function

      end interface
!
!---------------Executable code--------------------------------------
!
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      allocate(adiag(ladc))

!  .. read elements

      REWIND(iobase)
      read(iobase) serialindicator   !discarded, but should be less than zero!

      DO I = 1,NBUFS
        READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          IF(IROW.EQ.ICOL) adiag(irow) = A

        ENDDO
      ENDDO

      io_out = get_file_unit()
      open(unit=io_out,file=diagfile, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='NEW')
      write(io_out,ERR=887) (adiag(ixx),ixx=1,ladc)
      close(io_out)

      deallocate(adiag)
      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 887  CALL QUIT('Matrix write error in XDIAG_INI_R')
 889  CALL QUIT('Matrix read error in XDIAG_INI_R')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XDIAG_INI_C(iobase,intbuf,nbufs,ladc)

! complex extension of XDIAG_INI_R

      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase
      integer                        :: intbuf,nbufs,ladc
!
!------------------------ local variables -------------------------
!
      Character*8,parameter                     :: diagfile='ADCDGTMP'
      integer                                   :: serialindicator
      complex*16, allocatable, dimension(:)     :: BUF
      INTEGER, allocatable, dimension(:)        :: IOI,IOJ
      complex*16, allocatable, dimension(:)     :: adiag
      INTEGER                                   :: I,K,IXX
      INTEGER                                   :: NACT,JDUMMY,IROW,ICOL
      complex*16                                :: A
      INTEGER                                   :: io_out
!
!--------------- Interfaces --------------------------------------------
!
      interface

        Function get_file_unit()
          integer                      :: get_file_unit
        END Function

      end interface
!
!---------------Executable code--------------------------------------
!
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      allocate(adiag(ladc))

c  .. read elements

      REWIND(iobase)
      read(iobase) serialindicator   !discarded, but should be less than zero!

      DO I = 1,NBUFS
        READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          IF(IROW.EQ.ICOL) adiag(irow) = A

        ENDDO
      ENDDO

      io_out = get_file_unit()
      open(unit=io_out,file=diagfile, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='NEW')
       write(io_out,ERR=887) (adiag(ixx),ixx=1,ladc)
       close(io_out)

      deallocate(adiag)
      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 887  CALL QUIT('Matrix write error in XDIAG_INI_C')
 889  CALL QUIT('Matrix read error in XDIAG_INI_C')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE fulldia_r3(io,fname,intbuf,ladc,desrep)

      implicit none

      INTEGER                       :: io,intbuf,ladc,desrep
      CHARACTER(len=6)              :: fname
!
!---------------Description--------------------------------------------
!
!  For not too large symmetric matrices the full diagonalization is
!  performed for debugging purposes. This adapted version calculates
!  eigenvalues and eigenvectors and writes them to file "FULLEVC.XX"
!
!---------------Local variables--------------------------------------
!
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      REAL*8, allocatable             :: A(:,:),WORK(:),EVL(:)
      REAL*8, allocatable             :: BUF(:)
      INTEGER, allocatable            :: IOI(:),IOJ(:)
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO,LWORK
      integer                         :: ioresult,serialindicator
      character*10                    :: fullevcname
      Real*8                          :: D0
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!
      IF(LADC.gt.2000) THEN
        write(*,*) 'Full diagonalization skipped'
        return
      ENDIF

      CALL PST('Entering POLPRP REAL FULL diagonalizer+')

      D0 = 0.0d0

      open(unit=io,file=fname, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='unknown')

      REWIND(io)
      read(io) serialindicator


      allocate(a(ladc,ladc))
      allocate(evl(ladc))
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      lwork=10*ladc
      allocate(work(lwork))
      
      A=0.0d0
!___________________________________________________________
!|
!|   serial case, exact buffer number is not known
!|   we read stream blocks until end of file is reached.

      Do

        READ(io,iostat=ioresult) (BUF(IXX),IXX=1,INTBUF),
     &                           (IOI(IXX),IXX=1,INTBUF),
     &                           (IOJ(IXX),IXX=1,INTBUF),
     &                            NACT,JDUMMY
        if (ioresult > 0) then
           call quit('ADC matrix read error (real). Terminating')
        endif
        if (ioresult < 0) then  ! EOF detected
           exit    !leave infinite do loop
        endif

! useful data in buffer

        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A(irow,icol) = BUF(K)
        ENDDO

      Enddo

      close(io)
!|
!|
!|__________________________________________________________

!
!  call symmetric diagonalizer
!  calculate only eigenvalues being placed in A
!
      JOBZ='V'
      UPLO='L'

      CALL DSYEV( JOBZ, UPLO, LADC, A, LADC, EVL, WORK, LWORK, INFO )

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in DSYEV. Info:',INFO
           CALL QUIT ('*** Error in fulldia_r3  ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 20):'
      DO I=1,min(20,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO
!
!  write full set of eigenvectors to file
!  Keep hereby the *same* format as in the ADCEVECS file in order
!  to make comparison more readable
!
      IF(desrep.GT.9) THEN
        WRITE(fullevcname,'(A8,I2)') 'FULLEVC.',desrep
      ELSE
        WRITE(fullevcname,'(A9,I1)') 'FULLEVC.0',desrep
      ENDIF
      open(io,file=fullevcname,access='sequential',
     &     status='unknown',form='unformatted')
      do k=1,ladc
!       write(*,*) 'Wrinting full eigenvector #',k
        write(io) k
        write(io) evl(k)
        write(io) D0
        write(io) (A(i,k),i=1,ladc)
      enddo
      close(io)

!
!  release all memory
!
      deallocate(a)
      deallocate(evl)
      deallocate(buf)
      deallocate(ioi)
      deallocate(ioj)
      deallocate(work)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE fulldia_c3(io,fname,intbuf,ladc)

      implicit none

      INTEGER                       :: io,intbuf,ladc
      CHARACTER(len=6)              :: fname
!
!---------------Description--------------------------------------------
!
!  For not too large Hermitian matrices the full diagonalization is
!  performed for debugging purposes. The required arrays are allocated
!  dynamically in order to achieve maximal independence.
!
!  The only specific entities are two available file handles,
!  the length of the write chunks and the length of the ADC matrix.  
!
!---------------Local variables--------------------------------------
!
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      COMPLEX*16,   allocatable       :: A(:,:),CWORK(:),BUF(:)
      REAL*8,       allocatable       :: EVL(:),RWORK(:)
      INTEGER,      allocatable       :: IOI(:),IOJ(:)
      INTEGER                         :: LWORK
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO
      INTEGER                         :: ioresult,serialindicator
!
!---------------Executable code--------------------------------------
!
!  print info
!
!     write(*,*) 'file name: ',fname
!     write(*,*) 'intbuf:',intbuf
!     write(*,*) 'ladc:',ladc
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!


      IF(LADC.gt.2000) THEN                                
        write(*,*) 'Full diagonalization skipped'
        return                                             
      ENDIF

      CALL PST('Entering POLPRP FULL complex diagonalizer+')

      open(unit=io,file=fname, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='unknown')

      REWIND(io)
      read(io) serialindicator

      allocate(a(ladc,ladc))
      allocate(evl(ladc))
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      lwork=10*ladc
      allocate(rwork(lwork))
      allocate(cwork(lwork))
      
      A=(0.0d0,0.0d0)
!_________________________ reading section
!|
      Do

        READ(io,iostat=ioresult) (BUF(IXX),IXX=1,INTBUF),
     &                           (IOI(IXX),IXX=1,INTBUF),
     &                           (IOJ(IXX),IXX=1,INTBUF),
     &                            NACT,JDUMMY
        if (ioresult > 0) then
           call quit('ADC matrix read error (complex). Terminating')
        endif
        if (ioresult < 0) then  ! EOF detected
!          write(*,*) 'EOF detected.'
           exit    !leave infinite do loop
        endif

! useful data in buffer

        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A(irow,icol) = BUF(K)
        ENDDO

      Enddo

      close(io)

!|________________________ end reading section

!
!  call hermitian diagonalizer
!  calculate only eigenvalues being placed in A
!
      JOBZ='N'
      UPLO='L'

      CALL ZHEEV( JOBZ, UPLO, LADC, A, LADC, EVL, CWORK,
     &            LWORK, RWORK, INFO )

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in ZHEEV. Info:',INFO
           CALL QUIT ('*** Error in fulldia_c3 ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 20):'
      DO I=1,min(20,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO
!
!  release all memory
!
      deallocate(a)
      deallocate(evl)
      deallocate(buf)
      deallocate(ioi)
      deallocate(ioj)
      deallocate(rwork)
      deallocate(cwork)

      RETURN

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QMOOCCALC_R(ioda,qmu,aij,ladc,maxdavsp,locdim,iroot)

      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                      :: ioda,ladc,maxdavsp,locdim,iroot
      real*8,dimension(:)          :: qmu
      real*8,dimension(:,:)        :: aij
!
!------------------------ local variables -------------------------
!
      INTEGER                                 :: i,j
      real*8, allocatable, dimension(:)       :: aux
!
!---------------Description--------------------------------------------
!
!   This routine calculates the qm vector using the new b-vectors
!   from disk
!
      allocate(aux(ladc))
      qmu = 0.0d0
      do i=1,locdim
        read(ioda,rec=i,err=888) aux
        qmu(:) = qmu(:) + aij(i,iroot)*aux(:)
      enddo
      deallocate(aux)
      return

 888  call quit('Vector read error in QMOOCCALC_R!')

      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QMOOCCALC_C(ioda,qmu,aij,ladc,maxdavsp,locdim,iroot)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                      :: ioda,ladc,maxdavsp,locdim,iroot
      complex*16, dimension(:)     :: qmu
      complex*16, dimension(:,:)   :: aij
!
!------------------------ local variables -------------------------
!
      INTEGER                                     :: i,j
      complex*16, allocatable, dimension(:)       :: aux
!
!---------------Description--------------------------------------------
!
!   This routine calculates the qm vector using the new b-vectors
!   from disk
!
      allocate(aux(ladc))
      qmu = (0.0d0,0.0d0)
      do i=1,locdim
        read(ioda,rec=i,err=888) aux
        qmu(:) = qmu(:) + aij(i,iroot)*aux(:)
      enddo
      deallocate(aux)
      return

 888  call quit('Vector read error in QMOOCCALC_C!')

      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DISORT(n,arr,brr)
!
! adapted from numerical recipes in fortran 77
! to simultaneously shuffle the index array brr
! with the sorted array arr
!
      implicit none
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      INTEGER brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,atemp
      INTEGER b,btemp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        atemp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=atemp
        btemp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=btemp
        if(arr(l+1).gt.arr(ir))then
          atemp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=atemp
          btemp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=btemp
        endif
        if(arr(l).gt.arr(ir))then
          atemp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=atemp
          btemp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=btemp
        endif
        if(arr(l+1).gt.arr(l))then
          atemp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=atemp
          btemp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=btemp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        atemp=arr(i)
        arr(i)=arr(j)
        arr(j)=atemp
        btemp=brr(i)
        brr(i)=brr(j)
        brr(j)=btemp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)then
          stop 'increase NSTACK in DISORT'
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE restart_check(desrep,nroots,maxdavsp,ladc,
     &                         hp_name,ev_name,io_bold)

      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!   This routine checks the existence of some restart files and tries to fit
!   found restart vectors into the vector array. The order of importance is:
!
!   ADCEVECS.XX
!   DAVSTART.XX
!   Unity matrix
!
!------------------------ calling variables -------------------------
!
      integer                    :: desrep,nroots,maxdavsp,ladc
      character*8                :: hp_name,ev_name
      integer                    :: io_bold
!
!------------------------ local variables -------------------------
!
      logical                                :: isthere
      integer                                :: fileavail,io
      integer                                :: i,j,nsiz
      integer                                :: ida_iostat
      Character*11                           :: filename_d,filename_e
      Real*8                                 :: ortho,a1
      Real*8,allocatable,dimension(:)        :: qm
      Real*8,allocatable,dimension(:,:)      :: bstart
!
!------------------------ interface -------------------------
!
      Interface

        function get_file_unit()
           Integer    :: get_file_unit
        end function

      End interface
!
!------------------------ executable code -------------------------
!
      fileavail = 0
!
! check availability of DAVSTART.XX file
!
      IF(desrep.GT.9) THEN
        WRITE(filename_d,'(A8,A1,I2)') hp_name,'.',desrep
      ELSE
        WRITE(filename_d,'(A8,A2,I1)') hp_name,'.0',desrep
      ENDIF
      inquire(file=filename_d,exist=isthere)
      if(isthere) fileavail = fileavail + 1
!
! check availability of ADCEVECS.XX file
!
      IF(desrep.GT.9) THEN
        WRITE(filename_e,'(A8,A1,I2)') ev_name,'.',desrep
      ELSE
        WRITE(filename_e,'(A8,A2,I1)') ev_name,'.0',desrep
      ENDIF
      inquire(file=filename_e,exist=isthere)
      if(isthere) fileavail = fileavail + 1
!
!  read corresponding file or pack in unity matrix
!
      io = get_file_unit()
      Select case(fileavail)

!------------------------------------------------------------

      case(2)                 !  ADCEVECS is present

      write(*,*) 'Using ',filename_e,' as restart file.'

!------------------------------------------------------------

      case(1)                 !  DAVSTART is present

      write(*,*) 'Using ',filename_d,' as restart file.'
      OPEN(io,file=filename_d,access='sequential',
     &     status='unknown',form='unformatted')
      READ(io) nsiz
      WRITE(*,*) 'Found',nsiz,' Start vectors.'

      If(2*nroots.gt.maxdavsp) then
         Call Quit('Increase Davidson space!')
      Endif

      If(nroots.gt.nsiz) then
        Write(*,*) 'You request more states than there',
     &             ' are in the Main space.'
        Write(*,*) 'Adapting requested roots',nroots,' to',nsiz
        nroots = nsiz
      Endif

      allocate(bstart(nsiz,nroots))
      bstart = 0.0d0

      do i=1,nroots
        Read(io) (bstart(j,i),j=1,nsiz)
      enddo
      write(*,*) nroots,' start vectors read.'
      write(*,*) 'Length of start vectors:',nsiz

      Close(io)

! check orthogonality of the start vectors

      ortho = 0.0d0
      Do i=1,nroots
      Do j=1,nroots
        a1 = dot_product(bstart(:,i),bstart(:,j))
        ortho = ortho + abs(a1)
      Enddo
      Enddo
      ortho = ortho - dble(nroots)
      WRITE(*,*) 'Orthogonality:',ortho
!
! now make start vectors the initial Davidson vectors for the first 
! macro iteration
! 
      allocate(qm(ladc))
      Do i=1,nroots
        qm(1:nsiz) = bstart(1:nsiz,i)
        qm(nsiz+1:ladc) = 0.0d0
        write(io_bold,rec=i,iostat=ida_iostat) qm
        if(ida_iostat.ne.0) call
     &     quit('Error writing direct access file for bold')
      Enddo

      deallocate(qm)
      deallocate(bstart)

!------------------------------------------------------------

      case(0)                 !  no file. Unity matrix
      write(*,*) 'No restart file found. Using unity matrix.'

!------------------------------------------------------------

      End select


      write(*,*) 'Restart check finalized.'
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
