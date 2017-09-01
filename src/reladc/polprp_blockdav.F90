!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine mpi_block_davidson_r(iobase, &
                                      intbuf, &
                                      ladc, &
                                      nmain, &
                                      desrep, &
                                      nmbase)
!
!  real Davidson diagonalizer for serial and parallel runs.
!
#if defined (VAR_MPI)
      use interface_to_mpi
#endif
      use polprp_cfg   ! make davidson control variables accessible.

      Implicit none
!
!  formal parameters of the subroutine
!
      integer                               :: iobase
      integer                               :: intbuf
      integer                               :: ladc
      integer                               :: nmain
      integer                               :: desrep
      character*6                           :: nmbase
!
!  control variables for Davidson
!
      integer                               :: nroots
      integer                               :: maxdavsp
      integer                               :: maxdavit
      real*8                                :: convthr

#if defined (VAR_MPI)

#include  "polprp_servercodes.h"
#include  "../relccsd/ccpar.inc"

#endif
!
!  local variables
!
      real*8,parameter                      :: rzero = 0.0d0
      real*8,parameter                      :: prethr = 0.001d0
      real*8                                :: resnorm
      real*8                                :: prenorm1,prenorm2
      real*8                                :: ocheck
      real*8,allocatable,dimension(:)       :: lambda
      real*8,allocatable,dimension(:)       :: bv,bw
      real*8,allocatable,dimension(:)       :: prc
      real*8,allocatable,dimension(:)       :: pro
      real*8,allocatable,dimension(:)       :: work
      real*8,allocatable,dimension(:)       :: adiag
      real*8,allocatable,dimension(:,:)     :: res
      real*8,allocatable,dimension(:,:)     :: hbv
      real*8,allocatable,dimension(:,:)     :: hij,hijcal

      integer                               :: i,j,krydim
      integer                               :: io_bold
      integer                               :: lwork,info
      integer                               :: kroot
      integer                               :: ida_iostat
      integer                               :: iallocstat

      logical                               :: allconverged
      logical,allocatable,dimension(:)      :: converged

      character*10                          :: fileadc
      character*60                          :: listout
!
!  interfaces for subroutines
!
      interface

        SUBROUTINE init_krylov_r(i1,i2,i3,i4,i5)
          integer                          :: i1,i2,i3,i4,i5
        END SUBROUTINE

        SUBROUTINE xmvmul_sp_r(i1,i2,i3,ra1,ra2)
          integer                          :: i1,i2,i3
          real*8, dimension(:)             :: ra1,ra2
        END SUBROUTINE

        SUBROUTINE calc_resvec_r(i1,i2,i3,i4,ra1,ra2,ra3,ra4)
          integer                      :: i1,i2,i3,i4
          real*8, dimension(:,:)       :: ra1
          real*8, dimension(:,:)       :: ra2
          real*8, dimension(:)         :: ra3
          real*8, dimension(:,:)       :: ra4
        END SUBROUTINE

        SUBROUTINE calc_precon_r(i1,ra1,r1,ra2,ra3)
          integer                      :: i1
          real*8, dimension(:)         :: ra1
          real*8                       :: r1
          real*8, dimension(:)         :: ra2
          real*8, dimension(:)         :: ra3
        END SUBROUTINE

        SUBROUTINE calc_preort_r(i1,i2,i3,ra1,ra2)
          integer                      :: i1,i2,i3
          real*8, dimension(:)         :: ra1
          real*8, dimension(:)         :: ra2
        END SUBROUTINE

        SUBROUTINE write_davevecs_r(i1,i2,i3,i4,i5,i6,ra1,ra2)
          integer                      :: i1,i2,i3,i4,i5,i6
          real*8, dimension(:)         :: ra1
          real*8, dimension(:,:)       :: ra2
        END SUBROUTINE

      end interface
!
! say hello
!
      Call PST('Entering real (S/P) Davidson diagonalizer+')
!
! initialize davidson control variables from POLPRP input
!
      nroots   =   polprp_davroots   !number of requested Davidson roots
      maxdavsp =   polprp_davmaxsp   !maximum dimension of Davidson space
      maxdavit =   polprp_davmaxit   !maximum number of Dav. iterations
      convthr  =   polprp_davconv    !convergence of Davidson eigenvectors
!
! allocate and initialize arrays according to calculation dimensions
!
      lwork = 3*maxdavsp
      allocate(hij(maxdavsp,maxdavsp))
         hij = rzero
      allocate(hijcal(maxdavsp,maxdavsp))
         hijcal = rzero
      allocate(lambda(maxdavsp))
         lambda = 0.0d0
      allocate(bv(ladc))
         bv = rzero
      allocate(bw(ladc))
         bw = rzero
      allocate(res(ladc,nroots))
         res = rzero
      allocate(prc(ladc))
         prc = rzero
      allocate(pro(ladc))
         pro = rzero
      allocate(work(lwork))
         work = rzero
      allocate(adiag(ladc))
         adiag = rzero
      allocate(converged(nroots))
         converged = .false.

!*******************************************************
! print essential information to the user
!*******************************************************

      WRITE(*,*)
      WRITE(*,*) ' --------------------------------------'
      WRITE(*,*) ' ------ Davidson parameters -----------'
      WRITE(*,*) ' --------------------------------------'
      WRITE(*,*)
#if defined (VAR_MPI)
      WRITE (*,*) '*** parallel run on',NMPROC,' nodes.'
#else
      WRITE (*,*) '*** serial run'
#endif
      WRITE(*,*)
      WRITE(*,*) ' ------->> symmetry:              ',desrep
      WRITE(*,*) 'ADC matrix dimension (# of rows): ',ladc
      WRITE(*,*) 'Size of main space:               ',nmain
      WRITE(*,*)
      WRITE(*,*) 'Number of Davidson roots:         ',nroots
      WRITE(*,*) 'Maximum Davidson subspace:        ',maxdavsp
      WRITE(*,*) 'Number of Macro iterations:       ',maxdavit
      WRITE(*,'(A,E14.5)') ' Eigenvalue convergence up to:   ',convthr
      WRITE(*,*)
      WRITE(*,*) ' --------------------------------------'
      WRITE(*,*) ' ------ Progress of calculation:'
      WRITE(*,*) ' --------------------------------------'
      WRITE(*,*)


!*******************************************************
! Start Davidson run. First read ADC matrix diagonal
!*******************************************************

      open(unit=iobase,file='ADCDGTMP', access='SEQUENTIAL', &
           form='UNFORMATTED',status='unknown')
      read(iobase) (adiag(i),i=1,ladc)
      close(iobase,status='DELETE')
      write(*,*) 'Diagonal elements read'

!*******************************************************
! Open main ADC matrix file. Remains open during cycles
!*******************************************************

#if defined (VAR_MPI)
      write (fileadc,'(A6,A1,I1)') nmbase,'.',MASTER
#else
      write (fileadc,'(A6)') nmbase
#endif
      open(iobase,FILE=fileadc,FORM='UNFORMATTED',STATUS='UNKNOWN')
      write(*,*) 'ADC file ',fileadc,' opened.'

!*******************************************************
! Open direct access file for current Krylov space
!*******************************************************

      io_bold = iobase + 1
      open(io_bold,file='DAVXBO',access='direct', &
           recl=ladc*8,status='unknown')
      write(*,*) 'Direct access file for Krylov space opened.'

!*******************************************************
! Fill Krylov space with start vectors
! Initial dimension is equal to number of roots
!*******************************************************

      call init_krylov_r(io_bold,desrep,nroots,maxdavsp,ladc)
      krydim = nroots

      WRITE(*,*)
      WRITE(*,'(3X,A,4X,A,7X,A,7X,A,5X,A)') 'root','dim',  &
                 'exc.energy','error','status'
      WRITE(*,'(70A1)') ('-',i=1,70)
      allconverged = .false.
!_______________________________________________________________
!|
!| Microiteration Davidson loop until all roots converged or the
!| Maximum Davidson space is reached.
!|
!|
      DO WHILE (allconverged.eqv..false.)

      allocate(hbv(ladc,krydim),stat=iallocstat) ! start dimension only!
         hbv = rzero
      if(iallocstat.ne.0) stop 'allocation problem'

!*******************************************************
! Step X: Form subspace Hamiltonian,n current Krylov
! space dimension:   H * b_old is stored
! Will go away: orthogonality check of current space
!*******************************************************

      ocheck = 0.0d0
      do j=1,krydim

        read(io_bold,rec=j,iostat=ida_iostat) bw
        call xmvmul_sp_r(iobase,intbuf,ladc,bw(:),hbv(:,j))

        do i=j,krydim

          read(io_bold,rec=i) bv
          hij(i,j) = dot_product(bv,hbv(:,j))
          hij(j,i) = hij(i,j)

!         ocheck=ocheck + dot_product(bw,bv)

        enddo

      enddo
!     ocheck = ocheck - dble(krydim)
!     write(*,*) 'Orthogonality dev. of Krylov space:',ocheck

!*******************************************************
! Step X+Y: Diagonalize current subspace Hamiltonian
! hijcal is copy, contains eigenvectors after DSYEV call
!*******************************************************

      hijcal = hij
      call dsyev('V','L',krydim,hijcal,maxdavsp,lambda,work,lwork,info)
      if(info.ne.0) then
        call quit('Error in DSYEV (real Davidson!)')
      else
!       write(*,*) 'Eigenvalues of hij matrix:'
!       do i=1,krydim
!         write(*,*) i,lambda(i)
!       enddo
      endif

!*******************************************************
! Step X+Y: calculate nroots residual vectors. Use
! Current dimension of Krylov space!
!*******************************************************

      call calc_resvec_r(io_bold,ladc,nroots,krydim, &
                         hijcal,hbv,lambda,res)

!*******************************************************
! Step X+Y: Enter loop over the roots. Form residual vector
! on the fly once it is needed. 
! If convergence is not reached determine a preconditioner
! and enlarge space if preconditioner possesses "enough
! new dimensionality".
!*******************************************************

      do kroot = 1,nroots

        resnorm=sqrt(dot_product(res(:,kroot),res(:,kroot)))

        write(listout,'(A3,I3,A3,I5,6X,F10.7,6X,F10.7)') '   ',  &
            kroot,'   ',krydim,lambda(kroot),resnorm


        if(resnorm.le.convthr) then

          converged(kroot) = .true.
          listout(47:48) =  ' *'
          listout(52:54) = '   '
          write(*,*) listout
          cycle

        else

          converged(kroot) = .false.
          listout(47:48) =  ' o'

          call calc_precon_r(ladc,adiag,lambda(kroot),res(:,kroot),prc)
          prenorm1=sqrt(dot_product(prc,prc))

          call calc_preort_r(io_bold,ladc,krydim,prc,pro)
          prenorm2=sqrt(dot_product(pro,pro))

          if( (prenorm2/prenorm1).ge.prethr) then

            pro = pro/prenorm2
            krydim = krydim + 1
            write(io_bold,rec=krydim) pro  !extend trial space
            listout(52:54) = ' XT'

          else

            listout(52:54) = ' NX'

          endif
          write(*,*) listout

        endif

      enddo ! kroot
      write(*,*)
!
!  check if all roots converged.
!
      allconverged = .true.
      do kroot=1,nroots
        if(converged(kroot).eqv..false.) then
          allconverged = .false.
        endif
      enddo
!
!  deallocate large HBV array and
!  check if dimensionality extension is still in allowed range.
!
      deallocate(hbv)

      if( (krydim + nroots).gt.maxdavsp) then
        write(*,*) 'Maximum Davidson space reached. Exiting'
        exit
      endif


      END DO   !Microiteration
!|
!|
!|
!|       end of microiteration loop
!|
!|______________________________________________________

      if(allconverged) then
        write(*,*) 'All roots converged.'
        call write_davevecs_r(io_bold,desrep,nroots,nmain,ladc, &
                              krydim,lambda,hijcal)
      else
        write(*,*) 'Nonconverged roots remain.'
      endif
!
!  cleanup
!  deallocate arrays
!  close and delete ADC matrix file
!  close and delete direct access file
!
      deallocate(hij)
      deallocate(hijcal)
      deallocate(lambda)
      deallocate(bv)
      deallocate(bw)
      deallocate(res)
      deallocate(prc)
      deallocate(pro)
      deallocate(work)
      deallocate(adiag)
      deallocate(converged)

      close(iobase,status='delete')
      close(io_bold,status='delete')

      write(*,*)
      write(*,'(70A1)') ('-',i=1,70)
      write(*,*) 'Davidson finished.'
      write(*,'(70A1)') ('-',i=1,70)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE init_krylov_r(io_bold,desrep,nroots,maxdavsp,ladc)

      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!   This routine checks the existence of some restart files and tries to fit
!   found restart vectors into the vector array.
!
!------------------------ calling variables -------------------------
!
      integer                    :: io_bold
      integer                    :: desrep,nroots,maxdavsp,ladc
!
!------------------------ local variables -------------------------
!
      logical                                :: isthere
      integer                                :: io
      integer                                :: i,j,k,ix,nsiz
      integer                                :: ida_iostat
      Character*11                           :: filename_d,filename_e
      Real*8                                 :: ortho,a1
      Real*8,allocatable,dimension(:)        :: qm,ec
      Real*8,allocatable,dimension(:,:)      :: bstart
      integer                                :: xnroots,xnmain,xladc
      integer                                :: unroots,uladc
      integer                                :: jdummy
      Real*8                                 :: rdummy
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
!
! check availability of DAVSTART.XX file
!
      IF(desrep.GT.9) THEN
        WRITE(filename_d,'(A8,A1,I2)') 'DAVSTART','.',desrep
      ELSE
        WRITE(filename_d,'(A8,A2,I1)') 'DAVSTART','.0',desrep
      ENDIF
      inquire(file=filename_d,exist=isthere)
      if(.not.isthere) stop 'Missing Davidson start vectors!'
!
!  read file
!
      io = get_file_unit()

      write(*,*) 'Reading start vectors from ',filename_d
      OPEN(io,file=filename_d,access='sequential', &
           status='unknown',form='unformatted')
      READ(io) nsiz
      WRITE(*,*) 'Found',nsiz,' Start vectors.'

      If(2*nroots.gt.maxdavsp) then
         Call Quit('Increase Davidson space!')
      Endif

      If(nroots.gt.nsiz) then
        Write(*,*) 'You request more states than there', &
                   ' are in the Main space.'
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
      WRITE(*,*) 'Orthogonality of start vectors:',ortho
!
! now make start vectors the initial Davidson vectors for the first 
! macro iteration
! 
      allocate(qm(ladc))
      Do i=1,nroots
        qm(1:nsiz) = bstart(1:nsiz,i)
        qm(nsiz+1:ladc) = 0.0d0
        write(io_bold,rec=i,iostat=ida_iostat) qm
        if(ida_iostat.ne.0) call &
           quit('Error writing direct access file for bold')
      Enddo

      deallocate(qm)
      deallocate(bstart)


      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine calc_resvec_r(iu,ladc,nroots,krydim, &
                               cji,hbv,lambda,res)
!
!  compute residual vectors res, filled within this routine.
!
      implicit none
!
! formal arguments
!
      integer                               :: iu
      integer                               :: ladc
      integer                               :: nroots
      integer                               :: krydim
      real*8,dimension(:,:)                 :: cji
      real*8,dimension(:,:)                 :: hbv
      real*8,dimension(:)                   :: lambda
      real*8,dimension(:,:)                 :: res
!
! local variables
!
      integer                               :: j,kroot
      real*8,parameter                      :: rzero = 0.0d0
      real*8,allocatable,dimension(:)       :: aux,bv
!
! executable code
! form sum_j (C_ji*(H*b_j - lambda_i*b_j))
!
      allocate(bv(ladc))
      allocate(aux(ladc))
!
! first loop over Krylov space (minimize disk reads!)
! Inner loop: over nroots
!
      res = rzero
      Do j=1,krydim
        read(iu,rec=j) bv
        Do kroot=1,nroots
          aux(:) = hbv(:,j) - lambda(kroot)*bv(:)
          aux(:) = aux(:)*cji(j,kroot)
          res(:,kroot) = res(:,kroot) + aux(:)
        Enddo  !kroot
      Enddo  !j

      deallocate(aux)
      deallocate(bv)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine calc_precon_r(ladc,adiag,lambdai,res,prc)
!
! compute preconditioner
!
      implicit none
!
! formal arguments
!
      integer                               :: ladc
      real*8,dimension(:)                   :: adiag
      real*8                                :: lambdai
      real*8,dimension(:)                   :: res
      real*8,dimension(:)                   :: prc
!
! local variables
!
      integer                               :: j
      real*8                                :: a1
!
! executable code
!

      Do j=1,ladc
        a1=lambdai - adiag(j)
        if(dabs(a1).lt.1.0e-03) then
          prc(j) = 100.0d0 * res(j)
        else
          prc(j)=res(j)/a1
        endif
      Enddo

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine calc_preort_r(iu,ladc,krydim,prc,pro)
!
! orthogonalize preconditioner to current subspace
! prc goes in, unchanged
! pro goes out, assigned
!
      implicit none
!
! formal arguments
!
      integer                               :: iu
      integer                               :: ladc
      integer                               :: krydim
      real*8,dimension(:)                   :: prc
      real*8,dimension(:)                   :: pro
!
! local variables
!
      integer                               :: j
      real*8                                :: a1
      real*8,allocatable,dimension(:)       :: aux
!
! executable code
!
      allocate(aux(ladc))


      pro = prc
      Do j = 1,krydim
        read(iu,rec=j) aux
        a1=dot_product(aux,prc)
        pro(:) = pro(:) - a1*aux(:)
      Enddo

      deallocate(aux)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE write_davevecs_r(io_bold,desrep,nroots,nmain,ladc, &
                                  krydim,lambda,cji)

      implicit none
!
!---------------Description--------------------------------------------
!
! converged eigenvectors (real) are written out.
!
!------------------------ calling variables -------------------------
!
      integer                    :: io_bold
      integer                    :: desrep,nroots,nmain,ladc,krydim
      real*8,dimension(:)        :: lambda
      real*8,dimension(:,:)      :: cji
!
!------------------------ local variables -------------------------
!
      integer                                :: k,j
      integer                                :: iu
      real*8                                 :: a1,enorm
      real*8,parameter                       :: rzero=(0.0d0)
      real*8,allocatable,dimension(:)        :: aux,bv
      character*11                           :: davevcname
!
!------------------------ local variables -------------------------
!
      interface

        Function get_file_unit()
          integer                      :: get_file_unit
        END Function

      end interface
!
!------------------------ executable code -------------------------
!
      allocate(aux(ladc))
      allocate(bv(ladc))

! *** open file

      iu = get_file_unit()
!     write(*,*) 'Unit for writing eigenvectors:',iu
      if(desrep.GT.9) THEN
        write(davevcname,'(A8,A1,I2)') 'ADCEVECS','.',desrep
      else
        write(davevcname,'(A8,A2,I1)') 'ADCEVECS','.0',desrep
      endif
      write(*,*) 'Writing out eigenvectors to file ',davevcname
      write(*,*)
      open(iu,file=davevcname,access='sequential', &
           status='unknown',form='unformatted')

! *** write header information

      write(iu) nroots
      write(iu) nmain
      write(iu) ladc

! *** loop over the roots

      Do k = 1,nroots

        write(iu) k
        write(iu) lambda(k)

! *** build kth eigenvector

        aux = rzero
        Do j=1,krydim
          read(io_bold,rec=j) bv
          aux = aux + cji(j,k)*bv
        Enddo

! *** normalize kth eigenvector

        enorm = sqrt(dot_product(aux,aux))
        write(*,*) 'Norm of',k,'th eigenvector:',enorm

! *** compute pole strength of kth eigenvector

        a1 = 0.0d0
        do j=1,nmain
          a1 = a1 + aux(j)*aux(j)
        enddo
        write(iu) a1
        write(iu) (aux(j),j=1,ladc)

      Enddo   !k
!
      close(iu)

      deallocate(aux)
      deallocate(bv)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
