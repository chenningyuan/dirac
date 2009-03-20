module SymInfHSCC

   use IOmodule
   implicit none

   private ! ... specify variables as private by default

   integer, parameter :: maxirrep = 8
   character*4 point_group
   integer naxis, norder, numirrep
   integer irrep_product(maxirrep,maxirrep)
   integer totally_symmetric_irrep
   character*4 irrep_label(maxirrep)

!  public variables

   public :: maxirrep
   public :: point_group
   public :: naxis, norder, numirrep
   public :: irrep_product
   public :: totally_symmetric_irrep
   public :: irrep_label

 ! ... public routines ...

   public :: read_gamess_syminf
   public :: print_irrep_tables

!  private variables

   character*4 C1
   character*4 C2,  Cs,  Ci
   character*4 C2h, C2v, D2
   character*4 D2h

! ... nice symmetry specifications ... :)
   parameter (C1  = "C1 ")
   parameter (C2  = "C2 ")
   parameter (Cs  = "Cs ")
   parameter (Ci  = "Ci ")
   parameter (C2v = "C2v")
   parameter (C2h = "C2h")
   parameter (D2  = "D2 ")
   parameter (D2h = "D2h")

   character*4 A, B, AP, APP
   character*4 AG, AU
   character*4 A1, A2, B1, B2, B3
   character*4 BG, BU
   character*4 B1G, B1U
   character*4 B2G, B2U
   character*4 B3G, B3U

   parameter (A   = "A   ")
   parameter (B   = "B   ")
   parameter (AP  = "A'  ")
   parameter (APP = "A'' ")
   parameter (AG  = "AG  ")
   parameter (AU  = "AU  ")
   parameter (A1  = "A1  ")
   parameter (A2  = "A2  ")
   parameter (B1  = "B1  ")
   parameter (B2  = "B2  ")
   parameter (B3  = "B3  ")
   parameter (BG  = "BG  ")
   parameter (BU  = "BU  ")
   parameter (B1G = "B1G ")
   parameter (B1U = "B1U ")
   parameter (B2G = "B2G ")
   parameter (B2U = "B2U ")
   parameter (B3G = "B3G ")
   parameter (B3U = "B3U ")

   character*4 nolabel
   parameter (nolabel = "??? ")

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine read_gamess_syminf
!-----------------------------------------------------------------------
!
!     this routine reads point group symmetry information
!     from syminf file generated from gamess.
!
!     The data in syminf module are set in this call.
!
!     Called from modules/gmsmod.F90/InitGMSMod
!
!-----------------------------------------------------------------------

      implicit none
      character*5 group
      character*4 ptgrp
      integer nir
      integer ir, jr, irx, jrx, nrx
      integer kr, krxx(maxirrep)
      integer iabelian, irptype, iavail
      integer :: ierr=0 ! MI: important initialization for gfrortran 
                        ! which does not assign automatically zeroes 

      if (get_iprint() >= 5) then
       call Header('Output from "read_gamess_syminf"',-1)
      endif

      open(USYMINF,file=FileSYMINF)

      if (get_iprint() >= 5) then
        write(UOUT,'(2x,a,$)')  &
        'open file with Gamess symmetry information, FileSYMINF='
        write(UOUT,'(a)') FileSYMINF
      endif

!     read in point group

      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)

      read(USYMINF,9001) group, naxis, norder
 9001 format(3X,A5,3X,2X,1X,I2,2X,2X,1X,I3,1X)

      if (get_iprint() >= 5) then
         write(UOUT,'(2x,a,a)')  'group=',group
         write(UOUT,'(2x,a,i2)') 'naxis=',naxis
         write(UOUT,'(2x,a,i2)') 'norder=',norder
      endif

!    ... set up point group ...      
      ptgrp = " "

      if (group.eq."C1")                     ptgrp = C1
      if (group.eq."CS")                     ptgrp = Cs
      if (group.eq."CI")                     ptgrp = Ci
      if ((group.eq."CN" ).and.(naxis.eq.2)) ptgrp = C2
      if ((group.eq."S2N").and.(naxis.eq.1)) ptgrp = Ci
      if ((group.eq."CNH").and.(naxis.eq.2)) ptgrp = C2h
      if ((group.eq."CNV").and.(naxis.eq.2)) ptgrp = C2v
      if ((group.eq."DN" ).and.(naxis.eq.2)) ptgrp = D2
      if ((group.eq."DNH").and.(naxis.eq.2)) ptgrp = D2h

      if (ptgrp.eq." ") then
         write(UOUT,*)   &
         'read_gamess_syminf: unable to read point group information from syminf'
         call flushout
         call quit('unable to read point group information from syminf')
      endif

      if (get_iprint() >= 5) then
         write(UOUT,'(2x,a,a)') 'ptgrp=',ptgrp
      endif

      call setup_irrep_tables(ptgrp,ierr)

      if (ierr.ne.0) then
         write(UOUT,*) "read_gamess_syminf: error output :" 
         write(UOUT,*) "ptgrp=",ptgrp," ierr=",ierr
         write(UOUT,*) "trouble in setting up this point group !"
         call flushout
         call quit('read_gamess_syminf: trouble in setting up this point group !')
      endif

!     read type of group

      read(USYMINF,*)
      read(USYMINF,*) iabelian
      if (iabelian.ne.0) then
         write(UOUT,*) "syminf file contains non-abelian group!"
         call flushout
         call Quit('syminf file contains non-abelian group!')
      endif

      if (get_iprint() >= 5) then
        write(UOUT,'(2x,a,i3)') 'iabelian=',iabelian
      endif

!     read number of irreps

      read(USYMINF,*)
      read(USYMINF,*) nir
      if (nir.ne.numirrep) then
         write(UOUT,*) "syminf file contains invalid numirrep"
         call flushout
         call quit('syminf file contains invalid numirrep')
      endif

      if (get_iprint() >= 5) then
        write(UOUT,'(2x,a,2i3)') 'nir,numirrep:',nir,numirrep
      endif

!     read in type of irreps

      read(USYMINF,*)
      read(USYMINF,*) irptype
      if (irptype.ne.0) then
         write(UOUT,*) "syminf file contains a group with complex irreps"
         call flushout
         call quit('syminf file contains a group with complex irreps')
      endif

      if (get_iprint() >= 5) then
         write(UOUT,'(2x,a,i3)') 'irptype=',irptype
      endif

!     read in availability of irreps

      read(USYMINF,*)
      read(USYMINF,*) iavail
      if (iavail.ne.1) then
         write(UOUT,*)  &
         "syminf file contains a group whose irrep labels are not available !"
         call flushout
         call Quit('syminf file contains a group whose irrep labels are not available !')
      endif

      if (get_iprint() >= 5) then
         write(UOUT,'(2x,a,i3)') 'iavail=',iavail
      endif

!     skip ao/salc irrep information

      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      do ir = 1, nir
         read(USYMINF,*)
      enddo
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)

!     read in irrep multiplication table
!     and compare with the one we have

      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)
      read(USYMINF,*)

      do ir = 1, nir
      do jr = 1, nir

         read(USYMINF,*) irx, jrx, nrx

         if ((irx.ne.ir).or.(jrx.ne.jr).or.(nrx.ne.1 )) then
            write(UOUT,*) "Error reading irrep multiplication table from syminf file "
            call flushout
            call Quit('Error reading irrep multiplication table from syminf file')
         endif

         read(USYMINF,*) (krxx(kr),kr=1,nir)

         if (krxx(irrep_product(ir,jr)).ne.1) then
            write(UOUT,*) "irrep product table in syminf file does"
            write(UOUT,*) "not match the irrep table of the program"
            call flushout
            call Quit('irrep product table in syminf file does' &
                    //' not match the irrep table of the program')
         endif

      enddo
      enddo

!     skip further information in syminf file

      close(USYMINF)

      if (get_iprint() >= 5) then
        write(UOUT,'(/,2x,a)') '...control print-out of the multiplication table:'
        call print_irrep_tables
      endif

      return

      end subroutine read_gamess_syminf

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine setup_irrep_tables(group,ierr)
!-----------------------------------------------------------------------
!
!     this routine sets up irrep multiplication table for
!     abelian point groups C1, C2, Cs, Ci, C2h, C2v, D2, D2h
!
!     the order of irreps, irrep labels and irrep multiplication table
!     are the same as used in standard symmetry routines in gamess
!
!     called from "read_gamess_syminf"
!
!-----------------------------------------------------------------------

      implicit none
      character*(*) group
      integer ierr
      integer ir, jr

      numirrep = 0
      irrep_product = 0
      do ir = 1, maxirrep
        irrep_label(ir) = nolabel
      enddo

      if (group.eq.C1) then

         point_group = C1
         numirrep = 1
         irrep_product(1,1) = 1
         irrep_label(1) = A

      endif

      if ((group.eq.C2).or.(group.eq.Cs).or.(group.eq.Ci)) then

         numirrep = 2
         irrep_product(1,1) = 1
         irrep_product(2,1) = 2
         irrep_product(2,2) = 1

         if (group.eq.C2) then
            point_group = C2
            irrep_label(1) = A
            irrep_label(2) = B
         endif

         if (group.eq.Cs) then
            point_group = Cs
            irrep_label(1) = AP
            irrep_label(2) = APP
         endif

         if (group.eq.Ci) then
            point_group = Ci
            irrep_label(1) = AG
            irrep_label(2) = AU
         endif

      endif

      if ((group.eq.C2v).or.(group.eq.C2h).or.(group.eq.D2 )) then

         numirrep = 4

         irrep_product(1,1) = 1
         irrep_product(2,1) = 2
         irrep_product(2,2) = 1

         irrep_product(3,1) = 3
         irrep_product(4,1) = 4
         irrep_product(3,2) = 4
         irrep_product(4,2) = 3

         irrep_product(3,3) = 1
         irrep_product(4,3) = 2
         irrep_product(4,4) = 1

         if (group.eq.C2v) then
            point_group = C2v
            irrep_label(1) = A1
            irrep_label(2) = A2
            irrep_label(3) = B1
            irrep_label(4) = B2
         endif

         if (group.eq.C2h) then
            point_group = C2h
            irrep_label(1) = AG
            irrep_label(2) = AU
            irrep_label(3) = BU
            irrep_label(4) = BG
         endif

         if (group.eq.D2) then
            point_group = D2
            irrep_label(1) = A
            irrep_label(2) = B1
            irrep_label(3) = B3
            irrep_label(4) = B2
         endif

      endif

      if (group.eq.D2h) then

         numirrep = 8

         irrep_product(1,1) = 1
         irrep_product(2,1) = 2
         irrep_product(2,2) = 1

         irrep_product(3,1) = 3
         irrep_product(4,1) = 4
         irrep_product(3,2) = 4
         irrep_product(4,2) = 3

         irrep_product(3,3) = 1
         irrep_product(4,3) = 2
         irrep_product(4,4) = 1

         irrep_product(5,1) = 5
         irrep_product(6,1) = 6
         irrep_product(5,2) = 6
         irrep_product(6,2) = 5

         irrep_product(5,3) = 7
         irrep_product(6,3) = 8
         irrep_product(5,4) = 8
         irrep_product(6,4) = 7

         irrep_product(5,5) = 1
         irrep_product(6,5) = 2
         irrep_product(6,6) = 1

         irrep_product(7,1) = 7
         irrep_product(8,1) = 8
         irrep_product(7,2) = 8
         irrep_product(8,2) = 7
         
         irrep_product(7,3) = 5
         irrep_product(8,3) = 6
         irrep_product(7,4) = 6
         irrep_product(8,4) = 5

         irrep_product(7,5) = 3
         irrep_product(8,5) = 4
         irrep_product(7,6) = 4
         irrep_product(8,6) = 3

         irrep_product(7,7) = 1
         irrep_product(8,7) = 2
         irrep_product(8,8) = 1

         point_group = D2h

         irrep_label(1) = AG
         irrep_label(2) = AU
         irrep_label(3) = B3U
         irrep_label(4) = B3G
         irrep_label(5) = B1G
         irrep_label(6) = B1U
         irrep_label(7) = B2U
         irrep_label(8) = B2G

      endif

      if (numirrep.eq.0) then
        write(UOUT,*) 'setup_irrep_tables:  group=',group
        write(UOUT,*) 'setup_irrep_tables:  numirrep=',numirrep
        write(UOUT,*) 'setup_irrep_tables:  error exit !'
        call flushout
        ierr = 1
        return
      endif

!     specify totally symmetric irrep

      totally_symmetric_irrep = 1

!     symmetrize the irrep_product table

      do ir = 1, numirrep
      do jr = 1, (ir-1)
         irrep_product(jr,ir) = irrep_product(ir,jr)
      enddo
      enddo

      return

      end subroutine setup_irrep_tables

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine print_irrep_tables
!----------------------------------------------------------------------
!
!         print out irrep_product table from Gamess enviroment
!
!----------------------------------------------------------------------
      implicit none
      integer ir, jr

      write(UOUT,9001) point_group
      write(UOUT,9002) (irrep_label(ir),ir=1,numirrep)
      if (numirrep.eq.1) write(UOUT,9003)
      if (numirrep.eq.2) write(UOUT,9004)
      if (numirrep.eq.4) write(UOUT,9005)
      if (numirrep.eq.8) write(UOUT,9006)
      do ir = 1, numirrep
         write(UOUT,9007) irrep_label(ir),  &
         (irrep_label(irrep_product(ir,jr)),jr=1,numirrep)
         write(UOUT,9008)
      enddo
      write(UOUT,*)
  
      call flushout

      return

 9001 format("---------------------------------------",3(1H-),/,    &
             "IRREP MULTIPLICATION TABLE FOR GROUP = ",A4,/,        &
             "---------------------------------------",3(1H-),/)
 9002 format(6(1H ),'|',2(1H ),30(A4,2X))
 9003 format(6(1H-),'+',2(1H-), 6(1H-))
 9004 format(6(1H-),'+',2(1H-),12(1H-))
 9005 format(6(1H-),'+',2(1H-),24(1H-))
 9006 format(6(1H-),'+',2(1H-),48(1H-))
 9007 format(A4,2X, '|',2(1H ),10(A4,2X))
 9008 format(6(1H ),'|')

      end subroutine print_irrep_tables

end module SymInfHSCC
