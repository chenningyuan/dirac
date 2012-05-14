Known bugs
==========

 This page contains patches to known bugs.


Bugs
====

-  CAMB3LYP parameter change requires recompile - UE
-  Expectation values fail for AOC HF - UE
-  Warnings "ERROR(READIN): Error in reading molecule.MOL." printed for
   no reason(?). - UE
-  [STRIKEOUT:IPARCTL is not passed correctly to Dirac when 64 bit
   integers are used - UE] Fixed by moving main to fortran 90.

DIRAC10
-------

Automatic augmentation for larger basis sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Automatic augmentation using s-aug, d-aug, t-aug, ... prefixes.

Only affects larger basis sets and could lead to a linearly dependent
basis. What happens is this: larger basis sets get split into blocks to
save memory. Automatic augmentation incorrectly augmented each block.
After this patch it only augments the most diffuse block as it should.

Save the following to a file:

::

    diff --git a/abacus/herbas.F b/abacus/herbas.F
    index 1c82a11..d1c0a9d 100644
    --- a/abacus/herbas.F
    +++ b/abacus/herbas.F
    @@ -253,6 +253,8 @@ C for backward compatibility the old notation is retained
                 BASTMP = BASNAM(1:ISTART - 1)
              END IF
           END IF
    +
    +      iaug_save = iaug
     C
     C     Determine if there are any user added basis functions/polarization
     C     functions
    @@ -368,6 +370,15 @@ C                INTEXP = INTEXP_SAVE
                      ENDIF
                    END IF

    +               if (split) then
    +!                 block has been split and this is not the last part
    +!                 augment only last part
    +!                 do not augment this part
    +                  iaug = 0
    +               else
    +                  iaug = iaug_save
    +               end if
    +
                    INTEXP_SAVE = INTEXP
                    INTORB_SAVE = INTORB
                    CALL READNU(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU,ALPHA,

then apply the patch in the DIRAC root directory:

::

     patch -p1 < file

finally, recompile DIRAC:

::

     make

DFT geometry optimization is incorrect
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A bug in the DFT geometry optimization was unfortunately released
together with DIRAC10. The bug leads to non-convergent geometry steps,
shows incorrectly integrated number of electrons and incorrect XC
energies at displaced geometries. The following patch restores geometry
optimization at the DFT level.

Save the following to a file:

::

    diff --git a/abacus/abaopt.F b/abacus/abaopt.F
    index b961e5d..0574a43 100644
    --- a/abacus/abaopt.F
    +++ b/abacus/abaopt.F
    @@ -2260,6 +2260,7 @@ C     (by running HERMIT and SIRIUS).
     C
           use dirac_cfg
           use dft_cfg
    +      use xc_driver

     #include "implicit.h"
     #include "priunit.h"
    @@ -2290,11 +2291,6 @@ C
           END IF
           NEWGEO = .TRUE.

    -      if (dirac_cfg_dft_calculation) then
    -!       get a new grid at new geometry
    -        call reset_grid()
    -      end if
    -
           IF (.NOT. EXSIR) THEN
              CALL EXESIR(WORK,LWORK,IWOFF,WRKDLM)
              IF (DOCCSD) THEN
    @@ -2356,6 +2352,11 @@ ckr#endif
              ENERGY = EMCSCF
           END IF
     #else /* PRG_DIRAC */
    +
    +      if (dirac_cfg_dft_calculation) then
    +!       get a new grid at new geometry
    +        call reset_grid()
    +      end if
     C
     C     Geometry optimization for Dirac
     C
    diff --git a/dft/basis_info.F90 b/dft/basis_info.F90
    index b32c86a..499cffb 100644
    --- a/dft/basis_info.F90
    +++ b/dft/basis_info.F90
    @@ -320,17 +320,21 @@ contains
           inquire(file  = file_name, exist = file_exists)

           if (file_exists) then
    -         close(file_unit, status = 'keep')
    -         return
    -      else
              open(file_unit,            &
                   file   = file_name,   &
    -              status = 'new',       &
    +              status = 'old',       &
                   form   = 'formatted', &
                   access = 'sequential')
    -         rewind(file_unit)
    +         close(file_unit, status = 'delete')
           end if

    +      open(file_unit,            &
    +           file   = file_name,   &
    +           status = 'new',       &
    +           form   = 'formatted', &
    +           access = 'sequential')
    +      rewind(file_unit)
    +
           only_large = .false.
           if (mc == 1) then
              only_large = .true.
    diff --git a/dirac/dirscf.F b/dirac/dirscf.F
    index 869d69f..4d0dc9d 100644
    --- a/dirac/dirscf.F
    +++ b/dirac/dirscf.F
    @@ -1041,6 +1041,10 @@ C
     C   Called from PSIDHF
     C
     C*****************************************************************************
    +
    +      use basis_info
    +      use occupation_info
    +
     #include "implicit.h"
     #include "priunit.h"
     #include "dummy.h"
    @@ -1068,6 +1072,11 @@ C       in ONEFCK ... no need to prepare LUFCK1 !!!!
           CALL WRIFCK(LUFCK1,FOCK,1)
           CLOSE(LUFCK1,STATUS='KEEP')

    +!     write info to files
    +!     this is then used in dft and visualization
    +      call basis_info_write()
    +      call occupation_info_write()
    +
           IF(ONESYS) THEN
     C
     C     Solve one-electron system ...

then apply the patch in the DIRAC root directory:

::

     patch -p1 < file

finally, recompile DIRAC:

::

     make

Parallel DFT crashes with MPICH2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Save the following to a file:

::

    Index: dft/xc_driver.F90
    ===================================================================
    --- dft/xc_driver.F90   (revision 12324)
    +++ dft/xc_driver.F90   (revision 12325)
    @@ -1675,8 +1675,9 @@
     #ifdef VAR_MPI
        subroutine dft_collect()
           integer :: ierr
    -    real(8) :: sum_r
    -    integer :: sum_i
    +      real(8) :: sum_r
    +      real(8) :: dummy
    +      integer :: sum_i

           call mpi_reduce(nr_electrons,         &
                           sum_r,                &
    @@ -1729,43 +1730,51 @@
           end if

           if (nr_fmat > 0) then
    +!        suggestion for Radovan (no dcopy + no extra alloc needed) - 
    +!        a simple test with this code piece works fine for me.
              if (i_am_master) then
    -            call alloc(fmat_container, mat_dim, mat_dim, nz, nr_fmat)
    +           call mpi_reduce(MPI_IN_PLACE,               &
    +                           fmat,                       &
    +                           mat_dim*mat_dim*nz*nr_fmat, &
    +                           mpi_double_precision,       &
    +                           mpi_sum,                    &
    +                           0,                          &
    +                           mpi_comm_world,             &
    +                           ierr)
    +         else
    +!          declare dummy as --> real(8) :: dummy
    +!          alternatively pass MPI_IN_PLACE instead of dummy for the recvbuf.
    +           call mpi_reduce(fmat,                       &
    +                           dummy,                      &
    +                           mat_dim*mat_dim*nz*nr_fmat, &
    +                           mpi_double_precision,       &
    +                           mpi_sum,                    &
    +                           0,                          &
    +                           mpi_comm_world,             &
    +                           ierr)
              end if
    -
    -         call mpi_reduce(fmat,                       &
    -                         fmat_container,             &
    -                         mat_dim*mat_dim*nz*nr_fmat, &
    -                         mpi_double_precision,       &
    -                         mpi_sum,                    &
    -                         0,                          &
    -                         mpi_comm_world,             &
    -                         ierr)
    -
    -         if (i_am_master) then
    -            call dcopy(mat_dim*mat_dim*nz*nr_fmat, fmat_container, 1, fmat, 1)
    -            call dealloc(fmat_container)
    -         end if
           end if

           if (do_geo) then
              if (i_am_master) then
    -            call alloc(property_gradient_container, nr_atoms*3)
    +           call mpi_reduce(MPI_IN_PLACE,               &
    +                           property_gradient,          &
    +                           nr_atoms*3,                 &
    +                           mpi_double_precision,       &
    +                           mpi_sum,                    &
    +                           0,                          &
    +                           mpi_comm_world,             &
    +                           ierr)
    +         else
    +           call mpi_reduce(property_gradient,          &
    +                           dummy,                      &
    +                           nr_atoms*3,                 &
    +                           mpi_double_precision,       &
    +                           mpi_sum,                    &
    +                           0,                          &
    +                           mpi_comm_world,             &
    +                           ierr)
              end if
    -
    -         call mpi_reduce(property_gradient,           &
    -                         property_gradient_container, &
    -                         nr_atoms*3,                   &
    -                         mpi_double_precision,         &
    -                         mpi_sum,                      &
    -                         0,                            &
    -                         mpi_comm_world,               &
    -                         ierr)
    -
    -         if (i_am_master) then
    -            call dcopy(nr_atoms*3, property_gradient_container, 1, property_gradient, 1)
    -            call dealloc(property_gradient_container)
    -         end if
           end if

        end subroutine

then apply the patch in the DIRAC root directory:

::

     patch -p0 < file

finally, recompile DIRAC:

::

     make

User-defined CAM functional line is cut after character 80
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Save the following to a file:

::

    Index: dirac/dirrdn.F
    ===================================================================
    --- dirac/dirrdn.F      (revision 12323)
    +++ dirac/dirrdn.F      (revision 12324)
    @@ -4019,7 +4019,7 @@
           LOGICAL  SYMSTP,NOINPUT,SETFUN,SUCCESS,IOTC,GAUNT_SAVE,
          &         TWOCOMP_SAVE, TWOCOMPBSS_SAVE
           CHARACTER PROMPT*1, WORD*7, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
    -     &          WORD1*7, LINE*80
    +     &          WORD1*7, LINE*80, LINE100*100
           DIMENSION KSOP(0:7),WORK(LWORK)
     C
           SAVE SETFUN
    @@ -4332,12 +4332,13 @@
     !          .DFT
       117      CONTINUE
     C&&& DFT    : Kohn-Sham (DFT) Hamiltonian
    -               READ(LUCMD,'(A80)') LINE
    +!radovan: line can be longer than 80
    +               read(lucmd, '(a100)') line100
                    IF (.NOT. SETFUN) THEN
     C              ... hjaaj: DFTINPUT accumulates def. of functional
     C                  so if we call more than once (as in geom.opt.)
     C                  we get n*functional for n'th call !!!! /3-Jun-2003 hjaaj
    -                  CALL DFTINPUT(LINE, INPERR)
    +                  CALL DFTINPUT(line100, INPERR)
                       if (inperr .gt. 0) then
                         call quit('unknown functional after .DFT')
                       end if

then apply the patch in the DIRAC root directory:

::

     patch -p0 < file

finally, recompile DIRAC:

::

     make

Parallel MCSCF/CI calculations on one node with more than 2 processes are wrong (multiple nodes calculations are fine)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

to apply the patch, please save the following patch to file "file"

::

    2188c2188
    <               IF( IIOMOD_REL .eq. 0 )THEN
    ---
    > !             IF( IIOMOD_REL .eq. 0 )THEN
    2199c2199
    <               END IF
    ---
    > !             END IF
     

and type then in the luciarel/ subdirectory of DIRAC:

::

     patch transform_r.F < file

finally, recompile DIRAC:

::

     make

DIRAC08
-------

xlf compiler
^^^^^^^^^^^^

Necessary preprocessor statements are missing. Apply the following patch
to the **configure** script if you use the xlf compiler:

::

    Index: /Dirac/trunk/configure
    ===================================================================
    --- Dirac/configure (revision 7775)
    +++ Dirac/configure (revision 7790)
    @@ -2194,4 +2194,5 @@
     # IBM's xlf should be comma delimited, the others take spaces
          IBMMOD=1
    +     cpp_f90="$cpp_f90 `getmods`"
          ;;
        *)

Parallel 64-bit integer code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The released parallel DIRAC08 code is set up to always link with a
32-bit integer MPI library (tested for both 32-bit integer code and for
64-bit integer code with Open MPI and with Intel MPI (MPICH2 derivate)).
However, if you link a **parallel 64-bit integer code** with a **64-bit
integer MPI library** you must apply the following patches to get a
functioning parallel code:

::

    Index: Dirac/trunk/include/infpar.h
    ===================================================================
    --- Dirac/include/infpar.h (rev 7816)
    +++ Dirac/include/infpar.h (rev 7817)
    @@ -56,8 +56,12 @@
    -      INTEGER*4 ISTAT
    +      INTEGER   ISTAT

    Index: Dirac/trunk/gp/global_file.F
    ===================================================================
    --- Dirac/gp/global_file.F (revision 7817)
    +++ Dirac/gp/global_file.F (revision 7818)
    @@ -313,7 +313,7 @@

           INCLUDE "mpif.h"
     #include "global_file.h"
    -      INTEGER*4 STATUS(MPI_STATUS_SIZE)
    +      INTEGER   STATUS(MPI_STATUS_SIZE)

Parallel direct MP2 calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you run a parallel direct MP2 calculation using DIRAC08 with the
.CNVINT flag present on the \*DHFCAL deck the calculation will fail.
However, if you restart your MP2 calculation from a DFCOEF file and
deactivate the .CNVINT flag the calculation will run fine. Example:

::

    **WAVE FUNCTIONS
    .DHF
    .MP2
    *DHFCAL
    .CLOSED SHELL
    112
    .INTFLG
    1 1 1
    .ERGCNV
    1.0E-8
    .CNVINT
    0.01 0.0001
    *MP2CAL
    .OCCUP
    48..56
    .VIRTUAL
    all
    .VIRTHR
    1000.0
    **END OF

The above given example will fail if you run it in parallel.

Temporary solution: Stop after the SCF calculation and restart MP2:

::

    **WAVE FUNCTIONS
    .DHF
    .MP2
    *DHFCAL
    .CLOSED SHELL
    112
    .INTFLG
    1 1 1
    .ERGCNV
    1.0E-8
    #.CNVINT
    #0.01 0.0001
    *MP2CAL
    .OCCUP
    48..56
    .VIRTUAL
    all
    .VIRTHR
    1000.0
    **END OF

