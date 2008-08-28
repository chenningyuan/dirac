      subroutine scatter_vector(length,
     &                          orbital_rotation_vector,
     &                          h,
     &                          response_vector,
     &                          matrix,
     &                          pg_sym)

      implicit none

#include "priunit.h"
#include "dgroup.h"
#include "dcborb.h"

      integer          length
      integer          orbital_rotation_vector
      double precision h
      double precision response_vector
      double precision matrix
      integer          pg_sym

      integer          x, i, s, iz
      double precision f

      dimension        orbital_rotation_vector(2, length)
      dimension        response_vector(length, nz)
      dimension        matrix(norbt, norbt, nz)

      do x = 1, length

        i = orbital_rotation_vector(1, x)
        s = orbital_rotation_vector(2, x)

        do iz = 1, nz

          if (ipqtoq(iz, pg_sym-1) > 1) then
            f = -1.0
          else
            f =  1.0
          end if

!                row
!                ¦
!                ¦  column
!                ¦  ¦
          matrix(s, i, iz) = matrix(s, i, iz)
     &                     +     response_vector(x, iz)
          matrix(i, s, iz) = matrix(i, s, iz)
     &                     - f*h*response_vector(x, iz)
        end do
      end do

      end

      subroutine fill_work(do_pp,
     &                     do_pn,
     &                     length_pp,
     &                     length_pn,
     &                     orbital_rotation_vector_pp,
     &                     orbital_rotation_vector_pn,
     &                     jxop_pp,
     &                     jxop_pn)

      implicit none

      logical   do_pp
      logical   do_pn
      integer   length_pp
      integer   length_pn
      integer   orbital_rotation_vector_pp
      integer   orbital_rotation_vector_pn
      integer   jxop_pp
      integer   jxop_pn

      dimension orbital_rotation_vector_pp(2, *)
      dimension orbital_rotation_vector_pn(2, *)
      dimension jxop_pp(2, *)
      dimension jxop_pn(2, *)

      integer   i

      if (do_pp) then
        do i = 1, length_pp
          jxop_pp(1, i) = orbital_rotation_vector_pp(1, i)
          jxop_pp(2, i) = orbital_rotation_vector_pp(2, i)
        end do
      end if

      if (do_pn) then
        do i = 1, length_pn
          jxop_pn(1, i) = orbital_rotation_vector_pn(1, i)
          jxop_pn(2, i) = orbital_rotation_vector_pn(2, i)
        end do
      end if

      end

      subroutine general_density_block(density,
     &                                 block,
     &                                 get_density,
     &                                 irepop,
     &                                 gao_slice_left,
     &                                 gao_slice_right,
     &                                 buffer,
     &                                 component_left,
     &                                 component_right,
     &                                 scale_by)

      implicit none

      double precision density
      double precision block
      logical          get_density
      integer          irepop
      double precision gao_slice_left
      double precision gao_slice_right
      double precision buffer
      integer          component_left
      integer          component_right
      double precision scale_by

#include "dgroup.h"
#include "dcbham.h"
#include "dcbbas.h"

      dimension        block(ntbas(0), ntbas(0))
      dimension        gao_slice_left(ntbas(0))
      dimension        gao_slice_right(ntbas(0))
      dimension        buffer(ntbas(0))

      integer          i1, i2, n1, n2, j, irep1, irep2
      double precision f, ddot

      if (.not.((component_left + component_right > 2)
     &    .and. (mc == 1))) then

        do irep1 = 0, nbsym - 1
          n1 = nbbas(irep1, component_left)
          if (n1 > 0) then

            irep2 = ieor(irep1, irepop)
            n2 = nbbas(irep2, component_right)
            if (n2 > 0) then

              if (get_density) then
                i1 = ibbas(irep1, component_left)  + 1
                i2 = ibbas(irep2, component_right) + 1

                call dgemv('N',
     &                     nbbas(irep1, component_left),
     &                     nbbas(irep2, component_right),
     &                     1.0d0,
     &                     block(i1, i2),
     &                     ntbas(0),
     &                     gao_slice_right(i2),
     &                     1,
     &                     0.0d0,
     &                     buffer,
     &                     1)

                density = density
     &                  + scale_by*ddot(nbbas(irep1, component_left),
     &                                  gao_slice_left(i1),
     &                                  1,
     &                                  buffer,
     &                                  1)

              else
                i1 = ibbas(irep1, component_left)
                i2 = ibbas(irep2, component_right)

                do j = 1, n2
                  f = scale_by*density*gao_slice_right(i2 + j)
                  call daxpy(n1,
     &                       f,
     &                       gao_slice_left(i1 + 1),
     &                       1,
     &                       block(i1 + 1, i2 + j),
     &                       1)
                end do
              end if
            end if
          end if
        end do
      end if

      end

      subroutine aoresponse_input(word)

      use aoresponse_input_m

      implicit none

#include "priunit.h"

      integer, parameter :: ntable = 60
      integer            :: i
      character(1)       :: prompt
      character(7)       :: word, word_save
      character(7)       :: table(ntable)

!                   1         2         3         4
      data table /'.SKIPEE','.SKIPEP','.THRESH','.ALPHA ',
!                   5         6         7         8
     &            '.BETA  ','.GAMMA ','.FREQ  ','.BETA2N',
!                   9        10        11        12
     &            '.DC-SHG','.THG   ','.IDRI  ','.DC-KER',
!                  13        14        15        16
     &            '.FREQ_P','.0GAMMA','.0BETA ','.INTEGR',
!                  17        18        19        20
     &            '.VISUAL','.TESTXC','.xxxxxx','.xxxxxx',
!                  21        22        23        24
     &            '.xxxxxx','.xxxxxx','.xxxxxx','.xxxxxx',
!                  25        26        27        28
     &            '.xxxxxx','.xxxxxx','.xxxxxx','.xxxxxx',
!                  29        30        31        32
     &            '.xxxxxx','.xxxxxx','.xxxxxx','.xxxxxx',
!                  33        34        35        36
     &            '.xxxxxx','.xxxxxx','.xxxxxx','.xxxxxx',
!                  37        38        39        40
     &            '.xxxxxx','.xxxxxx','.xxxxxx','.xxxxxx',
!                  41        42        43        44
     &            '.xxxxxx','.xxxxxx','.xxxxxx','.xxxxxx',
!                  45        46        47        48
     &            '.xxxxxx','.xxxxxx','.xxxxxx','.xxxxxx',
!                  49        50        51        52
     &            '.xxxxxx','.xxxxxx','.xxxxxx','.xxxxxx',
!                  53        54        55        56
     &            '.xxxxxx','.xxxxxx','.xxxxxx','.xxxxxx',
!                  57        58        59        60
     &            '.xxxxxx','.xxxxxx','.xxxxxx','.xxxxxx'/


!     init
!     ====

      do i = 1, max_nr_op
        perturbation_op(i)%op_name = ''
      end do

      call add_op('XDIPLEN', 1, 1, (/1.0d0/), (/'XDIPLEN '/))
      call add_op('YDIPLEN', 1, 1, (/1.0d0/), (/'YDIPLEN '/))
      call add_op('ZDIPLEN', 1, 1, (/1.0d0/), (/'ZDIPLEN '/))

      aoresponse_skip_pp            = .false.
      aoresponse_skip_pn            = .false.

      aoresponse_threshold_response = 1.0d-8
      aoresponse_threshold_rhs      = 1.0d-8
      aoresponse_freq               = 0.0d0
      aoresponse_freq_process       = 0.0d0

      aoresponse_alpha              = .false.
      aoresponse_beta               = .false.
      aoresponse_beta_static        = .false.
      aoresponse_beta_2np1          = .false.
      aoresponse_gamma              = .false.

      aoresponse_gamma_static       = .false.
      aoresponse_dc_shg             = .false.
      aoresponse_thg                = .false.
      aoresponse_idri               = .false.
      aoresponse_dc_kerr            = .false.

      aoresponse_integrate          = .false.
      aoresponse_visual             = .false.
      aoresponse_test_xc            = .false.


!     process input
!     =============

      if (word .eq. '*AORESP') then
        word_save = word
  100   continue
        read(lucmd, '(a7)') word
        prompt = word(1:1)
        if (prompt .eq. '!' .or. prompt .eq. '#') then
          go to 100
        else if (prompt .eq. '.') then
          do 200 i = 1, ntable
            if (table(i) .eq. word) then
            go to ( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     &             11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     &             21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
     &             31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
     &             41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
     &             51, 52, 53, 54, 55, 56, 57, 58, 59, 60), i
            end if
  200     continue
          if (word .eq. '.OPTION') then
            call prtab(ntable,
     &                 table,
     &                 word_save//' input keywords',
     &                 lupri)
            go to 100
          end if
          write(lupri, '(/,3a,/)') ' keyword "',
     &                             word,
     &                             '" not recognized in *AORESPONSE'
          call prtab(ntable,
     &               table,
     &               word_save//' input keywords',
     &               lupri)
          call quit('illegal keyword in *AORESPONSE')

    1     continue
          aoresponse_skip_pp = .true.
          go to 100

    2     continue
          aoresponse_skip_pn = .true.
          go to 100

    3     continue
          read(lucmd, *) aoresponse_threshold_response,
     &                   aoresponse_threshold_rhs
          go to 100

    4     continue
          aoresponse_alpha = .true.
          go to 100

    5     continue
          aoresponse_beta = .true.
          go to 100

    6     continue
          aoresponse_gamma = .true.
          go to 100

    7     continue
          read(lucmd, *) aoresponse_freq(1)
          read(lucmd, *) aoresponse_freq(2)
          read(lucmd, *) aoresponse_freq(3)
          read(lucmd, *) aoresponse_freq(4)
          go to 100

    8     continue
          aoresponse_beta_2np1 = .true.
          go to 100

    9     continue
          aoresponse_dc_shg  = .true.
          go to 100

   10     continue
          aoresponse_thg = .true.
          go to 100

   11     continue
          aoresponse_idri = .true.
          go to 100

   12     continue
          aoresponse_dc_kerr = .true.
          go to 100

   13     continue
          read(lucmd, *) aoresponse_freq_process
          go to 100

   14     continue
          aoresponse_gamma_static = .true.
          go to 100

   15     continue
          aoresponse_beta_static = .true.
          go to 100

   16     continue
          aoresponse_integrate = .true.
          go to 100

   17     continue
          aoresponse_visual = .true.
          go to 100

   18     continue
          aoresponse_test_xc = .true.
          go to 100

   19     continue
   20     continue
   21     continue
   22     continue
   23     continue
   24     continue
   25     continue
   26     continue
   27     continue
   28     continue
   29     continue
   30     continue
   31     continue
   32     continue
   33     continue
   34     continue
   35     continue
   36     continue
   37     continue
   38     continue
   39     continue
   40     continue
   41     continue
   42     continue
   43     continue
   44     continue
   45     continue
   46     continue
   47     continue
   48     continue
   49     continue
   50     continue
   51     continue
   52     continue
   53     continue
   54     continue
   55     continue
   56     continue
   57     continue
   58     continue
   59     continue
   60     continue
          go to 100

        else if (prompt .eq. '*') then
          go to 300
        else
          write(lupri, '(/,3a,/)') ' keyword "',
     &                             word,
     &                             '" not recognized in *AORESPONSE'
          call prtab(ntable,
     &               table,
     &               word_save//' input keywords',
     &               lupri)
          call quit('illegal keyword in *AORESPONSE')
        end if
      end if
  300 continue

      end
