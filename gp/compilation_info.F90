subroutine print_compilation_info()
#include "priunit.h"
    write(lupri, '(a)') ' System                   | Linux-3.0.0-12-generic'
    write(lupri, '(a)') ' Processor                | x86_64'
    write(lupri, '(a)') ' Internal math            | OFF'
    write(lupri, '(a)') ' 64-bit integers          | OFF'
    write(lupri, '(a)') ' MPI                      | OFF'
    write(lupri, '(a)') ' XCFun                    | ON'
    write(lupri, '(a)') ' Fortran compiler         | /usr/bin/f95'
    write(lupri, '(a)') ' Fortran compiler version | GNU Fortran (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1'
    write(lupri, '(a)') ' Fortran flags            | -g -fcray-pointer -fbacktrace -DVAR_GFORTRAN -DVAR'
    write(lupri, '(a)') '                          | _MFDS -fno-range-check -O0 -g'
    write(lupri, '(a)') ' C compiler               | /usr/bin/gcc'
    write(lupri, '(a)') ' C compiler version       | gcc (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1'
    write(lupri, '(a)') ' C flags                  |   -O0 -g3'
    write(lupri, '(a)') ' BLAS                     | /usr/lib/libblas.so'
    write(lupri, '(a)') ' LAPACK                   | /usr/lib/liblapack.so'
    write(lupri, '(a)') ' Configuration time       | 2011-11-11 11:04:37.449961'
end subroutine
