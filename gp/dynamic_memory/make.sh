rm *.a *.o *.mod

#xlf -O3 -g  -pg -c memory_tracker.f90
#gfortran  -g  -c memory_tracker.f90
gfortran  -g -fdefault-integer-8   -c memory_tracker.f90

#xlf -O3 -g  -pg -c memory_allocator.f90
gfortran -g -c -fdefault-integer-8 memory_allocator.f90

#xlf -O3 -g  -pg -c test_allocator.f90
#gfortran -g -c test_allocator.f90
gfortran -g -c -fdefault-integer-8  test_allocator.f90

#xlc -O3 -g  -DSYS_LINUX -pg -c -I/usr/include tracing.c
#xlc -O3 -g  -DSYS_LINUX -pg -c -I/usr/include tracker_util.c

#xlc -O3 -g  -pg -c -I/usr/include tracing.c
gcc -g -c -I/usr/include tracing.c
#xlc -O3 -g  -pg -c -I/usr/include tracker_util.c
gcc -g  -c -I/usr/include tracker_util.c

ar cr allocator.a memory_allocator.o memory_tracker.o tracing.o tracker_util.o

#xlf -O3 -g -pg -Wl,-E -I. -o test.x test_allocator.o allocator.a -lc
#gfortran -g  -I. -o test.x test_allocator.o allocator.a -lc
gfortran -g -fdefault-integer-8  -I. -o test.x test_allocator.o allocator.a -lc


