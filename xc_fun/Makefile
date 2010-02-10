#This file probably requires GNU Make, but should be easily
#convertible to regular Make format.

# Compile time options:
# -DNDEBUG Turn off run time checks in production calculations to gain speed
# -DNO_STDCXX Do not use any part of the C++ standard 
#             library, to allow linking without it.
# -DWITH_QD Enable quad double precision functions. This requires
#           a _patched_ version of the QD library
# -DWITH_SINGLE Enable single precision interface. This doubles the
#               size of the binary, so be careful.

# C++ compiler and flags
CXX ?= g++
<<<<<<< .mine
CXXFLAGS ?= -g -Wall -O -DNDEBUG #-DWITH_QD
=======
CXXFLAGS ?= -g -Wall -O3 -DNDEBUG #-DWITH_QD
>>>>>>> .r10240
CXXFLAGS+=-Iinclude -Isrc/taylor -Isrc/functionals -Llib -fno-rtti -fno-exceptions
LIBS=-lxc_fun #-lqd

# F90 compiler and flags, used for example code
FC=gfortran
FFLAGS=-Wall -Llib -Jfortran # -DWITH_QD

lib/libxc_fun.a: src/xc_fun.o
	ar r $@ $<

src/xc_fun.o: $(wildcard src/functionals/*.h) src/taylor/taylor.h

stripped: lib/libxc_fun.a
	strip --strip-unneeded $<

benchmark: test/benchmark.cpp lib/libxc_fun.a
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS)

funtest: test/funtest.cpp lib/libxc_fun.a
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS)

funeval: test/funeval.cpp lib/libxc_fun.a
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS)

fortran/example.o: fortran/xc_fun_module.o

fortran_example: fortran/example.o fortran/xc_fun_module.o lib/libxc_fun.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) -lstdc++

clean:
	rm -f `find . -name '*.o' -o -name '*~' -o -name '*.a'` fortran_example benchmark testxc


.SUFFIXES: .F90

.F90.o:
	$(FC) $(FFLAGS) -c -o $*.o $*.F90
