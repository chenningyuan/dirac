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
CXXFLAGS ?= -g -Wall -O #-DWITH_QD
CXXFLAGS+=-Iinclude -Isrc -Isrc/taylor -Isrc/functionals -Llib -fno-rtti -fno-exceptions
LIBS= #-lqd

# F90 compiler and flags, used for example code
FC=gfortran
FFLAGS=-Wall -Llib -Jfortran # -DWITH_QD

FUN_OBJ:=$(patsubst %.cpp,%.o,$(wildcard src/functionals/*.cpp))


xcfun_test: src/xcfun.o src/functional.o src/settings.o src/evaluators.o $(FUN_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

lib: lib/libxcfun.so

lib/libxcfun.so: src/functional.o src/settings.o src/evaluators.o $(FUN_OBJ)
	$(CXX) $(CXXFLAGS) -fPIC -shared src/xcfun.cpp -o lib/libxcfun.so $(LIBS) $^ -DXCFUN_LIB	


stripped: lib/libxc_fun.a
	strip --strip-unneeded $<

funeval: test/funeval.cpp lib/libxcfun.a
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS) -lxcfun

example: test/example.cpp lib/libxcfun.so
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS) -lxcfun

fortran/example.o: fortran/xc_fun_module.o

fortran_example: fortran/example.o fortran/xc_fun_module.o lib/libxc_fun.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) -lstdc++

clean:
	rm -f `find . -name '*.o' -o -name '*~' -o -name '*.a'` fortran_example benchmark testxc


.SUFFIXES: .F90

.F90.o:
	$(FC) $(FFLAGS) -c -o $*.o $*.F90
