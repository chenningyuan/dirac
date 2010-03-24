#This file probably requires GNU Make, but should be easily
#convertible to regular Make format.

# Compile time options:
# -DNDEBUG Turn off run time checks in production calculations to gain speed
# To compile without -lstdc++ needed use (g++) -fno-rtti -fno-exceptions 
# -fno-threadsafe-statics. Note the last option.

# C++ compiler and flags
CXX ?= g++
CXXFLAGS ?= -g -Wall -O #-DWITH_QD
CXXFLAGS+=-Iinclude -Isrc -Isrc/taylor -Isrc/functionals -Llib -fno-rtti -fno-exceptions -fno-threadsafe-statics
LIBS=

# F90 compiler and flags, used for example code
FC=gfortran
FFLAGS=-Wall -Llib -Jfortran

BASE_OBJ:=src/xcfun.o src/functional.o src/settings.o src/evaluators.o src/fortran.o
FUN_OBJ:=$(patsubst %.cpp,%.o,$(wildcard src/functionals/*.cpp))

testall:

lib: lib/libxcfun.so

lib/libxcfun.so: $(BASE_OBJ) $(FUN_OBJ)
	$(CXX) $(CXXFLAGS) -fPIC -shared -o lib/libxcfun.so $(LIBS) $^	


stripped: lib/libxc_fun.a
	strip --strip-unneeded $<

funeval: test/funeval.cpp lib/libxcfun.so
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS) -lxcfun

example: test/example.cpp lib/libxcfun.so
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS) -lxcfun

testall: test/testall.cpp lib/libxcfun.so
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS) -lxcfun

fortran/example.o: fortran/xcfun_module.o

fortran_example: fortran/example.o fortran/xcfun_module.o lib/libxcfun.so
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) -lstdc++

clean:
	rm -f `find . -name '*.o' -o -name '*~' -o -name '*.a'` fortran_example benchmark testxc


.SUFFIXES: .F90

.F90.o:
	$(FC) $(FFLAGS) -c -o $*.o $*.F90
