#This file probably requires GNU Make, but should be easily
#convertible to regular Make format.

# Compile time options:
# -DNDEBUG Turn off run time checks to gain speed
# -DNO_STDCXX Do not use any part of the C++ standard 
#             library, to allow linking without it.
#  XC_MAX_ORDER (in include/xc_fun.h) Maximal derivative order.

CXX ?= g++
CXXFLAGS ?= -Wall -O3 -DNDEBUG -DNO_STDCXX
CXXFLAGS+=-Iinclude -Isrc/taylor -Isrc/functionals -Llib -fno-rtti -fno-exceptions

lib/libxc_fun.a: src/xc_fun.o
	ar r $@ $< 

stripped: lib/libxc_fun.a
	strip --strip-unneeded $<

benchmark: test/benchmark.cpp lib/libxc_fun.a
	$(CXX) $(CXXFLAGS) $< -o $@ -lxc_fun

clean:
	rm -f `find . -name '*.o' -o -name '*~' -o -name '*.a'`
