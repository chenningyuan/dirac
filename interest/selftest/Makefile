f90:=ifort -opt-report
normalflags:=-c -O2
highflags:=-c -O3
lowflags:=-c -O0
nopflags:=-c -O0
default:=-c  -O0
LIB_PATH:=../lib
LIB:=interest

all: libinterest_selftest.x
	@echo "Interest selftest part is compiled!"

clean:
	@set -e; \
	echo -n "Cleaning Interest selftest part ... "; \
	rm -f *.mod *.o *.a *.x mak.log 2>&1 && echo "Sucessfully done!" \
        || (echo "Error occured!" && exit 1);

ROUTINES = module_interest_selftest.o main_interest_selftest.o
           

libinterest_selftest.x: ${ROUTINES}
	${f90} -o libinterest_selftest.x  $(ROUTINES) -L$(LIB_PATH) -l$(LIB)
#
# ROUTINES
#
main_interest_selftest.o: main_interest_selftest.f90
	${f90} ${highflags} main_interest_selftest.f90 
module_interest_selftest.o: module_interest_selftest.f90
	${f90} ${highflags} module_interest_selftest.f90 
