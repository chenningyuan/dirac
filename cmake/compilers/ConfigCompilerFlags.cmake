include(SaveCompilerFlags)

if(CMAKE_C_COMPILER_WORKS)
    include(CFlags)
	# check for correct compiler version with Python script and CMake variable (if defined)
	include(FindCCompilerVersion)
endif()

if(CMAKE_CXX_COMPILER_WORKS)
    include(CXXFlags)
	# check for correct compiler version with Python script and CMake variable (if defined)
	include(FindCXXCompilerVersion)
endif()

if(CMAKE_Fortran_COMPILER_WORKS)
    include(FortranFlags)
	# check for correct compiler version with Python script and CMake variable (if defined)
	include(FindFortranCompilerVersion)
endif()
