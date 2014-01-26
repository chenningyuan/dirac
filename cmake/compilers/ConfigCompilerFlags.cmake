include(SaveCompilerFlags)

if(CMAKE_C_COMPILER_WORKS)
    include(CFlags)
endif()

if(CMAKE_CXX_COMPILER_WORKS)
    include(CXXFlags)
endif()

if(CMAKE_Fortran_COMPILER_WORKS)
    include(FortranFlags)
endif()

if(CMAKE_C_COMPILER_WORKS AND CMAKE_CXX_COMPILER_WORKS AND CMAKE_Fortran_COMPILER_WORKS)
	# check for correct compiler version with Python script and CMake variable (if defined)
	include(FindCompilerVersion)
endif()