include(SaveCompilerFlags)

if(CMAKE_C_COMPILER_WORKS AND CMAKE_CXX_COMPILER_WORKS AND CMAKE_Fortran_COMPILER_WORKS)
    include(CFlags)
	include(CXXFlags)
	include(FortranFlags)
	
	# check for correct compiler version with Python script and CMake variable (if defined)
	# radovan: please no FATAL_ERRORS until this is debugged
	#          with Python 2.6 it gives me version 0.0 with Intel 12.0 and quits
	# include(FindCompilerVersion)
endif()
