if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qzerosize -qextname -qsuppress=cmpmsg -qtbtable")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g")
# miro: needed for linking C/C++ objects with XLF compiled sources
    add_definitions(-DVAR_XLF)
endif()
