if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # switch off warnings due to incompatibility XCode 4 and Intel 11 on OsX 10.6
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Qoption,ld,-w")
    endif()
endif()

# workaround for strange problem with cilkrts lib
# intel does not ship the static version of cilkrts
# and for some reason it gets detected and included in these variables
# breaking the static linking buildup with intel compilers
foreach(_lang C CXX)
    if(CMAKE_${_lang}_IMPLICIT_LINK_LIBRARIES)
        list(REMOVE_ITEM CMAKE_${_lang}_IMPLICIT_LINK_LIBRARIES "cilkrts")
    endif()
endforeach()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    # remove -rdynamic flag offensive for PGI Fortran in static linking
    if (CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
        list(REMOVE_ITEM CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-rdynamic")
        message(STATUS "pgf90 flag -rdynamic fremoved from CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS")
    endif()
endif()

if ( ENABLE_STATIC_LINKING AND BLAS_TYPE MATCHES OPENBLAS AND (CMAKE_Fortran_COMPILER_ID MATCHES Intel OR CMAKE_Fortran_COMPILER_ID MATCHES GNU) AND (MPI_FOUND AND ENABLE_MPI) )
    #miro: workaround to get static linking with OpenMPI/Intel,GNU/OPENBLAS
    list(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES  "pthread")
    list(REMOVE_ITEM CMAKE_C_IMPLICIT_LINK_LIBRARIES  "pthread")
    list(REMOVE_ITEM CMAKE_CXX_IMPLICIT_LINK_LIBRARIES  "pthread")
    message(STATUS "Linking flag -lpthread removed from CMAKE_SHARED_LIBRARY_LINK_Fortran/CXX/C_FLAGS (OpenMPI/Intel,GNU/OPENBLAS/static)")
endif()

if ( ENABLE_STATIC_LINKING AND BLAS_TYPE MATCHES OPENBLAS AND (CMAKE_Fortran_COMPILER_ID MATCHES Intel) AND (MPI_FOUND AND ENABLE_MPI) )
    #miro: workaround to get dirac.x for OpenMPI/Intel/OPENBLAS/static linking
    list(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES  "dl")
    list(REMOVE_ITEM CMAKE_C_IMPLICIT_LINK_LIBRARIES  "dl")
    list(REMOVE_ITEM CMAKE_CXX_IMPLICIT_LINK_LIBRARIES  "dl")
    message(STATUS "Linking flag -ldl removed from CMAKE_SHARED_LIBRARY_LINK_Fortran/CXX/C_FLAGS (OpenMPI/Intel/OPENBLAS/static)")
endif()


