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
