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

# add stdc++ library to satisfy pgi compilers
if(CMAKE_CXX_COMPILER_ID MATCHES PGI)
    list(APPEND CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++")
endif()
