message("Report of configuration settings:")
message("CMAKE_SYSTEM_NAME          ${CMAKE_SYSTEM_NAME}          ")
message("ENABLE_INTERNAL_MATH      ${ENABLE_INTERNAL_MATH}      ")
message("ENABLE_64BIT_INTEGERS      ${ENABLE_64BIT_INTEGERS}      ")
message("ENABLE_MPI                 ${ENABLE_MPI}                 ")
message("ENABLE_XCFUN               ${ENABLE_XCFUN}               ")
message("ENABLE_STATIC_LINKING      ${ENABLE_STATIC_LINKING}      ")
message("CMAKE_BUILD_TYPE           ${CMAKE_BUILD_TYPE}           ")
message("CMAKE_Fortran_COMPILER     ${CMAKE_Fortran_COMPILER}     ")
message("FORTRAN_COMPILER_VERSION   ${FORTRAN_COMPILER_VERSION}   ")
message("CMAKE_Fortran_FLAGS        ${CMAKE_Fortran_FLAGS}        ")
message("CMAKE_Fortran_FLAGS_DEBUG  ${CMAKE_Fortran_FLAGS_DEBUG}  ")
message("CMAKE_Fortran_FLAGS_RELEAS ${CMAKE_Fortran_FLAGS_RELEASE}")
message("CMAKE_C_COMPILER           ${CMAKE_C_COMPILER}           ")
message("C_COMPILER_VERSION         ${C_COMPILER_VERSION}         ")
message("CMAKE_C_FLAGS              ${CMAKE_C_FLAGS}              ")
message("CMAKE_C_FLAGS_DEBUG        ${CMAKE_C_FLAGS_DEBUG}        ")
message("CMAKE_C_FLAGS_RELEASE      ${CMAKE_C_FLAGS_RELEASE}      ")
if (BLAS_LIBRARIES)
    message("BLAS_LIBRARIES             ${BLAS_LIBRARIES}             ")
endif()
if (LAPACK_LIBRARIES)
    message("LAPACK_LIBRARIES           ${LAPACK_LIBRARIES}           ")
endif()
message("CMAKE_Fortran_LINK_EXECUTABLE  ${CMAKE_Fortran_LINK_EXECUTABLE} ")
if (CMAKE_Fortran_LINK_FLAGS)
    message("CMAKE_Fortran_LINK_FLAGS   ${CMAKE_Fortran_LINK_FLAGS} ")
endif()
if (LINK_FLAGS)
    message("LINK_FLAGS                 ${LINK_FLAGS} ")
endif()
if (FLAGS)
    message("FLAGS                      ${FLAGS} ")
endif()
if (LINK_LIBRARIES)
    message("LINK_LIBRARIES             ${LINK_LIBRARIES} ")
endif()
get_directory_property( DirDefs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
foreach( d ${DirDefs} )
    message( STATUS "Found Define: " ${d} )
endforeach()

message("GIT_REVISION               ${GIT_REVISION}               ")

## check integer setting of OpenMPI Fortran compiler, mpif90
## TODO: do the real check, not only print warning message !
if(ENABLE_MPI)
     message("\nPlease check whether your MPI library library ")
    if(ENABLE_64BIT_INTEGERS)
        message("is configured to work with Fortran 64-bit integers (integer size 8):")
    else()
        message("is configured to work with Fortran 32-bit integers (integer size 4):")
    endif()
    message("example for OpenMPI: <OMPI path>/ompi_info -a | grep 'Fort integer size'")
# python script to perform the OpenMPI check of integer
message("results is:")
EXECUTE_PROCESS(COMMAND python -c "
import sys
import os
import string
mpi_path=os.path.dirname('${CMAKE_Fortran_COMPILER}'); 
ompi_info=os.path.join(mpi_path,'ompi_info'); 
if os.path.isfile(ompi_info):
    command=ompi_info+' -a | grep \"Fort integer size\" '
    print command
    sys.stdout.flush()
    os.system(command)
    sys.exit(0)
else:
    print >> stderr,'can not determine ompi_info !'
    sys.exit(-1)
"
RESULT_VARIABLE ret_var
)
#message("ret_var=${ret_var}")
endif()

