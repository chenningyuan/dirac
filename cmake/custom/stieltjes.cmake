option(ENABLE_STIELTJES "Enable stieltjes library" ON)

if(ENABLE_STIELTJES)
    message(STATUS "Stieltjes external module ENABLED")
    # we will be using f90 modules
 #  set(CMAKE_Fortran_MODULE_DIRECTORY
 #      ${PROJECT_BINARY_DIR}/external/stieltjes-build
 #      ${CMAKE_Fortran_MODULE_DIRECTORY}
 #      )
    include_directories (${PROJECT_BINARY_DIR}/external/stieltjes-build)

    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
    )

    add_external(stieltjes)

    set(EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libstieltjes.a
        ${EXTERNAL_LIBS})
else()
    message(STATUS "Stieltjes external module DISABLED")
endif()
