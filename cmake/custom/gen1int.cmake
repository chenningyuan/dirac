option(ENABLE_GEN1INT "Enable Gen1Int library" ON)

if(ENABLE_GEN1INT)
    include_directories(${PROJECT_BINARY_DIR}/external/gen1int-build)

    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    )

    add_external(gen1int)

    set(EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libgen1int.a
        ${EXTERNAL_LIBS})
endif()

message(STATUS "Gen1Int module: ${ENABLE_GEN1INT}")
