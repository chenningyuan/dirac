option(ENABLE_PCMSOLVER "Enable Polarizable Continuum Model through PCMSolver" ON)

if(ENABLE_PCMSOLVER)
    find_package(ZLIB)
    if(NOT ZLIB_FOUND)
        set(ENABLE_PCMSOLVER OFF)
        message(STATUS "Polarizable Continuum Model via PCMSolver DISABLED")
        message(STATUS " PCMSolver dependencies NOT satisfied:")
        if(NOT ZLIB_FOUND)
           IF (ENABLE_STATIC_LINKING)
              message(STATUS "  - install static Zlib development libraries")
           else()
              message(STATUS "  - install dynamic Zlib development libraries")
           endif()
        endif()
    else()
        add_definitions(-DHAS_PCMSOLVER)
        message(STATUS "Polarizable Continuum Model via PCMSolver ENABLED")
        set(PCMSOLVER_TESTS OFF)
        get_filename_component(ZLIB_ROOT ${ZLIB_LIBRARIES} PATH)

        # PCMSolver does not know profile
        if(CMAKE_BUILD_TYPE MATCHES "profile")
            set(PCM_BUILD_TYPE "release")
        else()
            set(PCM_BUILD_TYPE ${CMAKE_BUILD_TYPE})
        endif()

        set(ExternalProjectCMakeArgs
            -DCMAKE_BUILD_TYPE=${PCM_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
            -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
            -DEXTRA_Fortran_FLAGS=${PCM_EXTRA_Fortran_FLAGS}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DEXTRA_C_FLAGS=${PCM_EXTRA_C_FLAGS}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DEXTRA_CXX_FLAGS=${PCM_EXTRA_CXX_FLAGS}
            -DENABLE_EIGEN_MKL=${ENABLE_EIGEN_MKL}
            -DBoost_NO_BOOST_CMAKE=ON
            -DBOOST_INCLUDEDIR=${BOOST_INCLUDEDIR}
            -DBOOST_LIBRARYDIR=${BOOST_LIBRARYDIR}
            -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
            -DENABLE_TESTS=${PCMSOLVER_TESTS}
            -DZLIB_ROOT=${ZLIB_ROOT}
            -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
            -DSUBMODULES_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
            -DENABLE_CXX11_SUPPORT=OFF
            -DENABLE_LOGGER=OFF
            -DENABLE_TIMER=OFF
            -DBUILD_STANDALONE=OFF
            -DENABLE_FORTRAN_API=ON
            )

            if(PCMSOLVER_TESTS)
                add_external(pcmsolver ${CMAKE_CTEST_COMMAND})
            else()
                add_external(pcmsolver)
            endif()

            link_directories(${PROJECT_BINARY_DIR}/external/lib)
            set(EXTERNAL_LIBS
                ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}
                pcm
                ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}
                getkw
                ${ZLIB_LIBRARIES}
                ${EXTERNAL_LIBS})
            if (ENABLE_STATIC_LINKING)
                list(REMOVE_ITEM EXTERNAL_LIBS "c")
                message(STATUS "PCM accompanying libs cleared of 'libc' due to static linking")
            endif()
    endif()
    set(PCMSOLVER_PARSE_DIR ${PROJECT_BINARY_DIR}/external/bin)
    include_directories(${PROJECT_BINARY_DIR}/external/pcmsolver-build/modules)
else()
    message(STATUS "Polarizable Continuum Model via PCMSolver DISABLED")
endif()
