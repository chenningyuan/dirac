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
        set(PCMSOLVER_TESTS ON)
        get_filename_component(ZLIB_ROOT ${ZLIB_LIBRARIES} PATH)

        if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
            set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/external)
            set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
            set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
        else()
            set(EXTERNAL_PROJECT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/external)
            set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE INTERNAL "")
        endif()

        # PCMSolver does not know profile
        if(CMAKE_BUILD_TYPE MATCHES "profile")
            set(PCM_BUILD_TYPE "release")
        else()
            set(PCM_BUILD_TYPE ${CMAKE_BUILD_TYPE})
        endif()

        set(ExternalProjectCMakeArgs
            # PCMSolver does not know profile
            if(CMAKE_BUILD_TYPE MATCHES "profile")
                set(CMAKE_BUILD_TYPE "release")
                endif
                -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PROJECT_INSTALL_PREFIX} #${PROJECT_BINARY_DIR}/external
                -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                -DEXTRA_Fortran_FLAGS=${PCM_EXTRA_Fortran_FLAGS}
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DEXTRA_C_FLAGS=${PCM_EXTRA_C_FLAGS}
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DEXTRA_CXX_FLAGS=${PCM_EXTRA_CXX_FLAGS}
                -DENABLE_CXX11_SUPPORT=OFF
                -DPARENT_INCLUDE_DIR=${PROJECT_BINARY_DIR}/external/include
                -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
                -DEIGEN3_ROOT=${EIGEN3_ROOT}
                -DDISABLE_EIGEN_OWN=${DISABLE_EIGEN_OWN}
                -DENABLE_EIGEN_MKL=${ENABLE_EIGEN_MKL}
                -DPARENT_BINARY_DIR=${PROJECT_BINARY_DIR}
                -DBOOST_INCLUDEDIR=${BOOST_INCLUDEDIR}
                -DBOOST_LIBRARYDIR=${BOOST_LIBRARYDIR}
                -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
                -DENABLE_TESTS=${PCMSOLVER_TESTS}
                -DZLIB_ROOT=${ZLIB_ROOT}
                -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
                -DENABLE_LOGGER=OFF
                -DENABLE_TIMER=OFF
                )

            if(PCMSOLVER_TESTS)
                add_external(pcmsolver ${CMAKE_CTEST_COMMAND})
            else()
                add_external(pcmsolver)
            endif()

            link_directories(${EXTERNAL_PROJECT_INSTALL_PREFIX}/lib
                ${EXTERNAL_PROJECT_INSTALL_PREFIX})
            set(EXTERNAL_LIBS
                stdc++
                pcm
                stdc++
                getkw
                ${ZLIB_LIBRARIES}
                ${EXTERNAL_LIBS})
    endif()
    set(PCMSOLVER_PARSE_DIR ${EXTERNAL_PROJECT_INSTALL_PREFIX}/bin)
else()
    message(STATUS "Polarizable Continuum Model via PCMSolver DISABLED")
endif()
