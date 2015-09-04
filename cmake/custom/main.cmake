include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/cfun
    )

include(Definitions)
include(Sources)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

# fix for Windows with 64bit integers on and build type release
# when using GNU (gfortran) compilers
# there is one source file for which -O3/-O2 is not good
if((CMAKE_SYSTEM_NAME MATCHES "Windows") AND (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    AND (ENABLE_64BIT_INTEGERS)
    AND (CMAKE_BUILD_TYPE MATCHES "release"))
    set_source_files_properties(
        ${PROJECT_SOURCE_DIR}/src/abacus/hermodules.F90
        PROPERTIES COMPILE_FLAGS "-O1"
        )
endif()

# FIXME re-introduce
# set_source_files_properties(${PROJECT_BINARY_DIR}/config_info.F90 PROPERTIES GENERATED 1)

set(DIRAC_SOURCES
    ${C_SOURCES}
    ${FREE_FORTRAN_SOURCES}
    ${FIXED_FORTRAN_SOURCES}
    )

add_library(
    dirac
    ${DIRAC_SOURCES}
    )

add_dependencies(dirac run_pamadm)
add_dependencies(dirac git_info)

if(ENABLE_XCFUN)
    add_dependencies(dirac xcfun_fortran_bindings)
endif()

if(ENABLE_PCMSOLVER)
    add_dependencies(dirac pcmsolver)
endif()

set(_list_of_executables)
macro(compile_standalone _executable _source _list_of_executables)
    add_executable(
        ${_executable}
        ${_source}
    )
    set_property(TARGET ${_executable} PROPERTY LINKER_LANGUAGE Fortran)
    set(_list_of_executables
        ${_list_of_executables}
        ${_executable}
    )
endmacro()
compile_standalone(dirac.x src/main/main.F90 "${_list_of_executables}")

option(ENABLE_UTILITIES "Enable utilities" ON)
if(ENABLE_UTILITIES)
    compile_standalone(dirac_mointegral_export.x utils/dirac_mointegral_export.F90 "${_list_of_executables}")
    compile_standalone(rspread.x utils/rspread.F90 "${_list_of_executables}")
    compile_standalone(labread.x utils/labread.F90 "${_list_of_executables}")
    compile_standalone(twofit.x utils/twofit.F90 "${_list_of_executables}")
    compile_standalone(cfread.x utils/cfread.F90 "${_list_of_executables}")
    compile_standalone(vibcal.x utils/vibcal.F90 "${_list_of_executables}")
    compile_standalone(polfit.x utils/polfit.F90 "${_list_of_executables}")
    compile_standalone(mx2fit.x utils/mx2fit.F90 "${_list_of_executables}")
    compile_standalone(diag.x utils/diag.F90 "${_list_of_executables}")
    compile_standalone(test_allocator.x utils/test_allocator.F90 "${_list_of_executables}")
endif()

option(ENABLE_RELCCSD_STANDALONE "Enable compilation of standalone relccsd" OFF)
if(ENABLE_RELCCSD_STANDALONE)
    compile_standalone(relccsd.x src/relccsd/ccmain.F "${_list_of_executables}")
endif()

option(ENABLE_HSCC "Enable compilation of the standalone hsmrcc.x" OFF)
if(ENABLE_HSCC) # note: this functionality depends on the pam flag
    compile_standalone(hsmrcc.x src/cc_external/hsmrcc/hsmrmain.c "${_list_of_executables}")
    set_property(TARGET hsmrcc.x PROPERTY LINKER_LANGUAGE C)
endif()

option(ENABLE_HSFS "Enable compilation of the standalone hsfscc.x" OFF)
if(ENABLE_HSFS) # note: this functionality depends on the pam flag
    compile_standalone(hsfscc.x src/cc_external/hsfscc/hsfs_main.F90 "${_list_of_executables}")
endif()

foreach(
    _executable
    ${_list_of_executables}
    )
    if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
        SET_TARGET_PROPERTIES(${_executable} PROPERTIES LINK_FLAGS "-Wl,-bbigtoc")
    endif()
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        SET_TARGET_PROPERTIES(${_executable} PROPERTIES LINK_FLAGS "-Wl,-E")
    endif()
    target_link_libraries(
        ${_executable}
        dirac
        ${EXTERNAL_LIBS}
        )

    if(ENABLE_STATIC_LINKING)
        if(ENABLE_OLD_LINKER)
           set(_no_dyn_export "")
        else()
           # Miro: the "-Wl,--no-export-dynamic" flag is only for newer
           # ld-linkers, fix for older linkers
           set(_no_dyn_export "-Wl,--no-export-dynamic")
        endif()
        target_link_libraries(${_executable} ${_no_dyn_export})
        unset(_no_dyn_export)
        if (BLAS_TYPE MATCHES Atlas)
            # activate extra linking flags for dirac.x only if this bug occurs, please do not delete
            # http://gcc.gnu.org/bugzilla/show_bug.cgi?id=30471#c7
            # this bug is: signal 11 (SIGSEGV):Segmentation fault, function __restore_rt (0x1613DA0),from file sigaction.c
            target_link_libraries(${_executable} -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
        endif()
        if (BLAS_TYPE MATCHES MKL AND CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            # miro: GNU+MKL-static needs extra lib
            target_link_libraries(${_executable} -ldl)
        endif()
    endif()
endforeach()

unset(_list_of_executables)
