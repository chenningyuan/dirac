include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/krcc
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libkrcc.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FREE_KRCC_FORTRAN_SOURCES
        krcc_inp.F90
)

set(FIXED_KRCC_FORTRAN_SOURCES
        cc_vec_fnc2.F
        optimization.F
        cc_vec_fnc.F
        op_t_occ_r.F
        integrals.F
        contraction_analyze.F
        len_for_tf_r.F
        dim_t1t2_to_t12_map_r.F
        krcc_exc_e.F
        contraction_performer.F
        diagonal.F
        krcc_memory.F
        program.F
        contraction_solver.F
        io.F
        analyze.F
        krcc_wrt.F
        dim_t1t2_to_t12_s_r.F
        interface.F
        ccinfo.F
        strings_cc.F
        fock.F
        h_exp_t_mem_r.F
        utils.F
        get_int.F
        inum_for_occ2_r.F
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_KRCC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_KRCC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()


if(ENABLE_RUNTIMECHECK)
    message(STATUS "runtime-check flags activated for the 'krcc' module without exceptions")
    set_source_files_properties(${FREE_KRCC_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
    set_source_files_properties(${FIXED_KRCC_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_runtimecheck})
endif()

add_library(
    krcc
    STATIC
    ${FREE_KRCC_FORTRAN_SOURCES}
    ${FIXED_KRCC_FORTRAN_SOURCES}
)

# dependencies...
add_dependencies(krcc gp)

