cmake_dependent_option(ENABLE_KRCC "KRCC module" ON
                       "ENABLE_UNRELEASED" OFF)

if(ENABLE_KRCC)
    add_definitions(-DMOD_KRCC)
    set(FIXED_FORTRAN_SOURCES
        ${FIXED_FORTRAN_SOURCES}
        src/krcc/cc_vec_fnc2.F
        src/krcc/optimization.F
        src/krcc/cc_vec_fnc.F
        src/krcc/op_t_occ_r.F
        src/krcc/integrals.F
        src/krcc/contraction_analyze.F
        src/krcc/len_for_tf_r.F
        src/krcc/dim_t1t2_to_t12_map_r.F
        src/krcc/krcc_inp.F
        src/krcc/krcc_exc_e.F
        src/krcc/contraction_performer.F
        src/krcc/diagonal.F
        src/krcc/krcc_memory.F
        src/krcc/program.F
        src/krcc/contraction_solver.F
        src/krcc/io.F
        src/krcc/analyze.F
        src/krcc/krcc_wrt.F
        src/krcc/dim_t1t2_to_t12_s_r.F
        src/krcc/interface.F
        src/krcc/ccinfo.F
        src/krcc/strings_cc.F
        src/krcc/fock.F
        src/krcc/h_exp_t_mem_r.F
        src/krcc/utils.F
        src/krcc/get_int.F
        src/krcc/inum_for_occ2_r.F
        )
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_KRCC)
endif()

message(STATUS "KRCC module: ${ENABLE_KRCC}")
