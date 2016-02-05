cmake_dependent_option(ENABLE_ERI "ERI integral routines" ON
                       "ENABLE_UNRELEASED" OFF)

if(ENABLE_ERI)
    add_definitions(-DMOD_ERI)
    set(FIXED_FORTRAN_SOURCES
        ${FIXED_FORTRAN_SOURCES}
        src/eri/eri2car2.F
        src/eri/eri2odb.F
        src/eri/fck3mpi.F
        src/eri/eri2vec.F
        src/eri/fck3int.F
        src/eri/eri2her.F
        src/eri/eri2csi.F
        src/eri/eri2aob.F
        src/eri/eri2fck.F
        src/eri/eri2pro.F
        src/eri/eri2sph.F
        src/eri/eri2drv.F
        src/eri/eri2cft.F
        src/eri/eri2car1.F
        src/eri/fck3drv.F
        src/eri/fck3aob.F
        src/eri/fck3ftuv.F
        src/eri/eri2r12.F
        src/eri/eri2ave.F
        src/eri/eri2out.F
        src/eri/fck3odb.F
        src/eri/eri2ftuv.F
        src/eri/fck3box.F
        src/eri/eri2fck3.F
        src/eri/fck3mom.F
        src/eri/erisph.F
        src/eri/eri2dis.F
        src/eri/eri2sym.F
        src/eri/eri2exc.F
        )
    dirac_test(integrals_eri "")
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_ERI)
endif()

message(STATUS "ERI integral routines: ${ENABLE_ERI}")
