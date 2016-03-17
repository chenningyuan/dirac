include_directories(
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/davidson_lille/include/
    )

if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For librelccsd.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

set(FIXED_RELCCSD_FORTRAN_SOURCES
   ccdriv.F
   cceqn_diis.F
   cceqn_driver_amplitudes.F
   cceqn_driver_initialization.F
   cceqn_driver_lambda.F
   cceqn_mp2_amplitudes.F
   cceqn_t1_amplitudes.F
   cceqn_t1_utils.F
   cceqn_t2_amplitudes.F
   cceqn_t3_energy_correction.F
   cceqn_utils.F
   ccfspc.F
   ccgenp.F
   ccgetv.F
   ccgrad.F
   ccintm.F
   ccints.F
   ccmain.F
   ccpam.F
   ccrpa.F
   ccsort.F
   ccxblas.F
   contraction.F90
   intermediates.F90
   lambda_equation.F90
   relcc_cfg.F90
   waio.F
   wrapper_sorting.F90 
   ${PROJECT_SOURCE_DIR}/src/prp/mp2lag.F
   ${PROJECT_SOURCE_DIR}/src/main/dirac.F
# needs fix missing object pamgeo
   ${PROJECT_SOURCE_DIR}/src/dirac/dirgeo.F
)

set(FREE_RELCCSD_FORTRAN_SOURCES
    relcc_cfg.F90
    # needs dependency
    ${PROJECT_SOURCE_DIR}/src/main/mpi_framework.F90

    davidson_lille/mathlib/interfaces/blas95/blas95.F90
    davidson_lille/mathlib/interfaces/blas95/caxpy_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/cgemm_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/cgemv_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/cscal_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/daxpy_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/ddot_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/dgemm_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/dgemv_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/dnrm2_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/dscal_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/dznrm2_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/f95_precision.F90
    davidson_lille/mathlib/interfaces/blas95/saxpy_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/scnrm2_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/sdot_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/sgemm_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/sgemv_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/snrm2_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/sscal_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/zaxpy_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/zdotc_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/zgemm_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/zgemv_mkl95.F90
    davidson_lille/mathlib/interfaces/blas95/zscal_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/cgeev_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/cgesvd_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/cgetrf_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/cgetri_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/chpev_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/dgeev_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/dgesvd_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/dgetrf_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/dgetri_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/dlamc1.F90
    davidson_lille/mathlib/interfaces/lapack95/dlamc2.F90
    davidson_lille/mathlib/interfaces/lapack95/dlamc3.F90
    davidson_lille/mathlib/interfaces/lapack95/dlamc4.F90
    davidson_lille/mathlib/interfaces/lapack95/dlamc5.F90
    davidson_lille/mathlib/interfaces/lapack95/dlamch.F90
    davidson_lille/mathlib/interfaces/lapack95/dlansp.F90
    davidson_lille/mathlib/interfaces/lapack95/dlanst.F90
    davidson_lille/mathlib/interfaces/lapack95/dlapy2.F90
    davidson_lille/mathlib/interfaces/lapack95/dlarf.F90
    davidson_lille/mathlib/interfaces/lapack95/dlarfg.F90
    davidson_lille/mathlib/interfaces/lapack95/dlartg.F90
    davidson_lille/mathlib/interfaces/lapack95/dlascl.F90
    davidson_lille/mathlib/interfaces/lapack95/dlaset.F90
    davidson_lille/mathlib/interfaces/lapack95/dlasr.F90
    davidson_lille/mathlib/interfaces/lapack95/dlasrt.F90
    davidson_lille/mathlib/interfaces/lapack95/dlassq.F90
    davidson_lille/mathlib/interfaces/lapack95/dopgtr.F90
    davidson_lille/mathlib/interfaces/lapack95/dorg2l.F90
    davidson_lille/mathlib/interfaces/lapack95/dorg2r.F90
    davidson_lille/mathlib/interfaces/lapack95/dspev.F90
    davidson_lille/mathlib/interfaces/lapack95/dspev_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/dsptrd.F90
    davidson_lille/mathlib/interfaces/lapack95/dsteqr.F90
    davidson_lille/mathlib/interfaces/lapack95/dsterf.F90
    davidson_lille/mathlib/interfaces/lapack95/lapack95.F90
    davidson_lille/mathlib/interfaces/lapack95/sgeev_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/sgetrf_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/sgetri_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/sspev_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/xerbla.F90
    davidson_lille/mathlib/interfaces/lapack95/zgeev_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/zgesvd_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/zgetrf_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/zgetri_mkl95.F90
    davidson_lille/mathlib/interfaces/lapack95/zhpev_mkl95.F90
#   davidson_lille/include/cschmi_mod.F90
#   davidson_lille/include/dclchv_mod.F90
#   davidson_lille/include/dclvhv_mod.F90
#   davidson_lille/include/dschmi_mod.F90
#   davidson_lille/include/hred_mod.F90
#   davidson_lille/include/savvls_mod.F90
#   davidson_lille/include/savvso_mod.F90
#   davidson_lille/tools90/incord.F90
#   davidson_lille/tools90/lectur.F90
#   davidson_lille/tools90/macro9.F90
#   davidson_lille/tools90/mdavi.F90
#   davidson_lille/tools90/mdbg.F90
#   davidson_lille/tools90/mdterm.F90
#   davidson_lille/tools90/mfacto.F90
#   davidson_lille/tools90/mgenrl.F90
#   davidson_lille/tools90/mmulsy.F90
##  davidson_lille/tools90/mparam.F90
#   davidson_lille/tools90/mr3vec.F90
#   davidson_lille/tools90/mtrace.F90
#   davidson_lille/tools90/parmfil.F90
#   davidson_lille/tools90/rdmfil.F90
#   davidson_lille/tools90/timing.F90
    davidson_lille/include/fiflsh_mod.F90
    davidson_lille/tools90/caps.F90
    davidson_lille/tools90/jpsort.F90
    davidson_lille/tools90/kindvar.F90
    davidson_lille/tools90/matpri.F90
    davidson_lille/tools90/mconst.F90
    davidson_lille/tools90/mdebug.F90
    davidson_lille/tools90/mercod.F90
    davidson_lille/tools90/mphase.F90
    davidson_lille/david_schmidt.F90
    davidson_lille/dummy_sigma.F90
    davidson_lille/david_driver.F90
    davidson_lille/cdavid_driver.F90
    davidson_lille/calvhv.F90
    davidson_lille/david_trial_vectors_utils.F90
    davidson_lille/trf.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set_source_files_properties(${FREE_RELCCSD_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree=f90")
    set_source_files_properties(${FIXED_RELCCSD_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
endif()

#miro: all these manually picked up sources are necessary for building up the relccsd library, "make relccsd"
# - the quit.F file is not listed due to replacement by DIRAC QUIT routine
# - the ccunix.c is not included
add_library(
    relccsd
    ${FIXED_RELCCSD_FORTRAN_SOURCES}
    ${FREE_RELCCSD_FORTRAN_SOURCES}
)

# add intermodule dependencies
add_dependencies(relccsd gp)
add_dependencies(relccsd dirac)
add_dependencies(relccsd main)
add_dependencies(relccsd interface_mpi)
