if (ENABLE_64BIT_INTEGERS AND ${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND CMAKE_Fortran_COMPILER_ID MATCHES XL)
    SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS>")
    message(STATUS "For libinterest.a, objects mode set to 64 bit on IBM AIX with XL Fortran compiler")
endif()

add_library(
    interest
    src/git_revision_info.f90
    src/module_interest_eri.f90
    src/module_interest_hrr.f90
    src/module_interest_one.f90
    src/module_interest_osr.f90
    )
