#
# Check the of integer size compatibility of math libraries (BLAS).
#
# These DIRAC CMake internal variables are needed:
#
#  MATH_LIBS , EXPLICIT_LIBS,  ENABLE_BUILTIN_BLAS
#  ENABLE_64BIT_INTEGERS 
#  CMAKE_Fortran_FLAGS
#  CMAKE_Fortran_COMPILER
#

option(ENABLE_MATH_INT_TEST "Enable checking of math libs integer size compatibility" ON)

if (ENABLE_MATH_INT_TEST)

  if ( (MATH_LIBS OR EXPLICIT_LIBS) AND NOT ENABLE_BUILTIN_BLAS)

    message(STATUS "Configure-time Linked math libs integer compatibility test ENABLED")

    include(CheckFortranSourceCompiles) # this macro may be in the newest CMake
    include(CheckFortranSourceRuns)   # this macro may be in newest CMake
    include(CheckBlasIntegerType) # own macro

    set(CMAKE_REQUIRED_FLAGS ${CMAKE_Fortran_FLAGS})
    # you can not always combine libraries from different  sources (like MKL, SYSTEM_NATIVE)!
    if (MATH_LIBS AND EXPLICIT_LIBS)
      set(CMAKE_REQUIRED_LIBRARIES "${MATH_LIBS};${EXPLICIT_LIBS}")
    elseif(MATH_LIBS AND NOT EXPLICIT_LIBS)
      set(CMAKE_REQUIRED_LIBRARIES "${MATH_LIBS}")
    elseif(NOT MATH_LIBS AND EXPLICIT_LIBS)
      set(CMAKE_REQUIRED_LIBRARIES "${EXPLICIT_LIBS}")
    else()
      message(STATUS "Have you any blas library linked to Dirac ? Last linking step may crash ...")
    endif()

#   perform BLAS-integer test only, LAPACK-integer test is inside executable, dirac.x
    CHECK_BLAS_INTEGER_TYPE(blas_test_compiled blas_test_run blas_integer_type)

    if ("${blas_integer_type}" EQUAL 1)
      message(STATUS "We have blas integer*4")
      if (ENABLE_64BIT_INTEGERS)
        message("Your linked math libraries: ${CMAKE_REQUIRED_LIBRARIES}")
        message(WARNING "Incompatibitility between blas integer*4 and required integer*8 variables")
      endif()
    elseif ("${blas_integer_type}" EQUAL 2)
      message(STATUS "We have blas integer*8")
      if (NOT ENABLE_64BIT_INTEGERS)
        message("Your linked math libraries: ${CMAKE_REQUIRED_LIBRARIES}")
        message(WARNING "Incompatibitility between blas integer*8 and the default integer*4 variables")
      endif()
    else()
      message(STATUS "Your linked math libraries: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "Blas integer*4/8 status not determined.")
      message(STATUS "WATCH THAT YOU HAVE COMPABILE INTEGER SIZES IN LINKED BLAS,LAPACK LIBRARIES!")
    endif()
  endif()
else()
  message(STATUS "Linked math libs integer size compatibility test DISABLED")
  message(STATUS "WATCH THAT YOU HAVE COMPABILE INTEGER SIZES IN LINKED BLAS,LAPACK LIBRARIES!")
endif()
