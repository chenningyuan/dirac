add_library(xcfun)

target_sources(xcfun
  PRIVATE
    src/xcfun.cpp
    src/evaluators.cpp
    src/functional.cpp
    src/parameters.cpp
    src/potential.cpp
    src/fortran.c
    src/functionals/beckex.cpp
    src/functionals/brx.cpp
    src/functionals/ktx.cpp
#   src/functionals/lb94.cpp
    src/functionals/ldaerfc.cpp
    src/functionals/ldaerfc_jt.cpp
    src/functionals/ldaerfx.cpp
    src/functionals/list_of_functionals.cpp
    src/functionals/lypc.cpp
    src/functionals/m05c.cpp
    src/functionals/m05x.cpp
    src/functionals/m05x2c.cpp
    src/functionals/m05x2x.cpp
    src/functionals/m06c.cpp
    src/functionals/m06hfx.cpp
    src/functionals/m06lc.cpp
    src/functionals/m06lx.cpp
    src/functionals/m06x.cpp
    src/functionals/m06x2c.cpp
    src/functionals/m06x2x.cpp
    src/functionals/optx.cpp
    src/functionals/p86c.cpp
    src/functionals/pbec.cpp
    src/functionals/pbex.cpp
    src/functionals/pw86x.cpp
    src/functionals/pw91k.cpp
    src/functionals/pw91x.cpp
    src/functionals/pw92c.cpp
    src/functionals/revpbex.cpp
    src/functionals/revtpssc.cpp
    src/functionals/revtpssx.cpp
    src/functionals/rpbex.cpp
    src/functionals/slaterx.cpp
    src/functionals/spbec.cpp
    src/functionals/tfk.cpp
    src/functionals/tpssc.cpp
    src/functionals/tpssx.cpp
    src/functionals/vwn5c.cpp
  )

target_include_directories(xcfun
  PUBLIC
    include
  PRIVATE
    src
    src/taylor
    src/functionals
  )

#-------------------------------------------------------------------------------

option(XCFUN_ENABLE_TESTS "Enable tests" ON)

if(XCFUN_ENABLE_TESTS)
  add_executable(testall)

  target_sources(testall
    PRIVATE
      test/testall.c
    )

  target_link_libraries(testall
    PRIVATE
      xcfun
    )

  add_test(xcfun ${PROJECT_BINARY_DIR}/src/xcfun/testall)
endif()

#-------------------------------------------------------------------------------

add_library(xcfun_fortran_bindings)

target_sources(xcfun_fortran_bindings
  PRIVATE
    fortran/xcfun_autogen.F90
    fortran/xcfun_module.F90
  )

target_link_libraries(xcfun_fortran_bindings
  PRIVATE
    xcfun
  )
