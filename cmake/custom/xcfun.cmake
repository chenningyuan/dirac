option(ENABLE_XCFUN "Enable XCFun library" ON)

if(ENABLE_XCFUN)
    add_subdirectory(src/xcfun)
    set(EXTERNAL_LIBS xcfun ${EXTERNAL_LIBS})
    include_directories(${PROJECT_BINARY_DIR}/src/xcfun/fortran)
else()
    add_definitions(-DEXCLUDE_XCFUN)
endif()
