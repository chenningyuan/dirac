option(ENABLE_INTEREST "Enable Interest" ON)

if(ENABLE_INTEREST)
    message(STATUS "Interest library enabled")
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/interest)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} interest)
endif()
