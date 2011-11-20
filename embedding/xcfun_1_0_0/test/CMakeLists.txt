
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/src/include)
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(funeval_1_0_0 funeval.c)
target_link_libraries(funeval_1_0_0 xcfun_1_0_0 m)

add_executable(testall_1_0_0 testall.c)
target_link_libraries(testall_1_0_0 xcfun_1_0_0 m)

if (XCFUN_ENABLE_TESTS)
	add_test(xcfun testall_1_0_0)
endif()
