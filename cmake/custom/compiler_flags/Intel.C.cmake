if(CMAKE_C_COMPILER_ID MATCHES Intel)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
    set(CMAKE_C_FLAGS_RELEASE "-O2")
    set(CMAKE_C_FLAGS_DEBUG "-O0")
    add_compile_options(-march=skylake-avx512) # to handle size_t and similar C99 features
endif()
