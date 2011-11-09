if(CMAKE_C_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_C_FLAGS_RELEASE "-g -O2 -Wno-unused")
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -static -fpic"
            )
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES Intel)
    set(CMAKE_C_FLAGS         "-wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_C_FLAGS_RELEASE "-g -O2")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
endif()

if(CMAKE_C_COMPILER_ID MATCHES PGI)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_C_FLAGS_RELEASE "-O3")
endif()

if(CMAKE_C_COMPILER_ID MATCHES XL)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   " ")
    set(CMAKE_C_FLAGS_RELEASE " ")
endif()
