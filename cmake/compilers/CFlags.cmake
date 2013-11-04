if(NOT DEFINED DEFAULT_C_FLAGS_SET OR RESET_FLAGS)

if(CMAKE_C_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_FLAGS         "-g")
    set(CMAKE_C_FLAGS_DEBUG   "-O0")
    set(CMAKE_C_FLAGS_RELEASE "-O2 -Wno-unused")
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -static -fpic"
            )
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES Intel)
    set(CMAKE_C_FLAGS         "-g -wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
    set(CMAKE_C_FLAGS_DEBUG   "-O0")
    set(CMAKE_C_FLAGS_RELEASE "-O2 -xHost")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")

    if(DEFINED MKL_FLAG)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MKL_FLAG}")
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES PGI)
    set(CMAKE_C_FLAGS         "-g")
    set(CMAKE_C_FLAGS_DEBUG   "-O0")
    set(CMAKE_C_FLAGS_RELEASE "-O3")
endif()

if(CMAKE_C_COMPILER_ID MATCHES XL)
    set(CMAKE_C_FLAGS         "-qcpluscmt")
    set(CMAKE_C_FLAGS_DEBUG   " ")
    set(CMAKE_C_FLAGS_RELEASE "-O3")
endif()

if(DEFINED EXTRA_C_FLAGS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_C_FLAGS}")
endif()

save_compiler_flags(C)
endif()
