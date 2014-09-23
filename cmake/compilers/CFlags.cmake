if(NOT DEFINED DEFAULT_C_FLAGS_SET OR RESET_FLAGS)

# custom flags are defined by setup --custom-cc-flags
if(DEFINED CUSTOM_C_FLAGS)
    # set custom compiler flags
    set(CMAKE_C_FLAGS "${CUSTOM_C_FLAGS}")
    # unset flags for other build types
    # everything is defined in CMAKE_C_FLAGS variable
    # note: unset() does not work in this case as expected:
    # does not override variable in cache
    set(CMAKE_C_FLAGS_DEBUG "")
    set(CMAKE_C_FLAGS_RELEASE "")
    #set(CMAKE_C_FLAGS_PROFILE "")
else()
    # custom flags are not defined
    if(CMAKE_C_COMPILER_ID MATCHES GNU)
        set(CMAKE_C_FLAGS         "-g")
        set(CMAKE_C_FLAGS_DEBUG   "-O0")
        set(CMAKE_C_FLAGS_RELEASE "-O2 -Wno-unused")
        #set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -g -pg")        
	if(ENABLE_STATIC_LINKING)
            set(CMAKE_C_FLAGS
                "${CMAKE_C_FLAGS} -static -fpic"
                )
        endif()
    elseif(CMAKE_C_COMPILER_ID MATCHES Intel)
        set(CMAKE_C_FLAGS         "-g -wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
        set(CMAKE_C_FLAGS_DEBUG   "-O0")
        set(CMAKE_C_FLAGS_RELEASE "-O2")
        
        # intel compiler xHOST flag
        if(XHOST_AVAILABLE)
            set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -xHOST")
        endif(XHOST_AVAILABLE)
        
        #set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -g -pg")
        set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")

        if(DEFINED MKL_FLAG)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MKL_FLAG}")
        endif()
    elseif(CMAKE_C_COMPILER_ID MATCHES PGI)
        set(CMAKE_C_FLAGS         "-g")
        set(CMAKE_C_FLAGS_DEBUG   "-O0")
        set(CMAKE_C_FLAGS_RELEASE "-O3")
        #set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -g -pg")
    elseif(CMAKE_C_COMPILER_ID MATCHES XL)
        set(CMAKE_C_FLAGS         "-qcpluscmt")
        set(CMAKE_C_FLAGS_DEBUG   " ")
        set(CMAKE_C_FLAGS_RELEASE "-O3")
        #set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -g -pg")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
        set(CMAKE_C_FLAGS         "-g")
        set(CMAKE_C_FLAGS_DEBUG   "-O0")
        set(CMAKE_C_FLAGS_RELEASE "-O2 -Wno-unused")
        if(ENABLE_STATIC_LINKING)
            set(CMAKE_C_FLAGS
                "${CMAKE_C_FLAGS} -Bstatic -fpic"
                )
        endif()
    else()
        message(FATAL_ERROR "Vendor of your C compiler is not supported")
    endif()

    if(DEFINED EXTRA_C_FLAGS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_C_FLAGS}")
    endif()
endif()

save_compiler_flags(C)
endif()
