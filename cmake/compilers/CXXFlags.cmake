# This variable is not guaranteed to be defined for all compilers or languages
if(NOT DEFINED CMAKE_CXX_COMPILER_ID)
	message(FATAL_ERROR "CMAKE_CXX_COMPILER_ID variable is not defined! (CMake Error)")
endif()

if (NOT DEFINED DEFAULT_CXX_FLAGS_SET OR RESET_FLAGS)

# custom flags are defined by setup --custom-cxx-flags
if(DEFINED CUSTOM_CXX_FLAGS)
    # set custom compiler flags
    set(CMAKE_CXX_FLAGS "${CUSTOM_CXX_FLAGS}")
    # unset flags for other build types
    # everything is defined in CMAKE_CXX_FLAGS variable
    # note: unset() does not work in this case as expected:
    # does not override variable in cache
    set(CMAKE_CXX_FLAGS_DEBUG "")
    set(CMAKE_CXX_FLAGS_RELEASE "")
    #set(CMAKE_CXX_FLAGS_PROFILE "")
else()
    # custom flags are not defined
    if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
        set (CMAKE_CXX_FLAGS "-g -Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused")
        set (CMAKE_CXX_FLAGS_DEBUG "-O0 -DDEBUG")
        set (CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -DNDEBUG -Wno-unused")
        #set (CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g -pg")
        if (ENABLE_CODE_COVERAGE)
            set (CMAKE_CXX_FLAGS
                "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
            set (CMAKE_CXX_LINK_FLAGS "-fprofile-arcs -ftest-coverage")
        endif()
    elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
        set (CMAKE_CXX_FLAGS "-Wno-unknown-pragmas")
        set (CMAKE_CXX_FLAGS_DEBUG "-O0 -debug -DDEBUG")
        set (CMAKE_CXX_FLAGS_RELEASE "-debug -O3 -DNDEBUG")
        #set (CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g -pg")
        set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES PGI)
        set(CMAKE_CXX_FLAGS         "-g")
        set(CMAKE_CXX_FLAGS_DEBUG   "-O0")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3")
        #set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g -pg")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES XL)
        set(CMAKE_CXX_FLAGS         "-g")
        set(CMAKE_CXX_FLAGS_DEBUG   "-O0")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3")
        #set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g -pg")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES Clang)
        set (CMAKE_CXX_FLAGS "-g -Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused")
        set (CMAKE_CXX_FLAGS_DEBUG "-O0 -DDEBUG")
        set (CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -DNDEBUG -Wno-unused")
        if (ENABLE_CODE_COVERAGE)
            set (CMAKE_CXX_FLAGS
                "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
            set (CMAKE_CXX_LINK_FLAGS "-fprofile-arcs -ftest-coverage")
        endif()

    else()
        message(FATAL_ERROR "Vendor of your C++ compiler is not supported")
    endif()

    if(DEFINED EXTRA_CXX_FLAGS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
    endif()
endif()
    
save_compiler_flags(CXX)
endif ()
