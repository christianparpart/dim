cmake_minimum_required(VERSION 3.10)

project(smath VERSION "0.0.0" LANGUAGES CXX)

set(MASTER_PROJECT OFF)
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    set(MASTER_PROJECT ON)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ${MASTER_PROJECT})

option(SMATH_TESTING "smath: enables building of unit tests [default: ${MASTER_PROJECT}" ${MASTER_PROJECT})
option(SMATH_EMBEDDED_FMTLIB "smath: uses embedded fmtlib [default: ${MASTER_PROJECT}" ${MASTER_PROJECT})
option(SMATH_EMBEDDED_CATCH2 "smath: uses embedded catch2 for testing [default: ${MASTER_PROJECT}]" ${MASTER_PROJECT})

# setting defaults
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang"))
    add_compile_options(-Wall)
    add_compile_options(-Wextra)
elseif(DEFINED MSVC)
    add_definitions(-DNOMINMAX)
endif()

option(SMATH_ENABLE_TIDY "Enable clang-tidy [default: OFF]" OFF)
if(SMATH_ENABLE_TIDY)
    find_program(CLANG_TIDY_EXE
        NAMES clang-tidy-9 clang-tidy-8 clang-tidy-7 clang-tidy
        DOC "Path to clang-tidy executable")
    if(NOT CLANG_TIDY_EXE)
        message(STATUS "clang-tidy not found.")
    else()
        message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
        set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
    endif()
endif()

# ----------------------------------------------------------------------------
# dependencies

if(SMATH_EMBEDDED_FMTLIB)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/fmt" EXCLUDE_FROM_ALL)
    add_definitions(-DFMT_USE_WINDOWS_H=0)
else()
    # master project must provide its own fmtlib
endif()

if(SMATH_TESTING AND SMATH_EMBEDDED_CATCH2)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/catch2")
else()
    # master project must provide its own catch2
endif()

# ----------------------------------------------------------------------------
# LIBRARY

set(smath_SRCS
    smath/isqrt.h
    smath/vec.h
    smath/vec_ostream.h
    smath/mat.h
    smath/mat_ostream.h
    smath/value_traits.h
)
add_library(smath INTERFACE)
target_include_directories(smath INTERFACE
    $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(smath INTERFACE cxx_std_17)
target_link_libraries(smath INTERFACE fmt::fmt-header-only)

# ----------------------------------------------------------------------------
# TEST

if(SMATH_TESTING)
    enable_testing()
    add_executable(test_smath
        smath/test_smath.cpp
        # TODO: smath/vec_test.cpp
        smath/mat_test.cpp
        smath/util_test.cpp
    )
    target_link_libraries(test_smath smath fmt::fmt-header-only Catch2::Catch2)
    add_test(test_smath ./test_smath)
endif()
