cmake_minimum_required(VERSION 3.10)

project(dim VERSION "0.0.0" LANGUAGES CXX)

set(MASTER_PROJECT OFF)
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    set(MASTER_PROJECT ON)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ${MASTER_PROJECT})

option(DIM_TESTING "dim: enables building of unit tests [default: ${MASTER_PROJECT}" ${MASTER_PROJECT})
option(DIM_EMBEDDED_FMTLIB "dim: uses embedded fmtlib [default: ${MASTER_PROJECT}" ${MASTER_PROJECT})
option(DIM_EMBEDDED_CATCH2 "dim: uses embedded catch2 for testing [default: ${MASTER_PROJECT}]" ${MASTER_PROJECT})

# setting defaults
if (NOT("${CMAKE_CXX_STANDARD}"))
    set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang"))
    add_compile_options(-Wall)
    add_compile_options(-Wextra)
elseif(DEFINED MSVC)
    add_definitions(-DNOMINMAX)
endif()

option(DIM_ENABLE_TIDY "Enable clang-tidy [default: OFF]" OFF)
if(DIM_ENABLE_TIDY)
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

if(DIM_EMBEDDED_FMTLIB)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/fmt" EXCLUDE_FROM_ALL)
    add_definitions(-DFMT_USE_WINDOWS_H=0)
else()
    # master project must provide its own fmtlib
endif()

if(DIM_TESTING AND DIM_EMBEDDED_CATCH2)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/catch2")
else()
    # master project must provide its own catch2
endif()

# ----------------------------------------------------------------------------
# LIBRARY

set(dim_SRCS
    dim/isqrt.h
    dim/vec.h
    dim/vec_ostream.h
    dim/mat.h
    dim/mat_ostream.h
    dim/value_traits.h
)
add_library(dim INTERFACE)
target_include_directories(dim INTERFACE
    $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(dim INTERFACE cxx_std_17)
target_link_libraries(dim INTERFACE fmt::fmt-header-only)

# ----------------------------------------------------------------------------
# TEST

if(DIM_TESTING)
    enable_testing()
    add_executable(test_dim
        dim/test_dim.cpp
        dim/vec_test.cpp
        dim/mat_test.cpp
        dim/mat_solver_test.cpp
        dim/util_test.cpp
    )
    target_link_libraries(test_dim dim fmt::fmt-header-only Catch2::Catch2)
    add_test(test_dim ./test_dim)
endif()
