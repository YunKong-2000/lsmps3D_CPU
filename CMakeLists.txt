cmake_minimum_required(VERSION 3.20)

project(lsmps3d_cpu LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(lsmps3d_core
    src/io/logger.cpp
)

target_include_directories(lsmps3d_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_options(lsmps3d_core
    PRIVATE
        $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall -Wextra -Wpedantic>
)

add_executable(lsmps3d
    src/main.cpp
)

target_link_libraries(lsmps3d
    PRIVATE
        lsmps3d_core
)

enable_testing()

add_executable(smoke_test
    tests/smoke_test.cpp
)

target_link_libraries(smoke_test
    PRIVATE
        lsmps3d_core
)

add_test(NAME smoke_test COMMAND smoke_test)
