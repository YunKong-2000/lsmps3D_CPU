cmake_minimum_required(VERSION 3.20)

project(lsmps3d_cpu LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(lsmps3d_core
    src/config/config_reader.cpp
    src/core/particle_set.cpp
    src/core/simulation_config.cpp
    src/io/logger.cpp
    src/io/vtk_writer.cpp
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

add_executable(core_data_test
    tests/core_data_test.cpp
)

target_link_libraries(core_data_test
    PRIVATE
        lsmps3d_core
)

add_test(NAME core_data_test COMMAND core_data_test)

add_executable(config_reader_test
    tests/config_reader_test.cpp
)

target_link_libraries(config_reader_test
    PRIVATE
        lsmps3d_core
)

add_test(NAME config_reader_test COMMAND config_reader_test)

add_executable(vtk_writer_test
    tests/vtk_writer_test.cpp
)

target_link_libraries(vtk_writer_test
    PRIVATE
        lsmps3d_core
)

add_test(NAME vtk_writer_test COMMAND vtk_writer_test)
