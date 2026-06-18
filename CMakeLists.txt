cmake_minimum_required(VERSION 3.20)

project(lsmps3d_cpu LANGUAGES CXX)

find_package(Eigen3 REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(lsmps3d_core
    src/config/config_reader.cpp
    src/core/particle_set.cpp
    src/core/simulation_config.cpp
    src/io/logger.cpp
    src/io/vtk_writer.cpp
    src/free_surface/free_surface_detector.cpp
    src/lsmps/lsmps_basis.cpp
    src/lsmps/lsmps_matrices.cpp
    src/lsmps/weight_function.cpp
    src/neighbor/cell_linked_list.cpp
    src/neighbor/neighbor_search.cpp
)

target_include_directories(lsmps3d_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_options(lsmps3d_core
    PRIVATE
        $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall -Wextra -Wpedantic>
)

target_link_libraries(lsmps3d_core
    PUBLIC
        Eigen3::Eigen
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

add_executable(neighbor_search_test
    tests/neighbor_search_test.cpp
)

target_link_libraries(neighbor_search_test
    PRIVATE
        lsmps3d_core
)

add_test(NAME neighbor_search_test COMMAND neighbor_search_test)

add_executable(neighbor_search_large_case_test
    tests/neighbor_search_large_case_test.cpp
)

target_link_libraries(neighbor_search_large_case_test
    PRIVATE
        lsmps3d_core
)

add_test(NAME neighbor_search_large_case_test COMMAND neighbor_search_large_case_test)

add_executable(free_surface_detector_test
    tests/free_surface_detector_test.cpp
)

target_link_libraries(free_surface_detector_test
    PRIVATE
        lsmps3d_core
)

add_test(NAME free_surface_detector_test COMMAND free_surface_detector_test)

add_executable(lsmps_matrix_test
    tests/lsmps_matrix_test.cpp
)

target_link_libraries(lsmps_matrix_test
    PRIVATE
        lsmps3d_core
)

add_test(NAME lsmps_matrix_test COMMAND lsmps_matrix_test)

add_executable(lsmps_operator_validation_test
    tests/lsmps_operator_validation_test.cpp
)

target_link_libraries(lsmps_operator_validation_test
    PRIVATE
        lsmps3d_core
)

add_test(NAME lsmps_operator_validation_test COMMAND lsmps_operator_validation_test)
