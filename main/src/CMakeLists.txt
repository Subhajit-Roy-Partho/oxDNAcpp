cmake_minimum_required(VERSION 3.10)
project(terminal)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("/usr/local/include")
link_directories("/usr/local/lib")

# In mac, brew install openmp seperately so that gcc headers are not overwrittern.
if(APPLE)
    link_directories("/usr/local/opt/libomp/lib")
    include_directories("/usr/local/opt/libomp/include")
endif(APPLE)


find_package(yaml-cpp CONFIG REQUIRED)
message(STATUS "Used config file: ${yaml-cpp_CONFIG}")

add_library(libraries
    SimulationManager.cpp
    Analysis.cpp
    AnalysisMulticore.cpp
)

add_executable(terminal terminal.cpp)
# target_link_libraries(terminal yaml-cpp)
target_link_libraries(terminal PUBLIC libraries yaml-cpp)
target_include_directories(terminal PUBLIC
    "${PROJECT_BINARY_DIR}"
)

add_executable(test test.cpp)
target_link_libraries(test PUBLIC libraries)
target_include_directories(test PUBLIC
    "${PROJECT_BINARY_DIR}"
)