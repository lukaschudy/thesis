cmake_minimum_required(VERSION 3.16)
project(ogcode_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(ogcode_cpp
    src/main.cpp
    src/fortran_runtime.cpp
    src/generic_routines.cpp
    src/globals.cpp
    src/QL_routines.cpp
    src/PI_routines.cpp
    src/LearningSimulation.cpp
    src/ConvergenceResults.cpp
    src/EquilibriumCheck.cpp
    src/QGapToMaximum.cpp
    src/ImpulseResponse.cpp
    src/LearningTrajectory.cpp
    src/DetailedAnalysis.cpp
)

target_include_directories(ogcode_cpp PRIVATE include)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(ogcode_cpp PRIVATE OpenMP::OpenMP_CXX)
endif()
