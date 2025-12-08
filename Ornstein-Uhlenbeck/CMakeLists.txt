cmake_minimum_required(VERSION 3.2)

project (sde VERSION 0.0 DESCRIPTION "" LANGUAGES CXX)


set (CMAKE_CXX_STANDARD 17)

set (CMAKE_CXX_STANDARD_REQUIRED True)

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "Release")
    message ("Hello")
endif()

option (STATICALLY_LINK "Build for release purposes (statically link everything)" OFF)
if (STATICALLY_LINK)
    set (BUILD_SHARED_LIBS OFF)
else ()
    set (BUILD_SHARED_LIBS ON)
endif ()

set (CMAKE_CXX_FLAGS
     "${CMAKE_CXX_FLAGS} \
     -std=c++17 \
     -Wall \
     -Wformat \
     -Wextra \
     -Wmissing-include-dirs \
     -fdiagnostics-color=always \
     -mavx"
)

if (CMAKE_BUILD_TYPE MATCHES "Release")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()


# 'build' Eigen library first if non-system
set (EIGEN3_DIR /opt/eigen-3.4.0)
find_package (Eigen3 3.4.0 REQUIRED)
 

# Look for ROOT. As long as an installation is available in the environment, 
# this should succeed.
message(STATUS "Looking for ROOT")
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}")



# eigen_test
add_executable (eigen_test eigen_test.cpp)
target_link_libraries (eigen_test PRIVATE Eigen3::Eigen)


# sde
add_executable (sde sde.cpp)
target_link_libraries (sde ${ROOT_LIBRARIES})




# solve
#find_package (Ceres REQUIRED)
#add_executable (solve solve.cpp)
#target_link_libraries (solve Ceres::ceres)


