cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(project_inliers)

find_package(PCL 1.2 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (project_inliers project_inliers.cpp)
target_link_libraries (project_inliers ${PCL_LIBRARIES})
 

FIND_PACKAGE (QT4 REQUIRED)

INCLUDE (${QT_USE_FILE})

INCLUDE_DIRECTORIES (${CMAKE_CURRENT_BINARY_DIR})

ADD_DEFINITIONS (${QT_DEFINITIONS})



