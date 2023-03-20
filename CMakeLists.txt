cmake_minimum_required(VERSION 3.20)
project(2)

include_directories(/usr/local/include)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_COMPILER g++)  # from hugo
set(CMAKE_CXX_FLAGS "-Wall -std=c++20")  # from hugo
set(CMAKE_BUILD_TYPE "Debug")  # from hugo
set(CMAKE_CXX_FLAGS "-O2" )  # from hugo


# CGAL  # from Hugo
find_package( CGAL QUIET COMPONENTS )
if ( CGAL_FOUND )
    include( ${CGAL_USE_FILE} )
    message(STATUS ${CGAL_LIBRARIES})
    message(STATUS ${CGAL_3RD_PARTY_LIBRARIES})
else()
    message(SEND_ERROR "CGAL library is required")
    return()
endif()

find_package(Eigen3 REQUIRED)

include(CGAL_Eigen3_support)

include_directories( ${CMAKE_SOURCE_DIR}/include/ )

FILE(GLOB SRC_FILES src/*.cpp)
add_executable(2 ${SRC_FILES} src/test.cpp src/json_to_obj.cpp src/Hemisphericality.cpp src/understanding_city_json.cpp src/compilation_v1.cpp)


#add_executable(2 main.cpp)

find_package( CGAL REQUIRED )
#find_package(Eigen3 REQUIRED)
find_package(Eigen3 3.1.0 QUIET)
#find_package( CGAL Eigen3_support )

target_link_libraries(${PROJECT_NAME} CGAL::CGAL CGAL::Eigen3_support)

#target_link_libraries(${PROJECT_NAME} CGAL::Eigen3_support )

