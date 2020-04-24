cmake_minimum_required(VERSION 3.1)

project(Connection)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_BUILD_TYPE Release)

find_package(CGAL)

include(${CGAL_USE_FILE})
# include_directories(${Matlab_INCLUDE_DIRS})

include_directories(eigen)

file(GLOB Headers "src/*.h")
file(GLOB Sources "src/*.cpp")

add_executable(Connection 
        ${Headers} 
        ${Sources})

target_link_libraries(Connection 
        ${CGAL_LIBRARIES} 
        ${CGAL_3RD_PARTY_LIBRARIES})