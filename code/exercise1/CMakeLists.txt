CMAKE_MINIMUM_REQUIRED(VERSION 3.2) 

PROJECT(SPLINE)  

LINK_DIRECTORIES("/usr/local/lib")
LINK_DIRECTORIES("/usr/lib")
LINK_DIRECTORIES("/usr/include")
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14") # 添加c++14标准支持

#default compiler
SET(CMAKE_C_COMPILER "/usr/bin/gcc") 
SET(CMAKE_CXX_COMPILER "/usr/bin/g++") 

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin") # 可执行文件输出目录

add_executable(testDirichlet
    ./src/Dirichlettest.cpp
)

add_executable(testNeumann
    ./src/Neumanntest.cpp
)

add_executable(testRobin
    ./src/Robintest.cpp
)
target_link_libraries(testDirichlet jsoncpp)
target_link_libraries(testNeumann jsoncpp)
target_link_libraries(testRobin jsoncpp)