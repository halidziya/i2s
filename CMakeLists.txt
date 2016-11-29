cmake_minimum_required(VERSION 3.6)
project(I2Slice)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
    i2s/i2gmm.cpp
    i2s/Restaurant.cpp
    i2s/Restaurant.h
    i2s/SGlobal.cpp
    i2s/SGlobal.h
    i2s/Table.cpp
    i2s/Table.h)

INCLUDE_DIRECTORIES(../FastMat/)
LINK_DIRECTORIES(../FastMat)
add_executable(i2s ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(i2s libFastMat)