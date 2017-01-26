cmake_minimum_required(VERSION 2.6)
project(I2Slice)

macro(use_cxx11)
    if (CMAKE_VERSION VERSION_LESS "3.1")
        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            set (CMAKE_CXX_FLAGS "--std=gnu++11 -O3 ${CMAKE_CXX_FLAGS}")
        endif ()
    else ()
        set (CMAKE_CXX_STANDARD 11)
    endif ()
endmacro(use_cxx11)


set(CMAKE_CXX_STANDARD 11)
use_cxx11()


set(SOURCE_FILES
    i2s/i2gmm.cpp
    i2s/Restaurant.cpp
    i2s/Restaurant.h
    i2s/Table.cpp
    i2s/Table.h)

INCLUDE_DIRECTORIES(../GMMBase/)
LINK_DIRECTORIES(../GMMBase)
INCLUDE_DIRECTORIES(../FastMat/)
LINK_DIRECTORIES(../FastMat)
find_package (Threads)
add_executable(i2slice ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(i2slice FastMat)
TARGET_LINK_LIBRARIES (i2slice ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(i2slice GMMBase)
