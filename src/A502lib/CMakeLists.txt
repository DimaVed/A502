cmake_minimum_required(VERSION 3.21)

project(
	A502lib
    LANGUAGES CXX
)

include_directories(${Boost_INCLUDE_DIRS})
file (GLOB_RECURSE  Src  *.cpp  *.h )
add_library(A502lib STATIC  ${Src} )

