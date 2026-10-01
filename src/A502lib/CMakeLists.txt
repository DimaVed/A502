cmake_minimum_required(VERSION 3.21)

project(
	A502lib
    LANGUAGES CXX
)

include_directories(${Boost_INCLUDE_DIRS})

if(USE_VTK_SUPPORT)
    include_directories(${VTK_INCLUDE_DIRS})
    add_compile_definitions(WITH_VTK)
endif()

file (GLOB_RECURSE  Src  *.cpp  *.h )
add_library(A502lib STATIC  ${Src} )

# Consumers (Visualizer, Solver, Tests) include "AeroModel.h" / "Logger.h".
target_include_directories(A502lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

