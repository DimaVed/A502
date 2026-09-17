cmake_minimum_required(VERSION 3.21)

project(
	Visualizer
    LANGUAGES CXX
)

include_directories(${Boost_INCLUDE_DIRS})
file (GLOB_RECURSE  Src  *.cpp  *.h )
add_executable(Visualizer   ${Src} )

