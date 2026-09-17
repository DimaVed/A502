cmake_minimum_required(VERSION 3.21)

project(
	Solver
    LANGUAGES CXX
)

include_directories(${Boost_INCLUDE_DIRS})
file (GLOB_RECURSE  Src  *.cpp  *.h )
add_executable(Solver   ${Src} )

