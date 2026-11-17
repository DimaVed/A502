cmake_minimum_required(VERSION 3.21)

project(
	Tests
    LANGUAGES CXX
)

include_directories(${Boost_INCLUDE_DIRS})
file (GLOB_RECURSE  Src  *.cpp  *.h )
add_executable(Tests   ${Src} )

# Tests locate the sample decks from the source tree.
target_compile_definitions(Tests PRIVATE A502_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

