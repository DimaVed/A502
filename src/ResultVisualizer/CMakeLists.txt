cmake_minimum_required(VERSION 3.21)

project(
	ResultVisualizer
    LANGUAGES CXX
)

include_directories(${Boost_INCLUDE_DIRS})
file (GLOB_RECURSE  Src  *.cpp  *.h )
add_executable(ResultVisualizer   ${Src} )

if(USE_VTK_SUPPORT)
    target_compile_definitions(ResultVisualizer PRIVATE WITH_VTK)
endif()
