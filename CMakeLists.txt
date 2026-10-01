# Установка минимальной версии CMake (более старые не умеют находить библиотеки)

cmake_minimum_required (VERSION 3.20)
# Создание решения c именем Proj
project (A502 LANGUAGES CXX)

enable_testing()

option(USE_VTK_SUPPORT "Build with VTK support" ON)
# Создаем переменную в GUI как выпадающий список (CACHE STRING + STRINGS property)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build")
if(MSVC)
    add_compile_options(/utf-8)
endif()



if(WIN32)
    set(CURRENT_VCPKG_ROOT $ENV{VCPKG_ROOT})
    set(VCPKG_TARGET_TRIPLET "x64-windows-static")
elseif(UNIX) 
    set(CURRENT_VCPKG_ROOT "$ENV{HOME}/vcpkg")
    set(VCPKG_TARGET_TRIPLET "x64-linux-release") 
endif()

message( STATUS "vcpkg root:" ${CURRENT_VCPKG_ROOT})
message( STATUS "vcpkg triplet:" ${VCPKG_TARGET_TRIPLET})
set(VTKCompileTools_DIR "${CURRENT_VCPKG_ROOT}/installed/${VCPKG_TARGET_TRIPLET}/share/vtk-compile-tools" CACHE PATH "")

include( ${CURRENT_VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake)
include_directories( ${INCLUDE_DIRECTORIES}  ${CURRENT_VCPKG_ROOT}/installed/${VCPKG_TARGET_TRIPLET}/include/ )
find_package(Catch2 CONFIG REQUIRED)
include(Catch)
# windows.h определяет макросы min и max из-за этого проигрывают cтандартные функции
 add_compile_definitions(NOMINMAX)
# настройка статической линковки Boost
set(Boost_USE_STATIC_LIBS       ON)
set(Boost_USE_MULTITHREADED     ON)
set(Boost_USE_STATIC_RUNTIME    ON)
# Поиск компонентов boost

find_package(Boost  REQUIRED COMPONENTS json log program_options locale serialization)
find_package(nlohmann_json  REQUIRED)
if(USE_VTK_SUPPORT)

    add_compile_definitions(WITH_VTK)
    find_package(VTK REQUIRED)
	if(VTK_LIBRARIES MATCHES "static")
    set_property(TARGET ${PROJECT_NAME} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()
endif()

find_package(OpenCASCADE REQUIRED)
set(OCCT_LIBS 
   TKernel 
   TKMath 
   TKG2d 
   TKG3d 
   TKGeomBase 
   TKGeomAlgo 
   TKBRep 
   TKTopAlgo 
   TKPrim 
   TKBool 
   TKBO 
   TKHLR 
   TKMesh 
   TKOpenGl 
   TKService 
   TKCAF 
   TKLCAF 
   TKV3d 
   TKVCAF 
   TKXCAF
   TKXSBase 
   TKDE 
   TKDEIGES 
   TKDESTEP
   TKOffset 
)


set(MKL_INTERFACE_FULL "intel_ilp64")
set(MKL_THREADING "sequential") # или "tbb"/"omp", если используете
set(MKL_LINK "static") 
set(MKL_USE_STATIC_LIBS ON)
find_package(MKL CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(HDF5 CONFIG REQUIRED)
if(NOT TARGET HDF5::HDF5 AND TARGET hdf5::hdf5-static)
    add_library(HDF5::HDF5 ALIAS hdf5::hdf5-static)
elseif(NOT TARGET HDF5::HDF5 AND TARGET hdf5::hdf5-shared)
    add_library(HDF5::HDF5 ALIAS hdf5::hdf5-shared)
endif()
find_package(HighFive CONFIG REQUIRED)

# Рекомендуется включать профилирование только в определенных конфигурациях
option(TRACY_ENABLE "Enable Tracy profiling" OFF)

if(TRACY_ENABLE)
    find_package(Tracy REQUIRED)
    # Передает макрос TRACY_TIMER_FALLBACK всем компилируемым файлам
    add_compile_definitions(TRACY_TIMER_FALLBACK)
    message(STATUS "Tracy profiling enabled with timer fallback")
endif()





# Линковка boost 
if(Boost_FOUND)
   message(STATUS "Boost was find!")
   #include_directories(${Boost_INCLUDE_DIRS})
   message(STATUS "Boost include dir " ${Boost_INCLUDE_DIRS} )

   else()
    message(STATUS "Fail  Boost weren't find!")
endif() 

#Настройка автокомпиляции автогенерируемых файлов QT
set(CMAKE_AUTOMOC TRUE)
set(CMAKE_AUTOUIC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# поиск потребных библиотек QT
#find_package(Qt5 COMPONENTS Widgets Core REQUIRED)


add_subdirectory (src/A502lib)
add_subdirectory (src/Tests)
add_subdirectory (src/Visualizer)
add_subdirectory (src/Solver)


include_directories(${HDF5_INCLUDE_DIRS})
set (TargetLIst  A502lib Tests Visualizer Solver )

set (LibLinkList  Tests Visualizer Solver  )



foreach(exec  ${TargetLIst})
    message( STATUS ${exec})
	# Установка 20 стандарта С++
	set_target_properties (${exec} PROPERTIES CXX_STANDARD 20 CMAKE_CXX_STANDARD_REQUIRED ON)
	# Установка статической линковки для MSVC
	set_property (TARGET ${exec} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	# Установка флага многопоточной компиляции
	if (MSVC)
	target_compile_options (${exec} PRIVATE "/MP")
	target_compile_options(${exec} PRIVATE "/bigobj")
	endif()
	if(MSVC)
   target_compile_options(${exec} PRIVATE /W3 )
	else()
    #target_compile_options(${exec} PRIVATE -Wall -Wextra -Wpedantic )
	endif()
	# Прекомпилированный заголовок (добавляется, только если файл существует)
	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Console/pch.h)
	target_precompile_headers(${exec} PRIVATE Console/pch.h)
	endif()
	# линковака бибилиотек подтянутых менеджером пакетов
	target_link_libraries(${exec} PRIVATE Eigen3::Eigen)
	target_link_libraries(${exec} PRIVATE HighFive)
	target_link_libraries(${exec}   PRIVATE ${Boost_LIBRARIES})
	target_include_directories(${exec} PRIVATE ${Boost_INCLUDE_DIRS})
	target_include_directories(${exec} PRIVATE ${OpenCASCADE_INCLUDE_DIRS})
    target_link_libraries(${exec} PRIVATE ${OCCT_LIBS})
	target_link_libraries(${exec} PRIVATE MKL::MKL)
	# Catch2WithMain provides main(); only Tests should use it.
	if(exec STREQUAL "Tests")
		target_link_libraries(${exec} PRIVATE Catch2::Catch2 Catch2::Catch2WithMain)
	endif()
	target_link_libraries(${exec} PRIVATE nlohmann_json::nlohmann_json)
    
	#target_link_libraries(${exec} PRIVATE Tracy::TracyClient)
	if(USE_VTK_SUPPORT)
    target_link_libraries(${exec} PRIVATE ${VTK_LIBRARIES})
    
    endif()
	
	
endforeach()



foreach(exec  ${LibLinkList})
target_link_libraries (${exec}  PRIVATE A502lib)
endforeach()
message( STATUS ${CMAKE_CURRENT_SOURCE_DIR})

