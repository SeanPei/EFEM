cmake_minimum_required (VERSION 3.5)

project ("EFEM")
set(CMAKE_BUILD_TYPE "Release")
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
    message(STATUS "optional:-std=c++11")
endif(CMAKE_COMPILER_IS_GNUCXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

if (WIN32)
	add_definitions(-DWIN32)
	add_definitions(-D_ENABLE_EXTENDED_ALIGNED_STORAGE)
	set(CMAKE_CXX_FLAGS "/utf-8 ${CMAKE_CXX_FLAGS}")
endif(WIN32)

add_definitions(-DEIGEN_DONT_PARALLELIZE)
add_definitions(-DEIGEN_USE_MKL_ALL)

# 加载CMake函数
include(EFEMCmakeUtilities)

# Boost
set(Boost_USE_STATIC_LIBS ON) 
find_package(Boost REQUIRED COMPONENTS filesystem system)
if (Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
else ()
	message(FATAL_ERROR "Need boost library, but not found!")
endif()
# VTK
set(VTK_USE_STATIC_LIBS ON) 
find_package(VTK REQUIRED)

# external库
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external)
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIR})
# 附加包含目录
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# 可执行文件输出目录
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# 头文件和源文件
set(EFEM_SOURCE_DIRS model core io material math fe simulation main)
foreach(DIR IN LISTS EFEM_SOURCE_DIRS)
  add_subdirectory(src/${DIR})
endforeach()