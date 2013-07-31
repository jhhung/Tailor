cmake_minimum_required (VERSION 2.6)
project (Tailor)

# Initialize CXXFLAGS.
set(CMAKE_CXX_FLAGS "-Ofast -std=c++11")
#if your boost is installed in somewhere CMake cannot find, please specified as below:
#SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "$HOME/Downloads/boost_1_51_0/")
#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "$HOME/Downloads/boost_1_51_0/lib")
# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

include (${CMAKE_ROOT}/Modules/FindBoost.cmake)

#boost
set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package( Boost COMPONENTS date_time filesystem system serialization system filesystem regex thread iostreams program_options)
find_package( Threads )
find_package( ZLIB )

if(Boost_FOUND AND Threads_FOUND AND ZLIB_FOUND)
	include_directories(${Boost_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS} ) 
# add the executable
	add_executable(bin/tailor src/main.cpp)
	add_definitions(-DmyALPHABET)
	target_link_libraries( bin/tailor ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_SOURCE_DIR}/lib/libabwt_table.a)
endif()
