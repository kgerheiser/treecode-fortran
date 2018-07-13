cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(treecode)

enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone -check all -warn all")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

add_subdirectory(src)
add_subdirectory(tests)

add_executable(treecode src/main.f90)
target_include_directories(treecode PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_link_libraries(treecode libtreecode)
