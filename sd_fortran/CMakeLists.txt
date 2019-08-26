cmake_minimum_required (VERSION 2.8)

enable_language (Fortran)
project (SCE)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS ${CMAKE_Fortran_SOURCE_FILE_EXTENSIONS} f03)
set(CMAKE_Fortran_FLAGS "-std=f2008 -Wall -Wextra")

set(SOURCE_FILES
    util.f90
    collisions.f90
    constants.f90
    droplet.f90
    math.f90
    sce.f90)

include_directories( /usr/local/include/fgsl )
link_directories( /usr/local/lib/ )

add_executable(sce ${SOURCE_FILES})
target_link_libraries( sce fgsl )
set_target_properties(sce PROPERTIES LINKER_LANGUAGE Fortran)
