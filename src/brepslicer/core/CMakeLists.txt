cmake_minimum_required(VERSION 3.15)
project(BrepSlicer VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(BIN_3RD_LIST)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_SOURCE_DIR}/_bin/x64_Debug)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/_bin/x64_Release)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE})

if(MSVC)
    include(cmake/msvc/occt_bambu_deps.cmake)
    include(cmake/msvc/copy_bin.cmake)
    add_compile_options(/utf-8 /MP /EHsc /bigobj)
    add_compile_definitions(NOMINMAX _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES)
    add_compile_options(/wd4244 /wd4267 /wd4996)
else()
    message(FATAL_ERROR "L0/L1 bootstrap currently targets MSVC. Linux cmake can be added later.")
endif()

# Slice algorithm: no OCCT headers, no OCCT link. Talks only to IShape / IBrepKernel.
add_library(brepslicer_algo STATIC
    src/io/JsonWriter.cpp
    src/io/SvgWriter.cpp
    src/io/DxfWriter.cpp
    src/geom/GeomUtil.cpp
    src/topo/FaceIndex.cpp
    src/intersect/QuadricPlane.cpp
    src/intersect/AnalyticIntersect.cpp
    src/intersect/BSplineFit.cpp
    src/intersect/UvMarch.cpp
    src/intersect/UvMatch.cpp
    src/contour/ContourAssembler.cpp
    src/slicer/Slicer.cpp
    src/verify/SliceVerify.cpp
)
target_include_directories(brepslicer_algo
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# OCC adapter: IShape/IKernel/IBuilder/ISectionRef + ShapeEngine strategy (D:/OCC pattern).
add_library(brepslicer_occ STATIC
    src/occ/OccShape.cpp
    src/occ/OccKernel.cpp
    src/occ/OccBuilder.cpp
    src/occ/OccSectionRef.cpp
    src/occ/OccStrategy.cpp
)
target_include_directories(brepslicer_occ
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(brepslicer_occ PUBLIC brepslicer_algo OCCT750)

function(brepslicer_add_app target)
    add_executable(${target} ${ARGN})
    target_link_libraries(${target} PRIVATE brepslicer_occ)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    add_dependencies(${target} _ready_bin)
    set_property(TARGET ${target} PROPERTY VS_DEBUGGER_ENVIRONMENT "PATH=${OCCT750_BINARY_DIR};%PATH%")
    set_property(TARGET ${target} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/_bin/x64_Release")
endfunction()

brepslicer_add_app(brep_slice    apps/brep_slice/main.cpp)
brepslicer_add_app(slice_verify  apps/slice_verify/main.cpp)
brepslicer_add_app(slice_selftest apps/selftest/main.cpp)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT brep_slice)
