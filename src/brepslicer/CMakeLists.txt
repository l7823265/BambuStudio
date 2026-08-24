# BrepSlicer (L0/L1/L2-v3) — parametric B-rep slicing, linked into libslic3r.
# Default: in-tree core at src/brepslicer/core (vendored for fork builds).
# Override with -DBREPSLICER_CORE_DIR=... if needed.

get_filename_component(_brepslicer_in_tree "${CMAKE_CURRENT_SOURCE_DIR}/core" ABSOLUTE)
get_filename_component(_brepslicer_sibling
    "${CMAKE_SOURCE_DIR}/../BrepSlicer-core-L2-src-20260819" ABSOLUTE)

if(EXISTS "${_brepslicer_in_tree}/include/brepslicer/Slicer.h")
    set(_brepslicer_default "${_brepslicer_in_tree}")
elseif(EXISTS "${_brepslicer_sibling}/include/brepslicer/Slicer.h")
    set(_brepslicer_default "${_brepslicer_sibling}")
else()
    set(_brepslicer_default "${_brepslicer_in_tree}")
endif()

set(BREPSLICER_CORE_DIR "${_brepslicer_default}" CACHE PATH
    "Path to BrepSlicer-core (STEP B-rep slicer, L0/L1/L2)")

if(NOT EXISTS "${BREPSLICER_CORE_DIR}/include/brepslicer/Slicer.h")
    message(FATAL_ERROR
        "BrepSlicer-core not found at ${BREPSLICER_CORE_DIR}. "
        "Expected vendored tree at src/brepslicer/core, or set -DBREPSLICER_CORE_DIR=...")
endif()

message(STATUS "BrepSlicer-core (L0/L1/L2): ${BREPSLICER_CORE_DIR}")

set(OpenCASCADE_DIR "${CMAKE_PREFIX_PATH}/lib/cmake/occt")
find_package(OpenCASCADE REQUIRED)

set(brepslicer_sources
    ${BREPSLICER_CORE_DIR}/src/io/JsonWriter.cpp
    ${BREPSLICER_CORE_DIR}/src/io/SvgWriter.cpp
    ${BREPSLICER_CORE_DIR}/src/io/DxfWriter.cpp
    ${BREPSLICER_CORE_DIR}/src/geom/GeomUtil.cpp
    ${BREPSLICER_CORE_DIR}/src/topo/FaceIndex.cpp
    ${BREPSLICER_CORE_DIR}/src/intersect/QuadricPlane.cpp
    ${BREPSLICER_CORE_DIR}/src/intersect/AnalyticIntersect.cpp
    ${BREPSLICER_CORE_DIR}/src/intersect/BSplineFit.cpp
    ${BREPSLICER_CORE_DIR}/src/intersect/UvMarch.cpp
    ${BREPSLICER_CORE_DIR}/src/intersect/UvMatch.cpp
    ${BREPSLICER_CORE_DIR}/src/contour/ContourAssembler.cpp
    ${BREPSLICER_CORE_DIR}/src/slicer/Slicer.cpp
    ${BREPSLICER_CORE_DIR}/src/verify/SliceVerify.cpp
    ${BREPSLICER_CORE_DIR}/src/occ/OccShape.cpp
    ${BREPSLICER_CORE_DIR}/src/occ/OccKernel.cpp
    ${BREPSLICER_CORE_DIR}/src/occ/OccBuilder.cpp
    ${BREPSLICER_CORE_DIR}/src/occ/OccSectionRef.cpp
    ${BREPSLICER_CORE_DIR}/src/occ/OccStrategy.cpp
)

add_library(brepslicer STATIC ${brepslicer_sources})
target_include_directories(brepslicer
    PUBLIC
        ${BREPSLICER_CORE_DIR}/include
        ${OpenCASCADE_INCLUDE_DIR}
    PRIVATE
        ${BREPSLICER_CORE_DIR}/src
)
# Expose private headers (BSplineFit) to adapters that sample L2 segments.
target_include_directories(brepslicer INTERFACE ${BREPSLICER_CORE_DIR}/src)
target_compile_features(brepslicer PUBLIC cxx_std_17)
target_compile_definitions(brepslicer PRIVATE NOMINMAX _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES)
if(MSVC)
    target_compile_options(brepslicer PRIVATE /utf-8 /bigobj)
endif()
